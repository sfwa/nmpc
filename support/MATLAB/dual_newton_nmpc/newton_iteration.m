% Do a single unconstrained Newton iteration of the OCP.
function [x, u, lambda, epsilon, fStar, H, alpha] = newton_iteration(x, u, lambda, ...
        process_fcn, cost_fcn, lb, ub, constr_eq_fcn, constr_bound_fcn)
    % MATLAB's one-based indexing is a real pain here for trying to be
    % consistent with the notation used in the paper.
    N = size(x, 2)-1;

    nx = size(x, 1);
    nu = size(u, 1);
    nz = nx + nu;

    assert(size(u, 2) == N+1, 'Inconsistent horizon length');
    assert(size(lambda, 2) == N+2, 'Inconsistent horizon length');
    assert(size(lambda, 1) == nx, 'Incorrect dual variable size');

    % Newton Hessian.
    H = zeros(nx*N, nx*N);

    act_tol = 1e-6; % Tolerance for active constraints.

    % Set up stage QPs. This could be optionally parallelised.
    %
    % Ideally with the active-set solver, the stage Hessians would be
    % calculated once during initialisation using finite differences, and
    % then kept up to date using L-BFGS or a Gauss-Newton approximation for
    % least-squares objectives.
    E_k = cell(N+1, 1);
    C_k = cell(N+1, 1);
    c_k = cell(N+1, 1);
    lb = reshape(lb, nz, N+1);
    ub = reshape(ub, nz, N+1);
    
    % Set up first stage initial value constraint.
    lb(1:nx, 1) = x(:, 1);
    ub(1:nx, 1) = x(:, 1);
    
    for kk = 0:N
        ii = kk + 1;

        [E_k{ii}, C_k{ii}, c_k{ii}] = setup_stage_qp(x(:, ii), u(:, ii), process_fcn);
        
        % Special cases.
        if kk == 0
            E_k{ii} = zeros(nx, nz);
        end
        
        if kk == N
            C_k{ii} = zeros(nx, nz);
        end
    end

    % Solve stage QPs to get imcumbent objective function value and initial
    % primal estimate.
    [z, active_set, fStar_inc, H_k, g_k, D_k] = solve_all_stage_qps(x, u, lambda, ...
        E_k, C_k, c_k, lb, ub, cost_fcn, constr_eq_fcn, constr_bound_fcn, act_tol);

    % Calculate newton gradient.
    g = calculate_newton_gradient(nx, N, z, E_k, C_k, c_k);
    
    % Newton Hessian calculation.
    for kk = 0:N
        ii = kk + 1;

        % Set up the Newton Hessian. This involves calculating a null-space
        % basis matrix Z_k. When the active-set stage QP solver is
        % implemented, this can be retrieved directly from that at no
        % additional cost.
        %
        % Meanwhile, calculate the nullspace basis matrix using a QR
        % decomposition, and then calculate the Z_k and P_k matrices as per
        % III-C in the qpDUNES paper.
        %
        % Note: An optimisation is that the null-space basis only needs to be
        % re-calculated if there's an active set change. Also, if a constant
        % Hessian is used and no active set change occurred, this step can be
        % skipped entirely.
        %
        % Use the pseudo-inverse when calculating the reduced Hessian in
        % order to be able to handle a positive semidefinite objective.
        active_set_k = active_set{ii};
        n_act = sum(active_set_k);
        D_k_act = D_k{ii}(active_set_k, :);
        [Q_k, ~] = qr(D_k_act');
        Z_k = Q_k(:, n_act+1:end); % Might need to transpose Q_k?
        P_k = Z_k * pinv(transpose(Z_k) * H_k{ii} * Z_k) * transpose(Z_k);

        % Diagonal part.
        if kk > 0
            span = 1:nx;
            H((kk-1)*nx + span, (kk-1)*nx + span) = ...
                C_k{ii-1} * P_km1 * transpose(C_k{ii-1}) + ...
                E_k{ii} * P_k * transpose(E_k{ii});

            % Sub-diagonal part.
            if kk < N
                H(kk*nx + span, (kk-1)*nx + span) = -C_k{ii} * P_k * transpose(E_k{ii});
                H((kk-1)*nx + span, kk*nx + span) = transpose(...
                    -C_k{ii} * P_k * transpose(E_k{ii}));
            end
        end

        % Store P_k for next iteration.
        P_km1 = P_k;
    end

    % Solve the Newton system using banded Cholesky factorisation, to yield
    % the step direction. For the C++ implementation, want to consider adding
    % on-the-fly regularisation.
    %
    % Use the algorithm on page 13 of "A Parallel Quadratic Programming
    % Method for Dynamic Optimization Problems" by Frasch et al.
    %
    % Alternatively, consider using the conjugate gradient method as
    % described in "An Improved Distributed Dual Newton-CG Method for
    % Convex Quadratic Programming Problems" by Kozma et al.
%     dLambda = mldivide(H, g);
    dLambda = reverse_cholesky(nx, N, H, g);

    % Calculate the step size via backtracking line search followed by
    % bisection for refinement. Need to look at each stage QP and find the
    % distance to the closest active set boundary and use that as the minimum
    % step length.
    alpha = 1;
    alphaMax = 1;
    alphaMin = 0;

    % Backtracking line search. In each iteration, check whether the
    % candidate step length is shorter than the active-set boundary of the
    % initial dual estimate, in which case there's no point trying any
    % shorter step.
    nMaxLineSearch = 20;
    alphaScale = 0.8;
    for ii = 1:nMaxLineSearch
        % Calculate candidate objective function value for the current step
        % size.
        lambda_cand = reshape(lambda, [], 1);
        lambda_cand(nx+1:end-nx) = lambda_cand(nx+1:end-nx) + alphaMax * dLambda;
        lambda_cand = reshape(lambda_cand, nx, []);

        [~, ~, fStar_cand] = solve_all_stage_qps(x, u, lambda_cand, ...
            E_k, C_k, c_k, lb, ub, cost_fcn, constr_eq_fcn, constr_bound_fcn, act_tol);

        % Terminate line search at the appropriate point.
        if fStar_cand > fStar_inc
            alphaMax = min(alphaMax / alphaScale, 1);
            break;
        end

        alphaMax = alphaMax * alphaScale;
    end
    
    % Only do the bisection search if the line search concluded; otherwise
    % just take the full step.
    if fStar_cand > fStar_inc
        % Bisection interval search. If the backtracking line search has found
        % the interval where the optimal objective lies, this step refines the
        % step length further within that interval.
        nMaxIntervalSearch = 50;
        bisectionTolerance = 1e-6;
        for ii = 1:nMaxIntervalSearch
            alpha = 0.5 * (alphaMax + alphaMin);

            % Calculate candidate objective function value for the current step
            % size.
            lambda_cand = reshape(lambda, [], 1);
            lambda_cand(nx+1:end-nx) = lambda_cand(nx+1:end-nx) + alpha * dLambda;
            lambda_cand = reshape(lambda_cand, nx, []);

            [z, ~, fStar_cand] = solve_all_stage_qps(x, u, lambda_cand, ...
                E_k, C_k, c_k, lb, ub, cost_fcn, ...
                constr_eq_fcn, constr_bound_fcn, act_tol);

            % Calculate Newton gradient for current dual variables.
            g = calculate_newton_gradient(nx, N, z, E_k, C_k, c_k);

            % Work out the Newton objective gradient.
            fDash = transpose(dLambda) * g;

            if abs(fDash) < bisectionTolerance
                break;
            elseif fDash < 0
                alphaMax = alpha;
            elseif fDash > 0
                alphaMin = alpha;
            end
        end
    end

    % Make the step.
    lambda = reshape(lambda, [], 1);
    lambda(nx+1:end-nx) = lambda(nx+1:end-nx) + alpha * dLambda;
    lambda = reshape(lambda, nx, []);
    
    % Re-solve stage QPs for a final time.
    [z, ~, fStar] = solve_all_stage_qps(x, u, lambda, ...
        E_k, C_k, c_k, lb, ub, cost_fcn, constr_eq_fcn, constr_bound_fcn, act_tol);
    
    % Compute the primal infeasibility from the Newton gradient.
    g = calculate_newton_gradient(nx, N, z, E_k, C_k, c_k);
    epsilon = norm(g);

    x = z(1:nx, :);
    u = z(nx+1:end, :);
end

% Set up a stage QP.
function [E_k, C_k, c_k] = setup_stage_qp(x_k, u_k, process_fcn)
    z_k = [x_k; u_k];
    
    C_k = jacobianest(process_fcn, z_k);
    c_k = process_fcn(zeros(size(z_k, 1), 1));
    
    E_k = [eye(size(x_k, 1)) zeros(size(x_k, 1), size(u_k, 1))];
end

% Solve all stage QPs and return the objective value.
function [z, active_set, fStar, H_k, g_k, D_k] = solve_all_stage_qps(x, u, lambda, ...
        E_k, C_k, c_k, lb, ub, cost_fcn, constr_eq_fcn, constr_bound_fcn, act_tol)
    fStar = 0;
    N = size(x, 2)-1;
    H_k = cell(N+1, 1);
    g_k = cell(N+1, 1);
    D_k = cell(N+1, 1);
    z = zeros(size(x, 1) + size(u, 1), N+1);

    % Solve stage QPs. This could be optionally parallelised.
    %
    % One-based indexing messes everything up here, so care needs to be taken
    % to avoid off-by-one errors in the indices.
    active_set = cell(N+1, 1);
    for kk = 0:N
        ii = kk + 1;

        [z(:, ii), active_set{ii}, fStar_k, H_k{ii}, g_k{ii}, D_k{ii}] = ...
        solve_stage_qp(x(:, ii), u(:, ii), ...
            lambda(:, ii), lambda(:, ii+1), E_k{ii}, C_k{ii}, c_k{ii}, ...
            lb(:, ii), ub(:, ii), cost_fcn, constr_eq_fcn, constr_bound_fcn, act_tol);

        % Add stage objective value to total objective value.
        fStar = fStar + fStar_k;
    end
end

% Solve a stage QP.

% In the future, write an online active-set solver which can take advantage
% of warm-starting when doing RTI. For now, just use the MATLAB 'fmincon'
% function. Also, if the cost function is really going to be nonlinear,
% probably want to use a Quasi-Newton method for online estimation of the
% Hessian.
function [z_k, active_set, fStar_k, H_k, g_k, D_k] = solve_stage_qp(...
        x_k, u_k, lambda_k, lambda_k_1, ...
        E_k, C_k, c_k, lb, ub, cost_fcn, constr_eq_fcn, constr_bound_fcn, act_tol)
    z_k = [x_k; u_k];
    
    % Update p_k and q_k with latest dual estimate.
    p_k = transpose([-E_k; C_k]) * [lambda_k; lambda_k_1];
    q_k = transpose([zeros(size(x_k, 1), 1); c_k]) * [lambda_k; lambda_k_1];
    
    % Could solve unconstrained problem here to get a point to linearise
    % the constraints around.
    
    % Linearise nonlinear constraints.
    if ~isempty(constr_eq_fcn)
        Aeq = jacobianest(constr_eq_fcn, z_k);
        beq = -constr_eq_fcn(zeros(size(z_k)));
    else
        Aeq = [];
        beq = [];
    end
    
    if ~isempty(constr_bound_fcn)
        A = jacobianest(constr_bound_fcn, z_k);
        b = -constr_bound_fcn(zeros(size(z_k)));
    else
        A = [];
        b = [];
    end
    
    % Calculate linearised cost function.
    H_k = hessian(cost_fcn, z_k);
    g_k = transpose(gradest(cost_fcn, zeros(size(z_k))));
    
    % Solve the QP.
    opts = optimoptions('quadprog', ...
        'Algorithm', 'interior-point-convex', ...
        'ConstraintTolerance', act_tol, ...
        'Display', 'off');
    [z_k, ~, ~, ~, lagrange] = quadprog(H_k, g_k + p_k, ...
        A, b, Aeq, beq, lb, ub, z_k, opts);
    
    % Calculate mu_k from the Lagrange multipliers so that it is in the same
    % order as the constraints in D_k.
    active_bounds = (abs(lagrange.lower) > act_tol) | (abs(lagrange.upper) > act_tol);
    D_k_bounds = eye(numel(lagrange.lower));
    active_set = active_bounds;
    D_k = D_k_bounds;

    if ~isempty(Aeq)
        active_set = [active_set; abs(lagrange.eqlin) > act_tol];
        D_k = [D_k; Aeq];
    end
    
    if ~isempty(A)
        active_set = [active_set; abs(lagrange.ineqlin) > act_tol];
        D_k = [D_k; A];
    end

    fStar_k = transpose(z_k) * H_k * z_k + transpose(g_k + p_k) * z_k + q_k;
end

function [c, ceq] = combined_constraints(z, constr_eq_fcn, constr_bound_fcn)
    if ~isempty(constr_bound_fcn)
        c = constr_bound_fcn(z);
    else
        c = [];
    end
    
    if ~isempty(constr_eq_fcn)
        ceq = constr_eq_fcn(z);
    else
        ceq = [];
    end
end

function g = calculate_newton_gradient(nx, N, z, E_k, C_k, c_k)
    % Newton Hessian and gradient calculation.
    g = zeros(nx, N);
    
    for kk = 0:N
        ii = kk + 1;

        % Set up the Newton gradient for this stage. Note that the negative
        % sign in front of the right hand side of equation (6) in the
        % qpDUNES paper is erroneous.
        grad_block = [-E_k{ii}; C_k{ii}] * z(:, ii) + [zeros(nx, 1); c_k{ii}];

        if kk > 0
            g(:, kk) = g(:, kk) + grad_block(1:nx);
        end

        if kk < N
            g(:, kk+1) = g(:, kk+1) + grad_block((nx + 1):end);
        end
    end
    
    g = reshape(g, [], 1);
end

% Solve the Newton system using banded Cholesky factorisation, to yield
% the step direction. Consider adding on-the-fly regularisation.
%
% Use the algorithm on page 13 of "A Parallel Quadratic Programming
% Method for Dynamic Optimization Problems" by Frasch et al.
function [x, R] = reverse_cholesky(nx, N, H, g)
    R = zeros(size(H));
    
    for kk = fliplr(1:N)
        imax = max(1, (kk - 2)*nx + 1);
        
        for jj = fliplr(((kk - 1)*nx + 1):kk*nx)
            w = H(jj, jj);
            
            lmax = min(N*nx, (kk + 1)*nx);

            for ll = jj+1:lmax
                w = w - R(jj, ll)^2;
            end

            R(jj, jj) = sqrt(w);

            for ii = fliplr(imax:jj-1)
                w = H(ii, jj);

                if ii > (kk - 1) * nx
                    lmax = min(N*nx, (kk + 1)*nx);
                else
                    lmax = min(N*nx, kk*nx);
                end

                for ll = jj+1:lmax
                    w = w - R(jj, ll) * R(ii, ll);
                end

                % Apply on-the-fly regularisation here if necessary.
                R(ii, jj) = w / R(jj, jj);
            end
        end
    end
    
    x = zeros(size(g));
    y = zeros(size(g));
    
    % Solve the R * y = g part using back substitution.
    for kk = fliplr(1:N)
        lmax = min(N*nx, (kk + 1)*nx);
        
        for jj = fliplr(((kk - 1)*nx + 1):kk*nx)
            w = g(jj);

            for ll = jj+1:lmax
                w = w - R(jj, ll) * y(ll);
            end
            
            y(jj) = w / R(jj, jj);
        end
    end
    
    % Solve the R^T * x = y part using forward substitution.
    for kk = 1:N
        lmax = max(1, (kk - 2)*nx + 1);
        
        for jj = ((kk - 1)*nx + 1):kk*nx
            w = y(jj);

            for ll = lmax:jj-1
                w = w - R(ll, jj) * x(ll);
            end
            
            x(jj) = w / R(jj, jj);
        end
    end
end
