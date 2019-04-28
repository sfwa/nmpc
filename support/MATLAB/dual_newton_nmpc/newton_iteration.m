% Do a single unconstrained Newton iteration of the OCP.
function [x, u, lambda, epsilon, fStar, H, alpha] = newton_iteration(x, u, lambda, ...
        process_fcn, cost_fcn, lb, ub, constr_eq_fcn, constr_bound_fcn)
    % MATLAB's one-based indexing is a real pain here for trying to be
    % consistent with the notation used in the paper.
    N = size(x, 2)-1;

    nx = size(x, 1);

    assert(size(u, 2) == N+1, 'Inconsistent horizon length');
    assert(size(lambda, 2) == N+2, 'Inconsistent horizon length');
    assert(size(lambda, 1) == nx, 'Incorrect dual variable size');

    % Newton Hessian.
    H = zeros(nx*N, nx*N);

    act_tol = 1e-6; % Tolerance for active constraints.

    [H_k, g_k, A, b, Aeq, beq, E_k, C_k, c_k, lb, ub] = setup_all_stage_qps(x, u, ...
        lb, ub, process_fcn, cost_fcn, constr_eq_fcn, constr_bound_fcn);

    % Solve stage QPs to get incumbent objective function value and initial
    % primal estimate.
    [z, dz, p, q, active_set, D_k] = solve_all_stage_qps(x, u, lambda, H_k, g_k, ...
        E_k, C_k, c_k, lb, ub, A, b, Aeq, beq, act_tol);
    fStar_inc = cost_fcn(reshape(z, [], 1), 1:N+1) + dot(reshape(dz, [], 1), reshape(p, [], 1)) + sum(q);

    % Calculate newton gradient.
    g = calculate_newton_gradient(nx, N, dz, E_k, C_k, c_k);
    
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
        Z_k = Q_k(:, n_act+1:end);
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
    dLambda = reverse_cholesky(nx, N, H, -g, 1e-10, 1e-6);

    % Calculate the step size via an Armijo backtracking line search. Need
    % to look at each stage QP and find the distance to the closest active
    % set boundary and use that as the minimum step length.
    alphaMax = 1;
%     alphaMin = 0;
    
    alphaScale = 0.5;
    while true
        % Calculate candidate objective function value for the current step
        % size.
        lambda_cand = reshape(lambda, [], 1);
        lambda_cand(nx+1:end-nx) = lambda_cand(nx+1:end-nx) + alphaMax * dLambda;
        lambda_cand = reshape(lambda_cand, nx, []);

        [z, dz, p, q, ~] = solve_all_stage_qps(x, u, lambda_cand, H_k, g_k, ...
            E_k, C_k, c_k, lb, ub, A, b, Aeq, beq, act_tol);
        fStar_cand = cost_fcn(reshape(z, [], 1), 1:N+1) + dot(reshape(dz, [], 1), reshape(p, [], 1)) + sum(q);
        
        % Calculate Newton gradient for current dual variables.
        g = calculate_newton_gradient(nx, N, dz, E_k, C_k, c_k);
        
        sigma = 0.9;

        % Terminate line search at the appropriate point.
        if fStar_cand >= (fStar_inc + sigma * transpose(g) * dLambda)
            break;
        end

        alphaMax = alphaMax * alphaScale;
    end
    
    alpha = alphaMax;

    % Make the step.
    lambda = reshape(lambda, [], 1);
    lambda(nx+1:end-nx) = lambda(nx+1:end-nx) + alpha * dLambda;
    lambda = reshape(lambda, nx, []);
    
    fStar = fStar_cand;
    
    % Compute the primal infeasibility from the Newton gradient.
    epsilon = norm(g);

    x = z(1:nx, :);
    u = z(nx+1:end, :);
end

% Set up stage QPs. This could be optionally parallelised.
%
% Ideally with the active-set solver, the stage Hessians would be
% calculated once during initialisation using finite differences, and
% then kept up to date using L-BFGS or a Gauss-Newton approximation for
% least-squares objectives.
function [H_k, g_k, A, b, Aeq, beq, E_k, C_k, c_k, lb, ub] = setup_all_stage_qps(x, u, lb, ub, ...
        process_fcn, cost_fcn, constr_eq_fcn, constr_bound_fcn)
    N = size(x, 2)-1;
    nx = size(x, 1);
    nu = size(u, 1);
    nz = nx + nu;

    H_k = cell(N+1, 1);
    g_k = cell(N+1, 1);
    A = cell(N+1, 1);
    b = cell(N+1, 1);
    Aeq = cell(N+1, 1);
    beq = cell(N+1, 1);
    E_k = cell(N+1, 1);
    C_k = cell(N+1, 1);
    c_k = cell(N+1, 1);
    lb = reshape(lb, nz, N+1);
    ub = reshape(ub, nz, N+1);
    
    for kk = 0:N
        ii = kk + 1;
        
        if kk < N
            x_k_1 = x(:, ii+1);
        else
            x_k_1 = process_fcn([x(:, ii); u(:, ii)]);
        end

        [E_k{ii}, C_k{ii}, c_k{ii}, H_k{ii}, g_k{ii}, A{ii}, b{ii}, Aeq{ii}, beq{ii}, lb(:, ii), ub(:, ii)] = setup_stage_qp(...
            x(:, ii), x_k_1, u(:, ii), lb(:, ii), ub(:, ii), process_fcn, @(x) cost_fcn(x, ii), constr_eq_fcn, constr_bound_fcn);
        
        % Special cases.
        if kk == 0
            E_k{ii} = zeros(nx, nz);
        end
        
        if kk == N
            C_k{ii} = zeros(nx, nz);
        end
    end
end

% Set up a stage QP.
function [E_k, C_k, c_k, H_k, g_k, A, b, Aeq, beq, lb, ub] = setup_stage_qp(x_k, x_k_1, u_k, lb, ub, ...
        process_fcn, cost_fcn, constr_eq_fcn, constr_bound_fcn)
    z_k = [x_k; u_k];
    
    C_k = estimate_jacobian(process_fcn, z_k);
    c_k = process_fcn(z_k) - x_k_1;
    
    E_k = [eye(size(x_k, 1)) zeros(size(x_k, 1), size(u_k, 1))];

    % Could solve unconstrained problem here to get a point to linearise
    % the constraints around?
    
    % Linearise nonlinear constraints.
    if ~isempty(constr_eq_fcn)
        Aeq = estimate_jacobian(constr_eq_fcn, z_k);
        beq = -constr_eq_fcn(zeros(size(z_k)));
    else
        Aeq = [];
        beq = [];
    end
    
    if ~isempty(constr_bound_fcn)
        A = estimate_jacobian(constr_bound_fcn, z_k);
        b = -constr_bound_fcn(zeros(size(z_k)));
    else
        A = [];
        b = [];
    end
    
    % Calculate lower and upper bounds.
    lb = lb - z_k;
    ub = ub - z_k;
    
    % Calculate linearised cost function.
    H_k = estimate_hessian(cost_fcn, z_k);
    g_k = transpose(estimate_jacobian(cost_fcn, z_k));
end

% Solve all stage QPs and return the objective value.
function [z, dz, p, q, active_set, D_k] = solve_all_stage_qps(x, u, lambda, H_k, g_k, ...
        E_k, C_k, c_k, lb, ub, A, b, Aeq, beq, act_tol)
    N = size(x, 2)-1;
    D_k = cell(N+1, 1);
    z = zeros(size(x, 1) + size(u, 1), N+1);
    dz = zeros(size(x, 1) + size(u, 1), N+1);
    p = zeros(size(x, 1) + size(u, 1), N+1);
    q = zeros(1, N+1);

    % Solve stage QPs. This could be optionally parallelised.
    %
    % One-based indexing messes everything up here, so care needs to be taken
    % to avoid off-by-one errors in the indices.
    active_set = cell(N+1, 1);
    for kk = 0:N
        ii = kk + 1;

        [z(:, ii), dz(:, ii), p(:, ii), q(:, ii), active_set{ii}, D_k{ii}] = ...
        solve_stage_qp([x(:, ii); u(:, ii)], lambda(:, ii), lambda(:, ii+1), ...
            H_k{ii}, g_k{ii}, E_k{ii}, C_k{ii}, c_k{ii}, ...
            lb(:, ii), ub(:, ii), A{ii}, b{ii}, Aeq{ii}, beq{ii}, act_tol);
    end
end

% Solve a stage QP.

% In the future, write an online active-set solver which can take advantage
% of warm-starting when doing RTI. For now, just use the MATLAB 'fmincon'
% function. Also, if the cost function is really going to be nonlinear,
% probably want to use a Quasi-Newton method for online estimation of the
% Hessian.
function [z_k, dz_k, p_k, q_k, active_set, D_k] = solve_stage_qp(...
        z_k, lambda_k, lambda_k_1, H_k, g_k, ...
        E_k, C_k, c_k, lb, ub, A, b, Aeq, beq, act_tol)
    
    % Update p_k and q_k with latest dual estimate.
    p_k = transpose([-E_k; C_k]) * [lambda_k; lambda_k_1];
    q_k = transpose([zeros(size(lambda_k, 1), 1); c_k]) * [lambda_k; lambda_k_1];
    
    % Solve the QP.
%     dz_k = -pinv(H_k)*(g_k + p_k);
%     active_bounds = ((lb - z_k) > act_tol) | ((dz_k - ub) > act_tol);
%     dz_k = max(min(dz_k, ub), lb);
    warning('off', 'optim:quadprog:WillBeRemoved');
    opts = optimoptions('quadprog', ...
        'Algorithm', 'active-set', ...
        'ConstraintTolerance', act_tol, ...
        'Display', 'off');
    [dz_k, ~, ~, ~, lagrange] = quadprog(H_k, g_k + p_k, ...
        A, b, Aeq, beq, lb, ub, zeros(size(z_k)), opts);
    warning('on', 'optim:quadprog:WillBeRemoved');
    
    % Calculate mu_k from the Lagrange multipliers so that it is in the same
    % order as the constraints in D_k.
    active_bounds = (abs(lagrange.lower) > act_tol) | (abs(lagrange.upper) > act_tol);
    D_k_bounds = eye(numel(lb));
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
    
    z_k = z_k + dz_k;
end

function g = calculate_newton_gradient(nx, N, z, E_k, C_k, c_k)
    % Newton Hessian and gradient calculation.
    g = zeros(nx, N);
    
    for kk = 0:N
        ii = kk + 1;

        % Set up the Newton gradient for this stage. Note that the negative
        % sign in front of the right hand side of equation (6) in the
        % qpDUNES paper is erroneous.
        grad_block = -([-E_k{ii}; C_k{ii}] * z(:, ii) + [zeros(nx, 1); c_k{ii}]);

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
function [x, R] = reverse_cholesky(nx, N, H, g, regTol, regEps)
    R = zeros(size(H));
    
    for kk = fliplr(1:N)
        imax = max(1, (kk - 2)*nx + 1);
        
        for jj = fliplr(((kk - 1)*nx + 1):kk*nx)
            w = H(jj, jj);
            
            lmax = min(N*nx, (kk + 1)*nx);

            for ll = jj+1:lmax
                w = w - R(jj, ll)^2;
            end

            % Apply on-the-fly regularisation here if necessary.
            if w < regTol
                w = w + regEps;
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

% Simple function to estimate a Jacobian using central differences.
function J_est = estimate_jacobian(fcn, x, h)
    if nargin < 3
        h = 1e-6;
    end

    N = size(x, 1);
    J_est = zeros(size(fcn(x), 1), N);
    for ii = 1:N
        dx = zeros(size(x));
        dx(ii) = h;
        J_est(:, ii) = (fcn(x + dx) - fcn(x - dx)) / (2*h);
    end
end

% Simple function to estimate a Hessian using central differences.
function H_est = estimate_hessian(fcn, x, h)
    if nargin < 3
        h = 1e-3;
    end
    
    N = size(x, 1);
    H_est = zeros(N);
    for ii = 1:N
        dx = zeros(size(x));
        dx(ii) = h;
        H_est(ii, :) = (estimate_jacobian(fcn, x + dx, h) - estimate_jacobian(fcn, x - dx, h)) / (2*h);
    end
    
    % Make the estimate symmetric.
    H_est = 0.5 * (H_est + transpose(H_est));
end
