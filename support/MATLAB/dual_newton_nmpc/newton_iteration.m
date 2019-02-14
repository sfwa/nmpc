% Do a single unconstrained Newton iteration of the OCP.
function [x, u, lambda, H, epsilon] = newton_iteration(x, u, lambda, ...
        process_fcn, cost_fcn, lb, ub, constr_eq_fcn, constr_bound_fcn)
    assert(size(x, 2) == size(u, 2) == size(lambda, 2), 'Inconsistent horizon length');
    assert(size(x, 1) + size(u, 1) == size(lambda, 1));

    N = size(x, 2);
    nx = size(x, 1);
    nu = size(u, 1);
    nz = nx + nu;

    % Newton Hessian and gradient.
    H = zeros(nz*(N-1), nz*(N-1));
    g = zeros(nz, (N-1));

    act_tol = 1e-6; % Tolerance for active constraints.

    % Set up stage QPs. This could be optionally parallelised.
    %
    % Ideally with the active-set solver, the stage Hessians would be
    % calculated once during initialisation using finite differences, and
    % then kept up to date using L-BFGS.
    E_k = [eye(nx) zeros(nx, nu)]; % E_k is the same for all stages.
    C_k = cell(N, 1);
    c_k = cell(N, 1);
    H_k = cell(N, 1);
    g_k = cell(N, 1);
    Aeq_k = cell(N, 1);
    beq_k = cell(N, 1);
    A_k = cell(N, 1);
    b_k = cell(N, 1);
    D_k = cell(N, 1);
    for kk = 1:N
        [C_k{kk}, c_k{kk}, H_k{kk}, g_k{kk}, Aeq_k{kk}, beq_k{kk}, A_k{kk}, b_k{kk}, D_k{kk}] = setup_stage_qp(...
            x(:, kk), u(:, kk), process_fcn, cost_fcn, constr_eq_fcn, constr_bound_fcn);
    end

    % Solve stage QPs. This could be optionally parallelised.
    %
    % One-based indexing messes everything up here, so care needs to be taken
    % to avoid off-by-one errors in the indices.
    for kk = 1:N
        if kk == 1
            lambda_k = zeros(size(lambda, 1), 1);
        else
            lambda_k = lambda(:, kk);
        end

        if kk == N
            lambda_k_1 = lambda(:, kk+1);
        else
            lambda_k_1 = zeros(size(lambda, 1), 1);
        end

        [z(:, kk), mu_k, ~, ~] = solve_stage_qp(x(:, kk), u(:, kk), lambda_k, lambda_k_1, ...
            E_k, C_k{kk}, c_k{kk}, H_k{kk}, g_k{kk}, Aeq_k{kk}, beq_k{kk}, A_k{kk}, b_k{kk}, lb, ub, act_tol);
    end

    % Initial value embedding here – first stage equality constraint.

    % Newton Hessian and gradient calculation.
    for kk = 1:N
        % Set up the Newton gradient for this stage.
        grad_block = -([-E_k; C_k] * z(:, kk) + [zeros(nz, 1); c_k]);

        if kk > 1
            g(:, kk-1) = g(:, kk-1) + grad_block(1:(nz));
        end

        if kk < N
            g(:, kk) = g(:, kk) + grad_block((nz + 1):end);
        end

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
        active_set = (mu_k <= act_tol);
        n_act = sum(active_set);
        D_k_act = D_k(active_set, :);
        Q = qr(D_k_act);
        B = Q(:, 1:n_act);
        N = Q(:, n_act+1:end);
        Z_k = [-inv(B)*N; eye(nz - n_act)];
        P_k = Z_k * inv(transpose(Z_k) * H_k * Z_k) * transpose(Z_k);

        % Diagonal part.
        span = 1:nz;
        if kk > 1
            H((kk-1)*nz + span, (kk-1)*nz + span) = C_km1 * P_km1 * transpose(C_km1) + E_k * P_k * transpose(E_k);
        end

        % Sub-diagonal part.
        H(, ) = -C_k * P_k * transpose(E_k);

        % Store P_k and C_k for next iteration.
        C_km1 = C_k;
        P_km1 = P_k;
    end

    % Compute the primal infeasibility from the Newton gradient.
    epsilon = norm(reshape(g, [], 1));

    % Solve the Newton system using banded Cholesky factorisation, to yield
    % the step direction. For the C++ implementation, want to consider adding
    % on-the-fly regularisation.
    %
    % Instrumentation: check positive definiteness of the Newton Hessian.
    dLambda = mldivide(H, g);

    % Calculate the step size via backtracking line search followed by
    % bisection for refinement. Need to look at each stage QP and find the
    % distance to the closest active set boundary and use that as the minimum
    % step length.

    % Return shifted dual variables ready for the next iteration.
end

% Set up a stage QP.

function [C_k, c_k, H_k, g_k, Aeq_k, beq_k, A_k, b_k, D_k] = setup_stage_qp(x_k, u_k, ...
        process_fcn, cost_fcn, constr_eq_fcn, constr_bound_fcn)
    C_k = jacobianest(process_fcn, z_k);
    c_k = process_fcn(zeros(size(z_k, 1), 1));

    % Estimate unconstrained Hessian for each stage.
    %
    % Could also use 'hessdiag' to estimate a diagonal Hessian in order to
    % use the clipping solver proposed in the qpDUNES paper.
    H_k = hessian(cost_fcn, z_k);
    g_k = gradest(cost_fcn, z_k);

    % Linearise constraints for each stage.
    if ~isempty(constr_eq_fcn)
        Aeq_k = jacobianest(constr_eq_fcn, z_k);
        beq_k = constr_eq_fcn(zeros(size(z_k, 1), 1));
    else
        Aeq_k = zeros(0, size(z_k, 1));
        beq_k = [];
    end

    if ~isempty(constr_bound_fcn)
        A_k = jacobianest(constr_bound_fcn, z_k);
        b_k = constr_bound_fcn(zeros(size(z_k, 1), 1));
    else
        A_k = zeros(0, size(z_k, 1));
        b_k = [];
    end

    % Construct linearised constraint matrix D_k. Note that the form of the
    % constraints here is a bit different to that in the qpDUNES paper; the
    % affine constraints with upper and lower bounds are more general than
    % these.
    D_k = [eye(size(z_k, 1)); eye(size(z_k, 1)); A_k; Aeq_k];
end

% Solve a stage QP.

% In the future, write an online active-set solver which can take advantage
% of warm-starting when doing RTI. For now, just use the MATLAB 'fmincon'
% function. Also, if the cost function is really going to be nonlinear,
% probably want to use a Quasi-Newton method for online estimation of the
% Hessian.
function [z_k, mu_k, p_k, q_k] = solve_stage_qp(x_k, u_k, lambda_k, lambda_k_1, ..
        E_k, C_k, c_k, H_k, g_k, Aeq_k, beq_k, A_k, b_k, lb, ub, act_tol)
    % Calculate linearised continuity constraints for each stage.
    z_k = [x_k; u_k];

    % Update p_k and q_k with latest dual estimate.
    p_k = g_k + transpose([-E_k; C_k]) * [lambda_k; lambda_k_1];
    q_k = transpose([zeros(size(x_k, 1), 1); c_k]) * [lambda_k; lambda_k_1];

    % Solve the QP.
    opts = optimoptions('quadprog', ...
        'Algorithm', 'interior-point-convex', ...
        'MaxIterations', 100, ...
        'ConstraintTolerance', act_tol, ...
        'LinearSolver', 'dense');
    [z_k, ~, ~, ~, lagrange] = quadprog(H_k, p_k, A_k, b_k, Aeq_k, beq_k, lb, ub, z_k, opts);

    % Calculate mu_k from the Lagrange multipliers so that it is in the same
    % order as the constraints in D_k.
    mu_k = [lagrange.lower; lagrange.upper; lagrange.ineqlin; lagrange.eqlin];
end

function objValue = calculate_objective_value()

end
