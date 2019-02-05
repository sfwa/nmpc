% Solve the optimal control problem using MATLAB's fmincon().
function [x_out, u_out, fval, exitflag, output] = bench_ocp(dt, x, u, u_min, u_max, speed_max)
    assert(size(x, 1) == 5, 'Wrong state dimension');
    assert(size(u, 1) == 2, 'Wrong control dimension');
    assert(size(x, 2) == size(u, 2), 'Inconsistent horizon length');
    assert(iscolumn(u_min) & numel(u_min) == 2, ...
        'Wrong dimension for control lower bound');
    assert(iscolumn(u_max) & numel(u_max) == 2, ...
        'Wrong dimension for control upper bound');
    assert(isscalar(dt), 'Time step must be scalar');
    assert(isscalar(speed_max), 'Maximum speed must be scalar');

    % Set up optimiser options.
    opts = optimoptions('fmincon', ...
        'Algorithm', 'interior-point', ...
        'MaxFunctionEvaluations', inf, ...
        'MaxIterations', 100, ...
        'Diagnostics', 'on', ...
        'PlotFcn', {'optimplotfval', 'optimplotfunccount'});

    % Set state vector lower and upper bounds (inc. control).
    lb = reshape(repmat([-inf(5, 1); u_min], 1, size(x, 2)), [], 1);
    ub = reshape(repmat([inf(5, 1); u_max], 1, size(x, 2)), [], 1);
    
    % Pack state and control into augmented state vector.
    z = reshape([x; u], [], 1);
    
    % Set up initial condition equality constraint.
    Aeq = zeros(numel(z));
    Aeq(1:size(x, 1), 1:size(x, 1)) = eye(size(x, 1));
    beq = zeros(size(z));
    beq(1:size(x, 1)) = x(:, 1);
    
    [z_out, fval, exitflag, output] = fmincon(...
        @bench_objective, ...
        z, ...
        [], [], ...
        Aeq, beq, ...
        lb, ub, ...
        @(x) bench_constraints(dt, speed_max, x), ...
        opts);
    
    z_out = reshape(z_out, 7, []);
    x_out = z_out(1:5, :);
    u_out = z_out(6:7, :);
end

% Set up objective function to minimise integrated displacement (i.e. try to
% fly to (0, 0)).
function obj_val = bench_objective(z)
    x = reshape(z, 7, []);
    obj_val = sum(sqrt(x(1, :).^2 + x(2, :).^2));
end

% Constraints for the optimal control problem.
function [c, ceq] = bench_constraints(dt, speed_max, z)
    horizon = reshape(z, 7, []);
    x = horizon(1:5, :);
    u = horizon(6:7, :);

    % Set up non-linear continuity constraints.
    ceq = reshape(bench_ivp(dt, x(:, 1:end-1), u(:, 1:end-1)) - x(:, 2:end), [], 1);

    % Set up non-linear maximum speed constraint.
    c = sqrt(x(3, :).^2 + x(4, :).^2) - speed_max;
end