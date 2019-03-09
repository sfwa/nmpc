% Solve the optimal control problem using MATLAB's fmincon().
function [state_horizon, control_horizon, process_fcn, cost_fcn, lb, ub, constr_eq_fcn, constr_bound_fcn] = bench_ocp(...
        x0, u0, x_ref, u_ref, N, dt, x_min, x_max, u_min, u_max, speed_max)
    assert(iscolumn(x_min) & numel(x_min) == 5, ...
        'Wrong dimension for state lower bound');
    assert(iscolumn(x_max) & numel(x_max) == 5, ...
        'Wrong dimension for state upper bound');
    assert(iscolumn(u_min) & numel(u_min) == 2, ...
        'Wrong dimension for control lower bound');
    assert(iscolumn(u_max) & numel(u_max) == 2, ...
        'Wrong dimension for control upper bound');
    assert(isscalar(dt), 'Time step must be scalar');
    assert(isscalar(speed_max), 'Maximum speed must be scalar');

    state_horizon = zeros(5, N+1);
    state_horizon(:, 1) = x0;
    control_horizon = zeros(2, N+1);
    control_horizon(:, 1) = u0;

    % Propagate initial condition though state horizon using process
    % function.
    for ii = 2:N+1
        state_horizon(:, ii) = bench_process(dt, ...
            state_horizon(:, ii-1), control_horizon(:, ii-1));
    end

    % Set state vector lower and upper bounds (inc. control).
    lb = repmat([x_min; u_min], 1, N+1);
    ub = repmat([x_max; u_max], 1, N+1);
    
    % Set up initial state constraint.
    lb(1:size(state_horizon, 1), 1) = x0;
    ub(1:size(state_horizon, 1), 1) = x0;
    
    lb = reshape(lb, [], 1);
    ub = reshape(ub, [], 1);
    
    process_fcn = @(z) bench_process(dt, z(1:5, :), z(6:7, :));
    cost_fcn = @(z, ii) bench_objective(z - reshape([x_ref(:, ii); u_ref(:, ii)], [], 1));
    constr_eq_fcn = @(z) bench_eq_constraints(z);
    constr_bound_fcn = @(z) bench_bound_constraints(z, speed_max);
end

% Solve the set of IVPs over the horizon.
% For this, x_in is of dimension Nx by M, where Nx is the state dimension and
% M is the horizon length. Also, u is of dimension Nu by M, where Nu is the
% control dimension.
function x_out = bench_process(dt, x_in, u)
    assert(size(x_in, 1) == 5, 'Wrong state dimension');
    assert(size(u, 1) == 2, 'Wrong control dimension');
    assert(size(x_in, 2) == size(u, 2), 'Inconsistent horizon length');
    assert(isscalar(dt), 'Time step must be scalar');

    x_out = zeros(size(x_in));
    for ii = 1:size(x_in, 2)
        x_out(:, ii) = rk4_step(dt, x_in(:, ii), u(:, ii));
    end
end

function x_out = rk4_step(dt, x_in, u)
    k1 = bench_dynamics(x_in, u);
    k2 = bench_dynamics(x_in + k1*dt*0.5, u);
    k3 = bench_dynamics(x_in + k2*dt*0.5, u);
    k4 = bench_dynamics(x_in + k3*dt, u);

    x_out = x_in + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end

% System dynamics model. Takes a state vector and a control vector, and
% returns the state vector derivative.
function xdot = bench_dynamics(x, u)
    % State vector:
    %     - x(1): x-axis position
    %     - x(2): y-axis position
    %     - x(3): x-axis velocity
    %     - x(4): y-axis velocity
    %     - x(5): heading
    %
    % Control vector:
    %     - u(1): angular velocity
    %     - u(2): forward acceleration
    xdot = [...
        x(3:4);
        [cos(x(5)); sin(x(5))] * u(2);
        u(1);
    ];
end

% Set up objective function to minimise integrated displacement (i.e. try to
% fly to (0, 0)).
function obj_val = bench_objective(z)
    z = reshape(z, 7, []);
    obj_val = 0.5 * (...
        dot(z(1, :), z(1, :)) + ...
        dot(z(2, :), z(2, :)) + ...
        dot(z(3, :), z(3, :))*1e-3 + ...
        dot(z(4, :), z(4, :))*1e-3 + ...
        dot(z(5, :), z(5, :))*1e-6 + ...
        dot(z(6, :), z(6, :))*1e-6 + ...
        dot(z(7, :), z(7, :))*1e-6);
end

% Constraints for the optimal control problem.
function ceq = bench_eq_constraints(z)
    ceq = [];
end

function c = bench_bound_constraints(z, speed_max)
    x = z(1:5, :);

    % Set up non-linear maximum speed constraint.
    c = sqrt(x(3, :).^2 + x(4, :).^2) - speed_max;
end
