% Solve the optimal control problem using MATLAB's fmincon().
function [state_horizon, control_horizon, process_fcn, cost_fcn, lb, ub, constr_eq_fcn, constr_bound_fcn] = bench_ocp(...
        x0, u0, N, dt, u_min, u_max, speed_max)
    assert(iscolumn(u_min) & numel(u_min) == 2, ...
        'Wrong dimension for control lower bound');
    assert(iscolumn(u_max) & numel(u_max) == 2, ...
        'Wrong dimension for control upper bound');
    assert(isscalar(dt), 'Time step must be scalar');
    assert(isscalar(speed_max), 'Maximum speed must be scalar');

    state_horizon = zeros(5, N);
    state_horizon(:, 1) = x0;
    control_horizon = zeros(2, N);
    control_horizon(:, 1) = u0;

    % Propagate initial condition though state horizon using process
    % function.
    for ii = 2:N
        state_horizon(:, ii) = bench_process(dt, ...
            state_horizon(:, ii-1), control_horizon(:, ii-1));
    end

    % Set state vector lower and upper bounds (inc. control).
    lb = reshape(repmat([-inf(5, 1); u_min], 1, N), [], 1);
    ub = reshape(repmat([inf(5, 1); u_max], 1, N), [], 1);
    
    process_fcn = @(z) bench_process(dt, z(1:5, :), z(6:7, :));
    cost_fcn = @(z) bench_objective(z);
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
    x = reshape(z, 7, []);
    obj_val = sum(sqrt(x(1, :).^2 + x(2, :).^2));
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
