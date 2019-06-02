function [state_horizon, control_horizon, process_fcn, cost_fcn, lb, ub, constr_eq_fcn, constr_bound_fcn] = bench_ocp_inverted_pendulum()
    N = 60;
    dt = 0.05;
    x0 = [0 pi 0 0 0].';
    x_ref = [0 0 0 0 0].';
    u_ref = 0;
    x_min = [-1 -inf -inf -inf -15].';
    x_max = [ 1  inf  inf  inf  15].';
    u_min = -30;
    u_max =  30;
    state_horizon = zeros(5, N+1);
    state_horizon(:, 1) = x0;
    control_horizon = zeros(1, N+1);

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
    
    process_fcn = @(z) bench_process(dt, z(1:5, :), z(6, :));
    cost_fcn = @(z, ii) bench_objective(z - repmat([x_ref; u_ref], N+1, 1));
    constr_eq_fcn = @(z) bench_eq_constraints(z);
    constr_bound_fcn = @(z) bench_bound_constraints(z);
end

% Solve the set of IVPs over the horizon.
% For this, x_in is of dimension Nx by M, where Nx is the state dimension and
% M is the horizon length. Also, u is of dimension Nu by M, where Nu is the
% control dimension.
function x_out = bench_process(dt, x_in, u)
    assert(size(x_in, 1) == 5, 'Wrong state dimension');
    assert(size(u, 1) == 1, 'Wrong control dimension');
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
    %     - x(1): cart position
    %     - x(2): pendulum angle
    %     - x(3): cart velocity
    %     - x(4): pendulum angular velocity
    %     - x(5): force on cart
    %
    % Control vector:
    %     - u(1): time derivative of force on cart
    
    % Solve the implicit ODE given in part IV of
    % https://cdn.syscop.de/publications/Kouzoupis2015.pdf for acceleration
    % and angular acceleration.
    m = 0.1;
    M = 1;
    l = 0.7;
    g = 9.81;
    lhs = [
          -(M + m) m*l*cos(x(2))
        -cos(x(2))             l
    ];
    rhs = [
        m*l*(x(4)^2)*sin(x(2)) - x(5)
        g*sin(x(2))
    ];
    implicit = lhs \ rhs;
    xdot = [...
        x(3);
        x(4);
        implicit(1);
        implicit(2);
        u(1);
    ];
end

% Set up objective function to minimise integrated displacement (i.e. try to
% get to 0).
function obj_val = bench_objective(z)
    nx = 5;
    nu = 1;
    z_k = reshape(z, nx+nu, []);
    N = size(z_k, 2);
    
    % Weights.
    Q = diag([0.5 1 0.001 0.001 0.001]);
    P = diag([0.5 1 0.001 0.001 0.001]);
    R = 0.001;
    
    W_k = mat2cell(repmat(blkdiag(Q, R), 1, 1, N-1), nx+nu, nx+nu, ones(N-1, 1));
    W = blkdiag(W_k{:}, blkdiag(P, 0));
    
    obj_val = z.' * W * z;
end

% Constraints for the optimal control problem.
function ceq = bench_eq_constraints(z)
    ceq = [];
end

function c = bench_bound_constraints(z)
    c = [];
end
