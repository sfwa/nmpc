% Simple double-integrator optimal control problem.

% Set up OCP.
horizon_length = 20;
dt = 0.2;

% Initial conditions.
init_state = [
    10;         % Position
     5;         % Velocity
];

init_control = 0; % Acceleration

state_min = [0; -5];
state_max = [inf; 5];
control_min = -10;
control_max = 10;
speed_max = 10;

x_ref = zeros(size(init_state, 1), horizon_length+1);
u_ref = zeros(size(init_control, 1), horizon_length+1);

x_ref(1, :) = 1:21;

% Set up optimal control problem.
[state_horizon, control_horizon, process_fcn, cost_fcn, lb, ub, constr_eq_fcn, constr_bound_fcn] = bench_ocp_simple(...
    init_state, init_control, x_ref, u_ref, horizon_length, dt, state_min, state_max, ...
    control_min, control_max, speed_max);

% Solve the OCP using fmincon.
opts = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'MaxFunctionEvaluations', inf, ...
    'MaxIterations', 100, ...
    'PlotFcn', {'optimplotfval', 'optimplotfunccount'});

% Pack state and control into augmented state vector.
z = reshape([state_horizon; control_horizon], [], 1);

% Set up initial condition state equality constraint.
Aeq = zeros(numel(z));
Aeq(1:2, 1:2) = eye(2);
beq = zeros(numel(z), 1);
beq(1:2) = init_state(:, 1);

[z_out, fval, exitflag, output] = fmincon(...
    @(x) cost_fcn(x, 1:horizon_length+1), ...
    z, ...
    [], [], ...
    Aeq, beq, ...
    lb, ub, ...
    @(z) fmincon_constraint_fcn(z, 2, 1, process_fcn, constr_eq_fcn, constr_bound_fcn), ...
    opts);

z_out = reshape(z_out, 3, []);
x_out = z_out(1:2, :);
u_out = z_out(3, :);

subplot(3, 1, 1)
plot(x_out(1, :))
subplot(3, 1, 2)
plot(x_out(2, :))
subplot(3, 1, 3)
plot(u_out)
