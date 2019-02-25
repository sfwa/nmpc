% Set up asteroids-style test scenario where the ship has to navigate to a
% specific point while avoiding asteroids.

% Set up OCP.
horizon_length = 20;
dt = 0.2;

% Initial conditions.
init_state = [
    10;         % Position
     5;         % Velocity
];

init_control = 0; % Acceleration

state_min = [-inf; -inf];
state_max = [inf; inf];
control_min = -10;
control_max = 10;
speed_max = [];

% Set up optimal control problem.
[state_horizon, control_horizon, process_fcn, cost_fcn, lb, ub, constr_eq_fcn, constr_bound_fcn] = bench_ocp_simple(...
    init_state, init_control, horizon_length, dt, state_min, state_max, ...
    control_min, control_max, speed_max);

% Set up initial dual vector.
lambda = zeros(size(state_horizon, 1), horizon_length+2);

% Number of Newton iterations to execute.
numIterations = 200;
for ii = 1:numIterations
    [state_horizon, control_horizon, lambda, epsilon, fStar] = newton_iteration(...
        state_horizon, control_horizon, lambda, ...
        process_fcn, cost_fcn, lb, ub, constr_eq_fcn, constr_bound_fcn);

    % Store information and set up for next iteration.
    fprintf('Iteration %d, primal: %.1f, dual: %.1f\n', ii, epsilon, fStar);
    
    if epsilon < 1e-6
        break;
    end
end
