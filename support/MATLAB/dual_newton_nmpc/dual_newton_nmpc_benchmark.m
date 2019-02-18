% Set up asteroids-style test scenario where the ship has to navigate to a
% specific point while avoiding asteroids.

% Set up OCP.
horizon_length = 100;
dt = 0.1;

% Initial conditions.
init_state = [
    10;         % X-position
    -10;        % Y-position
    0;          % X-velocity
    3;          % Y-velocity
    30*pi/180;  % Heading
];

init_control = [
    0;          % Angular rate
    0;          % Acceleration
];

control_min = [-180*pi/180; 0];
control_max = [180*pi/180; 10];
speed_max = 10;

% Set up optimal control problem.
[state_horizon, control_horizon, process_fcn, cost_fcn, lb, ub, constr_eq_fcn, constr_bound_fcn] = bench_ocp(...
    init_state, init_control, horizon_length, dt, control_min, control_max, speed_max);

% Set up initial dual vector.
lambda = zeros(size(state_horizon, 1), horizon_length+2);

% Number of Newton iterations to execute.
numIterations = 5;
for ii = 1:numIterations
    [state_horizon, control_horizon, lambda, epsilon, fStar] = newton_iteration(...
        state_horizon, control_horizon, lambda, ...
        process_fcn, cost_fcn, lb, ub, constr_eq_fcn, constr_bound_fcn);

    % Store information and set up for next iteration.
    fprintf('Iteration %d\n', ii);
end

hFig = figure; axis equal; hold on; grid on; grid minor;
xlim([min(min(x_out(2, :)) - 2, 0) max(max(x_out(2, :)) + 2, 0)]);
ylim([min(min(x_out(1, :)) - 2, 0) max(max(x_out(1, :)) + 2, 0)]);
ax = gca;

for ii = 1:horizon_length
    plot(x_out(2, 1:ii), x_out(1, 1:ii), 'LineWidth', 2);
    quiver(x_out(2, ii), x_out(1, ii), sin(x_out(5, ii)), cos(x_out(5, ii)), ...
        'MaxHeadSize', 2);
    ax.ColorOrderIndex = 1;
    pause(dt);
end
