% Set up asteroids-style test scenario where the ship has to navigate to a
% specific point while avoiding asteroids.

% Set up OCP.
horizon_length = 100;
dt = 0.1;
control_horizon = zeros(2, horizon_length);
state_horizon = zeros(5, horizon_length);

% Initial conditions.
state_horizon(:, 1) = [
    10;         % x-position
    -10;        % y-position
    0;          % x-velocity
    3;          % y-velocity
    30*pi/180;  % heading
];

% Propagate initial condition though state horizon.
for ii = 2:horizon_length
    state_horizon(:, ii) = bench_ivp(dt, ...
        state_horizon(:, ii-1), control_horizon(:, ii-1));
end

control_min = [-180*pi/180; 0];
control_max = [180*pi/180; 10];
speed_max = 10;

% Solve optimal control problem for each iteration.
[x_out, u_out, fval, exitflag, output] = bench_ocp(dt, ...
    state_horizon, control_horizon, control_min, control_max, speed_max);

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
