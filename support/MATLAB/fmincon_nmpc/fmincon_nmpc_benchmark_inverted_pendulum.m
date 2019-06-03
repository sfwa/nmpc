% Inverted pendulum control problem

% Set up optimal control problem.
[state_horizon, control_horizon, process_fcn, cost_fcn, lb, ub, constr_eq_fcn, constr_bound_fcn] = bench_ocp_inverted_pendulum();

% Solve the OCP using fmincon.
opts = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'MaxFunctionEvaluations', inf, ...
    'MaxIterations', 200, ...
    'PlotFcn', {'optimplotfval', 'optimplotfunccount'});

% Pack state and control into augmented state vector.
z = reshape([state_horizon; control_horizon], [], 1);
horizon_length = size(state_horizon, 2)-1;

[z_out, fval, exitflag, output] = fmincon(...
    @(x) cost_fcn(x, 1:horizon_length+1), ...
    z, ...
    [], [], ...
    [], [], ...
    lb, ub, ...
    @(z) fmincon_constraint_fcn(z, 5, 1, process_fcn, constr_eq_fcn, constr_bound_fcn), ...
    opts);

z_out = reshape(z_out, 6, []);
x_out = z_out(1:5, :);
u_out = z_out(6, :);

hFig = figure;
ax = gca;

for ii = 1:horizon_length
    plot([x_out(1, ii) x_out(1, ii) - 0.7*sin(x_out(2, ii))], [0 0.7*cos(x_out(2, ii))], '-o', 'LineWidth', 2);
    xlim([-1 1]); ylim([-1 1]); axis equal; grid on; grid minor;
    pause(0.05);
end

subplot(3, 1, 1)
plot(x_out(1, :))
subplot(3, 1, 2)
plot(x_out(2, :))
subplot(3, 1, 3)
plot(u_out)
