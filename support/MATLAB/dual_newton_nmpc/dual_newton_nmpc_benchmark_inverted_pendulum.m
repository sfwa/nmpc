% Inverted pendulum control problem

% Set up optimal control problem.
[state_horizon, control_horizon, process_fcn, cost_fcn, lb, ub, constr_eq_fcn, constr_bound_fcn] = bench_ocp_inverted_pendulum();

horizon_length = size(state_horizon, 2)-1;

% Set up initial dual vector.
lambda = zeros(size(state_horizon, 1), horizon_length+2);

% Number of Newton iterations to execute.
numIterations = 200;
for ii = 1:numIterations
    [state_horizon, control_horizon, lambda, epsilon, fStar, H, alpha] = newton_iteration(...
        state_horizon, control_horizon, lambda, ...
        process_fcn, cost_fcn, lb, ub, constr_eq_fcn, constr_bound_fcn);

    % Store information and set up for next iteration.
    fprintf('Iteration %d, primal: %.1f, dual: %.1f, rcond: %e, alpha: %e\n', ii, epsilon, fStar, rcond(H), alpha);
    
    if epsilon < 1e-4
        break;
    end
end

x_out = state_horizon;
u_out = control_horizon;

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
