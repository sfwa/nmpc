% Set up nonlinear constraint function as expected by fmincon().
function [c, ceq] = fmincon_constraint_fcn(z, process_fcn, constr_eq_fcn, constr_bound_fcn)
    % Set up non-linear equality constraints including continuity
    % constraints.
    horizon = reshape(z, 7, []);
    N = size(horizon, 2);
    x = horizon(1:5, :);
    u = horizon(6:7, :);
    ceq = [
        reshape(process_fcn(horizon(:, 1:end-1)) - x(:, 2:end), [], 1);
        reshape(constr_eq_fcn(horizon), [], 1);
    ];

    % Set up non-linear inequality constraints.
    c = reshape(constr_bound_fcn(horizon), [], 1);
end
