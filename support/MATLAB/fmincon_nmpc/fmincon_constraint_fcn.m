% Set up nonlinear constraint function as expected by fmincon().
function [c, ceq] = fmincon_constraint_fcn(z, nx, nu, process_fcn, constr_eq_fcn, constr_bound_fcn)
    % Set up non-linear equality constraints including continuity
    % constraints.
    horizon = reshape(z, nx+nu, []);
    N = size(horizon, 2);
    x = horizon(1:nx, :);
    u = horizon(nx+1:nx+nu, :);
    ceq = [
        reshape(process_fcn(horizon(:, 1:end-1)) - x(:, 2:end), [], 1);
        reshape(constr_eq_fcn(horizon), [], 1);
    ];

    % Set up non-linear inequality constraints.
    c = reshape(constr_bound_fcn(horizon), [], 1);
end
