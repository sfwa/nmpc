% Solve the set of IVPs over the horizon.
% For this, x_in is of dimension Nx by M, where Nx is the state dimension and
% M is the horizon length. Also, u is of dimension Nu by M, where Nu is the
% control dimension.
function x_out = bench_ivp(dt, x_in, u)
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
