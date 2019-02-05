% System dynamics model. Takes a state vector and a control vector, and
% returns the state vector derivative.
function xdot = bench_dynamics(x, u)
    xdot = zeros(size(x));

    % State vector:
    %     - x(1): x-axis position
    %     - x(2): y-axis position
    %     - x(3): x-axis velocity
    %     - x(4): y-axis velocity
    %     - x(5): heading
    %
    % Control vector:
    %     - u(1): angular velocity
    %     - u(2): forward acceleration
    xdot(1:2) = x(3:4);
    xdot(3:4) = [cos(x(5)); sin(x(5))] * u(2);
    xdot(5) = u(1);
end
