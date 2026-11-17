/*
The smallest problem NMPC can solve: a point mass on a line that must come to
rest at the origin, with its acceleration limited. It is the example of the
README, and it shows every part of the problem concept in about fifty lines.

The state is position and velocity, the command is acceleration, and a stage
is 50 ms of exact double-integrator motion. The cost holds position, velocity
and acceleration small, and a terminal cost holds the state at the origin.

The program flies it closed loop the way a real controller does: once a
stage it shifts the plan, solves one iteration from the measured state, and
applies the first command.
*/

#include <cmath>
#include <cstdio>
#include "NMPC/Solver.h"

template <class T>
struct Stop_Problem {
    using Scalar = T;
    enum { NX = 2, NU = 1, NR = 3, NT = 2, Stages = 20 };
    using State = Eigen::Matrix<T, 2, 1>;      /* position (m), velocity (m/s) */
    using Control = Eigen::Matrix<T, 1, 1>;    /* acceleration (m/s^2) */

    T dt = T(0.05);                             /* stage length, s */
    T limit = T(4);                             /* acceleration limit, m/s^2 */

    static State retract(const State& x, const State& dx) { return x + dx; }
    static State difference(const State& x, const State& y) { return x - y; }

    State step(int, const State& x, const Control& u) const {
        return State(x(0) + dt * x(1) + T(0.5) * dt * dt * u(0), x(1) + dt * u(0));
    }
    template <class A, class B>
    void step_jacobian(int, const State&, const Control&, A& a, B& b) const {
        a << T(1), dt, T(0), T(1);
        b << T(0.5) * dt * dt, dt;
    }

    /* Residuals carry the root of the stage length, so the weights are per second. */
    template <class R>
    void stage_residual(int, const State& x, const Control& u, R& r) const {
        const T s = std::sqrt(dt);
        r << s * x(0), s * T(0.5) * x(1), s * T(0.1) * u(0);
    }
    template <class JX, class JU>
    void residual_jacobians(int, const State&, const Control&, JX& jx, JU& ju) const {
        const T s = std::sqrt(dt);
        jx << s, T(0), T(0), s * T(0.5), T(0), T(0);
        ju << T(0), T(0), s * T(0.1);
    }
    template <class R> void terminal_residual(const State& x, R& r) const { r = T(3) * x; }
    template <class JX> void terminal_jacobian(const State&, JX& jx) const { jx = T(3) * JX::Identity(); }

    void control_bounds(int, Control& lower, Control& upper) const {
        lower(0) = -limit;
        upper(0) = limit;
    }
};

int main() {
    using Problem = Stop_Problem<float>;
    using Solver = NMPC::Solver<Problem>;

    static Solver::Workspace workspace;         /* owned by the caller; nothing is allocated */
    Solver solver(workspace);
    Problem problem;
    NMPC::Options<float> options;               /* embedded defaults: 1 iteration, 5 candidates, 1 retry */

    Problem::State x(2.0f, 1.0f);
    solver.reset(x, Problem::Control::Zero());
    for (int stage = 0; stage < 100; stage++) {
        if (stage > 0) {
            solver.advance(problem, 1);         /* the plan moves on by the stage just flown */
        }
        const NMPC::Result<float> result = solver.solve(problem, x, options);
        if (!result.usable()) {
            std::printf("stage %d: no usable plan\n", stage);
            return 1;
        }
        const Problem::Control u = solver.policy().u[0];
        if (stage % 10 == 0) {
            std::printf("t %4.2f s  position %7.4f m  velocity %7.4f m/s  acceleration %6.3f m/s^2\n",
                double(stage) * 0.05, double(x(0)), double(x(1)), double(u(0)));
        }
        x = problem.step(0, x, u);
    }
    std::printf("final position %.5f m, velocity %.5f m/s\n", double(x(0)), double(x(1)));
    return std::abs(x(0)) < 1e-3f && std::abs(x(1)) < 1e-3f ? 0 : 1;
}
