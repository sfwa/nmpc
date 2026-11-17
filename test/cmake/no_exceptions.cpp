// Instantiate the solver and the derivative checks with exceptions and RTTI disabled by the build.
#include "NMPC/Solver.h"
#include "NMPC/Check.h"

#if defined(__cpp_exceptions) || defined(__EXCEPTIONS)
#if !defined(_MSC_VER)
#error "This consumer must be compiled with exceptions disabled"
#endif
#endif

struct Integrator {
    using Scalar = double;
    enum { NX = 2, NU = 1, NR = 3, NT = 2, Stages = 10 };
    using State = Eigen::Matrix<double, 2, 1>;
    using Control = Eigen::Matrix<double, 1, 1>;
    static State retract(const State& x, const State& d) { return x + d; }
    static State difference(const State& x, const State& y) { return x - y; }
    State step(int, const State& x, const Control& u) const { return State(x(0) + 0.1 * x(1), x(1) + 0.1 * u(0)); }
    template <class A, class B>
    void step_jacobian(int, const State&, const Control&, A& a, B& b) const {
        a << 1.0, 0.1, 0.0, 1.0;
        b << 0.0, 0.1;
    }
    template <class R>
    void stage_residual(int, const State& x, const Control& u, R& r) const { r << x(0), x(1), 0.1 * u(0); }
    template <class JX, class JU>
    void residual_jacobians(int, const State&, const Control&, JX& jx, JU& ju) const {
        jx << 1.0, 0.0, 0.0, 1.0, 0.0, 0.0;
        ju << 0.0, 0.0, 0.1;
    }
    template <class R> void terminal_residual(const State& x, R& r) const { r = x; }
    template <class J> void terminal_jacobian(const State&, J& j) const { j.setIdentity(); }
    void control_bounds(int, Control& lower, Control& upper) const {
        lower(0) = -2.0;
        upper(0) = 2.0;
    }
};

int main() {
    static NMPC::Solver<Integrator>::Workspace workspace;
    NMPC::Solver<Integrator> solver(workspace);
    const Integrator problem;
    const Integrator::State x(1.0, 0.0);
    solver.reset(x, Integrator::Control::Zero());
    NMPC::Options<double> options;
    options.max_iterations = 5;
    const bool solved = solver.solve(problem, x, options).usable();
    const bool checked = NMPC::check_step_jacobian(problem, 0, x, Integrator::Control::Zero(), 1e-6).passed;
    return solved && checked ? 0 : 1;
}
