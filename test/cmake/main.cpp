// No Eigen include: NMPC headers must supply their own dependencies.
#include "NMPC/Solver.h"
#include <type_traits>

#ifdef real_t
#error "NMPC must not define real_t"
#endif
#ifdef _MSVC_LANG
static_assert(_MSVC_LANG >= 201703L, "NMPC must propagate C++17 to consumers");
#else
static_assert(__cplusplus >= 201703L, "NMPC must propagate C++17 to consumers");
#endif

#ifdef EXPECT_EIGEN_UNROLLING_LIMIT
static_assert(EIGEN_UNROLLING_LIMIT == EXPECT_EIGEN_UNROLLING_LIMIT,
    "NMPC_EIGEN_UNROLLING_LIMIT must reach consumers");
#endif

/* The header's version must be the CMake package's. */
constexpr bool same_string(const char* a, const char* b) {
    return *a == *b && (*a == '\0' || same_string(a + 1, b + 1));
}
static_assert(NMPC_VERSION == NMPC_VERSION_MAJOR * 10000 + NMPC_VERSION_MINOR * 100 + NMPC_VERSION_PATCH,
    "NMPC_VERSION must be made from its parts");
#ifdef EXPECT_VERSION_STRING
static_assert(same_string(NMPC_VERSION_STRING, EXPECT_VERSION_STRING),
    "NMPC_VERSION_STRING must match the version of the nmpc CMake package");
#endif

#ifdef EXPECT_DOUBLE
typedef double Scalar;
#else
typedef float Scalar;
#endif

/* The smallest problem: hold a scalar at zero. */
struct Hold {
    using Scalar = ::Scalar;
    enum { NX = 1, NU = 1, NR = 2, NT = 1, Stages = 5 };
    using State = Eigen::Matrix<Scalar, 1, 1>;
    using Control = Eigen::Matrix<Scalar, 1, 1>;
    static State retract(const State& x, const State& d) { return x + d; }
    static State difference(const State& x, const State& y) { return x - y; }
    State step(int, const State& x, const Control& u) const { return x + u; }
    void step_jacobian(int, const State&, const Control&, Eigen::Matrix<Scalar, 1, 1>& a,
            Eigen::Matrix<Scalar, 1, 1>& b) const {
        a(0) = Scalar(1);
        b(0) = Scalar(1);
    }
    void stage_residual(int, const State& x, const Control& u, Eigen::Matrix<Scalar, 2, 1>& r) const {
        r << x(0), u(0);
    }
    void residual_jacobians(int, const State&, const Control&, Eigen::Matrix<Scalar, 2, 1>& jx,
            Eigen::Matrix<Scalar, 2, 1>& ju) const {
        jx << Scalar(1), Scalar(0);
        ju << Scalar(0), Scalar(1);
    }
    void terminal_residual(const State& x, Eigen::Matrix<Scalar, 1, 1>& r) const { r = x; }
    void terminal_jacobian(const State&, Eigen::Matrix<Scalar, 1, 1>& j) const { j(0) = Scalar(1); }
    void control_bounds(int, Control& lower, Control& upper) const {
        lower(0) = Scalar(-1);
        upper(0) = Scalar(1);
    }
};

// Deliberately unused: the library must not export its warning policy.
int main(int argc, char**) {
    static NMPC::Solver<Hold>::Workspace workspace;
    NMPC::Solver<Hold> solver(workspace);
    Hold::State x;
    x(0) = Scalar(2);
    solver.reset(x, Hold::Control::Zero());
    const NMPC::Result<Scalar> result = solver.solve(Hold(), x, NMPC::Options<Scalar>());
    return result.usable() && result.cost_after < result.cost_before ? 0 : 1;
}
