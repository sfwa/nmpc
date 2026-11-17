// UKF and NMPC in one program: UKF's real_t is NMPC's scalar, and both see one
// Eigen.
#include <type_traits>

#include "NMPC/Solver.h"
#include "UKF/Core.h"

static_assert(std::is_same<real_t, float>::value,
              "Expected UKF in single precision");

struct Level {
  using Scalar = real_t;
  enum { NX = 1, NU = 1, NR = 2, NT = 1, Stages = 4 };
  using State = Eigen::Matrix<real_t, 1, 1>;
  using Control = Eigen::Matrix<real_t, 1, 1>;
  static State retract(const State& x, const State& d) { return x + d; }
  static State difference(const State& x, const State& y) { return x - y; }
  State step(int, const State& x, const Control& u) const { return x + u; }
  template <class A, class B>
  void step_jacobian(int, const State&, const Control&, A& a, B& b) const {
    a(0) = real_t(1);
    b(0) = real_t(1);
  }
  template <class R>
  void stage_residual(int, const State& x, const Control& u, R& r) const {
    r << x(0), u(0);
  }
  template <class JX, class JU>
  void residual_jacobians(int, const State&, const Control&, JX& jx,
                          JU& ju) const {
    jx << real_t(1), real_t(0);
    ju << real_t(0), real_t(1);
  }
  template <class R>
  void terminal_residual(const State& x, R& r) const {
    r = x;
  }
  template <class J>
  void terminal_jacobian(const State&, J& j) const {
    j(0) = real_t(1);
  }
  void control_bounds(int, Control& lower, Control& upper) const {
    lower(0) = real_t(-1);
    upper(0) = real_t(1);
  }
};

int main() {
  const UKF::Quaternion attitude = UKF::Quaternion::Identity();
  static NMPC::Solver<Level>::Workspace workspace;
  NMPC::Solver<Level> solver(workspace);
  Level::State x;
  x(0) = real_t(0.5) + attitude.w() - real_t(1);
  solver.reset(x, Level::Control::Zero());
  return solver.solve(Level(), x, NMPC::Options<real_t>()).usable() ? 0 : 1;
}
