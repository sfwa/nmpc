/*
Swinging a pendulum up and holding it inverted with a torque too weak to lift
it directly. The motor can hold only about a third of the pendulum's weight
torque, so the controller must pump energy in over several swings, which a
linear controller cannot plan.

The state is the angle from hanging straight down and the angular rate. A
stage is 20 ms, made of two semi-implicit Euler steps so that the step stays
stable across the whole swing. The cost asks for the inverted position
(angle pi), and the horizon of 1.5 s is long enough to see one swing ahead.

It is flown the way a real controller runs: one Gauss-Newton iteration per
stage from the plan of the stage before. Between solves, nothing here uses
the plan's gains; Policy.h has them for a faster inner loop.
*/

#include <cmath>
#include <cstdio>

#include "NMPC/Solver.h"

/* M_PI is not standard C++. */
constexpr double Pendulum_Pi = 3.14159265358979323846;

/*
Everything a user would tune, in one place, with its units.
*/
struct Pendulum_Tuning {
  float length = 0.5f;        /* m */
  float mass = 0.2f;          /* kg, at the end of a massless rod */
  float damping = 0.005f;     /* N m s/rad */
  float torque_limit = 0.35f; /* N m: about a third of m g l */

  float angle_weight = 3.0f;    /* per rad, per second */
  float rate_weight = 0.3f;     /* per rad/s, per second */
  float torque_weight = 0.5f;   /* per unit of the torque limit, per second */
  float terminal_weight = 4.0f; /* per rad and per rad/s */
};

template <class T>
struct Pendulum_Problem {
  using Scalar = T;
  enum { NX = 2, NU = 1, NR = 3, NT = 2, Stages = 75 };
  using State =
      Eigen::Matrix<T, 2, 1>; /* angle from hanging down (rad), rate (rad/s) */
  using Control = Eigen::Matrix<T, 1, 1>; /* torque (N m) */

  static constexpr int Substeps = 2;
  T dt = T(0.02);   /* stage length, s */
  T g = T(9.80665); /* m/s^2 */
  Pendulum_Tuning tuning;

  static State retract(const State& x, const State& d) { return x + d; }
  static State difference(const State& x, const State& y) { return x - y; }

  T inertia() const {
    return T(tuning.mass) * T(tuning.length) * T(tuning.length);
  }

  /* One substep: the rate first, then the angle from the new rate. */
  State substep(const State& x, T u, T h) const {
    const T acceleration =
        (u - T(tuning.damping) * x(1) -
         T(tuning.mass) * g * T(tuning.length) * std::sin(x(0))) /
        inertia();
    const T rate = x(1) + h * acceleration;
    return State(x(0) + h * rate, rate);
  }

  State step(int, const State& x, const Control& u) const {
    const T h = dt / T(Substeps);
    State next = x;
    for (int i = 0; i < Substeps; i++) {
      next = substep(next, u(0), h);
    }
    return next;
  }

  template <class A, class B>
  void step_jacobian(int, const State& x, const Control& u, A& a, B& b) const {
    const T h = dt / T(Substeps);
    const T J = inertia();
    a.setIdentity();
    b.setZero();
    State s = x;
    for (int i = 0; i < Substeps; i++) {
      /* The substep's own Jacobian, chained onto what came before. */
      const T da_dangle =
          -T(tuning.mass) * g * T(tuning.length) * std::cos(s(0)) / J;
      const T da_drate = -T(tuning.damping) / J;
      Eigen::Matrix<T, 2, 2> m;
      m(1, 0) = h * da_dangle;
      m(1, 1) = T(1) + h * da_drate;
      m(0, 0) = T(1) + h * m(1, 0);
      m(0, 1) = h * m(1, 1);
      Eigen::Matrix<T, 2, 1> n;
      n(1) = h / J;
      n(0) = h * n(1);
      a = (m * a).eval();
      b = (m * b + n).eval();
      s = substep(s, u(0), h);
    }
  }

  template <class R>
  void stage_residual(int, const State& x, const Control& u, R& r) const {
    const T s = std::sqrt(dt);
    r << s * T(tuning.angle_weight) * (x(0) - T(Pendulum_Pi)),
        s * T(tuning.rate_weight) * x(1),
        s * T(tuning.torque_weight) * u(0) / T(tuning.torque_limit);
  }
  template <class JX, class JU>
  void residual_jacobians(int, const State&, const Control&, JX& jx,
                          JU& ju) const {
    const T s = std::sqrt(dt);
    jx.setZero();
    ju.setZero();
    jx(0, 0) = s * T(tuning.angle_weight);
    jx(1, 1) = s * T(tuning.rate_weight);
    ju(2, 0) = s * T(tuning.torque_weight) / T(tuning.torque_limit);
  }
  template <class R>
  void terminal_residual(const State& x, R& r) const {
    r << T(tuning.terminal_weight) * (x(0) - T(Pendulum_Pi)),
        T(tuning.terminal_weight) * x(1);
  }
  template <class JX>
  void terminal_jacobian(const State&, JX& jx) const {
    jx = T(tuning.terminal_weight) * JX::Identity();
  }
  void control_bounds(int, Control& lower, Control& upper) const {
    lower(0) = -T(tuning.torque_limit);
    upper(0) = T(tuning.torque_limit);
  }

  State state_scale = State(T(1), T(5));
  Control control_scale = Control::Constant(T(tuning.torque_limit));
};

/*
The controller as an application holds it: the problem, the solver and its
workspace, and what it did last.
*/
class Pendulum_Controller {
 public:
  using Problem = Pendulum_Problem<float>;
  using Solver = NMPC::Solver<Problem>;

  explicit Pendulum_Controller(
      const Pendulum_Tuning& tuning = Pendulum_Tuning())
      : solver(workspace) {
    problem.tuning = tuning;
    problem.control_scale(0) = tuning.torque_limit;
  }

  /* Start from a measured state with no torque planned. */
  void initialise(const Problem::State& x) {
    solver.reset(x, Problem::Control::Zero());
    started = false;
  }

  /*
  One stage: shift the plan by the stage just flown, solve once from the state,
  and return the torque to apply. Returns false if the solve was not usable, in
  which case the torque is the plan in force carried on.
  */
  bool update(const Problem::State& x, float& torque) {
    if (started) {
      solver.advance(problem, 1);
    }
    started = true;
    last = solver.solve(problem, x, options);
    torque = solver.policy().u[0](0);
    if (!last.usable()) {
      solver.reset(x, solver.policy().u[0]);
      return false;
    }
    return true;
  }

  const NMPC::Result<float>& result() const { return last; }
  Problem problem;

 private:
  Solver::Workspace workspace;
  Solver solver;
  NMPC::Options<float> options;
  NMPC::Result<float> last;
  bool started = false;
};

int main() {
  static Pendulum_Controller controller;
  Pendulum_Controller::Problem::State x(0.05f, 0.0f);
  controller.initialise(x);

  int failures = 0;
  for (int stage = 0; stage < 500; stage++) {
    float torque = 0.0f;
    if (!controller.update(x, torque)) {
      failures++;
    }
    if (stage % 25 == 0) {
      std::printf(
          "t %5.2f s  angle %7.3f rad  rate %7.3f rad/s  torque %6.3f N m  "
          "cost %8.3f\n",
          double(stage) * 0.02, double(x(0)), double(x(1)), double(torque),
          double(controller.result().cost_after));
    }
    Pendulum_Controller::Problem::Control u;
    u(0) = torque;
    x = controller.problem.step(0, x, u);
  }
  const float error = std::abs(x(0) - float(Pendulum_Pi));
  std::printf(
      "after 10 s: %.4f rad from inverted, rate %.4f rad/s, %d unusable "
      "solves\n",
      double(error), double(x(1)), failures);
  return error < 0.01f && std::abs(x(1)) < 0.05f && failures == 0 ? 0 : 1;
}
