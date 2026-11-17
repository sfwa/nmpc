/*
Slewing a rigid body to a target attitude with limited torque: the problem
with a state on a manifold.

The attitude is a Hamilton quaternion from body to world, and the solver
never sees it: the problem gives a 6-number tangent (a small rotation in the
world frame, then the body rate) through retract and difference, and NMPC's
Types.h supplies the exponential and logarithm of SO(3). The step integrates
the rate with Euler's equations and the attitude with the exponential map,
and its Jacobian is analytic.

The program slews through 150 degrees about a tilted axis, with the torque on
its bound for most of the manoeuvre, and reports the attitude error.
*/

#include <cmath>
#include <cstdio>

#include "NMPC/Solver.h"

struct Attitude_Tuning {
  float inertia_x = 0.02f, inertia_y = 0.025f, inertia_z = 0.04f; /* kg m^2 */
  float torque_limit = 0.05f;   /* N m on each axis */
  float attitude_weight = 4.0f; /* per rad, per second */
  float rate_weight = 0.5f;     /* per rad/s, per second */
  float torque_weight = 2.0f;   /* per unit of the torque limit, per second */
  float terminal_weight = 4.0f; /* per rad and per rad/s */
};

template <class T>
struct Attitude_State {
  Eigen::Quaternion<T> q;   /* body to world */
  Eigen::Matrix<T, 3, 1> w; /* body rate, rad/s */
};

template <class T>
struct Attitude_Problem {
  using Scalar = T;
  enum { NX = 6, NU = 3, NR = 9, NT = 6, Stages = 40 };
  using State = Attitude_State<T>;
  using Control = Eigen::Matrix<T, 3, 1>; /* body torque, N m */
  using Tangent = Eigen::Matrix<T, 6, 1>;
  using V3 = Eigen::Matrix<T, 3, 1>;
  using M3 = Eigen::Matrix<T, 3, 3>;

  T dt = T(0.05);
  Attitude_Tuning tuning;
  Eigen::Quaternion<T> target = Eigen::Quaternion<T>::Identity();

  V3 inertia() const {
    return V3(T(tuning.inertia_x), T(tuning.inertia_y), T(tuning.inertia_z));
  }

  static State retract(const State& x, const Tangent& d) {
    State o;
    o.q = NMPC::attitude_retract<T>(x.q, V3(d.template head<3>()));
    o.w = x.w + d.template tail<3>();
    return o;
  }
  static Tangent difference(const State& x, const State& y) {
    Tangent d;
    d << NMPC::attitude_difference<T>(x.q, y.q), x.w - y.w;
    return d;
  }
  static bool finite(const State& x) {
    return x.q.coeffs().allFinite() && x.w.allFinite();
  }

  /* exp(phi + e) = exp(phi) exp(Jr(phi) e) to first order. */
  static M3 right_jacobian(const V3& phi) {
    const T angle2 = phi.squaredNorm();
    const M3 h = NMPC::hat<T>(phi);
    T a, b;
    if (angle2 < T(1e-4)) {
      a = T(0.5) - angle2 / T(24);
      b = T(1) / T(6) - angle2 / T(120);
    } else {
      const T angle = std::sqrt(angle2);
      a = (T(1) - std::cos(angle)) / angle2;
      b = (angle - std::sin(angle)) / (angle2 * angle);
    }
    return M3::Identity() - a * h + b * h * h;
  }

  V3 rate_next(const State& x, const Control& u) const {
    const V3 J = inertia();
    const V3 momentum(J.cwiseProduct(x.w));
    return x.w + dt * (u - x.w.cross(momentum)).cwiseQuotient(J);
  }

  State step(int, const State& x, const Control& u) const {
    State o;
    o.w = rate_next(x, u);
    o.q = (x.q * NMPC::so3_exp<T>(V3(dt * o.w))).normalized();
    return o;
  }

  template <class A, class B>
  void step_jacobian(int, const State& x, const Control& u, A& a, B& b) const {
    const V3 Jd = inertia();
    const V3 w_next = rate_next(x, u);
    const M3 J = Jd.asDiagonal();
    const M3 Jinv = Jd.cwiseInverse().asDiagonal();
    const M3 dw_dw(M3::Identity() -
                   dt * Jinv *
                       (NMPC::hat<T>(x.w) * J - NMPC::hat<T>(V3(J * x.w))));
    const M3 dw_du(dt * Jinv);
    const Eigen::Quaternion<T> q_next =
        (x.q * NMPC::so3_exp<T>(V3(dt * w_next))).normalized();
    const M3 dq_dw(dt * q_next.toRotationMatrix() *
                   right_jacobian(V3(dt * w_next)));
    a.setZero();
    b.setZero();
    a.template block<3, 3>(0, 0).setIdentity();
    a.template block<3, 3>(0, 3) = dq_dw * dw_dw;
    a.template block<3, 3>(3, 3) = dw_dw;
    b.template block<3, 3>(0, 0) = dq_dw * dw_du;
    b.template block<3, 3>(3, 0) = dw_du;
  }

  template <class R>
  void stage_residual(int, const State& x, const Control& u, R& r) const {
    const T s = std::sqrt(dt);
    r << s * T(tuning.attitude_weight) *
             NMPC::attitude_difference<T>(x.q, target),
        s * T(tuning.rate_weight) * x.w,
        s * T(tuning.torque_weight) / T(tuning.torque_limit) * u;
  }
  template <class JX, class JU>
  void residual_jacobians(int, const State& x, const Control&, JX& jx,
                          JU& ju) const {
    const T s = std::sqrt(dt);
    jx.setZero();
    ju.setZero();
    jx.template block<3, 3>(0, 0) =
        s * T(tuning.attitude_weight) *
        NMPC::so3_left_jacobian_inverse<T>(
            NMPC::attitude_difference<T>(x.q, target));
    jx.template block<3, 3>(3, 3) = s * T(tuning.rate_weight) * M3::Identity();
    ju.template block<3, 3>(6, 0) =
        s * T(tuning.torque_weight) / T(tuning.torque_limit) * M3::Identity();
  }
  template <class R>
  void terminal_residual(const State& x, R& r) const {
    r << T(tuning.terminal_weight) * NMPC::attitude_difference<T>(x.q, target),
        T(tuning.terminal_weight) * x.w;
  }
  template <class JX>
  void terminal_jacobian(const State& x, JX& jx) const {
    jx.setZero();
    jx.template block<3, 3>(0, 0) =
        T(tuning.terminal_weight) *
        NMPC::so3_left_jacobian_inverse<T>(
            NMPC::attitude_difference<T>(x.q, target));
    jx.template block<3, 3>(3, 3) = T(tuning.terminal_weight) * M3::Identity();
  }
  void control_bounds(int, Control& lower, Control& upper) const {
    lower = Control::Constant(-T(tuning.torque_limit));
    upper = Control::Constant(T(tuning.torque_limit));
  }
};

int main() {
  using Problem = Attitude_Problem<float>;
  using Solver = NMPC::Solver<Problem>;
  using V3 = Problem::V3;
  static Solver::Workspace workspace;
  Solver solver(workspace);
  Problem problem;
  NMPC::Options<float> options;

  const float slew = 150.0f * 3.14159265f / 180.0f;
  Problem::State x;
  x.q = NMPC::so3_exp<float>(V3(V3(1.0f, -0.5f, 2.0f).normalized() * slew));
  x.w.setZero();
  solver.reset(x, Problem::Control::Zero());

  for (int stage = 0; stage < 200; stage++) {
    if (stage > 0) {
      solver.advance(problem, 1);
    }
    const NMPC::Result<float> result = solver.solve(problem, x, options);
    if (!result.usable()) {
      std::printf("stage %d: no usable plan\n", stage);
      return 1;
    }
    const Problem::Control u = solver.policy().u[0];
    if (stage % 20 == 0) {
      std::printf(
          "t %5.2f s  attitude error %6.2f deg  rate %6.3f rad/s  torque %6.3f "
          "%6.3f %6.3f N m\n",
          double(stage) * 0.05,
          double(NMPC::attitude_difference<float>(x.q, problem.target).norm()) *
              57.29578,
          double(x.w.norm()), double(u(0)), double(u(1)), double(u(2)));
    }
    x = problem.step(0, x, u);
  }
  const float error =
      NMPC::attitude_difference<float>(x.q, problem.target).norm();
  std::printf("after 10 s: %.4f deg from the target, rate %.5f rad/s\n",
              double(error) * 57.29578, double(x.w.norm()));
  return error < 1e-3f && x.w.norm() < 1e-3f ? 0 : 1;
}
