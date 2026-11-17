/*
Small problems for the unit tests. None of them is a vehicle: each exercises
one part of the problem concept.

  Stop_Problem      a point mass on a line that must come to rest at the origin
                    within an acceleration limit (the README's first example)
  Linear_Problem    linear dynamics and a quadratic cost, so that one iteration
                    can be compared with the solution of the condensed least
                    squares problem
  Attitude_Problem  a rigid body's attitude and body rate driven by torques, on
                    the quaternion manifold
  Tail_Problem      the stop problem with a tail_control
  Failing_Problem   the stop problem with a model that can be told to fail
  Pinned_Problem    a scalar integrator over two stages whose second command is
                    pinned by equal bounds
  Indefinite_Problem the stop problem with its own cost expansion, which can be
                    told to give the first stage a command Hessian that no
                    regularisation in budget makes positive definite
*/

#ifndef NMPC_TEST_PROBLEMS_H
#define NMPC_TEST_PROBLEMS_H

#include <cmath>
#include <limits>

#include "NMPC/Solver.h"

template <class T>
struct Stop_Problem {
  using Scalar = T;
  enum { NX = 2, NU = 1, NR = 3, NT = 2, Stages = 20 };
  using State = Eigen::Matrix<T, 2, 1>;   /* position (m), velocity (m/s) */
  using Control = Eigen::Matrix<T, 1, 1>; /* acceleration (m/s^2) */

  T dt = T(0.05); /* stage length, s */
  T limit = T(4); /* acceleration limit, m/s^2 */

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

  /* Residuals carry the root of the stage length, so the weights are per
     second. */
  template <class R>
  void stage_residual(int, const State& x, const Control& u, R& r) const {
    const T s = std::sqrt(dt);
    r << s * x(0), s * T(0.5) * x(1), s * T(0.1) * u(0);
  }
  template <class JX, class JU>
  void residual_jacobians(int, const State&, const Control&, JX& jx,
                          JU& ju) const {
    const T s = std::sqrt(dt);
    jx << s, T(0), T(0), s * T(0.5), T(0), T(0);
    ju << T(0), T(0), s * T(0.1);
  }
  template <class R>
  void terminal_residual(const State& x, R& r) const {
    r = T(3) * x;
  }
  template <class JX>
  void terminal_jacobian(const State&, JX& jx) const {
    jx = T(3) * JX::Identity();
  }

  void control_bounds(int, Control& lower, Control& upper) const {
    lower(0) = -limit;
    upper(0) = limit;
  }
};

template <class T>
struct Tail_Problem : Stop_Problem<T> {
  using typename Stop_Problem<T>::Control;
  T tail = T(0.25);
  void tail_control(int, Control& u) const { u(0) = tail; }
};

template <class T>
struct Failing_Problem : Stop_Problem<T> {
  using typename Stop_Problem<T>::State;
  using typename Stop_Problem<T>::Control;
  bool fail = false;
  State step(int k, const State& x, const Control& u) const {
    State next = Stop_Problem<T>::step(k, x, u);
    if (fail) {
      next(1) = std::numeric_limits<T>::quiet_NaN();
    }
    return next;
  }
};

/*
The stop problem's own Gauss-Newton expansion. With `indefinite`, every stage
after the first weighs its command four times as much (so that its gains
change), and the first stage's command Hessian is -1000, so that a backward pass
completes every stage but the first and then fails.
*/
template <class T>
struct Indefinite_Problem : Stop_Problem<T> {
  using typename Stop_Problem<T>::State;
  using typename Stop_Problem<T>::Control;
  bool indefinite = false;

  template <class R>
  void stage_residual(int k, const State& x, const Control& u, R& r) const {
    Stop_Problem<T>::stage_residual(k, x, u, r);
    if (indefinite && k > 0) {
      r(2) *= T(2);
    }
  }
  template <class TX, class TU, class XX, class UX, class UU>
  void cost_expansion(int k, const State& x, const Control& u, TX& lx, TU& lu,
                      XX& Lxx, UX& Lux, UU& Luu) const {
    const T h = this->dt, w = indefinite && k > 0 ? T(0.04) : T(0.01);
    lx << h * x(0), h * T(0.25) * x(1);
    lu << h * w * u(0);
    Lxx << h, T(0), T(0), h * T(0.25);
    Lux.setZero();
    Luu << (indefinite && k == 0 ? T(-1000) : h * w);
  }
};

/*
x[k+1] = x[k] + u[k], with cost (x^2 + u^2) / 2 at each stage and x^2 / 2 at the
end, over two stages. The second stage's command has equal bounds, so it cannot
move: its gain must be zero, and the first stage's gain is then -2/3 (the second
stage passes its state on unchanged, so the cost to go there is x^2).
*/
template <class T>
struct Pinned_Problem {
  using Scalar = T;
  enum { NX = 1, NU = 1, NR = 2, NT = 1, Stages = 2 };
  using State = Eigen::Matrix<T, 1, 1>;
  using Control = Eigen::Matrix<T, 1, 1>;

  T pinned = T(0);

  static State retract(const State& x, const State& d) { return x + d; }
  static State difference(const State& x, const State& y) { return x - y; }
  State step(int, const State& x, const Control& u) const { return x + u; }
  template <class A, class B>
  void step_jacobian(int, const State&, const Control&, A& a, B& b) const {
    a(0, 0) = T(1);
    b(0, 0) = T(1);
  }
  template <class R>
  void stage_residual(int, const State& x, const Control& u, R& r) const {
    r << x(0), u(0);
  }
  template <class JX, class JU>
  void residual_jacobians(int, const State&, const Control&, JX& jx,
                          JU& ju) const {
    jx << T(1), T(0);
    ju << T(0), T(1);
  }
  template <class R>
  void terminal_residual(const State& x, R& r) const {
    r = x;
  }
  template <class JX>
  void terminal_jacobian(const State&, JX& jx) const {
    jx(0, 0) = T(1);
  }
  void control_bounds(int k, Control& lower, Control& upper) const {
    if (k == 1) {
      lower(0) = pinned;
      upper(0) = pinned;
    } else {
      lower(0) = T(-10);
      upper(0) = T(10);
    }
  }
};

/*
x[k+1] = A x[k] + B u[k] with three states and two commands, a quadratic cost on
states and commands with cross terms, and bounds that can be tightened.
*/
template <class T>
struct Linear_Problem {
  using Scalar = T;
  enum { NX = 3, NU = 2, NR = 5, NT = 3, Stages = 8 };
  using State = Eigen::Matrix<T, 3, 1>;
  using Control = Eigen::Matrix<T, 2, 1>;

  Eigen::Matrix<T, 3, 3> A;
  Eigen::Matrix<T, 3, 2> B;
  Eigen::Matrix<T, 5, 3> Wx;
  Eigen::Matrix<T, 5, 2> Wu;
  Eigen::Matrix<T, 3, 3> Wt;
  State x_target;
  Control lower_bound, upper_bound;
  State state_scale;
  Control control_scale;

  Linear_Problem() {
    A << T(1), T(0.1), T(0), T(0), T(0.95), T(0.1), T(0.02), T(0), T(0.9);
    B << T(0), T(0.01), T(0.1), T(0), T(0.05), T(0.08);
    Wx << T(2), T(0), T(0.3), T(0), T(1), T(0), T(0), T(0.2), T(0.5), T(0),
        T(0), T(0), T(0), T(0), T(0);
    Wu << T(0), T(0), T(0), T(0), T(0), T(0), T(0.3), T(0.05), T(0), T(0.2);
    Wt << T(4), T(0), T(0), T(0), T(3), T(0), T(0.5), T(0), T(2);
    x_target << T(1), T(0), T(-0.5);
    lower_bound = Control::Constant(-std::numeric_limits<T>::infinity());
    upper_bound = Control::Constant(std::numeric_limits<T>::infinity());
    state_scale << T(1), T(0.5), T(2);
    control_scale << T(1), T(3);
  }

  static State retract(const State& x, const State& d) { return x + d; }
  static State difference(const State& x, const State& y) { return x - y; }

  State step(int, const State& x, const Control& u) const {
    return A * x + B * u;
  }
  template <class MA, class MB>
  void step_jacobian(int, const State&, const Control&, MA& a, MB& b) const {
    a = A;
    b = B;
  }
  template <class R>
  void stage_residual(int, const State& x, const Control& u, R& r) const {
    r = Wx * (x - x_target) + Wu * u;
  }
  template <class JX, class JU>
  void residual_jacobians(int, const State&, const Control&, JX& jx,
                          JU& ju) const {
    jx = Wx;
    ju = Wu;
  }
  template <class R>
  void terminal_residual(const State& x, R& r) const {
    r = Wt * (x - x_target);
  }
  template <class JX>
  void terminal_jacobian(const State&, JX& jx) const {
    jx = Wt;
  }
  void control_bounds(int, Control& lower, Control& upper) const {
    lower = lower_bound;
    upper = upper_bound;
  }
};

/* A rigid body's attitude (Hamilton quaternion, body to world) and body rate,
   driven by torques. */
template <class T>
struct Attitude_State {
  Eigen::Quaternion<T> q;
  Eigen::Matrix<T, 3, 1> w;
};

template <class T>
struct Attitude_Problem {
  using Scalar = T;
  enum { NX = 6, NU = 3, NR = 9, NT = 6, Stages = 25 };
  using State = Attitude_State<T>;
  using Control = Eigen::Matrix<T, 3, 1>;
  using Tangent = Eigen::Matrix<T, 6, 1>;
  using V3 = Eigen::Matrix<T, 3, 1>;
  using M3 = Eigen::Matrix<T, 3, 3>;

  T dt = T(0.04);
  V3 inertia = V3(T(0.02), T(0.025), T(0.04)); /* kg m^2 */
  T torque_limit = T(0.05);                    /* N m */
  Eigen::Quaternion<T> target = Eigen::Quaternion<T>::Identity();
  V3 attitude_weight = V3(T(4), T(4), T(2));
  T rate_weight = T(0.5), torque_weight = T(2), terminal_weight = T(4);

  static State retract(const State& x, const Tangent& d) {
    State o;
    o.q = NMPC::attitude_retract<T>(x.q, V3(d.template head<3>()));
    o.w = x.w + d.template tail<3>();
    return o;
  }
  static Tangent difference(const State& x, const State& y) {
    Tangent d;
    d.template head<3>() = NMPC::attitude_difference<T>(x.q, y.q);
    d.template tail<3>() = x.w - y.w;
    return d;
  }
  static bool finite(const State& x) {
    return x.q.coeffs().allFinite() && x.w.allFinite();
  }

  /* The right Jacobian of SO(3): exp(phi + e) = exp(phi) exp(Jr(phi) e) to
     first order. */
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
    const V3 momentum(inertia.cwiseProduct(x.w));
    return x.w + dt * (u - x.w.cross(momentum)).cwiseQuotient(inertia);
  }

  State step(int, const State& x, const Control& u) const {
    State o;
    o.w = rate_next(x, u);
    o.q = (x.q * NMPC::so3_exp<T>(V3(dt * o.w))).normalized();
    return o;
  }

  template <class MA, class MB>
  void step_jacobian(int, const State& x, const Control& u, MA& a,
                     MB& b) const {
    const V3 w_next = rate_next(x, u);
    const M3 J = inertia.asDiagonal();
    const M3 Jinv = inertia.cwiseInverse().asDiagonal();
    const M3 dw_dw =
        M3::Identity() -
        dt * Jinv * (NMPC::hat<T>(x.w) * J - NMPC::hat<T>(V3(J * x.w)));
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
    r.template segment<3>(0) =
        s *
        attitude_weight.cwiseProduct(NMPC::attitude_difference<T>(x.q, target));
    r.template segment<3>(3) = s * rate_weight * x.w;
    r.template segment<3>(6) = s * torque_weight * (u / torque_limit);
  }
  template <class JX, class JU>
  void residual_jacobians(int, const State& x, const Control&, JX& jx,
                          JU& ju) const {
    const T s = std::sqrt(dt);
    const V3 phi = NMPC::attitude_difference<T>(x.q, target);
    jx.setZero();
    ju.setZero();
    jx.template block<3, 3>(0, 0) = s * attitude_weight.asDiagonal() *
                                    NMPC::so3_left_jacobian_inverse<T>(phi);
    jx.template block<3, 3>(3, 3) = s * rate_weight * M3::Identity();
    ju.template block<3, 3>(6, 0) =
        (s * torque_weight / torque_limit) * M3::Identity();
  }
  template <class R>
  void terminal_residual(const State& x, R& r) const {
    r.template head<3>() =
        terminal_weight * NMPC::attitude_difference<T>(x.q, target);
    r.template tail<3>() = terminal_weight * x.w;
  }
  template <class JX>
  void terminal_jacobian(const State& x, JX& jx) const {
    jx.setZero();
    jx.template block<3, 3>(0, 0) =
        terminal_weight * NMPC::so3_left_jacobian_inverse<T>(
                              NMPC::attitude_difference<T>(x.q, target));
    jx.template block<3, 3>(3, 3) = terminal_weight * M3::Identity();
  }
  void control_bounds(int, Control& lower, Control& upper) const {
    lower = Control::Constant(-torque_limit);
    upper = Control::Constant(torque_limit);
  }
};

#endif
