/*
A synthetic problem with the dimensions of the multirotor reference
controller (17 tangent states, 5 commands, 22 stage and 12 terminal
residuals, 30 stages), so that NMPC's own benchmark can be timed and gated
without a vehicle model.

It is not the multirotor. The backward pass, the stage QP, the line search
and the bookkeeping depend only on the dimensions and on how many commands
sit on their bounds, so they do the same work as in the reference
controller. The model step and the residuals are cheaper than the
multirotor's, so a whole solve here is faster than the vehicle's; the
vehicle's own benchmark lives in sfwa/models.

The model is a stable linear system with a smooth nonlinearity in its first
twelve states, which makes its Jacobian change with the state, and a pure
integrator in the seventeenth (as the contouring clock is). The reference
asks for more than the bounds allow for part of the horizon, so that some
commands are held at a bound as in a hard manoeuvre.
*/

#ifndef NMPC_BENCHMARK_SYNTHETIC_H
#define NMPC_BENCHMARK_SYNTHETIC_H

#include <cmath>

#include "NMPC/Solver.h"

template <class T, int N = 30>
struct Synthetic_Problem {
  using Scalar = T;
  enum { NX = 17, NU = 5, NR = 22, NT = 12, Stages = N };
  using State = Eigen::Matrix<T, NX, 1>;
  using Control = Eigen::Matrix<T, NU, 1>;

  T dt = T(0.02);
  Eigen::Matrix<T, NX, NX> F; /* the linear part of a stage */
  Eigen::Matrix<T, NX, NU> G;
  Eigen::Matrix<T, 12, 1>
      bend; /* strength of the nonlinearity in each of the first 12 states */
  Eigen::Matrix<T, 16, 1> weight;
  Control control_weight;
  Control lower, upper;
  State state_scale;
  Control control_scale;
  T phase = T(0); /* moves the reference between solves */

  Synthetic_Problem() {
    /* A deterministic, mildly coupled, stable stage map. */
    unsigned seed = 2026u;
    auto next = [&seed]() {
      seed = seed * 1664525u + 1013904223u;
      return T(double((seed >> 9) & 0x3FFFu) / 16383.0 - 0.5);
    };
    F.setZero();
    for (int i = 0; i < NX; i++) {
      F(i, i) = T(1) - T(0.01) * T(i % 5);
      for (int j = 0; j < NX; j++) {
        if (j != i && ((i + 2 * j) % 3) == 0) {
          F(i, j) = T(0.02) * next();
        }
      }
    }
    F.row(16).setZero();
    F(16, 16) = T(1);
    for (int i = 0; i < NX; i++) {
      for (int j = 0; j < NU; j++) {
        G(i, j) = dt * next();
      }
    }
    G.row(16).setZero();
    G(16, 4) = dt;
    for (int i = 0; i < 12; i++) {
      bend(i) = T(0.05) * next();
    }
    weight << T(10), T(10), T(10), T(2), T(2), T(2), T(10), T(10), T(5), T(1),
        T(1), T(1), T(0.2), T(0.2), T(0.2), T(0.2);
    control_weight = Control::Constant(T(0.5));
    lower << T(0), T(0), T(0), T(0), T(0);
    upper << T(8), T(8), T(8), T(8), T(1);
    state_scale = State::Constant(T(0.5));
    state_scale(16) = T(0.02);
    control_scale = Control::Ones();
    control_scale(4) = T(0.2);
  }

  static State retract(const State& x, const State& d) { return x + d; }
  static State difference(const State& x, const State& y) { return x - y; }

  State step(int, const State& x, const Control& u) const {
    State next(F * x + G * u);
    for (int i = 0; i < 12; i++) {
      next(i) += dt * bend(i) * std::sin(x(i));
    }
    return next;
  }
  template <class A, class B>
  void step_jacobian(int, const State& x, const Control&, A& a, B& b) const {
    a = F;
    for (int i = 0; i < 12; i++) {
      a(i, i) += dt * bend(i) * std::cos(x(i));
    }
    b = G;
  }

  /* The reference: a slow wave in the first sixteen states, and commands that
     ask past the upper bounds mid-horizon. */
  T reference_state(int k, int i) const {
    return T(2) * std::sin(phase + T(0.1) * T(k) + T(0.3) * T(i));
  }
  T reference_command(int k, int j) const {
    return j == 4 ? T(1)
                  : T(4) + T(6) * std::sin(phase + T(0.15) * T(k) + T(j));
  }

  template <class R>
  void stage_residual(int k, const State& x, const Control& u, R& r) const {
    const T s = std::sqrt(dt);
    for (int i = 0; i < 16; i++) {
      r(i) = s * weight(i) * (x(i) - reference_state(k, i));
    }
    for (int j = 0; j < NU; j++) {
      r(16 + j) = s * control_weight(j) * (u(j) - reference_command(k, j));
    }
    r(21) = s * T(20) * (u(4) - T(1)) + s * T(0.1) * x(16);
  }
  template <class JX, class JU>
  void residual_jacobians(int, const State&, const Control&, JX& jx,
                          JU& ju) const {
    const T s = std::sqrt(dt);
    jx.setZero();
    ju.setZero();
    for (int i = 0; i < 16; i++) {
      jx(i, i) = s * weight(i);
    }
    for (int j = 0; j < NU; j++) {
      ju(16 + j, j) = s * control_weight(j);
    }
    ju(21, 4) = s * T(20);
    jx(21, 16) = s * T(0.1);
  }
  template <class R>
  void terminal_residual(const State& x, R& r) const {
    for (int i = 0; i < NT; i++) {
      r(i) = weight(i) * (x(i) - reference_state(N, i));
    }
  }
  template <class JX>
  void terminal_jacobian(const State&, JX& jx) const {
    jx.setZero();
    for (int i = 0; i < NT; i++) {
      jx(i, i) = weight(i);
    }
  }
  void control_bounds(int, Control& lo, Control& hi) const {
    lo = lower;
    hi = upper;
  }

  /* A first state, away from the reference. */
  State start() const {
    State x;
    for (int i = 0; i < NX; i++) {
      x(i) = T(0.5) * std::cos(T(0.7) * T(i));
    }
    x(16) = T(0);
    return x;
  }
  Control hover() const {
    Control u;
    u << T(4), T(4), T(4), T(4), T(1);
    return u;
  }
};

#endif
