/*
Part of NMPC. Distributed under the MIT licence: see LICENSE in the root of
the repository.
*/

#ifndef NMPC_CHECK_H
#define NMPC_CHECK_H

#include <algorithm>
#include <cmath>
#include <limits>

#include "NMPC/Problem.h"

/*
Derivative checks for tests. Each compares a problem's analytic derivatives
with central differences in the tangent, at one stage, state and command, and
reports the worst entry. Central differences are a test tool here and never
part of a solve.

The error of an entry is |analytic - numerical| / max(1, |numerical|), so it
is absolute for small entries and relative for large ones. A step in each
direction is `relative_step` times that state's or command's scale. The
default step is the cube root of the type's epsilon, which balances the
truncation and rounding errors of a central difference; tolerances must suit
the precision, with a float check looser than a double one.
*/

namespace NMPC {

template <class T>
struct CheckResult {
  T error = T(0); /* the worst entry's error */
  int row = -1;   /* where it is */
  int column = -1;
  bool passed = true;

  void consider(T analytic, T numerical, int i, int j, T tolerance) {
    const T e =
        std::abs(analytic - numerical) / std::max(T(1), std::abs(numerical));
    if (!(e <= error)) {
      error = e;
      row = i;
      column = j;
    }
    passed = passed && e <= tolerance;
  }

  void merge(const CheckResult& other) {
    if (other.error > error) {
      error = other.error;
      row = other.row;
      column = other.column;
    }
    passed = passed && other.passed;
  }
};

template <class T>
T default_relative_step() {
  return std::cbrt(std::numeric_limits<T>::epsilon());
}

/* The step Jacobians of stage k against central differences of step(). */
template <class Problem>
CheckResult<typename ProblemTraits<Problem>::T> check_step_jacobian(
    const Problem& problem, int k,
    const typename ProblemTraits<Problem>::State& x,
    const typename ProblemTraits<Problem>::Control& u,
    typename ProblemTraits<Problem>::T tolerance,
    typename ProblemTraits<Problem>::T relative_step =
        default_relative_step<typename ProblemTraits<Problem>::T>()) {
  using Tr = ProblemTraits<Problem>;
  using T = typename Tr::T;
  using State = typename Tr::State;
  const typename Tr::Tangent Dx = detail::state_scale(problem);
  const typename Tr::Control Du = detail::control_scale(problem);

  typename Tr::MatA A;
  typename Tr::MatB B;
  problem.step_jacobian(k, x, u, A, B);
  const State next = problem.step(k, x, u);

  CheckResult<T> result;
  for (int j = 0; j < Tr::NX; j++) {
    const T h = relative_step * Dx(j);
    typename Tr::Tangent e = Tr::Tangent::Zero();
    e(j) = h;
    const State plus = problem.step(k, Problem::retract(x, e), u);
    const State minus =
        problem.step(k, Problem::retract(x, typename Tr::Tangent(-e)), u);
    const typename Tr::Tangent column =
        (Problem::difference(plus, next) - Problem::difference(minus, next)) /
        (T(2) * h);
    for (int i = 0; i < Tr::NX; i++) {
      result.consider(A(i, j), column(i), i, j, tolerance);
    }
  }
  for (int j = 0; j < Tr::NU; j++) {
    const T h = relative_step * Du(j);
    typename Tr::Control e = Tr::Control::Zero();
    e(j) = h;
    const State plus = problem.step(k, x, typename Tr::Control(u + e));
    const State minus = problem.step(k, x, typename Tr::Control(u - e));
    const typename Tr::Tangent column =
        (Problem::difference(plus, next) - Problem::difference(minus, next)) /
        (T(2) * h);
    for (int i = 0; i < Tr::NX; i++) {
      result.consider(B(i, j), column(i), i, Tr::NX + j, tolerance);
    }
  }
  return result;
}

/* The stage residual Jacobians of stage k against central differences of
   stage_residual(). Columns past NX are commands. */
template <class Problem>
CheckResult<typename ProblemTraits<Problem>::T> check_residual_jacobians(
    const Problem& problem, int k,
    const typename ProblemTraits<Problem>::State& x,
    const typename ProblemTraits<Problem>::Control& u,
    typename ProblemTraits<Problem>::T tolerance,
    typename ProblemTraits<Problem>::T relative_step =
        default_relative_step<typename ProblemTraits<Problem>::T>()) {
  using Tr = ProblemTraits<Problem>;
  using T = typename Tr::T;
  const typename Tr::Tangent Dx = detail::state_scale(problem);
  const typename Tr::Control Du = detail::control_scale(problem);

  typename Tr::StageJx Jx;
  typename Tr::StageJu Ju;
  problem.residual_jacobians(k, x, u, Jx, Ju);

  CheckResult<T> result;
  typename Tr::StageResidual plus, minus;
  for (int j = 0; j < Tr::NX; j++) {
    const T h = relative_step * Dx(j);
    typename Tr::Tangent e = Tr::Tangent::Zero();
    e(j) = h;
    problem.stage_residual(k, Problem::retract(x, e), u, plus);
    problem.stage_residual(k, Problem::retract(x, typename Tr::Tangent(-e)), u,
                           minus);
    for (int i = 0; i < Tr::NR; i++) {
      result.consider(Jx(i, j), (plus(i) - minus(i)) / (T(2) * h), i, j,
                      tolerance);
    }
  }
  for (int j = 0; j < Tr::NU; j++) {
    const T h = relative_step * Du(j);
    typename Tr::Control e = Tr::Control::Zero();
    e(j) = h;
    problem.stage_residual(k, x, typename Tr::Control(u + e), plus);
    problem.stage_residual(k, x, typename Tr::Control(u - e), minus);
    for (int i = 0; i < Tr::NR; i++) {
      result.consider(Ju(i, j), (plus(i) - minus(i)) / (T(2) * h), i,
                      Tr::NX + j, tolerance);
    }
  }
  return result;
}

/* The terminal residual's Jacobian against central differences of
   terminal_residual(). */
template <class Problem>
CheckResult<typename ProblemTraits<Problem>::T> check_terminal_jacobian(
    const Problem& problem, const typename ProblemTraits<Problem>::State& x,
    typename ProblemTraits<Problem>::T tolerance,
    typename ProblemTraits<Problem>::T relative_step =
        default_relative_step<typename ProblemTraits<Problem>::T>()) {
  using Tr = ProblemTraits<Problem>;
  using T = typename Tr::T;
  const typename Tr::Tangent Dx = detail::state_scale(problem);

  typename Tr::TerminalJx J;
  problem.terminal_jacobian(x, J);

  CheckResult<T> result;
  typename Tr::TerminalResidual plus, minus;
  for (int j = 0; j < Tr::NX; j++) {
    const T h = relative_step * Dx(j);
    typename Tr::Tangent e = Tr::Tangent::Zero();
    e(j) = h;
    problem.terminal_residual(Problem::retract(x, e), plus);
    problem.terminal_residual(Problem::retract(x, typename Tr::Tangent(-e)),
                              minus);
    for (int i = 0; i < Tr::NT; i++) {
      result.consider(J(i, j), (plus(i) - minus(i)) / (T(2) * h), i, j,
                      tolerance);
    }
  }
  return result;
}

/* The gradient of the linear cost term against central differences of
   linear_cost(). Row 0; columns past NX are commands. */
template <class Problem>
CheckResult<typename ProblemTraits<Problem>::T> check_linear_cost_gradient(
    const Problem& problem, int k,
    const typename ProblemTraits<Problem>::State& x,
    const typename ProblemTraits<Problem>::Control& u,
    typename ProblemTraits<Problem>::T tolerance,
    typename ProblemTraits<Problem>::T relative_step =
        default_relative_step<typename ProblemTraits<Problem>::T>()) {
  using Tr = ProblemTraits<Problem>;
  using T = typename Tr::T;
  static_assert(detail::HasLinearCost<Problem>::value,
                "The problem has no linear cost term");
  const typename Tr::Tangent Dx = detail::state_scale(problem);
  const typename Tr::Control Du = detail::control_scale(problem);

  typename Tr::Tangent gx;
  typename Tr::Control gu;
  problem.linear_cost_gradient(k, x, u, gx, gu);

  CheckResult<T> result;
  for (int j = 0; j < Tr::NX; j++) {
    const T h = relative_step * Dx(j);
    typename Tr::Tangent e = Tr::Tangent::Zero();
    e(j) = h;
    const T plus = problem.linear_cost(k, Problem::retract(x, e), u);
    const T minus = problem.linear_cost(
        k, Problem::retract(x, typename Tr::Tangent(-e)), u);
    result.consider(gx(j), (plus - minus) / (T(2) * h), 0, j, tolerance);
  }
  for (int j = 0; j < Tr::NU; j++) {
    const T h = relative_step * Du(j);
    typename Tr::Control e = Tr::Control::Zero();
    e(j) = h;
    const T plus = problem.linear_cost(k, x, typename Tr::Control(u + e));
    const T minus = problem.linear_cost(k, x, typename Tr::Control(u - e));
    result.consider(gu(j), (plus - minus) / (T(2) * h), 0, Tr::NX + j,
                    tolerance);
  }
  return result;
}

/*
A problem's own cost_expansion against the Gauss-Newton expansion built from its
residual Jacobians and linear cost gradient, for a problem that provides both.
Rows and columns are the entries of [lx; lu] and of the upper triangle of Lxx,
then Lux and Luu, in that order.
*/
template <class Problem>
CheckResult<typename ProblemTraits<Problem>::T> check_cost_expansion(
    const Problem& problem, int k,
    const typename ProblemTraits<Problem>::State& x,
    const typename ProblemTraits<Problem>::Control& u,
    typename ProblemTraits<Problem>::T tolerance) {
  using Tr = ProblemTraits<Problem>;
  using T = typename Tr::T;
  static_assert(detail::HasCostExpansion<Problem>::value,
                "The problem has no cost_expansion of its own");

  typename Tr::Tangent lx;
  typename Tr::Control lu;
  typename Tr::MatA Lxx;
  typename Tr::Gain Lux;
  typename Tr::MatU Luu;
  problem.cost_expansion(k, x, u, lx, lu, Lxx, Lux, Luu);

  typename Tr::StageResidual r;
  typename Tr::StageJx Jx;
  typename Tr::StageJu Ju;
  problem.stage_residual(k, x, u, r);
  problem.residual_jacobians(k, x, u, Jx, Ju);
  typename Tr::Tangent gx(Jx.transpose() * r);
  typename Tr::Control gu(Ju.transpose() * r);
  if constexpr (detail::HasLinearCost<Problem>::value) {
    typename Tr::Tangent lgx;
    typename Tr::Control lgu;
    problem.linear_cost_gradient(k, x, u, lgx, lgu);
    gx += lgx;
    gu += lgu;
  }
  const typename Tr::MatA Hxx(Jx.transpose() * Jx);
  const typename Tr::Gain Hux(Ju.transpose() * Jx);
  const typename Tr::MatU Huu(Ju.transpose() * Ju);

  CheckResult<T> result;
  for (int i = 0; i < Tr::NX; i++) {
    result.consider(lx(i), gx(i), i, -1, tolerance);
  }
  for (int i = 0; i < Tr::NU; i++) {
    result.consider(lu(i), gu(i), Tr::NX + i, -1, tolerance);
  }
  for (int j = 0; j < Tr::NX; j++) {
    for (int i = 0; i <= j; i++) {
      result.consider(Lxx(i, j), Hxx(i, j), i, j, tolerance);
    }
    for (int i = 0; i < Tr::NU; i++) {
      result.consider(Lux(i, j), Hux(i, j), Tr::NX + i, j, tolerance);
    }
  }
  for (int j = 0; j < Tr::NU; j++) {
    for (int i = 0; i < Tr::NU; i++) {
      result.consider(Luu(i, j), Huu(i, j), Tr::NX + i, Tr::NX + j, tolerance);
    }
  }
  return result;
}

}  // namespace NMPC

#endif
