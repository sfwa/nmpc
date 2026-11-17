/*
Part of NMPC. Distributed under the MIT licence: see LICENSE in the root of
the repository.
*/

#ifndef NMPC_BOXQP_H
#define NMPC_BOXQP_H

#include <Eigen/Cholesky>
#include <cmath>
#include <limits>

#include "NMPC/Types.h"

namespace NMPC {

/*
The stage problem of the iLQR with bounds on its commands:

    minimise  g'd + d'Hd / 2   subject to   lower <= d <= upper

for a symmetric positive definite H, by a primal active-set method.

The start is the given d moved inside the bounds, with every variable held
that sits on a bound and is pushed against it by the gradient. Each pass finds
the minimiser over the variables that are free and moves toward it until the
first of them meets a bound, which is then held as well. At a minimiser over
the free variables, the held variable whose gradient most wants it off its
bound is let go, one at a time, which is what makes the passes finite.

It ends at the problem's solution with `held` marking the variables on a bound
and `factor` holding the Cholesky factor of H with their rows and columns
replaced by the identity, which is what the stage's feedback gain is solved
with. A held variable gets no feedback.

A projected Newton method with a line search, which an earlier solver used,
stops short when a variable starts a hair inside a bound that the minimiser
lies beyond: every trial is cut back to the bound, none gains what the line
search asks, and the variable is left free. Warm-started from the last solve's
step, the same stage then does it again at every solve. The contouring
problem's clock rate, which sits on its bound nearly all the time, met this
and the plan stopped improving. This method does not have that failure.

A variable whose bounds meet (lower == upper) is held from the start and never
let go, whatever its gradient: it cannot move, so it must get no feedback, and
its gain would otherwise reach into the stages before it through the cost to go.

When it is at the minimiser over the free variables is decided without an
absolute threshold, so that it does not depend on the scale of the problem:
either the last pass took the whole Newton step, whose free variables then take
the Newton target's values (the minimiser, to the accuracy of the solve),
or every free gradient is within `tolerance` of the size of the terms that make
it, |g(i)| + sum_j |H(i, j) d(j)|, which is where rounding leaves it. A held
variable is let go only if its gradient pulls it inward by more than the same
share of its terms, so that rounding cannot let go of a variable and hold it
again for ever. The tolerance is relative; Options::box_tolerance is a few
units of the last place of the type.

Returns the passes used, or -1 if a factorisation failed. `finished`, if given,
says whether the method reached the solution within the budget. If the pass
budget runs out, d is inside the bounds and the factor matches what is held,
but d need not be the solution and, from a poor start, can be worse than d = 0.
From d = 0, where 0 is inside the bounds, it never is, since no pass raises the
objective.
*/
template <class T, int M>
int solve_box_qp(const Matrix<T, M, M>& H, const Vector<T, M>& g,
                 const Vector<T, M>& lower, const Vector<T, M>& upper,
                 Vector<T, M>& d, Eigen::Array<bool, M, 1>& held,
                 Eigen::LLT<Matrix<T, M, M>>& factor, int max_passes,
                 T tolerance = T(16) * std::numeric_limits<T>::epsilon(),
                 bool* finished = nullptr) {
  if (finished) {
    *finished = false;
  }
  d = d.cwiseMax(lower).cwiseMin(upper);
  Vector<T, M> gradient = g + H * d;
  /* The size of the terms that make each gradient, against which it is
     compared. */
  Vector<T, M> size;
  auto measure = [&]() {
    for (int i = 0; i < M; i++) {
      T s = std::abs(g(i));
      for (int j = 0; j < M; j++) {
        s += std::abs(H(i, j) * d(j));
      }
      size(i) = tolerance * s;
    }
  };
  measure();
  for (int i = 0; i < M; i++) {
    held(i) = !(upper(i) > lower(i)) ||
              (d(i) <= lower(i) && gradient(i) > size(i)) ||
              (d(i) >= upper(i) && gradient(i) < -size(i));
  }

  Eigen::Array<bool, M, 1> factored_for = held;
  bool factored = false;
  /* The factorisation with the rows and columns of the held variables made the
     identity. */
  auto refactor = [&]() {
    Matrix<T, M, M> masked = H;
    for (int i = 0; i < M; i++) {
      if (held(i)) {
        masked.row(i).setZero();
        masked.col(i).setZero();
        masked(i, i) = T(1);
      }
    }
    factor.compute(masked);
    factored_for = held;
    factored = true;
    return factor.info() == Eigen::Success;
  };

  bool at_minimiser = false;
  for (int pass = 0; pass < max_passes; pass++) {
    if ((!factored || (held != factored_for).any()) && !refactor()) {
      return -1;
    }
    if (!at_minimiser) {
      at_minimiser = true;
      for (int i = 0; i < M; i++) {
        if (!held(i) && std::abs(gradient(i)) > size(i)) {
          at_minimiser = false;
          break;
        }
      }
    }
    if (at_minimiser || held.all()) {
      /* At the minimiser over the free variables: let go of the held one pulled
         hardest inward. */
      int release = -1;
      T strongest = T(0);
      for (int i = 0; i < M; i++) {
        if (!held(i) || !(upper(i) > lower(i))) {
          continue;
        }
        const T inward = d(i) <= lower(i) ? -gradient(i) : gradient(i);
        if (inward > size(i) && inward > strongest) {
          strongest = inward;
          release = i;
        }
      }
      if (release < 0) {
        if (finished) {
          *finished = true;
        }
        return pass + 1;
      }
      held(release) = false;
      at_minimiser = false;
      continue;
    }

    /* Newton step on the free variables with the held ones where they are, as
       far as the first bound. */
    Vector<T, M> fixed = d;
    for (int i = 0; i < M; i++) {
      if (!held(i)) {
        fixed(i) = T(0);
      }
    }
    Vector<T, M> rhs = -(g + H * fixed);
    for (int i = 0; i < M; i++) {
      if (held(i)) {
        rhs(i) = T(0);
      }
    }
    const Vector<T, M> target = factor.solve(rhs);
    Vector<T, M> search = target - d;
    for (int i = 0; i < M; i++) {
      if (held(i)) {
        search(i) = T(0);
      }
    }
    /*
    Which bound is met first is decided on the target itself, not on d + (target
    - d), which can lose it. A target past a bound always limits the step, even
    when the fraction of the step to that bound rounds to one or more (from -1e8
    to a bound of 0 with a target of 1 it is 1 in single precision): the step
    then stops on the bound and holds it rather than taking the target outside
    the bounds.
    */
    T nearest = T(0);
    int limit = -1;
    for (int i = 0; i < M; i++) {
      if (held(i)) {
        continue;
      }
      const bool above = target(i) > upper(i);
      if (above || target(i) < lower(i)) {
        const T reach = ((above ? upper(i) : lower(i)) - d(i)) / search(i);
        if (limit < 0 || reach < nearest) {
          nearest = reach;
          limit = i;
        }
      }
    }
    if (limit < 0) {
      /*
      The whole Newton step: the free variables take the target's values
      themselves, which is the minimiser over them to rounding. Adding the step
      to d instead could lose the target to cancellation when d is far from it
      (from 1e8 to 1 in single precision gives 0).
      */
      for (int i = 0; i < M; i++) {
        if (!held(i)) {
          d(i) = target(i);
        }
      }
      at_minimiser = true;
    } else {
      if (nearest < T(1)) {
        d = (d + nearest * search).cwiseMax(lower).cwiseMin(upper);
      } else {
        /* The bound is reached at the whole step to rounding: take the target,
           as above, within the bounds. */
        for (int i = 0; i < M; i++) {
          if (!held(i)) {
            d(i) = target(i) > upper(i)
                       ? upper(i)
                       : (target(i) < lower(i) ? lower(i) : target(i));
          }
        }
      }
      d(limit) = target(limit) > upper(limit) ? upper(limit) : lower(limit);
      held(limit) = true;
      at_minimiser = false;
    }
    gradient = g + H * d;
    measure();
  }
  return refactor() ? max_passes : -1;
}

}  // namespace NMPC

#endif
