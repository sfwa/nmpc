/*
Part of NMPC. Distributed under the MIT licence: see LICENSE in the root of
the repository.
*/

#ifndef NMPC_BOXQP_H
#define NMPC_BOXQP_H

#include "NMPC/Types.h"
#include <Eigen/Cholesky>
#include <cmath>
#include <limits>

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
either the last pass took the whole Newton step, which ends at that minimiser,
or every free gradient is within `tolerance` of the size of the terms that make
it, |g(i)| + sum_j |H(i, j) d(j)|, which is where rounding leaves it. A held
variable is let go only if its gradient pulls it inward by more than the same
share of its terms, so that rounding cannot let go of a variable and hold it
again for ever. The tolerance is relative; Options::box_tolerance is a few
units of the last place of the type.

Returns the passes used, or -1 if a factorisation failed. If the pass budget
runs out, d is inside the bounds and the factor matches what is held.
*/
template <class T, int M>
int solve_box_qp(const Matrix<T, M, M>& H, const Vector<T, M>& g,
        const Vector<T, M>& lower, const Vector<T, M>& upper,
        Vector<T, M>& d, Eigen::Array<bool, M, 1>& held,
        Eigen::LLT<Matrix<T, M, M>>& factor, int max_passes,
        T tolerance = T(16) * std::numeric_limits<T>::epsilon()) {
    d = d.cwiseMax(lower).cwiseMin(upper);
    Vector<T, M> gradient = g + H * d;
    /* The size of the terms that make each gradient, against which it is compared. */
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
            (d(i) <= lower(i) && gradient(i) > size(i)) || (d(i) >= upper(i) && gradient(i) < -size(i));
    }

    Eigen::Array<bool, M, 1> factored_for = held;
    bool factored = false;
    /* The factorisation with the rows and columns of the held variables made the identity. */
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
            /* At the minimiser over the free variables: let go of the held one pulled hardest inward. */
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
                return pass + 1;
            }
            held(release) = false;
            at_minimiser = false;
            continue;
        }

        /* Newton step on the free variables with the held ones where they are, as far as the first bound. */
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
        T step = T(1);
        int limit = -1;
        for (int i = 0; i < M; i++) {
            if (held(i)) {
                continue;
            }
            if (search(i) > T(0) && d(i) + search(i) > upper(i)) {
                const T reach = (upper(i) - d(i)) / search(i);
                if (reach < step) {
                    step = reach;
                    limit = i;
                }
            } else if (search(i) < T(0) && d(i) + search(i) < lower(i)) {
                const T reach = (lower(i) - d(i)) / search(i);
                if (reach < step) {
                    step = reach;
                    limit = i;
                }
            }
        }
        d = d + step * search;
        if (limit >= 0) {
            d(limit) = search(limit) > T(0) ? upper(limit) : lower(limit);
            held(limit) = true;
        }
        /* The whole Newton step ends at the minimiser over the free variables, to rounding. */
        at_minimiser = limit < 0;
        gradient = g + H * d;
        measure();
    }
    return refactor() ? max_passes : -1;
}

}

#endif
