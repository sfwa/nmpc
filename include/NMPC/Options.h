/*
Part of NMPC. Distributed under the MIT licence: see LICENSE in the root of
the repository.
*/

#ifndef NMPC_OPTIONS_H
#define NMPC_OPTIONS_H

#include "NMPC/Config.h"
#include <limits>

namespace NMPC {

/*
What a solve may spend, in units of work and never in seconds, and the
tolerances it works to. The defaults are the embedded defaults: one iteration
a solve, five line-search candidates and one regularisation retry. A restart
after the plan and the state have parted wants a few iterations, which is a
matter of passing other options to that one solve.

Regularisation and the stage QP work in the problem's normalised units, where
each state and command is divided by its scale (see Problem.h).
*/
template <class T>
struct Options {
    /* Gauss-Newton iterations a solve. One is a real-time iteration from the previous plan. */
    int max_iterations = 1;

    /* Step lengths tried by the line search: 1, 1/2, 1/4 and so on. */
    int max_line_search = 5;

    /* Backward passes repeated with a larger regularisation after one fails or its step is refused. */
    int max_regularisation_retries = 1;

    /* Passes of the active-set method in each stage's box QP. */
    int max_box_iterations = 32;

    /*
    Levenberg-Marquardt regularisation on each stage's command Hessian. It starts at
    regularisation_initial at every solve: a value carried from one solve to the next rose to its
    maximum in a quiet stretch where no step could be verified, after which the controller stopped
    responding (one of the three defects of the earlier solver).
    */
    T regularisation_initial = T(1e-6);
    T regularisation_min = T(1e-9);
    T regularisation_max = T(1e6);
    T regularisation_factor = T(10);

    /* A trial step is accepted when it achieves this share of its predicted decrease. */
    T accept_ratio = T(1e-2);

    /*
    The precision floor: a solve stops, as converged, when the decrease that the backward pass
    predicts is below this share of the cost. Below about this level a step cannot be told from
    rounding, so the line search would refuse it and raise the regularisation for nothing. The
    default is a few units of the last place of the scalar type. Zero turns the test off.
    */
    T precision_floor = T(8) * std::numeric_limits<T>::epsilon();

    /*
    An optional deadline, for a solve that must publish within its period whatever the work it meets. When
    out_of_time is set, it is asked (with time_context) before each backward pass and each line-search trial;
    when it answers true the solve stops there with what it has: the last accepted step's plan, or, if no step
    was accepted, the previous plan rolled out from the new state (NoProgress). Result::out_of_time says it
    happened. Null, the default, leaves the solve to its work budgets alone, so a solve does the same arithmetic
    wherever it runs.
    */
    bool (*out_of_time)(void* context) = nullptr;
    void* time_context = nullptr;
};

}

#endif
