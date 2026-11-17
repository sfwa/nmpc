/*
Part of NMPC. Distributed under the MIT licence: see LICENSE in the root of
the repository.
*/

#ifndef NMPC_RESULT_H
#define NMPC_RESULT_H

#include "NMPC/Config.h"
#include <limits>

namespace NMPC {

/* Why a solve ended. */
enum class Status {
    /* The predicted decrease fell below the precision floor. The plan is the best this precision can verify. */
    Converged,
    /* The iteration budget was spent. With one iteration a solve this is the normal outcome. */
    IterationLimit,
    /* No step was accepted within the line search and regularisation budgets. The plan is the rollout of the previous one. */
    NoProgress,
    /* The rollout from the given state was not finite. The workspace still holds the previous plan. */
    NumericalFailure
};

/*
What a solve did and what the application may act on. A solve that is not
usable left the plan in the workspace exactly as it was, so the plan in force
can fly on into its next stage while the application decides what to do.
*/
template <class T>
struct Result {
    Status status = Status::NumericalFailure;

    /* Costs of the plan rolled out from the new state before the first step, and of the plan returned. */
    T cost_before = std::numeric_limits<T>::quiet_NaN();
    T cost_after = std::numeric_limits<T>::quiet_NaN();

    /*
    The decrease the last backward pass predicted for a full step. Against cost_before this says how
    far the plan was from optimal for the new state: a large share means the plan and the state
    have parted company, which is the sign for a restart.
    */
    T predicted_decrease = T(0);

    /* The regularisation in force at the end of the solve. */
    T regularisation = T(0);

    /* Work done. */
    int iterations = 0;
    int accepted_steps = 0;
    int backward_passes = 0;
    int line_search_trials = 0;
    int box_iterations = 0;

    /* Stages of the returned plan with at least one command on a bound. */
    int stages_on_bound = 0;

    /* Options::out_of_time stopped the solve before its work budgets were spent. */
    bool out_of_time = false;

    bool usable() const { return status != Status::NumericalFailure; }
};

}

#endif
