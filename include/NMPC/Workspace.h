/*
Part of NMPC. Distributed under the MIT licence: see LICENSE in the root of
the repository.
*/

#ifndef NMPC_WORKSPACE_H
#define NMPC_WORKSPACE_H

#include "NMPC/Policy.h"

namespace NMPC {

/*
Everything a solver keeps between and during solves, in storage the caller
owns. Nothing is allocated: an embedded application declares it static, and
sizeof(Workspace) is the solver's memory. The plan is the only part that
means anything outside a solve; the rest is scratch.
*/
template <class Problem, int N = Problem::Stages>
struct Workspace {
    using Traits = ProblemTraits<Problem>;
    using State = typename Traits::State;
    using Control = typename Traits::Control;
    using MatA = typename Traits::MatA;
    using MatB = typename Traits::MatB;

    /* The plan: states, commands and gains in the problem's units. */
    Policy<Problem, N> plan;

    /* Whether the plan came from a solve (true) or from reset (false). */
    bool solved = false;

    /* A trial plan, which becomes the plan when it is accepted. */
    State trial_x[N + 1];
    Control trial_u[N];

    /* Per stage, in normalised units: the step's Jacobians and the stage QP's step, which warm-starts the next one. */
    MatA A[N];
    MatB B[N];
    Control step[N];

    /* The bounds in force for this solve, read from the problem once. */
    Control lower[N];
    Control upper[N];
};

}

#endif
