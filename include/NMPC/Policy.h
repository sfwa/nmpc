/*
Part of NMPC. Distributed under the MIT licence: see LICENSE in the root of
the repository.
*/

#ifndef NMPC_POLICY_H
#define NMPC_POLICY_H

#include <algorithm>

#include "NMPC/Problem.h"

namespace NMPC {

/*
The plan as a policy: for each stage the planned state at its start, the
command, and the gain that turns a departure from the planned state into a
change of command. It is the solver's output and the stored-gain inner loop's
input, and it is a plain value: an application that runs the solve below the
inner loop in priority copies it whole into a buffer that the inner loop
switches to on a tick boundary.

Gains are in the problem's own units: a change of command per unit of
tangent. Between stages the inner loop applies

    command = u(k) + K(k) * difference(x, x_plan)

where x_plan is the stage's planned state, ideally carried forward through the
stage by the model so that the gain acts about a moving point.

`origin` and `ticks_per_stage` are the application's: the tick at which stage 0
comes into force and the inner-loop ticks in a stage. The solver leaves them
alone. With them, at(tick) finds the stage in force at a tick, and holds the
last stage once the plan has run out.
*/
template <class Problem, int N = Problem::Stages>
struct Policy {
  using Traits = ProblemTraits<Problem>;
  using T = typename Traits::T;
  using State = typename Traits::State;
  using Control = typename Traits::Control;
  using Gain = typename Traits::Gain;

  enum { Stages = N };
  static_assert(N > 0, "A plan needs at least one stage");

  /* x[k] is the state at the start of stage k; x[N] is where the plan ends. */
  State x[N + 1];
  Control u[N];
  Gain K[N];

  long origin = 0;
  long ticks_per_stage = 1;

  /* A view of one stage. */
  struct Stage {
    int index;
    const State& x;
    const Control& u;
    const Gain& K;
  };

  int stage_index(long tick) const {
    const long elapsed = std::max(0L, tick - origin);
    return int(std::min(elapsed / std::max(1L, ticks_per_stage), long(N - 1)));
  }

  Stage at(long tick) const {
    const int k = stage_index(tick);
    return Stage{k, x[k], u[k], K[k]};
  }

  /* The command of stage k for the state x, before any bound. */
  Control feedback(int k, const State& state, const State& planned) const {
    NMPC_ASSERT(k >= 0 && k < N);
    return u[k] + K[k] * Problem::difference(state, planned);
  }

  Control feedback(int k, const State& state) const {
    return feedback(k, state, x[k]);
  }
};

}  // namespace NMPC

#endif
