/*
Part of NMPC. Distributed under the MIT licence: see LICENSE in the root of
the repository.
*/

#ifndef NMPC_WORKSPACE_H
#define NMPC_WORKSPACE_H

#include "NMPC/Policy.h"

namespace NMPC {

/*
Where the step's Jacobians live during a solve.

  Stored     every stage's A and B are made once an iteration, before the
             backward pass,
             and kept in the workspace. A backward pass repeated with a larger
             regularisation uses them again. This is the default.
  PerStage   each stage's A and B are made as the backward pass reaches that
             stage and
             are not kept. The workspace is smaller by N (NX^2 + NX NU) scalars,
             and a backward pass repeated with a larger regularisation makes
             them again.

The two give the same plans to the last bit: only where the Jacobians are kept
and how often they are made differ.
*/
enum class Jacobians { Stored, PerStage };

namespace detail {

template <class MatA, class MatB, int N, Jacobians J>
struct JacobianStore {
  MatA A[N];
  MatB B[N];
};

template <class MatA, class MatB, int N>
struct JacobianStore<MatA, MatB, N, Jacobians::PerStage> {};

}  // namespace detail

/*
Everything a solver keeps between and during solves, in storage the caller
owns. Nothing is allocated: an embedded application declares it static.
The solver's memory is sizeof(Workspace), sizeof(Solver) for the solver
object (the scales of the last solve and their products), and the call stack
of a solve. The plan is the only part that means anything outside a solve;
the rest is scratch.
*/
template <class Problem, int N = Problem::Stages,
          Jacobians J = Jacobians::Stored>
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

  /* Per stage, in normalised units: the step's Jacobians (with
     Jacobians::Stored) and the stage QP's step, which warm-starts the next one.
     */
  detail::JacobianStore<MatA, MatB, N, J> jacobians;
  Control step[N];

  /*
  The gains of the last backward pass, in the problem's units. They become the
  plan's only when the whole pass has completed with a usable prediction, so
  that a pass that fails part way never leaves the plan with new gains at its
  later stages and old ones at its earlier stages.
  */
  using Gain = typename Traits::Gain;
  Gain gain[N];

  /* The bounds in force for this solve, read from the problem once. */
  Control lower[N];
  Control upper[N];
};

}  // namespace NMPC

#endif
