/*
Copyright (C) 2013 Daniel Dyer

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef OCP_H
#define OCP_H

#include <stdint.h>
#include <cmath>
#include <qpDUNES.h>

#include "types.h"
#include "state.h"
#include "integrator.h"

/* OCP control and prediction horizon (number of steps). */
#define OCP_HORIZON_LENGTH 500

/* OCP control step length (seconds). */
#define OCP_STEP_LENGTH (1.0/50.0)

/*
Optimal Control Problem object.
*/
class OptimalControlProblem {
    /* Integrator object, depends on selection in `config.h`. */
#if defined(NMPC_INTEGRATOR_RK4)
    IntegratorRK4 integrator;
#elif defined(NMPC_INTEGRATOR_HEUN)
    IntegratorHeun integrator;
#elif defined(NMPC_INTEGRATOR_EULER)
    IntegratorEuler integrator;
#endif
    
    DynamicsModel *dynamics;

    ReferenceVector reference_trajectory[OCP_HORIZON_LENGTH];
    ControlVector control_horizon[OCP_HORIZON_LENGTH];
    StateVector state_horizon[OCP_HORIZON_LENGTH];
    StateVector integrated_state_horizon[OCP_HORIZON_LENGTH];

    /*
    Difference between predicted state and reference trajectory for each
    point on the horizon, weighted by the weight matrices.
    */
    GradientVector gradient[OCP_HORIZON_LENGTH];

    /*
    The "deltas" are the set of perturbations from point around which the
    system has been linearised, which arise due to the feedback step (actual
    measurement data being supplied in each iteration). These are the vectors
    which the QP solver is trying to optimise.
    */
    GradientVector deltas[OCP_HORIZON_LENGTH];

    /*
    Hessian matrices for each point on the horizon. Recalculated each
    iteration for use by the QP solver.
    */
    HessianMatrix hessians[OCP_HORIZON_LENGTH];

    /*
    Constraint matrix and bounding vectors. Note that we could easily have
    a constraint matrix, upper and lower bound vector for each horizon step,
    and they could be recalculated each iteration from a set of nonlinear
    constraints if desired.

    For now, only one constraint matrix and set of bounding vectors will be
    used to reduce memory usage.

    May need to generate constraints per-step to handle airspeed limits,
    though.
    */
    InequalityConstraintMatrix inequality_constraints;
    InequalityConstraintVector upper_bound, lower_bound;

    /*
    Matrices for continuity constraints. These are actually Jacobian
    matrices, which contain the linearised dynamics model for each point on
    the horizon.
    */
    ContinuityConstraintMatrix jacobians[OCP_HORIZON_LENGTH];
    DeltaVector integration_residuals[OCP_HORIZON_LENGTH];

    /* Weight matrices. */
    StateWeightMatrix state_weights;
    ControlWeightMatrix control_weights;
    StateWeightMatrix terminal_weights;

    /* Data structures for use by qpDUNES. */
    qpData_t qp_data;
    qpOptions_t qp_options;

    GradientVector state_to_delta(
        const StateVector &s, const ControlVector &c, const ReferenceVector &r);
    DeltaVector state_to_delta(const StateVector &s1, const StateVector &s2);
    void calculate_deltas();
    void solve_ivps();
    void calculate_hessians();
    void initialise_qp();
    void initial_constraint(ReferenceVector measurement);
    void solve_qp();
    void update_horizon();

public:
    OptimalControlProblem();
    void initialise();
    void preparation_step();
    void feedback_step(ReferenceVector measurement);
    void set_dynamics_model(DynamicsModel *in) { dynamics = in; }
};

#endif
