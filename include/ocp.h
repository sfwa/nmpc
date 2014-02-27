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

    bool relative_positions;

    DynamicsModel *dynamics;

    ControlVector control_reference[OCP_HORIZON_LENGTH];
    StateVector state_reference[OCP_HORIZON_LENGTH+1];
    ControlVector control_horizon[OCP_HORIZON_LENGTH];
    StateVector state_horizon[OCP_HORIZON_LENGTH+1];
    StateVector integrated_state_horizon[OCP_HORIZON_LENGTH];

    /*
    Affine constraint matrix and bounding vectors. These will be generated
    each iteration from the non-linear constraints.
    */
    InequalityConstraintMatrix affine_constraints[OCP_HORIZON_LENGTH];
    InequalityConstraintVector affine_upper_bound[OCP_HORIZON_LENGTH];
    InequalityConstraintVector affine_lower_bound[OCP_HORIZON_LENGTH];

    /*
    Simple inequality constraints.
    */
    StateConstraintVector lower_state_bound, upper_state_bound;
    ControlConstraintVector lower_control_bound, upper_control_bound;

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

    DeltaVector state_to_delta(
        const StateVector &s1,
        const StateVector &s2);
    void calculate_gradient();
    void solve_ivps(uint32_t i);
    void initialise_qp();
    void update_qp();
    void initial_constraint(StateVector measurement);
    void solve_qp();

public:
    OptimalControlProblem(DynamicsModel *d, bool use_relative_positions);
    void initialise();
    void set_state_weights(const DeltaVector &in) {
        state_weights.diagonal() = in;
    }
    void set_control_weights(const ControlVector &in) {
        control_weights.diagonal() = in;
    }
    void set_terminal_weights(const DeltaVector &in) {
        terminal_weights.diagonal() = in;
    }
    void set_lower_control_bound(const ControlConstraintVector &in) {
        lower_control_bound = in;
    }
    void set_upper_control_bound(const ControlConstraintVector &in) {
        upper_control_bound = in;
    }
    void set_reference_point(const ReferenceVector &in, uint32_t i);
    void preparation_step();
    void feedback_step(StateVector measurement);
    const ControlVector& get_controls() const { return control_horizon[0]; }
    void update_horizon(ReferenceVector new_reference);
    void set_dynamics_model(DynamicsModel *in) { dynamics = in; }
};

#endif
