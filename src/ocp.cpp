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

#include <cmath>
#include <qpDUNES.h>

#include "types.h"
#include "ocp.h"
#include "state.h"
#include "debug.h"

OptimalControlProblem::OptimalControlProblem(DynamicsModel *d) {
#if defined(NMPC_INTEGRATOR_RK4)
    integrator = IntegratorRK4();
#elif defined(NMPC_INTEGRATOR_HEUN)
    integrator = IntegratorHeun();
#elif defined(NMPC_INTEGRATOR_EULER)
    integrator = IntegratorEuler();
#endif

    dynamics = d;

    /* Initialise inequality constraints to +/-infinity. */
    lower_state_bound = StateConstraintVector::Ones() * -NMPC_INFTY;
    upper_state_bound = StateConstraintVector::Ones() * NMPC_INFTY;

    lower_control_bound = ControlConstraintVector::Ones() * -NMPC_INFTY;
    upper_control_bound = ControlConstraintVector::Ones() * NMPC_INFTY;

    /* Initialise weight matrices. */
    state_weights = StateWeightMatrix::Identity();
    control_weights = ControlWeightMatrix::Identity();
    terminal_weights = StateWeightMatrix::Identity();

    qp_options = qpDUNES_setupDefaultOptions();
    qp_options.maxIter = 100;
    qp_options.printLevel = 2;
    qp_options.stationarityTolerance = 1e-6;
}

GradientVector OptimalControlProblem::state_to_delta(
const StateVector &s, const ControlVector &c, const ReferenceVector &r) {
    GradientVector delta;

    delta.segment<6>(0) = s.segment<6>(0) - r.segment<6>(0);

    /*
    In order to increase the linearity of the problem and avoid quaternion
    normalisation issues, we calculate the difference between attitudes as a
    3-vector of Modified Rodrigues Parameters (MRP).
    */
    Quaternionr err_q = (Quaternionr(r.segment<4>(6)) *
        Quaternionr(s.segment<4>(6)).conjugate());

    delta.segment<3>(6) = NMPC_MRP_F *
        (err_q.vec() / (NMPC_MRP_A + err_q.w()));

    delta.segment<6>(9) = s.segment<6>(10) - r.segment<6>(10);
    delta.segment<NMPC_CONTROL_DIM>(NMPC_DELTA_DIM) = 
        c - r.segment<NMPC_CONTROL_DIM>(NMPC_STATE_DIM);

    return delta;
}

DeltaVector OptimalControlProblem::state_to_delta(
const StateVector &s1, const StateVector &s2) {
    DeltaVector delta;

    delta.segment<6>(0) = s1.segment<6>(0) - s2.segment<6>(0);

    /*
    In order to increase the linearity of the problem and avoid quaternion
    normalisation issues, we calculate the difference between attitudes as a
    3-vector of Modified Rodrigues Parameters (MRP).
    */
    Quaternionr err_q = (Quaternionr(s2.segment<4>(6)) *
        Quaternionr(s1.segment<4>(6)).conjugate());

    delta.segment<3>(6) = NMPC_MRP_F *
        (err_q.vec() / (NMPC_MRP_A + err_q.w()));

    delta.segment<6>(9) = s1.segment<6>(10) - s2.segment<6>(10);

    return delta;
}

/*
Linearises the state and control horizon around the previous solution, to
allow the QP solver to further optimise control values over the horizon.
Also calculates the gradient vector, which is just each delta multiplied by
the relevant weight matrix (the final one being multiplied by the terminal
weight).
*/
void OptimalControlProblem::calculate_deltas() {
    uint32_t i;

    for(i = 0; i < OCP_HORIZON_LENGTH; i++ ) {
        deltas[i] = state_to_delta(
            state_horizon[i],
            control_horizon[i],
            reference_trajectory[i]);

        /*
        Calculates the gradient vector, which is the difference between each
        point on the state/control horizon and the corresponding point in the
        reference trajectory multiplied by the weights. The terminal weights
        are applied to the last state delta.
        */
        StateWeightMatrix *weights;
        if(i != OCP_HORIZON_LENGTH-1) {
            weights = &state_weights;
        } else {
            weights = &terminal_weights;
        }

        gradient[i].segment<NMPC_DELTA_DIM>(0) =
            *weights * deltas[i].segment<NMPC_DELTA_DIM>(0);

        gradient[i].segment<NMPC_CONTROL_DIM>(NMPC_DELTA_DIM) =
            control_weights *
            deltas[i].segment<NMPC_CONTROL_DIM>(NMPC_DELTA_DIM);
    }
}

/*
Solve the initial value problems in order to set up continuity constraints,
which effectively store the system dynamics for this SQP iteration.
At the same time, compute the Jacobian function by applying perturbations to
each of the variables in turn, for use in the continuity constraints.
*/
void OptimalControlProblem::solve_ivps() {
    uint32_t i, j;

    for(i = 0; i < OCP_HORIZON_LENGTH-1; i++) {
        /* Solve the initial value problem at this horizon step. */
        integrated_state_horizon[i] = integrator.integrate(
            State(state_horizon[i]),
            control_horizon[i],
            dynamics,
            OCP_STEP_LENGTH);

        for(j = 0; j < NMPC_GRADIENT_DIM; j++) {
            ReferenceVector perturbed_state;
            perturbed_state.segment<NMPC_STATE_DIM>(0) = state_horizon[i];
            perturbed_state.segment<NMPC_CONTROL_DIM>(NMPC_STATE_DIM) =
                control_horizon[i];
            StateVector new_state;

            /* Need to calculate quaternion perturbations using MRPs. */
            if(j >= 6 || j <= 8) {
                Vector3r d_p;
                d_p << 0.0, 0.0, 0.0;
                d_p[j-6] = NMPC_EPS_4RT;
                real_t x_2 = d_p.squaredNorm();
                real_t delta_w = (-NMPC_MRP_A * x_2 + NMPC_MRP_F * std::sqrt(
                    NMPC_MRP_F_2 + ((real_t)1.0 - NMPC_MRP_A_2) * x_2)) /
                    (NMPC_MRP_F_2 + x_2);
                Vector3r delta_xyz = (((real_t)1.0 / NMPC_MRP_F) *
                    (NMPC_MRP_A + delta_w)) * d_p;
                Quaternionr delta_q;
                delta_q.vec() = delta_xyz;
                delta_q.w() = delta_w;
                Quaternionr temp = delta_q *
                    Quaternionr(perturbed_state.segment<4>(6));
                perturbed_state.segment<4>(6) << temp.vec(), temp.w();
            } else {
                perturbed_state[j] += NMPC_EPS_4RT;
            }

            new_state.segment<NMPC_STATE_DIM>(0) = integrator.integrate(
                State(perturbed_state.segment<NMPC_STATE_DIM>(0)),
                perturbed_state.segment<NMPC_CONTROL_DIM>(NMPC_STATE_DIM),
                dynamics,
                OCP_STEP_LENGTH);

            /*
            Calculate delta between perturbed state and original state, to
            yield a full column of the Jacobian matrix.
            */
            jacobians[i].col(j) =
                state_to_delta(new_state, integrated_state_horizon[i]) /
                NMPC_EPS_4RT;
        }

    /*
    Calculate integration residuals; these are needed for the continuity
    constraints.
    */
    integration_residuals[i] = state_to_delta(
        integrated_state_horizon[i],
        state_horizon[i+1]);
    }
}

/*
Linearises the objective function around the current horizon by calculating
second-order sensitivities at each step on the horizon.
The forward-difference method for calculating second-order derivatives from
gradient calls used in this function is due to "Numerical Methods for
Unconstrained Optimization and Nonlinear Equations" by J. E. Dennis, Jr., and
Robert B. Schnabel, page 103.
This could probably be improved a lot by following the methods outlined in
"Efficient Reduced SQP Methods for the Optimization of Chemical Processes
Described by Large Sparse DAE Models" by Daniel B. Leineweber.
*/
void OptimalControlProblem::calculate_hessians() {
    uint32_t i, j, k, l;
    real_t base_gradient;
    GradientVector cost_gradients, integrated_delta;

    for(i = 0; i < OCP_HORIZON_LENGTH; i++) {
        /*
        First, calculate cost function gradient with no perturbation. This is
        easy since we've already solved the IVPs with perturbation; just
        convert the integrated state vector to a delta.
        */
        StateWeightMatrix *weights;
        if(i != OCP_HORIZON_LENGTH-1) {
            GradientVector integrated_delta = state_to_delta(
            integrated_state_horizon[i],
            control_horizon[i],
            reference_trajectory[i+1]);

            weights = &state_weights;
        } else {
            GradientVector integrated_delta = state_to_delta(
            integrated_state_horizon[i],
            control_horizon[i],
            reference_trajectory[i]);

            weights = &terminal_weights;
        }

        real_t si, sf, ci, cf;

        si = integrated_delta.segment<NMPC_DELTA_DIM>(0).transpose() *
            *weights * integrated_delta.segment<NMPC_DELTA_DIM>(0);

        sf = deltas[i].segment<NMPC_DELTA_DIM>(0).transpose() *
            *weights * deltas[i].segment<NMPC_DELTA_DIM>(0);

        ci = integrated_delta.segment<NMPC_CONTROL_DIM>(
            NMPC_DELTA_DIM).transpose() *
            control_weights *
            integrated_delta.segment<NMPC_CONTROL_DIM>(NMPC_DELTA_DIM);

        cf = deltas[i].segment<NMPC_CONTROL_DIM>(
            NMPC_DELTA_DIM).transpose() *
            control_weights *
            deltas[i].segment<NMPC_CONTROL_DIM>(NMPC_DELTA_DIM);

        base_gradient = (sf + cf) - (si + ci);

        /* Next, calculate cost function gradient for each parameter. */
        for(j = 0; j < NMPC_GRADIENT_DIM; j++) {
            if(deltas[i](j) >= NMPC_EPS_SQRT) {
                /* Perturb the delta vector by an appropriate amount. */
                GradientVector perturbed_delta = deltas[i];
                deltas[i](j) += std::abs(deltas[i](j))*NMPC_EPS_4RT;

                /*
                Apply the perturbed delta vector to the reference trajectory
                vector (deltas are the difference between reference points
                and points on the state horizon).
                */
                StateVector perturbed_state;
                perturbed_state.segment<6>(0) =
                    reference_trajectory[i].segment<6>(0) +
                    perturbed_delta.segment<6>(0);

                /* Apply the delta MRP to the state vector quaternion. */
                Vector3r d_p = perturbed_delta.segment<3>(6);
                real_t x_2 = d_p.squaredNorm();
                real_t delta_w = (-NMPC_MRP_A * x_2 + NMPC_MRP_F * std::sqrt(
                    NMPC_MRP_F_2 + ((real_t)1.0 - NMPC_MRP_A_2) * x_2)) /
                    (NMPC_MRP_F_2 + x_2);
                Vector3r delta_xyz = (((real_t)1.0 / NMPC_MRP_F) *
                    (NMPC_MRP_A + delta_w)) * d_p;
                Quaternionr delta_q;
                delta_q.vec() = delta_xyz;
                delta_q.w() = delta_w;
                Quaternionr temp = delta_q *
                    Quaternionr(perturbed_state.segment<4>(6));
                perturbed_state.segment<4>(6) << temp.vec(), temp.w();

                perturbed_state.segment<6>(10) =
                    reference_trajectory[i].segment<6>(10) +
                    perturbed_delta.segment<6>(9);
                ControlVector perturbed_control = 
                    reference_trajectory[i].segment<NMPC_CONTROL_DIM>(
                        NMPC_STATE_DIM) +
                    perturbed_delta.segment<NMPC_CONTROL_DIM>(
                        NMPC_DELTA_DIM);

                /* Integrate the perturbed state vector. */
                StateVector integrated_state = integrator.integrate(
                    State(perturbed_state),
                    perturbed_control,
                    dynamics,
                    OCP_STEP_LENGTH);

                /* Convert the integrated result back into a delta. */
                if(i != OCP_HORIZON_LENGTH-1) {
                    integrated_delta = state_to_delta(
                    integrated_state,
                    perturbed_control,
                    reference_trajectory[i+1]);
                } else {
                    integrated_delta = state_to_delta(
                    integrated_state,
                    perturbed_control,
                    reference_trajectory[i]);
                }

                /*
                Calculate the cost function on both the initial and final
                deltas.
                */
                si = integrated_delta.segment<NMPC_DELTA_DIM>(
                    0).transpose() *
                    *weights * integrated_delta.segment<NMPC_DELTA_DIM>(0);

                sf = perturbed_delta.segment<NMPC_DELTA_DIM>(
                    0).transpose() *
                    *weights * perturbed_delta.segment<NMPC_DELTA_DIM>(0);

                ci = integrated_delta.segment<NMPC_CONTROL_DIM>(
                    NMPC_DELTA_DIM).transpose() *
                    control_weights *
                    integrated_delta.segment<NMPC_CONTROL_DIM>(
                        NMPC_DELTA_DIM);

                cf = perturbed_delta.segment<NMPC_CONTROL_DIM>(
                    NMPC_DELTA_DIM).transpose() *
                    control_weights *
                    perturbed_delta.segment<NMPC_CONTROL_DIM>(
                        NMPC_DELTA_DIM);

                cost_gradients[j] =
                    (((sf + cf) - (si + ci)) - base_gradient) /
                    (2.0*std::abs(deltas[i](j))*NMPC_EPS_4RT);
            } else {
                cost_gradients[j] = 0;
            }
        }

        /*
        Using the evaluated gradient calls, we can now assemble a finite-
        difference approximation of the Hessian matrix using forward
        differences.
        */
        hessians[i].diagonal() = cost_gradients * (real_t)2.0;
        for(k = 0; k < NMPC_GRADIENT_DIM-1; k++) {
            for(l = k+1; l < NMPC_GRADIENT_DIM; l++) {
                hessians[i](k, l) = cost_gradients[k] + cost_gradients[l];
            }
        }
    }
}

/*
Uses all of the information calculated so far to set up the various qpDUNES
datastructures in preparation for the feedback step.
This is really inefficient right now – there's heaps of copying and
transpositions because the Eigen matrices are all in column-major format
whereas qpDUNES expects row-major arrays.
Could do this much more efficiently using the Eigen Map class, and possibly
avoid having to copy data at all.
*/
void OptimalControlProblem::initialise_qp() {
    uint32_t i;

    /* Set up problem dimensions. */
    /* TODO: Determine number of affine constraints (D), and add them. */
    qpDUNES_setup(
        &qp_data,
        OCP_HORIZON_LENGTH,
        NMPC_DELTA_DIM,
        NMPC_CONTROL_DIM,
        0,
        &qp_options);

    return_t status_flag;

    for(i = 0; i < OCP_HORIZON_LENGTH; i++) {
//        status_flag = qpDUNES_setupRegularInterval(
//            &qp_data,
//            );
        AssertOK(status_flag);
    }
}

/*
This step is the first part of the feedback step; the very latest sensor
measurement should be provided in order to to set up the initial state for
the SQP iteration. This allows the feedback delay to be significantly less
than one time step.
*/
void OptimalControlProblem::initial_constraint(ReferenceVector measurement) {

}

/* Solves the QP uses qpDUNES. */
void OptimalControlProblem::solve_qp() {
    return_t status_flag = qpDUNES_solve(&qp_data);
    AssertSolutionFound(status_flag);
}

void OptimalControlProblem::update_horizon() {

}

/* Copies the reference trajectory into the state and control horizons. */
void OptimalControlProblem::initialise() {
    uint32_t i;

    for(i = 0; i < OCP_HORIZON_LENGTH; i++) {
        state_horizon[i] =
            reference_trajectory[i].segment<NMPC_STATE_DIM>(0);
        control_horizon[i] =
            reference_trajectory[i].segment<NMPC_CONTROL_DIM>(
                NMPC_STATE_DIM);
    }
}

/*
Completes the most computationally intense part of the NMPC iteration; this
step is independent of the lastest sensor measurements and so can be
executed as soon as possible after the previous iteration.
*/
void OptimalControlProblem::preparation_step() {
    calculate_deltas();
    solve_ivps();
    calculate_hessians();
    initialise_qp();
}

/*
This step feeds the latest sensor measurements into the problem and runs the
QP solver. This step uses takes much less time than the previous step, and so
should be executed at the last possible moment with the latest sensor data,
in order to make the control latency significantly less than the horizon step
length.
*/
void OptimalControlProblem::feedback_step(ReferenceVector measurement) {
    initial_constraint(measurement);
    solve_qp();
}
