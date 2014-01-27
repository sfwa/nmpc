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

extern "C"
{
    #include <qpDUNES.h>
}

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
    qp_options.printLevel = 0;
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
    Quaternionr err_q = (Quaternionr(s.segment<4>(6)) *
        Quaternionr(r.segment<4>(6)).conjugate());

    if(err_q.w() < 0) {
        err_q = Quaternionr(-err_q.w(), -err_q.x(), -err_q.y(), -err_q.z());
    }

    delta.segment<3>(6) = NMPC_MRP_F *
        (err_q.vec() / (NMPC_MRP_A + err_q.w()));

    delta.segment<3>(9) = s.segment<3>(10) - r.segment<3>(10);
    delta.segment<NMPC_CONTROL_DIM>(NMPC_DELTA_DIM) = 
        c - r.segment<NMPC_CONTROL_DIM>(NMPC_STATE_DIM);

    return delta;
}

DeltaVector OptimalControlProblem::state_to_delta(
const StateVector &s1, const StateVector &s2) {
    DeltaVector delta;

    delta.segment<6>(0) = s2.segment<6>(0) - s1.segment<6>(0);

    /*
    In order to increase the linearity of the problem and avoid quaternion
    normalisation issues, we calculate the difference between attitudes as a
    3-vector of Modified Rodrigues Parameters (MRP).
    */
    Quaternionr err_q = (Quaternionr(s2.segment<4>(6)) *
        Quaternionr(s1.segment<4>(6)).conjugate());

    if(err_q.w() < 0) {
        err_q = Quaternionr(-err_q.w(), -err_q.x(), -err_q.y(), -err_q.z());
    }

    delta.segment<3>(6) = NMPC_MRP_F *
        (err_q.vec() / (NMPC_MRP_A + err_q.w()));

    delta.segment<3>(9) = s2.segment<3>(10) - s1.segment<3>(10);

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

        // std::cout << deltas[i].transpose() << std::endl;

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

        gradients[i].segment<NMPC_DELTA_DIM>(0) =
            *weights * deltas[i].segment<NMPC_DELTA_DIM>(0);

        gradients[i].segment<NMPC_CONTROL_DIM>(NMPC_DELTA_DIM) =
            control_weights *
            deltas[i].segment<NMPC_CONTROL_DIM>(NMPC_DELTA_DIM);
    }

    // std::cout << "++++++++++++++++" << std::endl << std::endl;
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
            real_t perturbation = NMPC_EPS_4RT;

            /* Need to calculate quaternion perturbations using MRPs. */
            if(j < 6) {
                perturbed_state[j] += perturbation;
            }
            else if(j >= 6 && j <= 8) {
                Vector3r d_p;
                d_p << 0.0, 0.0, 0.0;
                d_p[j-6] = perturbation;
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
            } else if(j < NMPC_DELTA_DIM) {
                perturbed_state[j+1] += perturbation;
            } else {
                /*
                Perturbations for the control inputs should be proportional
                to the control range to make sure we don't lose too much
                precision.
                */
                perturbation *=
                    (upper_control_bound[j-NMPC_DELTA_DIM] -
                    lower_control_bound[j-NMPC_DELTA_DIM]);
                perturbed_state[j+1] += perturbation;
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
                state_to_delta(integrated_state_horizon[i], new_state) /
                perturbation;
        }

        /*
        Calculate integration residuals; these are needed for the continuity
        constraints.
        */
        integration_residuals[i] = state_to_delta(
            state_horizon[i+1],
            integrated_state_horizon[i]);
    }
}

/*
Uses all of the information calculated so far to set up the various qpDUNES
datastructures in preparation for the feedback step.
This is really inefficient right now – there's heaps of probably unnecessary
copying going on.
*/
void OptimalControlProblem::initialise_qp() {
    uint32_t i;
    real_t Q[NMPC_DELTA_DIM*NMPC_DELTA_DIM];
    Eigen::Map<StateWeightMatrix> Q_map(Q);
    real_t R[NMPC_CONTROL_DIM*NMPC_CONTROL_DIM];
    Eigen::Map<ControlWeightMatrix> R_map(R);
    real_t P[NMPC_DELTA_DIM*NMPC_DELTA_DIM];
    Eigen::Map<StateWeightMatrix> P_map(P);
    real_t g[NMPC_GRADIENT_DIM];
    Eigen::Map<GradientVector> g_map(g);
    real_t C[(NMPC_STATE_DIM-1)*NMPC_GRADIENT_DIM];
    Eigen::Map<ContinuityConstraintMatrix> C_map(C);
    real_t c[NMPC_DELTA_DIM];
    Eigen::Map<DeltaVector> c_map(c);
    real_t zLow[NMPC_GRADIENT_DIM];
    Eigen::Map<GradientVector> zLow_map(zLow);
    real_t zUpp[NMPC_GRADIENT_DIM];
    Eigen::Map<GradientVector> zUpp_map(zUpp);

    zLow_map.segment<NMPC_DELTA_DIM>(0) = lower_state_bound;
    zUpp_map.segment<NMPC_DELTA_DIM>(0) = upper_state_bound;

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
        /* Copy the relevant data into the qpDUNES arrays. */
        Q_map = state_weights;
        R_map = control_weights;
        g_map = gradients[i];
        C_map = jacobians[i];
        c_map = integration_residuals[i];
        zLow_map.segment<NMPC_CONTROL_DIM>(NMPC_DELTA_DIM) =
            lower_control_bound - control_horizon[i];
        zUpp_map.segment<NMPC_CONTROL_DIM>(NMPC_DELTA_DIM) =
            upper_control_bound - control_horizon[i];

        status_flag = qpDUNES_setupRegularInterval(
            &qp_data, qp_data.intervals[i],
            0, Q, R, 0, g, C, 0, 0, c, zLow, zUpp, 0, 0, 0, 0, 0, 0, 0);
        AssertOK(status_flag);
    }

    /* Set up final interval. */
    P_map = terminal_weights;
    status_flag = qpDUNES_setupFinalInterval(&qp_data, qp_data.intervals[i],
        P, g, zLow, zUpp, 0, 0, 0);
    AssertOK(status_flag);

    qpDUNES_setupAllLocalQPs(&qp_data, QPDUNES_FALSE);

    qpDUNES_indicateDataChange(&qp_data);
}

/*
This step is the first part of the feedback step; the very latest sensor
measurement should be provided in order to to set up the initial state for
the SQP iteration. This allows the feedback delay to be significantly less
than one time step.
*/
void OptimalControlProblem::initial_constraint(StateVector measurement) {
    real_t zLow[NMPC_GRADIENT_DIM];
    Eigen::Map<GradientVector> zLow_map(zLow);
    real_t zUpp[NMPC_GRADIENT_DIM];
    Eigen::Map<GradientVector> zUpp_map(zUpp);

    /* Control constraints are unchanged. */
    zLow_map.segment<NMPC_CONTROL_DIM>(NMPC_DELTA_DIM) =
        lower_control_bound - control_horizon[0];
    zUpp_map.segment<NMPC_CONTROL_DIM>(NMPC_DELTA_DIM) =
        upper_control_bound - control_horizon[0];

    /*
    Initial delta is constrained to be the difference between the measurement
    and the initial state horizon point.
    */
    DeltaVector initial_delta = state_to_delta(
        state_horizon[0],
        measurement);
    zLow_map.segment<NMPC_DELTA_DIM>(0) = initial_delta;
    zUpp_map.segment<NMPC_DELTA_DIM>(0) = initial_delta;

    return_t status_flag = qpDUNES_updateIntervalData(
        &qp_data, qp_data.intervals[0],
        0, 0, 0, 0, zLow, zUpp, 0, 0, 0, 0);
    AssertOK(status_flag);

    qpDUNES_indicateDataChange(&qp_data);
}

/* Solves the QP using qpDUNES. */
void OptimalControlProblem::solve_qp() {
    uint32_t i;
    real_t solution[NMPC_GRADIENT_DIM*(OCP_HORIZON_LENGTH+1)];

    return_t status_flag = qpDUNES_solve(&qp_data);
    AssertSolutionFound(status_flag);

    /* Get the solution. */
    qpDUNES_getPrimalSol(&qp_data, solution);

    /*
    Apply the deltas to the reference trajectory to generate the new state
    horizon.
    */
    for(i = 0; i < OCP_HORIZON_LENGTH; i++) {
        Eigen::Map<GradientVector> solution_map(
            &solution[i*NMPC_GRADIENT_DIM]);

        // std::cout << gradients[i].transpose() << std::endl;
        // std::cout << solution_map.transpose() << std::endl;

        state_horizon[i].segment<6>(0) += solution_map.segment<6>(0);

        Vector3r d_p = solution_map.segment<3>(6);
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
            Quaternionr(state_horizon[i].segment<4>(6));
        state_horizon[i].segment<4>(6) << temp.vec(), temp.w();

        state_horizon[i].segment<3>(10) += solution_map.segment<3>(9);

        control_horizon[i] +=
            solution_map.segment<NMPC_CONTROL_DIM>(NMPC_DELTA_DIM);

        // std::cout << reference_trajectory[i].transpose() << std::endl;
        // std::cout << state_horizon[i].transpose() << "\t" << control_horizon[i].transpose() << std::endl << std::endl;
        // std::cout << control_horizon[i].transpose() << std::endl;
    }

    // std::cout << "=========" << std::endl << std::endl;
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
    initialise_qp();
}

/*
This step feeds the latest sensor measurements into the problem and runs the
QP solver. This step uses takes much less time than the previous step, and so
should be executed at the last possible moment with the latest sensor data,
in order to make the control latency significantly less than the horizon step
length.
*/
void OptimalControlProblem::feedback_step(StateVector measurement) {
    initial_constraint(measurement);
    solve_qp();
}

/*
Shift the horizon across and add a new point to the end of the reference
trajectory.
*/
void OptimalControlProblem::update_horizon(ReferenceVector new_reference) {
    uint32_t i;

    /* TODO: Do this more intelligently. */
    for(i = 0; i < OCP_HORIZON_LENGTH-1; i++) {
        state_horizon[i] = state_horizon[i+1];
        control_horizon[i] = control_horizon[i+1];
        reference_trajectory[i] = reference_trajectory[i+1];
    }

    /* Insert new reference point. */
    reference_trajectory[i] = new_reference;
    state_horizon[i] = new_reference.segment<NMPC_STATE_DIM>(0);
    control_horizon[i] = new_reference.segment<NMPC_CONTROL_DIM>(
        NMPC_STATE_DIM);

    /* Prepare the QP for the next solution. */
    qpDUNES_shiftLambda(&qp_data);
    qpDUNES_shiftIntervals(&qp_data);
}
