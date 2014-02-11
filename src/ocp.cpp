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
    qp_options.maxIter = 5;
    qp_options.printLevel = 0;
    qp_options.stationarityTolerance = 1e-3;
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
Solve the initial value problems in order to set up continuity constraints,
which effectively store the system dynamics for this SQP iteration.
At the same time, compute the Jacobian function by applying perturbations to
each of the variables in turn, for use in the continuity constraints.
*/
void OptimalControlProblem::solve_ivps(uint32_t i) {
    uint32_t j;

    /* Solve the initial value problem at this horizon step. */
    integrated_state_horizon[i] = integrator.integrate(
        State(state_reference[i]),
        control_reference[i],
        dynamics,
        OCP_STEP_LENGTH);

    for(j = 0; j < NMPC_GRADIENT_DIM; j++) {
        ReferenceVector perturbed_state;
        perturbed_state.segment<NMPC_STATE_DIM>(0) = state_reference[i];
        perturbed_state.segment<NMPC_CONTROL_DIM>(NMPC_STATE_DIM) =
            control_reference[i];
        StateVector new_state;
        real_t perturbation = NMPC_EPS_4RT;

        /* Need to calculate quaternion perturbations using MRPs. */
        if(j < 6) {
            perturbed_state[j] += perturbation;
        } else if(j >= 6 && j <= 8) {
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

    /* Gradient vector fixed to zero. */
    g_map = GradientVector::Zero();

    /* Continuity constraint constant term fixed to zero. */
    c_map = DeltaVector::Zero();

    /* Zero Jacobians for now */
    C_map = ContinuityConstraintMatrix::Zero();

    Q_map = state_weights;
    R_map = control_weights;

    /* Copy the relevant data into the qpDUNES arrays. */
    zLow_map.segment<NMPC_CONTROL_DIM>(NMPC_DELTA_DIM) = lower_control_bound;
    zUpp_map.segment<NMPC_CONTROL_DIM>(NMPC_DELTA_DIM) = upper_control_bound;

    for(i = 0; i < OCP_HORIZON_LENGTH; i++) {
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
Updates the QP with the latest linearisations.
*/
void OptimalControlProblem::update_qp() {

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
        lower_control_bound - control_reference[0];
    zUpp_map.segment<NMPC_CONTROL_DIM>(NMPC_DELTA_DIM) =
        upper_control_bound - control_reference[0];

    /*
    Initial delta is constrained to be the difference between the measurement
    and the initial state horizon point.
    */
    DeltaVector initial_delta = state_to_delta(
        state_reference[0],
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

    if (status_flag == QPDUNES_SUCC_OPTIMAL_SOLUTION_FOUND) {
        /* Get the solution. */
        qpDUNES_getPrimalSol(&qp_data, solution);

        /* Extract the first set of control values */
        Eigen::Map<GradientVector> solution_map(
            &solution[0]);

        control_horizon[0] =
            control_reference[0] +
            solution_map.segment<NMPC_CONTROL_DIM>(NMPC_DELTA_DIM);
    }
}

/* Copies the reference trajectory into the state and control horizons. */
void OptimalControlProblem::initialise() {
    initialise_qp();
}

/*
Completes the most computationally intense part of the NMPC iteration; this
step is independent of the lastest sensor measurements and so can be
executed as soon as possible after the previous iteration.
*/
void OptimalControlProblem::preparation_step() {
    update_qp();
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
    memmove(state_reference, &state_reference[1],
            sizeof(StateVector) * (OCP_HORIZON_LENGTH - 1));
    memmove(control_reference, &control_reference[1],
            sizeof(ControlVector) * (OCP_HORIZON_LENGTH - 1));

    /* Prepare the QP for the next solution. */
    qpDUNES_shiftLambda(&qp_data);
    qpDUNES_shiftIntervals(&qp_data);

    set_reference_point(new_reference, OCP_HORIZON_LENGTH - 1);
}

void OptimalControlProblem::set_reference_point(const ReferenceVector &in,
uint32_t i) {
    state_reference[i] = in.segment<NMPC_STATE_DIM>(0);

    if(i < OCP_HORIZON_LENGTH) {
        control_reference[i] =
            in.segment<NMPC_CONTROL_DIM>(NMPC_STATE_DIM);

        solve_ivps(i);

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

        return_t status_flag;

        /* Gradient vector fixed to zero. */
        g_map = GradientVector::Zero();

        /* Continuity constraint constant term fixed to zero. */
        c_map = DeltaVector::Zero();

        /* Copy the relevant data into the qpDUNES arrays. */
        C_map = jacobians[i];
        zLow_map.segment<NMPC_CONTROL_DIM>(NMPC_DELTA_DIM) =
            lower_control_bound - control_reference[i];
        zUpp_map.segment<NMPC_CONTROL_DIM>(NMPC_DELTA_DIM) =
            upper_control_bound - control_reference[i];

        status_flag = qpDUNES_updateIntervalData(
            &qp_data, qp_data.intervals[i],
            0, g, C, 0, zLow, zUpp, 0, 0, 0, 0);
        AssertOK(status_flag);

        qpDUNES_indicateDataChange(&qp_data);
    }
}
