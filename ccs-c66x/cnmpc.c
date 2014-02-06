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

#include "types.h"
#include "qpDUNES.h"

#ifndef __TI_COMPILER_VERSION__
    #include "config.h"
    #include "../c/cnmpc.h"
#else
    #define UKF_USE_DSP_INTRINSICS

    #include "config.h"
    #include "cnmpc.h"
#endif

#define X 0
#define Y 1
#define Z 2
#define W 3

/* Non-TI compatibility */
#ifndef __TI_COMPILER_VERSION__
#define _nassert(x)
#endif

#ifndef M_PI
#define M_PI ((real_t)3.14159265358979323846)
#define M_PI_2 (M_PI * 0.5)
#define M_PI_4 (M_PI * 0.25)
#endif

#ifdef NMPC_SINGLE_PRECISION
#define NMPC_EPS_SQRT ((real_t)3.450e-4)
#define NMPC_EPS_4RT ((real_t)1.857e-2)
#endif

#ifdef NMPC_DOUBLE_PRECISION
#define NMPC_EPS_SQRT ((real_t)1.490e-8)
#define NMPC_EPS_4RT ((real_t)1.221e-4)
#endif

#define NMPC_INFTY ((real_t)1.0e12)

#define G_ACCEL ((real_t)9.80665)
#define RHO ((real_t)1.225)

/* Math routines -- assume float */
static inline void quaternion_multiply(real_t *restrict res,
const real_t q1[4], const real_t q2[4]) {
    assert(res && q1 && q2 && res != q1 && res != q2);
    _nassert((size_t)res % 4 == 0);
    _nassert((size_t)q1 % 4 == 0);
    _nassert((size_t)q2 % 4 == 0);

    res[W] = q1[W]*q2[W] - q1[X]*q2[X] - q1[Y]*q2[Y] - q1[Z]*q2[Z];
    res[X] = q1[W]*q2[X] + q1[X]*q2[W] + q1[Y]*q2[Z] - q1[Z]*q2[Y];
    res[Y] = q1[W]*q2[Y] - q1[X]*q2[Z] + q1[Y]*q2[W] + q1[Z]*q2[X];
    res[Z] = q1[W]*q2[Z] + q1[X]*q2[Y] - q1[Y]*q2[X] + q1[Z]*q2[W];
}

static inline void quaternion_conjugate(real_t *restrict res,
const real_t q[4]) {
    assert(res && q && res != q);
    _nassert((size_t)res % 4 == 0);
    _nassert((size_t)q % 4 == 0);

    res[X] = q[X];
    res[Y] = q[Y];
    res[Z] = q[Z];
    res[W] = -q[W];
}

static inline void vector3_add(real_t *restrict res, const real_t v1[3],
const real_t v2[3]) {
    assert(res && v1 && v2 && res != v1 && res != v2);
    _nassert((size_t)res % 4 == 0);
    _nassert((size_t)v1 % 4 == 0);
    _nassert((size_t)v2 % 4 == 0);

    res[X] = v1[X] + v2[X];
    res[Y] = v1[Y] + v2[Y];
    res[Z] = v1[Z] + v2[Z];
}

static inline void vector3_subtract(real_t *restrict res, const real_t v1[3],
const real_t v2[3]) {
    assert(res && v1 && v2 && res != v1 && res != v2);
    _nassert((size_t)res % 4 == 0);
    _nassert((size_t)v1 % 4 == 0);
    _nassert((size_t)v2 % 4 == 0);

    res[X] = v1[X] - v2[X];
    res[Y] = v1[Y] - v2[Y];
    res[Z] = v1[Z] - v2[Z];
}

static real_t wind_velocity[3];

/* 32000B */
static real_t ocp_state_reference[OCP_HORIZON_LENGTH * NMPC_REFERENCE_DIM];

/* 6000B */
static real_t ocp_control_reference[OCP_HORIZON_LENGTH * NMPC_CONTROL_DIM];

static real_t ocp_lower_state_bound[NMPC_DELTA_DIM];
static real_t ocp_upper_state_bound[NMPC_DELTA_DIM];
static real_t ocp_lower_control_bound[NMPC_CONTROL_DIM];
static real_t ocp_upper_control_bound[NMPC_CONTROL_DIM];
static real_t ocp_state_weights[NMPC_DELTA_DIM]; /* diagonal only */
static real_t ocp_terminal_weights[NMPC_DELTA_DIM]; /* diagonal only */
static real_t ocp_control_weights[NMPC_CONTROL_DIM]; /* diagonal only */

static qpData_t qp_data;
static qpOptions_t qp_options;


void _state_to_delta(real_t *delta, const real_t *restrict s1,
const real_t *restrict s2) {
    assert(delta && s1 && s2);
    assert(delta != s1);
    assert(delta != s2);
    assert(s1 != s2);

    size_t i;

    /* Calculate deltas for position, velocity, and angular velocity */
    #pragma MUST_ITERATE(3, 3)
    for (i = 0; i < 3; i++) {
        delta[i] = s2[i] - s1[i];
        delta[i + 3] = s2[i + 3] - s1[i + 3];
        delta[i + 9] = s2[i + 10] - s1[i + 10];
    }

    /*
    In order to increase the linearity of the problem and avoid quaternion
    normalisation issues, we calculate the difference between attitudes as a
    3-vector of Modified Rodrigues Parameters (MRP).
    */
    real_t err_q[4], s1_conj[4];
    quaternion_conjugate(s1_conj, &s1[6]);
    quaternion_multiply(err_q, &s2[6], s1_conj);

    /* Ensure real part stays positive */
    if (err_q[W] < 0) {
        err_q[X] = -err_q[X];
        err_q[Y] = -err_q[Y];
        err_q[Z] = -err_q[Z];
        err_q[W] = -err_q[W];
    }

    real_t d = NMPC_MRP_F / (NMPC_MRP_A + err_q[W]);
    delta[6] = err_q[X] * d;
    delta[7] = err_q[Y] * d;
    delta[8] = err_q[Z] * d;
}

#define IVP_PERTURBATION NMPC_EPS_4RT
#define IVP_PERTURBATION_RECIP (1.0 / NMPC_EPS_4RT)
void _solve_interval_ivp(qpData_t *qp, interval_t *interval,
real_t *control_ref, real_t *state_ref) {
    size_t i, j;
    real_t integrated_state[NMPC_STATE_DIM], new_state[NMPC_STATE_DIM];
    real_t jacobian[NMPC_DELTA_DIM * NMPC_GRADIENT_DIM]; /* 720B */

    /* Solve the initial value problem at this horizon step. */
    integrated_state = integrator.integrate(
        State(state_ref),
        control_ref,
        dynamics,
        OCP_STEP_LENGTH);

    for(i = 0; i < NMPC_GRADIENT_DIM; i++) {
        real_t perturbed_reference[NMPC_REFERENCE_DIM];
        real_t perturbation = IVP_PERTURBATION,
               perturbation_recip = IVP_PERTURBATION_RECIP;

        #pragma MUST_ITERATE(NMPC_STATE_DIM, NMPC_STATE_DIM)
        for (j = 0; j < NMPC_STATE_DIM; j++) {
            perturbed_reference[j] = state_ref[j];
        }

        perturbed_reference[NMPC_STATE_DIM] = control_ref[0];
        perturbed_reference[NMPC_STATE_DIM + 1u] = control_ref[1];
        perturbed_reference[NMPC_STATE_DIM + 2u] = control_ref[2];

        /* Need to calculate quaternion perturbations using MRPs. */
        if (i < 6u) {
            perturbed_state[i] += perturbation;
        } else if (i >= 6u && i <= 8u) {
            Vector3r d_p;
            d_p << 0.0, 0.0, 0.0;
            d_p[i - 6u] = perturbation;
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
        } else if (i < NMPC_DELTA_DIM) {
            perturbed_state[i + 1u] += perturbation;
        } else {
            /*
            Perturbations for the control inputs should be proportional
            to the control range to make sure we don't lose too much
            precision.
            */
            perturbation *=
                (upper_control_bound[i - NMPC_DELTA_DIM] -
                lower_control_bound[i - NMPC_DELTA_DIM]);
            perturbation_recip = 1.0 / perturbation;
            perturbed_state[i + 1u] += perturbation;
        }

        new_state = integrator.integrate(
            State(perturbed_reference),
            &perturbed_referece[NMPC_STATE_DIM],
            dynamics,
            OCP_STEP_LENGTH);

        /*
        Calculate delta between perturbed state and original state, to
        yield a full column of the Jacobian matrix.
        */
        jacobian.col(i) =
            _state_to_delta(integrated_state, new_state) * perturbation_recip;
    }
}

/*
This step is the first part of the feedback step; the very latest sensor
measurement should be provided in order to to set up the initial state for
the SQP iteration. This allows the feedback delay to be significantly less
than one time step.
*/
void _initial_constraint(StateVector measurement) {
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
    DeltaVector initial_delta = _state_to_delta(
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
void _solve_qp() {
    uint32_t i;
    real_t solution[NMPC_GRADIENT_DIM*(OCP_HORIZON_LENGTH+1)];

    return_t status_flag = qpDUNES_solve(&qp_data);
    AssertSolutionFound(status_flag);

    /* Get the solution. */
    qpDUNES_getPrimalSol(&qp_data, solution);

    /* Get the first control element */
    Eigen::Map<GradientVector> solution_map(&solution[0]);
    control_horizon[0] =
        control_reference[0] +
        solution_map.segment<NMPC_CONTROL_DIM>(NMPC_DELTA_DIM);
}


void nmpc_init() {
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

    for(i = 0; i < OCP_HORIZON_LENGTH; i++) {
        /* Copy the relevant data into the qpDUNES arrays. */
        Q_map = state_weights;
        R_map = control_weights;
        C_map = jacobians[i];
        zLow_map.segment<NMPC_CONTROL_DIM>(NMPC_DELTA_DIM) =
            lower_control_bound - control_reference[i];
        zUpp_map.segment<NMPC_CONTROL_DIM>(NMPC_DELTA_DIM) =
            upper_control_bound - control_reference[i];

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

    qpDUNES_setupAllLocalQPs(&qp_data, QPDUNES_TRUE);

    qpDUNES_indicateDataChange(&qp_data);
}

void nmpc_preparation_step() {
    uint32_t i;
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

    for(i = 0; i < OCP_HORIZON_LENGTH; i++) {
        /* Copy the relevant data into the qpDUNES arrays. */
        C_map = jacobians[i];
        zLow_map.segment<NMPC_CONTROL_DIM>(NMPC_DELTA_DIM) =
            lower_control_bound - control_reference[i];
        zUpp_map.segment<NMPC_CONTROL_DIM>(NMPC_DELTA_DIM) =
            upper_control_bound - control_reference[i];

        status_flag = qpDUNES_updateIntervalData(
            &qp_data, qp_data.intervals[i],
            0, g, C, c, zLow, zUpp, 0, 0, 0, 0);
        AssertOK(status_flag);
    }

    /* Set up final interval. */
    status_flag = qpDUNES_updateIntervalData(&qp_data, qp_data.intervals[i],
        0, g, 0, 0, zLow, zUpp, 0, 0, 0, 0);
    AssertOK(status_flag);

    qpDUNES_indicateDataChange(&qp_data);
}

void nmpc_feedback_step(const real_t measurement[NMPC_STATE_DIM]) {
    initial_constraint(measurement);
    solve_qp();
}

void nmpc_get_controls(real_t controls[NMPC_CONTROL_DIM]) {
    assert(controls);

    /* Return the next control state */
    memcpy(controls, &ocp_control_horizon[0],
           sizeof(real_t) * NMPC_CONTROL_DIM);
}

void nmpc_update_horizon(real_t new_reference[NMPC_REFERENCE_DIM]) {
    uint32_t i;

    /* TODO: Do this more intelligently. */
    for(i = 0; i < OCP_HORIZON_LENGTH-1; i++) {
        state_reference[i] = state_reference[i+1];
        control_reference[i] = control_reference[i+1];
    }

    state_reference[i] = state_reference[i+1];

    /* Insert new reference point. */
    state_reference[i+1] = new_reference.segment<NMPC_STATE_DIM>(0);
    control_reference[i] = new_reference.segment<NMPC_CONTROL_DIM>(
        NMPC_STATE_DIM);

    /* Prepare the QP for the next solution. */
    qpDUNES_shiftLambda(&qp_data);
    qpDUNES_shiftIntervals(&qp_data);

    /* Solve the IVP for the new reference point */
    _solve_interval_ivp(&qp_data, qp_data.intervals[i], state_reference,
                        control_reference);
}

void nmpc_set_state_weights(const real_t coeffs[NMPC_DELTA_DIM]) {
    assert(coeffs);

    /*
    We only store the diagonal of the weight matrices, so they're effectively
    vectors.
    */
    memcpy(ocp_state_weights, coeffs, sizeof(ocp_state_weights));
}

void nmpc_set_control_weights(const real_t coeffs[NMPC_CONTROL_DIM]) {
    assert(coeffs);

    memcpy(ocp_control_weights, coeffs, sizeof(ocp_control_weights));
}

void nmpc_set_terminal_weights(const real_t coeffs[NMPC_DELTA_DIM]) {
    assert(coeffs);

    memcpy(ocp_terminal_weights, coeffs, sizeof(ocp_terminal_weights));
}

void nmpc_set_lower_control_bound(const real_t coeffs[NMPC_CONTROL_DIM]) {
    assert(coeffs);

    memcpy(ocp_lower_control_bound, coeffs, sizeof(ocp_lower_control_bound));
}

void nmpc_set_upper_control_bound(const real_t coeffs[NMPC_CONTROL_DIM]) {
    assert(coeffs);

    memcpy(ocp_upper_control_bound, coeffs, sizeof(ocp_upper_control_bound));
}

void nmpc_set_reference_point(const real_t coeffs[NMPC_REFERENCE_DIM],
uint32_t i) {
    assert(coeffs);
    assert(i <= OCP_HORIZON_LENGTH);

    memcpy(&ocp_state_reference[i * sizeof(real_t) * NMPC_STATE_DIM],
           coeffs, sizeof(real_t) * NMPC_STATE_DIM);

    /*
    Only set control and solve IVPs for regular points, not the final one
    */
    if (i < OCP_HORIZON_LENGTH) {
        memcpy(&ocp_control_reference[i * sizeof(real_t) * NMPC_CONTROL_DIM],
               &coeffs[NMPC_STATE_DIM], sizeof(real_t) * NMPC_CONTROL_DIM);

        /* Solve the IVP for the new reference point */
        _solve_interval_ivp(&qp_data, qp_data.intervals[i], state_reference,
                            control_reference);
    }
}

void nmpc_set_wind_velocity(real_t x, real_t y, real_t z) {
    wind_velocity[0] = x;
    wind_velocity[1] = y;
    wind_velocity[2] = z;
}

uint32_t nmpc_config_get_state_dim() {
    return NMPC_STATE_DIM;
}

uint32_t nmpc_config_get_control_dim() {
    return NMPC_CONTROL_DIM;
}

uint32_t nmpc_config_get_horizon_length() {
    return OCP_HORIZON_LENGTH;
}

real_t nmpc_config_get_step_length() {
    return OCP_STEP_LENGTH;
}

enum nmpc_precision_t nmpc_config_get_precision() {
#ifdef NMPC_SINGLE_PRECISION
    return NMPC_PRECISION_FLOAT;
#else
    return NMPC_PRECISION_DOUBLE;
#endif
}
