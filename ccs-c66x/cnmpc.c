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

/* 26052B */
static real_t ocp_state_reference[(OCP_HORIZON_LENGTH + 1u) * NMPC_STATE_DIM];

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

/* Current control solution */
static real_t ocp_control_value[NMPC_CONTROL_DIM];


void _state_model(real_t *restrict out, const real_t *restrict state,
const real_t *restrict control) {
    assert(out && state && control);
    _nassert((size_t)out % 4 == 0);
    _nassert((size_t)state % 4 == 0);
    _nassert((size_t)control % 4 == 0);

    /* See src/state.cpp */
    real_t accel[6];
    _state_x8_dynamics(accel, state, control);

    /* Change in position */
    out[0] = state[3];
    out[1] = state[4];
    out[2] = state[5];

    /* Change in velocity */
    real_t a[4];
    a[X] = -state[6];
    a[Y] = -state[7];
    a[Z] = -state[8];
    a[W] = state[9];
    _mul_quat_vec3(&out[3], a, accel);

    /*
    Change in attitude (XYZW): delta_att = 0.5 * (omega_v.conj() * att)
    */
    a[X] = -a[X];
    a[Y] = -a[Y];
    a[Z] = -a[Z];
    real_t omega_q_conj[4] = {
        -state[10],
        -state[11],
        -state[12],
        0
    };
    _mul_quat_quat(&out[6], omega_q_conj, a);
    out[6] *= 0.5;
    out[7] *= 0.5;
    out[8] *= 0.5;
    out[9] *= 0.5;

    /* Change in angular velocity */
    out[10] = accel[4];
    out[11] = accel[5];
    out[12] = accel[6];
}

void _state_integrate_rk4(struct ukf_state_t *restrict in,
const real_t delta) {
    assert(in);
    _nassert((size_t)in % 8 == 0);

    /* See include/integrator.h */
    struct ukf_state_t a, b, c, d;

    /* a = in.model() */
    memcpy(&a, in, sizeof(a));
    _ukf_state_model(&a);

    /* b = (in + 0.5 * delta * a).model() */
    _mul_state_scalar_add_state(&b, &a, delta * 0.5, in);
    _ukf_state_model(&b);

    /* c = (in + 0.5 * delta * b).model() */
    _mul_state_scalar_add_state(&c, &b, delta * 0.5, in);
    _ukf_state_model(&c);

    /* d = (in + delta * c).model */
    _mul_state_scalar_add_state(&d, &c, delta, in);
    _ukf_state_model(&d);

    /* in = in + (delta / 6.0) * (a + (b * 2.0) + (c * 2.0) + d) */
    real_t *const restrict aptr = (real_t*)&a,
           *const restrict bptr = (real_t*)&b,
           *const restrict cptr = (real_t*)&c,
           *const restrict dptr = (real_t*)&d,
           *const restrict iptr = (real_t*)in;

    real_t delta_on_3 = delta * (1.0/3.0), delta_on_6 = delta * (1.0/6.0);

    uint32_t i;
    #pragma MUST_ITERATE(6)
    for (i = 0; i < 6; i++) {
        iptr[i] += delta_on_3 * (bptr[i] + cptr[i]);
        iptr[i] += delta_on_6 * (aptr[i] + dptr[i]);

        iptr[i + 9] += delta_on_3 * (bptr[i + 9] + cptr[i + 9]);
        iptr[i + 9] += delta_on_6 * (aptr[i + 9] + dptr[i + 9]);
    }
    iptr[15] += delta_on_3 * (bptr[15] + cptr[15]);
    iptr[15] += delta_on_6 * (aptr[15] + dptr[15]);
}

void _state_x8_dynamics(real_t *restrict out, const real_t *restrict state,
const real_t *restrict control) {
    assert(in && control);
    _nassert((size_t)in % 4 == 0);
    _nassert((size_t)control % 4 == 0);

    /* Work out airflow in NED, then transform to body frame */
    real_t ned_airflow[3], airflow[3];

    ned_airflow[X] = wind_velocity[X] - state[3];
    ned_airflow[Y] = wind_velocity[Y] - state[4];
    ned_airflow[Z] = wind_velocity[Z] - state[5];
    _mul_quat_vec3(airflow, &state[6], ned_airflow);

    /*
    Rotate G_ACCEL by current attitude, and set acceleration to that initially

    Extracted out of _mul_quat_vec3 for g = {0, 0, G_ACCEL}
    */
    real_t rx = 0, ry = 0, rz = 0, tx, ty;
    tx = state[6 + Y] * (G_ACCEL * 2.0);
    ty = -state[6 + X] * (G_ACCEL * 2.0);
    rx = state[6 + W] * tx;
    ry = state[6 + W] * ty;
    ry += state[6 + Z] * tx;
    rx -= state[6 + Z] * ty;
    rz = G_ACCEL;
    rz -= state[6 + Y] * tx;
    rz += state[6 + X] * ty;

    out[0 + X] = rx;
    out[0 + Y] = ry;
    out[0 + Z] = rz;

    /*
    Calculate axial airflow
    */
    real_t airflow_x2, airflow_y2, airflow_z2;
    airflow_x2 = airflow[X]*airflow[X];
    airflow_y2 = airflow[Y]*airflow[Y];
    airflow_z2 = airflow[Z]*airflow[Z];

    /*
    Determine airflow magnitude, and the magnitudes of the components in
    the vertical and horizontal planes
    */
    real_t thrust, ve2 = (0.0025 * 0.0025) * control[0] * control[0];
    /* 1 / 3.8kg times area * density of air */
    thrust = max(0.0, ve2 - airflow_x2) * (0.26315789473684 * 0.5 * RHO * 0.025);

    /*
    Calculate airflow in the horizontal and vertical planes, as well as
    pressure
    */
    real_t v_inv, horizontal_v2, vertical_v, vertical_v_inv, qbar;

    horizontal_v2 = airflow_x2 + airflow_y2;
    qbar = (RHO * 0.5) * horizontal_v2;
    v_inv = sqrt_inv(max(1.0, horizontal_v2 + airflow_z2));

    vertical_v = fsqrt(airflow_x2 + airflow_z2);
    vertical_v_inv = recip(max(1.0, vertical_v));

    /* Work out sin/cos of alpha and beta */
    real_t alpha, sin_alpha, cos_alpha, sin_beta, cos_beta, a2, sin_cos_alpha;

    sin_beta = airflow[Y] * v_inv;
    cos_beta = vertical_v * v_inv;

    alpha = fatan2(-airflow[Z], -airflow[X]);
    a2 = alpha * alpha;

    sin_alpha = -airflow[Z] * vertical_v_inv;
    cos_alpha = -airflow[X] * vertical_v_inv;

    /* Work out aerodynamic forces in wind frame */
    real_t lift, alt_lift, drag, side_force;

    lift = (-5 * alpha + 1) * a2 + 2.5 * alpha + 0.12;
    /* Generalize lift force for very high / very low alpha */
    sin_cos_alpha = sin_alpha * cos_alpha;
    alt_lift = 0.8 * sin_cos_alpha;
    if ((alpha < 0.0 && lift > alt_lift) ||
        (alpha > 0.0 && lift < alt_lift)) {
        lift = alt_lift;
    }

    /* 0.26315789473684 is the reciprocal of mass (3.8kg) */
    lift = (qbar * 0.26315789473684) * lift;
    drag = (qbar * 0.26315789473684) * (0.05 + 0.7 * sin_alpha * sin_alpha);
    side_force = (qbar * 0.26315789473684) * 0.3 * sin_beta * cos_beta;

    /* Convert aerodynamic forces from wind frame to body frame */
    real_t x_aero_f = lift * sin_alpha - drag * cos_alpha -
                             side_force * sin_beta,
           z_aero_f = lift * cos_alpha + drag * sin_alpha,
           y_aero_f = side_force * cos_beta;

    out[0 + Y] += y_aero_f;
    out[0 + X] += x_aero_f + thrust;
    out[0 + Z] -= z_aero_f;

    /* Determine moments */
    real_t pitch_moment, yaw_moment, roll_moment,
           yaw_rate = in->angular_velocity[Z],
           pitch_rate = in->angular_velocity[Y],
           roll_rate = in->angular_velocity[X],
           left_aileron = control[1], right_aileron = control[2];
    pitch_moment = 0.001 - 0.1 * sin_cos_alpha - 0.003 * pitch_rate -
                   0.01 * (left_aileron + right_aileron);
    roll_moment = -0.03 * sin_beta - 0.015 * roll_rate +
                  0.025 * (left_aileron - right_aileron);
    yaw_moment = -0.02 * sin_beta - 0.05 * yaw_rate -
                 0.01 * (absval(left_aileron) + absval(right_aileron));
    pitch_moment *= qbar;
    roll_moment *= qbar;
    yaw_moment *= qbar;

    /*
    Calculate angular acceleration (tau / inertia tensor).
    Inertia tensor is:
        0.3 0 -0.0334
        0 0.17 0
        -0.0334 0 0.405
    So inverse is:
        3.36422 0 0.277444
        0 5.88235 0
        0.277444 0 2.49202
    */
    out[3 + Y] = 5.8823528 * pitch_moment;
    out[3 + X] = (3.364222 * roll_moment + 0.27744448 * yaw_moment);
    out[3 + Z] = (0.27744448 * roll_moment + 2.4920163 * yaw_moment);
}

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
real_t *control_ref, real_t *state_ref, real_t *out_jacobian) {
    size_t i, j;
    real_t integrated_state[NMPC_STATE_DIM], new_state[NMPC_STATE_DIM];

    /* Solve the initial value problem at this horizon step. */
    integrated_state = integrator.integrate(
        State(state_ref),
        control_ref,
        dynamics,
        OCP_STEP_LENGTH);

    for (i = 0; i < NMPC_GRADIENT_DIM; i++) {
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
            perturbed_state[i] += IVP_PERTURBATION;
        } else if (i >= 6u && i <= 8u) {
            real_t d_p[3] = { 0.0, 0.0, 0.0 }, delta_q[4], temp[4];
            d_p[i - 6u] = IVP_PERTURBATION;
            /* x_2 = squared norm of d_p */
            delta_q[W] = -(IVP_PERTURBATION * IVP_PERTURBATION) +
                         16.0 / (16.0 + IVP_PERTURBATION * IVP_PERTURBATION);
            delta_q[X] = ((real_t)1.0 / NMPC_MRP_F) *
                         (NMPC_MRP_A + delta_q[W]) * d_p[X];
            delta_q[Y] = ((real_t)1.0 / NMPC_MRP_F) *
                         (NMPC_MRP_A + delta_q[W]) * d_p[Y];
            delta_q[Z] = ((real_t)1.0 / NMPC_MRP_F) *
                         (NMPC_MRP_A + delta_q[W]) * d_p[Z];
            quaternion_multiply(temp, delta_q, &perturbed_state[6]);
            perturbed_state[6] = temp[0];
            perturbed_state[7] = temp[1];
            perturbed_state[8] = temp[2];
            perturbed_state[9] = temp[3];
        } else if (i < NMPC_DELTA_DIM) {
            perturbed_state[i + 1u] += IVP_PERTURBATION;
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
        _state_to_delta(out_jacobian[NMPC_DELTA_DIM * i], integrated_state,
                        new_state);
        #pragma MUST_ITERATE(NMPC_DELTA_DIM, NMPC_DELTA_DIM)
        for (j = 0; j < NMPC_DELTA_DIM; j++) {
            out_jacobian[NMPC_DELTA_DIM * i + j] *= perturbation_recip;
        }
    }
}

/*
This step is the first part of the feedback step; the very latest sensor
measurement should be provided in order to to set up the initial state for
the SQP iteration. This allows the feedback delay to be significantly less
than one time step.
*/
void _initial_constraint(real_t measurement[NMPC_STATE_DIM]) {
    real_t z_low[NMPC_GRADIENT_DIM], z_upp[NMPC_GRADIENT_DIM];
    size_t i;

    /*
    Initial delta is constrained to be the difference between the measurement
    and the initial state horizon point.
    */
    _state_to_delta(z_low, state_reference, measurement);
    memcpy(z_upp, z_low, sizeof(real_t) * NMPC_DELTA_DIM);

    /* Control constraints are unchanged. */
    #pragma MUST_ITERATE(NMPC_CONTROL_DIM, NMPC_CONTROL_DIM)
    for (i = 0; i < NMPC_CONTROL_DIM; i++) {
        z_low[NMPC_DELTA_DIM + i] = lower_control_bound[i] -
                                    control_reference[i];
        z_upp[NMPC_DELTA_DIM + i] = upper_control_bound[i] -
                                    control_reference[i];
    }

    return_t status_flag = qpDUNES_updateIntervalData(
        &qp_data, qp_data.intervals[0], 0, 0, 0, 0, zLow, zUpp, 0, 0, 0, 0);
    AssertOK(status_flag);

    qpDUNES_indicateDataChange(&qp_data);
}

/* Solves the QP using qpDUNES. */
void _solve_qp() {
    return_t status_flag;

    status_flag = qpDUNES_solve(&qp_data);

    if (status_flag == QPDUNES_SUCC_OPTIMAL_SOLUTION_FOUND) {
        size_t i;
        real_t solution[NMPC_GRADIENT_DIM * (OCP_HORIZON_LENGTH + 1u)];

        /* Get the solution. */
        qpDUNES_getPrimalSol(&qp_data, solution);

        /* Get the first set of control values */
        for (i = 0; i < NMPC_CONTROL_DIM; i++) {
            ocp_control_value[i] = control_reference[i] +
                                   solution[NMPC_DELTA_DIM + i];
        }
    } else {
        /* Flag an error in some appropriate way */
    }
}


void nmpc_init() {
    real_t jacobian[NMPC_DELTA_DIM * NMPC_GRADIENT_DIM], /* 720B */
           z_low[NMPC_GRADIENT_DIM],
           z_upp[NMPC_GRADIENT_DIM],
           gradient[NMPC_GRADIENT_DIM],
           c[NMPC_DELTA_DIM],
           state_weight_mat[NMPC_DELTA_DIM * NMPC_DELTA_DIM], /* 576B */
           control_weight_mat[NMPC_CONTROL_DIM * NMPC_CONTROL_DIM];
    return_t status_flag;
    size_t i;

    /* Initialise state inequality constraints to +/-infinity. */
    for (i = 0; i < NMPC_DELTA_DIM; i++) {
        lower_state_bound[i] = -NMPC_INFTY;
        upper_state_bound[i] = NMPC_INFTY;
    }

    /* qpDUNES configuration */
    qp_options = qpDUNES_setupDefaultOptions();
    qp_options.maxIter = 5;
    qp_options.printLevel = 0;
    qp_options.stationarityTolerance = 1e-3;

    /* Set up problem dimensions. */
    qpDUNES_setup(
        &qp_data,
        OCP_HORIZON_LENGTH,
        NMPC_DELTA_DIM,
        NMPC_CONTROL_DIM,
        0,
        &qp_options);

    /* Convert state and control diagonals into full matrices */
    memset(state_weight_mat, 0, sizeof(state_weight_mat));
    for (i = 0; i < NMPC_DELTA_DIM; i++) {
        state_weight_mat[NMPC_DELTA_DIM * i + i] = ocp_state_weights[i];
    }

    memset(control_weight_mat, 0, sizeof(control_weight_mat));
    for (i = 0; i < NMPC_CONTROL_DIM; i++) {
        control_weight_mat[NMPC_CONTROL_DIM * i + i] = ocp_control_weights[i];
    }

    /* Gradient vector fixed to zero. */
    memset(gradient, 0, sizeof(gradient));

    /* Continuity constraint constant term fixed to zero. */
    memset(c, 0, sizeof(c));

    /* Zero Jacobian for now */
    memset(jacobian, 0, sizeof(jacobian));

    /* Global state and control constraints */
    memcpy(z_low, lower_state_bound, sizeof(real_t) * NMPC_DELTA_DIM);
    memcpy(&z_low[NMPC_DELTA_DIM], lower_control_bound,
           sizeof(real_t) * NMPC_CONTROL_DIM);
    memcpy(z_upp, upper_state_bound, sizeof(real_t) * NMPC_DELTA_DIM);
    memcpy(&z_upp[NMPC_DELTA_DIM], upper_control_bound,
           sizeof(real_t) * NMPC_CONTROL_DIM);

    for (i = 0; i < OCP_HORIZON_LENGTH; i++) {
        /* Copy the relevant data into the qpDUNES arrays. */
        status_flag = qpDUNES_setupRegularInterval(
            &qp_data, qp_data.intervals[i],
            0, state_weights, control_weights, 0, gradient, jacobian, 0, 0, c,
            z_low, z_upp, 0, 0, 0, 0, 0, 0, 0);
        AssertOK(status_flag);
    }

    /* Set up final interval. */
    for (i = 0; i < NMPC_DELTA_DIM; i++) {
        state_weight_mat[NMPC_DELTA_DIM * i + i] = ocp_terminal_weights[i];
    }
    status_flag = qpDUNES_setupFinalInterval(&qp_data, qp_data.intervals[i],
        state_weight_mat, g, z_low, z_upp, 0, 0, 0);
    AssertOK(status_flag);

    qpDUNES_setupAllLocalQPs(&qp_data, QPDUNES_TRUE);

    qpDUNES_indicateDataChange(&qp_data);
}

void nmpc_preparation_step() {
}

void nmpc_feedback_step(const real_t measurement[NMPC_STATE_DIM]) {
    _initial_constraint(measurement);
    _solve_qp();
}

void nmpc_get_controls(real_t controls[NMPC_CONTROL_DIM]) {
    assert(controls);

    /* Return the next control state */
    memcpy(controls, &ocp_control_horizon[0],
           sizeof(real_t) * NMPC_CONTROL_DIM);
}

void nmpc_update_horizon(real_t new_reference[NMPC_REFERENCE_DIM]) {
    /*
    Shift reference state and control -- we need to track all these values
    so we can calculate the appropriate delta in _initial_constraint
    */
    memmove(state_reference, &state_reference[NMPC_STATE_DIM],
            sizeof(real_t) * NMPC_STATE_DIM * (OCP_HORIZON_LENGTH - 1u));
    memmove(control_reference, &control_reference[NMPC_CONTROL_DIM],
            sizeof(real_t) * NMPC_CONTROL_DIM * (OCP_HORIZON_LENGTH - 1u));

    /* Prepare the QP for the next solution. */
    qpDUNES_shiftLambda(&qp_data);
    qpDUNES_shiftIntervals(&qp_data);

    nmpc_set_reference_point(new_reference, OCP_HORIZON_LENGTH - 1);
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
        real_t jacobian[NMPC_DELTA_DIM * NMPC_GRADIENT_DIM], /* 720B */
               z_low[NMPC_GRADIENT_DIM],
               z_upp[NMPC_GRADIENT_DIM],
               gradient[NMPC_GRADIENT_DIM];
        return_t status_flag;
        size_t j;

        /* Copy the control reference */
        memcpy(&ocp_control_reference[i * sizeof(real_t) * NMPC_CONTROL_DIM],
               &coeffs[NMPC_STATE_DIM], sizeof(real_t) * NMPC_CONTROL_DIM);

        /* Zero the gradient */
        memset(gradient, 0, sizeof(gradient));

        /* Update state and control constraints */
        memcpy(z_low, lower_state_bound, sizeof(real_t) * NMPC_DELTA_DIM);
        memcpy(z_upp, upper_state_bound, sizeof(real_t) * NMPC_DELTA_DIM);

        #pragma MUST_ITERATE(NMPC_CONTROL_DIM, NMPC_CONTROL_DIM)
        for (j = 0; j < NMPC_CONTROL_DIM; j++) {
            z_low[NMPC_DELTA_DIM + j] = lower_control_bound[j] -
                                        coeffs[NMPC_STATE_DIM + j];
            z_upp[NMPC_DELTA_DIM + j] = upper_control_bound[j] -
                                        coeffs[NMPC_STATE_DIM + j];
        }

        /*
        Solve the IVP for the new reference point to get the Jacobian (aka
        continuity constraint matrix, C).
        */
        _solve_interval_ivp(&qp_data, qp_data.intervals[i], coeffs,
                            &coeffs[NMPC_STATE_DIM], jacobian);

        /* Copy the relevant data into the qpDUNES arrays. */
        status_flag = qpDUNES_updateIntervalData(
            &qp_data, qp_data.intervals[i],
            0, gradient, jacobian, 0, z_low, z_upp, 0, 0, 0, 0);
        AssertOK(status_flag);
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
