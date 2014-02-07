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

#define NMPC_SINGLE_PRECISION

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>

#include "qpDUNES.h"

#include <stdio.h>

void _print_matrix(const char *label, real_t mat[], size_t rows,
size_t cols) {
    printf("%s", label);
    for (size_t i = 0; i < cols; i++) {
        for (size_t j = 0; j < rows; j++) {
            printf("%12.6g ", mat[j*cols + i]);
        }
        printf("\n");
    }
}

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

#ifndef absval
#define absval(x) ((x) < 0 ? -x : x)
#endif

#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif

/* Non-TI compatibility */
#ifndef __TI_COMPILER_VERSION__
#define _nassert(x)
#endif

#ifndef M_PI
#define M_PI ((real_t)3.14159265358979323846)
#define M_PI_2 (M_PI * 0.5)
#define M_PI_4 (M_PI * 0.25)
#endif

#if defined(NMPC_SINGLE_PRECISION)
#define NMPC_EPS_4RT ((real_t)1.857e-2)
#elif defined(NMPC_DOUBLE_PRECISION)
#define NMPC_EPS_4RT ((real_t)1.221e-4)
#endif

#define NMPC_INFTY ((real_t)1.0e12)

#define G_ACCEL ((real_t)9.80665)
#define RHO ((real_t)1.225)

#define sqrt_inv(x) (real_t)(1.0 / sqrt((x)))
#define divide(a, b) ((a) / (b))
#define recip(a) (real_t)(1.0 / (a))
#define fsqrt(a) (real_t)sqrt((a))

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

static inline void quaternion_vector3_multiply(real_t *restrict res,
const real_t *restrict q, const real_t *restrict v) {
    /*
    Multiply a quaternion by a vector (i.e. transform a vectory by a
    quaternion)

    v' = q * v * conjugate(q), or:
    t = 2 * cross(q.xyz, v)
    v' = v + q.w * t + cross(q.xyz, t)

    http://molecularmusings.wordpress.com/2013/05/24/a-faster-quaternion-vector-multiplication/
    */

    assert(res && q && v && res != v && res != q);
    _nassert((size_t)res % 4 == 0);
    _nassert((size_t)q % 4 == 0);
    _nassert((size_t)v % 4 == 0);

    register real_t rx, ry, rz, tx, ty, tz;

    tx = q[Y]*v[Z];
    ty = q[Z]*v[X];
    tx -= q[Z]*v[Y];
    ty -= q[X]*v[Z];
    tz = q[X]*v[Y];
    ty *= 2.0;
    tz -= q[Y]*v[X];
    tx *= 2.0;
    tz *= 2.0;

    rx = v[X];
    rx += q[W]*tx;
    rx += q[Y]*tz;
    rx -= q[Z]*ty;
    res[X] = rx;

    ry = v[Y];
    ry += q[W]*ty;
    ry += q[Z]*tx;
    ry -= q[X]*tz;
    res[Y] = ry;

    rz = v[Z];
    rz += q[W]*tz;
    rz -= q[Y]*tx;
    rz += q[X]*ty;
    res[Z] = rz;
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

static inline void state_scale_add(
real_t *restrict res, const real_t *restrict s1, const real_t a,
const real_t *restrict s2) {
    assert(res && s1 && s2 && s1 != s2 && res != s1 && res != s2);
    _nassert((size_t)res % 4 == 0);
    _nassert((size_t)s1 % 4 == 0);
    _nassert((size_t)s2 % 4 == 0);

    size_t i;
    #pragma MUST_ITERATE(NMPC_STATE_DIM, NMPC_STATE_DIM)
    for (i = 0; i < NMPC_STATE_DIM; i++) {
        res[i] = s2[i] + s1[i] * a;
    }
}

static inline float fatan2(float y, float x) {
#define PI_FLOAT     (real_t)3.14159265
#define PIBY2_FLOAT  (real_t)1.5707963
// |error| < 0.005
    if (x == (real_t)0.0) {
        if (y > (real_t)0.0) return PIBY2_FLOAT;
        if (y == (real_t)0.0) return (real_t)0.0;
        return -PIBY2_FLOAT;
    }
    float res;
    float z = divide(y, x);
    if (absval(z) < (real_t)1.0) {
        res = divide(z, ((real_t)1.0 + (real_t)0.28 * z * z));
        if (x < (real_t)0.0) {
            if (y < (real_t)0.0) return res - PI_FLOAT;
            return res + PI_FLOAT;
        }
    } else {
        res = PIBY2_FLOAT - z / (z * z + (real_t)0.28);
        if (y < (real_t)0.0) return res - PI_FLOAT;
    }
    return res;
#undef PI_FLOAT
#undef PIBY2_FLOAT
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

static void _state_model(real_t *restrict out, const real_t *restrict state,
const real_t *restrict control);
static void _state_integrate_rk4(real_t *restrict out,
const real_t *restrict state, const real_t *restrict control,
const real_t delta);
static void _state_x8_dynamics(real_t *restrict out,
const real_t *restrict state, const real_t *restrict control);
static void _state_to_delta(real_t *delta, const real_t *restrict s1,
const real_t *restrict s2);
static void _solve_interval_ivp(const real_t *restrict state_ref,
const real_t *restrict control_ref, real_t *out_jacobian);
static void _initial_constraint(const real_t measurement[NMPC_STATE_DIM]);
static void _solve_qp(void);


static void _state_model(real_t *restrict out, const real_t *restrict state,
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
    a[X] = state[6 + X];
    a[Y] = state[6 + Y];
    a[Z] = state[6 + Z];
    a[W] = -state[6 + W];
    quaternion_vector3_multiply(&out[3], a, accel);

    /*
    Change in attitude (XYZW): delta_att = 0.5 * (omega_v.conj() * att)
    */
    a[W] = -a[W];
    real_t omega_q_conj[4] = {
        -state[10 + X],
        -state[10 + Y],
        -state[10 + Z],
        0
    };
    quaternion_multiply(&out[6], omega_q_conj, a);
    out[6 + X] *= 0.5;
    out[6 + Y] *= 0.5;
    out[6 + Z] *= 0.5;
    out[6 + W] *= 0.5;

    /* Change in angular velocity */
    out[10] = accel[3];
    out[11] = accel[4];
    out[12] = accel[5];
}

static void _state_integrate_rk4(real_t *restrict out,
const real_t *restrict state, const real_t *restrict control,
const real_t delta) {
    assert(out && state && control);
    _nassert((size_t)out % 4 == 0);
    _nassert((size_t)state % 4 == 0);
    _nassert((size_t)control % 4 == 0);

    /* See include/integrator.h */
    real_t a[NMPC_STATE_DIM], b[NMPC_STATE_DIM], c[NMPC_STATE_DIM],
           d[NMPC_STATE_DIM], temp[NMPC_STATE_DIM];

    /* a = in.model() */
    _state_model(a, state, control);

    /* b = (in + 0.5 * delta * a).model() */
    state_scale_add(temp, a, delta * 0.5f, state);
    _state_model(b, temp, control);

    /* c = (in + 0.5 * delta * b).model() */
    state_scale_add(temp, b, delta * 0.5f, state);
    _state_model(c, temp, control);

    /* d = (in + delta * c).model */
    state_scale_add(temp, c, delta, state);
    _state_model(d, temp, control);

    /* in = in + (delta / 6.0) * (a + (b * 2.0) + (c * 2.0) + d) */
    real_t delta_on_3 = delta * (1.0f/3.0f), delta_on_6 = delta * (1.0f/6.0f);
    size_t i;
    #pragma MUST_ITERATE(NMPC_STATE_DIM, NMPC_STATE_DIM)
    for (i = 0; i < NMPC_STATE_DIM; i++) {
        out[i] = state[i] +
                 delta_on_3 * (b[i] + c[i]) + delta_on_6 * (a[i] + d[i]);
    }
}

static void _state_x8_dynamics(real_t *restrict out,
const real_t *restrict state, const real_t *restrict control) {
    assert(out && state && control);
    _nassert((size_t)out % 4 == 0);
    _nassert((size_t)state % 4 == 0);
    _nassert((size_t)control % 4 == 0);

    /* Work out airflow in NED, then transform to body frame */
    real_t ned_airflow[3], airflow[3];

    ned_airflow[X] = wind_velocity[X] - state[3];
    ned_airflow[Y] = wind_velocity[Y] - state[4];
    ned_airflow[Z] = wind_velocity[Z] - state[5];
    quaternion_vector3_multiply(airflow, &state[6], ned_airflow);

    /*
    Rotate G_ACCEL by current attitude, and set acceleration to that initially

    Extracted out of _mul_quat_vec3 for g = {0, 0, G_ACCEL}
    */
    real_t rx = 0, ry = 0, rz = 0, tx, ty;
    tx = state[6 + Y] * (G_ACCEL * 2.0f);
    ty = -state[6 + X] * (G_ACCEL * 2.0f);
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
    real_t thrust, ve2 = (0.0025f * 0.0025f) * control[0] * control[0];
    /* 1 / 3.8kg times area * density of air */
    thrust = max(0.0f, ve2 - airflow_x2) *
             (0.26315789473684f * 0.5f * RHO * 0.025f);

    /*
    Calculate airflow in the horizontal and vertical planes, as well as
    pressure
    */
    real_t v_inv, horizontal_v2, vertical_v, vertical_v_inv, qbar;

    horizontal_v2 = airflow_x2 + airflow_y2;
    qbar = (RHO * 0.5f) * horizontal_v2;
    v_inv = sqrt_inv(max(1.0f, horizontal_v2 + airflow_z2));

    vertical_v = fsqrt(airflow_x2 + airflow_z2);
    vertical_v_inv = recip(max(1.0f, vertical_v));

    /* Work out sin/cos of alpha and beta */
    real_t alpha, sin_alpha, cos_alpha, sin_beta, cos_beta, a2, sin_cos_alpha;

    sin_beta = airflow[Y] * v_inv;
    cos_beta = vertical_v * v_inv;

    alpha = (float)atan2(-airflow[Z], -airflow[X]);
    a2 = alpha * alpha;

    sin_alpha = -airflow[Z] * vertical_v_inv;
    cos_alpha = -airflow[X] * vertical_v_inv;

    /* Work out aerodynamic forces in wind frame */
    real_t lift, alt_lift, drag, side_force;

    lift = (-5.0f * alpha + 1.0f) * a2 + 2.0f * alpha + 0.3f;
    /* Generalize lift force for very high / very low alpha */
    sin_cos_alpha = sin_alpha * cos_alpha;
    alt_lift = 0.8f * sin_cos_alpha;
    if ((alpha < -0.25f && lift > alt_lift) ||
        (alpha > 0.0f && lift < alt_lift)) {
        lift = alt_lift;
    }

    /* 0.26315789473684 is the reciprocal of mass (3.8kg) */
    lift = (qbar * 0.26315789473684f) * lift;
    drag = (qbar * 0.26315789473684f) *
           (0.05f + 0.8f * sin_alpha * sin_alpha);
    side_force = (qbar * 0.26315789473684f) * 0.3f * sin_beta * cos_beta;

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
           yaw_rate = state[10 + Z],
           pitch_rate = state[10 + Y],
           roll_rate = state[10 + X],
           left_aileron = control[1], right_aileron = control[2];
    pitch_moment = 0.01f - 0.03f * sin_cos_alpha - 0.002f * pitch_rate -
                   0.3f * (left_aileron + right_aileron);
    roll_moment = -0.03f * sin_beta - 0.01f * roll_rate +
                  0.4f * (left_aileron - right_aileron);
    yaw_moment = -0.02f * sin_beta - 0.05f * yaw_rate -
                 0.01f * (absval(left_aileron) + absval(right_aileron));
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
    out[3 + Y] = 5.8823528f * pitch_moment;
    out[3 + X] = (3.364222f * roll_moment + 0.27744448f * yaw_moment);
    out[3 + Z] = (0.27744448f * roll_moment + 2.4920163f * yaw_moment);
}

static void _state_to_delta(real_t *delta, const real_t *restrict s1,
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
#define IVP_PERTURBATION_RECIP (real_t)(1.0 / NMPC_EPS_4RT)
static void _solve_interval_ivp(const real_t *restrict state_ref,
const real_t *restrict control_ref, real_t *out_jacobian) {
    size_t i, j;
    real_t integrated_state[NMPC_STATE_DIM], new_state[NMPC_STATE_DIM];

    /* Solve the initial value problem at this horizon step. */
    _state_integrate_rk4(integrated_state, state_ref, control_ref,
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
            perturbed_reference[i] += IVP_PERTURBATION;
        } else if (i >= 6u && i <= 8u) {
            real_t d_p[3] = { 0.0, 0.0, 0.0 }, delta_q[4], temp[4];
            d_p[i - 6u] = IVP_PERTURBATION;
            /* x_2 = squared norm of d_p */
            delta_q[W] = -(IVP_PERTURBATION * IVP_PERTURBATION) +
                         (real_t)16.0 /
                         (real_t)(16.0 + IVP_PERTURBATION * IVP_PERTURBATION);
            delta_q[X] = ((real_t)1.0 / NMPC_MRP_F) *
                         (NMPC_MRP_A + delta_q[W]) * d_p[X];
            delta_q[Y] = ((real_t)1.0 / NMPC_MRP_F) *
                         (NMPC_MRP_A + delta_q[W]) * d_p[Y];
            delta_q[Z] = ((real_t)1.0 / NMPC_MRP_F) *
                         (NMPC_MRP_A + delta_q[W]) * d_p[Z];
            quaternion_multiply(temp, delta_q, &perturbed_reference[6]);
            perturbed_reference[6] = temp[0];
            perturbed_reference[7] = temp[1];
            perturbed_reference[8] = temp[2];
            perturbed_reference[9] = temp[3];
        } else if (i < NMPC_DELTA_DIM) {
            perturbed_reference[i + 1u] += IVP_PERTURBATION;
        } else {
            /*
            Perturbations for the control inputs should be proportional
            to the control range to make sure we don't lose too much
            precision.
            */
            perturbation *=
                (ocp_upper_control_bound[i - NMPC_DELTA_DIM] -
                ocp_lower_control_bound[i - NMPC_DELTA_DIM]);
            perturbation_recip = (real_t)1.0 / perturbation;
            perturbed_reference[i + 1u] += perturbation;
        }

        _state_integrate_rk4(new_state, perturbed_reference,
                             &perturbed_reference[NMPC_STATE_DIM],
                             OCP_STEP_LENGTH);

        /*
        Calculate delta between perturbed state and original state, to
        yield a full column of the Jacobian matrix.
        */
        real_t jacobian_col[NMPC_GRADIENT_DIM];
        _state_to_delta(jacobian_col, integrated_state, new_state);
        #pragma MUST_ITERATE(NMPC_GRADIENT_DIM, NMPC_GRADIENT_DIM)
        for (j = 0; j < NMPC_GRADIENT_DIM; j++) {
            out_jacobian[NMPC_DELTA_DIM * i + j] = jacobian_col[j] *
                                                   perturbation_recip;
        }
    }
}

/*
This step is the first part of the feedback step; the very latest sensor
measurement should be provided in order to to set up the initial state for
the SQP iteration. This allows the feedback delay to be significantly less
than one time step.
*/
static void _initial_constraint(const real_t measurement[NMPC_STATE_DIM]) {
    real_t z_low[NMPC_GRADIENT_DIM], z_upp[NMPC_GRADIENT_DIM];
    size_t i;

    /*
    Initial delta is constrained to be the difference between the measurement
    and the initial state horizon point.
    */
    _state_to_delta(z_low, ocp_state_reference, measurement);
    memcpy(z_upp, z_low, sizeof(real_t) * NMPC_DELTA_DIM);

    /* Control constraints are unchanged. */
    #pragma MUST_ITERATE(NMPC_CONTROL_DIM, NMPC_CONTROL_DIM)
    for (i = 0; i < NMPC_CONTROL_DIM; i++) {
        z_low[NMPC_DELTA_DIM + i] = ocp_lower_control_bound[i] -
                                    ocp_control_reference[i];
        z_upp[NMPC_DELTA_DIM + i] = ocp_upper_control_bound[i] -
                                    ocp_control_reference[i];
    }

    return_t status_flag = qpDUNES_updateIntervalData(
        &qp_data, qp_data.intervals[0], 0, 0, 0, 0, z_low, z_upp, 0, 0, 0, 0);
    assert(status_flag == QPDUNES_OK);

    qpDUNES_indicateDataChange(&qp_data);
}

/* Solves the QP using qpDUNES. */
static void _solve_qp(void) {
    return_t status_flag;

    status_flag = qpDUNES_solve(&qp_data);
    if (status_flag == QPDUNES_SUCC_OPTIMAL_SOLUTION_FOUND) {
        size_t i;
        real_t solution[NMPC_GRADIENT_DIM * (OCP_HORIZON_LENGTH + 1u)];

        /* Get the solution. */
        qpDUNES_getPrimalSol(&qp_data, solution);

        /* Get the first set of control values */
        for (i = 0; i < NMPC_CONTROL_DIM; i++) {
            ocp_control_value[i] = ocp_control_reference[i] +
                                   solution[NMPC_DELTA_DIM + i];
        }
    } else {
        /* Flag an error in some appropriate way */
    }
}


void nmpc_init(void) {
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
        ocp_lower_state_bound[i] = -NMPC_INFTY;
        ocp_upper_state_bound[i] = NMPC_INFTY;
    }

    /* qpDUNES configuration */
    qp_options = qpDUNES_setupDefaultOptions();
    qp_options.maxIter = 5;
    qp_options.printLevel = 10;
    qp_options.stationarityTolerance = 1e-3f;

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

    /* Set Jacobian to 1 for now */
    memset(jacobian, 0, sizeof(jacobian));

    /* Global state and control constraints */
    memcpy(z_low, ocp_lower_state_bound, sizeof(real_t) * NMPC_DELTA_DIM);
    memcpy(&z_low[NMPC_DELTA_DIM], ocp_lower_control_bound,
           sizeof(real_t) * NMPC_CONTROL_DIM);
    memcpy(z_upp, ocp_upper_state_bound, sizeof(real_t) * NMPC_DELTA_DIM);
    memcpy(&z_upp[NMPC_DELTA_DIM], ocp_upper_control_bound,
           sizeof(real_t) * NMPC_CONTROL_DIM);

    for (i = 0; i < OCP_HORIZON_LENGTH; i++) {
        /* Copy the relevant data into the qpDUNES arrays. */
        status_flag = qpDUNES_setupRegularInterval(
            &qp_data, qp_data.intervals[i],
            0, state_weight_mat, control_weight_mat, 0, gradient, jacobian,
            0, 0, c, z_low, z_upp, 0, 0, 0, 0, 0, 0, 0);
        assert(status_flag == QPDUNES_OK);
    }

    /* Set up final interval. */
    for (i = 0; i < NMPC_DELTA_DIM; i++) {
        state_weight_mat[NMPC_DELTA_DIM * i + i] = ocp_terminal_weights[i];
    }

    status_flag = qpDUNES_setupFinalInterval(&qp_data, qp_data.intervals[i],
        state_weight_mat, gradient, z_low, z_upp, 0, 0, 0);
    assert(status_flag == QPDUNES_OK);


    qpDUNES_setupAllLocalQPs(&qp_data, QPDUNES_TRUE);

    qpDUNES_indicateDataChange(&qp_data);
}

void nmpc_preparation_step(void) {
}

void nmpc_feedback_step(real_t measurement[NMPC_STATE_DIM]) {
    _initial_constraint(measurement);
    _solve_qp();
}

void nmpc_get_controls(real_t controls[NMPC_CONTROL_DIM]) {
    assert(controls);

    /* Return the next control state */
    memcpy(controls, ocp_control_value, sizeof(real_t) * NMPC_CONTROL_DIM);
}

void nmpc_update_horizon(real_t new_reference[NMPC_REFERENCE_DIM]) {
    /*
    Shift reference state and control -- we need to track all these values
    so we can calculate the appropriate delta in _initial_constraint
    */
    memmove(ocp_state_reference, &ocp_state_reference[NMPC_STATE_DIM],
            sizeof(real_t) * NMPC_STATE_DIM * (OCP_HORIZON_LENGTH - 1u));
    memmove(ocp_control_reference, &ocp_control_reference[NMPC_CONTROL_DIM],
            sizeof(real_t) * NMPC_CONTROL_DIM * (OCP_HORIZON_LENGTH - 1u));

    /* Prepare the QP for the next solution. */
    qpDUNES_shiftLambda(&qp_data);
    qpDUNES_shiftIntervals(&qp_data);

    nmpc_set_reference_point(new_reference, OCP_HORIZON_LENGTH - 1);
}

void nmpc_set_state_weights(real_t coeffs[NMPC_DELTA_DIM]) {
    assert(coeffs);

    /*
    We only store the diagonal of the weight matrices, so they're effectively
    vectors.
    */
    memcpy(ocp_state_weights, coeffs, sizeof(ocp_state_weights));
}

void nmpc_set_control_weights(real_t coeffs[NMPC_CONTROL_DIM]) {
    assert(coeffs);

    memcpy(ocp_control_weights, coeffs, sizeof(ocp_control_weights));
}

void nmpc_set_terminal_weights(real_t coeffs[NMPC_DELTA_DIM]) {
    assert(coeffs);

    memcpy(ocp_terminal_weights, coeffs, sizeof(ocp_terminal_weights));
}

void nmpc_set_lower_control_bound(real_t coeffs[NMPC_CONTROL_DIM]) {
    assert(coeffs);

    memcpy(ocp_lower_control_bound, coeffs, sizeof(ocp_lower_control_bound));
}

void nmpc_set_upper_control_bound(real_t coeffs[NMPC_CONTROL_DIM]) {
    assert(coeffs);

    memcpy(ocp_upper_control_bound, coeffs, sizeof(ocp_upper_control_bound));
}

void nmpc_set_reference_point(real_t coeffs[NMPC_REFERENCE_DIM],
uint32_t i) {
    assert(coeffs);
    assert(i <= OCP_HORIZON_LENGTH);

    memcpy(&ocp_state_reference[i * NMPC_STATE_DIM], coeffs,
           sizeof(real_t) * NMPC_STATE_DIM);

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
        memcpy(&ocp_control_reference[i * NMPC_CONTROL_DIM],
               &coeffs[NMPC_STATE_DIM], sizeof(real_t) * NMPC_CONTROL_DIM);

        /* Zero the gradient */
        memset(gradient, 0, sizeof(gradient));

        /* Update state and control constraints */
        memcpy(z_low, ocp_lower_state_bound, sizeof(real_t) * NMPC_DELTA_DIM);
        memcpy(z_upp, ocp_upper_state_bound, sizeof(real_t) * NMPC_DELTA_DIM);

        #pragma MUST_ITERATE(NMPC_CONTROL_DIM, NMPC_CONTROL_DIM)
        for (j = 0; j < NMPC_CONTROL_DIM; j++) {
            z_low[NMPC_DELTA_DIM + j] = ocp_lower_control_bound[j] -
                                        coeffs[NMPC_STATE_DIM + j];
            z_upp[NMPC_DELTA_DIM + j] = ocp_upper_control_bound[j] -
                                        coeffs[NMPC_STATE_DIM + j];
        }

        /*
        Solve the IVP for the new reference point to get the Jacobian (aka
        continuity constraint matrix, C).
        */
        _solve_interval_ivp(coeffs, &coeffs[NMPC_STATE_DIM], jacobian);

        if (i == 11 || i == 12) {
            _print_matrix("g:\n", gradient, NMPC_GRADIENT_DIM, 1);
            _print_matrix("C:\n", jacobian, NMPC_STATE_DIM - 1, NMPC_GRADIENT_DIM);
            _print_matrix("zLow:\n", z_low, NMPC_GRADIENT_DIM, 1);
            _print_matrix("zUpp:\n", z_upp, NMPC_GRADIENT_DIM, 1);
        }

        /* Copy the relevant data into the qpDUNES arrays. */
        status_flag = qpDUNES_updateIntervalData(
            &qp_data, qp_data.intervals[i],
            0, gradient, jacobian, 0, z_low, z_upp, 0, 0, 0, 0);
        assert(status_flag == QPDUNES_OK);
    }
}

void nmpc_set_wind_velocity(real_t x, real_t y, real_t z) {
    wind_velocity[0] = x;
    wind_velocity[1] = y;
    wind_velocity[2] = z;
}

uint32_t nmpc_config_get_state_dim(void) {
    return NMPC_STATE_DIM;
}

uint32_t nmpc_config_get_control_dim(void) {
    return NMPC_CONTROL_DIM;
}

uint32_t nmpc_config_get_horizon_length(void) {
    return OCP_HORIZON_LENGTH;
}

real_t nmpc_config_get_step_length(void) {
    return OCP_STEP_LENGTH;
}

enum nmpc_precision_t nmpc_config_get_precision(void) {
#ifdef NMPC_SINGLE_PRECISION
    return NMPC_PRECISION_FLOAT;
#else
    return NMPC_PRECISION_DOUBLE;
#endif
}
