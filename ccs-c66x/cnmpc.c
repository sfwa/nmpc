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

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <stdbool.h>

#ifndef __TI_COMPILER_VERSION__
    #include "config.h"
    #include "../c/cnmpc.h"
#else
    #define USE_DSP_INTRINSICS

    #include "config.h"
    #include "cnmpc.h"
#endif

#include "c66math.h"
#include "qpDUNES/qpDUNES.h"

/*
Use static allocation for qpDUNES structures, since the sizes are all known at
compile time -- see qpDUNES/setup_qp.c:40-267
*/
#define nX NMPC_DELTA_DIM
#define nU NMPC_CONTROL_DIM
#define nI OCP_HORIZON_LENGTH
#define nZ (nX + nU)
#define nV nZ

struct static_interval_t { /* 1584B + sizeof(interval) */
    interval_t interval;

    /*
    Statically-allocated storage for the interval matrices. Since nD is zero,
    some of these are set to 1 to avoid non-standard zero-length arrays.

    H_data and cholH_data are meant to be nV * nV, but since they'll both
    always be diagonal we can just use one row.
    */
    real_t H_data[nV]; /* 60B */
    real_t cholH_data[nV]; /* 60B */
    real_t q_data[nV]; /* 60B */
    real_t C_data[nX * nV]; /* 720B */
    real_t c_data[nX]; /* 48B */
    real_t zLow_data[nV]; /* 60B */
    real_t zUpp_data[nV]; /* 60B */
    real_t z_data[nV]; /* 60B */
    real_t y_data[2u * nV]; /* 120B */
    real_t lambdaK_data[nX]; /* 48B */
    real_t lambdaK1_data[nX]; /* 48B */
    real_t clippingSolver_qStep_data[nV]; /* 60B */
    real_t clippingSolver_zUnconstrained_data[nV]; /* 60B */
    real_t clippingSolver_dz_data[nV]; /* 60B */

    /*
    These are allocated in qpDUNES_setup rather than qpDUNES_allocInterval,
    but they're still per-interval.
    */
    real_t zVecTmp_data[nZ]; /* 60B */
};

struct static_qpdata_t {
    qpData_t qpdata;

    /* QP solver data */
    interval_t *intervals_data[nI + 1u]; /* 404B */
    struct static_interval_t interval_recs[nI + 1u]; /* 287244B + 101*sizeof(interval) */
    real_t lambda_data[nX * nI]; /* 4800B */
    real_t deltaLambda_data[nX * nI]; /* 4800B */
    real_t hessian_data[nX * 2u * nX * nI]; /* 115200B */
    real_t cholHessian_data[nX * 2u * nX * nI]; /* 115200B */
    real_t gradient_data[nX * nI]; /* 4800B */
    real_t xVecTmp_data[nX]; /* 48B */
    real_t xxMatTmp_data[nX * nX]; /* 576B */
    real_t zxMatTmp_data[nZ * nX]; /* 720B */

    /* Logging data */
    itLog_t itLog_data;
    int_t *ieqStatus_data[nI + 1u]; /* 404B */
    int_t *prevIeqStatus_data[nI + 1u]; /* 404B */
    int_t ieqStatus_n_data[(nI + 1u) * nZ]; /* 6060B */
    int_t prevIeqStatus_n_data[(nI + 1u) * nZ]; /* 6060B */
};

#undef nX
#undef nU
#undef nD
#undef nI
#undef nZ
#undef nV

/* Math routines */
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

static real_t wind_velocity[3];

/* 26052B */
static real_t ocp_state_reference[(OCP_HORIZON_LENGTH + 1u) * NMPC_STATE_DIM];

/* 6000B */
static real_t ocp_control_reference[(OCP_HORIZON_LENGTH + 1u) *
                                    NMPC_CONTROL_DIM];

static real_t ocp_lower_state_bound[NMPC_DELTA_DIM];
static real_t ocp_upper_state_bound[NMPC_DELTA_DIM];
static real_t ocp_lower_control_bound[NMPC_CONTROL_DIM];
static real_t ocp_upper_control_bound[NMPC_CONTROL_DIM];
static real_t ocp_state_weights[NMPC_DELTA_DIM]; /* diagonal only */
static real_t ocp_terminal_weights[NMPC_DELTA_DIM]; /* diagonal only */
static real_t ocp_control_weights[NMPC_CONTROL_DIM]; /* diagonal only */

static struct static_qpdata_t ocp_qp_data;
static bool ocp_state_position_is_delta = false;

/* Current control solution */
static real_t ocp_control_value[NMPC_CONTROL_DIM];
static real_t ocp_last_result;


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
const real_t *restrict control_ref, real_t *restrict out_jacobian,
const real_t *restrict next_state_ref, real_t *restrict out_residuals);
static void _initial_constraint(const real_t measurement[NMPC_STATE_DIM]);
static real_t _solve_qp(void);


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
    real_t omega_q_conj[4];
    omega_q_conj[X] = -state[10 + X];
    omega_q_conj[Y] = -state[10 + Y];
    omega_q_conj[Z] = -state[10 + Z];
    omega_q_conj[W] = 0.0;
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
    real_t airflow_x2, airflow_y2, airflow_z2, airflow_v2;
    airflow_x2 = airflow[X]*airflow[X];
    airflow_y2 = airflow[Y]*airflow[Y];
    airflow_z2 = airflow[Z]*airflow[Z];
    airflow_v2 = airflow_x2 + airflow_y2 + airflow_z2;

    /*
    Determine airflow magnitude, and the magnitudes of the components in
    the vertical and horizontal planes
    */
    real_t rpm = (control[0] - 0.15f) * 12000.0f, thrust,
           ve2 = (0.012f * 0.012f) * rpm * rpm;
    /* 1 / 3.8kg times area * density of air */
    thrust = (ve2 - airflow_v2) *
             (0.26315789473684f * 0.5f * RHO * 0.02f);
    if (thrust < 0.0f) {
        thrust = 0.0f;
    }

    /*
    Calculate airflow in the horizontal and vertical planes, as well as
    pressure
    */
    real_t v_inv, vertical_v, vertical_v_inv, qbar;

    qbar = (RHO * 0.5f) * airflow_v2;
    v_inv = recip_sqrt_f(max(1.0f, airflow_v2));

    vertical_v = sqrt_f(airflow_x2 + airflow_z2);
    vertical_v_inv = recip_f(max(1.0f, vertical_v));

    /* Work out sin/cos of alpha and beta */
    real_t sin_alpha, cos_alpha, sin_beta, cos_beta, sin_cos_alpha;

    sin_beta = airflow[Y] * v_inv;
    cos_beta = vertical_v * v_inv;

    sin_alpha = -airflow[Z] * vertical_v_inv;
    cos_alpha = -airflow[X] * vertical_v_inv;

    sin_cos_alpha = sin_alpha * cos_alpha;

    /* Work out aerodynamic forces in wind frame */
    real_t lift, drag, side_force;

    /* 0.26315789473684 is the reciprocal of mass (3.8kg) */
    lift = (qbar * 0.26315789473684f) * (0.7f * sin_cos_alpha + 0.15f);
    drag = (qbar * 0.26315789473684f) *
           (0.05f + 0.7f * sin_alpha * sin_alpha);
    side_force = (qbar * 0.26315789473684f) * 0.05f * sin_beta * cos_beta;

    /* Convert aerodynamic forces from wind frame to body frame */
    real_t x_aero_f = lift * sin_alpha - drag * cos_alpha +
                             side_force * sin_beta,
           z_aero_f = lift * cos_alpha + drag * sin_alpha,
           y_aero_f = side_force * cos_beta;

    out[0 + Y] += y_aero_f;
    out[0 + X] += x_aero_f + thrust;
    out[0 + Z] -= z_aero_f;

    /* Limit the cubic portion of control surface moments to a sane value */
    if (vertical_v < 20.0f) {
        vertical_v = 20.0f;
    }

    /* Determine moments */
    real_t pitch_moment, yaw_moment, roll_moment,
           yaw_rate = state[10 + Z],
           pitch_rate = state[10 + Y],
           roll_rate = state[10 + X],
           left_aileron = control[1] - 0.5f,
           right_aileron = control[2] - 0.5f;
    pitch_moment = 0.0f - 0.0f * sin_alpha - 0.0f * pitch_rate -
                   0.15f * (left_aileron + right_aileron) * vertical_v * 0.1f;
    roll_moment = 0.05f * sin_beta - 0.05f * roll_rate +
                  0.35f * (left_aileron - right_aileron) * vertical_v * 0.1f;
    yaw_moment = 0.025f * sin_beta - 0.025f * yaw_rate -
                 0.15f * (absval(left_aileron) - absval(right_aileron)) *
                 vertical_v * 0.1f;
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
    out[3 + Y] = 15.8823528f * pitch_moment;
    out[3 + X] = (3.864222f * roll_moment + 0.27744448f * yaw_moment);
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

    real_t d = NMPC_MRP_F * recip_f(NMPC_MRP_A + err_q[W]);
    delta[6] = err_q[X] * d;
    delta[7] = err_q[Y] * d;
    delta[8] = err_q[Z] * d;
}

#define IVP_PERTURBATION NMPC_EPS_4RT
#define IVP_PERTURBATION_RECIP (real_t)(1.0 / IVP_PERTURBATION)
static void _solve_interval_ivp(const real_t *restrict state_ref,
const real_t *restrict control_ref, real_t *restrict out_jacobian,
const real_t *restrict next_state_ref, real_t *restrict out_residuals) {
    size_t i, j;
    real_t integrated_state[NMPC_STATE_DIM], new_state[NMPC_STATE_DIM];

    /* Solve the initial value problem at this horizon step. */
    _state_integrate_rk4(integrated_state, state_ref, control_ref,
                         OCP_STEP_LENGTH);

    /*
    Calculate integration residuals -- the difference between the integrated
    state and the next state.
    */
    if (ocp_state_position_is_delta) {
        integrated_state[0] -= state_ref[0];
        integrated_state[1] -= state_ref[1];
        integrated_state[2] -= state_ref[2];
    }
    _state_to_delta(out_residuals, next_state_ref, integrated_state);
    if (ocp_state_position_is_delta) {
        integrated_state[0] += state_ref[0];
        integrated_state[1] += state_ref[1];
        integrated_state[2] += state_ref[2];
    }

    /* Calculate the Jacobian */
    for (i = 0; i < NMPC_GRADIENT_DIM; i++) {
        real_t perturbed_reference[NMPC_REFERENCE_DIM];

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
        } else {
            perturbed_reference[i + 1u] += IVP_PERTURBATION;
        }

        _state_integrate_rk4(new_state, perturbed_reference,
                             &perturbed_reference[NMPC_STATE_DIM],
                             OCP_STEP_LENGTH);

        /*
        Calculate delta between perturbed state and original state, to
        yield a full column of the Jacobian matrix. Transpose during the copy
        to match the qpDUNES row-major convention.
        */
        real_t jacobian_col[NMPC_DELTA_DIM];
        _state_to_delta(jacobian_col, integrated_state, new_state);
        #pragma MUST_ITERATE(NMPC_DELTA_DIM, NMPC_DELTA_DIM)
        for (j = 0; j < NMPC_DELTA_DIM; j++) {
            out_jacobian[NMPC_GRADIENT_DIM * j + i] = jacobian_col[j] *
                                                      IVP_PERTURBATION_RECIP;
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

    return_t status_flag;
    status_flag = qpDUNES_updateIntervalConstraints(
        &ocp_qp_data.qpdata, ocp_qp_data.qpdata.intervals[0], z_low, z_upp);
    assert(status_flag == QPDUNES_OK);
}

/* Solves the QP using qpDUNES. */
static real_t _solve_qp(void) {
    return_t status_flag;
    real_t objective_value;
    size_t i;

    status_flag = qpDUNES_solve(&ocp_qp_data.qpdata);

    if (status_flag == QPDUNES_SUCC_OPTIMAL_SOLUTION_FOUND) {
        /* Check the objective value. */
        objective_value =
            abs_f(qpDUNES_computeObjectiveValue(&ocp_qp_data.qpdata));

        /* Get the first set of control values */
        for (i = 0; i < NMPC_CONTROL_DIM; i++) {
            ocp_control_value[i] = ocp_control_reference[i] +
                ocp_qp_data.qpdata.intervals[0]->z.data[NMPC_DELTA_DIM + i];
        }

        return objective_value;
    } else {
        /* Check the objective value. */
        objective_value =
            abs_f(qpDUNES_computeObjectiveValue(&ocp_qp_data.qpdata));

        /*
        Get the first set of control values, but scale closer to the reference
        because they're probably all over the place.
        */
        for (i = 0; i < NMPC_CONTROL_DIM; i++) {
            ocp_control_value[i] = ocp_control_reference[i] +
                ocp_qp_data.qpdata.intervals[0]->z.data[NMPC_DELTA_DIM + i];
        }

        /* Flag an error in some appropriate way */
        return -1.0;
    }
}

/*
Set up a qpDUNES interval with static allocation -- refer to
qpDUNES_allocInterval at qpDUNES/setup_qp.c:210
*/
static void _init_static_interval(struct static_interval_t *i, size_t nV) {
    assert(i);

    i->interval.nV = (uint32_t)nV;

    i->interval.H.data = i->H_data;
    i->interval.H.sparsityType = QPDUNES_MATRIX_UNDEFINED;
    i->interval.cholH.data = i->cholH_data;
    i->interval.cholH.sparsityType = QPDUNES_MATRIX_UNDEFINED;

    i->interval.q.data = i->q_data;

    i->interval.C.data = i->C_data;
    i->interval.C.sparsityType = QPDUNES_MATRIX_UNDEFINED;
    i->interval.c.data = i->c_data;

    i->interval.zLow.data = i->zLow_data;
    i->interval.zUpp.data = i->zUpp_data;

    i->interval.z.data = i->z_data;

    i->interval.y.data = i->y_data;

    i->interval.lambdaK.data = i->lambdaK_data;
    i->interval.lambdaK.isDefined = QPDUNES_TRUE;

    i->interval.lambdaK1.data = i->lambdaK1_data;
    i->interval.lambdaK1.isDefined = QPDUNES_TRUE;

    i->interval.qpSolverClipping.qStep.data = i->clippingSolver_qStep_data;
    i->interval.qpSolverClipping.zUnconstrained.data =
        i->clippingSolver_zUnconstrained_data;
    i->interval.qpSolverClipping.dz.data = i->clippingSolver_dz_data;

    /*
    Per-interval allocation within qpDUNES_setup, migrated here for
    convenience
    */
    i->interval.zVecTmp.data = i->zVecTmp_data;
}

/*
Set up the qpDunes qpData_t structure with static allocation -- refer to
qpDUNES_setup at qpDUNES/setup_qp.c:40.

Equivalent to:
    qpDUNES_setup(
        &qp->qpdata,
        OCP_HORIZON_LENGTH,
        NMPC_DELTA_DIM,
        NMPC_CONTROL_DIM,
        0,
        opts);
*/
static void _init_static_qp(struct static_qpdata_t *qp,
const qpOptions_t *opts) {
    assert(qp);
    assert(opts);

    size_t i, nV;

    qp->qpdata.options = *opts;
    qp->qpdata.nI = OCP_HORIZON_LENGTH;
    qp->qpdata.nX = NMPC_DELTA_DIM;
    qp->qpdata.nU = NMPC_CONTROL_DIM;
    qp->qpdata.nZ = NMPC_DELTA_DIM + NMPC_CONTROL_DIM;
    qp->qpdata.nDttl = 0;

    qp->qpdata.intervals = qp->intervals_data;

    for (i = 0; i < OCP_HORIZON_LENGTH + 1u; i++) {
        if (i < OCP_HORIZON_LENGTH) {
            nV = NMPC_DELTA_DIM + NMPC_CONTROL_DIM;
        } else {
            nV = NMPC_DELTA_DIM;
        }

        qp->qpdata.intervals[i] = &(qp->interval_recs[i].interval);
        qp->qpdata.intervals[i]->id = (uint32_t)i;

        _init_static_interval(&(qp->interval_recs[i]), nV);
    }

    /* Last interval doesn't need a Jacobian */
    qpDUNES_setMatrixNull(&(qp->qpdata.intervals[OCP_HORIZON_LENGTH]->C));

    qp->qpdata.intervals[0]->lambdaK.isDefined = QPDUNES_FALSE;
    qp->qpdata.intervals[OCP_HORIZON_LENGTH]->lambdaK1.isDefined =
        QPDUNES_FALSE;

    qp->qpdata.lambda.data = qp->lambda_data;
    qp->qpdata.deltaLambda.data = qp->deltaLambda_data;

    qp->qpdata.hessian.data = qp->hessian_data;
    qp->qpdata.cholHessian.data = qp->cholHessian_data;
    qp->qpdata.gradient.data = qp->gradient_data;

    qp->qpdata.xVecTmp.data = qp->xVecTmp_data;
    qp->qpdata.xxMatTmp.data = qp->xxMatTmp_data;
    qp->qpdata.zxMatTmp.data = qp->zxMatTmp_data;

    qp->qpdata.optObjVal = -qp->qpdata.options.QPDUNES_INFTY;

    qp->qpdata.log.itLog = &(qp->itLog_data);
    qp->qpdata.log.itLog[0].ieqStatus = qp->ieqStatus_data;
    qp->qpdata.log.itLog[0].prevIeqStatus = qp->prevIeqStatus_data;
    for (i = 0; i < OCP_HORIZON_LENGTH + 1u; i++) {
        qp->qpdata.log.itLog[0].ieqStatus[i] =
            &(qp->ieqStatus_n_data[i * qp->qpdata.nZ]);
        qp->qpdata.log.itLog[0].prevIeqStatus[i] =
            &(qp->prevIeqStatus_n_data[i * qp->qpdata.nZ]);
    }
}

#pragma FUNC_EXT_CALLED(nmpc_init);
void nmpc_init(bool use_relative_positions) {
    real_t C[NMPC_DELTA_DIM * NMPC_GRADIENT_DIM], /* 720B */
           z_low[NMPC_GRADIENT_DIM],
           z_upp[NMPC_GRADIENT_DIM],
           c[NMPC_DELTA_DIM];
    return_t status_flag;
    qpOptions_t qp_options;
    size_t i;

    memset(&ocp_qp_data, 0, sizeof(ocp_qp_data));

    ocp_state_position_is_delta = use_relative_positions;

    /* Initialise state inequality constraints to +/-infinity. */
    for (i = 0; i < NMPC_DELTA_DIM; i++) {
        ocp_lower_state_bound[i] = -NMPC_INFTY;
        ocp_upper_state_bound[i] = NMPC_INFTY;
    }

    /* qpDUNES configuration */
    qp_options = qpDUNES_setupDefaultOptions();
    qp_options.maxIter = 4;
    qp_options.printLevel = 0;
    qp_options.stationarityTolerance = 1e-1f;

    /* Set up problem dimensions. */
    _init_static_qp(&ocp_qp_data, &qp_options);

    /* Continuity constraint constant term fixed to zero. */
    memset(c, 0, sizeof(c));

    /* Set Jacobian to 0 for now */
    memset(C, 0, sizeof(C));

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
            &ocp_qp_data.qpdata, ocp_qp_data.qpdata.intervals[i],
            ocp_state_weights, ocp_control_weights, C, c, z_low, z_upp);
        assert(status_flag == QPDUNES_OK);
    }

    /* Set up final interval. */
    status_flag = qpDUNES_setupFinalInterval(
        &ocp_qp_data.qpdata, ocp_qp_data.qpdata.intervals[OCP_HORIZON_LENGTH],
        ocp_terminal_weights, z_low, z_upp);
    assert(status_flag == QPDUNES_OK);

    qpDUNES_setupAllLocalQPs(&ocp_qp_data.qpdata);
}

#pragma FUNC_EXT_CALLED(nmpc_preparation_step);
void nmpc_preparation_step(void) {
    qpDUNES_prepare(&ocp_qp_data.qpdata);
}

#pragma FUNC_EXT_CALLED(nmpc_feedback_step);
void nmpc_feedback_step(real_t measurement[NMPC_STATE_DIM]) {
    _initial_constraint(measurement);
    ocp_last_result = _solve_qp();
}

#pragma FUNC_EXT_CALLED(nmpc_get_controls);
enum nmpc_result_t nmpc_get_controls(real_t controls[NMPC_CONTROL_DIM]) {
    assert(controls);

    /* Return the next control state */
    memcpy(controls, ocp_control_value, sizeof(real_t) * NMPC_CONTROL_DIM);

    if (ocp_last_result >= 0.0f) {
        return NMPC_OK;
    } else {
        return NMPC_INFEASIBLE;
    }
}

#pragma FUNC_EXT_CALLED(nmpc_get_objective_value);
real_t nmpc_get_objective_value(void) {
    return ocp_last_result;
}

#pragma FUNC_EXT_CALLED(nmpc_update_horizon);
void nmpc_update_horizon(real_t new_reference[NMPC_REFERENCE_DIM]) {
    /*
    Shift reference state and control -- we need to track all these values
    so we can calculate the appropriate delta in _initial_constraint
    */
    memmove(ocp_state_reference, &ocp_state_reference[NMPC_STATE_DIM],
            sizeof(real_t) * NMPC_STATE_DIM * OCP_HORIZON_LENGTH);
    memmove(ocp_control_reference, &ocp_control_reference[NMPC_CONTROL_DIM],
            sizeof(real_t) * NMPC_CONTROL_DIM * OCP_HORIZON_LENGTH);

    /* Prepare the QP for the next solution. */
    qpDUNES_shiftLambda(&ocp_qp_data.qpdata);
    qpDUNES_shiftIntervals(&ocp_qp_data.qpdata);

    nmpc_set_reference_point(new_reference, OCP_HORIZON_LENGTH);

    size_t i;
    for (i = 0; i < OCP_HORIZON_LENGTH - 1u; i++) {
        qpDUNES_setupStageQP(&ocp_qp_data.qpdata,
                             ocp_qp_data.qpdata.intervals[i]);
    }
}

#pragma FUNC_EXT_CALLED(nmpc_set_state_weights);
void nmpc_set_state_weights(real_t coeffs[NMPC_DELTA_DIM]) {
    assert(coeffs);

    /*
    We only store the diagonal of the weight matrices, so they're effectively
    vectors.
    */
    memcpy(ocp_state_weights, coeffs, sizeof(ocp_state_weights));
}

#pragma FUNC_EXT_CALLED(nmpc_set_control_weights);
void nmpc_set_control_weights(real_t coeffs[NMPC_CONTROL_DIM]) {
    assert(coeffs);

    memcpy(ocp_control_weights, coeffs, sizeof(ocp_control_weights));
}

#pragma FUNC_EXT_CALLED(nmpc_set_terminal_weights);
void nmpc_set_terminal_weights(real_t coeffs[NMPC_DELTA_DIM]) {
    assert(coeffs);

    memcpy(ocp_terminal_weights, coeffs, sizeof(ocp_terminal_weights));
}

#pragma FUNC_EXT_CALLED(nmpc_set_lower_control_bound);
void nmpc_set_lower_control_bound(real_t coeffs[NMPC_CONTROL_DIM]) {
    assert(coeffs);

    memcpy(ocp_lower_control_bound, coeffs, sizeof(ocp_lower_control_bound));
}

#pragma FUNC_EXT_CALLED(nmpc_set_upper_control_bound);
void nmpc_set_upper_control_bound(real_t coeffs[NMPC_CONTROL_DIM]) {
    assert(coeffs);

    memcpy(ocp_upper_control_bound, coeffs, sizeof(ocp_upper_control_bound));
}

#pragma FUNC_EXT_CALLED(nmpc_set_reference_point);
void nmpc_set_reference_point(real_t coeffs[NMPC_REFERENCE_DIM],
uint32_t i) {
    assert(coeffs);
    assert(i <= OCP_HORIZON_LENGTH);

    memcpy(&ocp_state_reference[i * NMPC_STATE_DIM], coeffs,
           sizeof(real_t) * NMPC_STATE_DIM);
    memcpy(&ocp_control_reference[i * NMPC_CONTROL_DIM],
           &coeffs[NMPC_STATE_DIM], sizeof(real_t) * NMPC_CONTROL_DIM);

    /*
    Only set control and solve IVPs for regular points, not the final one
    */
    if (i > 0 && i <= OCP_HORIZON_LENGTH) {
        real_t jacobian[NMPC_DELTA_DIM * NMPC_GRADIENT_DIM], /* 720B */
               z_low[NMPC_GRADIENT_DIM],
               z_upp[NMPC_GRADIENT_DIM],
               residuals[NMPC_STATE_DIM];
        return_t status_flag;
        real_t *state_ref = &ocp_state_reference[(i - 1u) * NMPC_STATE_DIM];
        real_t *control_ref =
                    &ocp_control_reference[(i - 1u) * NMPC_CONTROL_DIM];
        size_t j;

        /* Update state and control constraints */
        memcpy(z_low, ocp_lower_state_bound, sizeof(real_t) * NMPC_DELTA_DIM);
        memcpy(z_upp, ocp_upper_state_bound, sizeof(real_t) * NMPC_DELTA_DIM);

        #pragma MUST_ITERATE(NMPC_CONTROL_DIM, NMPC_CONTROL_DIM)
        for (j = 0; j < NMPC_CONTROL_DIM; j++) {
            z_low[NMPC_DELTA_DIM + j] = ocp_lower_control_bound[j] -
                                        control_ref[j];
            z_upp[NMPC_DELTA_DIM + j] = ocp_upper_control_bound[j] -
                                        control_ref[j];
        }

        /*
        Solve the IVP for the previous reference point to get the Jacobian
        (aka continuity constraint matrix, C) and integration residuals (c).

        We do this for the previous point because we need the current point to
        work out the residuals.
        */
        _solve_interval_ivp(state_ref, control_ref, jacobian, coeffs,
                            residuals);

        /* Copy the relevant data into the qpDUNES arrays. */
        status_flag = qpDUNES_setupRegularInterval(
            &ocp_qp_data.qpdata, ocp_qp_data.qpdata.intervals[i - 1u],
            ocp_state_weights, ocp_control_weights, jacobian, residuals,
            z_low, z_upp);
        assert(status_flag == QPDUNES_OK);

        /*
        Clear out the lambda guess for this interval to prevent old results
        from biasing the step size. Also clear out the clipping QP solver data
        for that interval.
        */
        memset(&ocp_qp_data.qpdata.lambda.data[(i - 1u) * NMPC_DELTA_DIM], 0,
               sizeof(float) * NMPC_DELTA_DIM);
        memset(&ocp_qp_data.interval_recs[i - 1u].clippingSolver_qStep_data,
               0, sizeof(float) * (NMPC_DELTA_DIM + NMPC_CONTROL_DIM));
        memset(&ocp_qp_data.interval_recs[i - 1u].clippingSolver_dz_data,
               0, sizeof(float) * (NMPC_DELTA_DIM + NMPC_CONTROL_DIM));

        /* Set up the QP in preparation for solving */
        qpDUNES_setupStageQP(
                &ocp_qp_data.qpdata, ocp_qp_data.qpdata.intervals[i - 1u]);
    }
}

#pragma FUNC_EXT_CALLED(nmpc_get_reference_point);
void nmpc_get_reference_point(real_t coeffs[NMPC_REFERENCE_DIM], uint32_t i) {
    assert(coeffs);
    assert(i <= OCP_HORIZON_LENGTH);

    memcpy(coeffs, &ocp_state_reference[i * NMPC_STATE_DIM],
           sizeof(real_t) * NMPC_STATE_DIM);
}

#pragma FUNC_EXT_CALLED(nmpc_set_wind_velocity);
void nmpc_set_wind_velocity(real_t x, real_t y, real_t z) {
    wind_velocity[0] = x;
    wind_velocity[1] = y;
    wind_velocity[2] = z;
}

#pragma FUNC_EXT_CALLED(nmpc_config_get_state_dim);
uint32_t nmpc_config_get_state_dim(void) {
    return NMPC_STATE_DIM;
}

#pragma FUNC_EXT_CALLED(nmpc_config_get_control_dim);
uint32_t nmpc_config_get_control_dim(void) {
    return NMPC_CONTROL_DIM;
}

#pragma FUNC_EXT_CALLED(nmpc_config_get_horizon_length);
uint32_t nmpc_config_get_horizon_length(void) {
    return OCP_HORIZON_LENGTH;
}

#pragma FUNC_EXT_CALLED(nmpc_config_get_step_length);
real_t nmpc_config_get_step_length(void) {
    return OCP_STEP_LENGTH;
}

#pragma FUNC_EXT_CALLED(nmpc_config_get_precision);
enum nmpc_precision_t nmpc_config_get_precision(void) {
#ifdef NMPC_SINGLE_PRECISION
    return NMPC_PRECISION_FLOAT;
#else
    return NMPC_PRECISION_DOUBLE;
#endif
}
