/*
Part of NMPC. Distributed under the MIT licence: see LICENSE in the root of
the repository.
*/

#ifndef NMPC_INNER_LOOP_H
#define NMPC_INNER_LOOP_H

#include "NMPC/Types.h"
#include <algorithm>
#include <cmath>

namespace NMPC {

/*
Building blocks of an inner loop that runs between solves. The loop is the sum
of three terms: the plan's command at this tick (the feedforward), an error law
on the distance between the state and the plan, and a correction from the
measured acceleration. These are the pieces they share.
*/

/* A second-order Butterworth low-pass, by the bilinear transform. */
template <class T>
struct Butterworth2 {
    T b0 = T(1), b1 = T(0), b2 = T(0), a1 = T(0), a2 = T(0);
    T x1 = T(0), x2 = T(0), y1 = T(0), y2 = T(0);

    /* Cutoff and sample rate in Hz. Clears the history. */
    void design(T cutoff, T rate) {
        const T pi_ = T(3.14159265358979323846);
        const T k = std::tan(pi_ * cutoff / rate), q = T(1) / std::sqrt(T(2)), n = T(1) / (T(1) + k / q + k * k);
        b0 = k * k * n;
        b1 = T(2) * b0;
        b2 = b0;
        a1 = T(2) * (k * k - T(1)) * n;
        a2 = (T(1) - k / q + k * k) * n;
        reset();
    }

    void reset() { x1 = x2 = y1 = y2 = T(0); }

    T operator()(T x) {
        const T y = b0 * x + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
        x2 = x1;
        x1 = x;
        y2 = y1;
        y1 = y;
        return y;
    }
};

/* Each element of u within [lower, upper]. */
template <class V>
V clamp(const V& u, const V& lower, const V& upper) {
    return u.cwiseMax(lower).cwiseMin(upper);
}

/* Whether any element of u is at or beyond its bound. */
template <class V>
bool on_bound(const V& u, const V& lower, const V& upper) {
    return (u.array() <= lower.array()).any() || (u.array() >= upper.array()).any();
}

/*
The incremental correction (INDI): what the model did not expect of the last
tick, low-passed, as actuator increments. Each of NY channels is a measured
quantity less the model's expectation of it, such as an angular acceleration
from the gyro's derivative or a specific force from the accelerometer. Each
channel has its own second-order low-pass and a gain that turns it into the
wanted change of what the actuators make (for a multirotor: the collective
thrust and the three moments, with gains of the mass and minus each inertia).
The allocation turns that into actuator increments, which are held within
per-actuator limits.

The model and its expectation are the application's: this type holds only the
filters and the arithmetic, so that one inner loop can use the model it plans
with. A model filter's force and moment belong in the expectation, so that the
correction holds only what they leave. What the correction holds is taken out
of the bounds and the actuator state that the solver sees.
*/
template <class T, int NY, int NA>
struct Incremental {
    Butterworth2<T> filter[NY];
    Vector<T, NY> gain = Vector<T, NY>::Ones();
    Matrix<T, NA, NY> allocation = Matrix<T, NA, NY>::Zero();
    Vector<T, NA> limit = Vector<T, NA>::Constant(T(1));
    Vector<T, NA> held = Vector<T, NA>::Zero();     /* the increments of the last update */

    /* One cutoff for every channel, in Hz, at the inner loop's rate. */
    void design(T cutoff, T rate) {
        for (auto& f : filter) {
            f.design(cutoff, rate);
        }
    }

    void reset() {
        for (auto& f : filter) {
            f.reset();
        }
        held.setZero();
    }

    /* The residuals of this tick, measured less expected; returns the actuator increments. */
    const Vector<T, NA>& update(const Vector<T, NY>& residual) {
        Vector<T, NY> wanted;
        for (int i = 0; i < NY; i++) {
            wanted(i) = gain(i) * filter[i](residual(i));
        }
        held = (allocation * wanted).cwiseMax(-limit).cwiseMin(limit);
        return held;
    }
};

/* No correction: holds nothing. */
template <class T, int NA>
struct NoCorrection {
    Vector<T, NA> held = Vector<T, NA>::Zero();
    void reset() {}
};

/*
An attitude and rate cascade, as the error law of a rigid body in place of the
solver's stored gains: the attitude error to a rate command, and a PID on the
rate to a moment about each body axis. The attitude error is the NED rotation
log(q * conj(q_plan)) taken into body axes. The planned rate is in the planned
attitude's body axes, so it is turned into the vehicle's own body axes,
conj(q) q_plan w_plan, before it is compared with the measured rate; the two
frames part whenever the attitude is off the plan. The rate's derivative term is
low-passed; the integral is held within a limit and stops while the command is
clipped. The defaults are the review prototype's for the F450: attitude 6, 6
and 3 per second, rate 18, 18 and 6 per second with an integral corner at
4 rad/s (2 in yaw) and a lead of 0.025 s for the rotor lag in roll and pitch.
The caller turns the moment into actuator commands through its allocation.
*/
template <class T>
struct AttitudeCascade {
    Vector<T, 3> k_attitude = Vector<T, 3>(T(6), T(6), T(3));     /* 1/s */
    Vector<T, 3> kp = Vector<T, 3>(T(18), T(18), T(6));           /* 1/s */
    Vector<T, 3> ki = Vector<T, 3>(T(72), T(72), T(12));          /* 1/s^2 */
    Vector<T, 3> kd = Vector<T, 3>(T(0.45), T(0.45), T(0));       /* s */
    Vector<T, 3> integral_limit = Vector<T, 3>(T(0.6), T(0.6), T(0.1));  /* N m */
    T h = T(0.001);
    T derivative_decay = T(0);

    Vector<T, 3> integral = Vector<T, 3>::Zero();   /* N m */
    Vector<T, 3> error = Vector<T, 3>::Zero();      /* the rate error of the last command, rad/s */
    Vector<T, 3> error_previous = Vector<T, 3>::Zero();
    Vector<T, 3> derivative = Vector<T, 3>::Zero();

    /* The tick, s, and the cutoff of the derivative's low-pass, Hz. Clears the state. */
    void design(T tick, T derivative_cutoff = T(80)) {
        h = tick;
        derivative_decay = std::exp(-tick * T(2) * T(3.14159265358979323846) * derivative_cutoff);
        reset();
    }

    void reset() {
        integral.setZero();
        error.setZero();
        error_previous.setZero();
        derivative.setZero();
    }

    /*
    The moment about the body axes, N m, toward the planned attitude and rate, for the inertia's diagonal.
    w is the measured rate in the body axes of q, and w_plan the planned rate in the body axes of q_plan.
    */
    Vector<T, 3> moment(const Quaternion<T>& q, const Vector<T, 3>& w, const Quaternion<T>& q_plan,
            const Vector<T, 3>& w_plan, const Vector<T, 3>& inertia) {
        const Vector<T, 3> theta = so3_log<T>(q * q_plan.conjugate());
        const Quaternion<T> q_body = q.conjugate();
        const Vector<T, 3> planned(q_body * (q_plan * w_plan));
        const Vector<T, 3> rate_command = planned - k_attitude.cwiseProduct(Vector<T, 3>(q_body * theta));
        error = rate_command - w;
        derivative = derivative_decay * derivative + (T(1) - derivative_decay) * (error - error_previous) / h;
        error_previous = error;
        return inertia.cwiseProduct(kp.cwiseProduct(error) + kd.cwiseProduct(derivative)) + integral;
    }

    /* After the command is formed: the integral moves unless the command was clipped. */
    void integrate(bool clipped, const Vector<T, 3>& inertia) {
        if (!clipped) {
            integral += h * inertia.cwiseProduct(ki.cwiseProduct(error));
            integral = integral.cwiseMax(-integral_limit).cwiseMin(integral_limit);
        }
    }
};

}

#endif
