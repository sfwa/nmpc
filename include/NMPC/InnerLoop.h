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

}

#endif
