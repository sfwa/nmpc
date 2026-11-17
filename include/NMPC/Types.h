/*
Part of NMPC. Distributed under the MIT licence: see LICENSE in the root of
the repository.
*/

#ifndef NMPC_TYPES_H
#define NMPC_TYPES_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>

#include "NMPC/Config.h"

namespace NMPC {

/* Fixed-size aliases on a scalar T. Every size in NMPC is known at compile
   time. */
template <class T, int Rows>
using Vector = Eigen::Matrix<T, Rows, 1>;

template <class T, int Rows, int Columns>
using Matrix = Eigen::Matrix<T, Rows, Columns>;

template <class T>
using Quaternion = Eigen::Quaternion<T>;

/*
Manifold helpers for states that hold an attitude.

NMPC's words for the two operations on a state are retract and difference,
where UKF says apply_delta and delta:

    retract(x, d)       moves the state x by the tangent d
    difference(x, y)    the tangent d at y with retract(y, d) == x

An attitude is a Hamilton quaternion, and a tangent on it is a small rotation
applied on the left, in the frame the quaternion rotates into:

    retract(q, phi)     = exp(phi) * q
    difference(q, r)    = log(q * conj(r))

These are the conventions of the multirotor and fixed-wing reference models
(attitude from body to NED, perturbed in NED). A problem with other
conventions writes its own retract and difference and does not use these.
*/

/* The skew-symmetric matrix of a, so that hat(a) * b == a.cross(b). */
template <class T>
Matrix<T, 3, 3> hat(const Vector<T, 3>& a) {
  Matrix<T, 3, 3> m;
  m << T(0), -a.z(), a.y(), a.z(), T(0), -a.x(), -a.y(), a.x(), T(0);
  return m;
}

/* The unit quaternion of the rotation vector phi (rad). Series near zero keep
   it exact in single precision. */
template <class T>
Quaternion<T> so3_exp(const Vector<T, 3>& phi) {
  const T angle2 = phi.squaredNorm();
  T c, s;
  if (angle2 < T(1e-4)) {
    c = T(1) - angle2 / T(8) + angle2 * angle2 / T(384);
    s = T(0.5) - angle2 / T(48) + angle2 * angle2 / T(3840);
  } else {
    const T angle = std::sqrt(angle2);
    c = std::cos(angle / T(2));
    s = std::sin(angle / T(2)) / angle;
  }
  return Quaternion<T>(c, s * phi.x(), s * phi.y(), s * phi.z());
}

/* The rotation vector of a unit quaternion, taking the shorter way round (angle
   at most pi). */
template <class T>
Vector<T, 3> so3_log(const Quaternion<T>& q_in) {
  Quaternion<T> q = q_in;
  if (q.w() < T(0)) {
    q.coeffs() = -q.coeffs();
  }
  const T n2 = q.vec().squaredNorm();
  T k;
  if (n2 < T(1e-6)) {
    k = T(2) / q.w() * (T(1) - n2 / (T(3) * q.w() * q.w()));
  } else {
    const T n = std::sqrt(n2);
    k = T(2) * std::atan2(n, q.w()) / n;
  }
  return k * q.vec();
}

/*
The inverse of the left Jacobian of SO(3) at phi, to second order:
d log(exp(delta) * exp(phi)) / d delta at delta = 0. It is what turns a
change of attitude into a change of an attitude error log(q * conj(r)).
*/
template <class T>
Matrix<T, 3, 3> so3_left_jacobian_inverse(const Vector<T, 3>& phi) {
  const Matrix<T, 3, 3> h = hat<T>(phi);
  return Matrix<T, 3, 3>::Identity() - T(0.5) * h + h * h / T(12);
}

/* Attitude retract and difference with the conventions above. */
template <class T>
Quaternion<T> attitude_retract(const Quaternion<T>& q,
                               const Vector<T, 3>& phi) {
  return (so3_exp<T>(phi) * q).normalized();
}

template <class T>
Vector<T, 3> attitude_difference(const Quaternion<T>& q,
                                 const Quaternion<T>& r) {
  return so3_log<T>(q * r.conjugate());
}

}  // namespace NMPC

#endif
