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

#ifndef TYPES_H_
#define TYPES_H_

#include "config.h"

/*
Definitions for numerical precision. Also included is a definition for the
fourth root of the machine precision, which is used when calculating the
finite-difference interval for the Hessian calculation.
*/

#ifdef NMPC_SINGLE_PRECISION
typedef float real_t;
#define NMPC_EPS_SQRT ((real_t)3.450e-4)
#define NMPC_EPS_4RT ((real_t)1.857e-2)
#endif

#ifdef NMPC_DOUBLE_PRECISION
typedef double real_t;
#define NMPC_EPS_SQRT ((real_t)1.490e-8)
#define NMPC_EPS_4RT ((real_t)1.221e-4)
#endif

#define G_ACCEL ((real_t)9.80665)
#define RHO ((real_t)1.225)

#include <Eigen/Core>
#include <Eigen/Geometry>

typedef Eigen::Matrix<real_t, 1, 1> Vector1r;
typedef Eigen::Matrix<real_t, 2, 1> Vector2r;
typedef Eigen::Matrix<real_t, 3, 1> Vector3r;
typedef Eigen::Matrix<real_t, 4, 1> Vector4r;
typedef Eigen::Matrix<real_t, 5, 1> Vector5r;
typedef Eigen::Quaternion<real_t> Quaternionr;
typedef Eigen::Matrix<real_t, 3, 3> Matrix3x3r;
typedef Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> MatrixXr;
typedef Eigen::Matrix<real_t, Eigen::Dynamic, 1> VectorXr;

typedef Eigen::Matrix<real_t, NMPC_STATE_DIM, 1> StateVector;
typedef Eigen::Matrix<real_t, NMPC_STATE_DIM, 1> StateVectorDerivative;

/*
Typedef for a gradient vector. Dimension is smaller than a reference vector
by one because orientation is represented by MRP rather than quaternion.
 */
typedef Eigen::Matrix<real_t, NMPC_GRADIENT_DIM, 1> GradientVector;

/* Hessian matrices for control and state vectors. */
typedef Eigen::Matrix<
    real_t,
    NMPC_GRADIENT_DIM,
    NMPC_GRADIENT_DIM> HessianMatrix;

/* Matrices for state and control weights. */
typedef Eigen::Matrix<
    real_t,
    NMPC_STATE_DIM-1,
    NMPC_STATE_DIM-1> StateWeightMatrix;

typedef Eigen::Matrix<
    real_t,
    NMPC_CONTROL_DIM,
    NMPC_CONTROL_DIM> ControlWeightMatrix;

/* Typedef for control vector. */
typedef Eigen::Matrix<real_t, NMPC_CONTROL_DIM, 1> ControlVector;

/* Reference vector is for both state and control. */
typedef Eigen::Matrix<real_t, NMPC_REFERENCE_DIM, 1> ReferenceVector;

/* Typedef for delta vector. */
typedef Eigen::Matrix<real_t, NMPC_DELTA_DIM, 1> DeltaVector;

/* Typedef for constraint vector. */
typedef Eigen::Matrix<
    real_t,
    NMPC_GRADIENT_DIM,
    1> InequalityConstraintVector;

/*
Typedef for inequality constraint matrix. Same dimensions as the Hessian
Matrix.
*/
typedef Eigen::Matrix<
    real_t,
    NMPC_GRADIENT_DIM,
    NMPC_GRADIENT_DIM> InequalityConstraintMatrix;

/* Typedef for continuity constraint matrix.  */
typedef Eigen::Matrix<
    real_t,
    NMPC_STATE_DIM-1,
    NMPC_GRADIENT_DIM> ContinuityConstraintMatrix;

typedef Eigen::Matrix<real_t, 6, 1> AccelerationVector;

#endif
