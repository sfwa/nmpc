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

#ifndef STATE_H
#define STATE_H

#include "types.h"

/* DynamicsModel forward declaration. */
class DynamicsModel;

/*
Definition for filter state vector.
Contents are as follows:
    - Position (3-vector, m, NED frame)
    - Linear Velocity (3-vector, m/s, NED frame)
    - Attitude (quaternion (x, y, z, w), describes rotation from local NED
      frame to body frame.)
    - Angular Velocity (3-vector, rad/s, body frame)
*/
class State: public StateVector {
public:
    State() : StateVector() {}

    template<typename OtherDerived>
    State(const Eigen::MatrixBase<OtherDerived>& other) :
        StateVector(other) { }

    template<typename OtherDerived>
    State & operator= (const Eigen::MatrixBase<OtherDerived>& other)
    {
        StateVector::operator=(other);
        return *this;
    }

    const StateVectorDerivative model(ControlVector c, DynamicsModel *d);

    /* Read-only accessors */
    const Vector3r position() const {
        return segment<3>(0);
    }

    const Vector3r velocity() const {
        return segment<3>(3);
    }

    const Vector4r attitude() const {
        return segment<4>(6);
    }

    const Vector3r angular_velocity() const {
        return segment<3>(10);
    }

    /* Mutable accessors */
    Eigen::VectorBlock<StateVector, 3> position() {
        return segment<3>(0);
    }

    Eigen::VectorBlock<StateVector, 3> velocity() {
        return segment<3>(3);
    }

    Eigen::VectorBlock<StateVector, 4> attitude() {
        return segment<4>(6);
    }

    Eigen::VectorBlock<StateVector, 3> angular_velocity() {
        return segment<3>(10);
    }
};

#endif
