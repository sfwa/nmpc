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

#include "nmpc/types.h"
#include "nmpc/state.h"
#include "nmpc/integrator.h"
#include "nmpc/dynamics.h"

#include "cnmpc.h"

static FixedWingFlightDynamicsModel fixed_wing_model =
    FixedWingFlightDynamicsModel();
static State current;
#if defined(NMPC_INTEGRATOR_RK4)
    IntegratorRK4 integrator;
#elif defined(NMPC_INTEGRATOR_HEUN)
    IntegratorHeun integrator;
#elif defined(NMPC_INTEGRATOR_EULER)
    IntegratorEuler integrator;
#endif

void nmpc_init() {

}

void nmpc_set_position(real_t lat, real_t lon, real_t alt) {
    current.position() << lat, lon, alt;
}

void nmpc_set_velocity(real_t x, real_t y, real_t z) {
    current.velocity() << x, y, z;
}

void nmpc_set_acceleration(real_t x, real_t y, real_t z) {
    current.acceleration() << x, y, z;
}

void nmpc_set_attitude(real_t w, real_t x, real_t y, real_t z) {
    current.attitude() << x, y, z, w;
}

void nmpc_set_angular_velocity(real_t x, real_t y, real_t z) {
    current.angular_velocity() << x, y, z;
}

void nmpc_set_angular_acceleration(real_t x, real_t y, real_t z) {
    current.angular_acceleration() << x, y, z;
}

void nmpc_set_wind_velocity(real_t x, real_t y, real_t z) {
    current.wind_velocity() << x, y, z;
}

void nmpc_get_state(struct nmpc_state_t *in) {
    in->position[0] = current.position()[0];
    in->position[1] = current.position()[1];
    in->position[2] = current.position()[2];
    in->velocity[0] = current.velocity()[0];
    in->velocity[1] = current.velocity()[1];
    in->velocity[2] = current.velocity()[2];
    in->acceleration[0] = current.acceleration()[0];
    in->acceleration[1] = current.acceleration()[1];
    in->acceleration[2] = current.acceleration()[2];
    in->attitude[0] = current.attitude()[0];
    in->attitude[1] = current.attitude()[1];
    in->attitude[2] = current.attitude()[2];
    in->attitude[3] = current.attitude()[3];
    in->angular_velocity[0] = current.angular_velocity()[0];
    in->angular_velocity[1] = current.angular_velocity()[1];
    in->angular_velocity[2] = current.angular_velocity()[2];
    in->angular_acceleration[0] = current.angular_acceleration()[0];
    in->angular_acceleration[1] = current.angular_acceleration()[1];
    in->angular_acceleration[2] = current.angular_acceleration()[2];
    in->wind_velocity[0] = current.wind_velocity()[0];
    in->wind_velocity[1] = current.wind_velocity()[1];
    in->wind_velocity[2] = current.wind_velocity()[2];
}

void nmpc_set_state(struct nmpc_state_t *in) {
    current <<
        in->position[0],
        in->position[1],
        in->position[2],
        in->velocity[0],
        in->velocity[1],
        in->velocity[2],
        in->acceleration[0],
        in->acceleration[1],
        in->acceleration[2],
        in->attitude[0],
        in->attitude[1],
        in->attitude[2],
        in->attitude[3],
        in->angular_velocity[0],
        in->angular_velocity[1],
        in->angular_velocity[2],
        in->angular_acceleration[0],
        in->angular_acceleration[1],
        in->angular_acceleration[2],
        in->wind_velocity[0],
        in->wind_velocity[1],
        in->wind_velocity[2];
}

void nmpc_integrate(float dt, real_t control_vector[NMPC_CONTROL_DIM]) {
    AccelerationVector temp = fixed_wing_model.evaluate(current, 
        Eigen::Matrix<real_t, NMPC_CONTROL_DIM, 1>(control_vector));
    current.acceleration() << temp.segment<3>(0);
    current.angular_acceleration() << temp.segment<3>(3);
    current = integrator.integrate(current, dt);
}

void nmpc_fixedwingdynamics_set_mass(real_t mass) {
    fixed_wing_model.set_mass(mass);
}

void nmpc_fixedwingdynamics_set_inertia_tensor(real_t inertia_tensor[9]) {
    fixed_wing_model.set_inertia_tensor(Matrix3x3r(inertia_tensor));
}

void nmpc_fixedwingdynamics_set_prop_coeffs(real_t in_prop_area,
real_t in_prop_cve){
    fixed_wing_model.set_prop_coeffs(in_prop_area, in_prop_cve);
}

void nmpc_fixedwingdynamics_set_drag_coeffs(real_t coeffs[5]) {
    fixed_wing_model.set_drag_coeffs(Vector5r(coeffs));
}

void nmpc_fixedwingdynamics_set_lift_coeffs(real_t coeffs[5]) {
    fixed_wing_model.set_lift_coeffs(Vector5r(coeffs));
}

void nmpc_fixedwingdynamics_set_side_coeffs(real_t coeffs[4],
real_t control[NMPC_CONTROL_DIM]) {
    fixed_wing_model.set_side_coeffs(Vector4r(coeffs),
        Vector4r(control));
}

void nmpc_fixedwingdynamics_set_pitch_moment_coeffs(real_t coeffs[2],
real_t control[NMPC_CONTROL_DIM]) {
    fixed_wing_model.set_pitch_moment_coeffs(Vector2r(coeffs),
        Vector4r(control));
}

void nmpc_fixedwingdynamics_set_roll_moment_coeffs(real_t coeffs[1],
real_t control[NMPC_CONTROL_DIM]) {
    fixed_wing_model.set_roll_moment_coeffs(Vector1r(coeffs),
        Vector4r(control));
}

void nmpc_fixedwingdynamics_set_yaw_moment_coeffs(real_t coeffs[2],
real_t control[NMPC_CONTROL_DIM]) {
    fixed_wing_model.set_yaw_moment_coeffs(Vector2r(coeffs),
        Vector4r(control));
}

uint32_t nmpc_config_get_state_dim() {
    return NMPC_STATE_DIM;
}

uint32_t nmpc_config_get_control_dim() {
    return NMPC_CONTROL_DIM;
}

enum nmpc_precision_t nmpc_config_get_precision() {
#ifdef NMPC_SINGLE_PRECISION
    return NMPC_PRECISION_FLOAT;
#else
    return NMPC_PRECISION_DOUBLE;
#endif
}
