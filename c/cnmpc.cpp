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
#include "state.h"
#include "integrator.h"
#include "dynamics.h"
#include "ocp.h"

#include "cnmpc.h"

static X8DynamicsModel dynamics_model =
    X8DynamicsModel();
static State current;
#if defined(NMPC_INTEGRATOR_RK4)
    IntegratorRK4 integrator;
#elif defined(NMPC_INTEGRATOR_HEUN)
    IntegratorHeun integrator;
#elif defined(NMPC_INTEGRATOR_EULER)
    IntegratorEuler integrator;
#endif

static OptimalControlProblem ocp =
    OptimalControlProblem(&dynamics_model);

void nmpc_init() {
    ocp.initialise();
}

void nmpc_preparation_step() {
    ocp.preparation_step();
}

void nmpc_feedback_step(real_t measurement[NMPC_STATE_DIM]) {
    Eigen::Map<StateVector> measurement_map =
        Eigen::Map<StateVector>(measurement);
    StateVector m = measurement_map;
    ocp.feedback_step(m);
}

void nmpc_set_state_weights(real_t coeffs[NMPC_DELTA_DIM]) {
    Eigen::Map<DeltaVector> state_weight_map =
        Eigen::Map<DeltaVector>(coeffs);
    DeltaVector state_weight = state_weight_map;
    ocp.set_state_weights(state_weight);
}

void nmpc_set_control_weights(real_t coeffs[NMPC_CONTROL_DIM]) {
    Eigen::Map<ControlVector> control_weight_map =
        Eigen::Map<ControlVector>(coeffs);
    ControlVector control_weight = control_weight_map;
    ocp.set_control_weights(control_weight);
}

void nmpc_set_terminal_weights(real_t coeffs[NMPC_DELTA_DIM]) {
    Eigen::Map<DeltaVector> terminal_weight_map =
        Eigen::Map<DeltaVector>(coeffs);
    DeltaVector terminal_weight = terminal_weight_map;
    ocp.set_terminal_weights(terminal_weight);
}

void nmpc_set_lower_control_bound(real_t coeffs[NMPC_CONTROL_DIM]) {
    Eigen::Map<ControlConstraintVector> control_constraint_map =
        Eigen::Map<ControlConstraintVector>(coeffs);
    ControlConstraintVector control_constraint = control_constraint_map;
    ocp.set_lower_control_bound(control_constraint);
}

void nmpc_set_upper_control_bound(real_t coeffs[NMPC_CONTROL_DIM]) {
    Eigen::Map<ControlConstraintVector> control_constraint_map =
        Eigen::Map<ControlConstraintVector>(coeffs);
    ControlConstraintVector control_constraint = control_constraint_map;
    ocp.set_upper_control_bound(control_constraint);
}

void nmpc_set_reference_point(real_t coeffs[NMPC_REFERENCE_DIM],
uint32_t i) {
    Eigen::Map<ReferenceVector> reference_map =
        Eigen::Map<ReferenceVector>(coeffs);
    ReferenceVector reference = reference_map;
    ocp.set_reference_point(reference, i);
}

void nmpc_fixedwingdynamics_set_position(
real_t lat, real_t lon, real_t alt) {
    current.position() << lat, lon, alt;
}

void nmpc_fixedwingdynamics_set_velocity(
real_t x, real_t y, real_t z) {
    current.velocity() << x, y, z;
}

void nmpc_fixedwingdynamics_set_attitude(
real_t w, real_t x, real_t y, real_t z) {
    current.attitude() << x, y, z, w;
}

void nmpc_fixedwingdynamics_set_angular_velocity(
real_t x, real_t y, real_t z) {
    current.angular_velocity() << x, y, z;
}

void nmpc_fixedwingdynamics_get_state(struct nmpc_state_t *in) {
    in->position[0] = current.position()[0];
    in->position[1] = current.position()[1];
    in->position[2] = current.position()[2];
    in->velocity[0] = current.velocity()[0];
    in->velocity[1] = current.velocity()[1];
    in->velocity[2] = current.velocity()[2];
    in->attitude[0] = current.attitude()[0];
    in->attitude[1] = current.attitude()[1];
    in->attitude[2] = current.attitude()[2];
    in->attitude[3] = current.attitude()[3];
    in->angular_velocity[0] = current.angular_velocity()[0];
    in->angular_velocity[1] = current.angular_velocity()[1];
    in->angular_velocity[2] = current.angular_velocity()[2];
}

void nmpc_fixedwingdynamics_set_state(struct nmpc_state_t *in) {
    current <<
        in->position[0],
        in->position[1],
        in->position[2],
        in->velocity[0],
        in->velocity[1],
        in->velocity[2],
        in->attitude[0],
        in->attitude[1],
        in->attitude[2],
        in->attitude[3],
        in->angular_velocity[0],
        in->angular_velocity[1],
        in->angular_velocity[2];
}

void nmpc_fixedwingdynamics_integrate(
float dt, real_t control_vector[NMPC_CONTROL_DIM]) {
    current = integrator.integrate(
        current,
        ControlVector(control_vector),
        &dynamics_model,
        dt);
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
