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

#ifndef CNMPC_INTERFACE_H
#define CNMPC_INTERFACE_H

#include <stdint.h>
#include "config.h"

#if defined(NMPC_SINGLE_PRECISION) && !defined(__USE_SINGLE_PRECISION__)
typedef float real_t;
#elif defined(NMPC_DOUBLE_PRECISION) && !defined(__USE_SINGLE_PRECISION__)
typedef double real_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct nmpc_state_t {
    real_t position[3];
    real_t velocity[3];
    real_t attitude[4]; /* x, y, z, W */
    real_t angular_velocity[3];
};

enum nmpc_result_t {
    NMPC_OK,
    NMPC_INFEASIBLE,
    NMPC_ERROR
};

/*
If use_relative_positions is true, the position component of each state vector
in the reference trajectory is assumed to be relative to the previous state
vector.
*/
void nmpc_init(bool use_relative_positions);
void nmpc_preparation_step(void);
void nmpc_feedback_step(real_t measurement[NMPC_STATE_DIM]);
enum nmpc_result_t nmpc_get_controls(real_t controls[NMPC_CONTROL_DIM]);
void nmpc_update_horizon(real_t new_reference[NMPC_REFERENCE_DIM]);

/* Functions for setting weights and bounds for the OCP solver. */
void nmpc_set_state_weights(real_t coeffs[NMPC_DELTA_DIM]);
void nmpc_set_control_weights(real_t coeffs[NMPC_CONTROL_DIM]);
void nmpc_set_terminal_weights(real_t coeffs[NMPC_DELTA_DIM]);
void nmpc_set_lower_control_bound(real_t coeffs[NMPC_CONTROL_DIM]);
void nmpc_set_upper_control_bound(real_t coeffs[NMPC_CONTROL_DIM]);
void nmpc_set_reference_point(real_t coeffs[NMPC_REFERENCE_DIM],
uint32_t i);

real_t nmpc_get_objective_value(void);

/* Function to set the wind estimate for the dynamics model. */
void nmpc_set_wind_velocity(real_t x, real_t y, real_t z);

/* Functions for setting different parts of the state vector. */
void nmpc_fixedwingdynamics_set_position(
    real_t lat, real_t lon, real_t alt);
void nmpc_fixedwingdynamics_set_velocity(
    real_t x, real_t y, real_t z);
void nmpc_fixedwingdynamics_set_attitude(
    real_t w, real_t x, real_t y, real_t z);
void nmpc_fixedwingdynamics_set_angular_velocity(
    real_t x, real_t y, real_t z);

/* Functions for getting the state vector. */
void nmpc_fixedwingdynamics_set_state(struct nmpc_state_t *in);
void nmpc_fixedwingdynamics_get_state(struct nmpc_state_t *in);

void nmpc_fixedwingdynamics_integrate(float dt,
real_t control_vector[NMPC_CONTROL_DIM]);

/*
Functions to access the compiled configuration
*/
enum nmpc_precision_t {
    NMPC_PRECISION_FLOAT = 0,
    NMPC_PRECISION_DOUBLE = 1
};

uint32_t nmpc_config_get_state_dim(void);
uint32_t nmpc_config_get_control_dim(void);
uint32_t nmpc_config_get_horizon_length(void);
real_t nmpc_config_get_step_length(void);
enum nmpc_precision_t nmpc_config_get_precision(void);

#ifdef __cplusplus
}
#endif

#endif
