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

#ifndef INTERFACE_H
#define INTERFACE_H

#include <stdint.h>
#include "config.h"

#ifdef NMPC_SINGLE_PRECISION
typedef float real_t;
#endif

#ifdef NMPC_DOUBLE_PRECISION
typedef double real_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct nmpc_state_t {
    real_t position[3];
    real_t velocity[3];
    real_t attitude[4]; /* w, x, y, z */
    real_t angular_velocity[3];
};

void nmpc_init();
void nmpc_preparation_step();
void nmpc_feedback_step(real_t measurement[NMPC_STATE_DIM]);
void nmpc_get_controls(real_t controls[NMPC_CONTROL_DIM]);
void nmpc_update_horizon(real_t new_reference[NMPC_REFERENCE_DIM]);

/* Functions for setting weights and bounds for the OCP solver. */
void nmpc_set_state_weights(real_t coeffs[NMPC_DELTA_DIM]);
void nmpc_set_control_weights(real_t coeffs[NMPC_CONTROL_DIM]);
void nmpc_set_terminal_weights(real_t coeffs[NMPC_DELTA_DIM]);
void nmpc_set_lower_control_bound(real_t coeffs[NMPC_CONTROL_DIM]);
void nmpc_set_upper_control_bound(real_t coeffs[NMPC_CONTROL_DIM]);
void nmpc_set_reference_point(real_t coeffs[NMPC_REFERENCE_DIM], uint32_t i);

/* Functions for setting different parts of the state vector. */
void nmpc_fixedwingdynamics_set_position(
    real_t lat, real_t lon, real_t alt);
void nmpc_fixedwingdynamics_set_velocity(
    real_t x, real_t y, real_t z);
void nmpc_fixedwingdynamics_set_attitude(
    real_t w, real_t x, real_t y, real_t z);
void nmpc_fixedwingdynamics_set_angular_velocity(
    real_t x, real_t y, real_t z);

/* Functions for getting the state vector and covariance. */
void nmpc_fixedwingdynamics_set_state(struct nmpc_state_t *in);
void nmpc_fixedwingdynamics_get_state(struct nmpc_state_t *in);

void nmpc_fixedwingdynamics_integrate(
    float dt, real_t control_vector[NMPC_CONTROL_DIM]);

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
