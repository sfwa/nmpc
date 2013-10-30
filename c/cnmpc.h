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
#include "nmpc/config.h"

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
    real_t acceleration[3];
    real_t attitude[4]; /* w, x, y, z */
    real_t angular_velocity[3];
    real_t angular_acceleration[3];
    real_t wind_velocity[3];
    real_t gyro_bias[3];
};

/* Functions for setting different parts of the state vector. */
void nmpc_set_position(real_t lat, real_t lon, real_t alt);
void nmpc_set_velocity(real_t x, real_t y, real_t z);
void nmpc_set_acceleration(real_t x, real_t y, real_t z);
void nmpc_set_attitude(real_t w, real_t x, real_t y, real_t z);
void nmpc_set_angular_velocity(real_t x, real_t y, real_t z);
void nmpc_set_angular_acceleration(real_t x, real_t y, real_t z);
void nmpc_set_wind_velocity(real_t x, real_t y, real_t z);
void nmpc_set_gyro_bias(real_t x, real_t y, real_t z);

/* Functions for getting the state vector and covariance. */
void nmpc_set_state(struct nmpc_state_t *in);
void nmpc_get_state(struct nmpc_state_t *in);

void nmpc_integrate(float dt, real_t control_vector[NMPC_CONTROL_DIM]);

/*
Functions to set airframe properties and coefficients for the fixed-wing
dynamics model.
*/
void nmpc_fixedwingdynamics_set_mass(real_t mass);
void nmpc_fixedwingdynamics_set_inertia_tensor(real_t inertia_tensor[9]);
void nmpc_fixedwingdynamics_set_prop_coeffs(real_t in_prop_area,
    real_t in_prop_cve);
void nmpc_fixedwingdynamics_set_drag_coeffs(real_t coeffs[5]);
void nmpc_fixedwingdynamics_set_lift_coeffs(real_t coeffs[5]);
void nmpc_fixedwingdynamics_set_side_coeffs(real_t coeffs[4],
    real_t control[NMPC_CONTROL_DIM]);
void nmpc_fixedwingdynamics_set_pitch_moment_coeffs(real_t coeffs[2],
    real_t control[NMPC_CONTROL_DIM]);
void nmpc_fixedwingdynamics_set_roll_moment_coeffs(real_t coeffs[1],
    real_t control[NMPC_CONTROL_DIM]);
void nmpc_fixedwingdynamics_set_yaw_moment_coeffs(real_t coeffs[2],
    real_t control[NMPC_CONTROL_DIM]);

/*
Functions to access the compiled configuration
*/
enum nmpc_precision_t {
    NMPC_PRECISION_FLOAT = 0,
    NMPC_PRECISION_DOUBLE = 1
};

uint32_t nmpc_config_get_state_dim(void);
uint32_t nmpc_config_get_control_dim(void);
enum nmpc_precision_t nmpc_config_get_precision(void);

#ifdef __cplusplus
}
#endif

#endif
