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

#ifndef CONFIG_H_
#define CONFIG_H_

/* Uncomment one of the following to select the float-point precision */
/* #define NMPC_SINGLE_PRECISION */
#define NMPC_DOUBLE_PRECISION

/*
Choose the integration method here:
    - NMPC_INTEGRATOR_RK4: 4th-order integration method. Requires most CPU time.
    - NMPC_INTEGRATOR_HEUN: 2nd-order Heun's method. Requires middle CPU time.
    - NMPC_INTEGRATOR_EULER: 1st-orer Euler's method. Requires least CPU time.
*/

#define NMPC_INTEGRATOR_RK4
/* #define NMPC_INTEGRATOR_HEUN */
/* #define NMPC_INTEGRATOR_EULER */

/* NMPC vector dimensioning. */
#define NMPC_CONTROL_DIM 4
#define NMPC_STATE_DIM 25

#endif
