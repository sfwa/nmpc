/*

Some code derived from TI's MATHLIB; original copyright notice follows.

MATHLIB -- TI Floating-Point Math Function Library

Copyright (C) 2011 Texas Instruments Incorporated - http://www.ti.com/

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:

   Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

   Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the
   distribution.

   Neither the name of Texas Instruments Incorporated nor the names of
   its contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef C66MATH_H
#define C66MATH_H


#ifdef USE_C66_INTRINSICS

static inline float recip_sqrt_f(float a) {
    const float half = 0.5f;
    const float onep5 = 1.5f;
    float x0, x1, x2, x3;

    if (a < 1.17549435e-38f) {
        x2 = _itof(0x7F800000);
    } else {
        x0 = _rsqrsp(a);
        x1 = a * x0;
        x3 = onep5 - (x1 * x0 * half);
        x1 = x0 * x3;
        x2 = x1 * (onep5  - (a * x1 * x1 * half));
    }

    return x2;
}

static inline float sqrt_f(float a) {
    const float half = 0.5f;
    const float onep5 = 1.5f;
    const float maxe = 3.402823466E+38f;
    float x0, x1, x2, x3;

    if (a <= 0.0f) {
        x2 = 0.0f;
    } else if (a > maxe) {
        x2 = maxe;
    } else {
        x0 = _rsqrsp(a);

        x1 = a * x0;
        x3 = onep5 - (x1 * x0 * half);
        x1 = x0 * x3;
        x2 = x1 * (onep5 - (a * x1 * x1 * half));
        x2 = a * x2;
    }

    return x2;
}

static inline float divide_f(float a, float b) {
    const float two = 2.0f;
    const float maxe = 3.402823466E+38f;
    float x;

    if (a == 0.0f) {
        x = 0.0f;
    } else if (_fabsf(b) > maxe && _fabsf(a) <= maxe) {
        x = 0.0f;
    } else {
        x = _rcpsp(b);
        x = x * (two - b * x);
        x = x * (two - b * x);
        x = a * x;
    }

    return x;
}

static inline float recip_f(float b) {
    const float two = 2.0f;
    float x;

    if (_fabsf(b) > 3.402823466E+38f) {
        x = 0.0f;
    } else {
        x = _rcpsp(b);
        x = x * (two - b * x);
        x = x * (two - b * x);
    }

    return x;
}

static inline float abs_f(float x) {
	return _fabsf(x);
}

#else

#define recip_sqrt_f(x) ((float)(1.0 / sqrt((x))))
#define divide_f(a, b) ((float)((a) / (b)))
#define recip_f(a) ((float)(1.0 / (a)))
#define sqrt_f(a) ((float)sqrt((a)))
#define abs_f(x) ((float)fabs(x))

#endif

/* Non-TI compatibility */
#ifndef __TI_COMPILER_VERSION__
#define _nassert(x) assert(x)
#endif

#endif
