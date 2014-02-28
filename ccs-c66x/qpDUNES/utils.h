/*
 *  This file is part of qpDUNES.
 *
 *  qpDUNES -- A DUal NEwton Strategy for convex quadratic programming.
 *  Copyright (C) 2012 by Janick Frasch, Hans Joachim Ferreau et al.
 *  All rights reserved.
 *
 *  qpDUNES is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  qpDUNES is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with qp42; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#ifndef QPDUNES_UTILS_H
#define QPDUNES_UTILS_H

#include <stddef.h>
#include "types.h"
#include "../c66math.h"


sparsityType_t qpDUNES_detectMatrixSparsity(const real_t* const M,
size_t nRows, size_t nCols);

return_t qpDUNES_updateMatrixData(matrix_t* const to,
const real_t* const from, size_t nRows, size_t nCols);

return_t qpDUNES_setupZeroMatrix(size_t nRows, size_t nCols, matrix_t* to);

return_t qpDUNES_setMatrixNull(matrix_t* const matrix);

return_t qpDUNES_setupScaledIdentityMatrix(size_t nRows, real_t scalar,
matrix_t* to);

return_t qpDUNES_setupVector(vector_t* const to, const real_t* const from,
size_t n);

return_t qpDUNES_updateVector(vector_t* const to, const real_t* const from,
size_t n);

return_t qpDUNES_updateSimpleBoundVector(qpData_t* qpData, vector_t* const to,
const real_t* const dBnd, const real_t* const xBnd, const real_t* const uBnd);

return_t qpDUNES_copyMatrix(matrix_t* const to, const matrix_t* const from,
size_t dim0, size_t dim1);

return_t qpDUNES_makeMatrixDense(matrix_t* const M, size_t dim0, size_t dim1);

return_t qpDUNES_transposeMatrix(matrix_t* const to,
const matrix_t* const from, size_t dim0, size_t dim1);

static inline void qpDUNES_setupZeroVector(vector_t* restrict const to,
size_t n) {
    assert(to && to->data);
    _nassert((size_t)to->data % 4 == 0);

    size_t i;

    for (i = 0; i < n; i++) {
        to->data[i] = 0.0;
    }
}

static inline void qpDUNES_setupUniformVector(vector_t* restrict const to,
real_t value, size_t n) {
    assert(to && to->data);
    _nassert((size_t)to->data % 4 == 0);

    size_t i;

    for (i = 0; i < n; i++) {
        to->data[i] = value;
    }
}

static inline void qpDUNES_copyVector(vector_t* restrict const to,
const vector_t* const from, size_t n) {
    assert(to && from && to->data && from->data);
    _nassert((size_t)to->data % 4 == 0);
    _nassert((size_t)from->data % 4 == 0);

    size_t i;

    for (i = 0; i < n; i++) {
        to->data[i] = from->data[i];
    }
}

static inline void qpDUNES_copyArray(real_t* restrict const to,
const real_t* const from, size_t n) {
    assert(to && from);
    _nassert((size_t)to % 4 == 0);
    _nassert((size_t)from % 4 == 0);

    size_t i;

    for (i = 0; i < n; i++) {
        to[i] = from[i];
    }
}

static inline real_t qpDUNES_fmax(real_t a, real_t b) {
    return (a > b) ? a : b;
}

static inline real_t qpDUNES_fmin(real_t a, real_t b) {
    return (a < b) ? a : b;
}

#endif  /* QPDUNES_UTILS_H */

