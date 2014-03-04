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
 *  License along with qpDUNES; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include <stdint.h>
#include <stddef.h>
#include <math.h>
#include <assert.h>

#include "utils.h"
#include "../c66math.h"


return_t qpDUNES_updateMatrixData(matrix_t* const to,
const real_t* const from, size_t nRows, size_t nCols) {
    size_t i;

    if (!from) {
        return QPDUNES_OK;
    }

    assert(from && to && to->data);
    _nassert((size_t)to->data % 4 == 0);
    _nassert((size_t)from % 4 == 0);

    switch (to->sparsityType) {
        case QPDUNES_DENSE:
            for (i = 0; i < nRows * nCols; i++) {
                to->data[i] = from[i];
            }
            break;
        case QPDUNES_DIAGONAL:
            for (i = 0; i < nRows; i++) {
                to->data[i] = from[i * nCols + i];
            }
            break;
        case QPDUNES_IDENTITY:
            break;
        case QPDUNES_MATRIX_UNDEFINED:
        case QPDUNES_SPARSE:
        case QPDUNES_ALLZEROS:
            assert(0 && "Invalid sparsity type");
            return QPDUNES_ERR_UNKNOWN_MATRIX_SPARSITY_TYPE;
    }

    return QPDUNES_OK;
}


return_t qpDUNES_setMatrixNull(matrix_t* const matrix) {
    assert(matrix);

    matrix->sparsityType = QPDUNES_MATRIX_UNDEFINED;

    return QPDUNES_OK;
}


return_t qpDUNES_setupVector(vector_t* const to, const real_t* const from,
size_t n) {
    return qpDUNES_updateVector(to, from, n);
}


return_t qpDUNES_updateVector(vector_t* const to, const real_t* const from,
size_t n) {
    size_t i;

    if (!n || !from) {
        return QPDUNES_OK;
    }

    assert(from && to && to->data);
    _nassert((size_t)to->data % 4 == 0);
    _nassert((size_t)from % 4 == 0);

    /* copy data */
    for (i = 0; i < n; i++) {
        to->data[i] = from[i];
    }

    return QPDUNES_OK;
}


return_t qpDUNES_updateSimpleBoundVector(qpData_t* qpData, vector_t* const to,
const real_t* const dBnd, const real_t* const xBnd,
const real_t* const uBnd) {
    assert(qpData && to && to->data);
    _nassert((size_t)to->data % 4 == 0);
    _nassert((size_t)dBnd % 4 == 0);
    _nassert((size_t)xBnd % 4 == 0);
    _nassert((size_t)uBnd % 4 == 0);

    size_t i;

    if (dBnd) {
        for (i = 0; i < _NX_ + _NU_; i++) {
            to->data[i] = dBnd[i];
        }
    } else {
        if (xBnd) {
            for (i = 0; i < _NX_; i++) {
                to->data[i] = xBnd[i];
            }
        }
        if (uBnd) {
            for (i = 0; i < _NU_; i++) {
                to->data[_NX_ + i] = uBnd[i];
            }
        }
    }

    return QPDUNES_OK;
}
