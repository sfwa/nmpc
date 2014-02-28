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


sparsityType_t qpDUNES_detectMatrixSparsity(const real_t* const M,
size_t nRows, size_t nCols) {
    sparsityType_t sparsityM;
    size_t i, j;

    if (!nRows || !nCols || !M) {
        return (sparsityType_t)QPDUNES_OK;
    }

    if (nRows != nCols) {
        sparsityM = QPDUNES_DENSE;
    }

    /* check for sparsity */
    sparsityM = QPDUNES_DIAGONAL;

    for (i = 0; i < nRows && sparsityM != QPDUNES_DENSE; i++) { /* check if dense */
        for (j = 0; i && sparsityM != QPDUNES_DENSE && j < i - 1u; j++) {   /* lower triangle */
            if (abs_f(M[i * nCols + j]) > 1e-15f) { /* TODO: make threshold adjustable! */
                sparsityM = QPDUNES_DENSE;
            }
        }
        for (j = i + 1u; sparsityM != QPDUNES_DENSE && j < nCols; j++) {    /* upper triangle */
            if (abs_f(M[i * nCols + j]) > 1e-15f) {
                sparsityM = QPDUNES_DENSE;
            }
        }
    }

    /* check whether diagonal or identity */
    if (sparsityM != QPDUNES_DENSE) {
        sparsityM = QPDUNES_IDENTITY;

        for (i = 0; i < nRows; i++) {
            if (abs_f(M[i * nCols + i] - 1.0f) > 1e-15f) {
                sparsityM = QPDUNES_DIAGONAL;
                break;
            }
        }
    }

    return sparsityM;
}


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


return_t qpDUNES_setupZeroMatrix(size_t nRows, size_t nCols, matrix_t* to) {
    assert(to && to->data);
    _nassert((size_t)to->data % 4 == 0);

    size_t i;

    for (i = 0; i < nRows * nCols; i++) {
        to->data[i] = 0.0;
    }

    to->sparsityType = QPDUNES_ALLZEROS;

    return QPDUNES_OK;
}


return_t qpDUNES_setMatrixNull(matrix_t* const matrix) {
    assert(matrix);

    matrix->sparsityType = QPDUNES_MATRIX_UNDEFINED;

    return QPDUNES_OK;
}


return_t qpDUNES_setupScaledIdentityMatrix(size_t nRows, real_t scalar,
matrix_t* to) {
    assert(to && to->data);
    _nassert((size_t)to->data % 4 == 0);

    size_t i;

    qpDUNES_setupZeroMatrix(nRows, nRows, to);
    for (i = 0; i < nRows; i++) {
        to->data[i * nRows + i] = scalar;
    }

    to->sparsityType = QPDUNES_DIAGONAL;

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


/** deep matrix copy */
return_t qpDUNES_copyMatrix(matrix_t* const to, const matrix_t* const from,
size_t dim0, size_t dim1) {
    assert(to && from && to->data && from->data);
    _nassert((size_t)to->data % 4 == 0);
    _nassert((size_t)from->data % 4 == 0);

    size_t i;

    /* choose appropriate copy routine */
    switch (from->sparsityType) {
        case QPDUNES_DENSE:
        case QPDUNES_SPARSE:
            for (i = 0; i < dim0 * dim1; i++) {
                to->data[i] = from->data[i];
            }
            to->sparsityType = QPDUNES_DENSE;
            break;
        case QPDUNES_DIAGONAL:
            /* matrix diagonal is saved in first line */
            for (i = 0; i < dim1; i++) {
                to->data[i] = from->data[i];
            }
            to->sparsityType = QPDUNES_DIAGONAL;
            break;
        case QPDUNES_IDENTITY:
            to->sparsityType = QPDUNES_IDENTITY;
            break;
        case QPDUNES_MATRIX_UNDEFINED:
        case QPDUNES_ALLZEROS:
            assert(0 && "Invalid sparsity type");
            return QPDUNES_ERR_UNKNOWN_MATRIX_SPARSITY_TYPE;
    }

    return QPDUNES_OK;
}


return_t qpDUNES_makeMatrixDense(matrix_t* const M_ptr, size_t dim0,
size_t dim1) {
    assert(M_ptr && M_ptr->data);
    assert(dim0);
    assert(dim1);

    ptrdiff_t i, j;
    /* enable matrix access by preprocessor macro */
    real_t* restrict M = M_ptr->data;

    switch (M_ptr->sparsityType) {
        case QPDUNES_DENSE:
            break;
        case QPDUNES_SPARSE :
            M_ptr->sparsityType = QPDUNES_DENSE;
            break;
        case QPDUNES_DIAGONAL:
            /*
            matrix diagonal is saved in first line; go through matrix back to
            front
            */
            for (i = (ptrdiff_t)dim0 - 1; i >= 0; i--) {
                for (j = (ptrdiff_t)dim1 - 1; j > i; j--) {
                    accM((size_t)i, (size_t)j, dim1) = 0.0;
                }
                accM((size_t)i, (size_t)i, dim1) = accM(0, (size_t)i, dim1);
                for (j = i - 1; j >= 0; j--) {
                    accM((size_t)i, (size_t)j, dim1) = 0.0;
                }
            }
            M_ptr->sparsityType = QPDUNES_DENSE;
            break;
        case QPDUNES_IDENTITY:
            for (i = (ptrdiff_t)dim0 - 1; i >= 0; i--) {
                for (j = (ptrdiff_t)dim1 - 1; j > i; j--) {
                    accM((size_t)i, (size_t)j, dim1) = 0.0;
                }
                accM((size_t)i, (size_t)i, dim1) = 1.0;
                for (j = i - 1; j >= 0; j--) {
                    accM((size_t)i, (size_t)j, dim1) = 0.0;
                }
            }
            M_ptr->sparsityType = QPDUNES_DENSE;
            break;
        case QPDUNES_MATRIX_UNDEFINED:
        case QPDUNES_ALLZEROS:
            assert(0 && "Invalid sparsity type");
            return QPDUNES_ERR_UNKNOWN_MATRIX_SPARSITY_TYPE;
    }

    return QPDUNES_OK;
}


/* transpose matrix (deep copy) */
return_t qpDUNES_transposeMatrix(matrix_t* const to,
const matrix_t* const from, size_t dim0, size_t dim1) {
    assert(to && from && to->data && from->data);
    _nassert((size_t)to->data % 4 == 0);
    _nassert((size_t)from->data % 4 == 0);

    size_t i, j;

    /* choose appropriate copy routine */
    switch (from->sparsityType) {
        case QPDUNES_DENSE:
            to->sparsityType = QPDUNES_DENSE;
            for (i = 0; i < dim1; i++) {  /* go by columns of from matrix */
                for (j = 0; j < dim0; j++) {
                    to->data[i * dim0 + j] = from->data[j * dim1 + i];
                }
            }
            break;
        case QPDUNES_SPARSE:
            to->sparsityType = QPDUNES_DENSE;
            for (i = 0; i < dim1; i++) {  /* go by columns of from matrix */
                for (j = 0; j < dim0; j++) {
                    to->data[i * dim0 + j] = from->data[j * dim1 + i];
                }
            }
            break;
        case QPDUNES_DIAGONAL:
            break;
        case QPDUNES_IDENTITY:
            to->sparsityType = QPDUNES_IDENTITY;
            break;
        case QPDUNES_MATRIX_UNDEFINED:
        case QPDUNES_ALLZEROS:
            assert(0 && "Invalid sparsity type");
            return QPDUNES_ERR_UNKNOWN_MATRIX_SPARSITY_TYPE;
    }

    return QPDUNES_OK;
}
