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


#include <stddef.h>
#include <stdint.h>
#include <assert.h>

#include "matrix_vector.h"
#include "../c66math.h"

/* Matrix-vector product y = z'*H*z */
real_t multiplyzHz(const vv_matrix_t* const H, const z_vector_t* const z,
const size_t nV) {
    assert(H && z);
    _nassert((size_t)H->data % 4 == 0);
    _nassert((size_t)z->data % 4 == 0);

    size_t j;
    real_t result = 0.0f;

    /*
    Multiply vector with diagonal matrix saved in first line. nV will either
    be 12 (_NX_ -- for the last interval) or 15 (_NZ_ -- for all others), but
    the actual size of the matrices is the same for all intervals so just
    take it as 15.
    */
    if (nV == 15u) {
        #pragma MUST_ITERATE(15, 15)
        for (j = 0; j < 15u; j++) {
            result += H->data[j] * z->data[j] * z->data[j];
        }
    } else {
        #pragma MUST_ITERATE(12, 12)
        for (j = 0; j < 12u; j++) {
            result += H->data[j] * z->data[j] * z->data[j];
        }
    }

    return result;
}

/*
Matrix-vector product y = invH*z, using a Cholesky factorization H = L*L^T,
where L is a lower triangular matrix.

Solve L*L^T * z = x for z by
1) solving L*y = x for y
2) solving L^T*z = y for z

nV is the dimension of the symmetric matrix.
*/
return_t multiplyInvHz(qpData_t* const qpData, z_vector_t* const res,
const vv_matrix_t* const cholH, const z_vector_t* const z, const size_t nV) {
    assert(qpData && res && cholH && z);
    assert(nV);
    _nassert((size_t)res->data % 4 == 0);
    _nassert((size_t)cholH->data % 4 == 0);
    _nassert((size_t)z->data % 4 == 0);

    size_t j;

    /*
    Solve H*res = z -- H elements are always reciprocal. The last interval has
    nV = 12 instead of nV = 15.
    */
    if (nV == 15u) {
        #pragma MUST_ITERATE(15, 15)
        for (j = 0; j < 15u; j++) {
            res->data[j] = cholH->data[j] * z->data[j];
        }
    } else {
        #pragma MUST_ITERATE(12, 12)
        for (j = 0; j < 12u; j++) {
            res->data[j] = cholH->data[j] * z->data[j];
        }
    }

    return QPDUNES_OK;
}

/* Matrix-vector product res = C*z */
return_t multiplyCz(qpData_t* const qpData, x_vector_t* const res,
const xz_matrix_t* const C, const z_vector_t* const z) {
    assert(qpData && res && C && z && res->data && C->data && z->data);
    _nassert((size_t)res->data % 4 == 0);
    _nassert((size_t)C->data % 4 == 0);
    _nassert((size_t)z->data % 4 == 0);

    /* only dense multiplication */
    size_t i, j;

    #pragma MUST_ITERATE(_NX_, _NX_)
    for (i = 0; i < _NX_; i++) {
        res->data[i] = 0.0;
        #pragma MUST_ITERATE(_NZ_, _NZ_)
        for (j = 0; j < _NZ_; j++) {
            res->data[i] += accC(i, j) * z->data[j];
        }
    }

    return QPDUNES_OK;
}

/* Matrix-vector product z = C.T*y */
return_t multiplyCTy(qpData_t* const qpData, z_vector_t* const res,
const xz_matrix_t* const C, const x_vector_t* const y) {
    assert(qpData && res && C && y && res->data && C->data && y->data);
    _nassert((size_t)res->data % 4 == 0);
    _nassert((size_t)C->data % 4 == 0);
    _nassert((size_t)y->data % 4 == 0);

    /* only dense multiplication */
    size_t i, j;

    /* change multiplication order for more efficient memory access */
    #pragma MUST_ITERATE(_NZ_, _NZ_)
    for (j = 0; j < _NZ_; j++) {
        res->data[j] = 0.0;
    }

    #pragma MUST_ITERATE(_NX_, _NX_)
    for (i = 0; i < _NX_; i++) {
        #pragma MUST_ITERATE(_NZ_, _NZ_)
        for (j = 0; j < _NZ_; j++) {
            res->data[j] += accC(i, j) * y->data[i];
        }
    }

    return QPDUNES_OK;
}

/* Matrix times inverse matrix product res = A * Q^-1 */
return_t multiplyAInvQ(qpData_t* const qpData,
xx_matrix_t* restrict const res, const xx_matrix_t* const C,
const vv_matrix_t* const cholH) {
    assert(qpData && res && C && cholH && res->data && C->data &&
           cholH->data);
    assert(cholH->sparsityType == QPDUNES_DIAGONAL);
    _nassert((size_t)res->data % 4 == 0);
    _nassert((size_t)C->data % 4 == 0);
    _nassert((size_t)cholH->data % 4 == 0);

    size_t i, j;

    res->sparsityType = QPDUNES_DENSE;

    /* scale A part of C column-wise */
    #pragma MUST_ITERATE(_NX_, _NX_)
    for (i = 0; i < _NX_; i++) {
        #pragma MUST_ITERATE(_NX_, _NX_)
        for (j = 0; j < _NX_; j++) {
            /*
            cholH is the actual matrix in diagonal case -- elements
            are stored as reciprocal
            */
            res->data[i * _NX_ + j] = accC(i, j) * cholH->data[j];
        }
    }

    return QPDUNES_OK;
}

/* Inverse matrix times identity matrix product res = Q^-1 * I */
return_t getInvQ(qpData_t* const qpData, xx_matrix_t* const res,
const vv_matrix_t* const cholH, size_t nV) {
    assert(qpData && res && cholH && res->data && cholH->data);
    assert(cholH->sparsityType == QPDUNES_DIAGONAL);
    _nassert((size_t)res->data % 4 == 0);
    _nassert((size_t)cholH->data % 4 == 0);

    size_t j;

    /*
    Backsolve on diagonal matrix: res for cholH*res = I. All elements of cholH
    are reciprocal.
    */
    if (nV == 15u) {
        #pragma MUST_ITERATE(15, 15)
        for (j = 0; j < 15u; j++) {
            res->data[j] = cholH->data[j];
        }
    } else {
        #pragma MUST_ITERATE(12, 12)
        for (j = 0; j < 12u; j++) {
            res->data[j] = cholH->data[j];
        }
    }

    return QPDUNES_OK;
}

/* M2 * M1^-1 * M2.T -- result gets added to res, not overwritten */
return_t addCInvHCT(qpData_t* const qpData, xx_matrix_t* const restrict res,
const vv_matrix_t* const restrict cholM1, const xz_matrix_t* const restrict M2,
const d2_vector_t* const y, /* vector containing non-zeros for columns of M2 to be eliminated */
zx_matrix_t* const restrict zxMatTmp) { /* temporary matrix of shape dim1 x dim0 */
    real_t* restrict yd = y->data;
    size_t i, j, l;
    return_t result = QPDUNES_OK;

    assert(cholM1->sparsityType == QPDUNES_DIAGONAL);

    qpDUNES_makeMatrixDense(res, _NX_, _NX_);

    /* compute M1^-1/2 * M2.T */

    /*
    Z already contains H^-1 * M2^T, therefore only multiplication with M2
    from left is needed

    compute M2 * Z as dyadic products
    */
    for (l = 0; l < _NZ_; l++) {
        for (j = 0; j < _NX_; j++) {
            /*
            M1 is the actual matrix in diagonal case; M2 is untransposed
            All elements of M1 are reciprocal
            */
            zxMatTmp->data[l * _NX_ + j] =
                M2->data[j * _NZ_ + l] * cholM1->data[l];

            if (abs_f(cholM1->data[l]) < qpData->options.QPDUNES_ZERO *
                                         abs_f(M2->data[j * _NZ_ + l])) {
                result = QPDUNES_ERR_DIVISION_BY_ZERO;
            }
        }

        /*
        only add columns of variables with inactive upper and lower bounds
        */
        if ((yd[2u * l] <= qpData->options.equalityTolerance) &&
                (yd[2u * l + 1u] <= qpData->options.equalityTolerance)) {
            #pragma MUST_ITERATE(_NX_, _NX_)
            for (i = 0; i < _NX_; i++) {
                #pragma MUST_ITERATE(_NX_, _NX_)
                for (j = 0; j < _NX_; j++) {
                    /* since M2 is dense, so is Z */
                    res->data[i * _NX_ + j] +=
                        M2->data[i * _NZ_ + l] *
                        zxMatTmp->data[l * _NX_ + j];
                }
            }
        } /* end of dyadic addend */
    }

    return result;
}

return_t addScaledVector(vector_t* restrict const res, real_t scalingFactor,
const vector_t* restrict const update, size_t len) {
    assert(res && update && res->data && update->data);
    _nassert((size_t)res->data % 4 == 0);
    _nassert((size_t)update->data % 4 == 0);

    size_t i;

    for (i = 0; i < len; i++) {
        res->data[i] += scalingFactor * update->data[i];
    }

    return QPDUNES_OK;
}

/* res = x + a*y  */
return_t addVectorScaledVector(vector_t* restrict const res,
const vector_t* restrict const x, real_t scalingFactor,
const vector_t* restrict const y, size_t len) {
    assert(res && x && y && res->data && x->data && y->data);
    _nassert((size_t)res->data % 4 == 0);
    _nassert((size_t)x->data % 4 == 0);
    _nassert((size_t)y->data % 4 == 0);

    size_t i;

    for (i = 0; i < len; i++) {
        res->data[i] = x->data[i] + scalingFactor * y->data[i];
    }

    return QPDUNES_OK;
}

return_t addToVector(vector_t* restrict const res,
const vector_t* restrict const update, size_t len) {
    assert(res && update && res->data && update->data);
    _nassert((size_t)res->data % 4 == 0);
    _nassert((size_t)update->data % 4 == 0);

    size_t i;

    for (i = 0; i < len; i++) {
        res->data[i] += update->data[i];
    }

    return QPDUNES_OK;
}

return_t subtractFromVector(vector_t* restrict const res,
const vector_t* const update, size_t len) {
    assert(res && update && res->data && update->data);
    _nassert((size_t)res->data % 4 == 0);
    _nassert((size_t)update->data % 4 == 0);

    size_t i;

    for (i = 0; i < len; i++) {
        res->data[i] -= update->data[i];
    }

    return QPDUNES_OK;
}

return_t negateVector(vector_t* const res, size_t len) {
    assert(res && res->data);
    _nassert((size_t)res->data % 4 == 0);

    size_t i;

    for (i = 0; i < len; i++) {
        res->data[i] = -res->data[i];
    }

    return QPDUNES_OK;
}

/* Low level scalar product */
real_t scalarProd(const vector_t* const x, const vector_t* const y,
size_t len) {
    assert(x && y);
    _nassert((size_t)x % 4 == 0);
    _nassert((size_t)y % 4 == 0);

    size_t i;
    real_t res = 0.0f;

    for (i = 0; i < len; i++) {
        res += x->data[i] * y->data[i];
    }

    return res;
}

real_t vectorNorm(const vector_t* const vec, size_t len) {
    return sqrt_f(scalarProd(vec, vec, len));
}
