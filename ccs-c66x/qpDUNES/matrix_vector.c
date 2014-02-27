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
    return multiplyVectorMatrixVector((matrix_t*)H, (vector_t*)z, nV);
}


/* Matrix-vector product y = invH*z */
return_t multiplyInvHz(qpData_t* const qpData, z_vector_t* const res,
const vv_matrix_t* const cholH, const z_vector_t* const z, const size_t nV) {
    return multiplyInvMatrixVector(qpData, (vector_t*)res, (matrix_t*)cholH,
                                   (vector_t*)z, nV);
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

    for (i = 0; i < _NX_; i++) {
        res->data[i] = 0.0;
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
    for (j = 0; j < _NZ_; j++) {
        res->data[j] = 0.0;
    }

    for (i = 0; i < _NX_; i++) {
        for (j = 0; j < _NZ_; j++) {
            res->data[j] += accC(i, j) * y->data[i];
        }
    }

    return QPDUNES_OK;
}


/* Matrix times inverse matrix product res = A * Q^-1 */
return_t multiplyAInvQ(qpData_t* const qpData, xx_matrix_t* const res,
const xx_matrix_t* const C, const xx_matrix_t* const cholH) {
    assert(qpData && res && C && cholH && res->data && C->data &&
           cholH->data);
    _nassert((size_t)res->data % 4 == 0);
    _nassert((size_t)C->data % 4 == 0);
    _nassert((size_t)cholH->data % 4 == 0);

    size_t i, j;

    res->sparsityType = QPDUNES_DENSE;

    /* choose appropriate multiplication routine */
    switch (cholH->sparsityType) {
        /* cholH invalid types */
        case QPDUNES_DENSE:
        case QPDUNES_SPARSE:
        case QPDUNES_MATRIX_UNDEFINED:
        case QPDUNES_ALLZEROS:
            assert(0 && "Invalid cholH sparsity type");
            return QPDUNES_ERR_UNKNOWN_MATRIX_SPARSITY_TYPE;
        /* cholH diagonal */
        case QPDUNES_DIAGONAL:
            /* scale A part of C column-wise */
            for (i = 0; i < _NX_; i++) {
                for (j = 0; j < _NX_; j++) {
                    /* cholH is the actual matrix in diagonal case */
                    res->data[i * _NX_ + j] =
                        divide_f(accC(i, j), cholH->data[j]);
                }
            }
            break;
        /* cholH identity */
        case QPDUNES_IDENTITY:
            /* copy A block */
            for (i = 0; i < _NX_; i++) {
                for (j = 0; j < _NX_; j++) {
                    res->data[i * _NX_ + j] = accC(i, j);
                }
            }
            break;
    }

    return QPDUNES_OK;
}


/* Inverse matrix times identity matrix product res = Q^-1 * I */
return_t getInvQ(qpData_t* const qpData, xx_matrix_t* const res,
const vv_matrix_t* const cholH, size_t nV) {
    assert(qpData && res && cholH && res->data && cholH->data);
    _nassert((size_t)res->data % 4 == 0);
    _nassert((size_t)cholH->data % 4 == 0);

    switch (cholH->sparsityType){
        /* cholM1 dense */
        case QPDUNES_DENSE:
        case QPDUNES_SPARSE:
        case QPDUNES_MATRIX_UNDEFINED:
        case QPDUNES_ALLZEROS:
            assert(0 && "Invalid cholH sparsity type");
            return QPDUNES_ERR_UNKNOWN_MATRIX_SPARSITY_TYPE;
        /* cholM1 diagonal */
        case QPDUNES_DIAGONAL:
            res->sparsityType = QPDUNES_DIAGONAL;
            /* cholH is just save in first line; first _NX_ elements are Q part */
            return backsolveMatrixDiagonalIdentity(qpData, res->data,
                                                   cholH->data, nV);
        /* cholM1 identity */
        case QPDUNES_IDENTITY:
            res->sparsityType = QPDUNES_IDENTITY;
            return QPDUNES_OK;
    }
}

return_t factorizeH(qpData_t* const qpData, vv_matrix_t* const cholH,
const vv_matrix_t* const H, size_t nV) {
    return factorizePosDefMatrix(qpData, (matrix_t*)cholH, (matrix_t*)H, nV);
}


return_t addCInvHCT(qpData_t* const qpData, xx_matrix_t* const res,
const vv_matrix_t* const cholH, const xz_matrix_t* const C,
const d2_vector_t* const y, zx_matrix_t* const zxMatTmp) {
    /* TODO: summarize to one function */
    return addMultiplyMatrixInvMatrixMatrixT(qpData, res, cholH, C, y->data,
            zxMatTmp, &(qpData->xVecTmp), _NX_, _NZ_);
}


return_t addScaledVector(vector_t* restrict const res, real_t scalingFactor,
const vector_t* const update, size_t len) {
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
const vector_t* const x, real_t scalingFactor, const vector_t* const y,
size_t len) {
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


return_t addVectors(vector_t* restrict const res, const vector_t* const x,
const vector_t* const y, size_t len) {
    assert(res && x && y && res->data && x->data && y->data);
    _nassert((size_t)res->data % 4 == 0);
    _nassert((size_t)x->data % 4 == 0);
    _nassert((size_t)y->data % 4 == 0);

    size_t i;

    for (i = 0; i < len; i++) {
        res->data[i] = x->data[i] + y->data[i];
    }

    return QPDUNES_OK;
}


return_t addToVector(vector_t* restrict const res,
const vector_t* const update, size_t len) {
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


/* Compute a Cholesky factorization of M */
return_t factorizePosDefMatrix(qpData_t* const qpData, matrix_t* const cholM,
const matrix_t* const M, size_t dim0) {
    /* choose appropriate factorization routine */
    switch (M->sparsityType) {
        case QPDUNES_DENSE:
        case QPDUNES_SPARSE:
            cholM->sparsityType = QPDUNES_DENSE;
            return denseCholeskyFactorization(qpData, cholM, M, dim0);
        case QPDUNES_DIAGONAL:
            /*
            cholM in this case is defined to contain full diagonal matrix
            (not a factor)
            */
            return qpDUNES_copyMatrix(cholM, M, dim0, dim0);
        case QPDUNES_IDENTITY:
            cholM->sparsityType = QPDUNES_IDENTITY;
            return QPDUNES_OK;
        case QPDUNES_MATRIX_UNDEFINED:
        case QPDUNES_ALLZEROS:
            assert(0 && "Invalid M sparsity type");
            return QPDUNES_ERR_UNKNOWN_MATRIX_SPARSITY_TYPE;
    }
}


/* Compute a dense Cholesky factorization of M  */
return_t denseCholeskyFactorization(qpData_t* const qpData,
matrix_t* const cholM, const matrix_t* const M, size_t dim0) {
    assert(qpData && cholM && M && cholM->data && M->data);
    _nassert((size_t)cholM->data % 4 == 0);
    _nassert((size_t)M->data % 4 == 0);

    size_t i, j, k;
    real_t sum;

    /* go by columns */
    for (i = 0; i < dim0; i++) {
        /* write diagonal element: j == i */
        sum = M->data[i * dim0 + i];

        /* subtract squared forepart of this row */
        for (k = 0; k < i; k++) {
            sum -= cholM->data[i * dim0 + k] * cholM->data[i * dim0 + k];
        }

        if (sum > qpData->options.QPDUNES_ZERO) {
            cholM->data[i * dim0 + i] = sqrt_f(sum);
        } else {
            /* matrix not positive definite */
            return QPDUNES_ERR_DIVISION_BY_ZERO;
        }

        /* write remainder of i-th column */
        for (j = i + 1u; j < dim0; j++) {
            sum =  M->data[j * dim0 + i];

            /* subtract forepart of this row times forepart of ii-th row */
            for (k = 0; k < i; k++) {
                sum -= cholM->data[i * dim0 + k] * cholM->data[j * dim0 + k];
            }

            cholM->data[j * dim0 + i] =
                divide_f(sum, cholM->data[i * dim0 + i]);
        }
    }

    return QPDUNES_OK;
}


/* Backsolve for a dense L compute res for L*res = b and L^T*res = b */
return_t backsolveDenseL(qpData_t* const qpData, real_t* const res,
const real_t* const L, const real_t* const b, boolean_t transposed,
size_t n) {
    assert(qpData && res && L && b);
    assert(n);
    _nassert((size_t)res % 4 == 0);
    _nassert((size_t)L % 4 == 0);
    _nassert((size_t)b % 4 == 0);

    size_t i, j;
    real_t sum;

    /* Solve La = b, where L might be transposed. */
    if (transposed == QPDUNES_FALSE) {
        /* solve L*a = b */
        for (i = 0; i < n; i++) {
            sum = b[i];
            for (j = 0; j < i; j++) {
                sum -= accL(i, j, n) * res[j];
            }

            if (abs_f(accL(i, i, n)) >=
                    qpData->options.QPDUNES_ZERO * abs_f(sum)) {
                res[i] = divide_f(sum, accL(i, i, n));
            } else {
                return QPDUNES_ERR_DIVISION_BY_ZERO;
            }
        }
    } else {
        /* solve L^Ta = b */
        for (i = n - 1u; i != SIZE_MAX; i--) {
            sum = b[i];
            for (j = i + 1u; j < n; j++) {
                sum -= accL(j, i, n) * res[j];
            }

            if (abs_f(accL(i, i, n)) >=
                    qpData->options.QPDUNES_ZERO * abs_f(sum)) {
                res[i] = divide_f(sum, accL(i, i, n));
            } else {
                return QPDUNES_ERR_DIVISION_BY_ZERO;
            }
        }
    }

    return QPDUNES_OK;
}


/* Backsolve for a diagonal M compute res for M*res = b */
return_t backsolveDiagonal(qpData_t* const qpData, real_t* const res,
const real_t* const M, const real_t* const b, size_t n) {
    assert(qpData && res && M && b);
    assert(n);
    _nassert((size_t)res % 4 == 0);
    _nassert((size_t)M % 4 == 0);
    _nassert((size_t)b % 4 == 0);

    size_t i;

    /* Solve M*res = b */
    for (i = 0; i < n; i++) {
        res[i] = divide_f(b[i], accM(0, i, n));
    }

    return QPDUNES_OK;
}


/* Matrix backsolve for M1 diagonal, M2 identity compute res for M1*res = I */
return_t backsolveMatrixDiagonalIdentity(qpData_t* const qpData,
real_t* const res, const real_t* const M1, size_t dim0) {
    assert(qpData && res && M1);
    _nassert((size_t)res % 4 == 0);
    _nassert((size_t)M1 % 4 == 0);

    size_t i;

    /* backsolve on diagonal matrix: res for M1*res = I */
    for (i = 0; i < dim0; i++) {
        if (abs_f(M1[i]) >= qpData->options.QPDUNES_ZERO) {
            /* M1 is the actual matrix in diagonal case */
            res[i] = recip_f(M1[i]);
        } else {
            return QPDUNES_ERR_DIVISION_BY_ZERO;
        }
    }

    return QPDUNES_OK;
}


/* Matrix backsolve for L, M both dense compute res for L*res = M^T */
return_t backsolveMatrixTDenseDenseL(qpData_t* const qpData,
real_t* const res, const real_t* const L,
const real_t* const M, /* untransposed M */
real_t* const sums, /* memory for saving intermediate results (for speedup) */
boolean_t transposedL, size_t dim0, size_t dim1 /* dimensions of M */) {
    assert(qpData && res && L && M && sums);
    assert(dim1);
    _nassert((size_t)res % 4 == 0);
    _nassert((size_t)L % 4 == 0);
    _nassert((size_t)M % 4 == 0);
    _nassert((size_t)sums % 4 == 0);

    size_t i, j, k;

    /* Solve L*A = B^T, where L might be transposed. */
    if (transposedL == QPDUNES_FALSE) {
        /* solve L*A = B^T */
        for (i = 0; i < dim1; i++) { /* go by rows */
            for (k = 0; k < dim0; k++) {
                sums[k] = accM(k, i, dim0);
            }
            for (j = 0; j < i; j++) {
                for (k = 0; k < dim0; k++) {
                    sums[k] -= accL(i, j, dim1) * res[j * dim0 + k];
                }
            }
            for (k = 0; k < dim0; k++) {
                if (abs_f(accL(i, i, dim1)) >=
                        qpData->options.QPDUNES_ZERO * abs_f(sums[k])) {
                    res[i * dim0 + k] = divide_f(sums[k], accL(i, i, dim1));
                } else {
                    return QPDUNES_ERR_DIVISION_BY_ZERO;
                }
            }
        }
    } else {
        /* solve L^T*A = B^T */
        for (i = dim1 - 1u; i != SIZE_MAX; i--) { /* go by rows, bottom-up */
            for (k = 0; k < dim0; k++) {
                sums[k] = accM(k, i, dim0);
            }
            for (j= i + 1u; j < dim1; j++) {
                for (k = 0; k < dim0; k++) {
                    sums[k] -= accL(j, i, dim1) * res[j * dim0 + k];
                }
            }
            for(k = 0; k < dim0; k++) {
                if (abs_f(accL(i, i, dim1)) >=
                        qpData->options.QPDUNES_ZERO * abs_f(sums[k])) {
                    res[i * dim0 + k] = divide_f(sums[k], accL(i, i, dim1));
                } else {
                    return QPDUNES_ERR_DIVISION_BY_ZERO;
                }
            }
        }
    }

    return QPDUNES_OK;
}


/* Matrix backsolve for L diagonal and M dense compute res for L*res = M^T */
return_t backsolveMatrixTDiagonalDense(qpData_t* const qpData,
real_t* const res, const real_t* const M1,
const real_t* const M2, /* untransposed M2 */
size_t dim0, size_t dim1) { /* dimensions of M2 */
    assert(qpData && res && M1 && M2);
    _nassert((size_t)res % 4 == 0);
    _nassert((size_t)M1 % 4 == 0);
    _nassert((size_t)M2 % 4 == 0);

    size_t i, j;

    for (i = 0; i < dim1; i++) {
        for (j = 0; j < dim0; j++) {
            if (abs_f(M1[i]) >=
                    qpData->options.QPDUNES_ZERO * abs_f(M2[j * dim1 + i])) {
                /*
                M1 is the actual matrix in diagonal case; M2 is untransposed
                */
                res[i * dim0 + j] = divide_f(M2[j * dim1 + i], M1[i]);
            } else {
                return QPDUNES_ERR_DIVISION_BY_ZERO;
            }
        }
    }

    return QPDUNES_OK;
}


/* Generic matrix-vector product b = M*x */
return_t multiplyMatrixVector(vector_t* const res, const matrix_t* const M,
const vector_t* const x, size_t dim0, size_t dim1) {
    assert(res && M && x && res->data && M->data && x->data);

    /* choose appropriate multiplication routine */
    switch (M->sparsityType) {
        case QPDUNES_DENSE:
            return multiplyMatrixVectorDense(res->data, M->data, x->data,
                                             dim0, dim1);
        case QPDUNES_SPARSE:
            return multiplyMatrixVectorSparse(res->data, M->data, x->data,
                                              dim0, dim1);
        case QPDUNES_DIAGONAL:
            return multiplyMatrixVectorDiagonal(res->data, M->data, x->data,
                                                dim0);
        case QPDUNES_IDENTITY:
            /* just copy vector */
            return qpDUNES_copyVector(res, x, dim0);
        case QPDUNES_MATRIX_UNDEFINED:
        case QPDUNES_ALLZEROS:
            assert(0 && "Invalid M sparsity type");
            return QPDUNES_ERR_UNKNOWN_MATRIX_SPARSITY_TYPE;
    }
}


/*
Generic vector-matrix-vector product b = x'*M*x - M has to be square matrix
*/
real_t multiplyVectorMatrixVector(const matrix_t* const M,
const vector_t* const x, size_t dim0) {
    assert(M && x && M->data && x->data);

    /* choose appropriate multiplication routine */
    switch (M->sparsityType) {
        case QPDUNES_DENSE:
        case QPDUNES_SPARSE:
            return multiplyVectorMatrixVectorDense(M->data, x->data, dim0);
        case QPDUNES_DIAGONAL:
            return multiplyVectorMatrixVectorDiagonal(M->data, x->data, dim0);
        case QPDUNES_IDENTITY:
            /* just square vector */
            return scalarProd(x, x, dim0);
        case QPDUNES_MATRIX_UNDEFINED:
        case QPDUNES_ALLZEROS:
            assert(0 && "Invalid M sparsity type");
            return -1.0;
    }
}


/* Generic transposed matrix-vector product res = A.T*x */
return_t multiplyMatrixTVector(vector_t* const res, const matrix_t* const M,
const vector_t* const x, size_t dim0, size_t dim1) {
    assert(res && M && x && res->data && x->data);

    /* choose appropriate multiplication routine */
    switch (M->sparsityType) {
        case QPDUNES_DENSE:
            return multiplyMatrixTVectorDense(res->data, M->data, x->data,
                                              dim0, dim1);
        case QPDUNES_SPARSE:
            return multiplyMatrixTVectorSparse(res->data, M->data, x->data,
                                               dim0, dim1);
        case QPDUNES_DIAGONAL:
            return multiplyMatrixVectorDiagonal(res->data, M->data, x->data,
                                                dim0);
        case QPDUNES_IDENTITY:
            /* just copy vector */
            return qpDUNES_copyVector(res, x, dim0);
        case QPDUNES_MATRIX_UNDEFINED:
        case QPDUNES_ALLZEROS:
            assert(0 && "Invalid M sparsity type");
            return QPDUNES_ERR_UNKNOWN_MATRIX_SPARSITY_TYPE;
    }
}


/*
Matrix-vector product res = invM*x, using a Cholesky factorization H = L*L^T,
where L is a lower triangular matrix.

Solve L*L^T * res = x for res by
1) solving L*y = x for y
2) solving L^T*res = y for res

dim0 is the dimension of the symmetric matrix.
*/
return_t multiplyInvMatrixVector(qpData_t* const qpData, vector_t* const res,
const matrix_t* const cholH, const vector_t* const x, size_t dim0) {
    assert(qpData && res && cholH && x && res->data && cholH->data &&
           x->data);

    return_t statusFlag;

    /* choose appropriate multiplication routine */
    switch (cholH->sparsityType) {
        case QPDUNES_DENSE:
            /* first backsolve: L*y = z */
            statusFlag = backsolveDenseL(qpData, res->data, cholH->data,
                                         x->data, QPDUNES_FALSE, dim0);
            /* second backsolve: L^T*res = y */
            statusFlag += backsolveDenseL(qpData, res->data, cholH->data,
                                          res->data, QPDUNES_TRUE, dim0);
            return statusFlag;
        case QPDUNES_SPARSE:
            /* first backsolve: L*y = z */
            statusFlag = backsolveDenseL(qpData, res->data, cholH->data,
                                         x->data, QPDUNES_FALSE, dim0);
            /* second backsolve: L^T*res = y */
            statusFlag += backsolveDenseL(qpData, res->data, cholH->data,
                                          res->data, QPDUNES_TRUE, dim0);
            return statusFlag;
        case QPDUNES_DIAGONAL:
            /*
            cholH in this case contains full diagonal matrix (not a factor)
            */
            return backsolveDiagonal(qpData, res->data, cholH->data, x->data,
                                     dim0);
        case QPDUNES_IDENTITY:
            /* just copy vector */
            return qpDUNES_copyVector(res, x, dim0);
        case QPDUNES_MATRIX_UNDEFINED:
        case QPDUNES_ALLZEROS:
            assert(0 && "Invalid cholH sparsity type");
            return QPDUNES_ERR_UNKNOWN_MATRIX_SPARSITY_TYPE;
    }
}


/* Dense generic matrix-vector product b = A*x */
return_t multiplyMatrixVectorDense(real_t* const res, const real_t* const M,
const real_t* const x, size_t dim0, size_t dim1) {
    assert(res && M && x);
    _nassert((size_t)res % 4 == 0);
    _nassert((size_t)M % 4 == 0);
    _nassert((size_t)x % 4 == 0);

    size_t i, j;

    for (i = 0; i < dim0; i++) {
        res[i] = 0.0;
        for (j = 0; j < dim1; j++) {
            res[i] += accM(i, j, dim1) * x[j];
        }
    }

    return QPDUNES_OK;
}


/* Dense generic vector-matrix-vector product a = x'*M*x */
real_t multiplyVectorMatrixVectorDense(const real_t* const M,
const real_t* const x, size_t dim0) {
    assert(M && x);
    _nassert((size_t)M % 4 == 0);
    _nassert((size_t)x % 4 == 0);

    size_t i, j;
    real_t result = 0.0;

    for (i = 0; i < dim0; i++) {
        for (j = 0; j < dim0; j++) {
            result += x[i] * accM(i, j, dim0) * x[j];
        }
    }

    return result;
}


/*
Sparse generic matrix-vector product b = A*x

WARNING: sparse matrix-vector multiplication not yet implemented
*/
return_t multiplyMatrixVectorSparse(real_t* res, const real_t* const M,
const real_t* const x, size_t dim0, size_t dim1) {
    return multiplyMatrixVectorDense(res, M, x, dim0, dim1);
}


/* Generic diagonal matrix-vector product b = A*x */
return_t multiplyMatrixVectorDiagonal(real_t* const res,
const real_t* const M, const real_t* const x, size_t dim0) {
    assert(res && M && x);
    _nassert((size_t)res % 4 == 0);
    _nassert((size_t)M % 4 == 0);
    _nassert((size_t)x % 4 == 0);

    size_t i;

    /* multiply vector with diagonal matrix saved in first line */
    for (i = 0; i < dim0; i++) {
        res[i] = accM(0, i, dim0) * x[i];
    }

    return QPDUNES_OK;
}


/* Generic diagonal vector-matrix-vector product a = x'*M*x */
real_t multiplyVectorMatrixVectorDiagonal(const real_t* const M,
const real_t* const x, size_t dim0) {
    assert(M && x);
    _nassert((size_t)M % 4 == 0);
    _nassert((size_t)x % 4 == 0);

    size_t j;
    real_t result = 0.;

    /* multiply vector with diagonal matrix saved in first line */
    for (j = 0; j < dim0; j++) {
        result += accM(0, j, dim0) * x[j] * x[j];
    }

    return result;
}


/* Dense generic transposed matrix-vector product res = M.T*x */
return_t multiplyMatrixTVectorDense(real_t* const res, const real_t* const M,
const real_t* const x, size_t dim0, size_t dim1) {
    assert(res && M && x);
    _nassert((size_t)res % 4 == 0);
    _nassert((size_t)M % 4 == 0);
    _nassert((size_t)x % 4 == 0);

    size_t i, j;

    /* change multiplication order for more efficient memory access */
    for (j = 0; j < dim1; j++) {
        res[j] = 0.0;
    }
    for (i = 0; i < dim0; i++) {
        for (j = 0; j < dim1; j++) {
            res[j] += accM(i, j, dim1) * x[i];
        }
    }

    return QPDUNES_OK;
}


/*
Sparse generic transposed matrix-vector product b = A*x

WARNING: sparse matrix-vector multiplication not yet implemented
*/
return_t multiplyMatrixTVectorSparse(real_t* res, const real_t* const M,
const real_t* const x, size_t dim0, size_t dim1) {
    return multiplyMatrixTVectorDense(res, M, x, dim0, dim1);
}


/* M2 * M1^-1 * M2.T -- result gets added to res, not overwritten */
return_t addMultiplyMatrixInvMatrixMatrixT(qpData_t* const qpData,
matrix_t* const res, const matrix_t* const cholM1, const matrix_t* const M2,
const real_t* const y, /* vector containing non-zeros for columns of M2 to be eliminated */
matrix_t* const Ztmp, /* temporary matrix of shape dim1 x dim0 */
vector_t* const vecTmp,
size_t dim0, size_t dim1) {/* dimensions of M2 */
    size_t i, j, l;

    /* compute M1^-1/2 * M2.T */
    switch (cholM1->sparsityType) {
        case QPDUNES_DENSE:
        case QPDUNES_SPARSE:
            backsolveMatrixTDenseDenseL(qpData, Ztmp->data, cholM1->data,
                                        M2->data, vecTmp->data, QPDUNES_FALSE,
                                        dim0, dim1);
            break;
        case QPDUNES_DIAGONAL:
            /* computes full inverse times M2! */
            backsolveMatrixTDiagonalDense(qpData, Ztmp->data, cholM1->data,
                                          M2->data, dim0, dim1);
            break;
        case QPDUNES_IDENTITY:
            qpDUNES_transposeMatrix(Ztmp, M2, dim0, dim1);
            break;
        case QPDUNES_MATRIX_UNDEFINED:
        case QPDUNES_ALLZEROS:
            assert(0 && "Invalid cholM1 sparsity type");
            return QPDUNES_ERR_UNKNOWN_MATRIX_SPARSITY_TYPE;
    }

    qpDUNES_makeMatrixDense(res, dim0, dim0);

    if (cholM1->sparsityType != QPDUNES_DIAGONAL) {
        /* compute Z.T * Z as dyadic products */
        for (l = 0; l < dim1; l++) {
            /*
            only add columns of variables with inactive lower and upper bounds
            */
            if ((y[2u * l] <= qpData->options.equalityTolerance) &&
                    (y[2u * l + 1u] <= qpData->options.equalityTolerance)) {
                for (i = 0; i < dim0; i++) {
                    for (j = 0; j < dim0; j++) {
                        /* since M2 is dense, so is Z */
                        res->data[i * dim0 + j] += Ztmp->data[l * dim0 + i] *
                                                   Ztmp->data[l * dim0 + j];
                    }
                }
            }
        }
    } else { /* diagonal H */
        /*
        Z already contains H^-1 * M2^T, therefore only multiplication with M2
        from left is needed

        compute M2 * Z as dyadic products
        */
        for (l = 0; l < dim1; l++) {
            /*
            only add columns of variables with inactive upper and lower bounds
            */
            if ((y[2u * l] <= qpData->options.equalityTolerance) &&
                    (y[2u * l + 1u] <= qpData->options.equalityTolerance)) {
                for (i = 0; i < dim0; i++) {
                    for (j = 0; j < dim0; j++) {
                        /* since M2 is dense, so is Z */
                        res->data[i * dim0 + j] += M2->data[i * dim1 + l] *
                                                   Ztmp->data[l * dim0 + j];
                    }
                }
            } /* end of dyadic addend */
        }
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
    real_t res = 0.0;

    for (i = 0; i < len; i++) {
        res += x->data[i] * y->data[i];
    }

    return res;
}


real_t vectorNorm(const vector_t* const vec, size_t len) {
    return sqrt_f(scalarProd(vec, vec, len));
}
