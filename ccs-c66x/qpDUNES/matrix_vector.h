/*
 *	This file is part of qpDUNES.
 *
 *	qpDUNES -- A DUal NEwton Strategy for convex quadratic programming.
 *	Copyright (C) 2012 by Janick Frasch, Hans Joachim Ferreau et al.
 *	All rights reserved.
 *
 *	qpDUNES is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	qpDUNES is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with qpDUNES; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


/**
 *	\file include/qp/matrix_vector.h
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */


#ifndef QP42_MATRIX_VECTOR_H
#define QP42_MATRIX_VECTOR_H


#include "types.h"
#include "utils.h"
#include <math.h>
#include <stddef.h>


real_t multiplyzHz(const vv_matrix_t* const H, const z_vector_t* const z,
const size_t nV);

return_t multiplyInvHz(qpData_t* const qpData, z_vector_t* const res,
const vv_matrix_t* const cholH, const z_vector_t* const z, const size_t nV);

return_t multiplyCz(qpData_t* const qpData, x_vector_t* const res,
const xz_matrix_t* const C, const z_vector_t* const z);

return_t multiplyCTy(qpData_t* const qpData, z_vector_t* const res,
const xz_matrix_t* const C, const x_vector_t* const y);

/** Inverse matrix times matrix product res = Q^-1 * A */
return_t multiplyAInvQ(qpData_t* const qpData,
xx_matrix_t* restrict const res, const xx_matrix_t* const C,
const vv_matrix_t* const cholH);

return_t getInvQ(qpData_t* const qpData, xx_matrix_t* restrict const res,
const vv_matrix_t* restrict const cholH,
const d2_vector_t* restrict const y);

return_t addScaledVector(xn_vector_t* const res, real_t scalingFactor,
const xn_vector_t* restrict const update, size_t len);

return_t addVectorScaledVector(vector_t* restrict const res,
const vector_t* restrict const x, real_t scalingFactor,
const vector_t* restrict const y, size_t len);

/** Low-level scalar product */
real_t scalarProd(const vector_t* const x, const vector_t* const y,
size_t len);

return_t addToVector(vector_t* restrict const res,
const vector_t* restrict const update, size_t len);

return_t subtractFromVector(vector_t* const res, const vector_t* const update,
size_t len);

return_t negateVector(vector_t* const res, size_t len);

real_t vectorNorm(const vector_t* const vec, size_t len);

return_t addCInvHCT(qpData_t* const qpData, xx_matrix_t* const restrict res,
const vv_matrix_t* const restrict cholH, const xz_matrix_t* const restrict C,
const d2_vector_t* const y, /* vector containing non-zeros for columns of C to be eliminated */
zx_matrix_t* const restrict zxMatTmp);

#endif	/* QP42_MATRIX_VECTOR_H */
