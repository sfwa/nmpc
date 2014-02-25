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
 *	\file src/utils.c
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */


#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "utils.h"


/* ----------------------------------------------
 * safe array offset routine, avoids NULL
 * pointer offsetting
 *
>>>>>>                                            */
/*inline real_t* offsetArray(	real_t* data,*/
const real_t* offsetArray(	const real_t* const data,
									int_t offset
									)
{
	return (data != 0) ? &(data[offset]) : 0;
}
/*<<< offsetArray */



void qpDUNES_copyRealArray( int_t n,
						 real_t* const to,
						 const real_t* const from
						 )
{
	int_t i;

	if ( from != 0 )
	{
		for( i=0; i<n; ++i )
			to[i] = from[i];
	}
}


sparsityType_t qpDUNES_detectMatrixSparsity(	const real_t* const M,
											int_t nRows,
											int_t nCols
											)
{
	sparsityType_t sparsityM;
	int_t i,j;

	if ( ( nRows < 1 ) || ( nCols < 1 ) || ( M == 0 ) )
		return (sparsityType_t)QPDUNES_OK;

	if ( nRows != nCols )
	{
		sparsityM = QPDUNES_DENSE;
		return (sparsityType_t)QPDUNES_OK;
	}

	/* check for sparsity */
	sparsityM = QPDUNES_DIAGONAL;

	for( i=0; i<nRows; ++i ) {	/* check if dense */
		for( j=0; j<i-1; ++j ) {	/* lower triangle */
			if ( fabs( M[i*nCols+j] ) > 1e-15 ) {	/* TODO: make threshold adjustable! */
				sparsityM = QPDUNES_DENSE;
				break;
			}
		}
		for( j=i+1; j<nCols; ++j ) {	/* upper triangle */
			if ( fabs( M[i*nCols+j] ) > 1e-15 ) {
				sparsityM = QPDUNES_DENSE;
				break;
			}
		}
	}

	/* check whether diagonal or identity */
	if ( sparsityM != QPDUNES_DENSE )
	{
		sparsityM = QPDUNES_IDENTITY;

		for( i=0; i<nRows; ++i ) {
			if ( fabs( M[i*nCols+i] - 1.0 ) > 1e-15 ) {
				sparsityM = QPDUNES_DIAGONAL;
				break;
			}
		}
	}

	return sparsityM;

}


return_t qpDUNES_updateMatrixData(	matrix_t* const to,
								const real_t* const from,
								int_t nRows,
								int_t nCols
								)
{
	int_t i;

	if ( from == 0 )
		return QPDUNES_OK;

	if ( to == 0 )
		return QPDUNES_ERR_INVALID_ARGUMENT;

	switch ( to->sparsityType )
	{
		case QPDUNES_DENSE:
			for( i=0; i<nRows*nCols; ++i )
				to->data[i] = from[i];
			break;

		case QPDUNES_DIAGONAL:
			for( i=0; i<nRows; ++i )
				to->data[i] = from[i*nCols+i];
			break;

		case QPDUNES_IDENTITY:
			break;

		default:
			return QPDUNES_ERR_UNKNOWN_MATRIX_SPARSITY_TYPE;
	}

	return QPDUNES_OK;
}


return_t qpDUNES_setupZeroMatrix(	int_t nRows,
								int_t nCols,
								matrix_t* to
								)
{
	int_t i;

	for( i=0; i<nRows*nCols; ++i )
		to->data[i] = 0.0;

	to->sparsityType = QPDUNES_ALLZEROS;

	return QPDUNES_OK;
}


return_t qpDUNES_setMatrixNull(	matrix_t* const matrix
								)
{
	matrix->sparsityType = QPDUNES_MATRIX_UNDEFINED;

	return QPDUNES_OK;
}


return_t qpDUNES_existsMatrix(	matrix_t* matrix
							)
{
	if( matrix->data == 0 ) {
		return (return_t)QPDUNES_FALSE;
	}
	else {
		return (return_t)QPDUNES_TRUE;
	}
}


return_t qpDUNES_existsVector(	vector_t* vector
							)
{
	if( vector->data == 0 ) {
		return (return_t)QPDUNES_FALSE;
	}
	else {
		return (return_t)QPDUNES_TRUE;
	}
}


return_t qpDUNES_setupIdentityMatrix(	matrix_t* to
									)
{
	to->sparsityType = QPDUNES_IDENTITY;

	return QPDUNES_OK;
}


return_t qpDUNES_setupScaledIdentityMatrix(	int_t nRows,
											real_t scalar,
											matrix_t* to
											)
{
	int_t i;

	qpDUNES_setupZeroMatrix( nRows,nRows,to );

	for( i=0; i<nRows; ++i )
		to->data[i*nRows+i] = scalar;

	to->sparsityType = QPDUNES_DIAGONAL;

	return QPDUNES_OK;
}



return_t qpDUNES_setupVector(	vector_t* const to,
							const real_t* const from,
							int_t n
							)
{
	return qpDUNES_updateVector( to,from,n );
}


return_t qpDUNES_updateVector(	vector_t* const to,
							const real_t* const from,
							int_t n
							)
{
	int_t i;

	if ( ( n < 1 ) || ( from == 0 ) )
		return QPDUNES_OK;

	if ( to == 0 )
		return QPDUNES_ERR_INVALID_ARGUMENT;

	/* copy data */
	for( i=0; i<n; ++i ) {
		to->data[i] = from[i];
	}

	return QPDUNES_OK;
}


return_t qpDUNES_updateSimpleBoundVector(	qpData_t* qpData,
										vector_t* const to,
										const real_t* const dBnd,
										const real_t* const xBnd,
										const real_t* const uBnd
										)
{
	uint_t i;

	if ( dBnd != 0 ) {
		for( i=0; i<_NX_+_NU_; ++i ) {
			to->data[i] = dBnd[i];
		}
	}
	else {
		if ( xBnd != 0 ) {
			for( i=0; i<_NX_; ++i ) {
				to->data[i] = xBnd[i];
			}
		}
		if ( uBnd != 0 ) {
			for( i=0; i<_NU_; ++i ) {
				to->data[_NX_+i] = uBnd[i];
			}
		}
	}

	return QPDUNES_OK;
}


return_t qpDUNES_updateConstraintVector( 	vector_t* const to,
										const real_t* const dBnd,
										int_t nD
										)
{
	int_t i;

	if ( dBnd != 0 ) {
		for( i=0; i<nD; ++i ) {
			to->data[i] = dBnd[i];
		}
	}
	else {
		return QPDUNES_ERR_INVALID_ARGUMENT;
	}

	return QPDUNES_OK;
}


return_t qpDUNES_setupZeroVector(	vector_t* const to,
								int_t n
								)
{
	int_t i;

	for( i=0; i<n; ++i )
		to->data[i] = 0.0;

	return QPDUNES_OK;
}


return_t qpDUNES_setupUniformVector(	vector_t* const to,
									real_t value,
									int_t n
									)
{
	int_t i;

	for( i=0; i<n; ++i ) {
		to->data[i] = value;
	}

	return QPDUNES_OK;
}


return_t qpDUNES_copyVector(	vector_t* const to,
							const vector_t* const from,
							int_t n
							)
{
	int_t ii;

	for( ii=0; ii<n; ++ii )
		to->data[ii] = from->data[ii];

	return QPDUNES_OK;
}


/** deep matrix copy */
return_t qpDUNES_copyMatrix(	matrix_t* const to,
							const matrix_t* const from,
							int_t dim0,
							int_t dim1
							)
{
	int_t ii;

	/** choose appropriate copy routine */
	switch( from->sparsityType )
	{
		case QPDUNES_DENSE		:
		case QPDUNES_SPARSE	:
			for( ii=0; ii < dim0*dim1; ++ii ) {
				to->data[ii] = from->data[ii];
			}
			to->sparsityType = QPDUNES_DENSE;
			break;
		case QPDUNES_DIAGONAL	:
			/* matrix diagonal is saved in first line */
			for( ii=0; ii < dim1; ++ii ) {
				to->data[ii] = from->data[ii];
			}
			to->sparsityType = QPDUNES_DIAGONAL;
			break;
		case QPDUNES_IDENTITY	:
			to->sparsityType = QPDUNES_IDENTITY;
			break;
		default				:
			return QPDUNES_ERR_UNKNOWN_MATRIX_SPARSITY_TYPE;
	}

	return QPDUNES_OK;
}
/*<<<< END OF qp42_copyMatrix */


/** ... */
return_t qpDUNES_makeMatrixDense( matrix_t* const M_ptr,
								int_t dim0,
								int_t dim1
								)
{
	int_t ii, jj;

	real_t* M = M_ptr->data; /* enable matrix access by preprocessor macro */

	switch( M_ptr->sparsityType )
	{
		case QPDUNES_DENSE		:
			break;

		case QPDUNES_SPARSE	:
			M_ptr->sparsityType = QPDUNES_DENSE;
			break;

		case QPDUNES_DIAGONAL	:
			/* matrix diagonal is saved in first line */
			for( ii=dim0-1; ii >= 0; --ii ) {	/* go through matrix back to front */
				for( jj=dim1-1; jj > ii; --jj ) {
					accM( ii,jj,dim1 ) = 0.;
				}
				accM( ii,ii,dim1 ) = accM( 0,ii,dim1 );
				for( jj=ii-1; jj >= 0; --jj ) {
					accM( ii,jj,dim1 ) = 0.;
				}
			}
			M_ptr->sparsityType = QPDUNES_DENSE;
			break;

		case QPDUNES_IDENTITY	:
			for( ii=dim0-1; ii >= 0; --ii ) {	/* go through matrix back to front */
				for( jj=dim1-1; jj > ii; --jj ) {
					accM( ii,jj,dim1 ) = 0.;
				}
				accM( ii,ii,dim1 ) = 1.;
				for( jj=ii-1; jj >= 0; --jj ) {
					accM( ii,jj,dim1 ) = 0.;
				}
			}
			M_ptr->sparsityType = QPDUNES_DENSE;
			break;

		default				:
			return QPDUNES_ERR_UNKNOWN_MATRIX_SPARSITY_TYPE;
	}

	return QPDUNES_OK;
}
/*<<<< END OF qp42_makeMatrixDense */


/** transpose matrix (deep copy) */
return_t qpDUNES_transposeMatrix(	matrix_t* const to,
								const matrix_t* const from,
								int_t dim0,
								int_t dim1
								)
{
	int_t ii,jj;

	/** choose appropriate copy routine */
	switch( from->sparsityType )
	{
		case QPDUNES_DENSE		:
			to->sparsityType = QPDUNES_DENSE;
			for( ii=0; ii < dim1; ++ii ) {	/* go by columns of from matrix */
				for( jj=0; jj < dim0; ++jj ) {
					to->data[ii*dim0+jj] = from->data[jj*dim1+ii];
				}
			}
			break;

		case QPDUNES_SPARSE	:
/*			qp42_printWarning( __FILE__, __LINE__, "Sparse tranposeMatrix not implemented. Copying densely instead." );*/
			to->sparsityType = QPDUNES_DENSE;
			for( ii=0; ii < dim1; ++ii ) {	/* go by columns of from matrix */
				for( jj=0; jj < dim0; ++jj ) {
					to->data[ii*dim0+jj] = from->data[jj*dim1+ii];
				}
			}
			break;

		case QPDUNES_DIAGONAL	:
			break;

		case QPDUNES_IDENTITY	:
			to->sparsityType = QPDUNES_IDENTITY;
			break;

		default				:
			return QPDUNES_ERR_UNKNOWN_MATRIX_SPARSITY_TYPE;
	}

	return QPDUNES_OK;
}
/*<<<< END OF qp42_transposeMatrix */



/** selftranspose a square matrix */
return_t qpDUNES_selftransposeMatrix(	matrix_t* const Mptr,
									int_t dim			/**< leading and secondary dimension of M */
									)
{
	int_t ii,jj;

	real_t swap;
	real_t* M = Mptr->data;

	/** choose appropriate copy routine */
	switch( Mptr->sparsityType )
	{
		case QPDUNES_DENSE		:
		case QPDUNES_SPARSE	:
			Mptr->sparsityType = QPDUNES_DENSE;
			for( ii=0; ii < dim; ++ii ) {	/* go by rows of untransposed M in lower triangle */
				for( jj=0; jj < ii; ++jj ) {
					swap = accM(ii,jj,dim);
					accM(ii,jj,dim) = accMT(ii,jj,dim);
					accMT(ii,jj,dim) = swap;
				}
			}
			break;

		case QPDUNES_DIAGONAL	:
		case QPDUNES_IDENTITY	:
			break;

		default				:
			return QPDUNES_ERR_UNKNOWN_MATRIX_SPARSITY_TYPE;
	}

	return QPDUNES_OK;
}
/*<<<< END OF qp42_selftransposeMatrix */



return_t qpDUNES_copyArray(	real_t* const to,
							const real_t* const from,
							int_t n
							)
{
	int_t ii;

	for( ii=0; ii<n; ++ii )
		to[ii] = from[ii];

	return QPDUNES_OK;
}


/* ----------------------------------------------
 * max routine
 >>>>>                                            */
int_t qpDUNES_max(	int_t a,
						int_t b )
{
	return (a > b) ? a : b;
}
/*<<< END OF qp42_max */


/* ----------------------------------------------
 * min routine
 >>>>>                                            */
int_t qpDUNES_min(	int_t a,
						int_t b )
{
	return (a < b) ? a : b;
}
/*<<< END OF qp42_min */


/* ----------------------------------------------
 * max routine
 >>>>>                                            */
real_t qpDUNES_fmax(	real_t a,
							real_t b )
{
	return (a > b) ? a : b;
}
/*<<< END OF qp42_fmax */


/* ----------------------------------------------
 * min routine
 >>>>>                                            */
real_t qpDUNES_fmin(	real_t a,
							real_t b )
{
	return (a < b) ? a : b;
}
/*<<< END OF qp42_fmin */


/* ----------------------------------------------
 * sign routine
 >>>>>                                            */
int_t qpDUNES_sign(	const qpData_t* const qpData,
					real_t a
					)
{
	return (a < -qpData->options.equalityTolerance) ? -1 : ( (a > qpData->options.equalityTolerance) ? 1 : 0 );
}
/*<<< END OF qp42_sign */



/*extern inline void qp42_assertOK(	return_t statusFlag,*/
void qpDUNES_assertOK(	return_t statusFlag,
							char* fileName,
							int_t lineNumber,
							char* errString
							)
{
	assert(statusFlag != QPDUNES_OK);
}

/*
 *	end of file
 */
