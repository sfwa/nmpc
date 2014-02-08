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
 *	License along with qp42; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


/**
 *	\file include/qp/setup_qp.h
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */


#ifndef QP42_SETUP_QP_H
#define QP42_SETUP_QP_H


#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "types.h"
#include "matrix_vector.h"
#include "utils.h"


void qpDUNES_indicateDataChange(	qpData_t* const qpData
									);


return_t qpDUNES_updateData(	qpData_t* const qpData,
								const real_t* const H_,
								const real_t* const g_,
								const real_t* const C_,
								const real_t* const c_,
								const real_t* const zLow_,
								const real_t* const zUpp_,
								const real_t* const D_,
								const real_t* const dLow_,
								const real_t* const dUpp_
								);


return_t qpDUNES_setupRegularInterval(	qpData_t* const qpData,
									interval_t* interval,
									const real_t* const H_,
									const real_t* const Q_,
									const real_t* const R_,
									const real_t* const S_,
									const real_t* const g_,
									const real_t* const C_,
									const real_t* const A_,
									const real_t* const B_,
									const real_t* const c_,
									const real_t* const zLow_,
									const real_t* const zUpp_,
									const real_t* const xLow_,
									const real_t* const xUpp_,
									const real_t* const uLow_,
									const real_t* const uUpp_,
									const real_t* const D_,
									const real_t* const dLow_,
									const real_t* const dUpp_
									);


return_t qpDUNES_setupFinalInterval(	qpData_t* const qpData,
									interval_t* interval,
 									const real_t* const H_,
 									const real_t* const g_,
									const real_t* const zLow_,
									const real_t* const zUpp_,
									const real_t* const D_,
									const real_t* const dLow_,
									const real_t* const dUpp_
									);

return_t qpDUNES_updateIntervalData(	qpData_t* const qpData,
										interval_t* interval,
										const real_t* const H_,
										const real_t* const g_,
										const real_t* const C_,
										const real_t* const c_,
										const real_t* const zLow_,
										const real_t* const zUpp_,
										const real_t* const D_,
										const real_t* const dLow_,
										const real_t* const dUpp_,
										vv_matrix_t* const cholH
										);


/* ----------------------------------------------
 * first set up of local QP
 *
 >>>>>>                                           */
return_t qpDUNES_setupAllLocalQPs(	qpData_t* const qpData,
								boolean_t isLTI
								);


return_t qpDUNES_setupStageQP(	qpData_t* const qpData,
								interval_t* const interval,
								boolean_t copyCholH
								);


return_t qpDUNES_shiftIntervals(	qpData_t* const qpData
								);

return_t qpDUNES_shiftLambda(	qpData_t* const qpData
							);


qpOptions_t qpDUNES_setupDefaultOptions(	/*qpData_t* const qpData*/
										);

return_t qpDUNES_setupLog(	qpData_t* const qpData
						);


#endif	/* QP42_SETUP_QP_H */


/*
 *	end of file
 */
