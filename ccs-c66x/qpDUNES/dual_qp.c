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

#include "dual_qp.h"
#include "../c66math.h"

#include <stdio.h>

return_t qpDUNES_prepare(qpData_t* const qpData) {
    size_t ii;

    return_t statusFlag = QPDUNES_OK; /* generic status flag */
    interval_t* restrict interval;
    itLog_t* itLogPtr = &(qpData->log.itLog[0]);

    /** (1) todo: initialize local active sets (at least when using qpOASES) with initial guess from previous iteration */

    /** (2) solve local QP problems for initial guess of lambda */
    /* resolve initial QPs for possibly changed bounds (initial value embedding) */
    for (ii = 0; ii < _NI_ + 1; ++ii) {
        interval = qpData->intervals[ii];

        /* clip solution: */
        /* TODO: already clip all QPs except for the first one (initial value embedding); but take care for MHE!!!*/
        statusFlag = directQpSolver_doStep( qpData,
                                            interval,
                                            &(interval->qpSolverClipping.dz), 1,
                                            &(interval->qpSolverClipping.zUnconstrained),
                                            &(interval->z),
                                            &(interval->y),
                                            &(interval->q),
                                            &(interval->p)
                                            );
    }

    qpData->objValIncumbent = qpDUNES_computeObjectiveValue(qpData);
    if (statusFlag != QPDUNES_OK) {
        return statusFlag;
    }

    /** (1Ba) set up Newton system */
    statusFlag = qpDUNES_setupNewtonSystem(qpData);
    if (statusFlag != QPDUNES_OK) {
        return statusFlag;
    }

    /** (1Bb) factorize Newton system */
    itLogPtr->isHessianRegularized = QPDUNES_FALSE;
    statusFlag = qpDUNES_factorNewtonSystem(qpData, &(itLogPtr->isHessianRegularized));        /* TODO! can we get a problem with on-the-fly regularization in partial refactorization? might only be partially reg.*/
    if (statusFlag != QPDUNES_OK) {
        return statusFlag;
    }

    return QPDUNES_OK;
}

/* main solve function */
return_t qpDUNES_solve(qpData_t* const qpData) {
    return_t statusFlag = QPDUNES_OK; /* generic status flag */
    x_vector_t* restrict xVecTmp = &qpData->xVecTmp;

    uint_t* itCntr = &(qpData->log.numIter);
    itLog_t* itLogPtr = &(qpData->log.itLog[0]);

    *itCntr = 0;
    itLogPtr->itNbr = 0;

    qpDUNES_prepare(qpData);

    /** LOOP OF NONSMOOTH NEWTON ITERATIONS */
    /*  ----------------------------------- */
    for ((*itCntr) = 1; (*itCntr) <= qpData->options.maxIter; (*itCntr)++) {
        itLogPtr->itNbr = *itCntr;

        /** (1) get a step direction:
         *      switch between gradient and Newton steps */

        /** calculate gradient and check gradient norm for convergence */
        statusFlag = qpDUNES_computeNewtonGradient(qpData, &qpData->gradient,
                                                   xVecTmp);
        if (statusFlag != QPDUNES_OK) {
            return statusFlag;
        }

        /** (1Bc) compute step direction */
        /* Force QPDUNES_NH_FAC_BAND_REVERSE */
        statusFlag = qpDUNES_solveNewtonEquationBottomUp(qpData, &(qpData->deltaLambda), &(qpData->cholHessian), &(qpData->gradient));
        if (statusFlag != QPDUNES_OK) {
            return statusFlag;
        }

        /** (2) do QP solution for full step */
        qpDUNES_solveAllLocalQPs(qpData, &(qpData->deltaLambda));
        /* clipping solver: now unsaturated dz is available locally */

        /** (4) determine step length: do line search along the way of the full step
         *      and do the step */
        statusFlag = qpDUNES_determineStepLength(qpData, &(qpData->lambda),
                &(qpData->deltaLambda), &(itLogPtr->numLineSearchIter),
                &(qpData->alpha), &qpData->objValIncumbent,
                itLogPtr->isHessianRegularized);
        if (statusFlag != QPDUNES_ERR_NUMBER_OF_MAX_LINESEARCH_ITERATIONS_REACHED &&
                statusFlag != QPDUNES_ERR_EXCEEDED_MAX_LINESEARCH_STEPSIZE &&
                statusFlag != QPDUNES_OK &&
                statusFlag != QPDUNES_ERR_DECEEDED_MIN_LINESEARCH_STEPSIZE) {
            return statusFlag;
        } else if (statusFlag == QPDUNES_ERR_DECEEDED_MIN_LINESEARCH_STEPSIZE) {
            return QPDUNES_ERR_NEWTON_SYSTEM_NO_ASCENT_DIRECTION;
        }
    }


    /* get number of performed iterations right (itCntr is going one up before realizing it's too big) */
    qpData->log.numIter = qpData->options.maxIter;

    return QPDUNES_ERR_ITERATION_LIMIT_REACHED;
}

/* solve local QPs for a multiplier guess lambda */
return_t qpDUNES_solveAllLocalQPs(qpData_t* const qpData,
const xn_vector_t* const lambda) {
    size_t kk;
    size_t errCntr = 0;
    return_t statusFlag;

    /* 1) update local QP data */
    interval_t* interval;

    for (kk = 0; kk < _NI_ + 1u; kk++) {
        interval = qpData->intervals[kk];
        if (kk < _NI_) {
            qpDUNES_updateVector(
                &interval->lambdaK1, &lambda->data[kk * _NX_], _NX_ );
        }
        if (kk > 0) {
            qpDUNES_updateVector(
                &interval->lambdaK, &lambda->data[(kk - 1u) * _NX_], _NX_ );

            /* Solve local QPs for the previous interval */
            interval = qpData->intervals[kk - 1u];
            clippingQpSolver_updateStageData(
                qpData, interval, &(interval->lambdaK), &(interval->lambdaK1));
            statusFlag = directQpSolver_solveUnconstrained(
                qpData, interval, &interval->qpSolverClipping.qStep);
            if (statusFlag != QPDUNES_OK) { /* note that QPDUNES_OK == 0 */
                errCntr++;
            }
        }
    }

    /* Solve local QPs for the last interval */
    interval = qpData->intervals[_NI_];
    clippingQpSolver_updateStageData(
        qpData, interval, &(interval->lambdaK), &(interval->lambdaK1));
    statusFlag = directQpSolver_solveUnconstrained(
                qpData, interval, &interval->qpSolverClipping.qStep);
    if (statusFlag != QPDUNES_OK) { /* note that QPDUNES_OK == 0 */
        errCntr++;
    }

    if (errCntr) {
        return QPDUNES_ERR_STAGE_QP_INFEASIBLE;
    } else {
        return QPDUNES_OK;
    }
}

return_t qpDUNES_setupNewtonSystem(qpData_t* const qpData) {
    size_t ii, jj, kk;

    xx_matrix_t* restrict xxMatTmp = &qpData->xxMatTmp;
    zx_matrix_t* restrict zxMatTmp = &qpData->zxMatTmp;
    interval_t** restrict intervals = qpData->intervals;
    xn2x_matrix_t* restrict hessian = &qpData->hessian;

    /** calculate hessian */

    /* 1) diagonal blocks */
    /*    E_{k+1} P_{k+1}^-1 E_{k+1}' + C_{k} P_{k} C_{k}'  for projected Hessian  P = Z (Z'HZ)^-1 Z'  */
    for (kk = 0; kk < _NI_; kk++) {
        /* get EPE part -- getInvQ not supported with matrices other than diagonal... is this even possible?*/
        getInvQ(qpData, xxMatTmp, &intervals[kk + 1u]->cholH,
                &intervals[kk + 1u]->y);

        /* add CPC part */
        addCInvHCT(qpData, xxMatTmp, &intervals[kk]->cholH,
                   &intervals[kk]->C, &intervals[kk]->y, zxMatTmp);

        /* write Hessian part */
        #pragma MUST_ITERATE(_NX_, _NX_)
        for (ii = 0; ii < _NX_; ii++) {
            #pragma MUST_ITERATE(_NX_, _NX_)
            for (jj = 0; jj < _NX_; jj++) {
                accHessian(kk, 0, ii, jj) = xxMatTmp->data[ii * _NX_ + jj];
            }
        }

        if (kk == 0) {
            continue;
        }

        multiplyAInvQ(qpData, &qpData->xxMatTmp, &intervals[kk]->C,
                      &intervals[kk]->cholH);

        /* write Hessian part */
        #pragma MUST_ITERATE(_NX_, _NX_)
        for (ii = 0; ii < _NX_; ii++) {
            #pragma MUST_ITERATE(_NX_, _NX_)
            for (jj = 0; jj < _NX_; jj++) {
                /* cheap way of annihilating columns; TODO: make already in multiplication routine! */
                /* check if local constraint lb_x is inactive*/
                if ((intervals[kk]->y.data[2u * jj] <= qpData->options.equalityTolerance) &&
                        /* check if local constraint ub_x is inactive*/
                        (intervals[kk]->y.data[2u * jj + 1u] <= qpData->options.equalityTolerance)) {
                    accHessian(kk, -1, ii, jj) = -xxMatTmp->data[ii * _NX_ + jj];
                } else {
                    /* eliminate column if variable bound is active */
                    accHessian(kk, -1, ii, jj) = 0.0;
                }
            }
        }
    }   /* END OF sub-diagonal block for loop */

    return QPDUNES_OK;
}

return_t qpDUNES_computeNewtonGradient(qpData_t* const qpData,
xn_vector_t* restrict gradient, x_vector_t* restrict gradPiece) {
    size_t k, i;

    const interval_t* restrict interval;
    real_t norm2 = 0.0;

    /* d/(d lambda_ii) for kk=0.._NI_-1 */
    for (k = 0; k < _NI_; k++) {
        interval = qpData->intervals[k];
        /* ( C_kk*z_kk^opt + c_kk ) - x_(kk+1)^opt */
        multiplyCz(qpData, gradPiece, &(interval->C),
                   &(interval->z));

        /* write gradient part */
        #pragma MUST_ITERATE(_NX_, _NX_)
        for (i = 0; i < _NX_; i++) {
            gradPiece->data[i] += interval->c.data[i] -
                                  qpData->intervals[k + 1u]->z.data[i];
            gradient->data[k * _NX_ + i] = gradPiece->data[i];
            norm2 += gradPiece->data[i] * gradPiece->data[i];
        }
    }

    if (norm2 < qpData->options.stationarityTolerance *
                 qpData->options.stationarityTolerance) {
        return QPDUNES_SUCC_OPTIMAL_SOLUTION_FOUND;
    } else {
        return QPDUNES_OK;
    }
}


return_t qpDUNES_factorNewtonSystem(qpData_t* const qpData,
boolean_t* const isHessianRegularized) {
    size_t ii, jj, kk;

    return_t statusFlag;

    real_t minDiagElem = qpData->options.QPDUNES_INFTY;

    xn2x_matrix_t* hessian = &(qpData->hessian);
    xn2x_matrix_t* cholHessian = &(qpData->cholHessian);

    /* Try to factorize Newton Hessian, to check if positive definite */
    /* Force QPDUNES_NH_FAC_BAND_REVERSE */
    statusFlag = qpDUNES_factorizeNewtonHessianBottomUp(qpData, cholHessian,
                                                        hessian);

    /* check maximum diagonal element */
    if (statusFlag == QPDUNES_OK) {
        for (kk = 0; kk < _NI_; kk++) {
            for (ii = 0; ii < _NX_; ii++) {
                if (minDiagElem > accCholHessian(kk, 0, ii, ii)) {
                    minDiagElem = accCholHessian(kk, 0, ii, ii);
                }
            }
        }
    }

    /* regularize if Cholesky failed */
    if ((statusFlag == QPDUNES_ERR_DIVISION_BY_ZERO) ||
            /* or if diagonal elements are too small */
            (minDiagElem < qpData->options.newtonHessDiagRegTolerance)) {
        /* Force QPDUNES_REG_LEVENBERG_MARQUARDT */
        for (kk = 0; kk < _NI_; kk++) {
            for (jj = 0; jj < _NX_; jj++) {
                accHessian(kk, 0, jj, jj) += qpData->options.regParam;
            }
        }
        *isHessianRegularized = QPDUNES_TRUE;

        /* refactor Newton Hessian */
        /* Force QPDUNES_NH_FAC_BAND_REVERSE */
        statusFlag = qpDUNES_factorizeNewtonHessianBottomUp(
            qpData, cholHessian, hessian);  /* refactor full hessian */
    }

    return statusFlag;
}


/*
Bottom-up block-tridiagonal Cholesky for special storage format of Newton
matrix
*/
return_t qpDUNES_factorizeNewtonHessianBottomUp(qpData_t* const qpData,
xn2x_matrix_t* const restrict cholHessian,
xn2x_matrix_t* const restrict hessian) {
    int_t jj, ii, kk, ll;
    real_t sum;

    int_t blockIdxStart = _NI_ - 1;

    /* go by block columns */
    for (kk = blockIdxStart; kk >= 0; --kk) {
        /* go by in-block columns */
        #pragma MUST_ITERATE(_NX_ , _NX_)
        for (jj = _NX_ - 1; jj >= 0; --jj) {
            /* 1) compute diagonal element: ii == jj */
            /* take diagonal element of original */
            sum = accHessian(kk,0,jj,jj);

            /* subtract squared rearpart of corresponding row (transposed access, therefore rest of column): */
            /*  - this diagonal block */
            for( ll = jj+1; ll < _NX_; ++ll ) {
                /* TODO: get rid of transposed access...maybe start to save Hessian also in upper triangular format */
                sum -= accCholHessian(kk,0,ll,jj) * accCholHessian(kk,0,ll,jj); /* transposed access */
            }
            /*  - this row's subdiagonal block */
            if( kk < _NI_-1 ) { /* for all block columns but the last one */
                #pragma MUST_ITERATE(_NX_ , _NX_)
                for( ll = 0; ll < _NX_; ++ll ) {
                    sum -= accCholHessian(kk+1,-1,ll,jj) * accCholHessian(kk+1,-1,ll,jj);   /* transposed access */
                }
            }


            /* 2) check for too small diagonal elements */
            if ( sum < 1.e2f*qpData->options.equalityTolerance ) {  /* matrix not positive definite */
                return QPDUNES_ERR_DIVISION_BY_ZERO;
            }

            accCholHessian(kk,0,jj,jj) = sqrt_f( sum );


            /* 3) write remainder of jj-th column (upwards! via transposed access: jj-th row, leftwards): */
            /*  - this diagonal block */
            for( ii=jj-1; ii>=0; --ii )
            {
                sum = accHessian(kk,0,jj,ii);   /* transposed access */

                /* subtract rear part of this row times rear part of jj-th row */
                /*  - diagonal block */
                for( ll = jj+1; ll < _NX_; ++ll ) {
                    sum -= accCholHessian(kk,0,ll,ii) * accCholHessian(kk,0,ll,jj);     /* transposed access */
                }
                /*  - subdiagonal block */
                if( kk < _NI_-1 ) { /* for all block rows but the last one */
                    #pragma MUST_ITERATE(_NX_ , _NX_)
                    for( ll = 0; ll < _NX_; ++ll ) {
                        sum -= accCholHessian(kk+1,-1,ll,ii) * accCholHessian(kk+1,-1,ll,jj);   /* transposed access */
                    }
                }

                /* write transposed! (otherwise it's upper triangular matrix) */
                accCholHessian(kk,0,jj,ii) = divide_f(sum, accCholHessian(kk,0,jj,jj));
            }
            /*  - following row's subdiagonal block */
            if( kk > 0 ) {  /* for all block rows but the first one */
                #pragma MUST_ITERATE(_NX_ , _NX_)
                for( ii=_NX_-1; ii>=0; --ii )
                {
                    sum = accHessian(kk,-1,jj,ii);  /* transposed access */

                    /* subtract rear part of this row times rear part of jj-th row (only this block is non-zero) */
                    for( ll = jj+1; ll < _NX_; ++ll ) {
                        sum -= accCholHessian(kk,-1,ll,ii) * accCholHessian(kk,0,ll,jj);    /* transposed access */
                    }

                    /* write transposed! (otherwise it's upper triangular matrix) */
                    accCholHessian(kk,-1,jj,ii) = divide_f(sum, accCholHessian(kk,0,jj,jj));
                }
            }
        } /* next column */
    } /* next block column */


    return QPDUNES_OK;
}


/*
Special backsolve for backwards factorized block tridiagonal Newton matrix
*/
return_t qpDUNES_solveNewtonEquationBottomUp(qpData_t* const qpData,
xn_vector_t* const res,
const xn2x_matrix_t* const cholHessian, /**< lower triangular Newton Hessian factor */
const xn_vector_t* const gradient) {
    int_t ii, jj, kk;

    real_t tmp;

    /* solve L^T*x = g */
    #pragma MUST_ITERATE(_NI_, _NI_)
    for (kk = (_NI_ - 1); kk >= 0; kk--) /* go by block rows bottom up */
    {
        #pragma MUST_ITERATE(_NX_, _NX_)
        for (ii = (_NX_ - 1); ii >= 0; ii--) /* go by in-block rows top down */
        {
            tmp = 0.0;
            /* subtract all previously resolved unknowns ... */
            for (jj = ii + 1; jj < _NX_; ++jj) { /* ... of corresponding diagonal block */
                tmp -= accCholHessian(kk,0,jj,ii) * res->data[kk*_NX_+jj]; /* transposed access */
            }
            if (kk < _NI_ - 1) { /* ... of corresponding superdiagonal block, access via following row's subdiagonal block (if not first block row from bottom) */
                #pragma MUST_ITERATE(_NX_, _NX_)
                for (jj = 0; jj < _NX_; ++jj) {
                    tmp -= accCholHessian(kk+1,-1,jj,ii) * res->data[(kk+1)*_NX_+jj];/* TODO: maybe change access pattern, start with superdiag block, so cholH access is more continuous*/
                }
            }

            /* divide by diagonal element */
            res->data[kk * _NX_ + ii] = divide_f(
                gradient->data[kk * _NX_ + ii] + tmp,
                accCholHessian(kk, 0, ii, ii));
        }
    }

    /* solve L*res = x */
    #pragma MUST_ITERATE(_NI_, _NI_)
    for (kk = 0; kk < _NI_; ++kk) /* go by block rows top down */
    {
        #pragma MUST_ITERATE(_NX_, _NX_)
        for (ii = 0; ii < _NX_; ++ii) /* go by in-block rows top down */
        {
            tmp = 0.0;
            /* subtract all previously resolved unknowns ... */
            if (kk > 0) { /* ... of corresponding subdiagonal block (if not first block row) */
                #pragma MUST_ITERATE(_NX_ , _NX_)
                for (jj = 0; jj < _NX_; ++jj) {
                    tmp -= accCholHessian(kk,-1,ii,jj)* res->data[(kk-1)*_NX_+jj];
                }
            }
            for (jj = 0; jj < ii; ++jj) { /* ... of corresponding diagonal block */
                tmp -= accCholHessian(kk,0,ii,jj)* res->data[kk*_NX_+jj];
            }

            /* divide by diagonal element */
            res->data[kk * _NX_ + ii] = divide_f(
                 /* intermediate result of first backsolve is stored in res */
                res->data[kk * _NX_ + ii] + tmp,
                accCholHessian(kk, 0, ii, ii));
        }
    }

    return QPDUNES_OK;
}

return_t qpDUNES_determineStepLength(qpData_t* const qpData,
xn_vector_t* const lambda, xn_vector_t* const deltaLambdaFS,
uint_t* const itCntr, real_t* const alpha, real_t* const objValIncumbent,
boolean_t newtonHessianRegularized) {
    return_t statusFlag;

    size_t kk;

    interval_t* interval;

    size_t nV = _NX_ * _NI_;

    real_t alphaMin = 0.0;
    real_t alphaMax = 1.0;
    real_t alphaASChange = qpData->options.QPDUNES_INFTY;

    *itCntr = 0;

    /* compute minimum step size for active set change */
    /* WARNING: THIS ONLY WORKS IF ALL INTERVALS ARE OF THE SAME TYPE */
    alphaMin = qpData->options.QPDUNES_INFTY;
    for (kk = 0; kk < _NI_ + 1u; kk++) {
        directQpSolver_getMinStepsize(qpData->intervals[kk], &alphaASChange);
        if (alphaASChange < alphaMin) {
            alphaMin = alphaASChange;
        }
    }


    /* take full step and leave */
    if ((alphaMin > 1.0f - qpData->options.equalityTolerance) &&
            (newtonHessianRegularized == QPDUNES_FALSE)) {
        *alpha = 1.0f;

        addVectorScaledVector(lambda, lambda, *alpha, deltaLambdaFS, nV); /* temporary; TODO: move out to mother function */
        for (kk = 0; kk < _NI_ + 1u; kk++) {
            interval = qpData->intervals[kk];
            /* update primal, dual, and internal QP solver variables */
            directQpSolver_doStep(qpData, interval,
                        &(interval->qpSolverClipping.dz), *alpha,
                        &(interval->qpSolverClipping.zUnconstrained),
                        &(interval->z), &(interval->y), &(interval->q),
                        &(interval->p));
        }
        *objValIncumbent = qpDUNES_computeObjectiveValue(qpData);
        return QPDUNES_OK;
    }


    /* do a line search */
    /* force QPDUNES_LS_ACCELERATED_GRADIENT_BISECTION_LS */
    statusFlag = qpDUNES_backTrackingLineSearch(
        qpData, alpha, itCntr, deltaLambdaFS, nV, 0.0f, alphaMax,
        *objValIncumbent);
    if (statusFlag == QPDUNES_ERR_DECEEDED_MIN_LINESEARCH_STEPSIZE) { /* handle backtracking line search errors */
        return statusFlag;
    }
    /* take last alpha that did not yet lead to ascent */
    alphaMax = qpDUNES_fmin(
        alphaMax,
        divide_f(*alpha, qpData->options.lineSearchReductionFactor));
    statusFlag = qpDUNES_bisectionIntervalSearch(
        qpData, alpha, itCntr, deltaLambdaFS, nV, alphaMin, alphaMax);

    /* UPDATE VARIABLES */
    /* lambda */
    addScaledVector(lambda, *alpha, deltaLambdaFS, nV);
    /* stage QP variables */
    for (kk = 0; kk < _NI_ + 1u; kk++) {
        interval = qpData->intervals[kk];
        /* TODO: this might have already been done in line search; do not redo */
        /* update primal, dual, and internal QP solver variables */
        directQpSolver_doStep(qpData, interval,
                    &(interval->qpSolverClipping.dz), *alpha,
                    &(interval->qpSolverClipping.zUnconstrained),
                    &(interval->z), &(interval->y), &(interval->q),
                    &(interval->p));
    }
    *objValIncumbent = qpDUNES_computeObjectiveValue(qpData);

    /* return */
    return statusFlag;
}

return_t qpDUNES_backTrackingLineSearch(qpData_t* const qpData,
real_t* const alpha, uint_t* restrict const itCntr,
const xn_vector_t* const deltaLambdaFS, size_t nV, real_t alphaMin,
real_t alphaMax, real_t const objValIncumbent) {
    real_t objVal;
    real_t minimumProgress = qpData->options.lineSearchMinRelProgress *
        abs_f(objValIncumbent) + qpData->options.lineSearchMinAbsProgress;
    real_t normDeltaLambda = vectorNorm(deltaLambdaFS, nV);

    *alpha = alphaMax;

    /** perform line search */
    for (/*continuous itCntr*/;
            (*itCntr) < qpData->options.maxNumLineSearchIterations;
            (*itCntr)++) {
        /* get objective value */
        objVal = qpDUNES_computeParametricObjectiveValue(qpData, *alpha);

        /* check for progress */
        if (objVal > objValIncumbent + minimumProgress) {
            return QPDUNES_OK;
        } else { /* try smaller step size */
            *alpha = (*alpha) * qpData->options.lineSearchReductionFactor;
        }

        /* ensure minimum step size */
        if (normDeltaLambda * (*alpha - alphaMin)
                < qpData->options.equalityTolerance) {
            *alpha = alphaMin;
            return QPDUNES_ERR_DECEEDED_MIN_LINESEARCH_STEPSIZE;
        }
    }

    return QPDUNES_ERR_NUMBER_OF_MAX_LINESEARCH_ITERATIONS_REACHED;
}

return_t qpDUNES_bisectionIntervalSearch(qpData_t* const qpData,
real_t* const alpha, uint_t* restrict const itCntr,
const xn_vector_t* const deltaLambdaFS, size_t nV, real_t alphaMin,
real_t alphaMax) {
    size_t k, i;
    interval_t** restrict intervals = qpData->intervals;

    z_vector_t* restrict zTry;
    x_vector_t *restrict xVec = &qpData->xVecTmp;
    real_t alphaC = 1.0;
    real_t alphaSlope;
    real_t slopeNormalization;

    /* demand more stationarity for smaller steps */
    slopeNormalization = qpDUNES_fmin(
        1.0f, vectorNorm((vector_t*)deltaLambdaFS, nV));

    /*
    TODO: take line search iterations and maxNumLineSearchRefinementIterations
    together!
    */

    /** (1) check if full step is stationary or even still ascent direction */
    for ( /*continuous itCntr*/;
            (*itCntr) < qpData->options.maxNumLineSearchRefinementIterations;
            (*itCntr)++) {

        /* update z locally according to alpha guess */
        for (k = 0; k <= _NI_; k++) {
            zTry = &(intervals[k]->zVecTmp);
            /* get primal variables for trial step length */
            addVectorScaledVector(
                zTry, &(intervals[k]->qpSolverClipping.zUnconstrained),
                alphaMax, &(intervals[k]->qpSolverClipping.dz),
                intervals[k]->nV);
            directQpSolver_saturateVector(
                qpData, zTry, &(intervals[k]->y), &(intervals[k]->zLow),
                &(intervals[k]->zUpp), intervals[k]->nV);
        }

        /*
        manual gradient computation; TODO: use function, but watch out with
        z, dz, zTry, etc.
        */
        alphaSlope = 0.0f;
        for (k = 0; k < _NI_; k++) {
            /* ( A_kk*x_kk^opt + B_kk*u_kk^opt + c_kk ) - x_(kk+1)^opt */
            multiplyCz(qpData, xVec, &(intervals[k]->C),
                      &(intervals[k]->zVecTmp));

            /* write gradient part */
            #pragma MUST_ITERATE(_NX_, _NX_)
            for (i = 0; i < _NX_; i++) {
                xVec->data[i] += intervals[k]->c.data[i];
                xVec->data[i] -= intervals[k + 1]->zVecTmp.data[i];
                alphaSlope += xVec->data[i] *
                              deltaLambdaFS->data[k * _NX_ + i];
            }
        }

        /* take full step if stationary */
        if (abs_f(alphaSlope) <=
                abs_f(qpData->options.lineSearchStationarityTolerance *
                      slopeNormalization)) {
            *alpha = alphaMax;
            return QPDUNES_OK;
        }

        /* go into normal interval search if full step leads to descent */
        if ((alphaSlope < 0.0f) ^ (slopeNormalization < 0.0f)) {
            break;
        }

        /* increase step size otherwise (full step still leads to ascent) */
        alphaMin = alphaMax;
        alphaMax *= qpData->options.lineSearchIncreaseFactor;

        /* break if maximum step size reached */
        if (alphaMax > qpData->options.lineSearchMaxStepSize) {
            *alpha = alphaMin;
            return QPDUNES_ERR_EXCEEDED_MAX_LINESEARCH_STEPSIZE;
        }
    }


    /** (2) regular bisection interval search */
    for ( /*continuous itCntr*/;
            (*itCntr) < qpData->options.maxNumLineSearchRefinementIterations;
            (*itCntr)++) {
        alphaC = 0.5f * (alphaMin + alphaMax);

        /* update z locally according to alpha guess */
        for (k = 0; k <= _NI_; k++) {
            zTry = &(intervals[k]->zVecTmp);
            /* get primal variables for trial step length */
            addVectorScaledVector(
                zTry, &(intervals[k]->qpSolverClipping.zUnconstrained),
                alphaC, &(intervals[k]->qpSolverClipping.dz),
                intervals[k]->nV);
            directQpSolver_saturateVector(
                qpData, zTry, &(intervals[k]->y), &(intervals[k]->zLow),
                &(intervals[k]->zUpp), intervals[k]->nV);
        }

        /*
        manual gradient computation; TODO: use function, but watch out with
        z, dz, zTry, etc.
        */
        alphaSlope = 0.0f;
        for (k = 0; k < _NI_; k++) {
            /* ( A_kk*x_kk^opt + B_kk*u_kk^opt + c_kk ) - x_(kk+1)^opt */
            multiplyCz(qpData, xVec, &(intervals[k]->C),
                       &(intervals[k]->zVecTmp));

            /* write gradient part */
            #pragma MUST_ITERATE(_NX_, _NX_)
            for (i = 0; i < _NX_; i++) {
                xVec->data[i] += intervals[k]->c.data[i];
                xVec->data[i] -= intervals[k + 1]->zVecTmp.data[i];
                alphaSlope += xVec->data[i] *
                              deltaLambdaFS->data[k * _NX_ + i];
            }
        }

        /* check for stationarity in search direction */
        if (abs_f(alphaSlope) <=
                abs_f(qpData->options.lineSearchStationarityTolerance *
                      slopeNormalization)) {
            *alpha = alphaC;
            return QPDUNES_OK;
        } else {
            /* half interval */
            if (alphaSlope > 0) { /* ascent right of gradient */
                alphaMin = alphaC; /* throw out left interval */
            } else { /* ascent left of gradient */
                alphaMax = alphaC; /* throw out right interval */
            }
        }
    }

    *alpha = alphaC;

    return QPDUNES_ERR_NUMBER_OF_MAX_LINESEARCH_ITERATIONS_REACHED;
}

real_t qpDUNES_computeObjectiveValue(qpData_t* const qpData) {
    size_t k;
    interval_t* restrict interval;
    real_t objVal = 0.0f, qVal, lVal, cVal;

    for (k = 0; k <= _NI_; k++) {
        interval = qpData->intervals[k];

        /* quadratic objective part */
        qVal =
            0.5f * multiplyzHz(&(interval->H), &(interval->z), interval->nV);
        /* linear objective part */
        lVal = scalarProd(&(interval->q), &(interval->z), interval->nV);
        /* constant objective part */
        cVal = interval->p;

        /* sum up */
        objVal += qVal + lVal + cVal;
    }

    return objVal;
}

real_t qpDUNES_computeParametricObjectiveValue(qpData_t* const qpData,
const real_t alpha) {
    size_t k;
    z_vector_t* qTry;
    real_t pTry;
    real_t objVal = 0.0f;
    interval_t* restrict interval;

    for (k = 0; k <= _NI_; k++) {
        interval = qpData->intervals[k];
        qTry = &(interval->zVecTmp);

        /* get primal variables for trial step length */
        directQpSolver_doStep(
            qpData, interval, &(interval->qpSolverClipping.dz), alpha,
            &(interval->z), &(interval->z), &(interval->y), qTry, &pTry);

        /* quadratic objective part */
        interval->optObjVal =
            0.5f * multiplyzHz(&(interval->H), &(interval->z), interval->nV);
        /* linear objective part */
        interval->optObjVal += scalarProd(qTry, &(interval->z), interval->nV);
        /* constant objective part */
        interval->optObjVal += pTry;

        objVal += interval->optObjVal;
    }

    return objVal;
}
