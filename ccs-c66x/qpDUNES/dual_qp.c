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

/**
 *  \file src/dual_qp.c
 *  \author Janick Frasch, Hans Joachim Ferreau
 *  \version 1.0beta
 *  \date 2012
 */


#include "dual_qp.h"
#include "../c66math.h"


/* ----------------------------------------------
 * main solve function
 *
 >>>>>>                                           */
return_t qpDUNES_solve(qpData_t* const qpData) {
    uint_t ii, kk;

    return_t statusFlag = QPDUNES_OK; /* generic status flag */
    int_t lastActSetChangeIdx = _NI_;
    real_t objValIncumbent = qpData->options.QPDUNES_INFTY;

    int_t* itCntr = &(qpData->log.numIter);
    itLog_t* itLogPtr = &(qpData->log.itLog[0]);

    *itCntr = 0;
    itLogPtr->itNbr = 0;

    /** (1) todo: initialize local active sets (at least when using qpOASES) with initial guess from previous iteration */

    /** (2) solve local QP problems for initial guess of lambda */

    /* resolve initial QPs for possibly changed bounds (initial value embedding) */
    for (ii = 0; ii < _NI_ + 1; ++ii) {
        interval_t* interval = qpData->intervals[ii];

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
    objValIncumbent = qpDUNES_computeObjectiveValue(qpData);
    if (statusFlag != QPDUNES_OK) {
        return statusFlag;
    }
    /* get active set of local constraints */
    itLogPtr->nActConstr = qpDUNES_getActSet( qpData, itLogPtr->ieqStatus );
    itLogPtr->nChgdConstr = qpDUNES_compareActSets( qpData,
                                                    (const int_t * const * const ) itLogPtr->ieqStatus, /* explicit casting necessary due to gcc bug */
                                                    (const int_t * const * const ) itLogPtr->prevIeqStatus,
                                                    &lastActSetChangeIdx );


    /** LOOP OF NONSMOOTH NEWTON ITERATIONS */
    /*  ----------------------------------- */
    for ((*itCntr) = 1; (*itCntr) <= qpData->options.maxIter; ++(*itCntr)) {
        itLogPtr->itNbr = *itCntr;


        /** (1) get a step direction:
         *      switch between gradient and Newton steps */
        itLogPtr->isHessianRegularized = QPDUNES_FALSE;
        if ((*itCntr > 1) && (*itCntr - 1 <= qpData->options.nbrInitialGradientSteps)) { /* always do one Newton step first */
            /** (1Aa) get a gradient step */
            qpDUNES_computeNewtonGradient(qpData, &(qpData->gradient),
                    &(qpData->xVecTmp));

            /** (1Ab) do gradient step */
            qpDUNES_copyVector(&(qpData->deltaLambda), &(qpData->gradient),
                    _NI_ * _NX_);
            statusFlag = QPDUNES_OK;
        } else {
            /** (1Ba) set up Newton system */
            statusFlag = qpDUNES_setupNewtonSystem(qpData);
            switch (statusFlag) {
                case QPDUNES_OK:
                    break;
                case QPDUNES_SUCC_OPTIMAL_SOLUTION_FOUND: /* zero gradient norm detected */
                    /* ...and leave */
                    return QPDUNES_SUCC_OPTIMAL_SOLUTION_FOUND;
                default:
                    return statusFlag;
            }

            /** (1Bb) factorize Newton system */
            statusFlag = qpDUNES_factorNewtonSystem(qpData, &(itLogPtr->isHessianRegularized), lastActSetChangeIdx);        /* TODO! can we get a problem with on-the-fly regularization in partial refactorization? might only be partially reg.*/
            switch (statusFlag) {
                case QPDUNES_OK:
                    break;
                default:
                    return statusFlag;
            }

            /** (1Bc) compute step direction */
            /* Force QPDUNES_NH_FAC_BAND_REVERSE */
            statusFlag = qpDUNES_solveNewtonEquationBottomUp(qpData, &(qpData->deltaLambda), &(qpData->cholHessian), &(qpData->gradient));
            if (statusFlag != QPDUNES_OK) {
                return statusFlag;
            }
        }


        /** (2) do QP solution for full step */
        qpDUNES_solveAllLocalQPs(qpData, &(qpData->deltaLambda));
        /* clipping solver: now unsaturated dz is available locally */

        /** (4) determine step length: do line search along the way of the full step
         *      and do the step */
        statusFlag = qpDUNES_determineStepLength(qpData, &(qpData->lambda),
                &(qpData->deltaLambda), &(itLogPtr->numLineSearchIter),
                &(qpData->alpha), &objValIncumbent,
                itLogPtr->isHessianRegularized);
        switch (statusFlag) {
            case QPDUNES_OK:
            case QPDUNES_ERR_NUMBER_OF_MAX_LINESEARCH_ITERATIONS_REACHED:
            case QPDUNES_ERR_EXCEEDED_MAX_LINESEARCH_STEPSIZE:
                break;
            case QPDUNES_ERR_DECEEDED_MIN_LINESEARCH_STEPSIZE: /* deltaLambda is no ascent direction */
                return QPDUNES_ERR_NEWTON_SYSTEM_NO_ASCENT_DIRECTION;
            default:
                return statusFlag;
        }


        /** (5) regular log and display iteration */
        /* get active set of local constraints */
        /* - save old active set */
        if (qpData->options.logLevel >= QPDUNES_LOG_ITERATIONS) {
            itLogPtr->prevIeqStatus = qpData->log.itLog[(*itCntr) - 1].ieqStatus;
        } else {
            /* itLogPtr stays constant */
            /* copy prevIeqStatus */
            for (kk = 0; kk < _NI_ + 1; ++kk) {
                for (ii = 0; ii < _ND(kk)+_NV(kk); ++ii ) {
                    itLogPtr->prevIeqStatus[kk][ii] = itLogPtr->ieqStatus[kk][ii];
                }
            }
        }
        /* - get new active set */
        itLogPtr->nActConstr = qpDUNES_getActSet( qpData, itLogPtr->ieqStatus );
        itLogPtr->nChgdConstr = qpDUNES_compareActSets( qpData,
                                                     (const int_t * const * const ) itLogPtr->ieqStatus, /* explicit casting necessary due to gcc bug */
                                                     (const int_t * const * const ) itLogPtr->prevIeqStatus,
                                                     &lastActSetChangeIdx);
        qpDUNES_logIteration(qpData, itLogPtr, objValIncumbent, lastActSetChangeIdx);
    }


    /* get number of performed iterations right (itCntr is going one up before realizing it's too big) */
    qpData->log.numIter = qpData->options.maxIter;

    return QPDUNES_ERR_ITERATION_LIMIT_REACHED;
}
/*<<< END OF qpDUNES_solve */


/* ----------------------------------------------
 * log all data of this iteration
 *
 >>>>>>                                           */
void qpDUNES_logIteration(  qpData_t* qpData,
                        itLog_t* itLogPtr,
                        real_t objValIncumbent,
                        int_t lastActSetChangeIdx
                        )
{
    itLogPtr->gradNorm = vectorNorm(&(qpData->gradient), _NI_ * _NX_);
    itLogPtr->stepNorm = vectorNorm(&(qpData->deltaLambda), _NI_ * _NX_);
    itLogPtr->stepSize = qpData->alpha;
    itLogPtr->lambdaNorm = vectorNorm(&(qpData->lambda), _NI_ * _NX_);
    itLogPtr->objVal = objValIncumbent;
    itLogPtr->lastActSetChangeIdx = lastActSetChangeIdx;
}
/*<<< END OF qpDUNES_logIteration */


/* ----------------------------------------------
 * update all qSteps and pSteps (linear and constant objective function contribution) of the local QPs
 *
 >>>>>>                                           */
return_t qpDUNES_updateAllLocalQPs( qpData_t* const qpData,
                                    const xn_vector_t* const lambda
                                    )
{
    int_t kk;
    interval_t* interval;

    /* first interval: */
    interval = qpData->intervals[0];
    qpDUNES_updateVector( &(interval->lambdaK1), &(lambda->data[0]), _NX_ );
    /* intermediate intervals: */
    for (kk = 1; kk < _NI_; ++kk) {
        interval = qpData->intervals[kk];
        qpDUNES_updateVector( &(interval->lambdaK), &(lambda->data[(kk - 1) * _NX_]), _NX_ );
        qpDUNES_updateVector( &(interval->lambdaK1), &(lambda->data[kk * _NX_]), _NX_ );
    }
    /* last interval: */
    interval = qpData->intervals[_NI_];
    qpDUNES_updateVector( &(interval->lambdaK), &(lambda->data[(_NI_ - 1) * _NX_]), _NX_ );

    for (kk = 0; kk < _NI_ + 1; ++kk) {
        interval = qpData->intervals[kk];
        switch (interval->qpSolverSpecification) {
        case QPDUNES_STAGE_QP_SOLVER_CLIPPING:
            clippingQpSolver_updateStageData( qpData, interval, &(interval->lambdaK), &(interval->lambdaK1) );
            break;
        case QPDUNES_STAGE_QP_SOLVER_QPOASES:
            assert(0);
            break;
        default:
            return QPDUNES_ERR_UNKNOWN_ERROR;
        }
    }

    return QPDUNES_OK;
}
/*<<< END OF qpDUNES_updateAllLocalQPs */


/* ----------------------------------------------
 * solve local QPs for a multiplier guess lambda
 *
 >>>>>>                                           */
return_t qpDUNES_solveAllLocalQPs(  qpData_t* const qpData,
                                const xn_vector_t* const lambda
                                )
{
    int_t kk;
    int_t errCntr = 0;
    return_t statusFlag;

    /* 1) update local QP data */
    qpDUNES_updateAllLocalQPs(qpData, lambda);

    /* 2) solve local QPs */
    /* TODO: check what happens in case of errors (return)*/
    /* Note: const variables are predetermined shared (at least on apple)*/
    for (kk = 0; kk < _NI_ + 1; ++kk) {
        statusFlag = qpDUNES_solveLocalQP(qpData, qpData->intervals[kk]);
        if (statusFlag != QPDUNES_OK) { /* note that QPDUNES_OK == 0 */
            errCntr++;
        }
    }

    if (errCntr > 0) {
        return QPDUNES_ERR_STAGE_QP_INFEASIBLE;
    }

    return QPDUNES_OK;
}
/*<<< END OF qpDUNES_solveAllLocalQPs */


/* ----------------------------------------------
 * solve local QP
 *
 >>>>>>                                           */
return_t qpDUNES_solveLocalQP(  qpData_t* const qpData,
                            interval_t* const interval
                            )
{
    return directQpSolver_solveUnconstrained(qpData, interval, &(interval->qpSolverClipping.qStep)); /* solve QPs in first-order term updates only, to mimic homotopy */
}
/*<<< END OF qpDUNES_solveLocalQP */


/* ----------------------------------------------
 * ...
 *
 >>>>>>                                           */
return_t qpDUNES_setupNewtonSystem( qpData_t* const qpData
                                    )
{
    int_t ii, jj, kk;

    x_vector_t* restrict xVecTmp = &(qpData->xVecTmp);
    xx_matrix_t* restrict xxMatTmp = &(qpData->xxMatTmp);
    zx_matrix_t* restrict zxMatTmp = &(qpData->zxMatTmp);
    interval_t** restrict intervals = qpData->intervals;
    xn2x_matrix_t* restrict hessian = &(qpData->hessian);

    /** calculate gradient and check gradient norm for convergence */
    qpDUNES_computeNewtonGradient(qpData, &(qpData->gradient), xVecTmp);
    if ((vectorNorm(&(qpData->gradient), _NX_ * _NI_)
            < qpData->options.stationarityTolerance)) {
        return QPDUNES_SUCC_OPTIMAL_SOLUTION_FOUND;
    }


    /** calculate hessian */

    /* 1) diagonal blocks */
    /*    E_{k+1} P_{k+1}^-1 E_{k+1}' + C_{k} P_{k} C_{k}'  for projected Hessian  P = Z (Z'HZ)^-1 Z'  */
    for (kk = 0; kk < _NI_; ++kk) {
        /* check whether block needs to be recomputed */
        if ( (intervals[kk]->actSetHasChanged == QPDUNES_TRUE) || (intervals[kk+1]->actSetHasChanged == QPDUNES_TRUE) ) {
            /* get EPE part */
            getInvQ(qpData, xxMatTmp, &(intervals[kk + 1]->cholH), intervals[kk + 1]->nV); /* getInvQ not supported with matrices other than diagonal... is this even possible? */

            /* Annihilate columns in invQ; WARNING: this can really only be applied for diagonal matrices */
            qpDUNES_makeMatrixDense(xxMatTmp, _NX_, _NX_);
            #pragma MUST_ITERATE(_NX_, _NX_)
            for (ii = 0; ii < _NX_; ++ii) {
                if ((intervals[kk + 1]->y.data[2 * ii] >= qpData->options.equalityTolerance) ||     /* check if local constraint lb_x is active*/
                    (intervals[kk + 1]->y.data[2 * ii + 1] >= qpData->options.equalityTolerance))   /* check if local constraint ub_x is active*/   /* WARNING: weakly active constraints are excluded here!*/
                {
                    xxMatTmp->data[ii * _NX_ + ii] = 0.0f;
                }
            }

            /* add CPC part */
            addCInvHCT(qpData, xxMatTmp, &(intervals[kk]->cholH), &(intervals[kk]->C), &(intervals[kk]->y), zxMatTmp);

            /* write Hessian part */
            #pragma MUST_ITERATE(_NX_, _NX_)
            for (ii = 0; ii < _NX_; ++ii) {
                #pragma MUST_ITERATE(_NX_, _NX_)
                for (jj = 0; jj < _NX_; ++jj) {
                    accHessian( kk, 0, ii, jj ) = xxMatTmp->data[ii * _NX_ + jj];
                    /* clean xxMatTmp */
                    xxMatTmp->data[ii * _NX_ + jj] = 0.0f; /* TODO: this cleaning part is probably not needed, but we need to be very careful if we decide to leave it out! */
                }
            }
        }
    }   /* END OF diagonal block for loop */

    /* 2) sub-diagonal blocks */
    for (kk = 1; kk < _NI_; ++kk) {
        if (intervals[kk]->actSetHasChanged == QPDUNES_TRUE) {
            multiplyAInvQ( qpData, &(qpData->xxMatTmp), &(intervals[kk]->C), &(intervals[kk]->cholH) );

            /* write Hessian part */
            #pragma MUST_ITERATE(_NX_, _NX_)
            for (ii=0; ii<_NX_; ++ii) {
                #pragma MUST_ITERATE(_NX_, _NX_)
                for (jj=0; jj<_NX_; ++jj) {
                    /* cheap way of annihilating columns; TODO: make already in multiplication routine! */
                    if ( ( intervals[kk]->y.data[2*jj] <= qpData->options.equalityTolerance ) &&        /* check if local constraint lb_x is inactive*/
                         ( intervals[kk]->y.data[2*jj+1] <= qpData->options.equalityTolerance ) )       /* check if local constraint ub_x is inactive*/
                    {
                        accHessian( kk, -1, ii, jj ) = - xxMatTmp->data[ii * _NX_ + jj];
                    } else {
                        /* eliminate column if variable bound is active */
                        accHessian( kk, -1, ii, jj ) = 0.;
                    }
                }
            }
        }
    }   /* END OF sub-diagonal block for loop */

    return QPDUNES_OK;
}
/*<<< END OF qpDUNES_setupNewtonSystem */


return_t qpDUNES_computeNewtonGradient(qpData_t* const qpData,
xn_vector_t* restrict gradient, x_vector_t* restrict gradPiece) {
    size_t k, i;

    interval_t** intervals = qpData->intervals;

    /* d/(d lambda_ii) for kk=0.._NI_-1 */
    for (k = 0; k < _NI_; k++) {
        /* ( C_kk*z_kk^opt + c_kk ) - x_(kk+1)^opt */
        multiplyCz(qpData, gradPiece, &(intervals[k]->C),
                   &(intervals[k]->z));
        addToVector(gradPiece, &(intervals[k]->c), _NX_);
        subtractFromVector(gradPiece, &(intervals[k + 1u]->z), _NX_);

        /* write gradient part */
        #pragma MUST_ITERATE(_NX_, _NX_)
        for (i = 0; i < _NX_; i++) {
            gradient->data[k * _NX_ + i] = gradPiece->data[i];
        }
    }
    return QPDUNES_OK;
}


return_t qpDUNES_factorNewtonSystem( qpData_t* const qpData,
                                     boolean_t* const isHessianRegularized,
                                     int_t lastActSetChangeIdx
                                     )
{
    int_t ii, jj, kk;

    return_t statusFlag;

    real_t minDiagElem = qpData->options.QPDUNES_INFTY;

    xn2x_matrix_t* hessian = &(qpData->hessian);
    xn2x_matrix_t* cholHessian = &(qpData->cholHessian);

    /* Try to factorize Newton Hessian, to check if positive definite */
    /* Force QPDUNES_NH_FAC_BAND_REVERSE */
    statusFlag = qpDUNES_factorizeNewtonHessianBottomUp( qpData, cholHessian, hessian, lastActSetChangeIdx, isHessianRegularized );

    /* check maximum diagonal element */
    if (statusFlag == QPDUNES_OK) {
        for (kk = 0; kk < _NI_; ++kk) {
            for (ii = 0; ii < _NX_; ++ii) {
                if (minDiagElem > accCholHessian(kk, 0, ii, ii) ) {
                    minDiagElem = accCholHessian(kk, 0, ii, ii);
                }
            }
        }
    }


    if ( ( statusFlag == QPDUNES_ERR_DIVISION_BY_ZERO ) ||                  /* regularize if Cholesky failed */
         ( minDiagElem < qpData->options.newtonHessDiagRegTolerance ) )     /* or if diagonal elements are too small */
    {
        /* Force QPDUNES_REG_LEVENBERG_MARQUARDT */
        for (kk = 0; kk < _NI_; ++kk) {
                for (jj = 0; jj < _NX_; ++jj) {
                    accHessian( kk, 0, jj, jj )+= qpData->options.regParam;
                }
            }
        *isHessianRegularized = QPDUNES_TRUE;

        /* refactor Newton Hessian */
        /* Force QPDUNES_NH_FAC_BAND_REVERSE */
        statusFlag = qpDUNES_factorizeNewtonHessianBottomUp( qpData, cholHessian, hessian, _NI_+1, isHessianRegularized );  /* refactor full hessian */
        if ( statusFlag != QPDUNES_OK ) {
            return statusFlag;
        }
    }
    else {
        if ( statusFlag != QPDUNES_OK ) {
            return statusFlag;
        }
    }

    return QPDUNES_OK;
}
/*<<< END OF qpDUNES_factorNewtonSystem */


/* ----------------------------------------------
 * Bottom-up block-tridiagonal Cholesky for special storage format of Newton matrix
 *
 >>>>>>                                           */
return_t qpDUNES_factorizeNewtonHessianBottomUp( qpData_t* const qpData,
                                              xn2x_matrix_t* const restrict cholHessian,
                                              xn2x_matrix_t* const restrict hessian,
                                              int_t lastActSetChangeIdx,            /**< index from where the reverse factorization is restarted */
                                              boolean_t* restrict isHessianRegularized
                                              )
{
    int_t jj, ii, kk, ll;
    real_t sum;

    int_t blockIdxStart = (lastActSetChangeIdx >= 0)  ?
        (lastActSetChangeIdx < _NI_ - 1 ? lastActSetChangeIdx : _NI_ - 1)  :
        -1;

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
/*<<< END OF qpDUNES_factorizeNewtonHessianBottomUp */


/* ----------------------------------------------
 * special backsolve for backwards factorized block tridiagonal Newton matrix
 *
 >>>>>>                                           */
return_t qpDUNES_solveNewtonEquationBottomUp(   qpData_t* const qpData,
                                            xn_vector_t* const res,
                                            const xn2x_matrix_t* const cholHessian, /**< lower triangular Newton Hessian factor */
                                            const xn_vector_t* const gradient   )
{
    int_t ii, jj, kk;

    real_t sum;

    /* solve L^T*x = g */
    #pragma MUST_ITERATE(_NI_, _NI_)
    for (kk = (_NI_ - 1); kk >= 0; kk--) /* go by block rows bottom up */
    {
        #pragma MUST_ITERATE(_NX_, _NX_)
        for (ii = (_NX_ - 1); ii >= 0; ii--) /* go by in-block rows top down */
        {
            sum = gradient->data[kk * _NX_ + ii];
            /* subtract all previously resolved unknowns ... */
            for (jj = ii + 1; jj < _NX_; ++jj) { /* ... of corresponding diagonal block */
                sum -= accCholHessian(kk,0,jj,ii)* res->data[kk*_NX_+jj]; /* transposed access */
            }
            if (kk < _NI_ - 1) { /* ... of corresponding superdiagonal block, access via following row's subdiagonal block (if not first block row from bottom) */
                #pragma MUST_ITERATE(_NX_, _NX_)
                for (jj = 0; jj < _NX_; ++jj) {
                    sum -= accCholHessian(kk+1,-1,jj,ii)* res->data[(kk+1)*_NX_+jj];/* TODO: maybe change access pattern, start with superdiag block, so cholH access is more continuous*/
                }
            }

            /* divide by diagonal element */
            res->data[kk * _NX_ + ii] = divide_f(sum, accCholHessian(kk,0,ii,ii));
        }
    }

    /* solve L*res = x */
    #pragma MUST_ITERATE(_NI_, _NI_)
    for (kk = 0; kk < _NI_; ++kk) /* go by block rows top down */
    {
        #pragma MUST_ITERATE(_NX_, _NX_)
        for (ii = 0; ii < _NX_; ++ii) /* go by in-block rows top down */
        {
            sum = res->data[kk * _NX_ + ii]; /* intermediate result of first backsolve is stored in res */
            /* subtract all previously resolved unknowns ... */
            if (kk > 0) { /* ... of corresponding subdiagonal block (if not first block row) */
                #pragma MUST_ITERATE(_NX_ , _NX_)
                for (jj = 0; jj < _NX_; ++jj) {
                    sum -= accCholHessian(kk,-1,ii,jj)* res->data[(kk-1)*_NX_+jj];
                }
            }
            for (jj = 0; jj < ii; ++jj) { /* ... of corresponding diagonal block */
                sum -= accCholHessian(kk,0,ii,jj)* res->data[kk*_NX_+jj];
            }

            /* divide by diagonal element */
            res->data[kk * _NX_ + ii] = divide_f(sum, accCholHessian(kk,0,ii,ii));
        }
    }

    return QPDUNES_OK;
}
/*<<< END OF qpDUNES_solveNewtonEquationBottomUp */


/* ----------------------------------------------
 * ...
 *
 >>>>>>                                           */
return_t qpDUNES_determineStepLength(   qpData_t* const qpData,
                                    xn_vector_t* const lambda,
                                    xn_vector_t* const deltaLambdaFS,
                                    uint_t* const itCntr,
                                    real_t* const alpha,
                                    real_t* const objValIncumbent,
                                    boolean_t newtonHessianRegularized
                                    )
{
    return_t statusFlag;

    int_t kk;

    interval_t* interval;

    int_t nV = _NX_ * _NI_;

    real_t alphaMin = 0.;
    real_t alphaMax = 1.;
    real_t alphaASChange = qpData->options.QPDUNES_INFTY;

    xn_vector_t* lambdaTry = &(qpData->xnVecTmp);

    *itCntr = 0;

    /* compute minimum step size for active set change */
    /* WARNING: THIS ONLY WORKS IF ALL INTERVALS ARE OF THE SAME TYPE */
    if ( qpData->intervals[0]->qpSolverSpecification    == QPDUNES_STAGE_QP_SOLVER_CLIPPING )
    {
        alphaMin = qpData->options.QPDUNES_INFTY;
    }
    for ( kk = 0; kk < _NI_ + 1; ++kk )
    {
        if (qpData->intervals[kk]->qpSolverSpecification == QPDUNES_STAGE_QP_SOLVER_CLIPPING)
        {
            directQpSolver_getMinStepsize(qpData->intervals[kk], &alphaASChange);
            if (alphaASChange < alphaMin) {
                alphaMin = alphaASChange;
            }
        }
        /* TODO: compute minimum stepsize for qpOASES */
    }


    /* take full step and leave */
    if ( (alphaMin > 1.0f - qpData->options.equalityTolerance) && (newtonHessianRegularized == QPDUNES_FALSE) )
    {
        *alpha = 1.0f;

        addVectorScaledVector(lambda, lambda, *alpha, deltaLambdaFS, nV); /* temporary; TODO: move out to mother function */
        for (kk = 0; kk < _NI_ + 1; ++kk) {
            interval = qpData->intervals[kk];
            /* update primal, dual, and internal QP solver variables */
            switch (interval->qpSolverSpecification) {
            case QPDUNES_STAGE_QP_SOLVER_CLIPPING:
                directQpSolver_doStep(qpData, interval,
                        &(interval->qpSolverClipping.dz), *alpha,
                        &(interval->qpSolverClipping.zUnconstrained),
                        &(interval->z), &(interval->y), &(interval->q),
                        &(interval->p));
                break;

            case QPDUNES_STAGE_QP_SOLVER_QPOASES:
                assert(0);
                break;

            default:
                return QPDUNES_ERR_UNKNOWN_ERROR;
            }
        }
        *objValIncumbent = qpDUNES_computeObjectiveValue(qpData);
        return QPDUNES_OK;
    }


    /* do a line search */
    /* force QPDUNES_LS_ACCELERATED_GRADIENT_BISECTION_LS */
    statusFlag = qpDUNES_backTrackingLineSearch(qpData, alpha, itCntr, deltaLambdaFS, lambdaTry, nV, 0., alphaMax, *objValIncumbent);
    if (statusFlag == QPDUNES_ERR_DECEEDED_MIN_LINESEARCH_STEPSIZE) { /* handle backtracking line search errors */
        return statusFlag;
    }
    alphaMax = qpDUNES_fmin(alphaMax, divide_f(*alpha, qpData->options.lineSearchReductionFactor)); /* take last alpha that did not yet lead to ascent */
    statusFlag = qpDUNES_bisectionIntervalSearch( qpData, alpha, itCntr, deltaLambdaFS, lambdaTry, nV, alphaMin, alphaMax );

    /* UPDATE VARIABLES */
    /* lambda */
    addScaledVector(lambda, *alpha, deltaLambdaFS, nV);
    /* stage QP variables */
    for (kk = 0; kk < _NI_ + 1; ++kk) {
        interval = qpData->intervals[kk];
        /* TODO: this might have already been done in line search; do not redo */
        /* update primal, dual, and internal QP solver variables */
        switch (interval->qpSolverSpecification) {
        case QPDUNES_STAGE_QP_SOLVER_CLIPPING:
            directQpSolver_doStep(qpData, interval,
                    &(interval->qpSolverClipping.dz), *alpha,
                    &(interval->qpSolverClipping.zUnconstrained),
                    &(interval->z), &(interval->y), &(interval->q),
                    &(interval->p));
            break;

        case QPDUNES_STAGE_QP_SOLVER_QPOASES:
            assert(0);
            break;

        default:
            return QPDUNES_ERR_UNKNOWN_ERROR;
        }
    }
    *objValIncumbent = qpDUNES_computeObjectiveValue(qpData);

    /* return */
    return statusFlag;
}
/*<<< END OF qpDUNES_determineStepLength */


return_t qpDUNES_backTrackingLineSearch(qpData_t* const qpData,
real_t* const alpha, uint_t* const itCntr,
const xn_vector_t* const deltaLambdaFS, xn_vector_t* const lambdaTry,
size_t nV, real_t alphaMin, real_t alphaMax, real_t const objValIncumbent) {
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
real_t* const alpha, uint_t* const itCntr,
const xn_vector_t* const deltaLambdaFS, xn_vector_t* const lambdaTry,
size_t nV, real_t alphaMin, real_t alphaMax) {
    size_t k, i;
    interval_t* restrict interval;

    z_vector_t* restrict zTry;
    real_t alphaC = 1.0;
    real_t alphaSlope;
    real_t slopeNormalization;

    /* demand more stationarity for smaller steps */
    slopeNormalization = qpDUNES_fmin(
        1.0f, vectorNorm((vector_t*)deltaLambdaFS, nV));

    /* todo: get memory passed on from determine step length */
    xn_vector_t* restrict gradientTry = &(qpData->xnVecTmp2);
    /* todo: no need to recompute gradient in next Newton iteration! */

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
            interval = qpData->intervals[k];
            zTry = &(interval->zVecTmp);
            /* get primal variables for trial step length */
            addVectorScaledVector(
                zTry, &(interval->qpSolverClipping.zUnconstrained), alphaMax,
                &(interval->qpSolverClipping.dz), interval->nV);
            directQpSolver_saturateVector(
                qpData, zTry, &(interval->y), &(interval->zLow),
                &(interval->zUpp), interval->nV);
        }

        /*
        manual gradient computation; TODO: use function, but watch out with
        z, dz, zTry, etc.
        */
        for (k = 0; k < _NI_; k++) {
            /* ( A_kk*x_kk^opt + B_kk*u_kk^opt + c_kk ) - x_(kk+1)^opt */
            multiplyCz(qpData, &(qpData->xVecTmp),
                      &(qpData->intervals[k]->C),
                      &(qpData->intervals[k]->zVecTmp));
            addToVector(
                &(qpData->xVecTmp), &(qpData->intervals[k]->c), _NX_);

            /* write gradient part */
            #pragma MUST_ITERATE(_NX_, _NX_)
            for (i = 0; i < _NX_; i++) {
                qpData->xVecTmp.data[i] -=
                    qpData->intervals[k + 1]->zVecTmp.data[i];
                gradientTry->data[k * _NX_ + i] = qpData->xVecTmp.data[i];
            }
        }
        alphaSlope = scalarProd(gradientTry, deltaLambdaFS, nV);

        /* take full step if stationary */
        if (abs_f(divide_f(alphaSlope, slopeNormalization))
                <= qpData->options.lineSearchStationarityTolerance) {
            *alpha = alphaMax;
            return QPDUNES_OK;
        }

        /* go into normal interval search if full step leads to descent */
        if (divide_f(alphaSlope, slopeNormalization) < 0.0f) {
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
        alphaC = 0.5 * (alphaMin + alphaMax);

        /* update z locally according to alpha guess */
        for (k = 0; k <= _NI_; k++) {
            interval = qpData->intervals[k];
            zTry = &(interval->zVecTmp);
            /* get primal variables for trial step length */
            addVectorScaledVector(
                zTry, &(interval->qpSolverClipping.zUnconstrained), alphaC,
                &(interval->qpSolverClipping.dz), interval->nV);
            directQpSolver_saturateVector(
                qpData, zTry, &(interval->y), &(interval->zLow),
                &(interval->zUpp), interval->nV);
        }

        /*
        manual gradient computation; TODO: use function, but watch out with
        z, dz, zTry, etc.
        */
        for (k = 0; k < _NI_; k++) {
            /* ( A_kk*x_kk^opt + B_kk*u_kk^opt + c_kk ) - x_(kk+1)^opt */
            multiplyCz(qpData, &(qpData->xVecTmp),
                       &(qpData->intervals[k]->C),
                       &(qpData->intervals[k]->zVecTmp));
            addToVector(
                &(qpData->xVecTmp), &(qpData->intervals[k]->c), _NX_);

            /* write gradient part */
            #pragma MUST_ITERATE(_NX_, _NX_)
            for (i = 0; i < _NX_; i++) {
                qpData->xVecTmp.data[i] -=
                    qpData->intervals[k + 1]->zVecTmp.data[i];
                gradientTry->data[k * _NX_ + i] = qpData->xVecTmp.data[i];
            }
        }
        alphaSlope = scalarProd(gradientTry, deltaLambdaFS, nV);

        /* check for stationarity in search direction */
        if (abs_f(divide_f(alphaSlope, slopeNormalization))
                <= qpData->options.lineSearchStationarityTolerance) {
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


void qpDUNES_getPrimalSol(const qpData_t* const qpData, real_t* const z) {
    size_t k;

    for (k = 0; k <= _NI_; k++) {
        qpDUNES_copyArray(&(z[k * _NZ_]), qpData->intervals[k]->z.data,
                          qpData->intervals[k]->nV);
    }
}


real_t qpDUNES_computeObjectiveValue(qpData_t* const qpData) {
    size_t k;
    interval_t* restrict interval;
    real_t objVal = 0.0f;

    for (k = 0; k <= _NI_; k++) {
        interval = qpData->intervals[k];

        /* quadratic objective part */
        interval->optObjVal =
            0.5f * multiplyzHz(&(interval->H), &(interval->z), interval->nV);
        /* linear objective part */
        interval->optObjVal += scalarProd(&(interval->q), &(interval->z),
                                          interval->nV);
        /* constant objective part */
        interval->optObjVal += interval->p;

        /* sum up */
        objVal += interval->optObjVal;
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


/* ----------------------------------------------
 * Get number of active local constraints
 *
 *   Note: this overwrites AS (TODO: is this qpData->ieqStatus?)
 *
 >>>>>>                                           */
uint_t qpDUNES_getActSet( const qpData_t* const qpData,
                       int_t * const * const actSetStatus) {
    uint_t ii = 0;
    uint_t kk = 0;

    uint_t nActConstr = 0;

    for (kk = 0; kk < _NI_ + 1; ++kk) {
        for (ii = 0; ii < _ND(kk) + _NV(kk); ++ii ) {
            /* TODO: make this quick hack clean for general multiplier usage...! */
            /* go through multipliers in pairs by two */
            if ( qpData->intervals[kk]->y.data[2*ii] > qpData->options.equalityTolerance ) { /* lower bound active */
                actSetStatus[kk][ii] = -1;
                ++nActConstr;
            }
            else {
                if ( qpData->intervals[kk]->y.data[2*ii+1] > qpData->options.equalityTolerance ) { /* upper bound active */
                    actSetStatus[kk][ii] = 1;
                    ++nActConstr;
                }
                else {      /* no constraint bound active */
                    actSetStatus[kk][ii] = 0;
                }
            }
        }
    }

    return nActConstr;
}
/*<<< END OF qpDUNES_countActConstr */


/* ----------------------------------------------
 * Get number of differences between two active sets
 * TODO: do this based on alpha, no need for actually comparing active sets !!
*/
uint_t qpDUNES_compareActSets( const qpData_t* const qpData,
                            const int_t * const * const newActSetStatus,
                            const int_t * const * const oldActSetStatus,
                            int_t * const lastActSetChangeIdx) {
    uint_t ii, kk;
    uint_t nChgdConstr = 0;

    *lastActSetChangeIdx = -1;

    for (kk = 0; kk < _NI_+1; ++kk) {
        qpData->intervals[kk]->actSetHasChanged = QPDUNES_FALSE;
        for (ii = 0; ii < _ND(kk)+_NV(kk); ++ii ) {
            /* TODO: maybe include check whether lb = ub? Is a jump from lb to ub (or even to inactive, though unlikely) in this case really an active set change? */
            if( newActSetStatus[kk][ii] != oldActSetStatus[kk][ii] ) {
                ++nChgdConstr;
                qpData->intervals[kk]->actSetHasChanged = QPDUNES_TRUE;
                *lastActSetChangeIdx = kk;
            }
        }
    }

    return nChgdConstr;
}
