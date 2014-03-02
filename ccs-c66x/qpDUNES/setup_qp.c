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

#include <assert.h>

#include "dual_qp.h"
#include "setup_qp.h"
#include "../c66math.h"


void qpDUNES_indicateDataChange(qpData_t* const qpData) {
    assert(qpData);

    size_t k, i;

    /*
    initialize prevIeqStatus to safe values when data was changed to force
    Hessian refactorization
    */
    for (k = 0; k < _NI_ + 1u; k++) {
        for (i = 0; i < _ND(k) + _NV(k); i++) {
            /* some safe dummy value */
            qpData->log.itLog[0].prevIeqStatus[k][i] = -42;
        }
    }
}

return_t qpDUNES_setupRegularInterval(qpData_t* const qpData,
interval_t* interval, const real_t* const Q_, const real_t* const R_,
const real_t* const C_, const real_t* const c_, const real_t* const zLow_,
const real_t* const zUpp_) {
    size_t i;
    size_t nV = interval->nV;

    vv_matrix_t* H = &(interval->H);
    xz_matrix_t* C = &(interval->C);

    sparsityType_t sparsityQ;
    sparsityType_t sparsityR;

    /** (1) quadratic term of cost function */
    /* assemble Hessian */
    /* TODO: move Q, R out to MPC module */
    /* detect sparsity of Q, R */
    sparsityQ = (Q_ != 0) ?
        qpDUNES_detectMatrixSparsity(Q_, _NX_, _NX_) : QPDUNES_IDENTITY;
    sparsityR = (R_ != 0) ?
        qpDUNES_detectMatrixSparsity(R_, _NU_, _NU_) : QPDUNES_IDENTITY;

    assert(sparsityQ == QPDUNES_DIAGONAL);
    assert(sparsityR == QPDUNES_DIAGONAL);

    /* write Hessian blocks */
    H->sparsityType = QPDUNES_DIAGONAL;
    /* write diagonal in first line for cache efficiency */
    /* Q part */
    for (i = 0; i < _NX_; i++) {
        accH(0, i) = Q_[i * _NX_ + i];
    }
    /* R part */
    for (i = 0; i < _NU_; i++) {
        accH(0, _NX_ + i) = R_[i * _NU_ + i];
    }

    /** (2) linear term of cost function -- always 0 */

    /** (3) dynamic system */
    if (C->sparsityType == QPDUNES_MATRIX_UNDEFINED) {
        C->sparsityType = QPDUNES_DENSE;
    }
    if (C_) {
        /* set up C directly */
        qpDUNES_updateMatrixData((matrix_t*)C, C_, _NX_, _NZ_);
    } else {
        return QPDUNES_ERR_INVALID_ARGUMENT;
    }

    if (c_) {
        qpDUNES_setupVector((vector_t*)&(interval->c), c_, _NX_);
    }
    else {
        qpDUNES_setupZeroVector((vector_t*)&(interval->c), _NX_);
    }

    /** (4) bounds */
    qpDUNES_setupUniformVector(
        (vector_t*)&(interval->zLow), -qpData->options.QPDUNES_INFTY, nV);
    qpDUNES_updateSimpleBoundVector(
        qpData, (vector_t*)&(interval->zLow), zLow_, NULL, NULL);
    qpDUNES_setupUniformVector(
        (vector_t*)&(interval->zUpp), qpData->options.QPDUNES_INFTY, nV);
    qpDUNES_updateSimpleBoundVector(
        qpData, (vector_t*)&(interval->zUpp), zUpp_, NULL, NULL);

    return QPDUNES_OK;
}


return_t qpDUNES_setupFinalInterval(qpData_t* const qpData,
interval_t* interval, const real_t* const H_, const real_t* const zLow_,
const real_t* const zUpp_) {
    size_t nV = interval->nV;

    vv_matrix_t* H = &(interval->H);

    /** (1) quadratic term of cost function */
    if (H_) {    /* H given */
        H->sparsityType = qpDUNES_detectMatrixSparsity(H_, nV, nV);
        qpDUNES_updateMatrixData((matrix_t*)H, H_, nV, nV);
    } else {
        qpDUNES_setupScaledIdentityMatrix(_NX_, qpData->options.regParam,
                                          (matrix_t*)H);
    }

    /** (2) linear term of cost function -- unused */

    /** (3) local bounds */
    qpDUNES_setupUniformVector(
        (vector_t*)&(interval->zLow), -qpData->options.QPDUNES_INFTY, nV);
    qpDUNES_updateVector((vector_t*)&(interval->zLow), zLow_, nV);
    qpDUNES_setupUniformVector(
        (vector_t*)&(interval->zUpp), qpData->options.QPDUNES_INFTY, nV);
    qpDUNES_updateVector((vector_t*)&(interval->zUpp), zUpp_, nV);

    return QPDUNES_OK;
}


return_t qpDUNES_updateIntervalData(qpData_t* const qpData,
interval_t* interval, const real_t* const C_, const real_t* const c_,
const real_t* const zLow_, const real_t* const zUpp_) {
    size_t nV = interval->nV;

    qpDUNES_updateVector((vector_t*)&(interval->zLow), zLow_, nV);
    qpDUNES_updateVector((vector_t*)&(interval->zUpp), zUpp_, nV);

    /** re-factorize Hessian for direct QP solver if needed */
    if (C_ || c_) {
        /** copy data */
        qpDUNES_updateMatrixData((matrix_t*)&(interval->C), C_, _NX_, _NZ_);
        qpDUNES_updateVector((vector_t*)&(interval->c), c_, _NX_);
        qpDUNES_setupStageQP(qpData, interval, QPDUNES_TRUE);
    }

    return QPDUNES_OK;
}


return_t qpDUNES_setupAllLocalQPs(qpData_t* const qpData) {
    size_t k;
    interval_t* interval;

    /* (1) set up initial lambda guess */
    qpDUNES_updateVector(
        &(qpData->intervals[0]->lambdaK1), &(qpData->lambda.data[0]), _NX_);
    for (k = 1u; k < _NI_; k++) {
        qpDUNES_updateVector(&(qpData->intervals[k]->lambdaK),
                             &(qpData->lambda.data[(k - 1u) * _NX_]), _NX_);
        qpDUNES_updateVector(&(qpData->intervals[k]->lambdaK1),
                             &(qpData->lambda.data[k * _NX_]), _NX_ );
    }
    qpDUNES_updateVector(&(qpData->intervals[_NI_]->lambdaK),
                         &(qpData->lambda.data[(_NI_ - 1u) * _NX_]), _NX_);

    /* (2) decide which QP solver to use and set up */
    for (k = 0; k < _NI_ + 1u; k++) {
        interval = qpData->intervals[k];

        /* (c) prepare stage QP solvers */
        qpDUNES_setupStageQP(qpData, interval, QPDUNES_TRUE);
    }

    return QPDUNES_OK;
}


return_t qpDUNES_setupStageQP(qpData_t* const qpData,
interval_t* const interval, boolean_t refactorHessian) {
    assert(qpData && interval);

    return_t statusFlag;

    /* (b) prepare clipping QP solver */
    if (refactorHessian == QPDUNES_TRUE) {
        /*
        Only first Hessian needs to be factorized in LTI case, others can
        be copied; last one might still be different, due to terminal
        cost, even in LTI case
        */
        factorizeH(qpData, &(interval->cholH), &(interval->H),
                   interval->nV);
    }

    /* (c) solve unconstrained local QP for g and initial lambda guess: */
    /*     - get (possibly updated) lambda guess */
    if (interval->id > 0) {     /* lambdaK exists */
        qpDUNES_updateVector(
            &(interval->lambdaK),
            &(qpData->lambda.data[(interval->id - 1u) * _NX_]),
            _NX_);
    }
    if (interval->id < _NI_) {      /* lambdaK1 exists */
        qpDUNES_updateVector(
            &(interval->lambdaK1),
            &(qpData->lambda.data[interval->id * _NX_]),
            _NX_);
    }

    interval->p = 0.0f;

    /*     - update first order term */
    /* reset q; qStep is added in qpDUNES_solve, when bounds are known */
    qpDUNES_setupZeroVector(&(interval->q), interval->nV);
    clippingQpSolver_updateStageData(qpData, interval,
                                     &(interval->lambdaK),
                                     &(interval->lambdaK1));
    /* Note: qStep is rewritten in line before */
    /*     - solve */
    statusFlag = directQpSolver_solveUnconstrained(
        qpData, interval, &(interval->qpSolverClipping.qStep));
    if (statusFlag != QPDUNES_OK) {
        return statusFlag;
    }

    /* reset zUnconstrained */
    qpDUNES_setupZeroVector(&(interval->qpSolverClipping.zUnconstrained),
                            interval->nV);

    return statusFlag;
}


return_t qpDUNES_shiftIntervals(qpData_t* const qpData) {
    size_t k;

    /** (1) Shift Interval pointers */
    /*  save pointer to first interval */
    interval_t* freeInterval = qpData->intervals[0];

    /*  shift all but the last interval (different size) left */
    for (k = 0; k < _NI_ - 1u; k++) {
        qpData->intervals[k] = qpData->intervals[k + 1u];
        qpData->intervals[k]->id = (uint_t)k; /* correct stage index */
    }
    /*  hang the free interval on the second but last position */
    qpData->intervals[_NI_ - 1u] = freeInterval;
    qpData->intervals[_NI_ - 1u]->id = _NI_ - 1u; /* correct stage index */

    /* update definedness of lambda parts */
    qpData->intervals[0]->lambdaK.isDefined = QPDUNES_FALSE;
    qpData->intervals[_NI_ - 1u]->lambdaK.isDefined = QPDUNES_TRUE;

    return QPDUNES_OK;
}


return_t qpDUNES_shiftLambda(qpData_t* const qpData) {
    assert(qpData && qpData->lambda.data);
    _nassert((size_t)qpData->lambda.data % 4 == 0);

    size_t k, i;

    for (k = 0; k < _NI_ - 1u; k++) {
        for (i = 0; i < _NX_; i++) {
            qpData->lambda.data[k * _NX_ + i] =
                qpData->lambda.data[(k + 1u) * _NX_ + i];
        }
    }

    return QPDUNES_OK;
}


qpOptions_t qpDUNES_setupDefaultOptions(void) {
    qpOptions_t options;

    /* iteration limits */
    options.maxIter                     = 100;
    options.maxNumLineSearchIterations  = 5;               /* 0.3^19 = 1e-10 */
    options.maxNumLineSearchRefinementIterations    = 10;   /* 0.62^49 = 1e-10 */

    /* printing */
    options.printLevel                  = 2;
    options.printIntervalHeader         = 20;
    options.printIterationTiming        = QPDUNES_FALSE;
    options.printLineSearchTiming       = QPDUNES_FALSE;

    /* logging */
    options.logLevel                    = QPDUNES_LOG_OFF;

    /* numerical tolerances */
    options.stationarityTolerance       = 1e-6f;
    options.equalityTolerance           = 2.221e-14f;
    options.newtonHessDiagRegTolerance  = 1e-10f;
    options.activenessTolerance         = 1e4f * options.equalityTolerance;
    options.QPDUNES_ZERO                = 1e-16f;
    options.QPDUNES_INFTY               = 1e16f;
    options.ascentCurvatureTolerance    = 1e-6f;

    /* additional options */
    options.nbrInitialGradientSteps     = 0;
    options.checkForInfeasibility       = QPDUNES_FALSE;

    /* regularization option */
    options.regType                     = QPDUNES_REG_LEVENBERG_MARQUARDT;
    options.regParam                    = 1e-6f;    /**< the regularization parameter added on singular Hessian elements
                                                         - should be quite a bit bigger than regularization tolerance
                                                         - assumption: if regularization needed, than Hessian has a singular direction
                                                         - in this singular direction i want to do mostly a gradient step,
                                                           few Hessian information usable
                                                      */

    options.nwtnHssnFacAlg              = QPDUNES_NH_FAC_BAND_REVERSE;


    /* line search options */
    options.lsType                          = QPDUNES_LS_ACCELERATED_GRADIENT_BISECTION_LS;
    options.lineSearchReductionFactor       = 0.1f;  /**< needs to be between 0 and 1 */
    options.lineSearchIncreaseFactor        = 1.5f;  /**< needs to be greater than 1 */
    options.lineSearchMinAbsProgress        = options.equalityTolerance;
    options.lineSearchMinRelProgress        = 1e-10f;
    options.lineSearchStationarityTolerance = 1e-6f;
    options.lineSearchMaxStepSize           = 0.25f;
    options.lineSearchNbrGridPoints         = 5;

    /* qpOASES options */
    options.qpOASES_terminationTolerance    = 1e-12f;   /*< stationarity tolerance for qpOASES, see qpOASES::Options -> terminationTolerance */

    return options;
}
