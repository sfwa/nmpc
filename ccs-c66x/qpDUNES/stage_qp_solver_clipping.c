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
#include <float.h>

#include "stage_qp_solver_clipping.h"
#include "../c66math.h"


/*
update QP data

qStep = C.T*lambdaK1 - [lambdaK.T 0]
pStep = c*lambdaK1
*/
return_t clippingQpSolver_updateStageData(qpData_t* const qpData,
interval_t* const interval, const z_vector_t* const lambdaK,
const z_vector_t* const lambdaK1) {
    assert(qpData && interval && lambdaK && lambdaK1);

    size_t i;

    /*
    WARNING: working with qStep might lead to high cancellation errors (due to
    subtraction, result might be almost zero) but we can hope for this to
    appear only when lambda steps are identical over iterations (-> infeasible
    problem, see journal paper). Always recomputing q is less efficient.
    */

    if (lambdaK1->isDefined == QPDUNES_TRUE) {
        /* qStep = C.T*lambdaK1 */
        multiplyCTy(qpData, &(interval->qpSolverClipping.qStep),
                    &(interval->C), lambdaK1);
        /* pStep = c*lambdaK1 -- constant objective term */
        interval->qpSolverClipping.pStep =
            scalarProd(lambdaK1, &(interval->c), _NX_);
    } else {
        /* qStep = 0 */
        qpDUNES_setupZeroVector(&(interval->qpSolverClipping.qStep),
                                interval->nV);
        /* pStep = 0 --  constant objective term */
        interval->qpSolverClipping.pStep = 0.0f;
    }

    if (lambdaK->isDefined == QPDUNES_TRUE) {
        /* qStep -= [lambdaK.T 0] */
        for (i = 0; i < _NX_; i++) {
            interval->qpSolverClipping.qStep.data[i] -= lambdaK->data[i];
        }
    }

    return QPDUNES_OK;
}



/* solve unconstrained QP */
return_t directQpSolver_solveUnconstrained(qpData_t* const qpData,
interval_t* const interval, const z_vector_t* const qStep) {
    assert(qpData && interval && qStep);

    return_t statusFlag;

    /* solve unconstraint QP */
    statusFlag = multiplyInvHz(qpData, &(interval->qpSolverClipping.dz),
                               &(interval->cholH), qStep, interval->nV);
    negateVector(&(interval->qpSolverClipping.dz), interval->nV);

    return statusFlag;
}


/*
gets the step size to the first active set change if it is shorter than an
incumbent step size initially in alphaMin

compare gaps from z (i.e., multipliers from previous iteration) with dz
*/
return_t directQpSolver_getMinStepsize(const interval_t* const interval,
real_t* alphaMin) {
    assert(interval && alphaMin);

    size_t i;
    real_t alphaASChange, tmp1, tmp2;

    for (i = 0; i < interval->nV; i++) {
        /*
        The original qpDUNES implementation assumes compiler/chip support for
        1./0. == inf, and (2. < inf) == TRUE

        We can't really make that guarantee on CCS/C66x so encode it
        explicitly instead.
        */
        if (abs_f(interval->y.data[2u * i]) > 1e-15) {
            tmp1 = divide_f(interval->qpSolverClipping.dz.data[i],
                            interval->y.data[2u * i]);
        } else {
            tmp1 = FLT_MAX;
        }
        if (abs_f(interval->y.data[2 * i + 1u]) > 1e-15) {
            tmp2 = divide_f(interval->qpSolverClipping.dz.data[i],
                            -interval->y.data[2 * i + 1u]);
        } else {
            tmp2 = FLT_MAX;
        }

        alphaASChange = recip_f(qpDUNES_fmax(tmp1, tmp2));

        if (alphaASChange > 1e-15 && alphaASChange < *alphaMin) {
            *alphaMin = alphaASChange;
        }
    }

    return QPDUNES_OK;
}


/* do a step of length alpha */
return_t directQpSolver_doStep(qpData_t* const qpData,
interval_t* const interval, const z_vector_t* const stepDir, real_t alpha,
z_vector_t* const zUnconstrained, z_vector_t* const z, d2_vector_t* const mu,
z_vector_t* const q, real_t* const p) {
    assert(qpData && interval && stepDir && zUnconstrained && z && mu && q);

    size_t i;

    /* update primal solution and get dual solution */
    addVectorScaledVector(zUnconstrained,
                          &(interval->qpSolverClipping.zUnconstrained), alpha,
                          stepDir, interval->nV);

    /*
    skip copying if zUnconstrained and z are pointing to the same object
    (e.g., during line search, when just trying steps)
    */
    if (z != zUnconstrained) {
        qpDUNES_copyVector(z, zUnconstrained, interval->nV);
    }
    directQpSolver_saturateVector(qpData, z, mu, &(interval->zLow),
                                  &(interval->zUpp), interval->nV);

    /* update q */
    for (i = 0; i < interval->nV; i++) {
        q->data[i] = interval->q.data[i] +
                     alpha * interval->qpSolverClipping.qStep.data[i];
    }

    /* update p */
    *p = interval->p + alpha * interval->qpSolverClipping.pStep;

    return QPDUNES_OK;
}


return_t directQpSolver_saturateVector(qpData_t* const qpData,
d_vector_t* const vec,
d2_vector_t* const mu, /* pseudo multipliers, resembling the gaps to the bounds; + active, - inactive */
const d_vector_t* const lb, const d_vector_t* const ub, size_t nV) {
    assert(qpData && vec && mu && lb && ub);

    size_t i;

    for (i = 0; i < nV; i++) {
        /* feasibility gap to lower bound; negative value means inactive */
        mu->data[2u * i] = lb->data[i] - vec->data[i];
        /* feasibility gap to upper bound; negative value means inactive */
        mu->data[2u * i + 1u] = vec->data[i] - ub->data[i];
        if (mu->data[2u * i] >= -qpData->options.activenessTolerance) {
            /*
            TODO: check if this implementation of activeness Tolerance makes
            sense
            */
            vec->data[i] = lb->data[i];
        } else if (mu->data[2u * i + 1u] >=
                    -qpData->options.activenessTolerance) {
            vec->data[i] = ub->data[i];
        }
    }

    return QPDUNES_OK;
}
