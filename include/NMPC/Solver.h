/*
Part of NMPC. Distributed under the MIT licence: see LICENSE in the root of
the repository.
*/

#ifndef NMPC_SOLVER_H
#define NMPC_SOLVER_H

#include "NMPC/BoxQP.h"
#include "NMPC/Options.h"
#include "NMPC/Problem.h"
#include "NMPC/Result.h"
#include "NMPC/Workspace.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>

namespace NMPC {

/*
Box-constrained Gauss-Newton iLQR with a feasible rollout.

A solve linearises the problem about the plan, runs one backward pass that
solves a box QP at each stage, and one forward pass with a line search on the
step length. With the default options that is the whole solve: one iteration
from the previous plan shifted by a stage, which is a real-time iteration.

    solver.reset(x0, u0)                    a first plan: the command u0 at every stage
    solver.advance(problem, 1)              shift the plan one stage, once per solve period
    result = solver.solve(problem, x, options)
    solver.policy()                         the plan's states, commands and gains

The pieces that matter, each learned from a failed flight or a review finding:

  - Analytic step and residual Jacobians. Finite differences belong in tests (Check.h).
  - The backward pass forms one triangle of each symmetric product. It is the
    largest single part of a solve on a Cortex-M7.
  - Bounds through a box QP at each stage. The gain is zero in the direction of a
    command held at a bound.
  - Regularisation starts afresh at every solve.
  - The line search accepts a step against both terms of the predicted decrease.
  - A solve stops when the predicted decrease is below what the type resolves, and at once
    when it is zero.
  - Each stage's Jacobians are kept for the iteration by default; Solver<Problem, N,
    Jacobians::PerStage> makes them in the backward pass instead, for a smaller workspace.
  - A solve whose rollout is not finite leaves the plan as it was.

Two things that belong to the application are worth stating here. First, the
plan comes into force one solve period after the solve starts, so solve from
the state predicted for that moment, not the state measured when it starts.
At periods longer than about 20 ms that prediction must follow the stored
gains as well as the plan's commands. Second, if a solve is not usable, the
plan in force flies on into its next stage and the next solve should start
from a reset.
*/
template <class Problem, int N = Problem::Stages, Jacobians Storage = Jacobians::Stored>
class Solver {
public:
    using Traits = ProblemTraits<Problem>;
    using T = typename Traits::T;
    using State = typename Traits::State;
    using Control = typename Traits::Control;
    using Tangent = typename Traits::Tangent;
    using MatA = typename Traits::MatA;
    using MatB = typename Traits::MatB;
    using Gain = typename Traits::Gain;
    using MatU = typename Traits::MatU;
    using Workspace = NMPC::Workspace<Problem, N, Storage>;
    using Policy = NMPC::Policy<Problem, N>;

    enum { NX = Traits::NX, NU = Traits::NU, NR = Traits::NR, NT = Traits::NT, Stages = N };

    static_assert(detail::HasCostExpansion<Problem>::value || NR > 0,
        "A problem needs stage residuals or its own cost_expansion");

    explicit Solver(Workspace& workspace) : ws(workspace) {}

    /* A first plan: the command u0 at every stage from the state x0, with no gains. */
    void reset(const State& x0, const Control& u0) {
        for (int k = 0; k < N; k++) {
            ws.plan.x[k] = x0;
            ws.plan.u[k] = u0;
            ws.plan.K[k].setZero();
            ws.step[k].setZero();
        }
        ws.plan.x[N] = x0;
        ws.solved = false;
    }

    /*
    Move the plan on by a number of whole stages, once per solve period. A stage that enters at the
    end takes the problem's tail_control with no gain where the problem has one, and repeats the
    last stage's command and gain where it does not. Only the commands and gains matter to the
    next solve, which starts with a rollout; the states are shifted so that the plan stays readable.

    The stage QPs' steps, which start the next solve's stage QPs, stay with their stage numbers
    and are not shifted. That is what the review prototype did, and the reference flights were
    flown with it. The prototype's stage QP stopped on an absolute test of its gradient, so where
    a stage's Hessian was small its result depended on where it started, by up to 2e-4 N in a
    rotor command on the multirotor. The stop is now relative (Options::box_tolerance), which
    leaves only rounding to depend on the start.
    */
    void advance(const Problem& problem, int stages) {
        stages = std::min(std::max(stages, 0), N);
        if (stages == 0) {
            return;
        }
        /* The last stage's command and gain, kept before the shift, which overwrites them when every stage moves. */
        const Control last_u = ws.plan.u[N - 1];
        const Gain last_K = ws.plan.K[N - 1];
        for (int k = 0; k + stages < N; k++) {
            ws.plan.x[k] = ws.plan.x[k + stages];
            ws.plan.u[k] = ws.plan.u[k + stages];
            ws.plan.K[k] = ws.plan.K[k + stages];
        }
        for (int k = N - stages; k < N; k++) {
            ws.plan.x[k] = ws.plan.x[N];
            if constexpr (detail::HasTailControl<Problem>::value) {
                problem.tail_control(k, ws.plan.u[k]);
                ws.plan.K[k].setZero();
            } else {
                (void)problem;
                ws.plan.u[k] = last_u;
                ws.plan.K[k] = last_K;
            }
        }
    }

    /*
    For a solve period shorter than a stage: each stage's command moves that fraction of the way to the next
    stage's, and the states and gains stay. The next rollout makes the states again.
    */
    void advance_fraction(T fraction) {
        for (int k = 0; k + 1 < N; k++) {
            ws.plan.u[k] += fraction * (ws.plan.u[k + 1] - ws.plan.u[k]);
        }
    }

    Result<T> solve(const Problem& problem, const State& x0, const Options<T>& options) {
        Result<T> result;
        Dx = detail::state_scale(problem);
        Du = detail::control_scale(problem);
        inverse_Dx = Dx.cwiseInverse();
        scale_products();
        for (int k = 0; k < N; k++) {
            problem.control_bounds(k, ws.lower[k], ws.upper[k]);
        }
        regularisation = options.regularisation_initial;

        /* The previous plan's commands, rolled out from the new state. Nothing changes if that fails. */
        int bound_stages = 0;
        T cost = rollout(problem, x0, bound_stages);
        if (!std::isfinite(cost)) {
            result.regularisation = regularisation;
            return result;
        }
        accept();
        ws.solved = true;
        result.status = Status::IterationLimit;
        result.cost_before = cost;
        result.cost_after = cost;
        result.stages_on_bound = bound_stages;

        for (int iteration = 0; iteration < options.max_iterations; iteration++) {
            result.iterations = iteration + 1;
            linearise(problem);

            bool accepted = false;
            bool converged = false;
            for (int retry = 0; retry <= options.max_regularisation_retries && !accepted; retry++) {
                T dV1 = T(0);
                T dV2 = T(0);
                result.backward_passes++;
                if (!backward(problem, options, dV1, dV2, result.box_iterations)) {
                    regularisation = std::min(regularisation * options.regularisation_factor,
                        options.regularisation_max);
                    continue;
                }
                result.predicted_decrease = -(dV1 + dV2);
                /*
                A predicted decrease of zero (every stage's step is zero, as at an exact optimum of zero cost)
                is convergence whatever the floor: no step can be accepted, and retries and a line search
                would only spend time to report NoProgress.
                */
                if (!(result.predicted_decrease > T(0)) || (options.precision_floor > T(0) &&
                        result.predicted_decrease < options.precision_floor * std::abs(cost))) {
                    converged = true;
                    break;
                }
                T alpha = T(1);
                for (int trial = 0; trial < options.max_line_search; trial++, alpha *= T(0.5)) {
                    result.line_search_trials++;
                    const T expected = -(alpha * dV1 + alpha * alpha * dV2);
                    const T trial_cost = forward(problem, x0, alpha, bound_stages);
                    if (std::isfinite(trial_cost) && expected > T(0) &&
                            (cost - trial_cost) > options.accept_ratio * expected) {
                        accept();
                        cost = trial_cost;
                        result.stages_on_bound = bound_stages;
                        accepted = true;
                        break;
                    }
                }
                if (accepted) {
                    regularisation = std::max(regularisation / options.regularisation_factor,
                        options.regularisation_min);
                } else {
                    regularisation = std::min(regularisation * options.regularisation_factor,
                        options.regularisation_max);
                }
            }
            if (converged) {
                result.status = Status::Converged;
                break;
            }
            if (!accepted) {
                result.status = Status::NoProgress;
                break;
            }
            result.accepted_steps++;
            result.cost_after = cost;
        }
        result.regularisation = regularisation;
        return result;
    }

    const Policy& policy() const { return ws.plan; }
    const Workspace& workspace() const { return ws; }

    /*
    The parts of a solve, public so that a benchmark can time each one. They work on the plan in the
    workspace with the bounds and scales of the last solve.
    */

    /* Roll the plan's commands, clamped to the bounds, out from x0 into the trial plan. Returns its cost. */
    T rollout(const Problem& problem, const State& x0, int& bound_stages) {
        bound_stages = 0;
        ws.trial_x[0] = x0;
        T cost = T(0);
        for (int k = 0; k < N; k++) {
            const Control u = ws.plan.u[k].cwiseMax(ws.lower[k]).cwiseMin(ws.upper[k]);
            ws.trial_u[k] = u;
            bound_stages += on_bound(k, u);
            cost += detail::stage_cost(problem, k, ws.trial_x[k], u);
            ws.trial_x[k + 1] = problem.step(k, ws.trial_x[k], u);
            if (!detail::finite<Problem>(ws.trial_x[k + 1])) {
                return std::numeric_limits<T>::quiet_NaN();
            }
        }
        return cost + detail::terminal_cost(problem, ws.trial_x[N]);
    }

    /* The step's Jacobians at every stage of the plan, normalised. With Jacobians::PerStage the backward pass makes them instead. */
    void linearise(const Problem& problem) {
        if constexpr (Storage == Jacobians::Stored) {
            for (int k = 0; k < N; k++) {
                linearise_stage(problem, k, ws.jacobians.A[k], ws.jacobians.B[k]);
            }
        } else {
            (void)problem;
        }
    }

    /* The step's Jacobians at stage k of the plan, normalised. */
    void linearise_stage(const Problem& problem, int k, MatA& A, MatB& B) const {
        problem.step_jacobian(k, ws.plan.x[k], ws.plan.u[k], A, B);
        A = A.cwiseProduct(scale_A);
        B = B.cwiseProduct(scale_B);
    }

    /*
    The backward pass, with one triangle of each symmetric product. It stores each stage's step in
    the workspace and its gain in the plan, and returns false if a stage QP could not be factored.
    dV1 and dV2 are the linear and quadratic terms of the predicted change of cost for a full step.
    */
    bool backward(const Problem& problem, const Options<T>& options, T& dV1, T& dV2, int& box_iterations) {
        Tangent Vx;
        MatA Vxx;
        terminal_expansion(problem, Vx, Vxx);
        for (int k = N - 1; k >= 0; k--) {
            Linearisation here;
            if constexpr (Storage == Jacobians::PerStage) {
                linearise_stage(problem, k, here.A, here.B);
            }
            const MatA& A = stage_A(here, k);
            const MatB& B = stage_B(here, k);
            Tangent Qx;
            Control Qu;
            MatA Qxx;
            Gain Qux;
            MatU Quu;
            stage_expansion(problem, k, Qx, Qu, Qxx, Qux, Quu);

            MatA VA;
            MatB VB;
            for (int j = 0; j < NX; j++) {
                for (int i = 0; i < NX; i++) {
                    VA(i, j) = Vxx.col(i).dot(A.col(j));
                }
            }
            for (int j = 0; j < NU; j++) {
                for (int i = 0; i < NX; i++) {
                    VB(i, j) = Vxx.col(i).dot(B.col(j));
                }
            }
            for (int i = 0; i < NX; i++) {
                Qx(i) += A.col(i).dot(Vx);
            }
            for (int i = 0; i < NU; i++) {
                Qu(i) += B.col(i).dot(Vx);
            }
            for (int j = 0; j < NX; j++) {
                for (int i = 0; i <= j; i++) {
                    Qxx(i, j) += A.col(i).dot(VA.col(j));
                }
                for (int i = 0; i < NU; i++) {
                    Qux(i, j) += B.col(i).dot(VA.col(j));
                }
            }
            for (int j = 0; j < NU; j++) {
                for (int i = 0; i <= j; i++) {
                    Quu(i, j) += B.col(i).dot(VB.col(j));
                    Quu(j, i) = Quu(i, j);
                }
            }

            MatU H = Quu;
            H.diagonal().array() += regularisation;
            const Control lo = (ws.lower[k] - ws.plan.u[k]).cwiseQuotient(Du);
            const Control hi = (ws.upper[k] - ws.plan.u[k]).cwiseQuotient(Du);
            Control d = ws.step[k];
            Eigen::Array<bool, NU, 1> held;
            Eigen::LLT<MatU> factor;
            const int passes = solve_box_qp<T, NU>(H, Qu, lo, hi, d, held, factor, options.max_box_iterations,
                options.box_tolerance);
            if (passes < 0) {
                return false;
            }
            box_iterations += passes;

            Gain rhs = -Qux;
            for (int i = 0; i < NU; i++) {
                if (held(i)) {
                    rhs.row(i).setZero();
                }
            }
            /*
            The gain, column by column, by forward and back substitution with the stage QP's factor. Eigen 3.2.8's
            solve with a matrix right-hand side goes through its blocked kernels, which at these sizes cost a
            Cortex-M7 about 5 % more of the backward pass and 0.7 kB more stack; its own vector solve is slower
            still. The answer is the same to rounding. A held command's row of the factor is the identity and
            its right-hand side is zero, so its gain is zero.
            */
            Gain K;
            {
                const MatU& L = factor.matrixLLT();
                Control inverse_diagonal;
                for (int i = 0; i < NU; i++) {
                    inverse_diagonal(i) = T(1) / L(i, i);
                }
                for (int j = 0; j < NX; j++) {
                    Control c(rhs.col(j));
                    for (int i = 0; i < NU; i++) {
                        T v = c(i);
                        for (int m = 0; m < i; m++) {
                            v -= L(i, m) * c(m);
                        }
                        c(i) = v * inverse_diagonal(i);
                    }
                    for (int i = NU - 1; i >= 0; i--) {
                        T v = c(i);
                        for (int m = i + 1; m < NU; m++) {
                            v -= L(m, i) * c(m);
                        }
                        c(i) = v * inverse_diagonal(i);
                    }
                    K.col(j) = c;
                }
            }
            for (int i = 0; i < NU; i++) {
                if (held(i)) {
                    K.row(i).setZero();
                }
            }
            if (!K.allFinite() || !d.allFinite()) {
                return false;
            }
            ws.step[k] = d;
            for (int j = 0; j < NX; j++) {
                for (int i = 0; i < NU; i++) {
                    ws.plan.K[k](i, j) = Du(i) * K(i, j) * inverse_Dx(j);
                }
            }

            const Control Quud(Quu * d);
            dV1 += d.dot(Qu);
            dV2 += T(0.5) * d.dot(Quud);
            const Control t = Qu + Quud;
            for (int i = 0; i < NX; i++) {
                Vx(i) = Qx(i) + K.col(i).dot(t) + Qux.col(i).dot(d);
            }
            Gain M;
            for (int j = 0; j < NX; j++) {
                M.col(j) = Qux.col(j) + Quu * K.col(j);
            }
            for (int j = 0; j < NX; j++) {
                for (int i = 0; i <= j; i++) {
                    const T v = Qxx(i, j) + Qux.col(i).dot(K.col(j)) + K.col(i).dot(M.col(j));
                    Vxx(i, j) = v;
                    Vxx(j, i) = v;
                }
            }
        }
        return true;
    }

    /* A trial plan from x0 with a step of length alpha and the new gains. Returns its cost. */
    T forward(const Problem& problem, const State& x0, T alpha, int& bound_stages) {
        bound_stages = 0;
        ws.trial_x[0] = x0;
        T cost = T(0);
        for (int k = 0; k < N; k++) {
            const Tangent deviation = Problem::difference(ws.trial_x[k], ws.plan.x[k]);
            const Control u = (ws.plan.u[k] + alpha * Du.cwiseProduct(ws.step[k]) + ws.plan.K[k] * deviation)
                .cwiseMax(ws.lower[k]).cwiseMin(ws.upper[k]);
            ws.trial_u[k] = u;
            bound_stages += on_bound(k, u);
            cost += detail::stage_cost(problem, k, ws.trial_x[k], u);
            ws.trial_x[k + 1] = problem.step(k, ws.trial_x[k], u);
            if (!detail::finite<Problem>(ws.trial_x[k + 1])) {
                return std::numeric_limits<T>::quiet_NaN();
            }
        }
        return cost + detail::terminal_cost(problem, ws.trial_x[N]);
    }

    /* The regularisation in force, in normalised units. */
    T regularisation = T(0);

private:
    Workspace& ws;
    Tangent Dx = Tangent::Ones();
    Tangent inverse_Dx = Tangent::Ones();
    Control Du = Control::Ones();
    /*
    The products of the scales that normalise the Jacobians and the cost expansion, made once a solve: scale_A(i, j)
    is inverse_Dx(i) Dx(j), and so on. Each normalised coefficient is the same product as before, in the same order.
    */
    MatA scale_A, scale_xx;
    MatB scale_B;
    Gain scale_ux;
    MatU scale_uu;

    void scale_products() {
        for (int j = 0; j < NX; j++) {
            for (int i = 0; i < NX; i++) {
                scale_A(i, j) = inverse_Dx(i) * Dx(j);
                scale_xx(i, j) = Dx(i) * Dx(j);
            }
            for (int i = 0; i < NU; i++) {
                scale_ux(i, j) = Du(i) * Dx(j);
            }
        }
        for (int j = 0; j < NU; j++) {
            for (int i = 0; i < NX; i++) {
                scale_B(i, j) = inverse_Dx(i) * Du(j);
            }
            for (int i = 0; i < NU; i++) {
                scale_uu(i, j) = Du(i) * Du(j);
            }
        }
    }

    /* A stage's Jacobians, made in the backward pass, with Jacobians::PerStage; nothing with Jacobians::Stored. */
    struct Empty {};
    struct Made {
        MatA A;
        MatB B;
    };
    using Linearisation = typename std::conditional<Storage == Jacobians::PerStage, Made, Empty>::type;

    const MatA& stage_A(const Linearisation& here, int k) const {
        if constexpr (Storage == Jacobians::PerStage) {
            (void)k;
            return here.A;
        } else {
            (void)here;
            return ws.jacobians.A[k];
        }
    }
    const MatB& stage_B(const Linearisation& here, int k) const {
        if constexpr (Storage == Jacobians::PerStage) {
            (void)k;
            return here.B;
        } else {
            (void)here;
            return ws.jacobians.B[k];
        }
    }

    void accept() {
        for (int k = 0; k < N; k++) {
            ws.plan.x[k] = ws.trial_x[k];
            ws.plan.u[k] = ws.trial_u[k];
        }
        ws.plan.x[N] = ws.trial_x[N];
    }

    int on_bound(int k, const Control& u) const {
        return ((u.array() <= ws.lower[k].array()) || (u.array() >= ws.upper[k].array())).any() ? 1 : 0;
    }

    /* The Gauss-Newton expansion of stage k's cost about the plan, normalised. Lxx holds its upper triangle. */
    void stage_expansion(const Problem& problem, int k, Tangent& lx, Control& lu, MatA& Lxx, Gain& Lux,
            MatU& Luu) const {
        const State& x = ws.plan.x[k];
        const Control& u = ws.plan.u[k];
        if constexpr (detail::HasCostExpansion<Problem>::value) {
            problem.cost_expansion(k, x, u, lx, lu, Lxx, Lux, Luu);
            lx = lx.cwiseProduct(Dx);
            lu = lu.cwiseProduct(Du);
            /* Only the upper triangle of Lxx is the problem's to fill, and read. */
            Lxx.template triangularView<Eigen::Upper>() = Lxx.cwiseProduct(scale_xx);
            Lux = Lux.cwiseProduct(scale_ux);
            Luu = Luu.cwiseProduct(scale_uu);
        } else {
            typename Traits::StageResidual r;
            typename Traits::StageJx Jx;
            typename Traits::StageJu Ju;
            problem.stage_residual(k, x, u, r);
            problem.residual_jacobians(k, x, u, Jx, Ju);
            for (int j = 0; j < NX; j++) {
                Jx.col(j) *= Dx(j);
            }
            for (int j = 0; j < NU; j++) {
                Ju.col(j) *= Du(j);
            }
            for (int j = 0; j < NX; j++) {
                lx(j) = Jx.col(j).dot(r);
                for (int i = 0; i <= j; i++) {
                    Lxx(i, j) = Jx.col(i).dot(Jx.col(j));
                }
                for (int i = 0; i < NU; i++) {
                    Lux(i, j) = Ju.col(i).dot(Jx.col(j));
                }
            }
            for (int j = 0; j < NU; j++) {
                lu(j) = Ju.col(j).dot(r);
                for (int i = 0; i <= j; i++) {
                    Luu(i, j) = Ju.col(i).dot(Ju.col(j));
                    Luu(j, i) = Luu(i, j);
                }
            }
            if constexpr (detail::HasLinearCost<Problem>::value) {
                Tangent gx;
                Control gu;
                problem.linear_cost_gradient(k, x, u, gx, gu);
                lx += gx.cwiseProduct(Dx);
                lu += gu.cwiseProduct(Du);
            }
        }
    }

    /* The terminal cost's expansion at the end of the plan, normalised, both triangles. */
    void terminal_expansion(const Problem& problem, Tangent& Vx, MatA& Vxx) const {
        const State& x = ws.plan.x[N];
        if constexpr (detail::HasTerminalExpansion<Problem>::value) {
            problem.terminal_expansion(x, Vx, Vxx);
            Vx = Vx.cwiseProduct(Dx);
            for (int j = 0; j < NX; j++) {
                for (int i = 0; i <= j; i++) {
                    Vxx(i, j) *= Dx(i) * Dx(j);
                    Vxx(j, i) = Vxx(i, j);
                }
            }
        } else if constexpr (NT > 0) {
            typename Traits::TerminalResidual r;
            typename Traits::TerminalJx J;
            problem.terminal_residual(x, r);
            problem.terminal_jacobian(x, J);
            for (int j = 0; j < NX; j++) {
                J.col(j) *= Dx(j);
            }
            for (int j = 0; j < NX; j++) {
                Vx(j) = J.col(j).dot(r);
                for (int i = 0; i <= j; i++) {
                    Vxx(i, j) = J.col(i).dot(J.col(j));
                    Vxx(j, i) = Vxx(i, j);
                }
            }
        } else {
            (void)problem;
            Vx.setZero();
            Vxx.setZero();
        }
    }
};

}

#endif
