#include "Typed.h"
#include <cmath>
#include <limits>
#include "NMPC/Solver.h"
#include "Problems.h"

template <typename T>
class SolverTest : public ::testing::Test {};

typedef ::testing::Types<float, double> Precisions;
TYPED_TEST_SUITE(SolverTest, Precisions);

namespace {

template <class T> T loose() { return std::is_same<T, float>::value ? T(2e-3) : T(1e-8); }

/*
The linear problem condensed into one least-squares problem in the stacked commands U: the
residuals are G U + h. Its minimiser is the plan that one Gauss-Newton iteration must find.
*/
struct Condensed {
    enum { N = 8, NX = 3, NU = 2, Rows = 5 * N + 3, Vars = NU * N };
    Eigen::Matrix<double, Rows, Vars> G;
    Eigen::Matrix<double, Rows, 1> h;

    explicit Condensed(const Linear_Problem<double>& p, const Eigen::Matrix<double, 3, 1>& x0) {
        G.setZero();
        h.setZero();
        Eigen::Matrix<double, NX, NX> power = Eigen::Matrix<double, NX, NX>::Identity();
        Eigen::Matrix<double, NX, Vars> map = Eigen::Matrix<double, NX, Vars>::Zero();
        for (int k = 0; k <= N; k++) {
            const Eigen::Matrix<double, NX, 1> free_state(power * x0);
            if (k < N) {
                G.block<5, Vars>(5 * k, 0) = p.Wx * map;
                G.block<5, NU>(5 * k, NU * k) += p.Wu;
                h.segment<5>(5 * k) = p.Wx * (free_state - p.x_target);
                map = (p.A * map).eval();
                map.block<NX, NU>(0, NU * k) += p.B;
                power = (p.A * power).eval();
            } else {
                G.block<3, Vars>(5 * N, 0) = p.Wt * map;
                h.segment<3>(5 * N) = p.Wt * (free_state - p.x_target);
            }
        }
    }

    Eigen::Matrix<double, Vars, 1> unbounded() const {
        const Eigen::Matrix<double, Vars, Vars> H(G.transpose() * G);
        return H.llt().solve(-(G.transpose() * h));
    }

    Eigen::Matrix<double, Vars, 1> bounded(const Eigen::Matrix<double, NU, 1>& lower,
            const Eigen::Matrix<double, NU, 1>& upper) const {
        const Eigen::Matrix<double, Vars, Vars> H(G.transpose() * G);
        const Eigen::Matrix<double, Vars, 1> g(G.transpose() * h);
        Eigen::Matrix<double, Vars, 1> lo, hi, d = Eigen::Matrix<double, Vars, 1>::Zero();
        for (int k = 0; k < N; k++) {
            lo.segment<NU>(NU * k) = lower;
            hi.segment<NU>(NU * k) = upper;
        }
        Eigen::Array<bool, Vars, 1> held;
        Eigen::LLT<Eigen::Matrix<double, Vars, Vars>> factor;
        NMPC::solve_box_qp<double, Vars>(H, g, lo, hi, d, held, factor, 500);
        return d;
    }

    double cost(const Eigen::Matrix<double, Vars, 1>& U) const { return 0.5 * (G * U + h).squaredNorm(); }
};

}

TYPED_TEST(SolverTest, OneIterationOnALinearProblemIsTheRiccatiSolution) {
    typedef TypeParam T;
    typedef Linear_Problem<T> Problem;
    typedef NMPC::Solver<Problem> Solver;
    static typename Solver::Workspace workspace;
    Solver solver(workspace);
    Problem problem;
    const typename Problem::State x0(T(0.3), T(-0.2), T(0.8));

    NMPC::Options<T> options;
    options.regularisation_initial = T(0);
    options.precision_floor = T(0);
    solver.reset(x0, Problem::Control::Zero());
    const NMPC::Result<T> result = solver.solve(problem, x0, options);
    ASSERT_EQ(result.status, NMPC::Status::IterationLimit);
    EXPECT_EQ(result.accepted_steps, 1);
    EXPECT_EQ(result.line_search_trials, 1);

    const Condensed reference(Linear_Problem<double>(), x0.template cast<double>());
    const auto U = reference.unbounded();
    for (int k = 0; k < Condensed::N; k++) {
        for (int i = 0; i < 2; i++) {
            EXPECT_NEAR(double(solver.policy().u[k](i)), U(2 * k + i), 1e2 * double(loose<T>()))
                << "stage " << k << ", command " << i;
        }
    }
    EXPECT_NEAR(double(result.cost_after), reference.cost(U), double(loose<T>()));
    /* The predicted decrease is exact on a linear-quadratic problem. */
    EXPECT_NEAR(double(result.cost_before - result.cost_after), double(result.predicted_decrease),
        10.0 * double(loose<T>()));
}

TYPED_TEST(SolverTest, BoundsMatchAConvergedReference) {
    typedef TypeParam T;
    typedef Linear_Problem<T> Problem;
    typedef NMPC::Solver<Problem> Solver;
    static typename Solver::Workspace workspace;
    Solver solver(workspace);
    Problem problem;
    problem.lower_bound << T(-0.4), T(-0.2);
    problem.upper_bound << T(0.5), T(0.1);
    const typename Problem::State x0(T(0.3), T(-0.2), T(0.8));

    NMPC::Options<T> options;
    options.max_iterations = 200;
    options.max_line_search = 20;
    solver.reset(x0, Problem::Control::Zero());
    const NMPC::Result<T> result = solver.solve(problem, x0, options);
    ASSERT_TRUE(result.usable());
    EXPECT_GT(result.stages_on_bound, 0);

    Linear_Problem<double> problem_double;
    const Condensed reference(problem_double, x0.template cast<double>());
    const auto U = reference.bounded(problem.lower_bound.template cast<double>(),
        problem.upper_bound.template cast<double>());
    for (int k = 0; k < Condensed::N; k++) {
        for (int i = 0; i < 2; i++) {
            EXPECT_NEAR(double(solver.policy().u[k](i)), U(2 * k + i), 1e2 * double(loose<T>()))
                << "stage " << k << ", command " << i;
            EXPECT_GE(solver.policy().u[k](i), problem.lower_bound(i));
            EXPECT_LE(solver.policy().u[k](i), problem.upper_bound(i));
            /* No feedback in the direction of a held command. */
            if (solver.policy().u[k](i) == problem.lower_bound(i) || solver.policy().u[k](i) == problem.upper_bound(i)) {
                EXPECT_EQ(solver.policy().K[k].row(i).norm(), T(0));
            }
        }
    }
    EXPECT_NEAR(double(result.cost_after), reference.cost(U), 10.0 * double(loose<T>()));
}

TYPED_TEST(SolverTest, StopProblemRespectsItsLimitAndComesToRest) {
    typedef TypeParam T;
    typedef Stop_Problem<T> Problem;
    typedef NMPC::Solver<Problem> Solver;
    static typename Solver::Workspace workspace;
    Solver solver(workspace);
    Problem problem;
    NMPC::Options<T> options;

    typename Problem::State x(T(1), T(0));
    solver.reset(x, Problem::Control::Zero());
    /* Closed loop: one iteration a solve, the plan shifted one stage, the plant stepped with the first command. */
    for (int tick = 0; tick < 200; tick++) {
        if (tick > 0) {
            solver.advance(problem, 1);
        }
        const NMPC::Result<T> result = solver.solve(problem, x, options);
        ASSERT_TRUE(result.usable());
        const T u = solver.policy().u[0](0);
        EXPECT_LE(std::abs(u), problem.limit);
        typename Problem::Control command;
        command(0) = u;
        x = problem.step(0, x, command);
    }
    EXPECT_NEAR(double(x(0)), 0.0, 1e-3);
    EXPECT_NEAR(double(x(1)), 0.0, 1e-3);
}

TYPED_TEST(SolverTest, AttitudeOnTheManifoldConverges) {
    typedef TypeParam T;
    typedef Attitude_Problem<T> Problem;
    typedef NMPC::Solver<Problem> Solver;
    typedef Eigen::Matrix<T, 3, 1> V3;
    static typename Solver::Workspace workspace;
    Solver solver(workspace);
    Problem problem;
    NMPC::Options<T> options;

    typename Problem::State x;
    x.q = NMPC::so3_exp<T>(V3(T(0.6), T(-0.4), T(1.2)));
    x.w = V3(T(0.2), T(0), T(-0.3));
    solver.reset(x, Problem::Control::Zero());
    T first_cost = T(0);
    for (int tick = 0; tick < 150; tick++) {
        if (tick > 0) {
            solver.advance(problem, 1);
        }
        const NMPC::Result<T> result = solver.solve(problem, x, options);
        ASSERT_TRUE(result.usable());
        if (tick == 0) {
            first_cost = result.cost_before;
        }
        x = problem.step(0, x, solver.policy().u[0]);
    }
    EXPECT_LT(NMPC::attitude_difference<T>(x.q, problem.target).norm(), T(0.01));
    EXPECT_LT(x.w.norm(), T(0.01));
    EXPECT_GT(first_cost, T(0));
}

TYPED_TEST(SolverTest, AFailedSolveLeavesThePlan) {
    typedef TypeParam T;
    typedef Failing_Problem<T> Problem;
    typedef NMPC::Solver<Problem> Solver;
    static typename Solver::Workspace workspace;
    Solver solver(workspace);
    Problem problem;
    NMPC::Options<T> options;
    const typename Problem::State x0(T(1), T(0.5));
    solver.reset(x0, Problem::Control::Zero());
    ASSERT_TRUE(solver.solve(problem, x0, options).usable());

    typename Solver::Policy before = solver.policy();
    problem.fail = true;
    const NMPC::Result<T> result = solver.solve(problem, x0, options);
    EXPECT_FALSE(result.usable());
    EXPECT_EQ(result.status, NMPC::Status::NumericalFailure);
    for (int k = 0; k < Solver::Stages; k++) {
        EXPECT_EQ(solver.policy().x[k], before.x[k]);
        EXPECT_EQ(solver.policy().u[k], before.u[k]);
        EXPECT_EQ(solver.policy().K[k], before.K[k]);
    }
    EXPECT_EQ(solver.policy().x[Solver::Stages], before.x[Solver::Stages]);
}

TYPED_TEST(SolverTest, StopsAtThePrecisionFloor) {
    typedef TypeParam T;
    typedef Stop_Problem<T> Problem;
    typedef NMPC::Solver<Problem> Solver;
    static typename Solver::Workspace workspace;
    Solver solver(workspace);
    Problem problem;
    NMPC::Options<T> options;
    options.max_iterations = 100;
    const typename Problem::State x0(T(0.5), T(0));
    solver.reset(x0, Problem::Control::Zero());
    const NMPC::Result<T> result = solver.solve(problem, x0, options);
    EXPECT_EQ(result.status, NMPC::Status::Converged);
    EXPECT_LT(result.iterations, 100);
    /* Solved again from the same state, the plan is already as good as the type can tell. */
    const NMPC::Result<T> again = solver.solve(problem, x0, options);
    EXPECT_EQ(again.status, NMPC::Status::Converged);
    EXPECT_EQ(again.accepted_steps, 0);
}

TYPED_TEST(SolverTest, RegularisationStartsAfreshEverySolve) {
    typedef TypeParam T;
    typedef Stop_Problem<T> Problem;
    typedef NMPC::Solver<Problem> Solver;
    static typename Solver::Workspace workspace;
    Solver solver(workspace);
    Problem problem;
    NMPC::Options<T> options;
    options.precision_floor = T(0);
    options.max_iterations = 3;
    const typename Problem::State x0(T(0.5), T(0));
    solver.reset(x0, Problem::Control::Zero());
    solver.solve(problem, x0, options);
    /* A quiet stretch: the same state again and again, where steps cannot be verified and raise the regularisation. */
    for (int i = 0; i < 10; i++) {
        solver.solve(problem, x0, options);
    }
    const typename Problem::State moved(T(1.5), T(-1));
    solver.advance(problem, 1);
    const NMPC::Result<T> result = solver.solve(problem, moved, options);
    EXPECT_GT(result.accepted_steps, 0);
    EXPECT_LE(result.regularisation, options.regularisation_initial);
}

TYPED_TEST(SolverTest, AdvanceFractionMovesEachCommandTowardTheNext) {
    typedef TypeParam T;
    typedef Stop_Problem<T> Problem;
    typedef NMPC::Solver<Problem> Solver;
    static typename Solver::Workspace workspace;
    Solver solver(workspace);
    Problem problem;
    const typename Problem::State x0(T(1), T(0));
    solver.reset(x0, Problem::Control::Zero());
    solver.solve(problem, x0, NMPC::Options<T>());
    const typename Solver::Policy before = solver.policy();
    solver.advance_fraction(T(0.25));
    const int N = Solver::Stages;
    for (int k = 0; k < N; k++) {
        const typename Problem::Control expected = k + 1 < N ?
            typename Problem::Control(before.u[k] + T(0.25) * (before.u[k + 1] - before.u[k])) : before.u[k];
        EXPECT_TRUE(solver.policy().u[k].isApprox(expected, T(1e-5)) || (expected.isZero() && solver.policy().u[k].isZero()));
        EXPECT_TRUE(solver.policy().K[k] == before.K[k]);
        EXPECT_TRUE(solver.policy().x[k] == before.x[k]);
    }
}

TYPED_TEST(SolverTest, AdvanceShiftsThePlanAndFillsTheTail) {
    typedef TypeParam T;
    {
        typedef Stop_Problem<T> Problem;
        typedef NMPC::Solver<Problem> Solver;
        static typename Solver::Workspace workspace;
        Solver solver(workspace);
        Problem problem;
        const typename Problem::State x0(T(1), T(0));
        solver.reset(x0, Problem::Control::Zero());
        solver.solve(problem, x0, NMPC::Options<T>());
        const typename Solver::Policy before = solver.policy();
        solver.advance(problem, 3);
        const int N = Solver::Stages;
        for (int k = 0; k + 3 < N; k++) {
            EXPECT_EQ(solver.policy().u[k], before.u[k + 3]);
            EXPECT_EQ(solver.policy().K[k], before.K[k + 3]);
            EXPECT_EQ(solver.policy().x[k], before.x[k + 3]);
        }
        /* Without a tail_control the entering stages repeat the last command and gain. */
        for (int k = N - 3; k < N; k++) {
            EXPECT_EQ(solver.policy().u[k], before.u[N - 1]);
            EXPECT_EQ(solver.policy().K[k], before.K[N - 1]);
            EXPECT_EQ(solver.policy().x[k], before.x[N]);
        }
    }
    {
        typedef Tail_Problem<T> Problem;
        typedef NMPC::Solver<Problem> Solver;
        static typename Solver::Workspace workspace;
        Solver solver(workspace);
        Problem problem;
        const typename Problem::State x0(T(1), T(0));
        solver.reset(x0, Problem::Control::Zero());
        solver.solve(problem, x0, NMPC::Options<T>());
        solver.advance(problem, 2);
        const int N = Solver::Stages;
        for (int k = N - 2; k < N; k++) {
            EXPECT_EQ(solver.policy().u[k](0), problem.tail);
            EXPECT_EQ(solver.policy().K[k].norm(), T(0));
        }
    }
}

TYPED_TEST(SolverTest, AdvanceByTheWholeHorizonRepeatsTheLastCommand) {
    /* Every stage enters at the end: each takes the old last command and gain, not the old first. */
    typedef TypeParam T;
    typedef Stop_Problem<T> Problem;
    typedef NMPC::Solver<Problem> Solver;
    static typename Solver::Workspace workspace;
    Solver solver(workspace);
    Problem problem;
    const typename Problem::State x0(T(1), T(0));
    solver.reset(x0, Problem::Control::Zero());
    solver.solve(problem, x0, NMPC::Options<T>());
    const int N = Solver::Stages;
    for (int k = 0; k < N; k++) {
        workspace.plan.u[k](0) = T(k + 1);
    }
    const typename Solver::Policy before = solver.policy();
    ASSERT_NE(before.K[0], before.K[N - 1]);
    for (const int stages : {N, N + 5}) {
        solver.advance(problem, stages);
        for (int k = 0; k < N; k++) {
            EXPECT_EQ(solver.policy().u[k](0), T(N)) << stages << ", stage " << k;
            EXPECT_EQ(solver.policy().K[k], before.K[N - 1]) << stages << ", stage " << k;
            EXPECT_EQ(solver.policy().x[k], before.x[N]) << stages << ", stage " << k;
        }
    }
}

TYPED_TEST(SolverTest, PinnedCommandGetsNoGainAndLeavesTheEarlierGainsAlone) {
    /*
    A command with equal bounds gets a zero gain, and the stages before it see a cost to go that it does
    not act on. At the zero state its gradient is zero, which used to leave it free: the gains were
    -0.6 and -0.5 in place of -2/3 and 0.
    */
    typedef TypeParam T;
    typedef Pinned_Problem<T> Problem;
    typedef NMPC::Solver<Problem> Solver;
    static typename Solver::Workspace workspace;
    Solver solver(workspace);
    Problem problem;
    NMPC::Options<T> options;
    options.regularisation_initial = T(0);
    options.regularisation_min = T(0);
    for (const T x : {T(0), T(1)}) {
        const typename Problem::State x0 = Problem::State::Constant(x);
        solver.reset(x0, Problem::Control::Zero());
        const NMPC::Result<T> result = solver.solve(problem, x0, options);
        ASSERT_TRUE(result.usable());
        EXPECT_NEAR(double(solver.policy().K[0](0, 0)), -2.0 / 3.0, 1e-6) << x;
        EXPECT_EQ(solver.policy().K[1](0, 0), T(0)) << x;
        EXPECT_EQ(solver.policy().u[1](0), problem.pinned) << x;
        EXPECT_NEAR(double(solver.policy().u[0](0)), -2.0 / 3.0 * double(x), 1e-6) << x;
    }
}

TYPED_TEST(SolverTest, AnExactOptimumOfZeroCostStopsAtOnce) {
    /*
    At rest at the target the cost and its predicted decrease are both zero. The solve stops after one
    backward pass as converged, whatever the precision floor, and does not spend a retry and a line search
    to report NoProgress.
    */
    typedef TypeParam T;
    typedef Stop_Problem<T> Problem;
    typedef NMPC::Solver<Problem> Solver;
    static typename Solver::Workspace workspace;
    Solver solver(workspace);
    Problem problem;
    for (const T floor : {T(0), NMPC::Options<T>().precision_floor}) {
        NMPC::Options<T> options;
        options.precision_floor = floor;
        solver.reset(Problem::State::Zero(), Problem::Control::Zero());
        const NMPC::Result<T> result = solver.solve(problem, Problem::State::Zero(), options);
        EXPECT_EQ(result.status, NMPC::Status::Converged);
        EXPECT_EQ(result.backward_passes, 1);
        EXPECT_EQ(result.line_search_trials, 0);
        EXPECT_EQ(result.cost_after, T(0));
        EXPECT_EQ(result.regularisation, options.regularisation_initial);
        /* The gains are those of the problem, for the next tick. */
        EXPECT_LT(solver.policy().K[0](0, 0), T(0));
    }
}

TYPED_TEST(SolverTest, PolicyFeedbackIsWhatTheForwardPassApplies) {
    typedef TypeParam T;
    typedef Attitude_Problem<T> Problem;
    typedef NMPC::Solver<Problem> Solver;
    typedef Eigen::Matrix<T, 3, 1> V3;
    static typename Solver::Workspace workspace;
    Solver solver(workspace);
    Problem problem;
    problem.torque_limit = T(10);   /* far from the bounds, so that no command is clamped */

    typename Problem::State x;
    x.q = NMPC::so3_exp<T>(V3(T(0.3), T(-0.2), T(0.4)));
    x.w = V3(T(0.1), T(0), T(-0.1));
    solver.reset(x, Problem::Control::Zero());
    solver.solve(problem, x, NMPC::Options<T>());
    const typename Solver::Policy plan = solver.policy();

    /* A rollout from a disturbed state with no step applies the gains alone. */
    typename Problem::State disturbed = x;
    disturbed.q = NMPC::attitude_retract<T>(x.q, V3(T(0.01), T(0.02), T(-0.01)));
    disturbed.w += V3(T(0.02), T(-0.01), T(0));
    int bound_stages = 0;
    solver.forward(problem, disturbed, T(0), bound_stages);
    const typename Problem::Control expected = plan.feedback(0, disturbed);
    for (int i = 0; i < 3; i++) {
        EXPECT_NEAR(double(solver.workspace().trial_u[0](i)), double(expected(i)), 1e-5);
    }
    EXPECT_GT(plan.K[0].norm(), T(0));
}

TYPED_TEST(SolverTest, PolicyLooksUpTheStageInForce) {
    typedef TypeParam T;
    typedef Stop_Problem<T> Problem;
    NMPC::Policy<Problem> plan;
    plan.origin = 100;
    plan.ticks_per_stage = 20;
    EXPECT_EQ(plan.stage_index(50), 0);
    EXPECT_EQ(plan.stage_index(100), 0);
    EXPECT_EQ(plan.stage_index(119), 0);
    EXPECT_EQ(plan.stage_index(120), 1);
    EXPECT_EQ(plan.stage_index(100 + 20 * 19), 19);
    EXPECT_EQ(plan.stage_index(100000), 19);
    EXPECT_EQ(plan.at(145).index, 2);
}

TYPED_TEST(SolverTest, PerStageJacobiansGiveTheSamePlans) {
    /*
    Jacobians made stage by stage in the backward pass, and not kept, give the plans of the stored ones
    with a smaller workspace, including when a backward pass is repeated and must make them again. The
    arithmetic is the same; only a compiler that fuses multiply-adds differently where it inlines
    differently could move the last bit, so the plans are compared to a few units of it.
    */
    typedef TypeParam T;
    typedef Attitude_Problem<T> Problem;
    typedef NMPC::Solver<Problem> Stored;
    typedef NMPC::Solver<Problem, Problem::Stages, NMPC::Jacobians::PerStage> PerStage;
    typedef Eigen::Matrix<T, 3, 1> V3;
    static_assert(sizeof(typename PerStage::Workspace) < sizeof(typename Stored::Workspace),
        "per-stage Jacobians keep a smaller workspace");
    /* Smaller by the Jacobians, less the few bytes of padding that the empty store may take. */
    EXPECT_GE(sizeof(typename Stored::Workspace) - sizeof(typename PerStage::Workspace) + 16,
        sizeof(typename Stored::MatA[Problem::Stages]) + sizeof(typename Stored::MatB[Problem::Stages]));
    static typename Stored::Workspace stored_workspace;
    static typename PerStage::Workspace per_stage_workspace;
    Stored stored(stored_workspace);
    PerStage per_stage(per_stage_workspace);
    Problem problem;
    const T tolerance = std::is_same<T, float>::value ? T(1e-5) : T(1e-13);

    typename Problem::State x;
    x.q = NMPC::so3_exp<T>(V3(T(0.6), T(-0.4), T(1.2)));
    x.w = V3(T(0.2), T(0), T(-0.3));
    stored.reset(x, Problem::Control::Zero());
    per_stage.reset(x, Problem::Control::Zero());
    for (int tick = 0; tick < 40; tick++) {
        NMPC::Options<T> options;
        /* Every fourth solve has no line search, so every backward pass is repeated with more regularisation. */
        if (tick % 4 == 3) {
            options.max_line_search = 0;
            options.max_regularisation_retries = 2;
        }
        if (tick > 0) {
            stored.advance(problem, 1);
            per_stage.advance(problem, 1);
        }
        const NMPC::Result<T> a = stored.solve(problem, x, options);
        const NMPC::Result<T> b = per_stage.solve(problem, x, options);
        ASSERT_EQ(a.status, b.status) << tick;
        ASSERT_EQ(a.backward_passes, b.backward_passes) << tick;
        EXPECT_EQ(a.box_iterations, b.box_iterations) << tick;
        for (int k = 0; k < Problem::Stages; k++) {
            EXPECT_LE((stored.policy().u[k] - per_stage.policy().u[k]).cwiseAbs().maxCoeff(), tolerance) << tick;
            EXPECT_LE((stored.policy().K[k] - per_stage.policy().K[k]).cwiseAbs().maxCoeff(),
                tolerance * (T(1) + stored.policy().K[k].cwiseAbs().maxCoeff())) << tick;
        }
        x = problem.step(0, x, stored.policy().u[0]);
    }
}
