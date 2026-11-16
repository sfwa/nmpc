/*
Host timings of the solver on the synthetic problem with the multirotor's
dimensions, in both precisions: a whole solve as the reference controller
runs it (shift, then one iteration), and its parts. Host timings are relative:
compare them between builds on one machine. The Cortex-M7 counts are in
cortex-m7/.
*/

#include <benchmark/benchmark.h>
#include "Synthetic.h"

namespace {

template <class T>
struct Fixture {
    using Problem = Synthetic_Problem<T>;
    using Solver = NMPC::Solver<Problem>;
    typename Solver::Workspace workspace;
    Solver solver;
    Problem problem;
    NMPC::Options<T> options;
    typename Problem::State x;

    Fixture() : solver(workspace) {
        x = problem.start();
        solver.reset(x, problem.hover());
        /* Settle into the steady pattern of a running controller. */
        for (int i = 0; i < 20; i++) {
            next();
        }
    }

    void next() {
        solver.advance(problem, 1);
        problem.phase += T(0.1);
        solver.solve(problem, x, options);
        x = problem.step(0, x, solver.policy().u[0]);
    }
};

template <class T>
Fixture<T>& fixture() {
    static Fixture<T>* f = new Fixture<T>();
    return *f;
}

}

template <class T>
void Solver_Solve(benchmark::State& state) {
    Fixture<T>& f = fixture<T>();
    while (state.KeepRunning()) {
        f.next();
    }
}

template <class T>
void Solver_Linearise(benchmark::State& state) {
    Fixture<T>& f = fixture<T>();
    while (state.KeepRunning()) {
        f.solver.linearise(f.problem);
    }
}

template <class T>
void Solver_Backward(benchmark::State& state) {
    Fixture<T>& f = fixture<T>();
    f.solver.linearise(f.problem);
    while (state.KeepRunning()) {
        T dV1 = T(0), dV2 = T(0);
        int passes = 0;
        benchmark::DoNotOptimize(f.solver.backward(f.problem, f.options, dV1, dV2, passes));
    }
}

template <class T>
void Solver_LineSearchTrial(benchmark::State& state) {
    Fixture<T>& f = fixture<T>();
    while (state.KeepRunning()) {
        int bound = 0;
        benchmark::DoNotOptimize(f.solver.forward(f.problem, f.x, T(1), bound));
    }
}

template <class T>
void Solver_Rollout(benchmark::State& state) {
    Fixture<T>& f = fixture<T>();
    while (state.KeepRunning()) {
        int bound = 0;
        benchmark::DoNotOptimize(f.solver.rollout(f.problem, f.x, bound));
    }
}

BENCHMARK_TEMPLATE(Solver_Solve, float);
BENCHMARK_TEMPLATE(Solver_Solve, double);
BENCHMARK_TEMPLATE(Solver_Linearise, float);
BENCHMARK_TEMPLATE(Solver_Linearise, double);
BENCHMARK_TEMPLATE(Solver_Backward, float);
BENCHMARK_TEMPLATE(Solver_Backward, double);
BENCHMARK_TEMPLATE(Solver_LineSearchTrial, float);
BENCHMARK_TEMPLATE(Solver_LineSearchTrial, double);
BENCHMARK_TEMPLATE(Solver_Rollout, float);
BENCHMARK_TEMPLATE(Solver_Rollout, double);
