/*
Instruction counts for the solver on a Cortex-M7, in single precision, on the
synthetic problem with the multirotor's dimensions. It runs bare-metal under
QEMU's mps2-an500 machine with -icount shift=0, where virtual time advances
one nanosecond per instruction, so the counts are exact for one build.

Counts are instructions, not time. On the STM32H735 at 550 MHz whole NMPC
solves took 1.19 to 1.22 cycles per counted instruction. Counts move by up to
15 % between compiler versions, so the gate in check.py compares a build only
with a baseline made by the same compiler.

The same file builds on the host, where it checks only that every part runs.
*/

#include "Synthetic.h"
#include <cstdint>

#ifdef __arm__
extern "C" {
void print(const char*);
void print_u(uint64_t);
void timer_start(void);
uint64_t instructions(void);
void stack_paint(void);
uint32_t stack_used(void);
}
#else
#include <cstdio>
static void print(const char* s) { std::fputs(s, stdout); }
static void print_u(uint64_t v) { std::printf("%llu", static_cast<unsigned long long>(v)); }
static void timer_start() {}
static uint64_t instructions() { return 0; }
static void stack_paint() {}
static uint32_t stack_used() { return 0; }
#endif

using Problem = Synthetic_Problem<float>;
using Solver = NMPC::Solver<Problem>;

static Solver::Workspace workspace;
static Problem problem;

static void report(const char* name, uint64_t count, uint32_t stack) {
    print(name);
    print(": ");
    print_u(count);
    print(" instructions, stack ");
    print_u(stack);
    print(" bytes\n");
}

template <class F>
static void measure(const char* name, uint32_t repeats, F&& f) {
    stack_paint();
    const uint64_t start = instructions();
    for (uint32_t i = 0; i < repeats; i++) {
        f();
    }
    const uint64_t end = instructions();
    report(name, (end - start) / repeats, stack_used());
}

int main() {
    timer_start();
    print("sizeof: workspace ");
    print_u(sizeof(Solver::Workspace));
    print(", problem ");
    print_u(sizeof(Problem));
    print(" bytes\n");

    Solver solver(workspace);
    NMPC::Options<float> options;
    Problem::State x = problem.start();
    solver.reset(x, problem.hover());

    /* The first solve, from a reset plan, is the longest repeatable solve of a flight. */
    measure("first solve", 1, [&]() { solver.solve(problem, x, options); });
    x = problem.step(0, x, solver.policy().u[0]);

    /* The steady pattern of a running controller: shift, solve once, step the plant. */
    for (int i = 0; i < 10; i++) {
        solver.advance(problem, 1);
        problem.phase += 0.1f;
        solver.solve(problem, x, options);
        x = problem.step(0, x, solver.policy().u[0]);
    }
    int usable = 0;
    measure("solve", 10, [&]() {
        solver.advance(problem, 1);
        problem.phase += 0.1f;
        usable += solver.solve(problem, x, options).usable() ? 1 : 0;
        x = problem.step(0, x, solver.policy().u[0]);
    });

    measure("linearise", 5, [&]() { solver.linearise(problem); });
    measure("backward pass", 5, [&]() {
        float dV1 = 0.0f, dV2 = 0.0f;
        int passes = 0;
        solver.backward(problem, options, dV1, dV2, passes);
    });
    measure("line search trial", 5, [&]() {
        int bound = 0;
        solver.forward(problem, x, 1.0f, bound);
    });
    measure("rollout", 5, [&]() {
        int bound = 0;
        solver.rollout(problem, x, bound);
    });

    print(usable == 10 ? "all solves usable\n" : "UNUSABLE SOLVE\n");
    return usable == 10 ? 0 : 1;
}
