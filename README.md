# NMPC

Nonlinear model predictive control for embedded systems. NMPC is to control
what [UKF](https://github.com/sfwa/ukf) is to estimation: a header-only C++17
library on [Eigen](https://eigen.tuxfamily.org), with fixed memory, no heap,
no exceptions and no virtual calls, that knows nothing about vehicles.

It provides:

  * A box-constrained Gauss-Newton iLQR solver with analytic derivatives,
    built to run one iteration per control period from the previous plan (a
    real-time iteration) on MCU-class hardware.
  * The plan as a policy: the planned state, command and feedback gain of
    every stage, for an inner loop that runs faster than the solver.
  * Building blocks for that inner loop: an incremental (INDI) correction
    with its low-pass filters, clamping, a reference table, and a reference
    clock that keeps time in single precision.
  * Derivative checks against central differences, for tests.

A problem is any type that meets the concept in
[`Problem.h`](include/NMPC/Problem.h): a model step, residuals whose squares
represent the cost, their analytic Jacobians, and bounds on the commands. States 
may live on a manifold, such as one with a quaternion.

Every type is a template on its scalar `T`, so one program can hold solvers in
`float` and in `double`. NMPC never defines `real_t`, which belongs to UKF.

## Releases

[![Latest release](https://img.shields.io/github/v/release/sfwa/nmpc?include_prereleases&sort=semver)](https://github.com/sfwa/nmpc/releases)

NMPC uses [Semantic Versioning](https://semver.org), with the same release
process as UKF. Its version numbers start at 2, since the earlier
implementation (tag `legacy`) counts as version 1. The first release will be
`2.0.0`, after a release candidate `2.0.0-rc.1`; until then master is the
only version. Use the [latest release](https://github.com/sfwa/nmpc/releases/latest)
and [pin it by its tag](RELEASING.md#using-a-release).

- [RELEASING.md](RELEASING.md): the compatibility contract (which interfaces,
  conventions, build options and configurations a release supports, and
  which are experimental), what patch, minor and major releases may change,
  how numerical changes are reported, and how releases are made.
- [CHANGELOG.md](CHANGELOG.md): what has changed in each release.
- [CONTRIBUTING.md](CONTRIBUTING.md): what a pull request needs.

## Usage

### A problem

A point mass on a line that must come to rest at the origin within an
acceleration limit:

```C++
#include "NMPC/Solver.h"

template <class T>
struct Stop_Problem {
    using Scalar = T;
    enum { NX = 2, NU = 1, NR = 3, NT = 2, Stages = 20 };
    using State   = Eigen::Matrix<T, 2, 1>;     /* position (m), velocity (m/s) */
    using Control = Eigen::Matrix<T, 1, 1>;     /* acceleration (m/s^2) */

    T dt = T(0.05);                             /* stage length, s */
    T limit = T(4);                             /* acceleration limit, m/s^2 */

    static State retract(const State& x, const State& dx) { return x + dx; }
    static State difference(const State& x, const State& y) { return x - y; }

    State step(int, const State& x, const Control& u) const {
        return State(x(0) + dt * x(1) + T(0.5) * dt * dt * u(0), x(1) + dt * u(0));
    }
    template <class A, class B>
    void step_jacobian(int, const State&, const Control&, A& a, B& b) const {
        a << T(1), dt, T(0), T(1);
        b << T(0.5) * dt * dt, dt;
    }

    /* Residuals carry the root of the stage length, so the weights are per second. */
    template <class R>
    void stage_residual(int, const State& x, const Control& u, R& r) const {
        const T s = std::sqrt(dt);
        r << s * x(0), s * T(0.5) * x(1), s * T(0.1) * u(0);
    }
    template <class JX, class JU>
    void residual_jacobians(int, const State&, const Control&, JX& jx, JU& ju) const {
        const T s = std::sqrt(dt);
        jx << s, T(0), T(0), s * T(0.5), T(0), T(0);
        ju << T(0), T(0), s * T(0.1);
    }
    template <class R> void terminal_residual(const State& x, R& r) const { r = T(3) * x; }
    template <class JX> void terminal_jacobian(const State&, JX& jx) const { jx = T(3) * JX::Identity(); }

    void control_bounds(int, Control& lower, Control& upper) const { lower(0) = -limit; upper(0) = limit; }
};
```

The cost is half the sum of the squared residuals, so the solver forms the
Gauss-Newton Hessian with no second derivatives. A problem may also give:

  * `state_scale` and `control_scale`, typical sizes that condition the solver
    in single precision;
  * `linear_cost` and `linear_cost_gradient`, a term with no Hessian such as a
    reward for progress along a path;
  * `cost_expansion` and `terminal_expansion`, its own Gauss-Newton expansion
    for speed where it knows the structure of its residuals;
  * `tail_control`, the command for a stage that enters the plan when it is
    shifted;
  * `finite`, required when `State` is not an Eigen type.

[`Problem.h`](include/NMPC/Problem.h) documents each of them.

NMPC says *retract* and *difference* where UKF says `apply_delta` and *delta*:
`retract(x, d)` moves a state by a tangent, and `difference(x, y)` is the
tangent at `y` that reaches `x`. [`Types.h`](include/NMPC/Types.h) has the
attitude versions, with a perturbation applied on the left.

### Solving

```C++
using Problem = Stop_Problem<float>;
using Solver  = NMPC::Solver<Problem>;

static Solver::Workspace workspace;             /* owned by the caller; nothing is allocated */
Solver solver(workspace);
Problem problem;
NMPC::Options<float> options;                   /* embedded defaults: 1 iteration, 5 candidates, 1 retry */

solver.reset(x0, Problem::Control::Zero());     /* the first plan is a constant command */

/* Once per solve period. */
solver.advance(problem, 1);                     /* shift the plan one stage; advance_fraction(f) for a period shorter than a stage */
const NMPC::Result<float> result = solver.solve(problem, x_predicted, options);
if (result.usable()) {
    const auto& plan = solver.policy();         /* plan.x[k], plan.u[k] and plan.K[k] for every stage */
    publish(plan, tick_due + ticks_per_solve);  /* it comes into force one period after the solve began */
} else {
    solver.reset(x_predicted, u_safe);          /* the plan in force flies on, and the next solve starts afresh */
}
```

Two rules belong to the application. The plan comes into force one solve
period after the solve starts, so solve from the state predicted for that
moment; at periods longer than about 20 ms that prediction must follow the
stored gains as well as the planned commands. And a solve that is not usable
has left the plan exactly as it was, so the plan in force can carry on into
its next stage.

`Options` sets the work a solve may do, in units of work rather than seconds:
iterations, line-search candidates, regularisation retries and stage QP
passes, and the tolerances it works to. The stage QP's stop is relative to
the size of the terms in each gradient (`box_tolerance`), so a problem's
answer does not depend on how its cost is scaled. A command whose bounds meet
is held and gets no gain. `Result` says whether the plan is usable, why the solve ended, the
work done, the costs before and after, the predicted decrease, and how many
stages have a command on a bound.

### The policy

`solver.policy()` is the plan as a plain value: for each stage its planned
state, command and gain. An inner loop between solves applies

```C++
command = plan.u[k] + plan.K[k] * Problem::difference(x, x_planned);
```

which `plan.feedback(k, x, x_planned)` computes. Gains are specified in terms
of the problem's units, and a command held at a bound has no gain. An application 
that solves below its inner loop in priority copies the policy whole into a buffer
that the inner loop switches to on a tick boundary; `origin`,
`ticks_per_stage` and `at(tick)` find the stage in force.

### The inner loop

The inner-loop and reference headers are experimental in 2.0: they may change
in a minor release ([RELEASING.md](RELEASING.md#experimental-interfaces)).

The inner loop is the sum of three terms: the plan's command at this tick,
an error law on the distance from the plan (the stored gains above), and a
correction from the measured acceleration. [`InnerLoop.h`](include/NMPC/InnerLoop.h)
holds the pieces that do not depend on the vehicle:

  * `Incremental<T, NY, NA>`: NY residuals, each the measured value less the
    model's expectation, low-passed by a second-order Butterworth each, scaled
    by a gain, and allocated to NA actuator increments within per-actuator
    limits. The expectation is the application's, from the model it plans
    with. `NoCorrection` holds nothing.
  * `AttitudeCascade<T>`: an attitude and rate cascade, as the error law of a
    rigid body in place of the stored gains: the attitude error to a rate
    command, and a PID on the rate to a moment about each body axis, whose
    integral stops while the command is clipped. The planned rate is turned
    from the plan's body axes into the vehicle's before it is compared with the
    measured rate. The defaults are the review prototype's for the F450.
  * `Butterworth2`, `clamp` and `on_bound`.

What the correction holds is taken out of the solver's bounds and of the
actuator state the solver sees, so that the plan is made within what is left.

[`Reference.h`](include/NMPC/Reference.h) holds a reference as a table of rows
at a fixed spacing (`PathTable`, which may hold a window of a longer table;
a lookup outside the window holds the window's nearest end and is counted), and the reference's `Clock`, a sum compensated for
rounding: a plain single-precision sum at 1 kHz drifts 2 to 5 ms in 20 to 30 s.
What the columns mean is the vehicle's, so the vehicle's path type turns the
two rows of a lookup into a state and a command.

### Derivative checks

[`Check.h`](include/NMPC/Check.h) compares a problem's step Jacobians,
residual Jacobians, terminal Jacobian, linear cost gradient and own cost
expansion with central differences in the tangent. Use it in tests: a solve
never differentiates numerically.

```C++
const NMPC::CheckResult<double> check = NMPC::check_step_jacobian(problem, k, x, u, 1e-7);
EXPECT_TRUE(check.passed) << "worst " << check.error << " at " << check.row << ", " << check.column;
```

## Building

NMPC is header-only. Add it as a subdirectory, after UKF if the program uses
both, and link `nmpc::nmpc`:

```cmake
add_subdirectory(ukf)
add_subdirectory(nmpc)
target_link_libraries(controller PRIVATE nmpc::nmpc ukf::single)
```

To use an installed copy instead, build and install NMPC, then find its
package:

```sh
cmake -S path/to/nmpc -B build-nmpc -DCMAKE_INSTALL_PREFIX=/opt/nmpc
cmake --build build-nmpc
cmake --install build-nmpc
```

```cmake
find_package(nmpc 2.0 REQUIRED CONFIG) # With CMAKE_PREFIX_PATH including /opt/nmpc.
target_link_libraries(controller PRIVATE nmpc::nmpc)
```

The package gives the same `nmpc::nmpc` target, carrying the unrolling limit
the installed copy was configured with, and sets `nmpc_VERSION` and
`nmpc_VERSION_PRERELEASE`. It accepts a request for any earlier version with
the same major version. A pre-release, such as master's `2.0.0-dev` or
`2.0.0-rc.1`, comes before `2.0.0`, so it does not meet a request for `2.0`:
find it with `find_package(nmpc REQUIRED CONFIG)`. With the downloaded Eigen, its headers (patched, for
3.2.8) are installed under `include/nmpc/eigen3`; with
`NMPC_USE_SYSTEM_EIGEN`, the package finds Eigen again when it is loaded.
Install rules are generated when NMPC is built on its own, or with
`NMPC_INSTALL=ON` as a subdirectory. A program that finds installed copies of
both NMPC and UKF should have installed both with the same Eigen.

The options are UKF's, with `NMPC_` for `UKF_`:

| CMake option | Default | Purpose |
| --- | --- | --- |
| `NMPC_USE_SYSTEM_EIGEN` | UKF's choice beside UKF, else `OFF` | Use an installed Eigen 3.3 or newer (`Eigen3::Eigen`) instead of downloading one |
| `NMPC_EIGEN_VERSION` | UKF's beside UKF, else `3.2.8` | Eigen release to download when not using an installed one: 3.2.8 (patched) or 3.3 and newer (unpatched) |
| `NMPC_EIGEN_URL` | UKF's beside UKF, else empty (GitLab release) | Eigen source archive, as a URL or local file, to download instead; set `NMPC_EIGEN_VERSION` to its version |
| `NMPC_EIGEN_UNROLLING_LIMIT` | UKF's beside UKF, else empty (Eigen's 100) | Set `EIGEN_UNROLLING_LIMIT` for every target that links `nmpc` |
| `NMPC_WARNINGS_AS_ERRORS` | `ON` | Treat warnings in NMPC's own targets as errors |
| `NMPC_BUILD_TESTS` | `ON` standalone, `OFF` as a subdirectory | The unit tests (`unittest`, `unittest_no_heap`, and `check` to run them) |
| `NMPC_BUILD_BENCHMARKS` | `ON` standalone, `OFF` as a subdirectory | The host benchmark (`benchmark`) |
| `NMPC_BUILD_EXAMPLES` | `ON` standalone, `OFF` as a subdirectory | The examples (`stop`, `pendulum`, `attitude`) |
| `NMPC_INSTALL` | `ON` standalone, `OFF` as a subdirectory | Generate install rules and the `nmpc` CMake package |

NMPC and UKF must see one Eigen. By default both download Eigen 3.2.8 and
apply the patches in [`patch/`](patch), which keep its LLT free of heap
allocation and let it compile as C++17; the Cortex-M7 figures below use it,
and on an STM32H735 it ran these kernels 11 to 15 % faster than 3.4.0. When
UKF is configured first, NMPC takes UKF's download rather than fetch its own,
its Eigen options default to UKF's, and a configure that sets them
differently fails rather than build a program with two Eigens. To use another
Eigen, set the same choice for both libraries:

- download another release with `-DUKF_EIGEN_VERSION=3.4.0` (or `5.0.0`, or
  any release from 3.3), used unpatched, and `-DUKF_EIGEN_URL=` to build
  offline or from a mirror; or
- use an installed Eigen 3.3 or newer, including Eigen 5, with
  `-DUKF_USE_SYSTEM_EIGEN=ON`, choosing between installations with
  `-DEigen3_DIR=` or `CMAKE_PREFIX_PATH`.

Use the `NMPC_` names when NMPC stands alone. A solve stays free of heap
allocation with every supported Eigen.

`NMPC_WARNINGS_AS_ERRORS` (on) applies to NMPC's own targets, which build
with `-Wall -Wextra -Wshadow -Wconversion -Wdouble-promotion`, without
exceptions or RTTI. None of it reaches a consumer.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target unittest unittest_no_heap benchmark stop pendulum attitude
ctest --test-dir build
```

The tests use an installed GoogleTest where there is one and fetch
GoogleTest 1.7.0 otherwise, as UKF's do; the benchmark does the same with
Google Benchmark 1.3.0. Every test runs in `float` and in `double`. CI runs
them with GCC 13 and 15, Clang 18 and 22 on Linux and Clang on macOS, in
Debug and Release, with Eigen 3.2.8, 3.4.0 and 5.0.0: the same
configurations as UKF's. [`test/cmake`](test/cmake) holds the consumer
fixtures: NMPC's tools stay opt-in as a subdirectory, it builds without
exceptions, a sibling keeps the parent's language and flags, beside UKF it
shares UKF's Eigen (`-DNMPC_TEST_UKF_DIR=<ukf checkout>`), and the installed
package works alone and beside an installed UKF
([`test/cmake/installed`](test/cmake/installed)).

## Embedded use

A solve allocates nothing (a test counts allocations). Its memory is in
three parts, each of which an application provisions:

  * the workspace, `sizeof(Solver::Workspace)`, which the caller owns: the
    plan, the trial plan, the stage QPs' steps, the candidate gains of the
    backward pass (which become the plan's only when a whole pass succeeds)
    and, by default, every stage's Jacobians;
  * the solver object, `sizeof(Solver)`, which holds the scales of the last
    solve and their products and lives as long as the workspace;
  * the call stack of a solve.

In single precision NMPC's own code calls no double-precision routine.

By default the workspace keeps every stage's step Jacobians for the
iteration. `NMPC::Solver<Problem, Problem::Stages, NMPC::Jacobians::PerStage>`
makes each stage's Jacobians as the backward pass reaches it instead, and
keeps none: the plans are the same, the workspace is smaller by
N (NX² + NX NU) scalars, the stack is larger by one stage's Jacobians, and a
backward pass repeated with more regularisation makes them again.

[`benchmark/cortex-m7`](benchmark/cortex-m7) counts instructions on a
Cortex-M7 under QEMU for a synthetic problem with the dimensions of the
multirotor reference controller (17 states, 5 commands, 30 stages). Counts
are exact for one build; on an STM32H735 at 550 MHz whole solves have run at
1.19 to 1.22 cycles per counted instruction. `check.py` is the CI gate: it
fails if NMPC's code links a double-precision routine or a count moves more
than 3 % from the baseline accepted for the same compiler.

```sh
python3 benchmark/cortex-m7/check.py build/src/eigen3
```

| Synthetic problem, single precision, `-O2`, arm-none-eabi-gcc 13.2.1 | Instructions | With per-stage Jacobians |
| --- | --- | --- |
| Solve (shift and one iteration) | 2.69 million | 2.65 million |
| First solve, from a reset plan | 2.79 million | 2.75 million |
| Backward pass | 2.24 million | 2.32 million (it makes the Jacobians) |
| Linearisation | 0.13 million | none |
| Line-search trial | 0.12 million | 0.12 million |
| Rollout | 0.11 million | 0.11 million |

| Memory for that problem, single precision | Default | With per-stage Jacobians |
| --- | --- | --- |
| Workspace | 72.5 kB | 27.6 kB |
| Solver object | 3.3 kB | 3.3 kB |
| Call stack of a solve, from the stack's top | 8.7 kB | 10.2 kB |

The benchmark keeps the workspace and the solver object in static storage,
and measures the stack from its original top, so the three add up to what
the solver needs. Of the workspace, 10.2 kB is the backward pass's candidate
gains.

## Licence

MIT. See [LICENSE](LICENSE).
