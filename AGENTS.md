# AGENTS.md

Guidance for coding agents, and for people, working in this repository. The
[README](README.md) is for users of the library, [CONTRIBUTING.md](CONTRIBUTING.md)
says what a pull request needs, and [RELEASING.md](RELEASING.md) holds the
compatibility contract and the release process. Where this file and those
differ, they win.

## What this is

NMPC is a header-only C++17 library for nonlinear model predictive control on
small processors: an iLQR solver with box-constrained stage QPs, the problem
concept it solves, the policy and the inner loop that run between solves, and
derivative checks. It depends only on Eigen. It is version 2: version 1 is the
older C, C66x and Python tree, tagged `legacy`.

It has two siblings. [UKF](https://github.com/sfwa/ukf) shares NMPC's build
options, Eigen selection and CI matrix, and one program may use both.
[models](https://github.com/sfwa/models) builds the vehicles and controllers on
the two libraries.

## Layout

| Path | Contents |
| --- | --- |
| `include/NMPC/` | The library. `Version.h` is the only place the version is written. |
| `test/` | GoogleTest tests: `unittest` (typed, in `float` and `double`, see `Typed.h`) and `unittest_no_heap` (counts allocations). |
| `test/cmake/` | Consumer fixtures: as a subdirectory, beside UKF, installed (`installed/`), and a configuration that must fail (`mismatch/`). |
| `benchmark/` | Google Benchmark solves. `benchmark/cortex-m7/` is the Cortex-M7 instruction-count gate, with its baseline in `accepted.json`. |
| `examples/` | `stop`, `pendulum` and `attitude`. |
| `cmake/` | `ReadVersion.cmake` (identical in UKF), the package configuration, and `EigenGccWarnings.h`. |
| `patch/` | Patches applied to the bundled Eigen 3.2.8. |

## Build and test

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --target unittest unittest_no_heap benchmark stop pendulum attitude --parallel 2
ctest --test-dir build --output-on-failure
```

The default build downloads Eigen 3.2.8 and patches it. An installed
GoogleTest and Google Benchmark are used when found, and are downloaded
otherwise. The [README](README.md#building) lists the options, including
another Eigen (`NMPC_USE_SYSTEM_EIGEN`, `NMPC_EIGEN_VERSION`) and the
unrolling limit. Project targets build with `-Wall -Wextra -Wshadow
-Wconversion -Wdouble-promotion -Werror`, without exceptions or RTTI.

If you change CMake files, also run the consumer fixtures, as CI does:

```sh
cmake -S test/cmake -B build-consumer && cmake --build build-consumer && ctest --test-dir build-consumer
```

Add `-DNMPC_TEST_UKF_DIR=path/to/ukf` to build them beside UKF. The CI
workflow, `.github/workflows/cmake.yml`, has the installed and mismatch
checks.

If you change solver code, run the Cortex-M7 gate. It needs
`gcc-arm-none-eabi`, `libnewlib-arm-none-eabi` and `qemu-system-arm`, and the
bundled Eigen (`cmake --build build --target eigen3`):

```sh
python3 benchmark/cortex-m7/check.py build/src/eigen3
```

It fails if NMPC links a double-precision routine in single precision, or if
a count moves more than 3 % from the baseline. Do not update `accepted.json`
to make it pass. Explain the change in the pull request, and the maintainers
accept a new baseline when it is worth it.

## Formatting and lint

Both are standard tools with stock settings, checked in CI by the `lint` job.

- **clang-format 18**, Google style, unmodified (`.clang-format`). Format
  everything before committing:

  ```sh
  git ls-files '*.h' '*.cpp' '*.c' | xargs clang-format-18 -i
  ```

  Other clang-format versions format a little differently, so use 18. A
  commit that only reformats goes in `.git-blame-ignore-revs`, so that
  `git blame` skips it. Where
  include order matters (a vendor header that must come first), wrap the
  block in `// clang-format off` and `// clang-format on`, with a comment
  saying why.

- **clang-tidy 18**, with the `bugprone-*`, `performance-*` and
  `clang-analyzer-*` groups as errors (`.clang-tidy` lists the four checks
  that are off, and why). It runs over the tests, benchmark and examples, and
  through them over the headers:

  ```sh
  cmake -S . -B build-tidy -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_COMPILER=clang++-18 -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
  cmake --build build-tidy --target eigen3
  run-clang-tidy-18 -p build-tidy -quiet "$PWD/(test|benchmark|examples)/"
  ```

  Fix a finding rather than silencing it. If a `NOLINT` is unavoidable, name
  the check and give the reason on the same line.

## Rules of the code

- **No exceptions, RTTI or heap in a solve.** Failures are reported through
  return values. `unittest_no_heap` counts allocations, and project targets
  build with `-fno-exceptions -fno-rtti`.
- **Both precisions.** Library code is templated on the scalar type and must
  work in `float` and `double`. Tests use the typed fixtures in `Typed.h` to
  run in both. In single precision NMPC must call no double-precision routine:
  write literals as `T(0.5)` and use `std::` maths functions. The Cortex-M7
  gate checks this.
- **Eigen 3.2.8 to 5.0.** The bundled Eigen 3.2.8 is the oldest supported.
  It has no `Eigen::Index` (use `std::ptrdiff_t`), and under `-Wconversion`
  it rejects `Vec v = A * b;` (write `Vec v(A * b);`).
- **CMake 3.10** is the minimum. Keep new CMake within it.
- **Keep step with UKF.** The build options, the Eigen selection, the CI
  matrix, the lint setup, `cmake/ReadVersion.cmake` and `EigenGccWarnings.h`
  match UKF's. A change to one of them belongs in both repositories.
- **Comments** explain why. Write them as full sentences.

## Versions, compatibility and the changelog

- Read [RELEASING.md](RELEASING.md) before changing anything under
  `include/NMPC/`. Covered interfaces cannot change incompatibly outside a
  major release. `NMPC::detail` is not covered. `InnerLoop.h` and
  `Reference.h` are experimental and may change in a minor release.
- Say in the pull request if a change can move a plan or the work a solve
  does, even when it is a fix.
- Fill in the pull request template's changelog line and compatibility
  guess. The maintainers write `CHANGELOG.md` and set versions. Do not edit
  `Version.h` or a released section of the changelog unless you are making a
  release.
- Never create, move or delete tags or GitHub Releases. Release tags are
  immutable and made only by the maintainers, following RELEASING.md.

## Before you finish

1. Formatted with clang-format 18, and clang-tidy 18 reports nothing.
2. Tests, benchmark and examples build with warnings as errors, and `ctest`
   passes.
3. Consumer fixtures, if CMake files changed. Cortex-M7 gate, if solver code
   changed.
4. Documentation for anything a user can see: the README or the header
   comment.
