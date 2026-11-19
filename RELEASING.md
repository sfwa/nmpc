# Releases and compatibility

NMPC follows [Semantic Versioning](https://semver.org), as
[UKF](https://github.com/sfwa/ukf/blob/master/RELEASING.md) does, and the two
libraries share one release process. This file says what a version number
promises, what is and is not covered by that promise, and how a release is
made. [CONTRIBUTING.md](CONTRIBUTING.md) says what a contributor needs to do,
which is much less. User-facing changes are listed in
[CHANGELOG.md](CHANGELOG.md).

Version numbers start at 2. The earlier C, C66x and Python implementation,
kept at the tag `legacy`, counts as version 1; it had no numbered releases and
is not maintained.

## The compatibility contract

From `2.0.0`, the maintainers support the interfaces and behaviour listed
here. NMPC is header-only, so the promise is **source compatibility**: code
that compiles and behaves as documented against one release of a major
version compiles and behaves as documented against every later release of
that major version, in a supported configuration. There is no promise of
binary compatibility between releases, and no promise of bitwise-identical
plans (see [Numerical changes](#numerical-changes)).

### Covered

- **The problem concept** in [`Problem.h`](include/NMPC/Problem.h): the
  members a problem must provide (`Scalar`, `State`, `Control`, `NX`, `NU`,
  `NR`, `NT`, `Stages`, `retract`, `difference`, `step`, `step_jacobian`,
  `stage_residual`, `residual_jacobians`, `terminal_residual`,
  `terminal_jacobian`, `control_bounds`), the optional ones (`state_scale`,
  `control_scale`, `linear_cost`, `linear_cost_gradient`, `cost_expansion`,
  `terminal_expansion`, `tail_control`, `finite`), and what each means. A
  problem written against one release of a major version keeps working; a new
  optional member may be added in a minor release, a new required one only in
  a major release.
- **The solver.** `NMPC::Solver<Problem, N, Jacobians>`, its template
  parameters and defaults, its member types, and the members `Solver(workspace)`,
  `reset`, `advance`, `advance_fraction`, `solve`, `policy` and `workspace`.
  `NMPC::Workspace` and `NMPC::Jacobians` as types the caller declares and
  owns.
- **Options and results.** Every field of `NMPC::Options` and its default
  value; `NMPC::Result`, its fields and `usable()`; the values of
  `NMPC::Status` and what each says about the plan.
- **The policy.** `NMPC::Policy` and its stages (`x`, `u`, `K`), `feedback`,
  `origin`, `ticks_per_stage` and `at`.
- **Derivative checks** in [`Check.h`](include/NMPC/Check.h) and
  `NMPC::CheckResult`.
- **Attitude and matrix helpers** in [`Types.h`](include/NMPC/Types.h):
  `Vector`, `Matrix`, `Quaternion`, `hat`, `so3_exp`, `so3_log`,
  `so3_left_jacobian_inverse`, `attitude_retract` and `attitude_difference`.
- **Conventions.** The cost is half the sum of squared residuals plus the
  linear term; residuals carry the root of the stage length; `retract` and
  `difference` work in the tangent with a left perturbation for attitude; a
  plan comes into force one solve period after the solve begins; gains are in
  the problem's units and a command held at a bound has no gain; a solve that
  is not usable leaves the plan exactly as it was; `Options` counts work, not
  time.
- **Error handling and resources.** No exceptions, RTTI or virtual calls are
  needed; failures are reported in `Result`. A solve allocates nothing: its
  memory is the caller's workspace, the solver object and the call stack. In
  single precision NMPC's own code calls no double-precision routine.
  `NMPC_ASSERT` is the hook for checks on the caller's use of the library.
- **Build interface.** The CMake target `nmpc::nmpc` (and `nmpc`); the
  installed package `find_package(nmpc)`; the options `NMPC_USE_SYSTEM_EIGEN`,
  `NMPC_EIGEN_VERSION`, `NMPC_EIGEN_URL`, `NMPC_EIGEN_UNROLLING_LIMIT`,
  `NMPC_WARNINGS_AS_ERRORS`, `NMPC_BUILD_TESTS`, `NMPC_BUILD_BENCHMARKS`,
  `NMPC_BUILD_EXAMPLES` and `NMPC_INSTALL` and their defaults, including
  taking UKF's Eigen when UKF is configured first; and the version macros in
  [`NMPC/Version.h`](include/NMPC/Version.h).
- **Supported configurations**, below.

### Not covered

- Anything in `namespace NMPC::detail`; `NMPC::ProblemTraits`;
  `NMPC::solve_box_qp`, which is the solver's stage QP and is public for its
  tests; and the solver's `rollout`, `linearise`, `linearise_stage`,
  `backward` and `forward`, which are public so that a benchmark can time each
  part of a solve.
- The members and layout of `NMPC::Workspace`, and private members of any
  class. `sizeof(Solver::Workspace)` and `sizeof(Solver)` are reported, not
  promised, and may grow in a minor release.
- Header layout beyond the documented `#include "NMPC/..."` lines.
- Plans to the last bit, instruction counts, timing, code size and stack use.
  They are measured and reported (see [Numerical changes](#numerical-changes)
  and [Release gates](#release-gates)) but not promised.
- The examples, the benchmarks and their accepted baselines.
- Configurations outside the supported list. They may work, and fixes are
  welcome, but a release may break them.

### Experimental interfaces

An interface marked experimental is in a release so that users can try it, but
it may change in a minor release. Experimental interfaces are listed here and
say so in their header; nothing else is experimental. An interface leaves this
list in a minor release, after which the contract covers it.

| Interface | Since | Why it is experimental |
| --- | --- | --- |
| [`InnerLoop.h`](include/NMPC/InnerLoop.h): `Incremental`, `NoCorrection`, `AttitudeCascade` (including its defaults, which are one prototype airframe's), `Butterworth2`, `clamp`, `on_bound` | 2.0.0 | The inner-loop pieces are still being shaped by the multirotor and fixed-wing controllers in sfwa/models. |
| [`Reference.h`](include/NMPC/Reference.h): `PathTable`, `Clock` | 2.0.0 | As above, for the reference and its clock. |

### Supported configurations

What CI tests is what is supported. A newer compiler or Eigen is not supported
until CI tests it, and a release notes when one is added. These are the same
as UKF's, so that one program can use both libraries in any of them.

| | Supported |
| --- | --- |
| Language | C++17 |
| Precision | `float` and `double`, in one program if wanted |
| Compilers | GCC 13 and 15, Clang 18 and 22 on Linux; Apple Clang on the current macOS runner |
| Embedded | Cortex-M7 in single precision with Ubuntu 24.04's `arm-none-eabi-gcc` and Eigen 3.2.8, built and instruction-counted under QEMU in CI |
| Eigen | 3.2.8 with this repository's patches (the default download), 3.4.0 and 5.0.0, downloaded or installed |
| CMake | 3.10 or newer to configure; 3.14 or newer to install the package version file as architecture-independent |
| Generators | Unix Makefiles and Ninja |
| Build types | Debug and Release |

MSVC is not in CI. The maintainers also measure solves on an STM32H735 and
report what they find, but board timings are not part of the contract.
Raising a minimum in this table is a major change; adding a newer compiler or
Eigen is a minor change.

## What each part of the version means

| Release | Policy | Example |
| --- | --- | --- |
| **Patch: `2.0.1`** | Compatible corrections and implementation improvements | Fix the stage QP's stop while keeping the documented result; make the backward pass faster |
| **Minor: `2.1.0`** | Compatible features and deprecations | Add an optional problem member, an option with a default that keeps today's behaviour, or a new inner-loop piece; support a newer compiler or Eigen |
| **Major: `3.0.0`** | Changes that need downstream adaptation | Add a required problem member; change an `Options` default, a `Status` meaning or the cost convention; remove a deprecated interface; raise a supported minimum |

The maintainers classify every change and choose each version.

### Numerical changes

A program that still compiles does not show that its controller behaves as it
did. Every change that can move a plan (its states, commands or gains, the
status of a solve, or the work it does) gets its own entry under **Numerical
changes** in the changelog and the release notes, whatever kind of release it
is in. The entry says:

- which configurations are affected (precision, Eigen, problem features such
  as bounds that meet, a linear cost or a problem's own cost expansion);
- what the effect is expected to be, and its size;
- the regression evidence: tests, the Cortex-M7 instruction counts, and flights
  of the sfwa/models controllers, with links;
- whether downstream users should re-fly recorded or simulated cases, or
  retune.

A numerical bug fix can be a patch even when plans change: incorrect
mathematics is not kept to preserve old output. A deliberate change to a
documented convention or an existing default is a major change.

### Deprecation

1. A replacement arrives in a minor release, with migration notes in the
   changelog.
2. The old interface is marked `[[deprecated("use ...")]]` in the same or a
   later minor release, and keeps working.
3. It is removed only in the next major release, and no sooner than six
   months after it was first marked deprecated.

### Support

The latest minor release of the current major version gets fixes. A previous
major version gets fixes only for a maintenance window announced in the
release notes of the major release that replaces it; without such an
announcement there is none.

## Cadence and branches

- A patch release follows when useful fixes are ready. A minor release is due
  every six to eight weeks while development is active, and is skipped when
  there is nothing to release. A major release follows a compatibility change
  that justifies it, not the calendar.
- `master` holds the next release. Contributions arrive on short-lived
  branches. Between releases, `NMPC_VERSION_PRERELEASE` on master is `"dev"`
  and the version is the next one planned, so a build of master is never
  mistaken for a release.
- When a fix must ship without what is on master, a maintenance branch
  `release/X.Y` is cut from the tag `vX.Y.Z` and the fix is cherry-picked to it.
- Incompatible work waits on a branch until the maintainers open the next
  major version; it does not merge into master before then.
- A GitHub milestone named for each planned version shows what is meant to be
  in it.

## Release gates

A release candidate or release is made from a commit on which all of these
have passed:

- the [CMake workflow](.github/workflows/cmake.yml), every configuration in
  Debug and Release: the unit tests in `float` and `double`, with
  `EIGEN_NO_MALLOC`, and the allocation test; the runnable examples; the
  consumer fixtures alone, beside UKF, with Eigen unrolling disabled, and from
  the installed package alone and beside an installed UKF;
- the Cortex-M7 job: no double-precision routine linked, and instruction
  counts within 3 % of the accepted baseline, or a new baseline whose change
  is explained in the release notes;
- for changes to solver mathematics since the last release, the sfwa/models
  flights compared with the last release;
- sfwa/models built and tested against the candidate, with the UKF release it
  will be used with.

## Making a release

1. Check the gates above on the commit to release.
2. In one pull request:
   - set the version in [`include/NMPC/Version.h`](include/NMPC/Version.h):
     the numbers, `NMPC_VERSION_PRERELEASE` (`"rc.1"` for a candidate, `""`
     for a release) and `NMPC_VERSION_STRING` to match. CMake reads the
     version from this header and refuses to configure if the string and the
     numbers disagree, so the package version cannot drift from it;
   - rename `## [Unreleased]` in [CHANGELOG.md](CHANGELOG.md) to
     `## [X.Y.Z] - YYYY-MM-DD` (or `X.Y.Z-rc.N`) and start a new empty
     `## [Unreleased]`;
   - check that the README describes this release.
3. After it merges, tag that commit with an annotated tag and push the tag:

   ```sh
   git tag -a vX.Y.Z -m "NMPC X.Y.Z"
   git push origin vX.Y.Z
   ```

   The tag push runs the CMake workflow again, including a check that the tag
   names the version in the header and that the changelog has its section.
4. Create a GitHub Release from the tag. Mark a candidate as a pre-release.
   The notes are curated, not generated: start from the changelog section and
   add known limitations, migration guidance, the numerical-change entries,
   the UKF release it was tested with, and the contributors to credit.
5. Close the milestone. On master, set the version to the next planned
   version with `NMPC_VERSION_PRERELEASE "dev"`.

A published tag is never moved, deleted and re-created, or force-pushed. A
mistake in a release is fixed by the next release. The documentation for a
release is the repository as of its tag.

### Release notes template

```markdown
## NMPC X.Y.Z

One paragraph: what this release is for and who should upgrade. Tested with UKF A.B.C.

### Changes
(From CHANGELOG.md.)

### Numerical changes
(Each entry: affected configurations, expected effect, evidence, re-fly or retune advice. "None" if none.)

### Migration
(What a user must change, if anything.)

### Known limitations

### Contributors
Thanks to @... for ...
```

## Using a release

Pin a release by its tag. As a subdirectory or with FetchContent, after UKF
so that both use one Eigen:

```cmake
include(FetchContent)
FetchContent_Declare(ukf GIT_REPOSITORY https://github.com/sfwa/ukf.git GIT_TAG v1.0.0)
FetchContent_Declare(nmpc GIT_REPOSITORY https://github.com/sfwa/nmpc.git GIT_TAG v2.0.0)
FetchContent_MakeAvailable(ukf nmpc)
target_link_libraries(controller PRIVATE nmpc::nmpc ukf::single)
```

From an installed copy, ask for the major version you were written against;
the package accepts any later release with the same major version:

```cmake
find_package(nmpc 2.0 REQUIRED CONFIG)
```

As in SemVer, a pre-release comes before its release: `2.0.0-rc.1` does not
meet a request for `2.0` or `2.0.0`, even with `EXACT`. To try a release
candidate, or master, find the package without a version and check
`nmpc_VERSION_PRERELEASE`.

In code, `NMPC_VERSION` is `MAJOR * 10000 + MINOR * 100 + PATCH`, and
`NMPC_VERSION_STRING` includes any pre-release:

```cpp
#include "NMPC/Version.h"
#if NMPC_VERSION < 20100
/* Before 2.1.0. */
#endif
```

## The path to 2.0.0

2.0.0 has a fixed scope: the interfaces above, as they are on master. New
features go to 2.1.0.

1. **Interface review.** Check every covered interface above against its
   documentation, and decide which interfaces are experimental (the table
   above is the proposal).
2. **Installation and packaging**: install rules, the `nmpc` CMake package and
   its version file, and the installed-consumer tests, matching UKF's.
3. **sfwa/models builds against NMPC from a tag**, so that the candidate is
   tested by the controllers it was written for.
4. **`v2.0.0-rc.1`**, then two to four weeks for downstream testing. Release
   blockers found in that time are fixed and a further candidate is made.
5. **`v2.0.0`**, from the last candidate with no blockers.
