# Changelog

User-facing changes to NMPC, newest first. Each release has a section named
for its version; [RELEASING.md](RELEASING.md) says what a version promises and
how numerical changes are reported. Changes on master that are not yet
released are under **Unreleased**.

Sections within a release: **Added**, **Changed**, **Deprecated**,
**Removed**, **Fixed**, and **Numerical changes** (anything that can move a
plan, with the configurations affected, the expected effect, the evidence,
and whether to re-fly or retune).

## [Unreleased]

The first numbered release of NMPC version 2, planned as `2.0.0` after the
candidate `2.0.0-rc.1`. Version 2 is a rewrite: the earlier C, C66x and
Python implementation (version 1, tag `legacy`) is removed and shares no
interface with it.

### Added

- A header-only C++17 library on Eigen: the problem concept, a
  box-constrained Gauss-Newton iLQR solver with analytic derivatives built for
  real-time iteration, the plan as a policy with feedback gains, and
  derivative checks ([#8]).
- Optional per-stage Jacobians (`NMPC::Jacobians::PerStage`), which keep the
  plans and shrink the workspace ([#9]).
- Inner-loop building blocks (`Incremental`, `NoCorrection`,
  `AttitudeCascade`, `Butterworth2`) and the reference table and clock
  (`PathTable`, `Clock`), with a fractional plan advance ([#8]).
  *Experimental in 2.0: see [RELEASING.md](RELEASING.md#experimental-interfaces).*
- `NMPC/Version.h`, with `NMPC_VERSION_MAJOR`, `_MINOR`, `_PATCH`,
  `_PRERELEASE`, `_STRING` and `NMPC_VERSION`. The CMake package carries the
  same version, and `find_package(nmpc 2.0)` accepts any later 2.x.
- Install rules and a CMake package: `find_package(nmpc)` gives `nmpc::nmpc`.
- CMake options matching UKF's: `NMPC_EIGEN_VERSION` and `NMPC_EIGEN_URL` to
  download another Eigen (3.3 or newer, unpatched), `NMPC_USE_SYSTEM_EIGEN`
  for an installed Eigen 3.3 to 5, `NMPC_EIGEN_UNROLLING_LIMIT` and
  `NMPC_INSTALL`. Beside UKF, NMPC's Eigen options default to UKF's, and a
  configure that would give the program two Eigens fails.
- CI on the same configurations as UKF (GCC 13 and 15, Clang 18 and 22, macOS,
  Eigen 3.2.8, 3.4.0 and 5.0.0), with installed-package consumers, and a
  Cortex-M7 instruction-count gate under QEMU ([#8]).

### Fixed

Review findings fixed before the first release ([#9], [#10], [#11]): the box
QP holds commands whose bounds meet, stops on a relative tolerance and limits
a step at any bound its target crosses; `PathTable::span` clamps to its
window; `AttitudeCascade` compares rates in the vehicle's axes; `advance(N)`
keeps the last command and gain; an uphill or non-finite predicted decrease
is not reported as converged; gains reach the plan only after a whole valid
backward pass.

[#8]: https://github.com/sfwa/nmpc/pull/8
[#9]: https://github.com/sfwa/nmpc/pull/9
[#10]: https://github.com/sfwa/nmpc/pull/10
[#11]: https://github.com/sfwa/nmpc/pull/11
