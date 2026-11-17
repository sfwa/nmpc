# Contributing to NMPC

Thank you for helping. A contribution needs three things, and the maintainers
do the rest:

1. **A clear explanation** in the pull request: what changes, why, and how it
   was checked.
2. **Regression coverage** that fits the change. Every test runs in `float`
   and `double`; a bug fix should come with a test that failed before it, and
   a new problem member or derivative should go through [`Check.h`](include/NMPC/Check.h).
3. **Documentation** for anything a user can see: the README or the header
   comment, as fits.

The pull request template asks for a one-line changelog entry and what you
think the change does to compatibility. A best guess is fine. The maintainers
decide how a change is classified, which release it goes in and what version
that release gets ([RELEASING.md](RELEASING.md) says how), and write
[CHANGELOG.md](CHANGELOG.md) from the entries. You do not need to edit the
changelog, choose a version or tidy your commit history: pull requests are
usually squashed when merged.

## Build and test

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target unittest unittest_no_heap benchmark stop pendulum attitude
ctest --test-dir build
```

The [README](README.md#building) has the details, including how to use
another Eigen. CI builds every pull request with the compilers and Eigen
versions listed under [supported configurations](RELEASING.md#supported-configurations),
with warnings as errors, checks the consumer fixtures in
[test/cmake](test/cmake), and counts Cortex-M7 instructions against an
accepted baseline. It is fine to open a pull request and let CI find what
your machine did not.

## Formatting and lint

Code is formatted with clang-format 18 in Google style (`.clang-format`), and
clang-tidy 18 checks it for bugs and performance problems (`.clang-tidy`). CI
runs both. Before pushing, format with:

```sh
git ls-files '*.h' '*.cpp' '*.c' | xargs clang-format-18 -i
```

[AGENTS.md](AGENTS.md) has the clang-tidy command, and the repository's
conventions for anyone, person or coding agent, making a change.

## Changes that move plans or counts

If your change can move a plan or the work a solve does, say so in the pull
request, even when it is a fix: which configurations it affects and how you
checked the size of the effect. If the Cortex-M7 gate fails because counts
moved by more than 3 %, explain why in the pull request; the maintainers
accept a new baseline (`benchmark/cortex-m7/accepted.json`) when the change
is worth it.

## Interfaces

The interfaces covered by the [compatibility contract](RELEASING.md#the-compatibility-contract)
cannot change incompatibly outside a major release. If your change needs to,
open an issue first so that the maintainers can plan it; an optional problem
member or an option whose default keeps today's behaviour is usually
possible. Code in `NMPC::detail` is not covered and can change freely, and
the [experimental](RELEASING.md#experimental-interfaces) inner-loop and
reference headers may change in a minor release.
