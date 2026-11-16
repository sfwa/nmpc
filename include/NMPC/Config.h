/*
Part of NMPC. Distributed under the MIT licence: see LICENSE in the root of
the repository.
*/

#ifndef NMPC_CONFIG_H
#define NMPC_CONFIG_H

/*
The library version, as major, minor and patch numbers. The rewrite starts at
2: the earlier C, C66x and Python implementation is version 1.
*/
#define NMPC_VERSION_MAJOR 2
#define NMPC_VERSION_MINOR 0
#define NMPC_VERSION_PATCH 0

/*
Assertions on the caller's use of the library, such as a stage index out of
range. They are off unless the application defines NMPC_ASSERT before the
first NMPC header, for example as its own fault handler. The library never
asserts on numerical outcomes: those are reported in the solver's result.
*/
#ifndef NMPC_ASSERT
#define NMPC_ASSERT(condition) ((void)0)
#endif

/*
NMPC never defines real_t. That name belongs to UKF, which selects it with
UKF_SINGLE_PRECISION or UKF_DOUBLE_PRECISION. Every NMPC type is a template on
its scalar T instead, so one program can hold solvers in both precisions.
*/

#endif
