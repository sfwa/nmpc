/*
Part of NMPC. Distributed under the MIT licence: see LICENSE in the root of
the repository.
*/

#ifndef NMPC_VERSION_H
#define NMPC_VERSION_H

/*
The library version, under Semantic Versioning (https://semver.org). The
rewrite starts at 2: the earlier C, C66x and Python implementation, tagged
legacy, is version 1. RELEASING.md says what each part of the version
promises.

This header is the one place the version is written. CMake reads it for the
nmpc package version (find_package(nmpc 2.0)), and a configure fails if the
string below disagrees with the numbers. NMPC_VERSION_PRERELEASE is empty in a
release and otherwise names the pre-release, such as "rc.1". Between releases
master carries the next planned version with the pre-release "dev".
*/
#define NMPC_VERSION_MAJOR 2
#define NMPC_VERSION_MINOR 0
#define NMPC_VERSION_PATCH 0
#define NMPC_VERSION_PRERELEASE "dev"
#define NMPC_VERSION_STRING "2.0.0-dev"

/* One number for comparisons in the preprocessor: 20000 for 2.0.0, 20103 for 2.1.3. */
#define NMPC_VERSION (NMPC_VERSION_MAJOR * 10000 + NMPC_VERSION_MINOR * 100 + NMPC_VERSION_PATCH)

#endif
