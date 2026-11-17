cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# The version is written once, in include/NMPC/Version.h.
include(${CMAKE_CURRENT_LIST_DIR}/cmake/ReadVersion.cmake)
read_version(${CMAKE_CURRENT_LIST_DIR}/include/NMPC/Version.h NMPC nmpc_version nmpc_version_prerelease)
project(nmpc VERSION ${nmpc_version} LANGUAGES CXX)

# NMPC and UKF must see one Eigen. When UKF is configured first in the same
# build, NMPC's Eigen options default to UKF's, and a configure that sets them
# differently fails rather than build a program with two Eigens.
set(nmpc_ukf_configured OFF)
if(DEFINED UKF_USE_SYSTEM_EIGEN)
    set(nmpc_ukf_configured ON)
endif()
set(nmpc_default_system_eigen OFF)
set(nmpc_default_eigen_version "3.2.8")
set(nmpc_default_eigen_url "")
set(nmpc_default_unrolling_limit "")
if(nmpc_ukf_configured)
    set(nmpc_default_system_eigen ${UKF_USE_SYSTEM_EIGEN})
    if(DEFINED UKF_EIGEN_VERSION)
        set(nmpc_default_eigen_version ${UKF_EIGEN_VERSION})
    endif()
    if(DEFINED UKF_EIGEN_URL)
        set(nmpc_default_eigen_url "${UKF_EIGEN_URL}")
    endif()
    if(DEFINED UKF_EIGEN_UNROLLING_LIMIT)
        set(nmpc_default_unrolling_limit "${UKF_EIGEN_UNROLLING_LIMIT}")
    endif()
endif()

option(NMPC_USE_SYSTEM_EIGEN "Use an installed Eigen3 package instead of downloading one" ${nmpc_default_system_eigen})
option(NMPC_WARNINGS_AS_ERRORS "Treat warnings in NMPC's own targets as errors" ON)
set(nmpc_project_source_dir ${CMAKE_CURRENT_SOURCE_DIR})

# A consumer needs only the header target unless it opts into the project tools.
set(nmpc_standalone OFF)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(nmpc_standalone ON)
endif()
option(NMPC_BUILD_TESTS "Configure the NMPC unit tests" ${nmpc_standalone})
option(NMPC_BUILD_BENCHMARKS "Configure the NMPC benchmarks" ${nmpc_standalone})
option(NMPC_BUILD_EXAMPLES "Configure the NMPC examples" ${nmpc_standalone})

option(NMPC_INSTALL "Generate install rules and an nmpc CMake package" ${nmpc_standalone})
set(NMPC_EIGEN_VERSION "${nmpc_default_eigen_version}" CACHE STRING
    "Eigen release to download when NMPC_USE_SYSTEM_EIGEN is OFF (3.2.8 is patched; 3.3 or newer is used as released)")
set(NMPC_EIGEN_URL "${nmpc_default_eigen_url}" CACHE STRING
    "Eigen source archive (URL or local file) to use instead of the GitLab release of NMPC_EIGEN_VERSION")

include(GNUInstallDirs)
add_library(nmpc INTERFACE)
add_library(nmpc::nmpc ALIAS nmpc)
target_compile_features(nmpc INTERFACE cxx_std_17)
target_include_directories(nmpc INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

if(nmpc_ukf_configured)
    if((NMPC_USE_SYSTEM_EIGEN AND NOT UKF_USE_SYSTEM_EIGEN) OR
            (UKF_USE_SYSTEM_EIGEN AND NOT NMPC_USE_SYSTEM_EIGEN))
        message(FATAL_ERROR "NMPC_USE_SYSTEM_EIGEN (${NMPC_USE_SYSTEM_EIGEN}) must match UKF_USE_SYSTEM_EIGEN "
            "(${UKF_USE_SYSTEM_EIGEN}): NMPC and UKF must see one Eigen")
    endif()
    if(NOT NMPC_USE_SYSTEM_EIGEN AND DEFINED UKF_EIGEN_VERSION AND
            (NOT NMPC_EIGEN_VERSION STREQUAL UKF_EIGEN_VERSION OR NOT NMPC_EIGEN_URL STREQUAL UKF_EIGEN_URL))
        message(FATAL_ERROR "NMPC_EIGEN_VERSION and NMPC_EIGEN_URL (${NMPC_EIGEN_VERSION} ${NMPC_EIGEN_URL}) must match "
            "UKF_EIGEN_VERSION and UKF_EIGEN_URL (${UKF_EIGEN_VERSION} ${UKF_EIGEN_URL}): NMPC uses UKF's Eigen")
    endif()
endif()

if(NMPC_USE_SYSTEM_EIGEN)
    # Any installed Eigen 3.3 or newer, including Eigen 5. Point CMake at a
    # particular installation with Eigen3_DIR or CMAKE_PREFIX_PATH. The version
    # is checked here rather than in find_package because Eigen 5's package
    # accepts a single requested version only within the same minor release.
    find_package(Eigen3 REQUIRED NO_MODULE)
    if(Eigen3_VERSION VERSION_LESS 3.3)
        message(FATAL_ERROR "NMPC needs Eigen 3.3 or newer when NMPC_USE_SYSTEM_EIGEN is ON; found ${Eigen3_VERSION}")
    endif()
    message(STATUS "NMPC: using installed Eigen ${Eigen3_VERSION}")
    target_link_libraries(nmpc INTERFACE Eigen3::Eigen)
else()
    include(ExternalProject)
    if(TARGET eigen3)
        # UKF, configured first, has already downloaded the same Eigen.
        message(STATUS "NMPC: using UKF's downloaded Eigen ${NMPC_EIGEN_VERSION}")
    else()
        if(NOT NMPC_EIGEN_VERSION STREQUAL "3.2.8" AND
                (NOT NMPC_EIGEN_VERSION MATCHES "^[0-9]+\\.[0-9]+(\\.[0-9]+)?$" OR NMPC_EIGEN_VERSION VERSION_LESS 3.3))
            message(FATAL_ERROR "NMPC_EIGEN_VERSION must be 3.2.8 or a release of 3.3 or newer, such as 3.4.0 or 5.0.0")
        endif()
        message(STATUS "NMPC: using downloaded Eigen ${NMPC_EIGEN_VERSION}")
        if(NMPC_EIGEN_URL STREQUAL "")
            set(nmpc_eigen_url
                https://gitlab.com/libeigen/eigen/-/archive/${NMPC_EIGEN_VERSION}/eigen-${NMPC_EIGEN_VERSION}.tar.gz)
        else()
            set(nmpc_eigen_url ${NMPC_EIGEN_URL})
        endif()
        if(NMPC_EIGEN_VERSION STREQUAL "3.2.8")
            # Preserve the allocation and C++17/libc++ compatibility fixes.
            set(nmpc_eigen_patch
                PATCH_COMMAND patch -p1 -t -N < ${CMAKE_CURRENT_SOURCE_DIR}/patch/eigen-3.2.8-no-malloc-in-llt.patch
                COMMAND patch -p1 -t -N < ${CMAKE_CURRENT_SOURCE_DIR}/patch/eigen-3.2.8-cxx17.patch)
        else()
            # Newer releases need no patches: the solver's factorisations use
            # NMPC-owned storage, so EIGEN_NO_MALLOC checks still apply.
            set(nmpc_eigen_patch "")
        endif()
        set_directory_properties(PROPERTIES EP_PREFIX .)
        ExternalProject_Add(eigen3
            URL ${nmpc_eigen_url}
            TIMEOUT 30
            ${nmpc_eigen_patch}
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            INSTALL_COMMAND ""
            LOG_DOWNLOAD ON
            LOG_CONFIGURE ON
            LOG_BUILD ON)
    endif()
    ExternalProject_Get_Property(eigen3 source_dir)
    # An installed package carries the downloaded headers with it. CMake
    # does not prefix a relative SYSTEM install path, so do it here unless
    # the include directory is already absolute.
    if(IS_ABSOLUTE "${CMAKE_INSTALL_INCLUDEDIR}")
        set(nmpc_installed_eigen_dir "${CMAKE_INSTALL_INCLUDEDIR}/nmpc/eigen3")
    else()
        set(nmpc_installed_eigen_dir "$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}/nmpc/eigen3")
    endif()
    target_include_directories(nmpc SYSTEM INTERFACE
        $<BUILD_INTERFACE:${source_dir}>
        $<INSTALL_INTERFACE:${nmpc_installed_eigen_dir}>)
    # INTERFACE dependencies propagate to consumers, including parallel builds.
    add_dependencies(nmpc eigen3)
endif()

# Eigen unrolls fixed-size expressions up to EIGEN_UNROLLING_LIMIT (100 by
# default). A lower limit trades speed for flash. Leave it empty for Eigen's
# default. It is defined for every target that links nmpc, which must agree
# with any other code in the program that includes Eigen; beside UKF it
# defaults to, and must equal, UKF_EIGEN_UNROLLING_LIMIT.
set(NMPC_EIGEN_UNROLLING_LIMIT "${nmpc_default_unrolling_limit}" CACHE STRING
    "EIGEN_UNROLLING_LIMIT for targets that use NMPC (empty: Eigen's default of 100; 0 disables unrolling)")
if(nmpc_ukf_configured AND DEFINED UKF_EIGEN_UNROLLING_LIMIT AND
        NOT NMPC_EIGEN_UNROLLING_LIMIT STREQUAL UKF_EIGEN_UNROLLING_LIMIT)
    message(FATAL_ERROR "NMPC_EIGEN_UNROLLING_LIMIT (${NMPC_EIGEN_UNROLLING_LIMIT}) must match "
        "UKF_EIGEN_UNROLLING_LIMIT (${UKF_EIGEN_UNROLLING_LIMIT})")
endif()
if(NOT NMPC_EIGEN_UNROLLING_LIMIT STREQUAL "")
    if(NOT NMPC_EIGEN_UNROLLING_LIMIT MATCHES "^[0-9]+$")
        message(FATAL_ERROR "NMPC_EIGEN_UNROLLING_LIMIT must be empty or a nonnegative integer")
    endif()
    target_compile_definitions(nmpc INTERFACE EIGEN_UNROLLING_LIMIT=${NMPC_EIGEN_UNROLLING_LIMIT})
endif()

# An installed tree is used with find_package(nmpc), giving the same nmpc::nmpc
# target as add_subdirectory.
if(NMPC_INSTALL)
    include(CMakePackageConfigHelpers)
    set_target_properties(nmpc PROPERTIES EXPORT_NAME nmpc)
    install(TARGETS nmpc EXPORT nmpcTargets)
    install(DIRECTORY include/NMPC DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    if(NOT NMPC_USE_SYSTEM_EIGEN)
        # The downloaded (and, for 3.2.8, patched) Eigen headers. Build the
        # project before installing so that they have been downloaded.
        install(DIRECTORY ${source_dir}/Eigen ${source_dir}/unsupported
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nmpc/eigen3)
    endif()
    set(nmpc_cmake_dir ${CMAKE_INSTALL_LIBDIR}/cmake/nmpc)
    install(EXPORT nmpcTargets NAMESPACE nmpc:: DESTINATION ${nmpc_cmake_dir})
    configure_package_config_file(cmake/nmpcConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/nmpcConfig.cmake
        INSTALL_DESTINATION ${nmpc_cmake_dir})
    # Within a major version, later releases are source compatible
    # (RELEASING.md). The package is header-only, so any architecture will do.
    set(nmpc_version_file_options COMPATIBILITY SameMajorVersion)
    if(NOT CMAKE_VERSION VERSION_LESS 3.14)
        list(APPEND nmpc_version_file_options ARCH_INDEPENDENT)
    endif()
    write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/nmpcBasicConfigVersion.cmake
        VERSION ${PROJECT_VERSION} ${nmpc_version_file_options})
    # A pre-release does not pass for its release (cmake/ConfigVersion.cmake.in).
    set(package_name nmpc)
    set(package_version ${PROJECT_VERSION})
    set(package_prerelease "${nmpc_version_prerelease}")
    configure_file(cmake/ConfigVersion.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/nmpcConfigVersion.cmake @ONLY)
    install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/nmpcConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/nmpcConfigVersion.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/nmpcBasicConfigVersion.cmake
        DESTINATION ${nmpc_cmake_dir})
endif()

# Warning policy is private to project executables, never consumers. Project
# targets are also built without C++ exceptions or RTTI: the library reports
# failures through return values and needs neither.
function(nmpc_project_warnings target)
    set_target_properties(${target} PROPERTIES CXX_EXTENSIONS OFF)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4)
        if(NMPC_WARNINGS_AS_ERRORS)
            target_compile_options(${target} PRIVATE /WX)
        endif()
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wshadow -Wconversion
            -Wdouble-promotion -fno-exceptions -fno-rtti)
        if(NMPC_WARNINGS_AS_ERRORS)
            target_compile_options(${target} PRIVATE -Werror)
        endif()
        if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND
                CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 13)
            # GCC diagnoses unreachable packet/QR/LLT paths and possibly
            # uninitialized values inside Eigen even with SYSTEM includes.
            # Limit this to dependency headers.
            target_compile_options(${target} PRIVATE
                "-include${nmpc_project_source_dir}/cmake/EigenGccWarnings.h")
        endif()
    endif()
endfunction()

if(NMPC_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test EXCLUDE_FROM_ALL)
endif()
if(NMPC_BUILD_BENCHMARKS)
    add_subdirectory(benchmark EXCLUDE_FROM_ALL)
endif()
if(NMPC_BUILD_EXAMPLES)
    add_subdirectory(examples EXCLUDE_FROM_ALL)
endif()
