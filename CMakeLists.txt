cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(nmpc LANGUAGES CXX)

option(NMPC_USE_SYSTEM_EIGEN "Use an installed Eigen3 package instead of patched Eigen 3.2.8" OFF)
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

add_library(nmpc INTERFACE)
add_library(nmpc::nmpc ALIAS nmpc)
target_compile_features(nmpc INTERFACE cxx_std_17)
target_include_directories(nmpc INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# NMPC and UKF must see one Eigen. When UKF is already configured in this build
# it has fetched the patched Eigen 3.2.8 as the external project target eigen3,
# so NMPC reuses that copy and fetches its own only when it stands alone.
if(NMPC_USE_SYSTEM_EIGEN)
    find_package(Eigen3 3.3 REQUIRED NO_MODULE)
    target_link_libraries(nmpc INTERFACE Eigen3::Eigen)
elseif(TARGET eigen3)
    include(ExternalProject)
    ExternalProject_Get_Property(eigen3 source_dir)
    target_include_directories(nmpc SYSTEM INTERFACE ${source_dir})
    add_dependencies(nmpc eigen3)
else()
    include(ExternalProject)
    set_directory_properties(PROPERTIES EP_PREFIX .)
    ExternalProject_Add(eigen3
        URL https://gitlab.com/libeigen/eigen/-/archive/3.2.8/eigen-3.2.8.tar.gz
        TIMEOUT 30
        # Preserve the allocation and C++17/libc++ compatibility fixes.
        PATCH_COMMAND patch -p1 -t -N < ${CMAKE_CURRENT_SOURCE_DIR}/patch/eigen-3.2.8-no-malloc-in-llt.patch
        COMMAND patch -p1 -t -N < ${CMAKE_CURRENT_SOURCE_DIR}/patch/eigen-3.2.8-cxx17.patch
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
        LOG_DOWNLOAD ON
        LOG_CONFIGURE ON
        LOG_BUILD ON)
    ExternalProject_Get_Property(eigen3 source_dir)
    target_include_directories(nmpc SYSTEM INTERFACE ${source_dir})
    # INTERFACE dependencies propagate to consumers, including parallel builds.
    add_dependencies(nmpc eigen3)
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
                CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 13 AND
                (NOT NMPC_USE_SYSTEM_EIGEN OR Eigen3_VERSION VERSION_EQUAL 3.4.0))
            # GCC diagnoses unreachable packet/QR/LLT paths in these Eigen
            # versions even with SYSTEM includes. Limit this to dependency headers.
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
