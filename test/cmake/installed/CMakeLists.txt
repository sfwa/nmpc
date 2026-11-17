cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(nmpc_installed_consumer LANGUAGES CXX)

# Consume an installed NMPC through its CMake package, from a parent that has
# not enabled C++17. Configure with CMAKE_PREFIX_PATH (or nmpc_DIR) naming the
# install prefix. With NMPC_TEST_INSTALLED_UKF on, an installed UKF is found
# as well and one program uses both, as an application pinned to releases of
# both libraries would.
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
enable_testing()
option(NMPC_TEST_INSTALLED_UKF "Also find an installed UKF and use both in one program" OFF)

# A release meets a request for its major version, since the package is
# source compatible within a major version. A pre-release, such as master's
# 2.0.0-dev, comes before 2.0.0 and is found without a version.
find_package(nmpc REQUIRED CONFIG)
if(nmpc_VERSION_PRERELEASE STREQUAL "")
    find_package(nmpc ${nmpc_VERSION_MAJOR}.0 REQUIRED CONFIG)
endif()

# Ask the installed version file directly, as find_package does, whether it
# accepts a request (PACKAGE_VERSION_COMPATIBLE) and whether it is an exact
# match (PACKAGE_VERSION_EXACT, which find_package ... EXACT requires).
function(expect_request request expected_COMPATIBLE expected_EXACT)
    set(PACKAGE_FIND_VERSION ${request})
    string(REGEX MATCH "^[0-9]+" PACKAGE_FIND_VERSION_MAJOR ${request})
    set(PACKAGE_VERSION_COMPATIBLE FALSE)
    set(PACKAGE_VERSION_EXACT FALSE)
    include(${nmpc_DIR}/nmpcConfigVersion.cmake)
    foreach(result COMPATIBLE EXACT)
        set(got FALSE)
        if(PACKAGE_VERSION_${result})
            set(got TRUE)
        endif()
        if(NOT got STREQUAL expected_${result})
            message(FATAL_ERROR "Asking for nmpc ${request}: PACKAGE_VERSION_${result} is ${got}, expected ${expected_${result}}")
        endif()
    endforeach()
endfunction()
math(EXPR next_major "${nmpc_VERSION_MAJOR} + 1")
expect_request(${next_major}.0 FALSE FALSE)
if(nmpc_VERSION_PRERELEASE STREQUAL "")
    expect_request(${nmpc_VERSION} TRUE TRUE)
else()
    # 2.0.0-rc.1 must not pass for 2.0.0, exactly or otherwise.
    expect_request(${nmpc_VERSION} FALSE FALSE)
    expect_request(${nmpc_VERSION_MAJOR} FALSE FALSE)
endif()
if(NOT TARGET nmpc::nmpc)
    message(FATAL_ERROR "The installed nmpc package must provide nmpc::nmpc")
endif()
set(nmpc_version_string ${nmpc_VERSION})
if(NOT nmpc_VERSION_PRERELEASE STREQUAL "")
    set(nmpc_version_string ${nmpc_VERSION}-${nmpc_VERSION_PRERELEASE})
endif()
message(STATUS "Found nmpc ${nmpc_version_string} at ${nmpc_DIR}")

add_executable(installed_double ../main.cpp)
target_link_libraries(installed_double PRIVATE nmpc::nmpc)
target_compile_definitions(installed_double PRIVATE EXPECT_DOUBLE)
add_executable(installed_single ../main.cpp)
target_link_libraries(installed_single PRIVATE nmpc::nmpc)
add_executable(installed_no_exceptions ../no_exceptions.cpp)
target_link_libraries(installed_no_exceptions PRIVATE nmpc::nmpc)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(installed_no_exceptions PRIVATE -fno-exceptions -fno-rtti)
elseif(MSVC)
    target_compile_definitions(installed_no_exceptions PRIVATE _HAS_EXCEPTIONS=0)
endif()
set(consumers installed_double installed_single installed_no_exceptions)
foreach(target installed_double installed_single)
    target_compile_definitions(${target} PRIVATE "EXPECT_VERSION_STRING=\"${nmpc_version_string}\"")
endforeach()

if(NMPC_TEST_INSTALLED_UKF)
    find_package(ukf REQUIRED CONFIG)
    add_executable(installed_with_ukf ../with_ukf.cpp)
    target_link_libraries(installed_with_ukf PRIVATE nmpc::nmpc ukf::single)
    list(APPEND consumers installed_with_ukf)
endif()

foreach(target ${consumers})
    add_test(NAME ${target} COMMAND ${target})
endforeach()
