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

# A request for this major version must be met; the package is source
# compatible within a major version.
find_package(nmpc 2.0 REQUIRED CONFIG)
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
