cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(nmpc_consumer LANGUAGES CXX)

# Exercise consumption from a parent that has not already enabled C++17.
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
enable_testing()

# With NMPC_TEST_UKF_DIR set to a UKF checkout, UKF is configured first, as
# sfwa/models does, and NMPC must take UKF's patched Eigen rather than fetch
# a second copy.
set(NMPC_TEST_UKF_DIR "" CACHE PATH "A UKF checkout to configure beside NMPC")
if(NMPC_TEST_UKF_DIR)
    set(UKF_USE_SYSTEM_EIGEN ${NMPC_USE_SYSTEM_EIGEN} CACHE BOOL "" FORCE)
    add_subdirectory(${NMPC_TEST_UKF_DIR} ukf)
endif()
add_subdirectory(../.. nmpc)

if(TARGET unittest OR TARGET unittest_no_heap OR TARGET benchmark OR TARGET stop OR TARGET attitude)
    message(FATAL_ERROR "NMPC project tools must be opt-in when used as a subdirectory")
endif()
if(NMPC_USE_SYSTEM_EIGEN AND TARGET eigen3)
    message(FATAL_ERROR "System Eigen must not configure the bundled download")
endif()
if(NMPC_TEST_UKF_DIR AND NOT NMPC_USE_SYSTEM_EIGEN)
    get_target_property(nmpc_includes nmpc INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(ukf_includes ukf INTERFACE_INCLUDE_DIRECTORIES)
    list(GET ukf_includes 1 ukf_eigen)
    list(FIND nmpc_includes ${ukf_eigen} shared)
    if(shared LESS 0)
        message(FATAL_ERROR "NMPC must use UKF's Eigen (${ukf_eigen}), not ${nmpc_includes}")
    endif()
endif()

add_executable(consumer_double main.cpp)
target_link_libraries(consumer_double PRIVATE nmpc::nmpc)
target_compile_definitions(consumer_double PRIVATE EXPECT_DOUBLE)
add_executable(consumer_single main.cpp)
target_link_libraries(consumer_single PRIVATE nmpc::nmpc)

# Consumers may build the solver with exceptions disabled.
add_executable(consumer_no_exceptions no_exceptions.cpp)
target_link_libraries(consumer_no_exceptions PRIVATE nmpc::nmpc)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(consumer_no_exceptions PRIVATE -fno-exceptions -fno-rtti)
elseif(MSVC)
    target_compile_definitions(consumer_no_exceptions PRIVATE _HAS_EXCEPTIONS=0)
endif()

# A sibling that does not use NMPC must retain the parent's language and flags.
add_executable(consumer_sibling sibling.cpp)
set(consumers consumer_double consumer_single consumer_no_exceptions consumer_sibling)

# Beside UKF, one program uses both libraries in one precision.
if(NMPC_TEST_UKF_DIR)
    add_executable(consumer_with_ukf with_ukf.cpp)
    target_link_libraries(consumer_with_ukf PRIVATE nmpc::nmpc ukf::single)
    list(APPEND consumers consumer_with_ukf)
endif()

foreach(target ${consumers})
    add_test(NAME ${target} COMMAND ${target})
endforeach()
