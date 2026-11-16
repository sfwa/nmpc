cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(unittest)

# An installed GoogleTest is used where there is one. Otherwise the same
# GoogleTest as UKF's tests is fetched.
find_package(GTest CONFIG QUIET)
find_package(Threads REQUIRED)
if(TARGET GTest::gtest_main)
    set(nmpc_gtest_libraries GTest::gtest_main)
else()
    include(ExternalProject)
    set_directory_properties(PROPERTIES EP_PREFIX .)
    ExternalProject_Add(
        googletest
        URL https://github.com/google/googletest/archive/release-1.7.0.zip
        TIMEOUT 30
        CMAKE_ARGS
            -DCMAKE_BUILD_TYPE:STRING=Release
            -DCMAKE_POLICY_VERSION_MINIMUM:STRING=3.10
        INSTALL_COMMAND ""
        LOG_DOWNLOAD ON
        LOG_CONFIGURE ON
        LOG_BUILD ON)
    ExternalProject_Get_Property(googletest source_dir)
    ExternalProject_Get_Property(googletest binary_dir)
    set(nmpc_gtest_include ${source_dir}/include)
    set(nmpc_gtest_libraries
        ${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}gtest.a
        ${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main.a)
endif()

# Every test runs in float and in double (typed tests), so one binary covers both precisions.
add_executable(unittest
    TestBoxQP.cpp
    TestSolver.cpp
    TestCheck.cpp)
# The allocation test replaces the global operator new, so it is a program of its own.
add_executable(unittest_no_heap
    TestNoHeap.cpp)

foreach(target unittest unittest_no_heap)
    target_compile_definitions(${target} PRIVATE EIGEN_NO_MALLOC)
    nmpc_project_warnings(${target})
    if(nmpc_gtest_include)
        target_include_directories(${target} SYSTEM PRIVATE ${nmpc_gtest_include})
        add_dependencies(${target} googletest)
    endif()
    target_link_libraries(${target} PRIVATE nmpc::nmpc ${nmpc_gtest_libraries} Threads::Threads)
    add_test(${target} ${target})
endforeach()

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS unittest unittest_no_heap)
