cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(benchmark)

# An installed Google Benchmark is used where there is one. Otherwise the same
# version as UKF's benchmarks is fetched.
find_package(benchmark CONFIG QUIET)
find_package(Threads REQUIRED)
if(TARGET benchmark::benchmark)
    set(nmpc_benchmark_libraries benchmark::benchmark)
else()
    include(ExternalProject)
    set_directory_properties(PROPERTIES EP_PREFIX .)
    ExternalProject_Add(
        googlebenchmark
        URL https://github.com/google/benchmark/archive/v1.3.0.zip
        TIMEOUT 30
        CMAKE_ARGS
            -DCMAKE_BUILD_TYPE:STRING=Release
            -DCMAKE_POLICY_VERSION_MINIMUM:STRING=3.10
            -DBENCHMARK_ENABLE_TESTING:BOOL=OFF
            # This old dependency enables -Werror in Release independently of NMPC.
            -DHAVE_CXX_FLAG_WERROR:BOOL=OFF
        INSTALL_COMMAND ""
        LOG_DOWNLOAD ON
        LOG_CONFIGURE ON
        LOG_BUILD ON)
    ExternalProject_Get_Property(googlebenchmark source_dir)
    ExternalProject_Get_Property(googlebenchmark binary_dir)
    set(nmpc_benchmark_include ${source_dir}/include)
    set(nmpc_benchmark_libraries ${binary_dir}/src/${CMAKE_FIND_LIBRARY_PREFIXES}benchmark.a)
endif()

# The target is named benchmark, as UKF's is; the imported library is benchmark::benchmark.
add_executable(nmpc_benchmark
    SolverBenchmark.cpp
    BenchmarkMain.cpp)
set_target_properties(nmpc_benchmark PROPERTIES OUTPUT_NAME benchmark)
add_custom_target(benchmark DEPENDS nmpc_benchmark)
target_compile_definitions(nmpc_benchmark PRIVATE EIGEN_NO_MALLOC)
nmpc_project_warnings(nmpc_benchmark)
if(nmpc_benchmark_include)
    target_include_directories(nmpc_benchmark SYSTEM PRIVATE ${nmpc_benchmark_include})
    add_dependencies(nmpc_benchmark googlebenchmark)
endif()
target_link_libraries(nmpc_benchmark PRIVATE nmpc::nmpc ${nmpc_benchmark_libraries} Threads::Threads)
