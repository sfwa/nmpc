cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(pendulum)

add_executable(pendulum pendulum.cpp)
target_link_libraries(pendulum PRIVATE nmpc::nmpc)
target_compile_definitions(pendulum PRIVATE EIGEN_NO_MALLOC)
nmpc_project_warnings(pendulum)
