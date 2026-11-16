cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(attitude)

add_executable(attitude attitude.cpp)
target_link_libraries(attitude PRIVATE nmpc::nmpc)
target_compile_definitions(attitude PRIVATE EIGEN_NO_MALLOC)
nmpc_project_warnings(attitude)
