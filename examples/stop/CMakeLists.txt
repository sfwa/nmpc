cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(stop)

add_executable(stop stop.cpp)
target_link_libraries(stop PRIVATE nmpc::nmpc)
target_compile_definitions(stop PRIVATE EIGEN_NO_MALLOC)
nmpc_project_warnings(stop)
