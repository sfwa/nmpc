CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(acado_nmpc CXX)
INCLUDE(ExternalProject)

# CMake module(s) path
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR})

# Set default ExternalProject root directory
SET_DIRECTORY_PROPERTIES(PROPERTIES EP_PREFIX .)

SET(CMAKE_CXX_FLAGS "-O3")

# Add ukf
ExternalProject_Add(
ukf
URL https://github.com/sfwa/ukf/archive/master.zip
TIMEOUT 30
# Disable install step
INSTALL_COMMAND ""
# Wrap download, configure and build steps in a script to log output
LOG_DOWNLOAD ON
LOG_CONFIGURE ON
LOG_BUILD ON)

ExternalProject_Get_Property(ukf source_dir)
SET(ukf_dir ${source_dir})

ExternalProject_Get_Property(ukf binary_dir)
SET(ukf_lib_dir ${binary_dir})

#
# Prerequisites
#
FIND_PACKAGE(ACADO REQUIRED)

#
# Include directories
#
INCLUDE_DIRECTORIES(
	include
	${ACADO_INCLUDE_DIRS}
	${ukf_dir}/include
	${ukf_lib_dir}/src/eigen3)

ADD_LIBRARY(ukflib STATIC IMPORTED)

#
# Build an executable
#
ADD_EXECUTABLE(acado_nmpc
	src/main.cpp)
TARGET_LINK_LIBRARIES(acado_nmpc ${ACADO_SHARED_LIBRARIES})
SET_TARGET_PROPERTIES(acado_nmpc PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

ADD_DEPENDENCIES(acado_nmpc ukf)
TARGET_LINK_LIBRARIES(acado_nmpc ${ukf_lib_dir}/libukf.a)