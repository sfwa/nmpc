CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(acado_nmpc CXX)
INCLUDE(ExternalProject)

# CMake module(s) path
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR})

# Set default ExternalProject root directory
SET_DIRECTORY_PROPERTIES(PROPERTIES EP_PREFIX .)

# Add eigen3.2
ExternalProject_Add(
eigen3
URL https://bitbucket.org/eigen/eigen/get/3.2.0.tar.gz
TIMEOUT 30
# Disable install step
INSTALL_COMMAND ""
# Wrap download, configure and build steps in a script to log output
LOG_DOWNLOAD ON
LOG_CONFIGURE ON
LOG_BUILD ON)

ExternalProject_Get_Property(eigen3 source_dir)
SET(eigen_dir ${source_dir})

INCLUDE_DIRECTORIES(include ${eigen_dir})

#
# Prerequisites
#
FIND_PACKAGE(ACADO REQUIRED)

#
# Include directories
#
INCLUDE_DIRECTORIES(
	include
	${ACADO_INCLUDE_DIRS})

ADD_LIBRARY(nmpclib STATIC
	src/nmpc/state.cpp
	src/nmpc/dynamics.cpp)

ADD_DEPENDENCIES(nmpclib eigen3)

ADD_SUBDIRECTORY(c EXCLUDE_FROM_ALL)

#
# Build an executable
#
ADD_EXECUTABLE(acado_nmpc
	src/main.cpp)
TARGET_LINK_LIBRARIES(acado_nmpc ${ACADO_SHARED_LIBRARIES})
SET_TARGET_PROPERTIES(acado_nmpc PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

ADD_DEPENDENCIES(acado_nmpc nmpclib)
TARGET_LINK_LIBRARIES(acado_nmpc nmpclib)