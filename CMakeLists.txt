CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(nmpc CXX)
INCLUDE(ExternalProject)

# CMake module(s) path
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR})

# Set default ExternalProject root directory
SET_DIRECTORY_PROPERTIES(PROPERTIES EP_PREFIX .)

# Add eigen3.2
ExternalProject_Add(
eigen3
URL https://bitbucket.org/eigen/eigen/get/3.2.0.tar.gz
TIMEOUT 300
# Disable install step
INSTALL_COMMAND ""
# Wrap download, configure and build steps in a script to log output
LOG_DOWNLOAD ON
LOG_CONFIGURE ON
LOG_BUILD ON)

ExternalProject_Get_Property(eigen3 source_dir)
SET(eigen_dir ${source_dir})

# Add qpDUNES
ExternalProject_Add(
qpDUNES
URL https://github.com/jfrasch/qpDUNES/archive/master.zip
TIMEOUT 300
# Disable install step
INSTALL_COMMAND ""
# Wrap download, configure and build steps in a script to log output
LOG_DOWNLOAD ON
LOG_CONFIGURE ON
LOG_BUILD ON)

ExternalProject_Get_Property(qpDUNES source_dir)
SET(qpDUNES_dir ${source_dir})

INCLUDE_DIRECTORIES(include ${eigen_dir} ${qpDUNES_dir}/include)

#
# Include directories
#
INCLUDE_DIRECTORIES(include)

ADD_LIBRARY(nmpclib STATIC
	src/state.cpp
	src/dynamics.cpp
	src/nmpc.cpp
	src/ocp.cpp)

ADD_DEPENDENCIES(nmpclib eigen3 qpDUNES)

ADD_SUBDIRECTORY(c EXCLUDE_FROM_ALL)
