cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(nmpc_mismatch LANGUAGES CXX)

# This configure must fail: NMPC is asked for a different Eigen from the one
# UKF, configured first, has downloaded. CI checks that it does.
set(NMPC_TEST_UKF_DIR "" CACHE PATH "A UKF checkout to configure beside NMPC")
set(UKF_EIGEN_VERSION 3.2.8 CACHE STRING "" FORCE)
add_subdirectory(${NMPC_TEST_UKF_DIR} ukf)
set(NMPC_EIGEN_VERSION 3.4.0 CACHE STRING "" FORCE)
add_subdirectory(../../.. nmpc)
