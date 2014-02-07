CMAKE_MINIMUM_REQUIRED(VERSION 2.8.7 FATAL_ERROR)
PROJECT(c66nmpc)

INCLUDE_DIRECTORIES(../include ../c qpDUNES)

set(CMAKE_C_FLAGS "-O3 -Weverything -Wno-documentation -Wno-padded -Wno-unknown-pragmas -Wno-float-equal -fPIC")

ADD_LIBRARY(c66nmpc SHARED
    cnmpc.c
    qpDUNES/dual_qp.c
    qpDUNES/matrix_vector.c
    qpDUNES/setup_qp.c
    qpDUNES/stage_qp_solver_clipping.c
    qpDUNES/stage_qp_solver_qpoases.cpp
    qpDUNES/utils.c)
ADD_LIBRARY(fcsnmpc STATIC
    cnmpc.c
    qpDUNES/dual_qp.c
    qpDUNES/matrix_vector.c
    qpDUNES/setup_qp.c
    qpDUNES/stage_qp_solver_clipping.c
    qpDUNES/stage_qp_solver_qpoases.cpp
    qpDUNES/utils.c)
