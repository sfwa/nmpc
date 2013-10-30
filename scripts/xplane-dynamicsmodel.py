import math
import socket
import sys
import time
import ctypes

import nmpc
nmpc.init()
from nmpc import _cnmpc, state

sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
sock.connect(('127.0.0.1', 51000))
sock.sendall("sub sim/operation/override/override_planepath\n")
sock.sendall("set sim/operation/override/override_planepath [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]\n")
sock.sendall("sub sim/flightmodel/position/psi\n")  # yaw
sock.sendall("sub sim/flightmodel/position/theta\n")  # pitch
sock.sendall("sub sim/flightmodel/position/phi\n")  # roll
sock.sendall("sub sim/flightmodel/position/local_x\n")
sock.sendall("sub sim/flightmodel/position/local_y\n")
sock.sendall("sub sim/flightmodel/position/local_z\n")
sock.sendall("extplane-set update_interval 0.02\n")

nmpc.configure_airframe(
    mass=3.8,
    inertia_tensor=[
        2.59e-1, 0, -0.334e-1,
        0, 1.47e-1, 0,
        -0.334e-1, 0, 4.05e-1],
    prop_coeffs=[0.025, 0.00250],
    drag_coeffs=[0.11, 0.00075, 0.4, 0.025, 0.005],
    lift_coeffs=[-3.7, -5.4, 1.3, 1.7, 0.05],
    side_coeffs=[
        0, -2.35e-01, -1.87e-03, 4.53e-04,
        0.0, 1.1e-02, -1.1e-02, 0.0],
    pitch_moment_coeffs=[-0.001, -0.014, 0.0, -0.03, -0.03, 0.0],
    roll_moment_coeffs=[-0.002, 0.0, -0.03, 0.03, 0.0],
    yaw_moment_coeffs=[0, -0.005, 0.0, 0.0, 0.0, 0.0])

TIMESTEP = 1.0/50.0  # 50Hz updates.
_cnmpc.nmpc_set_position(math.radians(-37.8136), math.radians(144.9), 200)
_cnmpc.nmpc_set_velocity(20, 0, 0)
_cnmpc.nmpc_set_acceleration(0, 0, 0)
_cnmpc.nmpc_set_attitude(1, 0, 0, 0)
_cnmpc.nmpc_set_angular_velocity(0, 0, 0)
_cnmpc.nmpc_set_angular_acceleration(0, 0, 0)
_cnmpc.nmpc_set_wind_velocity(0, 0, 0)
_cnmpc.nmpc_get_state(state)

control_vec = [0, 0, 0, 0]

while 1:
    nmpc.integrate(TIMESTEP, (ctypes.c_double * 4)(*control_vec))

    update = "world-set %.9f %.9f %.9f\n" \
        % (math.degrees(state.position[0]),
           math.degrees(state.position[1]),
           state.position[2])

    q = (state.attitude[3], -state.attitude[0], -state.attitude[1], -state.attitude[2])
    yaw = math.atan2(2.0 * (q[0] * q[3] + q[1] * q[2]), 1.0 - 2.0 * (q[2] ** 2.0 + q[3] ** 2.0))
    pitch = math.asin(2.0 * (q[0] * q[2] - q[3] * q[1]))
    roll = math.atan2(2.0 * (q[0] * q[1] + q[2] * q[3]), 1.0 - 2.0 * (q[1] ** 2.0 + q[2] ** 2.0))

    update += "set sim/flightmodel/position/psi %.6f\n" % math.degrees(yaw)
    update += "set sim/flightmodel/position/theta %.6f\n" % math.degrees(pitch)
    update += "set sim/flightmodel/position/phi %.6f\n" % math.degrees(roll)

    sock.sendall(update)
    time.sleep(TIMESTEP)
