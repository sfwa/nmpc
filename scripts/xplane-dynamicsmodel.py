import math
import socket
import sys
import time

from cnmpc import State, cnmpc

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

cukf.fixedwingdynamics_set_mass(3.8)
cukf.fixedwingdynamics_set_inertia_tensor((ctypes.c_double * 9)(
    2.59e-1, 0, -0.334e-1,
    0, 1.47e-1, 0,
    -0.334e-1, 0, 4.05e-1
))
cukf.fixedwingdynamics_set_prop_coeffs(0.025, 0.00250)
cukf.fixedwingdynamics_set_lift_coeffs((ctypes.c_double * 5)(
    -3.7, -5.4, 1.3, 1.7, 0.05))
cukf.fixedwingdynamics_set_drag_coeffs((ctypes.c_double * 5)(
    0.11, 0.00075, 0.4, 0.025, 0.005))
cukf.fixedwingdynamics_set_side_coeffs((ctypes.c_double * 12)(
    0, -2.35e-01, -1.87e-03, 4.53e-04,
    0.0, 1.1e-02, -1.1e-02, 0.0
    ))
cukf.fixedwingdynamics_set_pitch_moment_coeffs((ctypes.c_double * 6)(
    -0.001, -0.014, 0.0, -0.03, -0.03, 0.0
    ))
cukf.fixedwingdynamics_set_roll_moment_coeffs((ctypes.c_double * 5)(
    -0.002, 0.0, -0.03, 0.03, 0.0
    ))
cukf.fixedwingdynamics_set_yaw_moment_coeffs((ctypes.c_double * 6)(
    0, -0.005, 0.0, 0.0, 0.0, 0.0
    ))

NMPC_STATE = State()
cnmpc.set_position(math.radians(-37.8136), math.radians(144.9), 200)
cnmpc.set_velocity(20, 0, 0)
cnmpc.set_acceleration(0, 0, 0)
cnmpc.set_attitude(1, 0, 0, 0)
cnmpc.set_angular_velocity(0, 0, 0)
cnmpc.set_angular_acceleration(0, 0, 0)
cnmpc.set_wind_velocity(0, 0, 0)

while 1:
    update = "world-set %.9f %.9f %.9f\n" % (math.degrees(readings["pos_lat_rad"]), math.degrees(readings["pos_lng_rad"]), readings["pos_alt"])

    q = (readings["att_w"], -readings["att_x"], -readings["att_y"], -readings["att_z"])
    yaw = math.atan2(2.0 * (q[0] * q[3] + q[1] * q[2]), 1.0 - 2.0 * (q[2] ** 2.0 + q[3] ** 2.0))
    pitch = math.asin(2.0 * (q[0] * q[2] - q[3] * q[1]))
    roll = math.atan2(2.0 * (q[0] * q[1] + q[2] * q[3]), 1.0 - 2.0 * (q[1] ** 2.0 + q[2] ** 2.0))

    update += "set sim/flightmodel/position/psi %.6f\n" % math.degrees(yaw)
    update += "set sim/flightmodel/position/theta %.6f\n" % math.degrees(pitch)
    update += "set sim/flightmodel/position/phi %.6f\n" % math.degrees(roll)

    time.sleep(1.0/50.0)  # 50Hz updates

    print line
    sock.sendall(update)
