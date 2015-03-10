import os
import math
import socket
import sys
import time
import ctypes
import vectors
import bisect
import datetime
import copy

import nmpc
nmpc.init(implementation="c")
from nmpc import _cnmpc, state

def socket_readlines(socket):
    buf = socket.recv(4096)
    done = False
    while 1:
        if "\n" in buf:
            (line, buf) = buf.split("\n", 1)
            yield line + "\n"
        else:
            break


def values_to_dict(fields, values):
    out = {}
    for i, field in enumerate(fields):
        if field not in out:
            out[field] = values[i]

    return out


def convert_metric(fields, values):
    metric_values = []
    for field, value in zip(fields, values):
        unit = field.rpartition(",")[2].strip("_")
        if unit == "ktas":
            metric_values.append(float(value) * 1.852 / 3.6)
        elif unit == "ftmsl":
            metric_values.append(float(value) * 0.3048)
        elif unit == "lb":
            metric_values.append(float(value) * 0.45359237 * 9.80665)
        elif unit == "ftlb":
            metric_values.append(float(value) * 1.3558179483314004)
        elif unit == "deg":
            metric_values.append(math.radians(float(value)))
        else:
            metric_values.append(float(value))

    return metric_values


def euler_to_q(yaw, pitch, roll):
    return (vectors.Q.rotate("X", -roll) *
            vectors.Q.rotate("Y", -pitch) *
            vectors.Q.rotate("Z", -yaw))


def interpolate_reference(sample_time, points):
    index = bisect.bisect([p[0] for p in points], sample_time)
    if index >= len(points):
        return points[-1]
    delta = [b - a for a, b in zip(points[index-1], points[index])]
    residual_time = ((sample_time - points[index-1][0]) / delta[0])
    new_point = [(a + b*residual_time) \
        for a, b in zip(points[index-1], delta)]

    q = vectors.Q(new_point[7], new_point[8], new_point[9], new_point[10])
    q_norm = q.normalize()
    new_point[7] = q_norm[0]
    new_point[8] = q_norm[1]
    new_point[9] = q_norm[2]
    new_point[10] = q_norm[3]

    return new_point

# Load the data
headers = None
initialised = False
initial_time = 0.0

MAX_THROTTLE = 25000.0

nmpc.setup(
    state_weights=[1, 1, 1, 1, 1, 1, 1, 1, 1e1, 7e-1, 7e-1, 1e1],
    control_weights=[1e-1, 1e3, 1e3],
    terminal_weights=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    upper_control_bound=[1.0, 1.0, 1.0],
    lower_control_bound=[0, 0, 0])
nmpc.initialise_horizon()

xplane_reference_points = []

fields = ["time",
    "pos_x", "pos_y", "pos_z"
    "vel_n", "vel_e", "vel_d",
    "att_x", "att_y", "att_z", "att_w",
    "angvel_p", "angvel_q", "angvel_r"]
# Load all the X-Plane data in one go
for line in sys.stdin:
    if line.strip() == "":
        continue
    fields = list(field.strip("\n ") \
        for field in line.split("|") if field.strip("\n "))
    if fields[0] == "_real,_time":
        headers = fields
    else:
        data = values_to_dict(headers, fields)
        if not initialised:
            initial_time = float(data["_real,_time"])
            position_offset = [
                float(data["____X,____m"]),
                float(data["____Y,____m"]),
                float(data["____Z,____m"])]
            initialised = True
        attitude = euler_to_q(
            math.radians(float(data["hding,_true"])),
            math.radians(float(data["pitch,__deg"])),
            math.radians(float(data["_roll,__deg"])))
        velocity = (
            -float(data["___vZ,__m/s"]),
            float(data["___vX,__m/s"]),
            -float(data["___vY,__m/s"]))
        out = [
            float(data["_real,_time"]) - initial_time,
            -(float(data["____Z,____m"]) - position_offset[2]),
            (float(data["____X,____m"]) - position_offset[0]),
            -(float(data["____Y,____m"]) - position_offset[1]),
            velocity[0],
            velocity[1],
            velocity[2],
            attitude[0],
            attitude[1],
            attitude[2],
            attitude[3],
            float(data["____P,rad/s"]),
            float(data["____Q,rad/s"]),
            float(data["____R,rad/s"])]

        #print "\t".join(str(o) for o in out)

        xplane_reference_points.append(map(float, out))

# Disable X-Plane simulation and set up initial position.
sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
sock.connect(('127.0.0.1', 51000))
sock.sendall("sub sim/operation/override/override_planepath\n")
sock.sendall("set sim/operation/override/override_planepath [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]\n")
sock.sendall("sub sim/operation/override/override_control_surfaces\n")
sock.sendall("set sim/operation/override/override_control_surfaces 1\n")
sock.sendall("sub sim/operation/override/override_throttles\n")
sock.sendall("set sim/operation/override/override_throttles 1\n")
sock.sendall("sub sim/flightmodel/position/q\n")
sock.sendall("sub sim/flightmodel/position/psi 0\n")  # yaw
sock.sendall("sub sim/flightmodel/position/theta 0\n")  # pitch
sock.sendall("sub sim/flightmodel/position/phi 0\n")  # roll
sock.sendall("sub sim/flightmodel/position/P 0\n")
sock.sendall("sub sim/flightmodel/position/Q 0\n")
sock.sendall("sub sim/flightmodel/position/R 0\n")
sock.sendall("sub sim/flightmodel/position/local_x 0\n")
sock.sendall("sub sim/flightmodel/position/local_y 0\n")
sock.sendall("sub sim/flightmodel/position/local_z 0\n")
sock.sendall("sub sim/flightmodel/position/local_vx 0\n")
sock.sendall("sub sim/flightmodel/position/local_vy 0\n")
sock.sendall("sub sim/flightmodel/position/local_vz 0\n")
sock.sendall("sub sim/flightmodel/engine/ENGN_thro_use\n")
sock.sendall("sub sim/flightmodel/controls/wing1l_ail1def\n")
sock.sendall("sub sim/flightmodel/controls/wing1r_ail1def\n")
sock.sendall("sub sim/weather/wind_now_x_msc\n")
sock.sendall("sub sim/weather/wind_now_y_msc\n")
sock.sendall("sub sim/weather/wind_now_z_msc\n")
sock.sendall("extplane-set update_interval 0.02\n")

sock.recv(1024)
sock.sendall("world-set -37.8136 144.9 200\n")
position_offset = [0, 0, 0]
time.sleep(1.0)

try:
    for line in socket_readlines(sock):
        if line.find("local_x") >= 0:
            position_offset[1] = float(line.split(" ")[-1])
        elif line.find("local_y") >= 0:
            position_offset[2] = -float(line.split(" ")[-1])
        elif line.find("local_z") >= 0:
            position_offset[0] = -float(line.split(" ")[-1])
except socket.error:
    pass

# Set up the NMPC reference trajectory using correct interpolation.
for i in xrange(0, nmpc.HORIZON_LENGTH+1):
    horizon_point = [a for a in interpolate_reference(
        i*nmpc.STEP_LENGTH, xplane_reference_points)]
    horizon_point.extend([0.5, 0.5, 0.5])
    nmpc.set_reference(horizon_point[1:], i)

# Set up initial attitude, velocity and angular velocity.
initial_point = interpolate_reference(0, xplane_reference_points)
update = ""
q = (initial_point[10], -initial_point[7], -initial_point[8], -initial_point[9])
yaw = math.atan2(2.0 * (q[0] * q[3] + q[1] * q[2]), 1.0 - 2.0 * (q[2] ** 2.0 + q[3] ** 2.0))
pitch = math.asin(2.0 * (q[0] * q[2] - q[3] * q[1]))
roll = math.atan2(2.0 * (q[0] * q[1] + q[2] * q[3]), 1.0 - 2.0 * (q[1] ** 2.0 + q[2] ** 2.0))

# Need to calculate the X-Plane quaternion to set the heading, pitch and roll
# as used by the physics model.
xplane_q = [0, 0, 0, 1]
psi = yaw / 2.0
theta = pitch / 2.0
phi = roll / 2.0
xplane_q[0] = math.cos(psi) * math.cos(theta) * math.cos(phi) + math.sin(psi) * math.sin(theta) * math.sin(phi)
xplane_q[1] = math.cos(psi) * math.cos(theta) * math.sin(phi) - math.sin(psi) * math.sin(theta) * math.cos(phi)
xplane_q[2] = math.cos(psi) * math.sin(theta) * math.cos(phi) + math.sin(psi) * math.cos(theta) * math.sin(phi)
xplane_q[3] = -math.cos(psi) * math.sin(theta) * math.sin(phi) + math.sin(psi) * math.cos(theta) * math.cos(phi)
update += "set sim/flightmodel/position/q [%.6f,%.6f,%.6f,%.6f]\n" % (xplane_q[0], xplane_q[1], xplane_q[2], xplane_q[3])
update += "set sim/flightmodel/position/psi %.6f\n" % math.degrees(yaw)
update += "set sim/flightmodel/position/theta %.6f\n" % math.degrees(pitch)
update += "set sim/flightmodel/position/phi %.6f\n" % math.degrees(roll)

update += "set sim/flightmodel/position/local_vx %.6f\n" % initial_point[5]
update += "set sim/flightmodel/position/local_vy %.6f\n" % -initial_point[6]
update += "set sim/flightmodel/position/local_vz %.6f\n" % -initial_point[4]

update += "set sim/flightmodel/position/P %.6f\n" % initial_point[11]
update += "set sim/flightmodel/position/Q %.6f\n" % initial_point[12]
update += "set sim/flightmodel/position/R %.6f\n" % initial_point[13]

# Zero controls.
update += "set sim/flightmodel/engine/ENGN_thro [%.6f,0,0,0,0,0,0,0]\n" % 0
update += "set sim/flightmodel/controls/wing1l_ail1def %.6f\n" % 0
update += "set sim/flightmodel/controls/wing1r_ail1def %.6f\n" % 0

sock.sendall(update)

time.sleep(1.0)

if initial_point[10] < 0:
    initial_point[7] = -initial_point[7]
    initial_point[8] = -initial_point[8]
    initial_point[9] = -initial_point[9]
    initial_point[10] = -initial_point[10]

measured_state = {
    "x": 0,
    "y": 0,
    "z": 0,
    "vx": initial_point[4],
    "vy": initial_point[5],
    "vz": initial_point[6],
    "qx": initial_point[7],
    "qy": initial_point[8],
    "qz": initial_point[9],
    "qw": initial_point[10],
    "wx": initial_point[11],
    "wy": initial_point[12],
    "wz": initial_point[13]
}

wind_velocity = [0, 0, 0]

# Prepare to re-enable the flight model.
update = "set sim/operation/override/override_planepath [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]\n"

sock.setblocking(0)

for i in xrange(1000):
    start_iteration = datetime.datetime.now()

    nmpc.prepare()

    # Get latest "measured" data.
    try:
        for line in socket_readlines(sock):
            fields = line.split(" ")
            if len(fields) != 3:
                continue
            if fields[1] == "sim/flightmodel/position/local_x":
                measured_state["y"] = float(fields[2]) - position_offset[1]
            elif fields[1] == "sim/flightmodel/position/local_y":
                measured_state["z"] = -float(fields[2]) - position_offset[2]
            elif fields[1] == "sim/flightmodel/position/local_z":
                measured_state["x"] = -float(fields[2]) - position_offset[0]
            elif fields[1] == "sim/flightmodel/position/local_vx":
                measured_state["vy"] = float(fields[2])
            elif fields[1] == "sim/flightmodel/position/local_vy":
                measured_state["vz"] = -float(fields[2])
            elif fields[1] == "sim/flightmodel/position/local_vz":
                measured_state["vx"] = -float(fields[2])
            elif fields[1] == "sim/flightmodel/position/psi":
                yaw = math.radians(float(fields[2]))
            elif fields[1] == "sim/flightmodel/position/theta":
                pitch = math.radians(float(fields[2]))
            elif fields[1] == "sim/flightmodel/position/phi":
                roll = math.radians(float(fields[2]))
            elif fields[1] == "sim/flightmodel/position/P":
                measured_state["wx"] = math.radians(float(fields[2]))
            elif fields[1] == "sim/flightmodel/position/Q":
                measured_state["wy"] = math.radians(float(fields[2]))
            elif fields[1] == "sim/flightmodel/position/R":
                measured_state["wz"] = math.radians(float(fields[2]))
            elif fields[1] == "sim/weather/wind_now_x_msc":
                wind_velocity[1] = float(fields[2])
            elif fields[1] == "sim/weather/wind_now_y_msc":
                wind_velocity[2] = -float(fields[2])
            elif fields[1] == "sim/weather/wind_now_z_msc":
                wind_velocity[0] = -float(fields[2])
    except socket.error:
        pass

    nmpc.set_wind_velocity(wind_velocity)

    # Recalculate quaternion in case euler angles have been updated.
    attitude = euler_to_q(yaw, pitch, roll)
    if attitude[3] < 0:
        measured_state["qx"] = -attitude[0]
        measured_state["qy"] = -attitude[1]
        measured_state["qz"] = -attitude[2]
        measured_state["qw"] = -attitude[3]
    else:
        measured_state["qx"] = attitude[0]
        measured_state["qy"] = attitude[1]
        measured_state["qz"] = attitude[2]
        measured_state["qw"] = attitude[3]

    state = [
        measured_state["x"],
        measured_state["y"],
        measured_state["z"],
        measured_state["vx"],
        measured_state["vy"],
        measured_state["vz"],
        measured_state["qx"],
        measured_state["qy"],
        measured_state["qz"],
        measured_state["qw"],
        measured_state["wx"],
        measured_state["wy"],
        measured_state["wz"]]
    print state

    nmpc.solve(state)

    control_vec = nmpc.get_controls()
    print ("t: %.2f " % (i*nmpc.STEP_LENGTH)) + repr(control_vec)

    update += "set sim/flightmodel/engine/ENGN_thro_use [%.6f,0,0,0,0,0,0,0]\n" % control_vec[0]
    update += "set sim/flightmodel/controls/wing1l_ail1def %.6f\n" % math.degrees(control_vec[1] - 0.5)
    update += "set sim/flightmodel/controls/wing1r_ail1def %.6f\n" % math.degrees(control_vec[2] - 0.5)

    # Add one to the index because of the terminal point.
    horizon_point = [a for a in interpolate_reference(
        (i+1+nmpc.HORIZON_LENGTH)*nmpc.STEP_LENGTH, xplane_reference_points)]
    horizon_point.extend([0.5, 0.5, 0.5])
    nmpc.update_horizon(horizon_point[1:])

    sock.sendall(update)

    update = ""

    s = nmpc.STEP_LENGTH - (datetime.datetime.now() - start_iteration).total_seconds()
    if s > 0:
        time.sleep(s)

sock.sendall("set sim/operation/override/override_planepath [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]\n")
