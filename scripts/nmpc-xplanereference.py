import os
import math
import socket
import sys
import time
import ctypes
import vectors

import nmpc
nmpc.init()
from nmpc import _cnmpc, state

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


# Load the data
headers = None
initialised = False
initial_time = 0.0

nmpc.configure_airframe(
    mass=3.8,
    inertia_tensor=[
        2.59e-1, 0, -0.334e-1,
        0, 1.47e-1, 0,
        -0.334e-1, 0, 4.05e-1],
    prop_coeffs=[0.025, 0.00250],
    drag_coeffs=[0.0, 0.0, 0.2, 0.0, 0.05],
    lift_coeffs=[-3.7, -5.4, 1.3, 1.7, 0.05],
    side_coeffs=[
        0, 2.35e-01, -1.87e-03, 4.53e-04,
        0.0, 1.1e-02, -1.1e-02],
    pitch_moment_coeffs=[-0.01, -0.0018, 0.0, -0.001, -0.001],
    roll_moment_coeffs=[-0.002, 0.0, -0.003, 0.003],
    yaw_moment_coeffs=[0, -0.005, 0.0, 0.0, 0.0])

nmpc.setup(
    state_weights=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0],
    control_weights=[0, 0, 0],
    terminal_weights=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0],
    upper_control_bound=[15000, 1.0, 1.0],
    lower_control_bound=[0, -1.0, -1.0])

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

        print "\t".join(str(o) for o in out)

        readings = dict(zip(fields, map(float, out)))