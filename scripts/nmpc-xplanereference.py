import os
import math
import socket
import sys
import time
import ctypes
import vectors
import bisect

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


def interpolate_reference(sample_time, points):
    index = bisect.bisect([p[0] for p in points], sample_time)
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
    state_weights=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    control_weights=[1e-10, 1, 1],
    terminal_weights=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    upper_control_bound=[18000, 1.0, 1.0],
    lower_control_bound=[0, -1.0, -1.0])

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

# Set up the NMPC reference trajectory using correct interpolation.
for i in xrange(0, nmpc.HORIZON_LENGTH):
    horizon_point = interpolate_reference(
        i*nmpc.STEP_LENGTH, xplane_reference_points)
    horizon_point.extend([0, 0, 0])
    nmpc.set_reference(horizon_point[1:], i)

nmpc.initialise_horizon()
for i in range(8):
    nmpc.prepare()
    nmpc.solve(interpolate_reference(0, xplane_reference_points)[1:])