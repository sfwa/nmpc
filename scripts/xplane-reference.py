#!/usr/bin/env python
# encoding=utf-8

import os
import sys
import math
import time
import ctypes
import geomag
import vectors
import argparse
import datetime
import functools
import texttable
import collections


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

# Load all the X-Plane data in one go
for line in sys.stdin:
    if line.strip() == "":
        continue
    fields = list(field.strip("\n ") for field in line.split("|") if field.strip("\n "))
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
