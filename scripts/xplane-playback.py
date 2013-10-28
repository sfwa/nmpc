import math
import socket
import sys
import time

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

fields = ["time", "pos_lat_rad", "pos_lng_rad", "pos_alt", "att_x", "att_y", "att_z", "att_w"]
for line in sys.stdin:
    readings = dict(zip(fields, map(float, line.strip("\n").split("\t"))))

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
