import numpy as np


def rad2deg(rad):
    return rad / np.pi * 180

def deg2rad(deg):
    return deg / 180 * np.pi

def getRotMat(q):
    q1, q2, q3, q4 = q[:]
    A = np.zeros((3,3))
    A[0,0] = q1**2-q2**2-q3**2+q4**2
    A[0,1] = 2*(q1*q2 + q3*q4)
    A[0,2] = 2*(q1*q3 - q2*q4)
    A[1,0] = 2*(q1*q2 - q3*q4)
    A[1,1] = -q1**2+q2**2-q3**2+q4**2
    A[1,2] = 2*(q2*q3 + q1*q4)
    A[2,0] = 2*(q1*q3 + q2*q4)
    A[2,1] = 2*(q2*q3 - q1*q4)
    A[2,2] = -q1**2-q2**2+q3**2+q4**2
    return A

def getEulerAngles(q):
    m = getRotMat(q)
    test = -m[2, 0]
    if test > 0.99999:
        yaw = 0
        pitch = np.pi / 2
        roll = np.arctan2(m[0, 1], m[0, 2])
    elif test < -0.99999:
        yaw = 0
        pitch = -np.pi / 2
        roll = np.arctan2(-m[0, 1], -m[0, 2])
    else:
        yaw = np.arctan2(m[0, 1], m[1, 1])
        pitch = np.arctan2(m[2, 0], m[2, 2])
        roll = np.arcsin(-m[2, 1])

    yaw = rad2deg(yaw)
    pitch = rad2deg(pitch)
    roll = rad2deg(roll)

    return yaw, pitch, roll
