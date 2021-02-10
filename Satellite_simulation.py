import numpy as np
import math
import scipy as sp
from scipy import special
from scipy import signal
from Satellite_model import Satelite_model
from Parameters import SET_PARAMS
import Quaternion_functions
import Sensors
import Actuators
from Controller import Control
import Earth_model
import Disturbances
import matplotlib
from matplotlib import pyplot as plt
import control

pi = math.pi

class Simulation:
    def __init__(self):
        print('Start Simulation')
        self.altitude = SET_PARAMS().Semi_major_axis

    class ADCS:
        def init(self):
            print('ADCS initiated')

        def step(action):
            print('action taken')

        def reset(self):
            print('reset environment')

if __name__ == "__main__":

    Satellite_control_model, Satellite_disturbance_model = Control().Reaction_Wheel(SET_PARAMS.Inertia, SET_PARAMS.ws1, SET_PARAMS.hs1, SET_PARAMS.tau_w, SET_PARAMS.Ts)
    A = Satellite_control_model.A
    B = Satellite_control_model.B
    C = Satellite_control_model.C
    D = Satellite_control_model.D

    om_bw = 0.2
    Kpd = 2**om_bw*np.array((SET_PARAMS.Inertia[0], SET_PARAMS.Inertia[1], SET_PARAMS.Inertia[2]))
    Kp = om_bw**2*np.array((SET_PARAMS.Inertia[0], SET_PARAMS.Inertia[1], SET_PARAMS.Inertia[2]))

    Ky = np.concatenate((Kpd, 2*Kp), axis = 1)
    Kr = -2*np.array((Kp))
    D_cl = np.zeros((np.shape(C)[0], np.shape(Kr)[1]))
    Closed_loop_satellite_control_model = control.StateSpace(A - np.matmul(np.matmul(B,Ky),C), np.matmul(B,Kr), C, D_cl)

    print(control.pole(Satellite_control_model))
    print(control.zero(Satellite_control_model))

    mag, phase, omega = control.freqresp(Satellite_control_model, [0.1, 1., 10.])

    t, y = step_response(Closed_loop_satellite_control_model)
    T, yout, xout = control.forced_response(Closed_loop_satellite_control_model, u, t)
    print(T)
    print(yout)
    print(xout)

    for j in range(6):
        plt.plot(T, yout[j])
        plt.show()

    T, yout = control.step_response(Satellite_control_model)

    for j in range(6):
        plt.plot(T, yout[j])
        plt.show()