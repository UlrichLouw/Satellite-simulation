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

    Satellite_control_model, Satellite_disturbance_model = Control().Reaction_Wheel(SET_PARAMS.Inertia, SET_PARAMS.ws0, SET_PARAMS.hs0, SET_PARAMS.tau_w, SET_PARAMS.Ts)

    #A, B, C, D = Satellite_control_model
    om_bw = 0.2
    Kpd = 2*SET_PARAMS.Inertia*om_bw
    Kp = SET_PARAMS.Inertia*om_bw**2

    Ky = np.array(([Kpd, 2*Kp]))
    Kr = -2*np.array(([Kp]))
    #D_cl = np.zeros(())