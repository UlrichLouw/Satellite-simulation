import numpy as np
import math
from Satellite_model import Satelite_model
from Parameters import SET_PARAMS
import Quaternion_functions
import Sensors
import Actuators
from Controller import Control
import Earth_model
import Disturbances

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
    pass