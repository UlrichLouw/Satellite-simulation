import numpy as np
from Earth_model import orbit
from Parameters import SET_PARAMS
from sgp4.api import Satrec, WGS72
from sgp4.api import jday

sat = Satrec()

class Sensors:
    def __init__(self):
        self.orbit = orbit
        self.EFC_to_EIC = self.orbit.EFC_to_EIC(SET_PARAMS.time)
        self.EIC_to_ORC = self.orbit.EIC_to_ORC()
        sat.sgp4init(WGS72, 'i', 5, SET_PARAMS.epoch, 
        SET_PARAMS.Drag_term, 0.0, 0.0, SET_PARAMS.eccentricicity, 
        SET_PARAMS.Argument_of_perigee, SET_PARAMS.inclination,
        SET_PARAMS.Mean_anomaly, SET_PARAMS.we*60, SET_PARAMS.RAAN
        )      

    def noise(self):
        pass

    def sun(self):
        pass

    def nadir(self):
        vector = np.array(([0,0,0]))
        return True, vector

    def magnetometer(self):
        pass

    def current_height_above_earth(self):
        h = 500e3
        return h

    def satellite_vector(self, t):
        r_sat = sat.sgp4(SET_PARAMS.J_t, SET_PARAMS.fr + t * 3.168808781403e-8)
        return r_sat