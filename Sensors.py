import numpy as np
from Earth_model import orbit
from Parameters import SET_PARAMS
from sgp4.api import Satrec, WGS72
from sgp4.api import jday
import math
pi = math.pi

class Sensors:
    def __init__(self):
        self.sat = Satrec()
        self.orbit = orbit()
        self.sat.sgp4init(WGS72, 'i', 5, SET_PARAMS.epoch, 
        SET_PARAMS.Drag_term, 0.0, 0.0, SET_PARAMS.eccentricity, 
        SET_PARAMS.Argument_of_perigee, SET_PARAMS.inclination,
        SET_PARAMS.Mean_anomaly, SET_PARAMS.wo*60, SET_PARAMS.RAAN
        )  
        self.satellite = self.sat.twoline2rv(SET_PARAMS.s_list, SET_PARAMS.t_list)
        e, self.r_sat, self.v_sat = self.satellite.sgp4(SET_PARAMS.J_t, SET_PARAMS.fr)  

    def noise(self):
        pass

    def sun(self, t):
        T_jc = (SET_PARAMS.J_t + SET_PARAMS.fr + t * 3.168808781403e-8 - 2452545)/36525
        M_o = 357.527723300 + 35999.050340*T_jc     #in degrees
        lambda_Mo = 280.460618400 + 36000.770053610*T_jc        #degrees
        lambda_e = lambda_Mo + 1.914666471*np.sin(M_o*pi/180) + 0.019994643*math.sin(2*M_o*pi/180)      #degrees
        epsilon =  23.439291 - 0.013004200*T_jc                 #degrees
        r_o = 1.00140612 - 0.016708617*np.cos(M_o*pi/180) - 0.000139589*np.cos(2*M_o*pi/180)        #degrees
        rsun = r_o * np.array(([np.cos(lambda_e*pi/180)],[np.cos(epsilon*pi/180)*np.sin(lambda_e*pi/180)],[np.sin(epsilon*pi/180)*np.sin(lambda_e*pi/180)]))
        S_EIC = rsun*(149597871)*1000 - np.reshape(self.r_sat, (3,1))
        
        """
        if D_sun < Radius_earth:
            self.in_sun_view = False
        else:
            self.in_sun_view = True
        """
        return S_EIC #, self.in_sun_view     #in m

    def nadir(self):
        vector = np.array(([0,0,0]))
        return True, vector

    def magnetometer(self):
        pass

    def current_height_above_earth(self):
        h = 500e3
        return h

    def satellite_vector(self, t):
        e, r_sat, v_sat = self.satellite.sgp4(SET_PARAMS.J_t, SET_PARAMS.fr + t * 3.168808781403e-8)
        self.r_sat = np.array((r_sat))*1000 # convert r_sat to m
        self.v_sat = np.array((v_sat))*1000 # v_sat to m/s
        self.A_EFC_to_EIC = self.orbit.EFC_to_EIC(t)
        self.r_sat_EIC = np.matmul(self.A_EFC_to_EIC, self.r_sat)
        self.v_sat_EIC = np.matmul(self.A_EFC_to_EIC, self.v_sat)
        self.A_EIC_to_ORC = self.orbit.EIC_to_ORC(self.r_sat_EIC, self.v_sat_EIC)
        self.r_sat = np.matmul(self.A_EIC_to_ORC, self.r_sat_EIC)
        self.v_sat = np.matmul(self.A_EIC_to_ORC, self.v_sat_EIC)
        return self.r_sat, self.v_sat, self.A_EIC_to_ORC
