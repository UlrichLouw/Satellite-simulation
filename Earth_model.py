from Parameters import SET_PARAMS
import math
import numpy as np

pi = math.pi
class orbit:
    def __init__(self):
        self.w_earth = SET_PARAMS.w_earth
        self.a_G0 = SET_PARAMS.a_G0

    def EFC_to_EIC(self, t):
        a_G = self.w_earth * t + self.a_G0   # angle in radians form the greenwich
        A = np.array(([[np.cos(a_G), -np.sin(a_G), 0], [np.sin(a_G), np.cos(a_G), 0], [0,0,1]]))
        return A

    def EIC_to_ORC(self, position_vector, velocity_vector):
        # position vector - Height from center of earth, 
        position_vector = position_vector/np.linalg.norm(position_vector)
        velocity_vector = velocity_vector/np.linalg.norm(velocity_vector)
        c = -position_vector   # position vector must be measured by sensors
        b = np.cross(velocity_vector, position_vector)/(np.linalg.norm(np.cross(velocity_vector, position_vector)))
        a = np.cross(b,c)
        A = np.reshape(np.array(([a],[b],[c])),(3,3))
        return A

class Earth:   
    def __init__(self):
        J_t = SET_PARAMS.J_t
        self.V = np.zeros((2))
        self.first = 1      # This is required for the geomagnetic_field_strength_func to initiate correctly

    def geomagnetic_field_strength_func(self, V):
        if self.first:
            self.V = [V,V]
            self.first = 0
        self.V[0] = V
        delta_V = -(self.V[0]-self.V[1])
        self.V[1] = V

        return delta_V

    def scalar_potential_function(self, orbit_radius, coevelation, longitude):
        rs = orbit_radius
        theta = coevelation
        lambda_ = longitude
        k = SET_PARAMS.k
        g = SET_PARAMS.g
        h = SET_PARAMS.h
        Re = SET_PARAMS.Re
        P = special.lpmn(k,k,theta)[0]
        V = 0
        for n in range(1,k):
            sum_ = 0
            for m in range(n):
                print(P)
                sum_ = sum_ + (g*np.cos(m*lambda_) + h*np.sin(m*lambda_)) * P[m,n]
            
            V = V + (Re/rs)**(n+1) * sum_

        return V