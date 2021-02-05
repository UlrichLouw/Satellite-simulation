import numpy as np
import math
import scipy as sp
from scipy import special

pi = math.pi

def euler_to_quaternion(roll, pitch, yaw):
    qx = np.sin(roll/2) * np.cos(pitch/2) * np.cos(yaw/2) - np.cos(roll/2) * np.sin(pitch/2) * np.sin(yaw/2)
    qy = np.cos(roll/2) * np.sin(pitch/2) * np.cos(yaw/2) + np.sin(roll/2) * np.cos(pitch/2) * np.sin(yaw/2)
    qz = np.cos(roll/2) * np.cos(pitch/2) * np.sin(yaw/2) - np.sin(roll/2) * np.sin(pitch/2) * np.cos(yaw/2)
    qw = np.cos(roll/2) * np.cos(pitch/2) * np.cos(yaw/2) + np.sin(roll/2) * np.sin(pitch/2) * np.sin(yaw/2)

    return [qx, qy, qz, qw]

def quaternion_error(current_quaternion, command_quaternion):
    qc1, qc2, qc3, qc4 = command_quaternion
    q_c = np.array([qc4, qc3, -qc2, -qc1],[-qc3, qc4, qc1, -qc2],[qc2, -qc1, qc4, -qc3], [qc1, qc2, qc3, qc4])
    q_error = np.matmul(q_c, current_quaternion)
    return q_error


class SET_PARAMS:
    eccentricicity = 0.0002704
    inclination = 97.2927 #degrees
    Semi_major_axis = 6879.55 #km The distance from the satellite to the earth + the earth radius
    RAAN = 0 #Right ascension of the ascending node
    AP = 0 #argument of perigee
    Re = 6371.2 #km magnetic reference radius
    Mean_motion = 15.215 #rev/day
    Period = 5678.7 #seconds
    Ix = 0.0071 #kg.m^2
    Iy = 0.035 #kg.m^2
    Iz = 0.337 #kg.m^2
    wo = (Mean_motion/(3600*24)) / (2*pi*Semi_major_axis * 1000) #rad/s
    k = 2 #order of expansion
    g = -2.2 #defined by IGRF
    h = -8 #
    m = 10
    n = 10


class Sensors:
    def sun(self):
        pass

    def nadir(self):
        pass

    def magnetometer(self):
        pass

class Actuators:
    def active_aerodynamic_roll_control_paddles(self):
        self.AARCP_num = 2

    def magnetic_torgue_rods(self):
        self.MTR_num = 3

    def nano_reaction_wheels(self):
        self.RW_num = 3

class Satelite_model:
    def __init__(self):
        self.Ix = SET_PARAMS.Ix
        self.Iy = SET_PARAMS.Iy
        self.Iz = SET_PARAMS.Iz
        self.V = np.zeros((2))
        self.first = 1
    
    def Inertia_func(self):
        Inertia = np.identity(3)*np.array((self.Ix, self.Iy, self.Iz))
        return Inertia

    def AOB_func(self, current_quaternion):
        q1, q2, q3, q4 = current_quaternion
        AOB = np.zeros((3,3))
        AOB[0,0] = q1**2-q2**2-q3**2+q4**2
        AOB[0,1] = 2*(q1*q2 + q3*q4)
        AOB[0,2] = 2*(q1*q3 - q2*q4)
        AOB[1,0] = 2*(q1*q2 - q3*q4)
        AOB[1,1] = -q1**2+q2**2-q3**2+q4**2
        AOB[1,2] = 2*(q2*q3 + q1*q4)
        AOB[2,0] = 2*(q1*q3 + q2*q4)
        AOB[2,1] = 2*(q2*q3 - q1*q4)
        AOB[2,2] = -q1**2-q2**2+q3**2+q4**2
        return AOB

    def q_derived_func(self, current_quaternion, wx, wy, wz, A, wo):
        q1, q2, q3, q4 = current_quaternion
        Q = np.array([[q4, -q3, q2],[q3, q4, -q1],[-q2, q1, q4],[-q1, -q2, -q3]])
        Aw = np.array([[wx + A[0,1]*wo],[wy + A[1,1]*wo],[wz + A[2,1]*wo]])
        q_derived = 0.5 * np.matmul(Q, Aw)
        return q_derived

    def wbo_func(self, wbi, A):
        wbo = wbi - np.matmul(A, np.array([0],[-SET_PARAMS.wo],[0]))
        return wbo

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

    def geomagnetic_field_strength_func(self, V):
        if self.first:
            self.V = [V,V]
            self.first = 0
        self.V[0] = V
        delta_V = -(self.V[0]-self.V[1])
        self.V[1] = V
        return delta_V

    def Sun_Position(self):
        

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

class Control:
    def __init__(self):
        pass

    def Magnetic(Select, Wbo, Qe, Magentic_model, CON, H, Wref, h_wheel):
        Nm = 1
        return Nm

    def Reaction_Wheel(H, Wref, EWbo, EQe, CON, Select):
        hw = 1
        Nw = 1
        return hw, Nw

    def Paddle(Select, Wbo, Qe, CON):
        Paddle_angle = 1
        return paddle_angle

if __name__ == "__main__":
    sim = Simulation()
    sat = Satelite_model()
    current_quaternion = euler_to_quaternion(0.2, 0.3, 0.2)
    AOB = sat.AOB_func(current_quaternion)
    print(sat.scalar_potential_function(SET_PARAMS.Semi_major_axis, SET_PARAMS.inclination-90, 1))
 