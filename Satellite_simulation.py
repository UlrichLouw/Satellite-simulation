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
    eccentricicity = 0.000092
    inclination = 97.4 #degrees
    Semi_major_axis = 6879.55 #km The distance from the satellite to the earth + the earth radius
    RAAN = 275 #Right ascension of the ascending node in degrees
    AP = 0 #argument of perigee
    Re = 6371.2 #km magnetic reference radius
    Mean_motion = 15.2355 #rev/day
    Mean_anomaly = 29.3 #degrees
    Argument_of_perigee = 57.4 #in degrees
    Drag_term = 0.000194
    omega = Argument_of_perigee
    Mass = 20 #kg
    Dimensions = np.array(([0.3, 0.3, 0.4])) # Lx, Ly, Lz
    Period = 5678.7 #seconds
    Ix = 0.4 #kg.m^2
    Iy = 0.45 #kg.m^2
    Iz = 0.3 #kg.m^2
    Iw = 88.1e-6 #kgm^2 Inertia of the RW-06 wheel
    wo = (Mean_motion/(3600*24)) / (2*pi*Semi_major_axis * 1000) #rad/s
    k = 2 #order of expansion
    g = -2.2 #defined by IGRF
    h = -8 #
    m = 10
    n = 10
    J_t = 2500000
    wbi = np.array([1],[1],[1])
    theta_d_max = 2 #degrees per second (theta derived), angular velocity
    theta_d_d = 0.133 # degrees per second^2 (rotation speed derived), angular acceleration
    h_ws_max = 15.7 # mNms
    N_ws_max = 1.05 #mNm


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

    def nano_reaction_wheels(self):     #RW-0.06 Sinclair Interplanetary
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

    def wbo_func(self, A):
        wbo = SET_PARAMS.wbi - np.matmul(A, np.array([0],[-SET_PARAMS.wo],[0]))
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
        T_jc = (SET_PARAMS.J_t - 2452545)/36525
        M_o = 357.527723300 + 35999.050340*T_jc     #in degrees
        lambda_Mo = 280.460618400 + 36000.770053610*T_jc        #degrees
        lambda_e = lambda_Mo + 1.914666471*np.sin(M_o*pi/180) + 0.019994643*math.sin(2*M_o*pi/180)      #degrees
        epsilon =  23.439291 - 0.013004200*T_jc                 #degrees
        r_o = 1.00140612 - 0.016708617*np.cos(M_o*pi/180) - 0.000139589*np.cos(2*M_o*pi/180)        #degrees
        rsun = r_o * np.array(([np.cos(lambda_e*pi/180)],[np.cos(epsilon*pi/180)*np.sin(lambda_e*pi/180)],[np.sin(epsilon*pi/180)*np.sin(lambda_e*pi/180)]))
        return rsun*(149597871)     #in km

class Disturbances:
    def __init__(self):
        pass

    def Gravity_gradient_func(self, AOB):
        zoB = AOB * np.array(([0],[0],[1]))
        kgx = 3 * SET_PARAMS.wo**2 * (SET_PARAMS.Iz - SET_PARAMS.Iy)
        kgy = 3 * SET_PARAMS.wo**2 * (SET_PARAMS.Ix - SET_PARAMS.Iz)
        kgz = 3 * SET_PARAMS.wo**2 * (SET_PARAMS.Iy - SET_PARAMS.Ix)
        Ngg = np.array(([kgx*AOB[1][2]*AOB[2][2]],[kgy*AOB[0][2]*AOB[2][2]],[kgz*AOB[0][2]*AOB[1][2]]))
        return Ngg

    def Aerodynamic(self, v_ab, Surface_area_i, incidence_angle,
                        normal_accommodation, ratio_of_molecular_exit, offset_vector, 
                        unit_inward_normal, h, h_o, H):
        Ai = Surface_area_i
        alpha_i = incidence_angle
        sigma_t = 0.8       # tangential_accommodation
        sigma_n = 0.8  # normal_accommodation
        S = 0.05        # ratio_of_molecular_exit
        r_pi = offset_vector
        n_i = unit_inward_normal
        p = 0.5 * (p_o * np.exp(-(h-h_o)/H))
        va = np.arraty(([0],[0],[w_e])) * r_sat - v_sat
        N_aero = 0
        for i in range(n):
            N_aero = N_aero + p * np.linalg.norm(va)**2 * Ai * np.heaviside(np.cos(alpha_i))*np.cos(alpha_i)*sigma_t*(r_pi) + 
                        (sigma_n * S + (2-sigma_n - sigma_t)*np.cos(alpha_i)*(r_pi))

        return N_aero

    def Wheel_Imbalance(self, rotation_rate):
        # Static imbalance with one mass a distance r from the centre of the flywheel
        # m - mass of imbalance in kg
        # r - distance from centre of flywheel in m
        Us = 0.0208 / 1e-5 #For the RW-0.06 wheels in kg/m; Us = m*r
        omega = rotation_rate #rad/s rotation rate of wheel
        phi_s = 0 #arbitrary phase
        w_z = np.array(([0,0,0])) #np.array 

        F_zs = Us * omega**2 * np.array(([np.sin(omega*t + phi_s)],[np.cos(omega*t + phi_s)], [0]))

        N_zs = np.matmul(w_z, F_zs)

        # Dynamic imbalance with two masses seperated by 180 degrees and distance, d
        # d - width of flywheel in m
        # r - distance from flywheel
        # m - mass in kg
        Ud = 0.0208 / 1e-7 #For the RW-0.06 wheels in kg/m^2; Ud = m*r*d
        N_zd = Ud*omega**2 * np.array(([np.sin(omega*t + phi_s)],[np.cos(omega*t + phi_s)], [0]))

        return N_zs

        
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
    print(Disturbances().Gravity_gradient_func())