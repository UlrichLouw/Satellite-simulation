import math
import numpy as np
import Quaternion_functions
from sgp4.api import jday

pi = math.pi

class SET_PARAMS:
    # All the parameters specific to the satellite and the mission
    """
    Orbit parameters
    """
    eccentricicity = 0.000092
    inclination = 97.4 #degrees
    Semi_major_axis = 6879.55 #km The distance from the satellite to the earth + the earth radius
    Height_above_earth_surface = 500e3 #distance above earth surface
    Scale_height = 8500 #scale height of earth atmosphere
    RAAN = 275*pi/180 #Right ascension of the ascending node in radians
    AP = 0 #argument of perigee
    Re = 6371.2 #km magnetic reference radius
    Mean_motion = 15.2355 #rev/day
    Mean_anomaly = 29.3 #degrees
    Argument_of_perigee = 57.4 #in degrees
    omega = Argument_of_perigee
    Period = 5678.7 #seconds
    J_t,fr = jday(2021,2,16,15,30,0) #current julian date
    epoch = J_t - 2433281.5
    we = (Mean_motion/(3600*24)) / (2*pi*Semi_major_axis * 1000) #rad/s
    """
    position parameters
    """
    v_sat = np.array(([0,we,0]))
    a_G0 = 0
    """
    Atmosphere (Aerodynamic)
    """
    Drag_term = 0.000194
    normal_accommodation = 0.8
    tangential_accommodation = 0.8
    ratio_of_molecular_exit = 0.05
    offset_vector = np.array(([0.01,0.01,0.01]))
    """
    Earth effects (geomagnetic)
    """
    k = 2 #order of expansion
    g = -2.2 #defined by IGRF
    h = -8 #
    m = 10
    n = 10
    """
    Satellite body
    """
    Mass = 20 #kg
    Dimensions = np.array(([0.3, 0.3, 0.4])) # Lx, Ly, Lz
    Ix = 14.3 #kg.m^2
    Iy = 13.6 #kg.m^2
    Iz = 4.6 #kg.m^2
    Iw = 88.1e-6 #kgm^2 Inertia of the RW-06 wheel
    """
    Satellite initial position
    """
    quaternion_initial = Quaternion_functions.euler_to_quaternion(0,0,0) #roll, pitch, yaw
    wbi = np.array(([0.25],[0.1],[0.3]))
    initial_angular_momentum = np.zeros((3,1))
    position_vector = np.array([0,0,0]) # position vector - Distance from center of earth (R), 
    velocity_vector = np.array([0, we, 0])
    """
    Max parameters of actuaters
    """
    theta_d_max = 2 #degrees per second (theta derived), angular velocity
    theta_d_d = 0.133 # degrees per second^2 (rotation speed derived), angular acceleration
    h_ws_max = 15.7 # mNms
    N_ws_max = 1.05 #mNm
    """
    Control parameters
    """
    Kp = 5
    Kd = 2
    w_ref = np.zeros((3,1))
    q_ref = Quaternion_functions.euler_to_quaternion(0,-45,90) #roll, pitch, yaw
    time = 0
    Ts = 0.01 # Discrete time for control

    


