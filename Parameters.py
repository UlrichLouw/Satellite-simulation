import math
import numpy as np
pi = math.pi

class Control_Parameters:
    hw_ref = 1

class SET_PARAMS:
    # All the parameters specific to the satellite and the mission
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
    Ix = 14.3 #kg.m^2
    Iy = 13.6 #kg.m^2
    Iz = 4.6 #kg.m^2
    Iw = 88.1e-6 #kgm^2 Inertia of the RW-06 wheel
    wo = (Mean_motion/(3600*24)) / (2*pi*Semi_major_axis * 1000) #rad/s
    k = 2 #order of expansion
    g = -2.2 #defined by IGRF
    h = -8 #
    m = 10
    n = 10
    J_t = 2500000
    wbi = np.array(([1],[1],[1]))
    theta_d_max = 2 #degrees per second (theta derived), angular velocity
    theta_d_d = 0.133 # degrees per second^2 (rotation speed derived), angular acceleration
    h_ws_max = 15.7 # mNms
    N_ws_max = 1.05 #mNm

    # This function provides the Inertia matrix for the satellite
    Inertia = np.identity(3)*np.array((Ix, Iy, Iz))
    hs0 = np.zeros((3,1))
    hs1 = np.array(([0],[0.5],[0.2]))
    ws0 = np.zeros((3,1))
    ws1 = np.array(([0],[-0.003],[0]))

    # Discrete time for control
    Ts = 0.1
    tau_w = 0.25 # wheel inertia


