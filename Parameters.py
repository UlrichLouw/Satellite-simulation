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
    eccentricity = 0.000092  # Update eccentricity list
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
    J_t,fr = jday(2020,2,16,15,30,0) #current julian date
    epoch = J_t - 2433281.5 + fr
    Drag_term = 0.000194 # Remember to update the list term
    wo = (Mean_motion/(3600*24)) / (2*pi*Semi_major_axis * 1000) #rad/s
    """
    TLE data
    """
    # s list
    satellite_number_list = '1 25544U'
    international_list = ' 98067A   '
    epoch_list = str("{:.8f}".format(epoch))
    mean_motion_derivative_first_list = '  .00000000'
    mean_motion_derivative_second_list = '  00000-0'
    Drag_term_list = '  19400-4' # B-star
    Ephereris_list = ' 0'
    element_num_checksum_list = '   703'
    s_list = satellite_number_list + international_list + epoch_list + mean_motion_derivative_first_list + mean_motion_derivative_second_list + Drag_term_list + Ephereris_list + element_num_checksum_list
    # r list
    line_and_satellite_number_list = '2 27843  '
    inclination_list = str("{:.4f}".format(inclination))
    intermediate_list = '  '
    RAAN_list = str("{:.4f}".format(RAAN*180/pi))
    intermediate_list_2 = ' '
    eccentricity_list = '0000920 '
    perigree_list = str("{:.5f}".format(Argument_of_perigee))
    intermediate_list_3 = intermediate_list_2
    mean_anomaly_list = str("{:.5f}".format(Mean_anomaly))
    intermediate_list_4 = intermediate_list_3
    mean_motion_list = str("{:8f}".format(Mean_motion))
    Epoch_rev_list = '447869'
    t_list = line_and_satellite_number_list + inclination_list + intermediate_list + RAAN_list + intermediate_list_2 + eccentricity_list + perigree_list + intermediate_list_3 + mean_anomaly_list + intermediate_list_4 + mean_motion_list + Epoch_rev_list
    """
    position parameters
    """
    a_G0 = 0
    """
    Atmosphere (Aerodynamic)
    """
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
    w_earth = 7.2921150e-5 #rad/s
    """
    Satellite body
    """
    Mass = 20 #kg
    Dimensions = np.array(([0.3, 0.3, 0.4])) # Lx, Ly, Lz
    Ix = 0.4 #kg.m^2
    Iy = 0.45 #kg.m^2
    Iz = 0.3 #kg.m^2
    Iw = 88.1e-6 #kgm^2 Inertia of the RW-06 wheel
    Surface_area_i = Dimensions[0] * Dimensions[1]
    """
    Satellite initial position
    """
    quaternion_initial = Quaternion_functions.euler_to_quaternion(0,0,0) #roll, pitch, yaw
    wbi = np.array(([0.1],[0.05],[0.01]))
    initial_angular_momentum = np.zeros((3,1))
    """
    Max parameters of actuaters
    """
    theta_d_max = 2 #degrees per second (theta derived), angular velocity
    theta_d_d = 0.133 # degrees per second^2 (rotation speed derived), angular acceleration
    h_ws_max = 60e-3 # Nms
    N_ws_max = 5e-3 #Nm
    """
    Control parameters
    """
    w_ref = np.zeros((3,1))
    q_ref = np.array(([0, 0, 1, 0]))        #Quaternion_functions.euler_to_quaternion(0,-45,90) #roll, pitch, yaw
    time = 0
    Ts = 0.1 # Time_step
    Kp = Ts*Ix * 30
    Kd = Ts*Ix * 125

    


