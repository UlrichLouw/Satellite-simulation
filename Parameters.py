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
    Mean_motion = 15.2355000000 #rev/day
    Mean_anomaly = 29.3 #degrees
    Argument_of_perigee = 57.4 #in degrees
    omega = Argument_of_perigee
    Period = 86400/Mean_motion #seconds
    J_t,fr = jday(2020,2,16,15,30,0) #current julian date
    epoch = J_t - 2433281.5 + fr
    Drag_term = 0.000194 # Remember to update the list term
    wo = (Mean_motion/(3600*24)) / (2*pi*Semi_major_axis * 1000) #rad/s
    """
    Sun synchronous parameters
    """
    """
    eccentricity = 0.0001313349 # Update eccentricity list
    inclination = 98.000640878#degrees
    Semi_major_axis = 7008.4715082 #km The distance from the satellite to the earth + the earth radius
    Height_above_earth_surface = 700e3 #distance above earth surface
    Scale_height = 8500 #scale height of earth atmosphere
    RAAN = 71.497449685*pi/180 #Right ascension of the ascending node in radians
    AP = 69.532933059 #argument of perigee
    Re = 6371.2 #km magnetic reference radius
    Mean_motion = 15.2355000000 #rev/day
    Mean_anomaly = 290.42758119 #degrees
    Argument_of_perigee = 69.532933059 #in degrees
    omega = Argument_of_perigee
    Period = 5834.805 #seconds
    J_t,fr = jday(2020,2,16,15,30,0) #current julian date
    epoch = J_t - 2433281.5 + fr
    Drag_term = 0.000194 # Remember to update the list term
    wo = (Mean_motion/(3600*24)) / (2*pi*Semi_major_axis * 1000) #rad/s
    """
    """
    TLE data
    """
    # s list
    satellite_number_list = '1 25544U'
    international_list = ' 98067A   '
    epoch_list = str("{:.8f}".format(epoch))
    mean_motion_derivative_first_list = '  .00001764'
    mean_motion_derivative_second_list = '  00000-0'
    Drag_term_list = '  19400-4' # B-star
    Ephereris_list = ' 0'
    element_num_checksum_list = '  7030'
    s_list = satellite_number_list + international_list + epoch_list + mean_motion_derivative_first_list + mean_motion_derivative_second_list + Drag_term_list + Ephereris_list + element_num_checksum_list
    # t list
    line_and_satellite_number_list = '2 27843  '
    inclination_list = str("{:.4f}".format(inclination))
    intermediate_list = ' '
    RAAN_list = str("{:.4f}".format(RAAN*180/pi))
    intermediate_list_2 = ' '
    eccentricity_list = '0000920  '
    perigree_list = str("{:.4f}".format(Argument_of_perigee))
    intermediate_list_3 = intermediate_list_2 + ' '
    mean_anomaly_list = str("{:.4f}".format(Mean_anomaly))
    intermediate_list_4 = intermediate_list_2
    mean_motion_list = str("{:8f}".format(Mean_motion)) + '00'
    Epoch_rev_list = '000009'
    t_list = line_and_satellite_number_list + inclination_list + intermediate_list + RAAN_list + intermediate_list_2 + eccentricity_list + perigree_list + intermediate_list_3 + mean_anomaly_list + intermediate_list_4 + mean_motion_list + Epoch_rev_list
    """
    Overwrite Jansen vuuren se waardes met sgp4 example
    """
    #s_list = '1 25544U 98067A   19343.69339541  .00001764  00000-0  38792-4 0  9991'
    #t_list = '2 25544  51.6439 211.2001 0007417  17.6667  85.6398 15.50103472202482'
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
    h = -8.0 #
    m = 10
    n = 10
    Radius_earth = 6371e3 # in m
    w_earth = 7.2921150e-5 #rad/s
    """
    Sun parameters
    """
    Radius_sun = 696340e3 # in m
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
    M_magnetic_max = 1 #Am^2
    """
    Control parameters
    """
    w_ref = np.zeros((3,1)) # desired angular velocity of satellite
    q_ref = np.array(([0, 0, 1, 0])) # initial position of satellite
    time = 0
    Ts = 1 # Time_step
    Kp = Ts*Ix * 30
    Kd = Ts*Ix * 125
    Kd_magnet = 1
    """
    Display parameters
    """
    faster_than_control = 1.0 # how much faster the satellite will move around the earth in simulation than the control
    Display = False  # if display is desired or not
    skip = 20   # the number of iterations before display
    """
    Sensor Parameters
    """
    magnetometer_noise = 0.5e-6 # standard deviation of magnetometer noise in Tesla
    sun_noise = 9 # degrees away from where the actual sun is
    """
    CSV file parameters
    """
    Save_file = False
    """
    Fault types
    """
    Fault_names = {
        1: "Reaction_wheel",
        2: "Sun sensor",
        3: "Magnetometer",
        4: "Earth sensor",
    }
    


