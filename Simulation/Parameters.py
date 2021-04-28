import math
import numpy as np
import Quaternion_functions
from sgp4.api import jday
from struct import *
from scipy import special

pi = math.pi

class SET_PARAMS:
    # All the parameters specific to the satellite and the mission
    """
    Orbit parameters
    """
    eccentricity = 0.000092             # Update eccentricity list
    inclination = 97.4                  # degrees
    Semi_major_axis = 6879.55           # km The distance from the satellite to the earth + the earth radius
    Height_above_earth_surface = 500e3  # distance above earth surface
    Scale_height = 8500                 # scale height of earth atmosphere
    RAAN = 275*pi/180                   # Right ascension of the ascending node in radians
    AP = 0                              # argument of perigee
    Re = 6371.2                         # km magnetic reference radius
    Mean_motion = 15.2355000000         # rev/day
    Mean_anomaly = 29.3                 # degrees
    Argument_of_perigee = 57.4          # in degrees
    omega = Argument_of_perigee
    Period = 86400/Mean_motion          # seconds
    J_t,fr = jday(2020,2,16,15,30,0)    # current julian date
    epoch = J_t - 2433281.5 + fr
    Drag_term = 0.000194                # Remember to update the list term
    wo = (Mean_motion/(3600*24)) / (2*pi*Semi_major_axis * 1000) #rad/s
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
    """
    s_list = '1 25544U 98067A   19343.69339541  .00001764  00000-0  38792-4 0  9991'
    t_list = '2 25544  51.6439 211.2001 0007417  17.6667  85.6398 15.50103472202482'
    """
    """
    position parameters
    """
    a_G0 = 0        # Angle from the greenwhich
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
    k = 10 #order of expansion
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
    quaternion_initial = np.array(([0, 0, 1, 0])) #Quaternion_functions.euler_to_quaternion(0,0,0) #roll, pitch, yaw
    wbi = np.array(([0.0],[0.0],[0.0]))
    initial_angular_wheels = np.zeros((3,1))
    """
    Max parameters of actuaters
    """
    wheel_angular_d_max = 2.0 #degrees per second (theta derived), angular velocity
    wheel_angular_d_d = 0.133 # degrees per second^2 (rotation speed derived), angular acceleration
    h_ws_max = 15.7e-3 # Nms
    N_ws_max = 1.05e-3 #Nm
    M_magnetic_max = 25e-6 #Nm
    """
    Control parameters
    """
    w_ref = np.zeros((3,1)) # desired angular velocity of satellite
    q_ref = np.array(([0, 0, 1, 0])) # initial position of satellite
    time = 1
    Ts = 1 # Time_step
    wn = 5e-6/Ts
    Kp = Ts*Ix * 12.5 #2 * wn**2 
    Kd = Ts*Ix * 125 #2 * wn * 0.707
    Kd_magnet = 1e-7
    Ks_magnet = 1e-7
    """
    Display parameters
    """
    faster_than_control = 1.0 # how much faster the satellite will move around the earth in simulation than the control
    Display = False # if display is desired or not
    skip = 20  # the number of iterations before display
    Number_of_orbits = 1
    Number_of_multiple_orbits = 1
    """
    Visualize measurements
    """
    Visualize = False
    """
    Sensor Parameters
    """
    Magnetometer_noise = 0.0001         #standard deviation of magnetometer noise in Tesla
    Sun_noise = 0.0001                    #standard deviation away from where the actual sun is
    Earth_noise = 0.0001                  #standard deviation away from where the actual earth is
    """
    CSV file parameters
    """
    save_as = ".csv"
    load_as = ".csv"
    """
    Storage of data for prediction
    """
    data_mode = "_buffer"
    buffer_mode = True
    buffer_size = 20

    # File names for the storage of the data attained during the simulation
    filename = "Data_files/Faults" + data_mode
    """
    Mode of operation
    """
    Mode = "Nominal"  
    """
    Fault types and fault parameters
    """
    number_of_faults = 17
    Fault_names = {
    "None": 1,
    "Electronics": 2,
    "Overheated": 3, 
    "Catastrophic_RW": 4,
    "Catastrophic_sun": 5, 
    "Erroneous": 6, 
    "Inverted_polarities": 7,
    "Interference": 8, 
    "Stop": 9, 
    "Closed_shutter": 10,
    "Increasing": 11, 
    "Decrease": 12, 
    "Oscillates": 13,
    "Bit_flip": 14,
    "Sign_flip": 15,
    "Insertion": 16,
    "High_noise": 17,
    }
    likelyhood_multiplier = 1000
    Fault_simulation_mode = 1 # Continued failure, a mistake that does not go back to normal
    Fault_simulation_mode = 0 # Failure is based on specified class failure rate. Multiple failures can occure simultaneously

Min_high_noise = 5
Max_high_noise = 10

Min_high_speed_percentage = 0.9
Max_high_speed_percentage = 1

min_inteference = 3
max_interference = 5

Min_low_speed_percentage = 0.9
Max_low_speed_percentage = 1

def bitflip(x,pos):
    fs = pack('f',x)
    bval = list(unpack('BBBB',fs))
    [q,r] = divmod(pos,8)
    bval[q] ^= 1 << r
    fs = pack('BBBB', *bval)
    fnew=unpack('f',fs)
    return fnew[0]

def random_size(minimum, maximum):
    return np.clip(np.random.normal((minimum+maximum)/2,(maximum-minimum)/2),minimum,maximum)

def random_bit_flip(input_var):
    position = np.random.randint(0, 32)
    input_var = bitflip(input_var, position)
    return input_var

def Reliability(t,n,B):
    return np.exp(-(t / n)**B)

def weibull(t,n,B):
    return (B / n) * (t / n)**(B - 1) * np.exp(-(t / n)**B)

class Fault_parameters:
    def __init__(self, Fault_per_hour, number_of_failures, failures, seed):
        self.np_random = np.random
        self.np_random.seed(seed)
        self.Fault_rate_per_hour = Fault_per_hour
        self.failure = "None"
        self.failures = failures
        self.number_of_failures = number_of_failures
        self.Fault_rate_per_second = self.Fault_rate_per_hour/3600
        self.n = 1/(self.Fault_rate_per_second)
        self.Beta = 0.4287
        self.time = int(SET_PARAMS.Number_of_orbits*SET_PARAMS.Period/(SET_PARAMS.faster_than_control*SET_PARAMS.Ts))
        gamma = special.gammainc(1/self.Beta, (self.time/self.n)**self.Beta)
        self.Reliability_area = self.n * ((self.time/self.n)*gamma)/(self.Beta*((self.time/self.n)**self.Beta)**(1/self.Beta))
        self.Reliability_area_per_time_step = 0

    def Failure_Reliability_area(self, t):
        self.Reliability_area_per_time_step += weibull(t, self.n, self.Beta)*SET_PARAMS.Ts
        Failed = True if self.np_random.uniform(0,1) < self.Reliability_area_per_time_step/self.Reliability_area else False
        if Failed:
            ind = self.np_random.randint(0,self.number_of_failures) 
            self.failure = self.failures[ind]

        return self.failure

    def Failure(self,t):
        mean = weibull(t, self.n, self.Beta)
        self.np_random.normal(mean, 1)

    def random_(self):
        return self.np_random.normal(self.weibull_mean, self.weibull_std) 

class Reaction_wheels(Fault_parameters):
    def __init__(self, seed):
        self.angular_wheels = SET_PARAMS.wheel_angular_d_max
        self.angular_wheels_max = SET_PARAMS.h_ws_max*random_size(minimum = Min_high_speed_percentage, maximum = Max_high_speed_percentage)
        self.Fault_rate_per_hour = 2.5e-7 * SET_PARAMS.likelyhood_multiplier
        self.number_of_failures = 3
        self.failures = {
            0: "Electronics",
            1: "Overheated",
            2: "Catastrophic_RW"
        }
        super().__init__(self.Fault_rate_per_hour, self.number_of_failures, self.failures, seed)


    def Electronics_failure(self, angular_wheels):
        self.angular_wheels = (angular_wheels - angular_wheels/100) if self.failure == "Electronics" else angular_wheels
        return self.angular_wheels

    def Overheated(self, angular_wheels):
        self.angular_wheels = (angular_wheels - angular_wheels/100) if self.failure == "Overheated" else angular_wheels
        return self.angular_wheels

    def Catastrophic(self, angular_wheels):
        self.angular_wheels = 0 if self.failure == "Catastrophic_RW" else angular_wheels
        return self.angular_wheels

class Sun_sensor(Fault_parameters):
    def __init__(self, seed):
        self.Fault_rate_per_hour = 8.15e-9 * SET_PARAMS.likelyhood_multiplier
        self.number_of_failures = 2
        self.failures = {
            0: "Catastrophic_sun",
            1: "Erroneous"
        }
        super().__init__(self.Fault_rate_per_hour, self.number_of_failures, self.failures, seed)

    def Catastrophic(self, sun_sensor):
        return np.zeros(sun_sensor.shape) if self.failure == "Catastrophic_sun" else sun_sensor

    def Erroneous(self, sun_sensor):
        # Sun_sensor must be provided as a unit vector
        return self.np_random.uniform(-1,1,sun_sensor.shape) if self.failure == "Erroneous" else sun_sensor

class Magnetorquers(Fault_parameters):
    def __init__(self, seed):
        self.Fault_rate_per_hour = 8.15e-9 * SET_PARAMS.likelyhood_multiplier
        self.number_of_failures = 2
        self.failures = {
            0: "Inverted_polarities",
            1: "Interference"
        }
        super().__init__(self.Fault_rate_per_hour, self.number_of_failures, self.failures, seed)
    
    def Inverted_polarities(self, magnetic_torquers):
        # Inverted polarities meand that the magnetic torquers will move in the oppositie direction (thus multiplied by -1)
        return -magnetic_torquers if self.failure == "Inverted_polarities" else magnetic_torquers

    def Interference(self, Magnetorquers):
        return Magnetorquers*random_size(min_inteference, max_interference) if self.failure == "Interference" else magnetic_torquers

class Magnetometers(Fault_parameters):
    def __init__(self, seed):
        self.Fault_rate_per_hour = 8.15e-9 * SET_PARAMS.likelyhood_multiplier
        self.number_of_failures = 2
        self.failures = {
            0: "Stop",
            1: "Interference"
        }
        super().__init__(self.Fault_rate_per_hour, self.number_of_failures, self.failures, seed)

    def Stop(self, magnetometer):
        # All of the magnetometers are zero
        return np.zeros(magnetometer.shape) if self.failure == "Stop" else magnetometer

    def Interference(self, magnetometers):
        return magnetometers*random_size(min_inteference, max_interference) if self.failure == "Interference" else magnetometers

class Earth_Sensor(Fault_parameters):
    def __init__(self, seed):
        self.Fault_rate_per_hour = 8.15e-9 * SET_PARAMS.likelyhood_multiplier
        self.number_of_failures = 1
        self.failures = {
            0: "None"
        }
        super().__init__(self.Fault_rate_per_hour, self.number_of_failures, self.failures, seed)

class Star_tracker(Fault_parameters):
    def __init__(self, seed):
        self.Fault_rate_per_hour = 8.15e-9 * SET_PARAMS.likelyhood_multiplier
        self.number_of_failures = 1
        self.failures = {
            0: "Closed_shutter"
        }
        super().__init__(self.Fault_rate_per_hour, self.number_of_failures, self.failures, seed)

    def Closed_shutter(self, Star_tracker):
        return np.zeros(Star_tracker.shape) if self.failure == "Closed_shutter" else Star_tracker

class Overall_control(Fault_parameters):
    def __init__(self, seed):
        self.Fault_rate_per_hour = 8.15e-9 * SET_PARAMS.likelyhood_multiplier
        self.number_of_failures = 3
        self.failures = {
            0: "Increasing",
            1: "Decrease",
            2: "Oscillates"
        }
        self.previous_mul = -1
        self.oscillation_magnitude = 0.1
        self.angular_wheels_max = SET_PARAMS.wheel_angular_d_max*random_size(minimum = Min_high_speed_percentage*0.75, maximum = Max_high_speed_percentage)
        self.angular_wheels_min = SET_PARAMS.wheel_angular_d_max*random_size(minimum = Min_low_speed_percentage, maximum = Max_low_speed_percentage)
        super().__init__(self.Fault_rate_per_hour, self.number_of_failures, self.failures, seed)

    def Increasing(self, angular_wheels):
        self.angular_wheels = np.minimum((angular_wheels - angular_wheels/100), self.angular_wheels_max*np.ones(angular_wheels.shape)) if self.failure == "Increasing" else angular_wheels
        return self.angular_wheels

    def Decrease(self, angular_wheels):
        self.angular_wheels = np.maximum(angular_wheels - angular_wheels/100, self.angular_wheels_min*np.ones(angular_wheels.shape)) if self.failure == "Decreasing" else angular_wheels
        return self.angular_wheels

    def Oscillates(self, angular_wheels):
        angular_wheels = (angular_wheels + angular_wheels*self.oscillation_magnitude*self.previous_mul) if self.failure == "Oscillates" else angular_wheels
        self.previous_mul = self.previous_mul*(-1)
        return angular_wheels


class Common_data_transmission(Fault_parameters):
    def __init__(self, seed):
        self.Fault_rate_per_hour = 8.15e-9 * SET_PARAMS.likelyhood_multiplier
        self.number_of_failures = 3
        self.failures = {
            0: "Bit_flip",
            1: "Sign_flip",
            2: "Insertion"
        }
        super().__init__(self.Fault_rate_per_hour, self.number_of_failures, self.failures, seed)

    def Bit_flip(self, value_to_change):
        if self.failure == "Bit_flip":
            ind = self.np_random.randint(0, len(value_to_change))
            value_to_change[ind] = random_bit_flip(value_to_change[ind])

        return value_to_change

    def Sign_flip(self, value_to_change):
        return -value_to_change if self.failure == "Sign_flip" else value_to_change

    def Insertion(self, value_to_change):
        return 0 if self.failure == "Insertion" else value_to_change

class Sensors_general(Fault_parameters):
    def __init__(self, seed):
        self.Fault_rate_per_hour = 8.15e-9 * SET_PARAMS.likelyhood_multiplier
        self.number_of_failures = 1
        self.failures = {
            0: "High_noise"
        }
        super().__init__(self.Fault_rate_per_hour, self.number_of_failures, self.failures, seed)

    def High_noise(self, sensor):
        return sensor*random_size(minimum = Min_high_noise, maximum = Max_high_noise) if self.failure == "High_noise" else sensor

