import numpy as np
import Controller
from Disturbances import Disturbances
import Parameters
from Parameters import SET_PARAMS, Fault_parameters
from Sensors import Sensors
import matplotlib.pyplot as plt
import Quaternion_functions
import mpld3
import time
import csv
import pandas as pd
from threading import Thread
import concurrent.futures
import multiprocessing
from pathlib import Path
import seaborn as sns
from matplotlib.ticker import EngFormatter
from decimal import Decimal
import math
from plotly.subplots import make_subplots
import plotly.graph_objects as go

Fault_names_to_num = SET_PARAMS.Fault_names

# ! The matplotlib cannot display plots while visual simulation runs.
# ! Consequently the Display and visualize parameters in Parameters 
# ! must be set as desired

if SET_PARAMS.Display:
    import Satellite_display as view

# The DCM must be calculated depending on the current quaternions
def Transformation_matrix(q):
    q1, q2, q3, q4 = q[:]
    A = np.zeros((3,3))
    A[0,0] = q1**2-q2**2-q3**2+q4**2
    A[0,1] = 2*(q1*q2 + q3*q4)
    A[0,2] = 2*(q1*q3 - q2*q4)
    A[1,0] = 2*(q1*q2 - q3*q4)
    A[1,1] = -q1**2+q2**2-q3**2+q4**2
    A[1,2] = 2*(q2*q3 + q1*q4)
    A[2,0] = 2*(q1*q3 + q2*q4)
    A[2,1] = 2*(q2*q3 - q1*q4)
    A[2,2] = -q1**2-q2**2+q3**2+q4**2
    return A

# Function to save a csv file of simulation data
def save_as_excel(Data, sheetnames):
    with pd.ExcelWriter(SET_PARAMS.filename + ".xlsx") as writer:
        i = 0
        for data in Data:
            df = pd.DataFrame(data, columns = data.keys())
            sheetname = sheetnames[i]
            df.to_excel(writer, sheet_name = sheetname, index = False)
            i += 1

####################
# SAVE AS CSV FILE #
####################

def save_as_csv(Data, orbit):
    df = pd.DataFrame(Data, columns = Data.keys())
    df.to_csv(SET_PARAMS.filename + str(orbit) + ".csv")

#######################################################
# FUNCTION TO SAVE A PICKLE FILE OF SIMULATION DATA   #
#######################################################
def save_as_pickle(Data, orbit):
    df = pd.DataFrame(Data, columns = Data.keys())
    df.to_pickle(SET_PARAMS.filename + str(orbit) + ".pkl")

##########################################
# FUNCTION TO VISUALIZE DATA AS GRAPHS   #
##########################################
def visualize_data(D, fault):
    for i in D:
        if i == "Current fault" or i == "Current fault binary" or i == "Current fault numeric":
            pass
        elif i == "Sun in view":
            pass
        else:
            y = np.array((D[i]))
            fig = make_subplots(rows=3, cols=1)
            x = y.shape[0]
            x = np.arange(0,x,1)
            y_min = np.amin(y)
            y_max = np.amax(y)

            fig.append_trace(go.Scatter(
                x=x,
                y=y[:,0],
                name = "x"
            ), row=1, col=1)

            fig.append_trace(go.Scatter(
                x=x,
                y=y[:,1],
                name = "y"
            ), row=2, col=1)

            fig.append_trace(go.Scatter(
                x=x,
                y=y[:,2],
                name = 'z'
            ), row=3, col=1)

            fig.update_yaxes(range=[y_min, y_max], row=1, col=1)
            fig.update_yaxes(range=[y_min, y_max], row=2, col=1)
            fig.update_yaxes(range=[y_min, y_max], row=3, col=1)
            fig.update_layout(height=600, width=600, title_text=str(i))
            fig.write_html("Plots/" + str(fault) +"/"+ str(i)+".html")

##############################################################################
# FUNCTION TO CALCULATE THE ANGULAR MOMENTUM BASED ON THE DERIVATIVE THEREOF #
##############################################################################
def rungeKutta_h(x0, angular, x, h, N_control):
    angular_momentum_derived = N_control
    n = int(np.round((x - x0)/h))

    y = angular
    for _ in range(n):
        k1 = h*(angular_momentum_derived) 
        k2 = h*((angular_momentum_derived) + 0.5*k1) 
        k3 = h*((angular_momentum_derived) + 0.5*k2) 
        k4 = h*((angular_momentum_derived) + 0.5*k3) 

        y = y + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)

        x0 = x0 + h; 
    
    return y

class Dynamics:
    # Initiate initial parameters for the beginning of each orbit set (fault)
    def __init__(self, seed):
        self.seed = seed
        self.np_random = np.random
        self.np_random.seed(seed)                   # Ensures that every fault parameters are implemented with different random seeds
        self.dist = Disturbances()                  # Disturbances of the simulation
        self.w_bi = SET_PARAMS.wbi                  # Angular velocity in ORC
        self.wo = SET_PARAMS.wo                     # Angular velocity of satellite around the earth
        self.angular_wheels = SET_PARAMS.initial_angular_wheels 
        self.q = SET_PARAMS.quaternion_initial      # Quaternion position
        self.t = SET_PARAMS.time                    # Beginning time
        self.dt = SET_PARAMS.Ts                     # Time step
        self.dh = self.dt/10                        # Size of increments for Runga-kutta method
        self.Ix = SET_PARAMS.Ix                     # Ixx inertia
        self.Iy = SET_PARAMS.Iy                     # Iyy inertia
        self.Iz = SET_PARAMS.Iz                     # Izz inertia
        self.Inertia = np.identity(3)*np.array(([self.Ix, self.Iy, self.Iz]))
        self.Iw = SET_PARAMS.Iw                     # Inertia of a reaction wheel
        self.angular_momentum = SET_PARAMS.initial_angular_wheels # Angular momentum of satellite wheels
        self.faster_than_control = SET_PARAMS.faster_than_control   # If it is required that satellite must move faster around the earth than Ts
        self.control = Controller.Control()         # Controller.py is used for control of satellite    
        self.initiate_fault_parameters()
        """
        ! self.Orbit_Data = {
        !    "Sun x": [], "Sun y": [], "Sun z": [],              #S_o measurement (vector of sun in ORC)
        !    "Magnetometer x": [], "Magnetometer y": [], "Magnetometer z": [],     #B vector in SBC
        !    "Earth x": [],"Earth y": [], "Earth z": [],            #Satellite position vector in ORC
        !    "Angular momentum of wheels x": [], "Angular momentum of wheels y": [], "Angular momentum of wheels z": [],      #Wheel angular velocity of each reaction wheel
        !    "Angular momentum of satellite x": [], "Angular momentum of satellite y": [], "Angular momentum of satellite z": [],
        !    "Sun in view": [],                              #True or False values depending on whether the sun is in view of the satellite
        !    "Current fault": [],                            #What the fault is that the system is currently experiencing
        !    "Current fault numeric": [],
        !    "Current fault binary": []
        !}
        """
        ####################################################
        #  THE ORBIT_DATA DICTIONARY IS USED TO STORE ALL  #
        #     THE MEASUREMENTS FOR EACH TIMESTEP (TS)      #
        # EACH ORBIT HAS AN INDUCED FAULT WITHIN THE ADCS. #
        ####################################################

        self.Orbit_Data = {
            "Sun": [],            #S_o measurement (vector of sun in ORC)
            "Magnetometer": [],    #B vector in SBC
            "Earth": [],           #Satellite position vector in ORC
            "Angular momentum of wheels": [],    #Wheel angular velocity of each reaction wheel
            "Angular velocity of satellite": [],
            "Sun in view": [],                              #True or False values depending on whether the sun is in view of the satellite
            "Current fault": [],                            #What the fault is that the system is currently experiencing
            "Current fault numeric": [],
            "Current fault binary": []
        }

        self.zeros = np.zeros((SET_PARAMS.number_of_faults,), dtype = int)

        self.fault = "None"                      # Current fault of the system

    def determine_earth_vision(self):
        #################################################################
        #      FOR THIS SPECIFIC SATELLITE MODEL, THE EARTH SENSOR      #
        #                    IS FIXED TO THE -Z FACE                    #
        # THIS IS ACCORDING TO THE ORBIT AS DEFINED BY JANSE VAN VUUREN #
        #             THIS IS DETERMINED WITH THE SBC FRAME             #
        #################################################################

        angle_difference = Quaternion_functions.rad2deg(np.arccos(np.clip(np.dot(self.r_sat_sbc, SET_PARAMS.Earth_sensor_position),-1,1)))
        if angle_difference < SET_PARAMS.Earth_sensor_angle:
            self.r_sat = self.Earth_sensor_fault.normal_noise(self.r_sat, SET_PARAMS.Earth_noise)
            self.r_sat = self.Earth_sensor_fault.General_sensor_high_noise(self.r_sat)
            self.r_sat = self.Common_data_transmission_fault.Bit_flip(self.r_sat)
            self.r_sat = self.Common_data_transmission_fault.Sign_flip(self.r_sat)
            self.r_sat = self.Common_data_transmission_fault.Insertion_of_zero_bit(self.r_sat) 
        else:
            self.r_sat_sbc = np.zeros(self.r_sat_sbc.shape)

    def determine_sun_vision(self):
        #################################################################
        #    FOR THIS SPECIFIC SATELLITE MODEL, THE FINE SUN SENSOR     #
        #       IS FIXED TO THE +X FACE AND THE COARSE SUN SENSOR       #
        #                   IS FIXED TO THE -X FACE.                    #
        # THIS IS ACCORDING TO THE ORBIT AS DEFINED BY JANSE VAN VUUREN #
        #             THIS IS DETERMINED WITH THE SBC FRAME             #
        #################################################################

        if self.sun_in_view:
            angle_difference_fine = Quaternion_functions.rad2deg(np.arccos(np.dot(self.S_b[:,0], SET_PARAMS.Fine_sun_sensor_position)))
            angle_difference_coarse = Quaternion_functions.rad2deg(np.arccos(np.dot(self.S_b[:,0], SET_PARAMS.Coarse_sun_sensor_position)))
            if angle_difference_fine < SET_PARAMS.Fine_sun_sensor_angle:
                self.S_b = self.Sun_sensor_fault.normal_noise(self.S_b, SET_PARAMS.Fine_sun_noise)

                ######################################################
                # IMPLEMENT ERROR OR FAILURE OF SENSOR IF APPLICABLE #
                ######################################################

                self.S_b = self.Sun_sensor_fault.Catastrophic_sun(self.S_b, "Fine")
                self.S_b = self.Sun_sensor_fault.Erroneous_sun(self.S_b, "Fine")
                self.S_b = self.Common_data_transmission_fault.Bit_flip(self.S_b)
                self.S_b = self.Common_data_transmission_fault.Sign_flip(self.S_b)
                self.S_b = self.Common_data_transmission_fault.Insertion_of_zero_bit(self.S_b)  
            

            elif angle_difference_coarse < SET_PARAMS.Coarse_sun_sensor_angle:
                self.S_b = self.Sun_sensor_fault.normal_noise(self.S_b, SET_PARAMS.Coarse_sun_noise)

                ######################################################
                # IMPLEMENT ERROR OR FAILURE OF SENSOR IF APPLICABLE #
                ######################################################

                self.S_b = self.Sun_sensor_fault.Catastrophic_sun(self.S_b, "Coarse")
                self.S_b = self.Sun_sensor_fault.Erroneous_sun(self.S_b, "Coarse")
                self.S_b = self.Common_data_transmission_fault.Bit_flip(self.S_b)
                self.S_b = self.Common_data_transmission_fault.Sign_flip(self.S_b)
                self.S_b = self.Common_data_transmission_fault.Insertion_of_zero_bit(self.S_b)  
            else:
                self.S_b = np.zeros(self.S_b.shape)
        

    def initiate_fault_parameters(self):
        #################################
        # ALL THE CURRENT FAULT CLASSES #
        #################################
        
        self.Reaction_wheel_fault = Parameters.Reaction_wheels(self.seed)
        self.Earth_sensor_fault = Parameters.Earth_Sensor(self.seed)    
        self.Sun_sensor_fault = Parameters.Sun_sensor(self.seed)
        self.Magnetometer_fault = Parameters.Magnetometers(self.seed)
        self.Magnetorquers_fault = Parameters.Magnetorquers(self.seed)
        self.Control_fault = Parameters.Overall_control(self.seed)
        self.Common_data_transmission_fault = Parameters.Common_data_transmission(self.seed)
    
    def initiate_purposed_fault(self, fault):
        self.fault = fault
        self.Reaction_wheel_fault.failure = self.fault
        self.Earth_sensor_fault.failure = self.fault
        self.Magnetometer_fault.failure = self.fault
        self.Sun_sensor_fault.failure = self.fault
        self.Magnetorquers_fault.failure = self.fault
        self.Control_fault.failure = self.fault
        self.Common_data_transmission_fault.failure = self.fault
        # ! self.Star_tracker_fault = Parameters.Star_tracker(self.seed)

    ########################################################################################
    # FUNCTION TO CALCULATE THE SATELLITE ANGULAR VELOCITY BASED ON THE DERIVATIVE THEREOF #
    ########################################################################################
    def rungeKutta_w(self, x0, w, x, h):      
        ######################################################
        # CONTROL TORQUES IMPLEMENTED DUE TO THE CONTROL LAW #
        ######################################################

        N_control_magnetic, N_control_wheel = self.control.control(w, self.q, self.Inertia, self.B, self.Control_fault)

        if "RW" in self.fault:
            N_control_wheel = self.Reaction_wheel_fault.Electronics_of_RW_failure(N_control_wheel)
            N_control_wheel = self.Reaction_wheel_fault.Overheated_RW(N_control_wheel)
            N_control_wheel = self.Reaction_wheel_fault.Catastrophic_RW(N_control_wheel)
            N_control_wheel = self.Control_fault.Increasing_angular_RW_momentum(N_control_wheel)
            N_control_wheel = self.Control_fault.Decreasing_angular_RW_momentum(N_control_wheel)
            N_control_wheel = self.Control_fault.Oscillating_angular_RW_momentum(N_control_wheel)

        self.angular_momentum = np.clip(rungeKutta_h(x0, self.angular_momentum, x, h, N_control_wheel), -SET_PARAMS.h_ws_max, SET_PARAMS.h_ws_max)

        N_aero = 0 # ! self.dist.Aerodynamic(self.A, self.A_EIC_to_ORC)

        ###################################
        # DISTURBANCE OF GRAVITY GRADIENT #
        ###################################

        Ngg = self.dist.Gravity_gradient_func(self.A_ORC_to_SBC) 

        n = int(np.round((x - x0)/h))
        y = w

        for _ in range(n):
            N_gyro = y * (np.matmul(self.Inertia,y) + self.angular_momentum)

            #############################################
            # DISTURBANCE OF A REACTION WHEEL IMBALANCE #
            #############################################

            N_rw = np.reshape(self.dist.Wheel_Imbalance(self.angular_momentum/self.Iw, h),(3,1))   

            ######################################################
            # ALL THE DISTURBANCE TORQUES ADDED TO THE SATELLITE #
            ######################################################

            N_disturbance = Ngg + N_aero + N_rw - N_gyro                
            N_control = N_control_magnetic - N_control_wheel
            N = N_control + N_disturbance
            k1 = h*((np.matmul(np.linalg.inv(self.Inertia), N))) 
            k2 = h*((np.matmul(np.linalg.inv(self.Inertia), N)) + 0.5*k1) 
            k3 = h*((np.matmul(np.linalg.inv(self.Inertia), N)) + 0.5*k2) 
            k4 = h*((np.matmul(np.linalg.inv(self.Inertia), N)) + 0.5*k3) 
            y = y + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)
    
            x0 = x0 + h; 

        y = np.clip(y, -SET_PARAMS.wheel_angular_d_max, SET_PARAMS.wheel_angular_d_max)

        return y

    ###########################################################################################
    # FUNCTION TO CALCULATE THE SATELLITE QUATERNION POSITION BASED ON THE DERIVATIVE THEREOF #
    ###########################################################################################
    def rungeKutta_q(self, x0, y0, x, h):      
        wx, wy, wz = self.w_bo[:,0]
        n = int(np.round((x - x0)/h))

        y = y0

        for _ in range(n):
            k1 = h*(0.5 * np.matmul(np.array(([0, wz, -wy, wx], [-wz, 0, wx, wy], [wy, -wx, 0, wz], [-wx, -wy, -wz, 0])),y))
            k2 = h*(0.5 * np.matmul(np.array(([0, wz, -wy, wx], [-wz, 0, wx, wy], [wy, -wx, 0, wz], [-wx, -wy, -wz, 0])),y + 0.5*k1))
            k3 = h*(0.5 * np.matmul(np.array(([0, wz, -wy, wx], [-wz, 0, wx, wy], [wy, -wx, 0, wz], [-wx, -wy, -wz, 0])),y + 0.5*k2))
            k4 = h*(0.5 * np.matmul(np.array(([0, wz, -wy, wx], [-wz, 0, wx, wy], [wy, -wx, 0, wz], [-wx, -wy, -wz, 0])),y + 0.5*k3))

            y = y + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)
    
            x0 = x0 + h; 
        
        y = y/np.linalg.norm(y)

        return y

    def Fault_implementation(self):
        if self.fault == "None":
            Faults = []
            True_faults = []
            Faults.append(self.Reaction_wheel_fault.Failure_Reliability_area(self.t)) 
            Faults.append(self.Earth_sensor_fault.Failure_Reliability_area(self.t))   
            Faults.append(self.Sun_sensor_fault.Failure_Reliability_area(self.t))
            Faults.append(self.Magnetometer_fault.Failure_Reliability_area(self.t))
            Faults.append(self.Magnetorquers_fault.Failure_Reliability_area(self.t))
            Faults.append(self.Control_fault.Failure_Reliability_area(self.t))
            Faults.append(self.Common_data_transmission_fault.Failure_Reliability_area(self.t))
            # ! Faults.append(self.Star_tracker_fault.Failure_Reliability_area(self.t))
            for fault in Faults:
                if fault != "None":
                    True_faults.append(fault)
                
            if True_faults:
                self.fault = True_faults[self.np_random.randint(0,len(True_faults))]
                print(self.fault)
                

    def rotation(self):
        ##############################################################
        #    DETERMINE WHETHER A FAUKLT OCCURED WITHIN THE SYSTEM    #
        # BASED ON THE STATISTICAL PROBABILITY DEFINED IN PARAMETERS #
        ##############################################################

        self.Fault_implementation()

        ########################################
        # DETERMINE THE ACTUAL POSITION OF THE #
        # SATELLITE FROM THE EARTH AND THE SUN #
        ########################################

        self.r_sat, self.v_sat_EIC, self.A_EIC_to_ORC, r_EIC = sense.satellite_vector(self.t*self.faster_than_control)
        self.r_EIC = self.Earth_sensor_fault.General_sensor_high_noise(np.asarray(r_EIC).astype(np.float64))

        self.S_EIC, self.sun_in_view = sense.sun(self.t*self.faster_than_control)

        ######################################
        # DETERMINE THE DCM OF THE SATELLITE #
        ######################################

        self.A_ORC_to_SBC = Transformation_matrix(self.q)
        self.A_EIC_to_SBC = np.matmul(self.A_EIC_to_ORC, self.A_ORC_to_SBC)
        
        self.r_sat_sbc = np.matmul(self.A_ORC_to_SBC, self.r_sat)
        self.r_sat_sbc = self.r_sat_sbc/np.linalg.norm(self.r_sat_sbc)
        
        self.S_b = np.matmul(self.A_EIC_to_SBC, self.S_EIC)
        
        if self.sun_in_view:
            self.S_b = self.S_b/np.linalg.norm(self.S_b)

        ##################################################
        # DETERMINE WHETHER THE SUN AND THE EARTH SENSOR #
        #   IS IN VIEW OF THE VECTOR ON THE SATELLITE    #
        ##################################################

        self.determine_sun_vision()
        self.determine_earth_vision()
    
        self.Beta = sense.magnetometer(self.t) 

        self.B = np.matmul(self.A_EIC_to_SBC,self.Beta)

        ######################################################
        # IMPLEMENT ERROR OR FAILURE OF SENSOR IF APPLICABLE #
        ######################################################

        self.B = self.Magnetometer_fault.Stop_magnetometers (self.B)
        self.B = self.Magnetometer_fault.Interference_magnetic(self.B)
        self.B = self.Magnetometer_fault.General_sensor_high_noise(self.B)
        self.B = self.Common_data_transmission_fault.Bit_flip(self.B)
        self.B = self.Common_data_transmission_fault.Sign_flip(self.B)
        self.B = self.Common_data_transmission_fault.Insertion_of_zero_bit(self.B)
        
        ########################################################
        # THE ERROR FOR W_BI IS WITHIN THE RUNGEKUTTA FUNCTION #
        ########################################################

        self.w_bi = self.rungeKutta_w(self.t, self.w_bi, self.t+self.dt, self.dh)
        
        self.w_bo = self.w_bi - np.matmul(self.A_ORC_to_SBC, np.array(([0],[self.wo],[0])))
        
        self.q = self.rungeKutta_q(self.t, self.q, self.t+self.dt, self.dh)
        
        self.t += self.dt
        
        self.update()
        
        return self.w_bi, self.q, self.A_ORC_to_SBC, self.r_EIC, self.sun_in_view

    def update(self):
        """
        ! self.Orbit_Data["Magnetometer x"].append(self.B[0])
        ! self.Orbit_Data["Magnetometer z"].append(self.B[1])
        ! self.Orbit_Data["Magnetometer y"].append(self.B[2])
        ! self.Orbit_Data["Sun x"].append(self.S_b[0][0])
        ! self.Orbit_Data["Sun y"].append(self.S_b[1][0])
        ! self.Orbit_Data["Sun z"].append(self.S_b[2][0])
        ! self.Orbit_Data["Earth x"].append(self.r_sat_sbc[0])
        ! self.Orbit_Data["Earth y"].append(self.r_sat_sbc[1])
        ! self.Orbit_Data["Earth z"].append(self.r_sat_sbc[2])
        ! self.Orbit_Data["Angular momentum of wheels x"].append(self.angular_momentum[0][0])
        ! self.Orbit_Data["Angular momentum of wheels y"].append(self.angular_momentum[1][0])
        ! self.Orbit_Data["Angular momentum of wheels z"].append(self.angular_momentum[2][0])
        ! self.Orbit_Data["Angular momentum of satellite x"].append(self.w_bi[0][0])
        ! self.Orbit_Data["Angular momentum of satellite y"].append(self.w_bi[1][0])
        ! self.Orbit_Data["Angular momentum of satellite z"].append(self.w_bi[2][0])
        """
        self.Orbit_Data["Magnetometer"].append(self.B)
        self.Orbit_Data["Sun"].append(self.S_b[:,0])
        self.Orbit_Data["Earth"].append(self.r_sat_sbc)
        self.Orbit_Data["Angular momentum of wheels"].append(self.angular_momentum[:,0])
        self.Orbit_Data["Angular velocity of satellite"].append(self.w_bi[:,0])
        self.Orbit_Data["Sun in view"].append(self.sun_in_view)
        if self.sun_in_view == False and (self.fault == "Catastrophic_sun" or self.fault == "Erroneous"):
            self.Orbit_Data["Current fault"].append("None")
            temp = list(self.zeros)
            temp[Fault_names_to_num["None"] - 1] = 1
            self.Orbit_Data["Current fault numeric"].append(temp)
            self.Orbit_Data["Current fault binary"].append(0)
        else:
            self.Orbit_Data["Current fault"].append(self.fault)
            temp = list(self.zeros)
            temp[Fault_names_to_num[self.fault] - 1] = 1
            self.Orbit_Data["Current fault numeric"].append(temp)
            self.Orbit_Data["Current fault binary"].append(0 if self.fault == "None" else 1)

def loop(index, Data, orbit_descriptions):
    print(SET_PARAMS.Fault_names_values[index])

    if SET_PARAMS.Display:
        satellite = view.initializeCube(SET_PARAMS.Dimensions)
        pv = view.ProjectionViewer(1920, 1080, satellite)
    
    for j in range(1, int(SET_PARAMS.Number_of_orbits*SET_PARAMS.Period/(SET_PARAMS.faster_than_control*SET_PARAMS.Ts)+1)):
        w, q, A, r, sun_in_view = D.rotation()
        if SET_PARAMS.Display and j%SET_PARAMS.skip == 0:
            pv.run(w, q, A, r, sun_in_view)
        
        if j%(int(SET_PARAMS.Number_of_orbits*SET_PARAMS.Period/(SET_PARAMS.faster_than_control*SET_PARAMS.Ts)/10)) == 0:
            print("Number of time steps for orbit loop number", index, " = ", "%.2f" % float(j/int(SET_PARAMS.Number_of_orbits*SET_PARAMS.Period/(SET_PARAMS.faster_than_control*SET_PARAMS.Ts))))

        if SET_PARAMS.Fault_simulation_mode == 2 and j%(int(SET_PARAMS.Number_of_orbits*SET_PARAMS.Period/(SET_PARAMS.faster_than_control*SET_PARAMS.Ts)/SET_PARAMS.fixed_orbit_failure)) == 0:
            D.initiate_purposed_fault(SET_PARAMS.Fault_names_values[index])
            if SET_PARAMS.Display:
                pv.fault = D.fault

    if SET_PARAMS.Visualize and SET_PARAMS.Display == False:
        path_to_folder = Path("Plots/" + str(D.fault))
        path_to_folder.mkdir(exist_ok=True)
        visualize_data(D.Orbit_Data, D.fault)
    
    elif SET_PARAMS.Display == True:
        pv.save_plot(D.fault)

    print("Number of multiple orbits", index)  
    Data[index] = D.Orbit_Data
    orbit_descriptions[index] = D.fault

    if SET_PARAMS.save_as == ".csv":
        save_as_csv(D.Orbit_Data, index)
    else:
        save_as_pickle(D.Orbit_Data, index)

################################################################
# FOR ALL OF THE FAULTS RUN A NUMBER OF ORBITS TO COLLECT DATA #
################################################################
if __name__ == "__main__":
    sense = Sensors()

    #########################################################
    # IF THE SAVE AS IS EQUAL TO XLSX, THE THREADING CANNOT #
    #           BE USED TO SAVE CSV FILES                   #     
    #########################################################
    if SET_PARAMS.save_as == ".xlsx":
        Data = []
        orbit_descriptions = []
        for i in range(4,SET_PARAMS.Number_of_multiple_orbits):
            D = Dynamics(i)

            print(SET_PARAMS.Fault_names_values[i+1])

            if SET_PARAMS.Display:
                satellite = view.initializeCube(SET_PARAMS.Dimensions)
                pv = view.ProjectionViewer(1920, 1080, satellite)
            
            for j in range(int(SET_PARAMS.Number_of_orbits*SET_PARAMS.Period/(SET_PARAMS.faster_than_control*SET_PARAMS.Ts)+1)):
                w, q, A, r, sun_in_view = D.rotation()
                if SET_PARAMS.Display and j%SET_PARAMS.skip == 0:
                    pv.run(w, q, A, r, sun_in_view)
                
                if j%(int(SET_PARAMS.Number_of_orbits*SET_PARAMS.Period/(SET_PARAMS.faster_than_control*SET_PARAMS.Ts)/10)) == 0:
                    print("Number of time steps for orbit loop number", i, " = ", "%.2f" % float(j/int(SET_PARAMS.Number_of_orbits*SET_PARAMS.Period/(SET_PARAMS.faster_than_control*SET_PARAMS.Ts))))

                if SET_PARAMS.Fault_simulation_mode == 2 and (j+1)%(int(SET_PARAMS.Number_of_orbits*SET_PARAMS.Period/(SET_PARAMS.faster_than_control*SET_PARAMS.Ts)/SET_PARAMS.fixed_orbit_failure)) == 0:
                    D.initiate_purposed_fault(SET_PARAMS.Fault_names_values[i+1])
                    if SET_PARAMS.Display:
                        pv.fault = D.fault

            if SET_PARAMS.Visualize and SET_PARAMS.Display == False:
                path_to_folder = Path("Plots/" + str(D.fault))
                path_to_folder.mkdir(exist_ok=True)
                visualize_data(D.Orbit_Data, D.fault)
            
            elif SET_PARAMS.Display == True:
                pv.save_plot(D.fault)

            Data.append(D.Orbit_Data)
            orbit_descriptions.append(str(D.fault))

        save_as_excel(Data, orbit_descriptions)

    ######################################################
    # IF THE SAVE AS IS NOT EQUAL TO XLSX, THE THREADING #
    #           CAN BE USED TO SAVE CSV FILES            #
    ######################################################
    else:
        threads = []

        manager = multiprocessing.Manager()
        Data = manager.dict()
        orbit_descriptions = manager.dict()

        for i in range(1, SET_PARAMS.Number_of_multiple_orbits+1):
            D = Dynamics(i)

            t = multiprocessing.Process(target=loop, args=(i,Data, orbit_descriptions))
            threads.append(t)
            t.start()
            print("Beginning of", i)
            if i%15 == 0 or i == SET_PARAMS.Number_of_multiple_orbits:
                for process in threads:     
                    process.join()

                threads = []





