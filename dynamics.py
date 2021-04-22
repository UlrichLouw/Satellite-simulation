import numpy as np
import Controller
from Disturbances import Disturbances
import Parameters
from Parameters import SET_PARAMS
from Sensors import Sensors
import matplotlib.pyplot as plt
import time
import csv
import pandas as pd
from threading import Thread
import concurrent.futures
import multiprocessing

Fault_names_to_num = SET_PARAMS.Fault_names
# The matplotlib cannot display plots while visual simulation runs.
# Consequently the Display and visualize parameters in Parameters 
# must be set as desired
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
    with pd.ExcelWriter(SET_PARAMS.excel_filename) as writer:
        i = 0
        for data in Data:
            df = pd.DataFrame(data, columns = data.keys())
            sheetname = sheetnames[i]
            df.to_excel(writer, sheet_name = sheetname, index = False)
            i += 1

# Save as csv file
def save_as_csv(Data, sheetnames):
    i = 0
    for data in Data:
        df = pd.DataFrame(data, columns = data.keys())
        df.to_csv(SET_PARAMS.csv_filename + str(i) + ".csv")
        i += 1

# Function to save a pickle file of simulation data
def save_as_pickle(Data):
    i = 0
    for data in Data:
        df = pd.DataFrame(data, columns = data.keys())
        df.to_pickle(SET_PARAMS.pickle_filename + str(i) + ".pkl")
        i += 1

# Function to visualize data as graphs
def visualize_data(D):
    for i in D.Orbit_Data:
        if i == "Current fault":
            pass
        elif i == "Sun in view":
            fig = plt.figure()
            y = np.array((D.Orbit_Data[i]))
            plt.plot(y)
            plt.title(i)
            plt.show() #plt.pause(0.001)
        else:
            fig = plt.figure()
            y = np.array((D.Orbit_Data[i]))
            plt.plot(y)
            plt.title(i)
            plt.show() #plt.pause(1)

# Function to calculate the angular momentum based on the derivative thereof
def rungeKutta_h(self, x0, angular, x, h, N_control):
    angular_momentum_derived = N_control
    n = int(np.round((x - x0)/h))

    y = angular
    for i in range(n):
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
        self.np_random.seed(seed)                            # Ensures that every fault parameters are implemented with different random seeds
        self.dist = Disturbances()                  # Disturbances of the simulation
        self.w_bi = SET_PARAMS.wbi                  # Angular velocity in ORC
        self.wo = SET_PARAMS.wo                     # Angular velocity of satellite around the earth
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

        self.Orbit_Data = {
            "Sun x": [], "Sun y": [], "Sun z": [],              #S_o measurement (vector of sun in ORC)
            "Magnetometer x": [], "Magnetometer y": [], "Magnetometer z": [],     #B vector in SBC
            "Earth x": [],"Earth y": [], "Earth z": [],            #Satellite position vector in ORC
            "Angular momentum of wheels x": [],
            "Angular momentum of wheels y": [], 
            "Angular momentum of wheels z": [],      #Wheel angular velocity of each reaction wheel
            "Sun in view": [],                              #True or False values depending on whether the sun is in view of the satellite
            "Current fault": [],                            #What the fault is that the system is currently experiencing
            "Current fault numeric": [],
            "Current fault binary": []
        }

        self.zeros = np.zeros((17,), dtype = int)

        self.fault = "None"                      # Current fault of the system
        # The Orbit_Data dictionary is used to store all the measurements for each timestep (Ts)
        # Each Orbit has an induced fault within the ADCS. 

    def initiate_fault_parameters(self):
        """
        All the current fault classes
        """
        self.Reaction_wheel_fault = Parameters.Reaction_wheels(self.seed)
        self.Earth_sensor_fault = Parameters.Earth_Sensor(self.seed)    
        self.Sun_sensor_fault = Parameters.Sun_sensor(self.seed)
        self.Magnetometer_fault = Parameters.Magnetometers(self.seed)
        self.Magnetorquers_fault = Parameters.Magnetorquers(self.seed)
        self.Control_fault = Parameters.Overall_control(self.seed)
        self.Common_data_transmission_fault = Parameters.Common_data_transmission(self.seed)
        self.Sensors_general_fault = Parameters.Sensors_general(self.seed)
        #self.Star_tracker_fault = Parameters.Star_tracker(self.seed)

    # Function to calculate the satellite angular velocity based on the derivative thereof
    def rungeKutta_w(self, x0, w, x, h, r_sat):
        if self.fault == "None":
            N_aero = 0 #self.dist.Aerodynamic(self.A, self.A_EIC_to_ORC)
            N_rw = np.reshape(self.dist.Wheel_Imbalance(w, x),(3,1))    # Disturbance of a reaction wheel imbalance
            Ngg = self.dist.Gravity_gradient_func(self.A)               # Disturbance of gravity gradient
            
            # Control torques implemented due to the control law
            N_control_magnetic, N_control_wheel = self.control.control(w, self.q, self.Inertia, self.B, self.Control_fault)

            self.angular_momentum = np.clip(rungeKutta_h(self, x0, self.angular_momentum, x, h, N_control_wheel), -SET_PARAMS.h_ws_max, SET_PARAMS.h_ws_max)
            
            N_gyro = -w * (np.matmul(self.Inertia,w) + self.angular_momentum)

            N_disturbance = Ngg + N_aero + N_rw + N_gyro                # All then disturbance torques added to the satellite

            n = int(np.round((x - x0)/h))
            y = w

            for i in range(n):
                k1 = h*((N_gyro + (-N_control_magnetic  + N_disturbance))) 
                k2 = h*((N_gyro + (-N_control_magnetic  + N_disturbance)) + 0.5*k1) 
                k3 = h*((N_gyro + (-N_control_magnetic  + N_disturbance)) + 0.5*k2) 
                k4 = h*((N_gyro + (-N_control_magnetic  + N_disturbance)) + 0.5*k3) 
                y = y + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)
        
                x0 = x0 + h; 

            np.clip(y, -SET_PARAMS.wheel_angular_d_max, SET_PARAMS.wheel_angular_d_max)
        else:
            y = self.Reaction_wheel_fault.Electronics_failure(w)
            y = self.Reaction_wheel_fault.Overheated(y)
            y = self.Reaction_wheel_fault.Catastrophic(y)
            y = self.Control_fault.Increasing(y)
            y = self.Control_fault.Decrease(y)
            y = self.Control_fault.Oscillates(y)

        return y

    # Function to calculate the satellite quaternion position based on the derivative thereof
    def rungeKutta_q(self, x0, y0, x, h):
        
        wx, wy, wz = self.w_bo[:,0]
        n = int(np.round((x - x0)/h))

        y = y0

        for i in range(n):
            k1 = h*(0.5 * np.matmul(np.array(([0, wz, -wy, wx], [-wz, 0, wx, wy], [wy, -wx, 0, wz], [-wx, -wy, -wz, 0])),y))
            k2 = h*(0.5 * np.matmul(np.array(([0, wz, -wy, wx], [-wz, 0, wx, wy], [wy, -wx, 0, wz], [-wx, -wy, -wz, 0])),y + 0.5*k1))
            k3 = h*(0.5 * np.matmul(np.array(([0, wz, -wy, wx], [-wz, 0, wx, wy], [wy, -wx, 0, wz], [-wx, -wy, -wz, 0])),y + 0.5*k2))
            k4 = h*(0.5 * np.matmul(np.array(([0, wz, -wy, wx], [-wz, 0, wx, wy], [wy, -wx, 0, wz], [-wx, -wy, -wz, 0])),y + 0.5*k3))

            y = y + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)
    
            x0 = x0 + h; 
        
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
            Faults.append(self.Sensors_general_fault.Failure_Reliability_area(self.t))
            #Faults.append(self.Star_tracker_fault.Failure_Reliability_area(self.t))
            for fault in Faults:
                if fault != "None":
                    True_faults.append(fault)
                
            if True_faults:
                self.fault = True_faults[self.np_random.randint(0,len(True_faults))]
                print(self.fault)

    def rotation(self):
        if self.t == SET_PARAMS.time:
            self.initiate_fault_parameters()

        self.Fault_implementation()

        self.r_sat, v_sat, self.A_EIC_to_ORC, r_EIC = sense.satellite_vector(self.t*self.faster_than_control, error=self.Earth_sensor_fault)

        self.S_EIC, self.sun_in_view = sense.sun(self.t*self.faster_than_control)       

        self.A = np.matmul(self.A_EIC_to_ORC, Transformation_matrix(self.q))
        
        self.r_sat_sbc = np.matmul(self.A, self.r_sat)
        
        self.S_o = np.matmul(self.A_EIC_to_ORC, self.S_EIC)
        
        self.S_b = np.matmul(self.A, self.S_o)
        
        self.Beta = sense.magnetometer(self.t) 

        self.B = np.matmul(self.A,np.matmul(self.A_EIC_to_ORC,self.Beta))
        
        # The error for w_bi is within the rungeKutta function
        self.w_bi = self.rungeKutta_w(self.t, self.w_bi, self.t+self.dt, self.dh, self.r_sat)
        
        self.w_bo = self.w_bi - np.matmul(self.A, np.array(([0],[self.wo],[0])))
        
        self.q = self.rungeKutta_q(self.t, self.q, self.t+self.dt, self.dh)
        
        self.t += self.dt

        if self.fault != "None":
            self.r_sat = self.Sensors_general_fault.High_noise(self.r_sat)
            self.S_EIC = self.Sensors_general_fault.High_noise(self.S_EIC)
            self.S_EIC = self.Sun_sensor_fault.Catastrophic(self.S_EIC)
            self.S_EIC = self.Sun_sensor_fault.Erroneous(self.S_EIC)
            self.S_EIC = self.Common_data_transmission_fault.Bit_flip(self.S_EIC)
            self.S_EIC = self.Common_data_transmission_fault.Sign_flip(self.S_EIC)
            self.S_EIC = self.Common_data_transmission_fault.Insertion(self.S_EIC)
            self.Beta = self.Magnetometer_fault.Stop(self.Beta)
            self.Beta = self.Magnetometer_fault.Interference(self.Beta)
            self.Beta = self.Common_data_transmission_fault.Bit_flip(self.Beta)
            self.Beta = self.Common_data_transmission_fault.Sign_flip(self.Beta)
            self.Beta = self.Common_data_transmission_fault.Insertion(self.Beta)
        
        self.update()
        
        return self.w_bi, self.q, self.A, r_EIC, self.sun_in_view

    def update(self):
        self.Orbit_Data["Magnetometer x"].append(self.B[0])
        self.Orbit_Data["Magnetometer z"].append(self.B[1])
        self.Orbit_Data["Magnetometer y"].append(self.B[2])
        self.Orbit_Data["Sun x"].append(self.S_b[0][0])
        self.Orbit_Data["Sun y"].append(self.S_b[1][0])
        self.Orbit_Data["Sun z"].append(self.S_b[2][0])
        self.Orbit_Data["Earth x"].append(self.r_sat_sbc[0])
        self.Orbit_Data["Earth y"].append(self.r_sat_sbc[1])
        self.Orbit_Data["Earth z"].append(self.r_sat_sbc[2])
        self.Orbit_Data["Angular momentum of wheels x"].append(self.angular_momentum[0][0])
        self.Orbit_Data["Angular momentum of wheels y"].append(self.angular_momentum[1][0])
        self.Orbit_Data["Angular momentum of wheels z"].append(self.angular_momentum[2][0])
        self.Orbit_Data["Sun in view"].append(self.sun_in_view)
        if self.sun_in_view == False and self.fault == "Catastrophic_sun" or self.fault == "Erroneous":
            self.Orbit_Data["Current fault"].append("None")
        else:
            self.Orbit_Data["Current fault"].append(self.fault)
        temp = list(self.zeros)
        temp[Fault_names_to_num[self.fault] - 1] = 1
        self.Orbit_Data["Current fault numeric"].append(temp)
        self.Orbit_Data["Current fault binary"].append(0 if self.fault == "NoneNone" else 1)

def loop(index, Data, orbit_descriptions):
    if SET_PARAMS.Display:
        satellite = view.initializeCube(SET_PARAMS.Dimensions)
        pv = view.ProjectionViewer(1920, 1080, satellite)

    """
    for i in range(int(SET_PARAMS.Number_of_orbits*SET_PARAMS.Period/(SET_PARAMS.faster_than_control*SET_PARAMS.Ts))):
        w, q, A, r, sun_in_view = D.rotation()
        if SET_PARAMS.Display and i%SET_PARAMS.skip == 0:
            pv.run(w, q, A, r, sun_in_view)
    """
    while D.fault == "None":
        w, q, A, r, sun_in_view = D.rotation()
        if SET_PARAMS.Display and i%SET_PARAMS.skip == 0:
            pv.run(w, q, A, r, sun_in_view)
    
    """
    for i in range(int(SET_PARAMS.Number_of_orbits*SET_PARAMS.Period/(SET_PARAMS.faster_than_control*SET_PARAMS.Ts))):
        w, q, A, r, sun_in_view = D.rotation()
        if SET_PARAMS.Display and i%SET_PARAMS.skip == 0:
            pv.run(w, q, A, r, sun_in_view)
    """

    if SET_PARAMS.Visualize and SET_PARAMS.Display == False:
        visualize_data(D)
        
    print("Number of multiple orbits", index)  
    Data[index] = D.Orbit_Data
    orbit_descriptions[index] =str(index)
    df = pd.DataFrame(D.Orbit_Data, columns = D.Orbit_Data.keys())
    df.to_csv(SET_PARAMS.csv_filename + str(i) + ".csv")


if __name__ == "__main__":
    # FOR ALL OF THE FAULTS RUN A NUMBER OF ORBITS TO COLLECT DATA
    sense = Sensors()
    threads = []
    """
    for i in range(SET_PARAMS.Number_of_multiple_orbits):
        threads.append(i)

    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = executor.map(loop, threads)
    """
    manager = multiprocessing.Manager()
    Data = manager.dict()
    orbit_descriptions = manager.dict()

    for i in range(SET_PARAMS.Number_of_multiple_orbits):
        D = Dynamics(i)
        t = multiprocessing.Process(target=loop, args=(i,Data, orbit_descriptions))
        threads.append(t)
        t.start()
    
    dataframe = []
    i = 0
    for process in threads:
        process.join()
        dataframe.append(pd.DataFrame.from_dict(Data[i]))
        i += 1

    if SET_PARAMS.Save_excel_file:
        save_as_excel(Data.values(), orbit_descriptions.values())
    elif SET_PARAMS.Save_csv_file:
        save_as_csv(Data.values(), orbit_descriptions.values())
    else:
        save_as_pickle(Data.values())
"""
if __name__ == "__main__":
    D = Dynamics()
    sense = Sensors()
    print(D.Reaction_wheel_fault.Reliability_area)
    for j in range(SET_PARAMS.Number_of_multiple_orbits):
        for i in range(1,int(SET_PARAMS.Number_of_orbits*SET_PARAMS.Period/(SET_PARAMS.faster_than_control*SET_PARAMS.Ts))):
            if D.Reaction_wheel_fault.Failure_Reliability_area(i) == True:
                print("Failed")
"""

