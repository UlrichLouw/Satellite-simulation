import numpy as np
import Controller
from Disturbances import Disturbances
from Parameters import SET_PARAMS
from Sensors import Sensors
import matplotlib.pyplot as plt
import time
import csv
import pandas as pd

Fault_names_to_num = {
    "NoneNone": 1,
    "Sun sensorx": 2,
    "Sun sensory": 3, 
    "Sun sensorz": 4,
    "Magnetometerx": 5, 
    "Magnetometery": 6, 
    "Magnetometerz": 7,
    "Earth sensorx": 8, 
    "Earth sensory": 9, 
    "Earth sensorz": 10,
    "Reaction wheelx": 11, 
    "Reaction wheely": 12, 
    "Reaction wheelz": 13,
    "Controlall": 14
}
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

# Function to save a pickle file of simulation data
def save_as_pickle(Data):
    Data_frame = pd.DataFrame(Data, columns = Data.keys())
    Data_frame.to_pickle(SET_PARAMS.pickle_filename)

# Function to visualize data as graphs
def visualize_data(D):
    for i in D.Orbit_Data:
        if i == "Current fault":
            pass
        else:
            if i == "Sun in view":
                fig = plt.figure()
                y = np.array((D.Orbit_Data[i]))
                plt.plot(y)
                plt.title(i)
                plt.show() #plt.pause(0.001)
            else:
                fig = plt.figure()
                for j in D.Orbit_Data[i]:
                    y = np.array((D.Orbit_Data[i][j]))
                    plt.plot(y)

                plt.legend(D.Orbit_Data[i])
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
    def __init__(self):
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
        self.angular_momentum = SET_PARAMS.initial_angular_momentum # Angular momentum of satellite wheels
        self.faster_than_control = SET_PARAMS.faster_than_control   # If it is required that satellite must move faster around the earth than Ts
        self.control = Controller.Control()         # Controller.py is used for control of satellite    
        self.fault = "NoneNone"                         # Current fault of the system
        # All of the faults can be implemented on a single sensor in a specified axis (3)
        self.Earth_sensor_fault = [False, False, False]     
        self.Reaction_wheel_fault = [False, False, False]
        self.Sun_sensor_fault = [False, False, False]
        self.Magnetometer_fault = [False, False, False]
        self.Control_fault = False
        # The Orbit_Data dictionary is used to store all the measurements for each timestep (Ts)
        # Each Orbit has an induced fault within the ADCS. 
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
        self.zeros = np.zeros((14,), dtype = int)

    # Function to calculate the satellite angular velocity based on the derivative thereof
    def rungeKutta_w(self, x0, w, x, h, r_sat):
        N_aero = 0 #self.dist.Aerodynamic(self.A, self.A_EIC_to_ORC)
        N_rw = np.reshape(self.dist.Wheel_Imbalance(w, x),(3,1))    # Disturbance of a reaction wheel imbalance
        Ngg = self.dist.Gravity_gradient_func(self.A)               # Disturbance of gravity gradient
        
        # Control torques implemented due to the control law
        N_control_magnetic, N_control_wheel = self.control.control(w, self.q, self.Inertia, self.B, self.Control_fault)

        # Change the reaction wheel and reaction torque to 0 (depending on the fault)
        if any(self.Reaction_wheel_fault):
            self.angular_momentum[np.where(self.Reaction_wheel_fault)[0]] = 0
            N_control_wheel[np.where(self.Reaction_wheel_fault)[0]] = 0
        else:
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
        
        return np.clip(y, -SET_PARAMS.wheel_angular_d_max, SET_PARAMS.wheel_angular_d_max)

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

    def Fault_implementation(self, index, direction):
        if int(self.t) == SET_PARAMS.Fault_time:
            if direction == "x":
                i = 0
            elif direction == "y":
                i = 1
            elif direction == "z":
                i = 2
            
            if index == "Sun sensor":
                self.Sun_sensor_fault[i] = True

            elif index == "Magnetometer":
                self.Magnetometer_fault[i] = True
            
            elif index == "Earth sensor":
                self.Earth_sensor_fault[i] = True

            elif index == "Reaction wheel":
                self.Reaction_wheel_fault[i] = True
            
            elif index == "Control":
                self.Control_fault = True
            
            self.fault = str(index) + str(direction)

    def rotation(self, index, direction):
        self.Fault_implementation(index, direction)
        self.r_sat, v_sat, self.A_EIC_to_ORC, r_EIC = sense.satellite_vector(self.t*self.faster_than_control, error=self.Earth_sensor_fault)
        self.S_EIC, self.sun_in_view = sense.sun(self.t*self.faster_than_control, self.Sun_sensor_fault)    
        self.A = np.matmul(self.A_EIC_to_ORC, Transformation_matrix(self.q))
        self.r_sat_sbc = np.matmul(self.A, self.r_sat)
        self.S_o = np.matmul(self.A_EIC_to_ORC, self.S_EIC)
        self.S_b = np.matmul(self.A, self.S_o)
        self.Beta = sense.magnetometer(self.t, error = self.Magnetometer_fault) 
        self.B = np.matmul(self.A,np.matmul(self.A_EIC_to_ORC,self.Beta))
        self.w_bi = self.rungeKutta_w(self.t, self.w_bi, self.t+self.dt, self.dh, self.r_sat)
        self.w_bo = self.w_bi - np.matmul(self.A, np.array(([0],[self.wo],[0])))
        self.q = self.rungeKutta_q(self.t, self.q, self.t+self.dt, self.dh)
        self.t += self.dt
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
        self.Orbit_Data["Current fault"].append(self.fault)
        temp = list(self.zeros)
        temp[Fault_names_to_num[self.fault] - 1] = 1
        self.Orbit_Data["Current fault numeric"].append(temp)
        self.Orbit_Data["Current fault binary"].append(0 if self.fault == "NoneNone" else 1)

if __name__ == "__main__":
    # FOR ALL OF THE FAULTS RUN A NUMBER OF ORBITS TO COLLECT DATA
    Data = []
    orbit_descriptions = []
    for index in SET_PARAMS.Fault_names:  
        for direction in SET_PARAMS.Fault_names[index]:
            if SET_PARAMS.Display:
                satellite = view.initializeCube(SET_PARAMS.Dimensions)
                pv = view.ProjectionViewer(1920, 1080, satellite)

            D = Dynamics()
            sense = Sensors()   
            for i in range(int(SET_PARAMS.Number_of_orbits*SET_PARAMS.Period/(SET_PARAMS.faster_than_control*SET_PARAMS.Ts))):
                w, q, A, r, sun_in_view = D.rotation(index, direction)
                if SET_PARAMS.Display and i%SET_PARAMS.skip == 0:
                    pv.run(w, q, A, r, sun_in_view)

            if SET_PARAMS.Visualize and SET_PARAMS.Display == False:
                visualize_data(D)

            Data.append(D.Orbit_Data)
            orbit_descriptions.append(str(index) + str(direction))

    if SET_PARAMS.Save_excel_file:
        save_as_excel(Data, orbit_descriptions)
    else:
        save_as_pickle(Data)

