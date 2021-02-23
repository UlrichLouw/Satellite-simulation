import numpy as np
import Controller
from Disturbances import Disturbances
from Parameters import SET_PARAMS
from Sensors import Sensors
import matplotlib.pyplot as plt
import time
import csv
import pandas as pd

if SET_PARAMS.Display:
    import Satellite_display as view

csv_columns = ["Sun", "Magnetometer", "Earth","Wheel speed", "Sun in view"]

# The Orbit_Data dictionary is used to store all the measurements for each timestep (Ts)
# Each Orbit has an induced fault within the ADCS. 
Orbit_Data = {
    "Sun": {"x": [],"y": [], "z": []},              #S_o measurement (vector of sun in ORC)
    "Magnetometer": {"x": [],"y": [], "z": []},     #B vector in SBC
    "Earth": {"x": [],"y": [], "z": []},            #Satellite position vector in ORC
    "Wheel speed": {"x": [],"y": [], "z": []},      #Wheel angular velocity of each reaction wheel
    "Sun in view": []       #True or False values depending on whether the sun is in view of the satellite
}

# The Data dictionary is used to store each orbit
Data = {
    "Orbit": [],            #The data of a single orbit
    "Fault": []             #The fault induced during the orbit
}

# File names for the storage of the data attained during the simulation
csv_file = "Fault-prediction/Data_files/Faults.csv"
pickle_filename = "Fault-prediction/Data_files/Faults.pkl"

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

def save_as_csv():
    df = pd.DataFrame.from_dict(Data)
    df.to_csv(index=False)
    with open(csv_file, 'w') as csvfile:
        writer = csv.writer(csvfile)
        for key, value in Data.items():
            writer.writerow([key, value])

def save_as_pickle():
    Data_frame = pd.DataFrame.from_dict(Data)
    Data_frame.to_pickle(pickle_filename)

def visualize_data():
    for i in Orbit_Data:
        if i == "Sun in view":
            fig = plt.figure()
            y = np.array((Orbit_Data[i]))
            plt.plot(y)
            plt.title(i)
            plt.pause(1)
        else:
            fig = plt.figure()
            for j in Orbit_Data[i]:
                y = np.array((Orbit_Data[i][j]))
                plt.plot(y)

            plt.legend(Orbit_Data[i])
            plt.title(i)
            plt.pause(1)

class Dynamics:
    def __init__(self):
        self.dist = Disturbances()
        self.w_bi = SET_PARAMS.wbi
        self.wo = SET_PARAMS.wo
        self.q = SET_PARAMS.quaternion_initial
        self.t = SET_PARAMS.time
        self.dt = SET_PARAMS.Ts
        self.dh = self.dt/10
        self.Ix = SET_PARAMS.Ix
        self.Iy = SET_PARAMS.Iy
        self.Iz = SET_PARAMS.Iz
        self.Inertia = np.identity(3)*np.array(([self.Ix, self.Iy, self.Iz]))
        self.Iw = SET_PARAMS.Iw
        self.angular_momentum = SET_PARAMS.initial_angular_momentum
        self.faster_than_control = SET_PARAMS.faster_than_control
        self.control = Controller.Control()
        self.angular_momentum = np.zeros((3,1))

    def rungeKutta_h(self, x0, angular, x, h, N_control):
        angular_momentum_derived = N_control
        n = int(np.round((x - x0)/h))

        y = angular
        for i in range(n):
            k1 = h*(y) 
            k2 = h*((y) + 0.5*k1) 
            k3 = h*((y) + 0.5*k2) 
            k4 = h*((y) + 0.5*k3) 

            y = y + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)
    
            x0 = x0 + h; 
        
        return y

    def rungeKutta_w(self, x0, w, x, h, r_sat, v_sat):
        N_aero = 0 #self.dist.Aerodynamic(self.A, self.A_EIC_to_ORC)
        N_rw = np.reshape(self.dist.Wheel_Imbalance(w, x),(3,1))
        Ngg = self.dist.Gravity_gradient_func(self.A)
        N_disturbance = Ngg + N_aero + N_rw #Ignore gyroscope
        N_control = self.control.control(w, self.q, self.Inertia, self.B)
        self.angular_momentum = self.angular_momentum + N_control*self.dt
        """        
        self.angular_momentum = self.rungeKutta_h(x0, self.angular_momentum, x, h, N_control) 
        """
        n = int(np.round((x - x0)/h))

        y = w
        for i in range(n):
            k1 = h*((-w * np.matmul(self.Inertia,w)) + (-N_control+N_disturbance)) 
            k2 = h*(((-w * np.matmul(self.Inertia,w)) + (-N_control+N_disturbance)) + 0.5*k1) 
            k3 = h*(((-w * np.matmul(self.Inertia,w)) + (-N_control+N_disturbance)) + 0.5*k2) 
            k4 = h*(((-w * np.matmul(self.Inertia,w)) + (-N_control+N_disturbance)) + 0.5*k3) 

            y = y + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)
    
            x0 = x0 + h; 
        
        return y


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

    def rotation(self):
        self.r_sat, v_sat, self.A_EIC_to_ORC, r_EIC = sense.satellite_vector(self.t*self.faster_than_control)
        self.S_EIC, self.sun_in_view = sense.sun(self.t*self.faster_than_control)
        self.S_o = np.matmul(self.A_EIC_to_ORC, self.S_EIC)
        self.A = np.matmul(self.A_EIC_to_ORC, Transformation_matrix(self.q))
        self.Beta = sense.magnetometer(self.t) 
        self.Beta += np.random.normal(0,SET_PARAMS.magnetometer_noise,self.Beta.shape)
        self.B = np.matmul(self.A,np.matmul(self.A_EIC_to_ORC,self.Beta))
        self.w_bi = self.rungeKutta_w(self.t, self.w_bi, self.t+self.dt, self.dh, self.r_sat, v_sat)
        self.w_bo = self.w_bi - np.matmul(self.A, np.array(([0],[self.wo],[0])))
        self.q = self.rungeKutta_q(self.t, self.q, self.t+self.dt, self.dh)
        self.t += self.dt
        self.update()
        return self.w_bi, self.q, self.A, r_EIC, self.sun_in_view

    def update(self):
        i = 0
        for index in Orbit_Data["Sun"]:
            Orbit_Data["Magnetometer"][index].append(self.B[i])
            Orbit_Data["Sun"][index].append(self.S_o[i][0])
            Orbit_Data["Earth"][index].append(self.r_sat[i])
            Orbit_Data["Wheel speed"][index].append(self.angular_momentum[i][0])
            i += 1
        Orbit_Data["Sun in view"].append(self.sun_in_view)


if __name__ == "__main__":
    D = Dynamics()
    sense = Sensors()

    if SET_PARAMS.Display:
        satellite = view.initializeCube(SET_PARAMS.Dimensions)
        pv = view.ProjectionViewer(1920, 1080, satellite)

    for index in SET_PARAMS.Fault_names:     
        for i in range(int(SET_PARAMS.Period/(SET_PARAMS.faster_than_control*SET_PARAMS.Ts))):
            w, q, A, r, sun_in_view = D.rotation()
            if SET_PARAMS.Display and i%SET_PARAMS.skip == 0:
                pv.run(w, q, A, r, sun_in_view)

        Data["Orbit"].append(Orbit_Data)
        Data["Fault"].append(SET_PARAMS.Fault_names[index])

        if SET_PARAMS.visualize:
            visualize_data()

    if SET_PARAMS.Save_csv_file:
        save_as_csv()
    else:
        save_as_pickle()

