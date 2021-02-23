import numpy as np 
import math
import Quaternion_functions
from Parameters import SET_PARAMS

pi = math.pi

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

class Control:
    def __init__(self):
        self.Kp = SET_PARAMS.Kp
        self.Kd = SET_PARAMS.Kd
        self.w_ref = SET_PARAMS.w_ref
        self.q_ref = SET_PARAMS.q_ref
        self.N_max = SET_PARAMS.N_ws_max
        self.first = False
        self.angular_momentum = np.zeros((3,1))

    def control(self, w, q, Inertia, B, error = False):
        if error:
            N = self.magnetic_torquers(B, w)       #control_wheels do not respond
            self.angular_momentum = np.zeros((3,1))
        else:
            N = self.control_wheel(w, q, Inertia)
            self.angular_momentum += N*SET_PARAMS.Ts
            self.angular_momentum = np.clip(self.angular_momentum, -SET_PARAMS.h_ws_max, SET_PARAMS.h_ws_max)
            N = self.magnetic_torquers(B, w)
        return N, self.angular_momentum

    def control_wheel(self, w, q, Inertia):
        q_error = Quaternion_functions.quaternion_error(q, self.q_ref)
        w_error = self.w_ref - w
        N = np.clip(np.reshape((self.Kp * np.matmul(Inertia, q_error[0:3])),(3,1)) - self.Kd * np.matmul(Inertia, w_error), -self.N_max,self.N_max)
        return N
    
    def magnetic_torquers(self, B, w):
        if self.first == 0:
            self.first = True
            My = 0.0
            Mx = 0.0
            Mz = 0.0
            Beta = 0.0
        else:
            Beta = np.arccos(B[1]/np.linalg.norm(B))
            My = SET_PARAMS.Kd_magnet * (Beta - self.Beta)/SET_PARAMS.Ts
            if B[2] > B[0]:
                Mx = SET_PARAMS.Ks_magnet * (w[1][0] - self.w_ref[1][0])*np.sign(B[2])
                Mz = 0
            else:
                Mz = SET_PARAMS.Ks_magnet * (self.w_ref[1][0] - w[1][0])*np.sign(B[0])
                Mx = 0

        M = np.array(([[Mx],[My],[Mz]]))
        self.Beta = Beta
        N = np.reshape(np.matmul(M,np.reshape(B,(1,3)))[1,:],(3,1))
        N = np.clip(N, -SET_PARAMS.M_magnetic_max, SET_PARAMS.M_magnetic_max)
        return N