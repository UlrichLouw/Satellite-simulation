import numpy as np 
import math
import Simulation.Quaternion_functions as Quaternion_functions
from Simulation.Parameters import SET_PARAMS

pi = math.pi

class Control:
    def __init__(self):
        self.Kp = SET_PARAMS.Kp
        self.Kd = SET_PARAMS.Kd
        self.w_ref = SET_PARAMS.w_ref
        self.q_ref = SET_PARAMS.q_ref
        self.N_max = SET_PARAMS.N_ws_max
        self.first = True
        self.mode = SET_PARAMS.Mode

    def control(self, w, q, Inertia, B, error = False):
        if self.mode == "Nominal":   # Normal operation
            N_magnet = np.zeros((3,1))
            N_wheel = self.control_wheel(w, q, Inertia)

        elif self.mode == "Safe":    # Detumbling mode
            N_magnet = self.magnetic_torquers(B, w)
            N_wheel = np.zeros((3,1))
    
        return N_magnet, N_wheel

    def control_wheel(self, w, q, Inertia):
        q_error = Quaternion_functions.quaternion_error(q, self.q_ref)
        w_error = self.w_ref - w
        N = np.clip(np.reshape((self.Kp * np.matmul(Inertia, q_error[0:3])),(3,1)) - self.Kd * np.matmul(Inertia, w_error), -self.N_max,self.N_max)
        return N
    
    def magnetic_torquers(self, B, w):
        if self.first == True:
            self.first = False
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