import numpy as np 
import math
import Quaternion_functions
from Parameters import SET_PARAMS

pi = math.pi

class Control:
    def __init__(self):
        self.Kp = SET_PARAMS.Kp
        self.Kd = SET_PARAMS.Kd
        self.w_ref = SET_PARAMS.w_ref
        self.q_ref = SET_PARAMS.q_ref
        self.N_max = SET_PARAMS.N_ws_max
        self.first = False

    def control(self, w, q, Inertia, B):
        N = self.control_wheel(w, q, Inertia) + self.magnetic_torquers(B)
        return N

    def control_wheel(self, w, q, Inertia):
        q_error = Quaternion_functions.quaternion_error(q, self.q_ref)
        w_error = self.w_ref - w
        N = np.clip(np.reshape((self.Kp * np.matmul(Inertia, q_error[0:3])),(3,1)) - self.Kd * np.matmul(Inertia, w_error), -self.N_max,self.N_max)
        return N
    
    def magnetic_torquers(self, B):
        if self.first == 0:
            self.first = True
            My = 0.0
            Beta = 0.0
        else:
            Beta = np.arccos(B[1]/np.linalg.norm(B))
            My = SET_PARAMS.Kd_magnet * (Beta - self.Beta)/SET_PARAMS.Ts

        M = np.array(([[0],[My],[0]]))
        self.Beta = Beta
        N = np.reshape(np.matmul(M,np.reshape(B,(1,3)))[1,:],(3,1))
        return N