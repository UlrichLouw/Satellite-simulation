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

    def control(self, w, q, Inertia):
        q_error = Quaternion_functions.quaternion_error(q, self.q_ref)
        w_error = self.w_ref - w
        N = np.clip(np.reshape((self.Kp * np.matmul(Inertia, q_error[0:3])),(3,1)) - self.Kd * np.matmul(Inertia, w_error), -self.N_max,self.N_max)
        return N