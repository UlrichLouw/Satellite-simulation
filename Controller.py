import numpy as np 
import math
import Quaternion_functions

pi = math.pi

class Control:
    def __init__(self):
        self.Kp = 0.3
        self.Kd = 0.1
        self.w_ref = np.zeros((3,1))
        self.q_ref = Quaternion_functions.euler_to_quaternion(90,90,0)
        self.control_max = 0.01

    def control(self, w, q, Inertia):
        q_error = Quaternion_functions.quaternion_error(q, self.q_ref)
        w_error = self.w_ref - w
        N = self.Kp * np.matmul(Inertia, q_error[0:3]) + self.Kd * np.matmul(Inertia, w_error)
        return N