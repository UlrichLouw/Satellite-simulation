import Parameters
import math
import numpy as np 
from scipy import signal
import control

pi = math.pi

class Control:
    def __init__(self):
        self.zero3 = np.zeros((3,3))
        self.zero63 = np.zeros((6,3))
        self.zero66 = np.zeros((6,6))
        self.eye3 = np.identity(3)
        self.eye6 = np.identity(6)

    def Magnetic(self, Select, Wbo, Qe, Magentic_model, CON, H, Wref, h_wheel):
        Nm = 1
        return Nm

    def Reaction_Wheel(self, Is, ws, hw, tau_w, Ts):
        Is1 , Is2, Is3 = Is

        hw1, hw2, hw3 = hw[:,0]
        ws1, ws2, ws3 = ws[:,0]

        Is11, Is12, Is13 = Is[0,:]
        Is21, Is22, Is23 = Is[1,:]
        Is31, Is32, Is33 = Is[2,:]

        Aw1 = np.matmul(np.array(([0, Is31, Is21], [-2*Is31, -Is32, -Is33 + Is11], [2*Is21, Is22 - Is11, Is23])), ws)
        Aw2 = np.matmul(np.array(([Is31, 2*Is32, Is33 - Is22], [-Is32, 0, Is12], [Is22 - Is11, 2*Is12, - Is13])), ws)
        Aw3 = np.matmul(np.array(([-Is21, Is33 - Is22, -2*Is23], [Is11 - Is33, Is12, 2*Is13], [Is23, - Is13, 0])), ws)
        Aw = np.array(([Aw1, Aw2, Aw3]))[0]

        Is_invert = np.linalg.inv(Is)

        Ahw = np.array(([0, -hw3, hw2], [hw3, 0, -hw1], [-hw2, hw1, 0]))
        Awh = np.array(([0, ws3, -ws2], [-ws3, 0, ws1], [ws2, -ws1, 0]))

        Aww = Aw + Ahw

        Iinvaww = np.matmul(Is_invert, Aww)

        A = np.concatenate((Iinvaww, self.zero3, Iinvaww, -Is_invert))
        A = np.concatenate((A, np.concatenate((0.5 * self.eye3, self.zero3, self.zero3, self.zero3))), axis = 1)
        A = np.concatenate((A, np.concatenate((self.zero3, self.zero3, self.zero3, self.eye3))), axis = 1)
        A = np.concatenate((A, np.concatenate((self.zero3, self.zero3, self.zero3, -self.eye3 / tau_w))), axis = 1)

        Bu = np.concatenate((self.zero3,self.zero3,self.zero3,-self.eye3 / tau_w))

        Bd = np.concatenate((Is_invert, self.zero3, self.zero3, self.zero3))

        C = np.concatenate((self.eye6, self.zero66), axis =1)

        D = self.zero63

        Satellite_control_model = control.StateSpace(A, Bu, C, D)

        #Satellite_control_model = signal.StateSpace(A, Bu, C, D, dt = Ts)

        Satellite_disturbance_model = signal.StateSpace(A, Bd, C, D, dt = Ts)

        #h_ref = Control_Parameters.h_ref

        #K = h_ref / SET_PARAMS.Iw

        #K = K * 60/(2*pi)


        return Satellite_control_model, Satellite_disturbance_model

    def Paddle(Select, Wbo, Qe, CON):
        Paddle_angle = 1
        return paddle_angle