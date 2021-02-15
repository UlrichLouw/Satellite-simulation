import numpy as np
import Satellite_display as view
import Controller

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

class Dynamics:
    def __init__(self):
        self.w_bi = np.array(([0.1],[0],[0]))
        self.q = np.array(([0,0,0,1]))
        self.t = 0
        self.dt = 0.01
        self.dh = 0.001
        self.N = np.array(([0.1], [0.1],[0.1]))
        self.Ix = 2
        self.Iy = 2
        self.Iz = 2
        self.Inertia = np.array(([self.Ix, self.Iy, self.Iz]))
        self.Iw = 0.05
        self.angular_momentum = np.zeros((3,1))

    def rungeKutta_h(self, x0, angular, x, h, N_control):
        angular_momentum_derived = N_control/self.Iw
        n = int(np.round((x - x0)/h))

        y = angular
        for i in range(n):
            k1 = h*(y + angular_momentum_derived) 
            k2 = h*(y + angular_momentum_derived + 0.5*k1) 
            k3 = h*(y + angular_momentum_derived + 0.5*k2) 
            k4 = h*(y + angular_momentum_derived + 0.5*k3) 

            y = y + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)
    
            x0 = x0 + h; 
        
        return y

    def rungeKutta_w(self, x0, w, x, h):
        Ngg = 0
        N_aero = 0
        N_rw = 0
        N_disturbance = Ngg + N_aero + N_rw #Ignore gyroscope
        N_control = Controller.Control().control(w, self.q, self.Inertia)
        self.angular_momentum = self.rungeKutta_h(x0, self.angular_momentum, x, h, N_control)
        n = int(np.round((x - x0)/h))

        y = w

        for i in range(n):
            k1 = h*((-w * np.matmul(self.Inertia,w) + self.angular_momentum) + (-N_control+N_disturbance) - self.angular_momentum) 
            k2 = h*(((-w * np.matmul(self.Inertia,w) + self.angular_momentum) + (-N_control+N_disturbance) - self.angular_momentum) + 0.5*k1) 
            k3 = h*(((-w * np.matmul(self.Inertia,w) + self.angular_momentum) + (-N_control+N_disturbance) - self.angular_momentum) + 0.5*k2) 
            k4 = h*(((-w * np.matmul(self.Inertia,w) + self.angular_momentum) + (-N_control+N_disturbance) - self.angular_momentum) + 0.5*k3) 

            y = y + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)
    
            x0 = x0 + h; 
        
        return y


    def rungeKutta_q(self, x0, y0, x, h):
        wx, wy, wz = self.w_bi[:,0]
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
        self.w_bi = self.rungeKutta_w(self.t, self.w_bi, self.t+self.dt, self.dh)
        self.q = self.rungeKutta_q(self.t, self.q, self.t+self.dt, self.dh)
        #print(self.q)
        print(self.w_bi)
        A = Transformation_matrix(self.q)
        self.t += self.dt
        return self.w_bi, self.q, A


if __name__ == "__main__":
    D = Dynamics()
    satellite = view.initializeCube()
    pv = view.ProjectionViewer(640, 480, satellite)
    for i in range(100000):
        w, q, A = D.rotation()
        pv.run(w, q, A)
