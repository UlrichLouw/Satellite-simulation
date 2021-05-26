import numpy as np
from Parameters import SET_PARAMS

class model_update():
    def __init__(self):
        self.w_prev = SET_PARAMS.w
        SET_PARAMS.wo 
        # ! self.w_Bo = self.w_prev - self.A_ORC_to_SBC @ np.array(([0],[self.w_prev],[0]))
        self.dw_prev= np.zeros((3,1))
        self.q_prev = SET_PARAMS.quaternion_initial
        self.Ix = SET_PARAMS.Ix                     # Ixx inertia
        self.Iy = SET_PARAMS.Iy                     # Iyy inertia
        self.Iz = SET_PARAMS.Iz                     # Izz inertia
        self.Nm_prev = 0
        self.Ngg_prev = 0
        self.Ngyro_prev = 0
        self.Nw_prev = 0
        self.Ts = SET_PARAMS.Ts
        self.Inertia = np.identity(3)*np.array(([self.Ix, self.Iy, self.Iz]))
        self.covariance_matrix_P_prev = np.zeros((7,7))
        

    ###################################
    # PERFORM THE ANGULAR RATE UPDATE #
    ###################################
    def model_w_update(self):
        self.dw =np.linalg.inv(self.Inertia) @ (self.Nm_prev - self.Nw_prev - self.Ngg_prev - self.Ngyro_prev)

        self.w = self.w + (self.Ts/2)*(3*(self.dw)-self.dw_prev)

    def model_quaternion_update(self, A):
        wx, wy, wz = self.w_Bo[:,0]
        self.omega_k = np.array(([0, wz, -wy, wx], [-wz, 0, wx, wy], [wy, -wx, 0, wz], [-wx, -wy, -wz, 0]))

        norm_w_Bo = np.sqrt((wx**2 + wy**2 + wz**2))

        # ? w_Bo = self.w - A*np.array(([0],[self.wo],[0]))

        kq = self.Ts * norm_w_Bo
        
        temp = np.cos(kq) * np.identity(4) + (1/norm_w_Bo) * np.sin(kq)*self.omega_k

        self.q = temp @ self.q_prev

    def measurement_update(self):
        pass

    def model_update(self, Nm, Nw, Ngg, Ngyro, A, hw):
        self.model_w_update()

        self.model_quaternion_update(A)
        
        self.state_vector = np.transpose(np.array((self.w, self.q)))

        self.F_t()

        self.phi = np.identity(7) + self.Ts*self.Ft + 0.5*self.Ts*self.Ft**2

        self.covariance_matrix_Q = np.array(([[SET_PARAMS.RW_sigma_x,0,0,0,0,0,0],[0,SET_PARAMS.RW_sigma_y,0,0,0,0,0],[0,0,SET_PARAMS.RW_sigma_z,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0]],[0,0,0,0,0,0,0]]))

        self.covariance_matrix_P = self.phi * self.covariance_matrix_P_prev * self.phi.T + self.covariance_matrix_Q

        left = self.covariance_matrix_P @ self.Hk.T

        # ! Change R_k since it is not yet implemented
        # ! R_k is defined in notes on Appendix B.45
        # ! Thereafter update P_k 
        # ! Self.Hk is also not defined

        

        right = self.Hk @ self.covariance_matrix_P @ self.Hk.T + self.Rk

        self.Kk = left @ np.linalg.inv(right)

        self.measurement_update()

        self.update(Nm, Nw, Ngg, Ngyro)

    def F_t(self, h):
        wx, wy, wz = self.w[:,0]
        hx, hy, hz = h[:,0]

        Ix = self.Ix
        Iy = self.Iy
        Iz = self.Iz

        a11 = 0
        a12 = (wz*(Iy - Iz) - hz)/Ix
        a13 = (wy*(Iy - Iz) + hy)/Ix

        a21 = (wz*(Iz - Ix) + hz)/Iy
        a22 = 0
        a23 = (wx*(Iz - Ix) - hx)/Iy

        a31 = (wy*(Ix - Iy) - hy)/Iz
        a32 = (wx*(Ix - Iy) + hx)/Iz
        a33 = 0
        
        TL = np.array(([[a11, a12, a13], [a21, a22, a23], [a31, a32, a33]]))

        kgx = SET_PARAMS.kgx
        kgy = SET_PARAMS.kgy
        kgz = SET_PARAMS.kgz
        
        K = np.array(([[2*kgx, 0, 0],[0 , 2*kgy, 0], [0 , 0 , 2*kgz]]))

        q1 = self.q[0]
        q2 = self.q[1]
        q3 = self.q[2]
        q4 = self.q[3]

        A13 = self.A[0, 2]
        A23 = self.A[1, 2]
        A33 = self.A[2, 2]

        d1 = np.array(([[(-q1*A23 + q4*A33)/Ix], [(-q1*A13 + q3*A33)/Iy], [(q3*A23 + q4*A13)/Iz]]))

        d2 = np.array(([[(-q2*A23 + q3*A33)/Ix], [(-q2*A13 - q4*A33)/Iy], [(-q4*A23 + q3*A13)/Iz]]))

        d3 = np.array(([[(q3*A23 + q2*A33)/Ix], [(q3*A13 + q1*A33)/Iy], [(q1*A23 + q2*A13)/Iz]]))

        d4 = np.array(([[(q4*A23 + q1*A33)/Ix], [(q4*A13 - q2*A33)/Iy], [(-q2*A23 + q1*A13)/Iz]]))

        D = np.concatenate((d1, d2, d3, d4))

        TR = K @ D
        
        BL = 0.5 * np.array(([[q4, -q3, q2],[q3, q4, -q1],[-q2, q1, q4],[-q1, -q2, -q3]]))

        BR = 0.5 * self.omega_k + self.wo * np.array(([[q1*q3, q1*q4, 1 - q1**2, -q1*q2],[q2*q3, q2*q4, -q1*q2, 1-q2**2],[-(1-q3**2), q3*q4, -q1*q3, -q2*q3],[q3*q4, -(1-q4**2), -q1*q4, -q2*q4]]))

        T = np.concatenate((TL, TR))

        B = np.concatenate((BL, BR))

        self.Ft = np.concatenate(T, B, axis=1)

    def measured_w_update(self):
        pass

    def update(self, Nm, Nw, Ngg, Ngyro):
        self.Nw_prev = Nw
        self.dw_prev = self.dw
        self.Nm_prev = Nm
        self.Ngg_prev = Ngg
        self.Ngyro_prev = Ngyro