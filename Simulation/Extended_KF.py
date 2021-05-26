import numpy as np
from Parameters import SET_PARAMS
Ts = SET_PARAMS.Ts

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

# ! Update Q_k, R_k
# ! Qk - System Noise Covariance Matrix page 120. Look at B.27
# ! Rk - Measurement Noise Covariance Matrix

class EKF():
    def __init__(self):
        self.angular_noise = SET_PARAMS.RW_sigma

        self.measurement_noise =  0.002

        self.P_k = np.eye(7)

        self.sigma_k = np.eye(7)

        self.x_k = np.concatenate(SET_PARAMS.wbi, SET_PARAMS.quaternion_initial)

        self.Inertia = SET_PARAMS.Inertia

        self.R_k, self.m_k = measurement_noise_covariance_matrix(self.measurement_noise)
        self.Q_k = system_noise_covariance_matrix(self.angular_noise)
        self.wo = SET_PARAMS.wo


    def Kalman_update(self, vmeas_k, v_ORC_k, Nm, Nw, Ngyro, h):
        # Model update
        H_k = Jacobian_H(v_k)
        dw_bi = delta_angular(self.Inertia, Nm, Nw, Ngyro)
        self.w_bi = state_model_update(dw_bi, self.w_bi)
        self.A_ORC_to_SBC = Transformation_matrix(self.q)
        vmodel_k = self.A_ORC_to_SBC @ v_ORC_k
        self.w_bo = self.w_bi - self.A_ORC_to_SBC @ np.array(([0],[self.wo],[0]))
        omega_k = omega_k_function(self.w_bo)
        kq = kq_function(self.w_bo)
        self.q = state_model_update_quaternion(self.q, kq, omega_k, self.w_bo)

        F_t = F_t_function(h, self.wi, self.Inertia, self.q, omega_k)
        self.sigma_k = sigma_k_function(F_t)

        x_k_update = np.concatenate((self.w_bi, self.q), axis = 1)
        P_k_update = state_covariance_matrix(self.Q_k, self.P_k, self.sigma_k)

        # Measurement update
        e_k = e_k_function(vmeas_k, self.A_ORC_to_SBC, vmodel_k)
        K_k = Jacobian_K(P_k_update, H_k, self.R_k)
        self.P_k = update_state_covariance_matrix(K_k, H_k, P_k_update, self.R_k)
        self.x_k = state_measurement_update(x_k_update, K_k, e_k)

        # Normalize the quaternion matrix
        self.q = self.x_k[3:]
        self.q = np.linalg.norm(self.q)
        self.x_k[3:] = self.q

        return self.x_k


def omega_k_function(w_bo):
    wx, wy, wz = w_bo[:,0]

    W = np.array(([0, wz, -wy, wx], 
                  [-wz, 0, wx, wy], 
                  [wy, -wx, 0, wz], 
                  [-wx, -wy, -wz, 0]))
    return W


def kq_function(w_bo):
    wx, wy, wz = w_bo[:,0]
    kq = Ts/2 * np.sqrt(wx**2 + wy**2 + wz**2)
    return kq


def measurement_noise_covariance_matrix(measurement_noise):
    m_k = np.array(([[measurement_noise], [measurement_noise], [measurement_noise]]))
    R_k = np.diag([measurement_noise, measurement_noise, measurement_noise]) ** 2
    return R_k, m_k


def system_noise_covariance_matrix(angular_noise):
    Q_k = np.diag([angular_noise,angular_noise,angular_noise]) ** 2
    return Q_k


def Jacobian_H(q, vmodel_k):
    q1, q2, q3, q4 = q[:,0]

    zero3 = np.zeros((3,3))
    h1 = 2 * np.array(([[q1, q2, q3], [q2, -q1, q4], [q3, -q4, -q1]])) @ vmodel_k
    h2 = 2 * np.array(([[-q2, q1, -q4], [q1, q2, q3], [q4, q3, -q2]])) @ vmodel_k
    h3 = 2 * np.array(([[-q3, q4, q1],[-q4, -q3, q2],[q1, q2, q3]])) @ vmodel_k
    h4 = 2 * np.array(([[q4, q3, -q2],[-q3, q4, q1],[q2, -q1, q4]])) @ vmodel_k
    H_k = np.concatenate((zero3, h1, h2, h3, h4), axis = 1)
    
    return H_k


def delta_angular(Inertia, Nm, Nw, Ngyro, Ngg):
    return Inertia @ (Nm - Nw - Ngyro - Ngg)


def state_model_update(delta_angular, x_prev):
    return x_prev + (Ts/2) * (3 * delta_angular - delta_angular)


def state_model_update_quaternion(q, kq, sigma_k, w_ob):
    return (np.cos(kq) * np.eye(4) + (1/w_ob)*np.sin(kq)*sigma_k) @ q


def state_covariance_matrix(Q_k, P_k, sigma_k):
    P_k = sigma_k @ P_k @ sigma_k.T + Q_k
    return P_k


def Jacobian_K(P_k, H_k, R_k):
    K_k = P_k @ H_k.T @ np.linalg.inv(H_k @ P_k @ H_k.T + R_k)
    return K_k


def update_state_covariance_matrix(K_k, H_k, P_k, R_k):
    P_k = (np.eye(7) - K_k @ H_k) @ P_k @ (np.eye(7) - K_k @ H_k) + K_k @ R_k @ K_k.T
    return P_k


def state_measurement_update(x_k, K_k, e_k):
    x_k = x_k + K_k @ e_k
    return x_k


def e_k_function(vmeas_k, A, vmodel_k):
    e_k = vmeas_k - A @ vmodel_k
    return e_k


def sigma_k_function(F_t):
    sigma_k = np.eye(7) + Ts*F_t + 0.5 * Ts**2 * F_t ** 2 
    return sigma_k


def F_t_function(h, wi, Inertia, q, omega_k):
    wx, wy, wz = wi[:,0]
    hx, hy, hz = h[:,0]

    Ix, Iy, Iz = Inertia[:,0]

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

    q1, q2, q3, q4 = q[:,0]

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

    BR = 0.5 * omega_k + SET_PARAMS.wo * np.array(([[q1*q3, q1*q4, 1 - q1**2, -q1*q2],[q2*q3, q2*q4, -q1*q2, 1-q2**2],[-(1-q3**2), q3*q4, -q1*q3, -q2*q3],[q3*q4, -(1-q4**2), -q1*q4, -q2*q4]]))

    T = np.concatenate((TL, TR))

    B = np.concatenate((BL, BR))

    Ft = np.concatenate(T, B, axis=1)
    return Ft

if __name__ == "__main__":
    rkf = EKF()
    v_k = np.zeros((3,))
    Nm = np.array(([0, 1, 0])).T
    Nw = np.array(([0.1, 0.1, -0.1])).T
    Ngyro = np.array(([-0.3, 0.15, 0])).T
    Ngg = np.array(([-0.3, 0.15, 0])).T
    for i in range(100):
        rkf.Kalman_update(v_k, Nm, Nw, Ngyro, Ngg)