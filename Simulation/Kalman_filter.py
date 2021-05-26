import numpy as np
from Parameters import SET_PARAMS
Ts = 1

class RKF():
    def __init__(self):
        self.angular_noise = SET_PARAMS.RW_sigma

        self.measurement_noise =  0.002

        self.P_k = np.eye(3)

        self.sigma_k = np.eye(3)

        self.x_k = SET_PARAMS.wbi

        self.Inertia = SET_PARAMS.Inertia

        self.R_k, self.m_k = measurement_noise_covariance_matrix(self.measurement_noise)
        self.Q_k = system_noise_covariance_matrix(self.angular_noise)

    def Kalman_update(self, v_k, Nm, Nw, Ngyro):
        # Model update
        H_k = Jacobian_H(v_k)
        w_b = delta_angular(self.Inertia, Nm, Nw, Ngyro)
        x_k_update = state_model_update(w_b, self.x_k)
        P_k_update = state_covariance_matrix(self.Q_k, self.P_k, self.sigma_k)

        # Measurement update
        y_k = measurement_state_y(H_k, w_b, self.m_k)
        K_k = Jacobian_K(P_k_update, H_k, self.R_k)
        self.P_k = update_state_covariance_matrix(K_k, H_k, P_k_update)
        self.x_k = state_measurement_update(x_k_update, K_k, y_k, H_k)
        return self.x_k

def measurement_noise_covariance_matrix(measurement_noise):
    m_k = np.array(([[measurement_noise], [measurement_noise], [measurement_noise]]))
    R_k = np.diag([measurement_noise, measurement_noise, measurement_noise]) ** 2
    return R_k, m_k

def system_noise_covariance_matrix(angular_noise):
    Q_k = np.diag([angular_noise,angular_noise,angular_noise]) ** 2
    return Q_k

def Jacobian_H(v_k):
    vx, vy, vz = v_k
    H_k = np.array(([[0, Ts*vz, -Ts*vy], [-Ts*vz, 0, Ts*vx], [Ts*vy, -Ts*vx, 0]]))
    return H_k

def delta_angular(Inertia, Nm, Nw, Ngyro):
    return Inertia @ (Nm - Nw - Ngyro)

def state_model_update(delta_angular, x_prev):
    return x_prev + (Ts/2) * (3 * delta_angular - delta_angular)

def state_covariance_matrix(Q_k, P_k, sigma_k):
    P_k = sigma_k @ P_k @ sigma_k.T + Q_k
    return P_k

def Jacobian_K(P_k, H_k, R_k):
    K_k = P_k @ H_k.T @ np.linalg.inv(H_k @ P_k @ H_k.T + R_k)
    return K_k

def update_state_covariance_matrix(K_k, H_k, P_k):
    P_k = (np.eye(3) - K_k @ H_k) @ P_k
    return P_k

def state_measurement_update(x_k, K_k, y_k, H_k):
    x_k = x_k + K_k @ (y_k - H_k @ x_k)
    return x_k

def measurement_state_y(H_k, w_b, m_k):
    y_k = H_k @ w_b + m_k
    return y_k


if __name__ == "__main__":
    rkf = RKF()
    v_k = np.zeros((3,))
    Nm = np.array(([0, 1, 0])).T
    Nw = np.array(([0.1, 0.1, -0.1])).T
    Ngyro = np.array(([-0.3, 0.15, 0])).T
    for i in range(100):
        rkf.Kalman_update(v_k, Nm, Nw, Ngyro)