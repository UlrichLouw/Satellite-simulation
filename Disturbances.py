from Parameters import SET_PARAMS

class Disturbances:
    def __init__(self):
        pass

    def Gravity_gradient_func(self, AOB):
        zoB = AOB * np.array(([0],[0],[1]))
        kgx = 3 * SET_PARAMS.wo**2 * (SET_PARAMS.Iz - SET_PARAMS.Iy)
        kgy = 3 * SET_PARAMS.wo**2 * (SET_PARAMS.Ix - SET_PARAMS.Iz)
        kgz = 3 * SET_PARAMS.wo**2 * (SET_PARAMS.Iy - SET_PARAMS.Ix)
        Ngg = np.array(([kgx*AOB[1][2]*AOB[2][2]],[kgy*AOB[0][2]*AOB[2][2]],[kgz*AOB[0][2]*AOB[1][2]]))
        
        return Ngg

    def Aerodynamic(self, v_ab, Surface_area_i, incidence_angle,
                        normal_accommodation, ratio_of_molecular_exit, offset_vector, 
                        unit_inward_normal, h, h_o, H):
        Ai = Surface_area_i
        alpha_i = incidence_angle
        sigma_t = 0.8       # tangential_accommodation
        sigma_n = 0.8  # normal_accommodation
        S = 0.05        # ratio_of_molecular_exit
        r_pi = offset_vector
        n_i = unit_inward_normal
        p = 0.5 * (p_o * np.exp(-(h-h_o)/H))
        va = np.arraty(([0],[0],[w_e])) * r_sat - v_sat
        N_aero = 0
        for i in range(n):
            N_aero = N_aero + p * np.linalg.norm(va)**2 * Ai * np.heaviside(np.cos(alpha_i))*np.cos(alpha_i)*sigma_t*(r_pi) + (sigma_n * S + (2-sigma_n - sigma_t)*np.cos(alpha_i)*(r_pi))

        return N_aero

    def Wheel_Imbalance(self, rotation_rate):
        # Static imbalance with one mass a distance r from the centre of the flywheel
        # m - mass of imbalance in kg
        # r - distance from centre of flywheel in m
        Us = 0.0208 / 1e-5 #For the RW-0.06 wheels in kg/m; Us = m*r
        omega = rotation_rate #rad/s rotation rate of wheel
        phi_s = 0 #arbitrary phase
        w_z = np.array(([0,0,0])) #np.array 

        F_zs = Us * omega**2 * np.array(([np.sin(omega*t + phi_s)],[np.cos(omega*t + phi_s)], [0]))

        N_zs = np.matmul(w_z, F_zs)

        # Dynamic imbalance with two masses seperated by 180 degrees and distance, d
        # d - width of flywheel in m
        # r - distance from flywheel
        # m - mass in kg
        Ud = 0.0208 / 1e-7 #For the RW-0.06 wheels in kg/m^2; Ud = m*r*d
        N_zd = Ud*omega**2 * np.array(([np.sin(omega*t + phi_s)],[np.cos(omega*t + phi_s)], [0]))

        return N_zs, N_zd
