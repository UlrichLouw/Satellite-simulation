from Parameters import SET_PARAMS
import math

pi = math.pi

class Earth:   
    def __init__(self):
        J_t = SET_PARAMS.J_t
        self.V = np.zeros((2))
        self.first = 1      # This is required for the geomagnetic_field_strength_func to initiate correctly

    def geomagnetic_field_strength_func(self, V):
        if self.first:
            self.V = [V,V]
            self.first = 0
        self.V[0] = V
        delta_V = -(self.V[0]-self.V[1])
        self.V[1] = V

        return delta_V

    def scalar_potential_function(self, orbit_radius, coevelation, longitude):
        rs = orbit_radius
        theta = coevelation
        lambda_ = longitude
        k = SET_PARAMS.k
        g = SET_PARAMS.g
        h = SET_PARAMS.h
        Re = SET_PARAMS.Re
        P = special.lpmn(k,k,theta)[0]
        V = 0
        for n in range(1,k):
            sum_ = 0
            for m in range(n):
                print(P)
                sum_ = sum_ + (g*np.cos(m*lambda_) + h*np.sin(m*lambda_)) * P[m,n]
            
            V = V + (Re/rs)**(n+1) * sum_

        return V

    def Sun_Position(self):
        T_jc = (J_T - 2452545)/36525
        M_o = 357.527723300 + 35999.050340*T_jc     #in degrees
        lambda_Mo = 280.460618400 + 36000.770053610*T_jc        #degrees
        lambda_e = lambda_Mo + 1.914666471*np.sin(M_o*pi/180) + 0.019994643*math.sin(2*M_o*pi/180)      #degrees
        epsilon =  23.439291 - 0.013004200*T_jc                 #degrees
        r_o = 1.00140612 - 0.016708617*np.cos(M_o*pi/180) - 0.000139589*np.cos(2*M_o*pi/180)        #degrees
        rsun = r_o * np.array(([np.cos(lambda_e*pi/180)],[np.cos(epsilon*pi/180)*np.sin(lambda_e*pi/180)],[np.sin(epsilon*pi/180)*np.sin(lambda_e*pi/180)]))

        return rsun*(149597871)     #in km