from Parameters import SET_PARAMS

class Satelite_model:
    # The satellite model refers to the forces and phenomena experience by the satellite during 
    # the mission

    def __init__(self):
        self.Ix = SET_PARAMS.Ix
        self.Iy = SET_PARAMS.Iy
        self.Iz = SET_PARAMS.Iz

    def AOB_func(self, current_quaternion):
        # This function provides a translation matrix from the ORC to SBC

        q1, q2, q3, q4 = current_quaternion
        AOB = np.zeros((3,3))
        AOB[0,0] = q1**2-q2**2-q3**2+q4**2
        AOB[0,1] = 2*(q1*q2 + q3*q4)
        AOB[0,2] = 2*(q1*q3 - q2*q4)
        AOB[1,0] = 2*(q1*q2 - q3*q4)
        AOB[1,1] = -q1**2+q2**2-q3**2+q4**2
        AOB[1,2] = 2*(q2*q3 + q1*q4)
        AOB[2,0] = 2*(q1*q3 + q2*q4)
        AOB[2,1] = 2*(q2*q3 - q1*q4)
        AOB[2,2] = -q1**2-q2**2+q3**2+q4**2

        return AOB

    def q_derived_func(self, current_quaternion, wx, wy, wz, A, wo):
        # For the control systems the derivative of the quaternions is required.

        q1, q2, q3, q4 = current_quaternion
        Q = np.array([[q4, -q3, q2],[q3, q4, -q1],[-q2, q1, q4],[-q1, -q2, -q3]])
        Aw = np.array([[wx + A[0,1]*wo],[wy + A[1,1]*wo],[wz + A[2,1]*wo]])
        q_derived = 0.5 * np.matmul(Q, Aw)

        return q_derived

    def wbo_func(self, A):
        wbo = SET_PARAMS.wbi - np.matmul(A, np.array([0],[-SET_PARAMS.wo],[0]))
        return wbo

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
        