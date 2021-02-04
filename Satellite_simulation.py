class SET_PARAMS:
    eccentricicity = 0.0002704
    inclination = 97.2927 #degrees
    Semi_major_axis = 6879.55 #km
    Mean_motion = 15.215 #rev/day
    Period = 5678.7 #seconds

class Sensors:
    def sun(self):
        pass

    def nadir(self):
        pass

    def magnetometer(self):
        pass

class Actuators:
    def active_aerodynamic_roll_control_paddles(self):
        self.AARCP_num = 2

    def magnetic_torgue_rods(self):
        self.MTR_num = 3

    def nano_reaction_wheels(self):
        self.RW_num = 3

class Simulation:
    def __init__(self):
        print('Start Simulation')
        self.altitude = SET_PARAMS().Semi_major_axis

    class ADCS:
        def init(self):
            print('ADCS initiated')

        def step(action):
            print('action taken')

        def reset(self):
            print('reset environment')

if __name__ == "__main__":
    sim = Simulation()