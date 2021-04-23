class Actuators:
    def active_aerodynamic_roll_control_paddles(self):
        self.AARCP_num = 2

    def magnetic_torgue_rods(self):
        self.MTR_num = 3

    def nano_reaction_wheels(self):     #RW-0.06 Sinclair Interplanetary
        self.RW_num = 4
        a = np.sqrt(1/3)
        b = np.sqrt(2,3)
        Awar = np.array(([a, b, 0],[a, -b, 0],[-a, 0, -b],[-a, 0, b]))
        w_w_ini = [100, 100, 100, 100]