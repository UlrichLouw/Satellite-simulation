import pandas as pd
import numpy as np

class Detection:
    def __init__(self):
        sun_threshold = 0.15
        earth_threshold = 0.1
        star_threshold = 0.25
        angular_threshold = 0.1 
        magnetometer_threshold = 0.2

    def Per_Timestep(self, args):
        sun, star, earth, angular_momentum = args

        self.sun_fault(sun)
        self.star_fault(star)
        self.earth_fault(earth)
        self.angular_momentum_fault(angular_momentum)

    #####################################################
    # IF THE THRESHOLD OF THE SUN VECTOR IS LARGER THAT #
    #     A SPECIFIED VALUE THEN IT RETURN AN ERROR     #
    #####################################################
    def sun_fault(self, sun):
        pass

    
    #####################################
    # IF THE STAR VECTOR IS LARGER THAN #
    # A GIVEN THRESHOLD RETURN AN ERROR #
    #####################################
    def star_fault(self, star):
        pass

    
    ########################################
    # IF THE EARTH VECTOR IS LARGER THAN A #
    # GIVEN THRESHOLD THEN RETURN AN ERROR #
    ########################################
    def earth_fault(self, earth):
        pass

    
    ######################################################
    # IF THE ANGULAR MOMENTUM IS LARGER THAN A SPECIFIED #
    #    VALUE OR REMAINS LARGER THAN RETURN AN ERROR    #
    ######################################################
    def angular_momentum_fault(self, angular_moment):
        pass

    
    ########################################################
    # IF THE MAGNETOMETER IS LARGER THAN A SPECIFIED VALUE #
    #        OR REMAINS LARGER THAN RETURN AN ERROR        #
    ########################################################
    def magnetometer_fault(self, magnetometer):
        pass