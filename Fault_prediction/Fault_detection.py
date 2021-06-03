import pandas as pd
import numpy as np
import collections
from Simulation.Parameters import SET_PARAMS

######################################################
# BASIC DETECTION WITH THRESHOLDS AND IF STATEMENTS. #
# THE CURRENT METHOD OF FDIR SYSTEMNS ON SATELLITES  #
######################################################

class Basic_detection:
    def __init__(self):
        sun_threshold = 0.15                
        earth_threshold = 0.1
        star_threshold = 0.25
        angular_threshold = 0.1 
        magnetometer_threshold = 0.2
        self.sun_buffer = collections.deque(maxlen = SET_PARAMS.buffer_size)
        self.earth_buffer = collections.deque(maxlen = SET_PARAMS.buffer_size)
        self.star_buffer = collections.deque(maxlen = SET_PARAMS.buffer_size)
        self.magnetometer_buffer = collections.deque(maxlen = SET_PARAMS.buffer_size)
        self.angular_threshold = collections.deque(maxlen = SET_PARAMS.buffer_size)
        self.sensors = {"sun": self.sun_buffer, 
                "earth": self.earth_buffer, 
                "star": self.star_buffer, 
                "Angular momentum of wheels": self.angular_threshold,
                "magnetometer": self.magnetometer_buffer}

    def Per_Timestep(self, Data):
        for sensor in self.sensors:
            self.sensors[sensor].append(Data[sensor][0])

        self.sun_fault(self.sensors['sun'])
        self.star_fault(self.sensors['star'])
        self.earth_fault(self.sensors['earth'])
        self.angular_momentum_fault(self.sensors['Angular momentum of wheels'])


    #####################################################
    # IF THE THRESHOLD OF THE SUN VECTOR IS LARGER THAT #
    #     A SPECIFIED VALUE THEN IT RETURN AN ERROR     #
    #####################################################

    def sun_fault(self, sun):
        sun = np.array((sun))
        current_sun = sun[0]
        mean_sun = np.mean(sun)
        var_sun = np.var(sun)
        norm_sun = np.linalg.norm(current_sun)

        if norm_sun != 1:
            Error = "SUN_BROKEN"

        
    ######################################################
    # IF THE THRESHOLD OF THE STAR VECTOR IS LARGER THAT #
    #     A SPECIFIED VALUE THEN IT RETURN AN ERROR      #
    ######################################################

    def star_fault(self, star):
        star = np.array((star))
        current_star = star[0]
        mean_star = np.mean(star)
        var_star = np.var(star)
        norm_star = np.linalg.norm(current_star)

        if norm_star != 1:
            Error = "STAR_BROKEN"

    
    ########################################
    # IF THE EARTH VECTOR IS LARGER THAN A #
    # GIVEN THRESHOLD THEN RETURN AN ERROR #
    ########################################

    def earth_fault(self, earth):
        earth = np.array((earth))
        current_earth = earth[0]
        mean_earth = np.mean(earth)
        var_earth = np.var(earth)
        norm_earth = np.linalg.norm(current_earth)

        if norm_earth != 1:
            Error = "EARTH_BROKEN"

    
    ######################################################
    # IF THE ANGULAR MOMENTUM IS LARGER THAN A SPECIFIED #
    #    VALUE OR REMAINS LARGER THAN RETURN AN ERROR    #
    ######################################################

    def angular_momentum_fault(self, angular_moment):
        angular = np.array((angular_moment))
        current_angular = angular[0]
        mean_angular = np.mean(angular)
        var_angular = np.var(angular)
        norm_angular = np.linalg.norm(current_angular)

        if norm_angular != 1:
            Error = "ANGULAR_BROKEN"

    
    ########################################################
    # IF THE MAGNETOMETER IS LARGER THAN A SPECIFIED VALUE #
    #        OR REMAINS LARGER THAN RETURN AN ERROR        #
    ########################################################

    def magnetometer_fault(self, magnetometer):
        magnetometer = np.array((magnetometer))
        current_magnetometer = magnetometer[0]
        mean_magnetometer = np.mean(magnetometer)
        var_magnetometer = np.var(magnetometer)
        norm_magnetometer = np.linalg.norm(current_magnetometer)

        if norm_magnetometer != 1:
            Error = "MAGNETOMETER_BROKEN"


#####################################################################
#   THIS CLASS IS FOR MORE COMPLEX OPERATIONS THAN IF STATEMENTS    #
# SUCH AS CORRELATION WITH THE MATRIX AND OTHER STATISTICAL METHODS #
#####################################################################

class Correlation_detection:
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


###################################################
# THIS CLASS IS THE MOST NOVEL DETECTION METHODS  #
# THIS CLASS WILL HAVE MULTIPLE METHODS TO CHOOSE #
# AND THESE METHODS WILL BE COMPARED TO SEE WHICH #
#       WORKS METHODS BEST FOR EACH SENSOR        #
###################################################

class Encompassing_detection:
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