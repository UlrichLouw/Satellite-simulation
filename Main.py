import sys
import numpy as np
from Simulation.Parameters import SET_PARAMS
import pandas as pd
import multiprocessing
from pathlib import Path
from Simulation.dynamics import Single_Satellite
from Simulation.Save_display import visualize_data, save_as_csv, save_as_pickle, save_as_excel
import Fault_prediction.Fault_detection as Fault_detection
import Simulation.Constellation as Constellation

# ! The matplotlib cannot display plots while visual simulation runs.
# ! Consequently the Display and visualize parameters in Parameters 
# ! must be set as desired

if SET_PARAMS.Display:
    import Simulation.Satellite_display as view

#####################################
# LOOP THROUGH DYNAMICS IF MULTIPLE #
#       THREADS ARE REQUIRED        #
#####################################

def loop(index, D, Data, orbit_descriptions):
    print(SET_PARAMS.Fault_names_values[index])

    if SET_PARAMS.Display:
        satellite = view.initializeCube(SET_PARAMS.Dimensions)
        pv = view.ProjectionViewer(1920, 1080, satellite)
    
    for j in range(1, int(SET_PARAMS.Number_of_orbits*SET_PARAMS.Period/(SET_PARAMS.faster_than_control*SET_PARAMS.Ts)+1)):
        w, q, A, r, sun_in_view = D.rotation()
        if SET_PARAMS.Display and j%SET_PARAMS.skip == 0:
            pv.run(w, q, A, r, sun_in_view)
        
        if j%(int(SET_PARAMS.Number_of_orbits*SET_PARAMS.Period/(SET_PARAMS.faster_than_control*SET_PARAMS.Ts)/10)) == 0:
            print("Number of time steps for orbit loop number", index, " = ", "%.2f" % float(j/int(SET_PARAMS.Number_of_orbits*SET_PARAMS.Period/(SET_PARAMS.faster_than_control*SET_PARAMS.Ts))))

        if SET_PARAMS.Fault_simulation_mode == 2 and j%(int(SET_PARAMS.Number_of_orbits*SET_PARAMS.Period/(SET_PARAMS.faster_than_control*SET_PARAMS.Ts)/SET_PARAMS.fixed_orbit_failure)) == 0:
            D.initiate_purposed_fault(SET_PARAMS.Fault_names_values[index])
            if SET_PARAMS.Display:
                pv.fault = D.fault

    if SET_PARAMS.Visualize and SET_PARAMS.Display == False:
        path_to_folder = Path("Plots/" + str(D.fault))
        path_to_folder.mkdir(exist_ok=True)
        visualize_data(D.Orbit_Data, D.fault)
    
    elif SET_PARAMS.Display == True:
        pv.save_plot(D.fault)

    print("Number of multiple orbits", index)  
    Data[index] = D.Orbit_Data
    orbit_descriptions[index] = D.fault

    if SET_PARAMS.save_as == ".csv":
        save_as_csv(D.Orbit_Data, index)
    else:
        save_as_pickle(D.Orbit_Data, index)

################################################################
# FOR ALL OF THE FAULTS RUN A NUMBER OF ORBITS TO COLLECT DATA #
################################################################
if __name__ == "__main__":
    #########################################################
    # IF THE SAVE AS IS EQUAL TO XLSX, THE THREADING CANNOT #
    #           BE USED TO SAVE CSV FILES                   #     
    #########################################################
    SET_PARAMS.Display = True
    SET_PARAMS.save_as = ".xlsx"
    SET_PARAMS.Kalman_filter_use = "EKF"
    SET_PARAMS.Number_of_orbits = 10
    SET_PARAMS.skip = 200
    SET_PARAMS.Number_of_satellites = 100
    SET_PARAMS.Constellation = False
    SET_PARAMS.k_nearest_satellites = 5


    #####################################
    # PARAMETERS FOR SATELLITE DYNAMICS #
    #####################################
    s_list, t_list, J_t, fr = SET_PARAMS.s_list, SET_PARAMS.t_list, SET_PARAMS.J_t, SET_PARAMS.fr
    #########################################################
    #   TO ENABLE A CONSTELLATION A CLASS IS CREATED THAT   #
    #     CONTAINS THE DATA OF THE ENTIRE CONSTELLATION     #
    #  THAT DATA IS TRANSFERED TO EACH SATELLITE DEPENDING  #
    # ON THE SATELLITES ID AND THE SATELLITES CLOSEST TO IT #
    #########################################################

    if SET_PARAMS.Constellation:
        Stellar = Constellation.Constellation(SET_PARAMS.Number_of_satellites)

        for sat_num in range(SET_PARAMS.Number_of_satellites):
            Stellar.initiate_satellite(sat_num)

        for j in range(int(SET_PARAMS.Number_of_orbits*SET_PARAMS.Period/(SET_PARAMS.faster_than_control*SET_PARAMS.Ts)+1)):
            for sat_num in range(SET_PARAMS.Number_of_satellites):
                Stellar.satellites[sat_num].run()

    elif SET_PARAMS.save_as == ".xlsx":
        FD = Fault_detection.Basic_detection()
        Data = []
        orbit_descriptions = []
        for i in range(SET_PARAMS.Number_of_multiple_orbits):
            D = Single_Satellite(i, s_list, t_list, J_t, fr)

            print(SET_PARAMS.Fault_names_values[i+1])

            if SET_PARAMS.Display:
                satellite = view.initializeCube(SET_PARAMS.Dimensions)
                pv = view.ProjectionViewer(1920, 1080, satellite)
            
            for j in range(int(SET_PARAMS.Number_of_orbits*SET_PARAMS.Period/(SET_PARAMS.faster_than_control*SET_PARAMS.Ts)+1)):
                w, q, A, r, sun_in_view = D.rotation()

                # Detect faults based on data from Dynamics (D):
                Fault = FD.Per_Timestep(D.Orbit_Data)
                if Fault != "None":
                    print(Fault)

                if SET_PARAMS.Display and j%SET_PARAMS.skip == 0:
                    pv.run(w, q, A, r, sun_in_view)
                
                if j%(int(SET_PARAMS.Number_of_orbits*SET_PARAMS.Period/(SET_PARAMS.faster_than_control*SET_PARAMS.Ts)/10)) == 0:
                    print("Number of time steps for orbit loop number", i, " = ", "%.2f" % float(j/int(SET_PARAMS.Number_of_orbits*SET_PARAMS.Period/(SET_PARAMS.faster_than_control*SET_PARAMS.Ts))))

                if SET_PARAMS.Fault_simulation_mode == 2 and (j+1)%(int(SET_PARAMS.Number_of_orbits*SET_PARAMS.Period/(SET_PARAMS.faster_than_control*SET_PARAMS.Ts)/SET_PARAMS.fixed_orbit_failure)) == 0:
                    D.initiate_purposed_fault(SET_PARAMS.Fault_names_values[i+1])
                    if SET_PARAMS.Display:
                        pv.fault = D.fault

            if SET_PARAMS.Visualize and SET_PARAMS.Display == False:
                path_to_folder = Path("Plots/" + str(D.fault))
                path_to_folder.mkdir(exist_ok=True)
                visualize_data(D.Orbit_Data, D.fault)
            
            elif SET_PARAMS.Display == True:
                pv.save_plot(D.fault)

            Data.append(D.Orbit_Data)
            orbit_descriptions.append(str(D.fault))

        save_as_excel(Data, orbit_descriptions)

    ######################################################
    # IF THE SAVE AS IS NOT EQUAL TO XLSX, THE THREADING #
    #           CAN BE USED TO SAVE CSV FILES            #
    ######################################################
    else:
        threads = []

        manager = multiprocessing.Manager()
        Data = manager.dict()
        orbit_descriptions = manager.dict()

        for i in range(1, SET_PARAMS.Number_of_multiple_orbits+1):
            D = Single_Satellite(i, s_list, t_list, J_t, fr)

            t = multiprocessing.Process(target=loop, args=(i, D, Data, orbit_descriptions))
            threads.append(t)
            t.start()
            print("Beginning of", i)
            if i%15 == 0 or i == SET_PARAMS.Number_of_multiple_orbits:
                for process in threads:     
                    process.join()

                threads = []
