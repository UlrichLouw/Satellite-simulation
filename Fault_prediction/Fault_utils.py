import pandas as pd
import numpy as np
from Parameters import SET_PARAMS
import collections

X_buffer = []
Y_buffer = []

excel_file = SET_PARAMS.excel_filename
xls = pd.ExcelFile(excel_file)
def Binary_split(classified_data):
    for index in range(len(classified_data.index)):
        if classified_data["Current fault"][index] == "None":
            classified_data["Current fault"][index] = 0
        else:
            classified_data["Current fault"][index] = 1

    return classified_data

def Dataset_order(index, binary_set, buffer, categorical_num, use_previously_saved_models = False, columns_compare = None, columns_compare_to = None):
    X_buffer_replaced = []
    if SET_PARAMS.Save_excel_file == True:
        Data = pd.read_excel(xls, str(index))

    else:
        pickle_file = SET_PARAMS.pickle_filename
        Data = pd.read_pickle(pickle_file)
    
    if binary_set and use_previously_saved_models == False:
        Orbit = Data.drop(columns = ['Current fault', 'Current fault numeric'])
    elif categorical_num == True:
        Orbit = Data.drop(columns = ['Current fault', 'Current fault binary'])
    else:
        Orbit = Binary_split(Data)

    if columns_compare != None:
        columns_to_keep = columns_compare + columns_compare_to
        for i in columns_to_keep:
            temp = Orbit[i]
        
        Orbit = Orbit[columns_to_keep]
        X = Orbit[columns_compare]
        Y = Orbit[columns_compare_to]
    else:
        Orbit.drop(columns = ['Sun in view'], inplace = True)
        X = Orbit.iloc[:,0:-1].values
        X_correlation_sun_earth_magnetometer = Orbit.iloc[:,0:9].values
        Y = Orbit.iloc[:,-1].values

    buffer_x = collections.deque(maxlen = SET_PARAMS.buffer_size)
    buffer_correlation_sun_earth_magnetometer = collections.deque(maxlen = SET_PARAMS.buffer_size)
    y = Y[SET_PARAMS.buffer_size - 1:]
    buffer_y = []

    if buffer == True:
        for i in range(SET_PARAMS.buffer_size - 1):
            buffer_x.append(X[i,:])
            buffer_correlation_sun_earth_magnetometer.append(X_correlation_sun_earth_magnetometer[i,:])

        for i in range(SET_PARAMS.buffer_size, X.shape[0]):
            buffer_x.append(X[i,:])
            buffer_correlation_sun_earth_magnetometer.append(X_correlation_sun_earth_magnetometer[i,:])
            if use_previously_saved_models == True:
                buffer_y.append(np.fromstring(y[i-SET_PARAMS.buffer_size][1:-1], dtype = float, sep=','))
            #Binary_stat_fault(buffer_correlation_sun_earth_magnetometer)
            X_buffer.append(np.asarray(buffer_x).flatten())
            X_buffer_replaced.append(np.asarray(X_buffer).flatten())

    X = np.asarray(X_buffer_replaced)
    if use_previously_saved_models == True:
        Y = np.asarray(buffer_y)
        Y = Y.reshape(X.shape[0], Y.shape[1])
        Y_buffer.append(Y)
    else:
        Y = np.asarray(Y[SET_PARAMS.buffer_size:]).reshape(X.shape[0],1)
        Y_buffer.append(Y)

    return Y, Y_buffer, X, X_buffer, Orbit