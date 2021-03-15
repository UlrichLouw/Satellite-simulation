import pandas as pd
from Parameters import SET_PARAMS
import numpy as np
import tensorflow as tf
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from tensorflow import feature_column
from tensorflow.keras import layers
from sklearn.metrics import confusion_matrix
import collections
from tensorflow.keras.models import model_from_json
import os

sc = StandardScaler()

excel_file = SET_PARAMS.excel_filename
xls = pd.ExcelFile(excel_file)
RANDOM_SEED = 0

Fault_names_to_num = {
    "NoneNone": 1,
    "Sun sensorx": 2,
    "Sun sensory": 3, 
    "Sun sensorz": 4,
    "Magnetometerx": 5, 
    "Magnetometery": 6, 
    "Magnetometerz": 7,
    "Earth sensorx": 8, 
    "Earth sensory": 9, 
    "Earth sensorz": 10,
    "Reaction wheelx": 11, 
    "Reaction wheely": 12, 
    "Reaction wheelz": 13,
    "Controlall": 14
}

loaded_model_1 = None
loaded_model_2 = None
loaded_model_3 = None
loaded_model_4 = None
loaded_model_5 = None
loaded_model_6 = None
loaded_model_7 = None
loaded_model_8 = None
loaded_model_9 = None
loaded_model_10 = None
loaded_model_11 = None
loaded_model_12 = None
loaded_model_13 = None
loaded_model_14 = None


model_names = {
    1: loaded_model_1,
    2: loaded_model_2,
    3: loaded_model_3,
    4: loaded_model_4,
    5: loaded_model_5,
    6: loaded_model_6,
    7: loaded_model_7,
    8: loaded_model_8,
    9: loaded_model_9,
    10: loaded_model_10,
    11: loaded_model_11,
    12: loaded_model_12,
    13: loaded_model_13,
    14: loaded_model_14
}

model_data_lists = {
    1: [],
    2: [],
    3: [],
    4: [],
    5: [],
    6: [],
    7: [],
    8: [],
    9: [],
    10: [],
    11: [],
    12: [],
    13: [],
    14: []
}

testing_data_lists = {
    1: [],
    2: [],
    3: [],
    4: [],
    5: [],
    6: [],
    7: [],
    8: [],
    9: [],
    10: [],
    11: [],
    12: [],
    13: [],
    14: []
}

def Binary_split(classified_data):
    for index in range(len(classified_data.index)):
        if classified_data["Current fault"][index] == "None":
            classified_data["Current fault"][index] = 0
        else:
            classified_data["Current fault"][index] = 1

    return classified_data

"""
def numerical_split(classified_data):
    zeros = np.zeros((14,), dtype = int)
    for index in range(len(classified_data.index)):
        zeros[int(Fault_names_to_num[classified_data["Current fault"][index]]) - 1] = 1
        classified_data["Current fault numeric"][index] = [zeros]
        zeros[int(Fault_names_to_num[classified_data["Current fault"][index]]) - 1] = 0

    return classified_data
"""

# A utility method to create a tf.data dataset from a Pandas Dataframe
def df_to_dataset(dataframe, shuffle=True, batch_size=32):
  dataframe = dataframe.copy()
  labels = dataframe.pop('Current fault')
  ds = tf.data.Dataset.from_tensor_slices((dict(dataframe), labels))
  if shuffle:
    ds = ds.shuffle(buffer_size=len(dataframe))
  ds = ds.batch(batch_size)
  return ds

def prediction_NN(X, Y, index, direction):
    X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size = 0.2)
    X_train = np.asarray(sc.fit_transform(X_train)).astype(np.float32)
    X_test = np.asarray(sc.transform(X_test)).astype(np.float32)

    y_train = np.asarray(y_train).astype(np.float32)
    y_test = np.asarray(y_test).astype(np.int)

    model = tf.keras.models.Sequential([
            tf.keras.layers.Dense(units=X.shape[1], activation = 'relu'),
            tf.keras.layers.Dense(units=128, activation='relu'),
            tf.keras.layers.Dropout(rate=0.2),
            tf.keras.layers.Dense(units=128, activation='relu'),
            tf.keras.layers.Dense(units=1, activation='sigmoid')
            ])

    model.compile(optimizer='adam',
            loss='binary_crossentropy',
            metrics=['Precision'])

    batch_size = 32 # A small batch sized is used for demonstration purposes

    model.fit(X_train, y_train, epochs=10, batch_size = batch_size, use_multiprocessing=True, verbose=1)

    y_pred = model.predict(X_test)

    model_json = model.to_json()
    with open("models/" + str(index) + str(direction) + ".json", "w") as json_file:
        json_file.write(model_json)

    model.save_weights("models/" + str(index) + str(direction) + ".h5")

    cm = confusion_matrix(y_test, y_pred.round())
    return cm

def prediction_NN_determine_other_NN(X, Y):
    X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size = 0.2)
    X_train = np.asarray(sc.fit_transform(X_train)).astype(np.float32)
    X_test = np.asarray(sc.transform(X_test)).astype(np.float32)

    y_train = np.asarray(y_train).astype(np.float32)
    y_test = np.asarray(y_test).astype(np.int)

    model = tf.keras.models.Sequential([
            tf.keras.layers.Dense(units=X.shape[1], activation = 'relu'),
            tf.keras.layers.Dense(units=128, activation='relu'),
            tf.keras.layers.Dropout(rate=0.2),
            tf.keras.layers.Dense(units=128, activation='relu'),
            tf.keras.layers.Dense(units=Y.shape[1], activation='softmax')
            ])

    model.compile(optimizer='adam',
            loss='categorical_crossentropy',
            metrics=['accuracy'])

    batch_size = 32 # A small batch sized is used for demonstration purposes

    model.fit(X_train, y_train, epochs=10, batch_size = batch_size, use_multiprocessing=True, verbose=1)

    y_pred = model.predict(X_test)

    ind = 1
    for index in SET_PARAMS.Fault_names:
        for direction in SET_PARAMS.Fault_names[index]:
            json_file = open("models/" + str(index) + str(direction) + ".json", 'r')
            loaded_model_json = json_file.read()
            json_file.close()
            model_names[ind] = model_from_json(loaded_model_json)
            model_names[ind].load_weights("models/" + str(index) + str(direction) + ".h5")
            model_names[ind].compile(optimizer='adam',
            loss='binary_crossentropy',
            metrics=['Precision'])
            ind += 1

    ind = []

    for i in range(y_pred.shape[0]):
        model_data_lists[np.argmax(y_pred[i])+1].append(X_test[i,:])
        testing_data_lists[np.argmax(y_pred[i])+1].append(0 if y_test[i][0] == 1 else 1)
        ind.append(np.argmax(y_pred[i])+1)

    for i in range(1,y_pred.shape[1]+1):
        res = True in (ele == i for ele in ind)
        if res:
            y_predicted = model_names[i].predict(np.asarray(model_data_lists[i]).astype(np.float32))
            cm = confusion_matrix(np.asarray(testing_data_lists[i]), np.asarray(y_predicted).round())
            print(cm, i)

    return cm

if __name__ == "__main__":
    confusion_matrices = []
    All_orbits = []
    buffer = True
    binary_set = True
    X_buffer = []
    Y_buffer = []
    ind = 0
    use_previously_saved_models = False
    categorical_num = True
    
    for index in SET_PARAMS.Fault_names:
        for direction in SET_PARAMS.Fault_names[index]:
            X_buffer_replaced = []
            if SET_PARAMS.Save_excel_file == True:
                Data = pd.read_excel(xls, str(index) + str(direction))

            else:
                pickle_file = SET_PARAMS.pickle_filename
                Data = pd.read_pickle(pickle_file)
            
            if binary_set == True and use_previously_saved_models == False:
                Orbit = Data.drop(columns = ['Current fault', 'Current fault numeric'])
            elif categorical_num == True:
                Orbit = Data.drop(columns = ['Current fault', 'Current fault binary'])
            else:
                Orbit = Binary_split(Data)

            Orbit.drop(columns = ['Sun in view'], inplace = True)
            X = Orbit.iloc[:,1:-1].values
            Y = Orbit.iloc[:,-1].values

            buffer_x = collections.deque(maxlen = SET_PARAMS.buffer_size)
            y = Y[SET_PARAMS.buffer_size - 1:]
            buffer_y = []

            for i in range(SET_PARAMS.buffer_size - 1):
                buffer_x.append(X[i,:])

            for i in range(SET_PARAMS.buffer_size, X.shape[0]):
                buffer_x.append(X[i,:])
                if use_previously_saved_models == True:
                    buffer_y.append(np.fromstring(y[i-SET_PARAMS.buffer_size][1:-1], dtype = float, sep=','))
                X_buffer.append(np.asarray(buffer_x).flatten())
                X_buffer_replaced.append(np.asarray(buffer_x).flatten())
            
            All_orbits.append(Orbit)
            X = np.asarray(X_buffer_replaced)
            if use_previously_saved_models == True:
                Y = np.asarray(buffer_y)
                Y = Y.reshape(X.shape[0], Y.shape[1])
            else:
                Y = np.asarray(Y[SET_PARAMS.buffer_size:]).reshape(X.shape[0],1)
                Y_buffer.append(Y)

            if use_previously_saved_models == False:
                cm = prediction_NN(X, Y, index, direction)
                print(cm, str(index) + str(direction))
    
    if buffer == False:
        All_orbits = pd.concat(All_orbits)
        X = All_orbits.iloc[:,1:-1].values
        Y = All_orbits.iloc[:,-1].values
    else:
        X = np.asarray(X_buffer)
        Y = np.asarray(Y_buffer).reshape(X.shape[0], Y.shape[1])

    if use_previously_saved_models == False:
        index = "all samples"
        cm = prediction_NN(X, Y, index, None)
        print(cm, index)

    else:
        cm = prediction_NN_determine_other_NN(X, Y)
        print(cm)
    