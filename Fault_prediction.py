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

sc = StandardScaler()

excel_file = SET_PARAMS.excel_filename
xls = pd.ExcelFile(excel_file)
RANDOM_SEED = 0

def Binary_split_with_buffer(classified_data):
    for index in range(len(classified_data.index)):
        if classified_data["Current fault"][index] == "None":
            classified_data["Current fault"][index] = 0
        else:
            classified_data["Current fault"][index] = 1

    return classified_data

# A utility method to create a tf.data dataset from a Pandas Dataframe
def df_to_dataset(dataframe, shuffle=True, batch_size=32):
  dataframe = dataframe.copy()
  labels = dataframe.pop('Current fault')
  ds = tf.data.Dataset.from_tensor_slices((dict(dataframe), labels))
  if shuffle:
    ds = ds.shuffle(buffer_size=len(dataframe))
  ds = ds.batch(batch_size)
  return ds

if __name__ == "__main__":
    confusion_matrices = []
    All_orbits = []
    buffer = True
    binary_set = True
    X_buffer = []
    Y_buffer = []
    
    for index in SET_PARAMS.Fault_names:
        for direction in SET_PARAMS.Fault_names[index]:
            if SET_PARAMS.Save_excel_file == True:
                Data = pd.read_excel(xls, str(index) + str(direction))

            else:
                pickle_file = SET_PARAMS.pickle_filename
                Data = pd.read_pickle(pickle_file)
            
            if binary_set == True:
                Orbit = Data.drop(columns = ['Current fault'], inplace = True)
            else:
                Orbit = Binary_split_with_buffer(Data)

            Orbit.drop(columns = ['Sun in view'], inplace = True)
            X = Orbit.iloc[:,1:-1].values
            Y = Orbit.iloc[:,-1].values

            buffer_x = collections.deque(maxlen = SET_PARAMS.buffer_size)
            buffer_y = Y[SET_PARAMS.buffer_size:]

            for i in range(SET_PARAMS.buffer_size - 1):
                buffer_x.append(X[i,:])

            for i in range(SET_PARAMS.buffer_size, X.shape[0]):
                buffer_x.append(X[i,:])
                X_buffer.append(np.asarray(buffer_x).flatten())
            
            Y_buffer.append(buffer_y)
            All_orbits.append(Orbit)
    
    if buffer == False:
        All_orbits = pd.concat(All_orbits)
        X = All_orbits.iloc[:,1:-1].values
        Y = All_orbits.iloc[:,-1].values
    else:
        X = np.asarray(X_buffer)
        Y = np.asarray(Y_buffer).reshape(X.shape[0], 1)

    X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size = 0.2)
    X_train = sc.fit_transform(X_train)
    X_test = sc.transform(X_test)

    X_train = np.asarray(X_train).astype(np.float32)
    y_train = np.asarray(y_train).astype(np.float32)
    y_test = np.asarray(y_test).astype(np.int)

    feature_columns = []

    # numeric cols
    for header in Orbit.columns.to_list():
        feature_columns.append(feature_column.numeric_column(header))

    classifier = tf.estimator.DNNClassifier(feature_columns= feature_columns, hidden_units = [30,10], n_classes = 2)

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

    model.fit(X_train, y_train, epochs=10, batch_size = batch_size, use_multiprocessing=True)

    y_pred = model.predict(X_test)

    cm = confusion_matrix(y_test, y_pred.round())
    print(cm)