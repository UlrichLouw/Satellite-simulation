import pandas as pd
from Parameters import SET_PARAMS
import numpy as np
import tensorflow as tf
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from tensorflow import feature_column
from tensorflow.keras import layers
from sklearn.metrics import confusion_matrix
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
    for index in SET_PARAMS.Fault_names:
        for direction in SET_PARAMS.Fault_names[index]:
            if SET_PARAMS.Save_excel_file == True:
                Data = pd.read_excel(xls, str(index) + str(direction))

            else:
                pickle_file = SET_PARAMS.pickle_filename
                Data = pd.read_pickle(pickle_file)
            
            Orbit = Binary_split_with_buffer(Data)
            Orbit.drop(columns = ['Sun in view'], inplace = True)
            All_orbits.append(Orbit)
    
    All_orbits = pd.concat(All_orbits)
    X = All_orbits.iloc[:,1:-1].values
    Y = All_orbits.iloc[:,-1].values
        
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
            tf.keras.layers.Dense(units=12, activation = 'relu'),
            tf.keras.layers.Dense(units=128, activation='relu'),
            tf.keras.layers.Dropout(rate=0.2),
            tf.keras.layers.Dense(units=128, activation='relu'),
            tf.keras.layers.Dense(units=1, activation='sigmoid')
            ])

    model.compile(optimizer='adam',
            loss='binary_crossentropy',
            metrics=['accuracy'])

    batch_size = 32 # A small batch sized is used for demonstration purposes

    model.fit(X_train, y_train, epochs=10, batch_size = batch_size, use_multiprocessing=True)

    y_pred = model.predict(X_test)

    cm = confusion_matrix(y_test, y_pred.round())
    print(cm)