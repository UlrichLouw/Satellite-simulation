import pandas as pd
from Parameters import SET_PARAMS
import numpy as np

pickle_filename = "Data_files/Faults.pkl"

df = pd.read_pickle(pickle_filename)
Data = pd.DataFrame.to_dict(df)

def Binary_split(classified_data):
    for index in range(len(classified_data["Current fault"])):
        if classified_data["Current fault"][index] == "None":
            classified_data["Current fault"][index] = 0
        else:
            classified_data["Current fault"][index] = 1
    return classified_data


if __name__ == "__main__":
    for index in SET_PARAMS.Fault_names:
        Orbit = Binary_split(Data["Orbit"][SET_PARAMS.Fault_names[index]])