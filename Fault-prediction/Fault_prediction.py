import pandas as pd
from Parameters import SET_PARAMS
import numpy as np

pickle_filename = "Faults.pkl"

df = pd.read_pickle(pickle_filename)
Data = pd.DataFrame.to_dict(df)

for i in SET_PARAMS.Fault_names:
    Orbit = Data["Orbit"][i]
    Orbit_fault = Data["Fault"][i]