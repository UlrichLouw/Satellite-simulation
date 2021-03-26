import pandas as pd
from Parameters import SET_PARAMS
import numpy as np
import collections
import itertools

MIN_CORRELATION = 0.6
MAX_VARIANCE = 10

def Binary_stat_fault(Data):
    Fault = False
    var = 0

    corr = np.min(np.corrcoef(Data))
    #var = np.var(Data)

    if corr < MIN_CORRELATION:
        Fault = True
        print(Fault)
    """
    elif var > MAX_VARIANCE:
        Fault = True
    """
    return Fault