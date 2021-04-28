import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import multivariate_normal
import random as rn
import eif as iso
import seaborn as sb
sb.set_style(style="whitegrid")
sb.set_color_codes()
import scipy.ndimage
from scipy.interpolate import griddata
import numpy.ma as ma
from numpy.random import uniform, seed

def getDepth(x, root, d):
    n = root.n
    p = root.p
    if root.ntype == 'exNode':
        return d
    else:
        if (x-p).dot(n) < 0:
            return getDepth(x,root.left,d+1)
        else:
            return getDepth(x,root.right,d+1)
        
def getVals(forest,x,sorted=True):
    theta = np.linspace(0,2*np.pi, forest.ntrees)
    r = []
    for i in range(forest.ntrees):
        temp = forest.compute_paths_single_tree(np.array([x]),i)
        r.append(temp[0])
    if sorted:
        r = np.sort(np.array(r))
    return r, theta

if __name__ == "__main__":
    confusion_matrices = []
    All_orbits = []
    X_buffer = []
    Y_buffer = []
    buffer = True
    binary_set = True
    use_previously_saved_models = False
    categorical_num = True
    
    for index in range(SET_PARAMS.Number_of_multiple_orbits):
        Y, Y_buffer, X, X_buffer, Orbit = Dataset_order(index, binary_set, buffer, categorical_num, use_previously_saved_models)
        All_orbits.append(Orbit)

        F1 = iso.iForest(X, ntrees = 200, sample_size = 256, Extension_Level=1)
        # Split each dataset into two halves: training set and test set
        train1 = Y[:int(nSamples/2)]
        train2 = X[:int(nSamples/2)]
        test1 = Y[int(nSamples/2):]
        test2 = X[int(nSamples/2):]

        # Create a cca object as an instantiation of the CCA object class. 
        cca = rcca.CCA(kernelcca = False, reg = 0., numCC = 2)

        # Use the train() method to find a CCA mapping between the two training sets.
        cca.train([train1, train2])

        # Use the validate() method to test how well the CCA mapping generalizes to the test data.
        # For each dimension in the test data, correlations between predicted and actual data are computed.
        testcorrs = cca.validate([test1, test2])
        print(testcorrs)