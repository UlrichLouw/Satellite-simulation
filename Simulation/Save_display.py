from Simulation.Parameters import SET_PARAMS
import pandas as pd
import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objects as go

# Function to save a csv file of simulation data
def save_as_excel(Data, sheetnames):
    with pd.ExcelWriter(SET_PARAMS.filename + ".xlsx") as writer:
        i = 0
        for data in Data:
            df = pd.DataFrame(data, columns = data.keys())
            sheetname = sheetnames[i]
            df.to_excel(writer, sheet_name = sheetname, index = False)
            i += 1

####################
# SAVE AS CSV FILE #
####################

def save_as_csv(Data, orbit):
    df = pd.DataFrame(Data, columns = Data.keys())
    df.to_csv(SET_PARAMS.filename + str(orbit) + ".csv")

#######################################################
# FUNCTION TO SAVE A PICKLE FILE OF SIMULATION DATA   #
#######################################################
def save_as_pickle(Data, orbit):
    df = pd.DataFrame(Data, columns = Data.keys())
    df.to_pickle(SET_PARAMS.filename + str(orbit) + ".pkl")

##########################################

##########################################
# FUNCTION TO VISUALIZE DATA AS GRAPHS   #
##########################################
def visualize_data(D, fault):
    for i in D:
        if i == "Current fault" or i == "Current fault binary" or i == "Current fault numeric":
            pass
        elif i == "Sun in view":
            pass
        else:
            y = np.array((D[i]))
            fig = make_subplots(rows=3, cols=1)
            x = y.shape[0]
            x = np.arange(0,x,1)
            y_min = np.amin(y)
            y_max = np.amax(y)

            fig.append_trace(go.Scatter(
                x=x,
                y=y[:,0],
                name = "x"
            ), row=1, col=1)

            fig.append_trace(go.Scatter(
                x=x,
                y=y[:,1],
                name = "y"
            ), row=2, col=1)

            fig.append_trace(go.Scatter(
                x=x,
                y=y[:,2],
                name = 'z'
            ), row=3, col=1)

            fig.update_yaxes(range=[y_min, y_max], row=1, col=1)
            fig.update_yaxes(range=[y_min, y_max], row=2, col=1)
            fig.update_yaxes(range=[y_min, y_max], row=3, col=1)
            fig.update_layout(height=600, width=600, title_text=str(i))
            fig.write_html("Plots/" + str(fault) +"/"+ str(i)+".html")