import json
import numpy as np
from matplotlib import pyplot as plt
import math
import time
import argparse

def Scatter_plotter(thrshld, xvar, yvar, stat=2, plot="Scatter_plot"): 

    with open("run_uniques/essentials.json", "r") as f:
        essentials = json.load(f)
    A_LEN = essentials["essentials"][0]
    F_LEN = essentials["essentials"][1]
    G_LEN = essentials["essentials"][2]
    F_RANGE = essentials["essentials"][3]
    A_RANGE = essentials["essentials"][4]
    G_RANGE = essentials["essentials"][5]
    AFG_PAIR = essentials["essentials"][6]
    trials = essentials["essentials"][7]

    with open("Merged_jsons/Merged_output.json", "r") as f:
        output = json.load(f)

    with open("Merged_jsons/Merged_Peaks.json", "r") as f:
        RHO_MOD = json.load(f)

    label=["amplitude","frequency", "gamma"]
    PSPACE_LEN=len(AFG_PAIR)
    SUCC_PAIRS=([],[])
    FAIL_PAIRS=([],[])

    #loops from trial/ parameter space pairs 
    for i in range(trials):

        succ_count_thrhld=0  #test condition for adding to fail array

        for j in range(PSPACE_LEN):

            RM_ij=np.array(RHO_MOD[str(i)][0][stat][j])  #get moded rho_ij background values
            FG_ij=RHO_MOD[str(i)][1][stat][j]            #get moded rho_ij foreground value
            N_BG_ij=len(RM_ij[ RM_ij > FG_ij])     # get background vals > froeground vals

            #tests if the given pair passes the threshold test
            if ((N_BG_ij==0) and (FG_ij > thrshld)):

                parameter_x=output[str(i)][0][xvar]  #get parameter 'xvar' of trial i
                parameter_y=output[str(i)][0][yvar]  #get parameter 'yvar' of trial i
                SUCC_PAIRS[0].append(parameter_x) #add to successes
                SUCC_PAIRS[1].append(parameter_y)
                succ_count_thrhld+=1
                break

        #sees if trial i failed the test; adds parameters to fail array
        if succ_count_thrhld==0:
            parameter_x=output[str(i)][0][xvar]
            parameter_y=output[str(i)][0][yvar]
            FAIL_PAIRS[0].append(parameter_x)
            FAIL_PAIRS[1].append(parameter_y)

    colors = ("red", "blue")
    groups = ("pass", "fail")
    marks = ("o", "*")
    data = ( SUCC_PAIRS, FAIL_PAIRS )

    # Create plot
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    for data, marks, color, group in zip(data, marks, colors, groups):
        x , y = data
        ax.scatter(x, y, marker=marks, c=color, edgecolors='none', s=30, label=group)

    plt.title('found-missed plot')
    plt.xlabel(label[xvar])
    plt.ylabel(label[yvar])
    plt.legend(loc=2)
    plt.show()
    plt.savefig("plots/{}.png".format(plot))

# def Scatter_plotter(thrshld, xvar, yvar, stat=2, plot="Scatter_plot"): 

if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--thrshld', type=int)
    parser.add_argument('--xvar', type=int)
    parser.add_argument('--yvar', type=int)
    parser.add_argument('--stat', nargs='?', const=1, type=int, default=2)
    parser.add_argument('--plot', nargs='?', const=1, type=str, default="Scatter_plot")
    args = parser.parse_args()

    Scatter_plotter(args.thrshld,args.xvar,args.yvar,args.stat,args.plot)
