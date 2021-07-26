import json
import numpy as np
from matplotlib import pyplot as plt
import math
import time
import argparse

def ROC_Data(T0, Tf, N, stat):

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
    
    with open("Merged_jsons/Merged_Peaks.json", "r") as f:
        RHO_MOD = json.load(f)

    PSPACE_LEN=len(AFG_PAIR)

    # Stats per threshold
    Detection_Prob = []
    New_False_Prob = []

    #sets current threshold value
    for thrshld in np.linspace(T0,Tf,N):

        # Detection/ False Alarm probability counters
        Detect_count = 0
        False_count = 0

        #loops from trial/ parameter space pairs 
        for i in range(trials):

            # Detection Probability
            Max_FG_ij = max(RHO_MOD[str(i)][1][stat]) # max of onsources per trial
            if Max_FG_ij > thrshld:
                Detect_count += 1

            for j in range(PSPACE_LEN):

                RM_ij=np.array(RHO_MOD[str(i)][0][stat][j])

                # False Alarm probability
                falses_ij = len(RM_ij[RM_ij > thrshld])
                False_count += falses_ij

        # Detection/False Alarm probability stats
        Detect_stat = Detect_count / trials
        False_stat = False_count / (len(RHO_MOD[str(i)][0][stat][0]) * PSPACE_LEN * trials)

        # Appending stat per threshold
        Detection_Prob.append(Detect_stat)
        New_False_Prob.append(False_stat)
        
    return(Detection_Prob, New_False_Prob)

def ROC_Curve(N, outputfile="ROC_test"):
    
    with open("run_uniques/essentials.json", "r") as f:
        essentials = json.load(f)
    stat_list = essentials["essentials"][8]

    
    with open("Merged_jsons/Merged_thresholds.json", "r") as f:
        Thresholds = json.load(f)
    
    tempn = len(Thresholds)
    
    for s in range(tempn):
        
        thresholds = Thresholds[str(s)]

        Detection_Prob, New_False_Prob = ROC_Data(min(thresholds), max(thresholds), N, s)
        plt.plot(New_False_Prob, Detection_Prob, label=stat_list[s])

    plt.xlabel("New_False_Probs")
    plt.ylabel("Detection_Probs")
    plt.title("ROC Curve:N={}".format(N))
    plt.legend()
    plt.savefig("plots/{}.png".format(outputfile))
    plt.show()

# def ROC_Curve(N, outputfile="ROC_test"):

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--N', type=int)
    parser.add_argument('--outputfile', nargs='?', const=1, type=str, default="ROC_test")
    args = parser.parse_args()

    ROC_Curve(args.N,args.outputfile)
