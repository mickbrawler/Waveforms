import json
import numpy as np
from matplotlib import pyplot as plt
import math
import time
import argparse

def test_plotter(T, N, var1=0, var2=0, stat=2, bg_test=True, plot1="test_plot1", plot2="test_plot2"):

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

    # Calculates the index of the attribute variable
    # (A, F, or G) which the heat map should marginalize
    # over
    if ((var1 <= 2) and (var2 <= 2)):    #checks if range valid
        if ((0<= var1) and (0<= var2)):  #checks if range valid
            if (var2 < var1):            #properly orders
                var1, var2 = var2, var1
                index=(1/4)*(var1-1)*(var1-2)*(var2)*(7-3*var2) #maps to needed index

            elif (var2==var1):
                index=var2  # if the values are the same, set as index

            else:
                index=(1/4)*(var1-1)*(var1-2)*(var2)*(7-3*var2) #maps to needed index

        else:
            print("improper index inputs: both must be >= 0")
            index=0
    else:
        print("improper index inputs: both must be <= 2")
        index=0

    mapping={}   #initalizes mapping for heat map

    holder=np.full((A_LEN, F_LEN, G_LEN), 0, dtype=list)

    # builds link between keys and np-array 
    for j in range(F_LEN):
        for i in range(A_LEN):
            for k in range(G_LEN):
                F=F_RANGE[j]
                A=A_RANGE[i]
                G=G_RANGE[k]

                holder[(i,j,k)]=[0]
                mapping.update({ (A,F,G) : holder[(i,j,k)]})

    #initalizes threshold axis points, counts for each threshold
    #and space of parameter space
    THRS_AXIS=[]
    COUNT_AXIS=[]
    PSPACE_LEN=len(AFG_PAIR)
    PSACE_PAIRS=[]

    #sets current threshold value
    # how many templates succeeded per threshold (succ_count_thrhld), how many trials succeeded per threshold (heatcount)
    for thrshld in np.linspace(0,T,N):

        THRS_AXIS.append(thrshld)
        succ_count_thrhld=0

        #loops from trial/ parameter space pairs 
        for i in range(trials):
            heat_count=0
            for j in range(PSPACE_LEN):
                
                RM_ij=np.array(RHO_MOD[str(i)][0][stat][j])
                FG_ij=RHO_MOD[str(i)][1][stat][j]
                N_BG_ij=len(RM_ij[ RM_ij > FG_ij])

                #tests if the given pair passes the threshold test
                if (bg_test==True):
                    if ((N_BG_ij==0) and (FG_ij > thrshld)): # (per template) offsource peaks < onsource peak, onsource peak > threshold

                        if heat_count==0:
                            succ_count_thrhld+=1
                            mapping[tuple(output[str(i)][0][:3])][0]+=1 # adds 1 to that paramater combo's value, default is 0 
                            heat_count+=1
                else:
                    if (FG_ij > thrshld): # (per template) onsource peak > threshold

                        if heat_count==0:
                            succ_count_thrhld+=1
                            mapping[tuple(output[str(i)][0][:3])][0]+=1 # adds 1 to that paramater combo's value, default is 0 
                            heat_count+=1

        COUNT_AXIS.append(succ_count_thrhld)

    plt.plot(THRS_AXIS,COUNT_AXIS)
    plt.xlabel("$Threshold$")
    plt.ylabel("Counts")
    plt.show()
    plt.savefig("plots/{}.png".format(plot1))
    plt.clf()

    #redefines things to be normal arrays inside
    for j in range(F_LEN):
        for i in range(A_LEN):
            for k in range(G_LEN):
                holder[(i,j,k)]=holder[(i,j,k)][0] # I think [0] is making the list value for each parameter combo just its value

    #marginalizes array in direction of index
    index=int(index) # index is still weird
    w=holder.sum(index)

    #builds a copy of the original array, but fixes things to be integers 
    cop=np.full(w.shape,0) # array in shape of marginalized array of holder values
    z=[(i,j) for i in range(w.shape[0]) for j in range(w.shape[1])] 
    for tup in z: # z holds the index for cop, these indexes are each of two variables in use
        cop[tup]=int(w[tup])

    label=np.array(["amplitude","frequency", "gamma"])
    x=np.array([0, 1, 2])
    x=x[ x != index] # parameter index that aren't index are used

    plt.imshow( cop, cmap=plt.cm.hot)
    plt.xlabel(label[x[1]])
    plt.ylabel(label[x[0]])
    plt.savefig("plots/{}.png".format(plot2))

    #plt.imshow( heat_array, cmap=plt.cm.hot) 

#def test_plotter(T, N, var1=0, var2=0, stat=2, bg_test=True, plot1="test_plot1", plot2="test_plot2"):

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--T', type=int)
    parser.add_argument('--N', type=int)
    parser.add_argument('--var1', nargs='?', const=1, type=int, default=0)
    parser.add_argument('--var2', nargs='?', const=1, type=int, default=0)
    parser.add_argument('--stat', nargs='?', const=1, type=int, default=2)
    parser.add_argument('--bg_test', nargs='?', const=1, type=bool, default=True)
    parser.add_argument('--plot1', nargs='?', const=1, type=str, default="test_plot1")
    parser.add_argument('--plot2', nargs='?', const=1, type=str, default="test_plot2")

    args = parser.parse_args()

    test_plotter(args.T,args.N,args.var1,args.var2,args.stat,args.bg_test,args.plot1,args.plot2)
