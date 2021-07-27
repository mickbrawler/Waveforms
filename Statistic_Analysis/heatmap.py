import json
import numpy as np
from matplotlib import pyplot as plt
import math
import time
import argparse

def index_counter(var1,var2):
    
    # Calculates the index of the attribute variable
    # (A, F, or G) which the heat map should marginaliz 
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

    return index


def heatmap(var1, var2, T, stat=2, plot1="heat_plot1", plot2 = "heat_plot2", max_OS=-1, max_BG=-1):

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
    
    with open("Merged_jsons/Merged_Max_OS.json", "r") as f:
        MAX_OS = json.load(f)

    with open("Merged_jsons/Merged_Max_BG_TEMP.json", "r") as f:
        MAX_BG_TEMP = json.load(f)
    
    index=index_counter(var1,var2)

    mapping={}   #initalizes mapping for heat map
    count_mapping={} #initalies array to count occurances of certain parameters
    BG_mapping={}

    holder=np.full((A_LEN, F_LEN, G_LEN), 0, dtype=list)
    count_holder=np.full((A_LEN, F_LEN, G_LEN), 0, dtype=list)
    BG_holder=np.full((A_LEN, F_LEN, G_LEN), 0, dtype=list)

    # builds link between keys and np-array 
    for j in range(F_LEN):
        for i in range(A_LEN):
            for k in range(G_LEN):
                F=F_RANGE[j]
                A=A_RANGE[i]
                G=G_RANGE[k]

                holder[(i,j,k)]=[0]
                count_holder[(i,j,k)]=[0]
                BG_holder[(i,j,k)]=[0]

                mapping.update({ (A,F,G) : holder[(i,j,k)]})
                count_mapping.update({ (A,F,G) : count_holder[(i,j,k)]})
                BG_mapping.update({ (A,F,G) : BG_holder[(i,j,k)]})

    #form BG mappings
    num_j=0
    for j in AFG_PAIR:
        BG_mapping[tuple(j)][0]=MAX_BG_TEMP[str(num_j)][stat]
        num_j+=1

    # form OS mappings
    for OS in MAX_OS.values():
        if (T<=OS[str(stat)][1]):
            mapping[tuple(OS[str(stat)][0])][0]+=float(OS[str(stat)][1])
            count_mapping[tuple(OS[str(stat)][0])][0]+=1

    #redefines things to be normal arrays inside 
    for j in range(F_LEN):
        for i in range(A_LEN):
            for k in range(G_LEN):
                holder[(i,j,k)]=holder[(i,j,k)][0]
                count_holder[(i,j,k)]=count_holder[(i,j,k)][0]
                BG_holder[(i,j,k)]=BG_holder[(i,j,k)][0]

    #marginalizes array in direction of index
    index=int(index)
    plot_array=holder.sum(index)
    normalizer=count_holder.sum(index)
    BG_plot=BG_holder.sum(index)

    #builds a copy of the original array, but fixes things to be integers 
    pa_cop=np.full(plot_array.shape,0.0)
    n_cop=np.full(normalizer.shape,0.0)
    BG_cop=np.full(BG_plot.shape,0.0)
    z=[(i,j) for i in range(plot_array.shape[0]) for j in range(plot_array.shape[1])]
    for tup in z:
        pa_cop[tup]=plot_array[tup]
        n_cop[tup]=normalizer[tup]
        BG_cop[tup]=BG_plot[tup]

    #normalize
    pa_cop[n_cop>0]=pa_cop[n_cop>0]/n_cop[n_cop>0]

    if (max_OS==-1):
        vmax_val=pa_cop.max()
    else:
        vmax_val=max_OS

    label=np.array(["A", "F", "G"])
    IL=np.array([0, 1, 2])
    IL=IL[ IL != index]

#     XV="self."+label[IL[1]]+"_RANGE"
    XV=label[IL[1]]+"_RANGE"
    XA=eval(XV)

#     YV="self."+label[IL[0]]+"_RANGE"
    YV=label[IL[0]]+"_RANGE"
    YA=eval(YV)

    contours = plt.contour(XA, YA, pa_cop, 10, colors='blue')
    plt.clabel(contours, inline=True, fontsize=8)

    plt.contourf(XA, YA, pa_cop, 100, cmap='hot', alpha=1, vmin=T, vmax=vmax_val);
    plt.colorbar();

    plt.xlabel(label[IL[1]])
    plt.ylabel(label[IL[0]])
    plt.title('Max On Source; t='+str(T))
    plt.show()
    plt.savefig("plots/{}.png".format(plot1))
    
    if (max_BG==-1):
        vmax_val=BG_cop.max()
    else:
        vmax_val=max_BG

    contours = plt.contour(XA, YA, BG_cop, 10, colors='blue')
    plt.clabel(contours, inline=True, fontsize=8)

    plt.contourf(XA, YA, BG_cop, 100, cmap='hot', alpha=1, vmin=T, vmax=vmax_val);
    plt.colorbar();

    plt.xlabel(label[IL[1]])
    plt.ylabel(label[IL[0]])
    plt.title('Max Background; t='+str(T))
    plt.show()
    plt.savefig("plots/{}.png".format(plot2))

# def heatmap(var1, var2, T, max_OS=-1, max_BG=-1, stat=2, plot1="heat_plot1", plot2 = "heat_plot2"):

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--var1', type=int)
    parser.add_argument('--var2', type=int)
    parser.add_argument('--T', type=int)
    parser.add_argument('--stat', nargs='?', const=1, type=int, default=2)
    parser.add_argument('--plot1', nargs='?', const=1, type=str, default="heat_plot1")
    parser.add_argument('--plot2', nargs='?', const=1, type=str, default="heat_plot2")
    args = parser.parse_args() 

    heatmap(args.var1,args.var2,args.T,args.stat,args.plot1,args.plot2)
