import json
import numpy as np
from matplotlib import pyplot as plt
import math
import time
import argparse

# Produces a template given a position in parameter space
def template(A, f, gamma, duration, dt):
    t = np.arange(0, duration + dt, dt)
    w = 2 * np.pi * f
    return A * np.sin(w*t)*np.exp(-gamma*t)

# Produces a cross corelation function given a data input and a template in parameter space
def CrossCorrelation(data, template, dt):
    ii = 0
    M = []

    while len(data[ii:]) >= len(template):
        M.append(np.sum((data[ii: len(template) + ii] * template)))
        ii+=1

    return M

# Produces chi square at each "slide"
def ChiSquare(data, template, dt):
    ii = 0
    C = []

    while len(data[ii:]) >= len(template):
        C.append(-1 * np.sum((data[ii: len(template) + ii] - template) ** 2))
        ii += 1

    return C

def modulator(rho_ij, D, dt):

    rho_mod_D, RHO_ij = [] , rho_ij[:]

    dn , L = math.floor(2*D/dt) , len(RHO_ij)

    for i in range(0,L-(L%dn),dn):
        rho_mod_D.append(max(RHO_ij[i:i+dn]))

    if (L-(L%dn)) != L :
        rho_mod_D.append(max(RHO_ij[L-(L%dn):L]))

    return rho_mod_D

def statudio(trialn, D, N_A, N_g, N_f, t0_tf, T, trials, run1=True, seedn=1,
             N_t=10000, inputfile="input", A0=1, Af=50, g0=0, gf=2, F0=90,
             Ff=110):

    # initalizes the arrays which span parameter space, and their lengths
    A_RANGE=np.linspace(A0,Af,N_A)
    G_RANGE=np.linspace(g0,gf,N_g)
    F_RANGE=np.linspace(F0,Ff,N_f)

    A_LEN, G_LEN, F_LEN = N_A, N_g, N_f

    # constructs timestep resolution, and saves N and t0/tf internally 
    N, dt = N_t, T/N_t

    # constructs time range to pick injected signal start time from/ corresponding length 
    t_RANGE=np.linspace(0,T-(t0_tf),int(N_t*(1-((t0_tf)/(T)))))
    t_LEN=len(t_RANGE)

    # initialize arrays for various data/cross-correlations/chi-squares 
    noise = []

    # constructs all templates which correspond to points in the parameter space
    TEMPLATES_AFG=[ template( A, f, g, t0_tf, dt) for A in A_RANGE 
                   for g in G_RANGE for f in F_RANGE]

    AFG_PAIR=[ [A, f, g] for A in A_RANGE 
                   for g in G_RANGE for f in F_RANGE]

    # Reads waveform data file 
    with open("run_uniques/{}-waveform_data.json".format(inputfile),"r") as f: 
        waveform_data = json.load(f)

    waveform_data = waveform_data[str(trialn)] # trialn's parameters and differently seeded data
    
    output={}

    output.update({trialn:[[],[],[],[]]})

    # isolates random a-g-f pair / data set 
    temp_AGFT, data = waveform_data[0], waveform_data[1]

    noise.append(data) 

    output[trialn][0], output[trialn][1] = temp_AGFT, data  # stores random a-g-f pair / data set 

    Quad_CRS = []
    Quad_CHI = []
    
    # performs base static calculation across parameter space
    # Quadrature Sum
    for n in range(seedn): # Use seedn as index for data

        CRS_COR, CHI_SQR = [[],[]]

        for temp in TEMPLATES_AFG:

            CC_dh = list(CrossCorrelation(data[n], temp, dt))
            CRS_COR.append(CC_dh)

            CS_dh = list(ChiSquare(data[n], temp, dt))
            CHI_SQR.append(CS_dh)

        Quad_CRS.append(CRS_COR) # now a 3d list of seedn statistics, with 2d list statistics per waveform
        Quad_CHI.append(CHI_SQR)

    CRS_COR = np.sum(np.array(Quad_CRS) ** 2, axis = 0) ** .5 # Quadrature sum of each seed's statistic
    CHI_SQR = np.sum(np.array(Quad_CHI) ** 2, axis = 0) ** .5

    # stores base statistics to attribute
    cross_cor = CRS_COR
    chi = CHI_SQR
    output[trialn][2], output[trialn][3] = CRS_COR.tolist(), CHI_SQR.tolist() # store quadrature summed based statistics
    
#calculates test statistic, stores it internally,
#and returns a copy of it as a dictionary 
# trying connect stat and anlysis function
    
    # counts number of tempates in parameter space
    PSPACE_LEN = len(AFG_PAIR)
    
    # String to equation!
    stats = ["CC_IJ","CS_IJ","CC_IJ/abs((1+CS_IJ))"]

    # initalizes rho statistic dictionary/stat outputs collector
    RHO = {}
    stat_outputs = []
    
    for stat in stats: # objective is that the appropriate jsons holds multiple stat results

        # indexed to loops through parameter space templates and
        # calculates each rho_ij given template j
        rho_i = []
        for j in range(PSPACE_LEN):
            CC_IJ = np.array(cross_cor[j][:])
            CS_IJ = np.array(chi[j][:])

            # Evaluates string (Exec gave issues... eval is the same concept though)
            p = eval(stat)
            rho_i.append(list(p))
        stat_outputs.append(rho_i) # stat_outputs is 3d list holding a 2d list of a stat's outputs for each template
    RHO.update({ trialn : stat_outputs })

    if (2*D >= dt):
        
        # gets the length of linearized template space
        TEMP_LEN=len(cross_cor) # number of templates 

        RHO_MOD={}
        MAX_OS={}
        MAX_BG_TEMP={}
        BG_VALS_IJ, FG_VAL_IJ = [], []
        pot_thresholds = {}
        
        # seperates fg value from bg value
        T0_2D = math.floor(output[trialn][0][3]/(2*D))

        for j in range(TEMP_LEN):
            MAX_BG_TEMP.update({ j : list(np.zeros(len(stats)))})

        for stat in range(len(stats)):
            
            pot_thresholds.update({stat:[]})
            
            BG_VALS_ij, FG_VAL_ij = [], []
            for j in range(TEMP_LEN):

                # calculates bg values + fg values
                BG_VALS_ij.append(modulator(RHO[trialn][stat][j][:],D,dt))

                FG_VAL_ij.append(BG_VALS_ij[j].pop(T0_2D))
                    
            BG_VALS_IJ.append(BG_VALS_ij) # 3d list
            FG_VAL_IJ.append(FG_VAL_ij) # 2d list
            
            pot_thresholds[stat] += FG_VAL_ij # for each stat key, forground values of trialn are added
        RHO_MOD.update({ trialn: [ BG_VALS_IJ, FG_VAL_IJ ] }) # these are the peaks we look for
        
        stat_dict = {}
        for stat in range(len(stats)):
            stat_dict.update({stat:[tuple(output[trialn][0][0:3]),max(RHO_MOD[trialn][1][stat])] })
        MAX_OS.update({trialn: stat_dict}) # stores t0 and forground peaks for each trial

        for stat in range(len(stats)):
            for j in range(TEMP_LEN):
                new_max=max(RHO_MOD[trialn][0][stat][j]) # max ofsource peak for each template

                MAX_BG_TEMP[j][stat] = new_max # every value in MAX_BG_TEMP dictionary changes from 0 to that templates max
    else: 
        print("invalid D; it is required that 2*D >= T/N")

    # All jsons below will serve post processing

    # output holds trialn's parameters, data, cross-correlation, chi-square
    with open("output_folder/output_{}.json".format(trialn), "w") as f: 
        json.dump(output, f, indent=2, sort_keys=True)

    # RHO_MOD holds trialn's background and forground values
    with open("Peaks_folder/Peaks_{}.json".format(trialn), "w") as f:
        json.dump(RHO_MOD, f, indent=2, sort_keys=True)

    # MAX_OS holds trialn's max onsource peak
    with open("Max_OS_folder/Max_OS_{}.json".format(trialn), "w") as f:
        json.dump(MAX_OS, f, indent=2, sort_keys=True)
    
    # MAX_BG_TEMP and pot_thresholds were devised wrong. They used to rewrite the same file with
    # changed or more value in every run. This won't work when parallelized as some runs may finish
    # faster than even the first. They are now made for each trial, but because of their different form
    # compared to the simple first three, weirder json combiners must be used on them

    with open("Max_BG_TEMP_folder/Max_BG_TEMP{}.json".format(trialn), "w") as f: # Merged
        json.dump(MAX_BG_TEMP, f, indent=2, sort_keys=True)
        
    with open("thresholds_folder/thresholds{}.json".format(trialn), "w") as f:
        json.dump(pot_thresholds, f, indent=2, sort_keys=True)
    
    # Now only the essentials file is made when run1 is True
    # essentials holds all the values that don't change for each trialn run
    if run1 == True:
        essent = {"essentials":[A_LEN,F_LEN,G_LEN,list(F_RANGE),list(A_RANGE),list(G_RANGE),AFG_PAIR,trials,stats]}
        with open("run_uniques/essentials.json", "w") as f:
            json.dump(essent, f, indent=2, sort_keys=True)

#def statudio(trialn, D, N_A, N_g, N_f, t0_tf, T, trials, run1=True, seedn=1,
#             N_t=10000, inputfile="input", A0=1, Af=50, g0=0, gf=2, F0=90,
#             Ff=110):

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--trialn', type=int)
    parser.add_argument('--D', type=float)
    parser.add_argument('--N_A', type=int)
    parser.add_argument('--N_g', type=int)
    parser.add_argument('--N_f', type=int)
    parser.add_argument('--t0_tf', type=int)
    parser.add_argument('--T', type=int)
    parser.add_argument('--trials', type=int)
    parser.add_argument('--run1', nargs='?', const=1, type=bool, default=True)
    parser.add_argument('--seedn', nargs='?', const=1, type=int, default=1)
    parser.add_argument('--N_t', nargs='?', const=1, type=int, default=10000)
    parser.add_argument('--inputfile', nargs='?', const=1, type=str, default='input')

    args = parser.parse_args()

    statudio(args.trialn,args.D,args.N_A,args.N_g,args.N_f,args.t0_tf,args.T,args.trials,args.run1,args.seedn,args.N_t,args.inputfile)

