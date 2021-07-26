import numpy as np
import math
import json
import random
import argparse

def waveforms(N_A, N_g, N_f, t0_tf, T, B, trials, seedn=1, N=10000,
        inputfile="input", phi0=0, A0=1, Af=50, g0=0, gf=2, F0=90, Ff=110):

    # initalizes the arrays which span parameter space, and their lengths
    A_RANGE=np.linspace(A0,Af,N_A)
    G_RANGE=np.linspace(g0,gf,N_g)
    F_RANGE=np.linspace(F0,Ff,N_f)

    # number of parameters available
    A_LEN, G_LEN, F_LEN = N_A, N_g, N_f
    
    waveform_data={}
    for j in range(trials):
        waveform_data.update({j:[[],[]]})
        
        # calculates random indice for each parameter (A, f, g)
        A_RAN=random.randint(0,A_LEN-1)
        G_RAN=random.randint(0,G_LEN-1)
        F_RAN=random.randint(0,F_LEN-1)
        
        # calculates random parameters A, f, g
        A, gamma, f = A_RANGE[A_RAN], G_RANGE[G_RAN], F_RANGE[F_RAN]
        
        dt=T/N # time resolution

        t0=(T-t0_tf)*np.random.random(1)[0]  # randomly generate start time
        START_INDEX=math.floor(t0/dt)        # find index associated with time

        ##NOTE: using 't0' instead of some multiple of dt may cause issues later
        
        SIG_LEN = (math.floor(t0_tf/dt)+1 if (t0_tf != T) else N) # calculate # of indexes signal takes
        INJECTED = np.zeros(N)                 # initalize injected signal, with N size array of zeroes
        for i in range(SIG_LEN):
            INJECTED[START_INDEX + i]=t0+i*dt       # fill in injected signal with its time stamps

        w = 2 * np.pi * f
        
        # replace timestamps with their displacement values
        SR = INJECTED[START_INDEX : START_INDEX+SIG_LEN][:]
        INJECTED[START_INDEX : START_INDEX+SIG_LEN] = A*np.sin(w*(SR-t0) + phi0)*np.exp(-gamma*(SR-t0))
        
        # Purposed for Quadrature Sum
        D_i = [] # list of each differently seeded waveform
        for n in range(seedn):
            np.random.seed(seed = n)
            NOISE = np.random.normal(scale=(B/(math.sqrt(3)*2)), size=N)  # Noise!
            D_i.append(list(NOISE + INJECTED))  # complete data!
        
        # gets parameters and data for each trial, stuffs it into dictionary
        parameters = [A, f, gamma, t0]
        waveform_data[j][0], waveform_data[j][1] = parameters, D_i
        
    # each trial has list of parameters used and list of data values
    with open("run_uniques/{}-waveform_data.json".format(inputfile) , "w") as f:
        json.dump(waveform_data, f, indent=2, sort_keys=True)

#def waveforms(N_A, N_g, N_f, t0_tf, T, B, trials, seedn=1, N=10000,
#        inputfile="input", phi0=0, A0=1, Af=50, g0=0, gf=2, F0=90, Ff=110):

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--N_A', type=int)
    parser.add_argument('--N_g', type=int)
    parser.add_argument('--N_f', type=int)
    parser.add_argument('--t0_tf', type=int)
    parser.add_argument('--T', type=int)
    parser.add_argument('--B', type=int)
    parser.add_argument('--trials', type=int)
    parser.add_argument('--seedn', nargs='?', const=1, type=int, default=1)
    parser.add_argument('--N', nargs='?', const=1, type=int, default=10000)
    parser.add_argument('--inputfile', nargs='?', const=1, type=str, default='input')

    args = parser.parse_args()

    waveforms(args.N_A,args.N_g,args.N_f,args.t0_tf,args.T,args.B,args.trials,args.seedn,args.N,args.inputfile)
