import numpy as np
import pylab as pl

def waveform(f,A,t,t0,tend,gamma,phi0):
    
    if tend > t: # conditional for noise duration
        t = tend + 2

    # Method of finding dt
    i = 100 # Extent of n's tried to get dt

    n = np.arange(1,i+1)

    dt = (tend - t0)/ n

    n, dt = np.meshgrid(n, dt) # creates matrix of n (x axis), dt (y axis)

    multipledt = dt * n # multiply the two matrices returns dt multiplyed by 1 to i (i being how far you want to go to see if dt reaches t0 from 0)

    correctdt = dt[multipledt == t0] # make values in multipledt that are t0 True, return value at position of the Trues

    useddt = min(correctdt)
    
    # Obtain number of values from 0 to t0 to tend to 
    paddingL = (t0-0)/useddt # How many elements between 0 and t0
    
    signal = (tend - t0)/useddt # How many elements between t0 and tend
    
    paddingR = (t-tend)/useddt

    totlen = paddingL + signal + paddingR # number upto which dt will be multiplied by to get time stamps

    T = useddt * np.arange(totlen+1) # timestamps!
    
    # Must make Y have 0s for values outside t0 to tend.
    # This way d and Y graph are correct
    # Create boolean version of T array where Falses are the padding
    # Place on top of Y, makes Falses 0s

    
    BoolT = (T>t0) & (T<tend)
    BoolT = BoolT.astype(int)
    
    #Endgame
    w = 2*np.pi*f
    
    Y = A*np.sin(w*T + phi0)*np.exp(-gamma*T)
    Y = Y * BoolT
    
    np.random.seed(seed = 1)
    y = np.random.random(len(Y))
    
    d = Y + y # Complete Data
    
    # Graphing (Why not)
    pl.rcParams.update({'font.size': 18})
    pl.figure(figsize=(12,10))
    
    pl.plot(T, y, color = 'green', linewidth=2) # Noise
    pl.plot(T, d, color = 'black', linewidth=2) # Combined
    pl.plot(T, Y, color = 'orange', linewidth=2) # Signal
    pl.savefig('waveform.png')
    
waveform(1,1,10,2,5,0,0)
