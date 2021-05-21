import numpy as np
import pylab as pl

class Crosscor:

def template(f,gamma,t0,tend):
    
    global t
    global y
    global tzero
    global tEnd
    
    tzero = t0
    tEnd = tend
    
    t = np.arange(t0,tend,dt) # Template time series is made in same dt as waveform
    
    if t[-1] != tend:         # The dt may not allow t to reach tend (Rounding issue?), add it then
        t = np.append(t,tend)

    w = 2 * np.pi * f
    y = np.sin(w*t)*np.exp(-gamma*t)
    
# f,gamma, t0, tend
template(1,0,2,3)

def match():
    
    t_in_T = (T_full>=tzero) & (T_full<=tEnd) #t values in T_full are made TRUE (Remember T_full is same size as d)
    
    usedY = d[t_in_T] # Use boolean masking on d, now d values associated with values of template times are used

    M = np.sum(usedY * y) #Multiply result with template's y, Add em
    
    return(M)

match()
# Sliding functionality could involve adding usedY by dt in each loop when starting from 0
