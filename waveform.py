import numpy as np
import pylab as pl

def get_time_series(f, A, t, gamma, t0, tend):
    
    if t0 > t:                                           #Check if t0/tend are greater than t
        return "Invalid signal start indice"
    elif tend > t:
        return "Invalid signal end indice"
    
    phi0 = np.pi/6  
    w = 2*np.pi*f
    
    # Signal time series
    T = np.linspace(t0,tend, int(((tend-t0)/t)*10000))   # Array from t0 to tend,spaced out a fraction of the # of times t is
    # T is a numpy array that starts at t0 and ends at tend, with linearly spaced times 

    Y = A*np.sin(w*T + phi0)*np.exp(-gamma*T)
    # Y is the a sinusoidal displacement that lasts for the entire duration of time (t0, tend)

    # Carrier time series
    t = np.linspace(0, t, 10000)
    # t is a numpy array that starts at 0 and ends at t.

    y = A*np.sin(w*t + phi0)*np.exp(-gamma*t)
    # y is the displacement during this time duration t
    
    # Plotting
    pl.rcParams.update({'font.size': 18})
    pl.figure(figsize=(12,10))
    
    pl.plot(t, y, linewidth=2)
    pl.plot(T, Y, linewidth=2)                           # Plot of "signal" overlaps carrier wave
    
    pl.plot(t, A*np.exp(-gamma*t), 'k--', linewidth=2)
    pl.plot(t, -A*np.exp(-gamma*t), 'k--', linewidth=2)
    pl.xlabel('Time (s)')
    pl.ylabel('$\\sin(\\omega t)$')
    
    return [T,Y]
