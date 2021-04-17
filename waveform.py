import numpy as np
import pylab as pl

def get_time_series(f, A, t, gamma, t0, tend):

    if t0 > t:
        return "Invalid signal start indice"
    elif tend > t:
        return "Invalid signal end indice"

    phi0 = np.pi/6
    w = 2*np.pi*f

    #Signal time series
    T = np.linspace(t0, tend, int(((tend-t0)/t)*10000))
    Y = A*np.sin(w*t + phi0)*np.exp(-gamma*t)

    #Carrier time series
    t = np.linspace(0,t,10000)
    y = A*np.sin(w*t +phi0)*np.exp(-gamma*t)

    #Plotting
    pl.rcParams.update({'font.size':18})
    pl.figure(figsize=(12,10))

    pl.plot(t, y, linewidth=2)
    pl.plot(T, Y, linewidth=2)

    pl.plot(t, A*np.exp(-gamma*t), 'k--', linewidth=2)
    pl.plot(t, -A*np.exp(-gamma*t), 'k--', linewidth=2)
    pl.xlabel('Time (s)')
    pl.ylabel('$\\sin(\\omega t)$')

    return[T,Y]
