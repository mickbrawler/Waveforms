import numpy as np
import pylab as pl

def waveform(f,A,b,t0,tend,d_end_t,gamma,phi0,N=10000):
    
    T = np.linspace(t0,tend,N)
    dt = np.mean(np.diff(T))
    
    t = t0
    t_minus = []
    while t >= 0:
        t = t - dt
        t_minus.append(t)
    t_minus = np.array(t_minus)[::-1]
    t_minus = t_minus[t_minus >= 0]
    t_plus = np.arange(tend+dt,d_end_t,dt)
    T_full = np.hstack((t_minus,T,t_plus))
    print(T_full)
    print(np.std(np.diff(T_full)))
    
    w = 2 * np.pi * f
    y = A*np.sin(w*T + phi0)*np.exp(-gamma*(T-t0))
    y_minus = np.zeros_like(t_minus)
    y_plus = np.zeros_like(t_plus)
    y_full = np.hstack((y_minus, y, y_plus))
    
    noise = -b+2*b*np.random.random(len(T_full))
    
    return(T_full,y_full+noise,T,y)

T,y,t,signal = waveform(20,10,100,2,3,4,1,0)
pl.plot(T,y)
pl.plot(t,signal,'r')