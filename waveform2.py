import numpy as np
import pylab as pl

def waveform(f,A,b,t0,tend,d_end_t,gamma,phi0,N=1000):
    
    # Conditional for noise duration
    if tend > d_end_t: 
        d_end_t = tend + 2
    
    T = np.linspace(t0,tend,N) # Time stamps of signal
    dt = np.mean(np.diff(T))
    
    t = t0
    t_minus = []
    while t >= 0:
        t = t - dt
        t_minus.append(t)  # Create time spamps from (t0-dt) to 0

    t_minus = np.array(t_minus)[::-1]  # Reverse to be from 0 to t0
    t_minus = t_minus[t_minus >= 0]  # Eliminate numbers less than 0
    
    t_plus = np.arange(tend+dt,d_end_t,dt)  # Time stamps from (tend+dt) to d_end_t, in dt's
    
    T_full = np.hstack((t_minus,T,t_plus))  # Connect time stamps
    
    dev = np.std(np.diff(T_full))  # Standard deviation in dt's of T_full
    
    w = 2 * np.pi * f  
    y = A*np.sin(w*T + phi0)*np.exp(-gamma*(T-t0))
    
    # Padding of signal data
    y_minus = np.zeros_like(t_minus)
    y_plus = np.zeros_like(t_plus)
    y_full = np.hstack((y_minus, y, y_plus))
    
    np.random.seed(seed = 1)
    noise = -b+2*b*np.random.random(len(T_full))  # Noise!
    
    d = noise + y_full  # Complete Data!
    
    # Graphing   
    pl.rcParams.update({'font.size': 18})
    pl.figure(figsize=(12,10))
    pl.plot(T_full, noise, color = 'green', linewidth=2)  # Noise
    pl.plot(T_full, d, color = 'black', linewidth=2)  # Combined
    pl.plot(T_full, y_full, color = 'orange', linewidth=2)  # Signal
    pl.savefig('waveform.png')
    
    return(dt,T_full,d)

# f,A,b,t0,tend,d_end_t,gamma,phi0

def template(f,gamma,t0,tend):
    
    dt = waveform(1,1,1,1,3,5,0,0)[0]
    t = np.arange(t0,tend,dt)
    w = 2 * np.pi * f  
    y = np.sin(w*t)*np.exp(-gamma*t0)
    
    pl.rcParams.update({'font.size': 18})
    pl.figure(figsize=(12,10))
    pl.plot(t, y, color = 'green', linewidth=2)
    
    return(t,y)

# f,gamma, t0, tend                             
template(1,0,1,3)
