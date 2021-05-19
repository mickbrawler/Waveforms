import numpy as np
import pylab as pl

def waveform(f,A,b,t0,tend,d_end_t,gamma,phi0,N=1000):
    
    global dt      
    global T_full 
    global d
    
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
    
# f,A,b,t0,tend,d_end_t,gamma,phi0
waveform(1,1,1,2,3,5,0,0)

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
