import numpy as np
import pylab as pl

def waveform(f,A,b,t0,tend,t,gamma,phi0,n = 1000):
    
    # Conditional for noise duration
    if tend > t: 
        t = tend + 2

    # Method of finding dt (Generally it uses the next best n if it doesn't work)
    n = np.arange(1,n+1) # 1 to given stepcount will be used to find stepcount that works best

    dt = (tend - t0)/n # dt is now the dt's from t0 to tend from the stepcounts from 1 to given

    n, dt = np.meshgrid(n, dt) # Creates matrixes of n (x axis), dt (y axis)

    tcandidates = dt * n # Find all potential times for each dt

    correctdt = min(dt[(tcandidates == t0) | (tcandidates == (tend-t0))]) # Match potential times with t0 and (t-tend)
         
    # Obtain n value for all three intervals: 0_t0, t0_tend, tend_t
    leftn = (t0-0)/correctdt 
    signaln = (tend - t0)/correctdt 
    rightn = (t-tend)/correctdt
    
    totlen = leftn + signaln + rightn # Number upto which dt will be multiplied by to get time stamps

    T = correctdt * np.arange(totlen+1) # Timestamps!
 
    #Endgame
    w = 2*np.pi*f
    
    BoolT = (T>t0) & (T<tend)    # Create boolean version of T array where Falses are the padding
    BoolT = BoolT.astype(int)
    Y = A*np.sin(w*T + phi0)*np.exp(-gamma*T)
    Y = Y * BoolT    # Isolate signal data, make the rest 0s
    
    np.random.seed(seed = 1)
    y = b * np.random.uniform(-1,1,len(Y)) # Noise!
    
    d = Y + y # Complete Data!
    
    # Graphing
    pl.rcParams.update({'font.size': 18})
    pl.figure(figsize=(12,10))
    pl.plot(T, y, color = 'green', linewidth=2) # Noise
    pl.plot(T, d, color = 'black', linewidth=2) # Combined
    pl.plot(T, Y, color = 'orange', linewidth=2) # Signal
    pl.savefig('waveform.png')
        
    return (T)

# f,A,b,t0,tend,t,gamma,phi0
waveform(1,1,1,1,3,5,0,0)
