import json
import numpy as np
import pylab as pl

def waveform(f, A, b, t0, tend, d_end_t=None, gamma=0.0, phi0=0.0, 
             N=1000, verbose=False, seed_number=None, project_name=None):
    """
    METHOD
    ======
    Takes input parameters of a wave and the strength and duration of
    noise, and returns the data.

    PARAMETERS
    ==========
    f : (Float) Frequency of the signal
    A : (Float) Amplitude of the signal
    b : (Float) Amplitude of the noise
    t0 : (Float) Timestamp of the beginning of the signal
    tend : (Float) Time stamp of the end of the signal
    d_end_t : (Float) Time stamp of the end time of the data. Default = None
    gamma : (Float) Attenuation factor of the signal. Default = 0.0
    phi0 : (Float) Initial phase of the signal. Default = 0.0
    N : (Int) Total number of time stamps. Default = 1000
    verbose: (Bool) Set True to get diagnostic stdout. Default = False
    seed_number: (Int) Number set to seed noise. Default = None
    project_name: (String) Name given to png and json file created. Default = None

    OUTPUT
    ======
    A tuple of a float and two numpy arrays (dt, T_full, d), where dt is 
    the resolution of the time series. T_full is the full list of time stamps
    of the data starting at 0 and ending and d_end_t, and d is the 
    corresponding displacement values in the data.

    """
    
    # Conditional for noise duration
    # If the data-end time is supplied to be too small:
    if verbose:
        print("Making sure that the stretch of data is longer than signal")
    assert t0 > 0, "Signal should start later than t=0"
    if (d_end_t is None) or (tend > d_end_t - 10):
        d_end_t = tend + 10
        if verbose:
            print("data end time is set at {}".format(d_end_t))
    
    T = np.linspace(t0, tend, N) # Time stamps of signal
    dt = np.mean(np.diff(T)) # figuring out the resolution of the series
    if verbose:
        print("Mean value of timing resolution = {}".format(dt))
    
    t = t0 # Initializing the time series at the start time
    t_minus = [] # To populate time stamps prior to the signal start
    while t >= 0: # Making sure that we reach all the way back to zero.
        t = t - dt
        t_minus.append(t)  # Create time spamps from (t0-dt) to 0

    t_minus = np.array(t_minus)[::-1]  # Reverse to be from 0 to t0
    t_minus = t_minus[t_minus >= 0]  # Eliminate numbers less than 0
    
    t_plus = np.arange(tend+dt, d_end_t, dt)  # Time stamps from (tend+dt) to d_end_t, in dt's
    
    T_full = np.hstack((t_minus, T, t_plus))  # Connect time stamps
    
    dev = np.std(np.diff(T_full))  # Standard deviation in dt's of T_full
    if verbose:
        print("Standard deviation of the resolution of time = {}".format(dev))

    if verbose:
        print("Creating time series of the signal...")
    w = 2 * np.pi * f  
    y = A*np.sin(w*T + phi0)*np.exp(-gamma*(T-t0))

    
    # Padding of signal data
    if verbose:
        print("Creating the zero-padded signal...")
    y_minus = np.zeros_like(t_minus)
    y_plus = np.zeros_like(t_plus)
    y_full = np.hstack((y_minus, y, y_plus))
    
    if verbose:
        print("Creating random noise...")
    if seed_number is None:
        seed_number = 1
    np.random.seed(seed = seed_number)
    noise = -b+2*b*np.random.random(len(T_full))  # Noise!
    
    if verbose:
        print("Creating final data")
    d = noise + y_full  # Complete Data!
    
    # Graphing   
    pl.rcParams.update({'font.size': 18})
    pl.figure(figsize=(20,15))
    pl.plot(T_full, noise, color = 'green', linewidth=2)  # Noise
    pl.plot(T_full, d, color = 'black', linewidth=2)  # Combined
    pl.plot(T, y, color = 'orange', linewidth=2)  # Signal
    pl.xlabel("Time")
    pl.ylabel("displacement")
    text = "f={}; A={}; b={}; t0={}; tend={}; gamma={}; N={}"
    pl.title(text.format(f, A, b, t0, tend, gamma, N))
    #if project_name is None:
    #    project_name = "test"
    #pl.savefig("figures/{}-waveform_plot-f_{}-A_{}-b_{}-t0_{}-tend_{}-gamma_{}-seed_{}.png".format(project_name, f, A, b, t0, tend, gamma, seed_number))
    
    T_full = T_full
    d = d
    #data = {"dt" : dt, "t_full" : T_full, "d" : d}
    #outputfile = "data/{}-waveform_data-f_{}-A_{}-b_{}-t0_{}-tend_{}-gamma_{}-seed_{}.json".format(project_name, f, A, b, t0, tend, gamma, seed_number)
    #with open(outputfile, "w") as f:
    #    json.dump(data, f, indent=2, sort_keys=True)
    return(dt, T_full, d)
