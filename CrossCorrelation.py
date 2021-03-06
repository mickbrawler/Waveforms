import json
import numpy as np
import pylab as pl
#from numba import jit

#@jit(nopython=True)
def match(data, template, dt):
        
    ii = 0
    time_slides = []
    M = []
        
    while len(data[ii:]) >= len(template):
        time_slides.append(ii*dt)
        
        M.append(np.sum((data[ii: len(template) + ii] * template)))
        #M.append(np.sum((data[ii: len(template) + ii] * template) / (1 + ((data[ii:len(template) + ii] - template) ** 2))))
        #M.append(np.sum(((data[ii: len(template) + ii] * template) / (1 + (((data[ii:len(template) + ii] - template) ** 2) / template)))))
        ii += 1
        
    M = np.array(M)
    time_slides = np.array(time_slides)
        
    return(time_slides, M)

def ChiSquare(data, template, dt):

    ii = 0
    time_slides = []
    C = []
        
    while len(data[ii:]) >= len(template):
        time_slides.append(ii*dt)

        C.append(np.sum((data[ii:len(template) + ii] - template) ** 2))
        ii += 1
        
    return(time_slides, C)

class Crosscor:
    def __init__(self, filename):
        with open(filename, "r") as f:
            data = json.load(f)
        self.dt = data["dt"]
        self.tfull = np.array(data["t_full"])
        self.d = np.array(data["d"])

    def template(self, f, gamma, duration):
        
        t = np.arange(0, duration + self.dt, self.dt)
        self.t = t
        w = 2 * np.pi * f
        self.y = np.sin(w*t)*np.exp(-gamma*t)
    
#@jit(nopython=True)
def matcharrays(f_low, f_hi, gamma_low, gamma_hi, datafile,
           tmplt_dur, matchfile, df=1.0, dg=0.1):
    
    """
    METHOD: Takes as input the upper and lower values of frequency 
    and gammas, constructs a bank of templates using this range of values, 
    and then computes the matches for each set.

    PARAMETERS:
    -----------
    f_low: Lower bound of the frequency grid
    f_hi: Upper bound of the frequency grid
    gamma_low: Lower bound of the gamma grid
    gamma_hi: Upper bound of the gamma grid
    datafile: The JSON file with the data time series
    tmplt_dur: The duration of the templates
    df: Step-size in frequency (default = 1.0)
    dg: Step-size in gamma (default = 0.1)
    matchfile: Name for json file with matches
    
    OUTPUT:  Creates json file with frequency in a list, gamma in a list,
    time in a 2d list, and matches in a 2d list.
    """
    
    f = np.arange(f_low, f_hi+df, df)
    g = np.arange(gamma_low, gamma_hi + dg, dg)

    F = []
    G = []
    M = []
    T = []

    Obj = Crosscor(datafile)
    for i in f:
        for j in g:
            Obj.template(i, j, tmplt_dur)
            t, m = match(Obj.d, Obj.y, Obj.dt)
            T.append(list(t))
            M.append(list(m))
            F.append(i)
            G.append(j)

    data = {"f" : F, "g" : G, "t" : T, "m" : M}
    
    outputfile = "results/{}.json".format(matchfile)
    with open(outputfile, "w") as f:
        json.dump(data, f, indent=2, sort_keys=True)

def combine(file1,file2,file3,file4,file5,outputfile):
    
    """
    METHOD: Takes two json files (can be adjusted to accept more) and performs
    the combined m operation to lower effect of noise. 

    PARAMETERS:
    -----------
    file1: (json) First waveform json file
    file2: (json) Second waveform json file
    file3: (json) Third waveform json file
    file4: (json) Four waveform json file
    file5: (json) Five waveform json file
    outputfile: (String) Name for txt file that will store data for analysis

    OUTPUT: Creates txt file with lists of each run of frequency, gamma, time
    slides, and matches. Return global maximum values
    """

    with open(file1, "r") as f:
        data = json.load(f)
    f1 = data["f"]
    g = data["g"]
    t = data["t"]
    m1 = data["m"]

    with open(file2, "r") as f:
        data = json.load(f)
    m2 = data["m"]
    
    with open(file3, "r") as f:
        data = json.load(f)
    m3 = data["m"]
    
    with open(file4, "r") as f:
        data = json.load(f)
    m4 = data["m"]
    
    with open(file5, "r") as f:
        data = json.load(f)
    m5 = data["m"]

    m1 = np.array(m1)
    m2 = np.array(m2)
    m3 = np.array(m3)
    m4 = np.array(m4)
    m5 = np.array(m5)
    t = np.array(t)

    combinedM = ((m1 ** 2) + (m2 ** 2) + (m3 ** 2) + (m4 ** 2) + (m5 ** 2)) ** .5
   
    M = []
    T = []
    step = 0

    for i in combinedM:
        x = np.argmax(i)
        T.append(t[step,x])
        M.append(combinedM[step,x])
        step += 1
    
    max_Match = np.argmax(np.array(M))

    globF = f1[max_Match]
    globG = g[max_Match]
    globT = T[max_Match]
    globM = M[max_Match]

    output = np.vstack((f1,g,T,M)).T
    outputfile = "results/{}.txt".format(outputfile)
    np.savetxt(outputfile, output, fmt="%f\t%f\t%f\t%f")

    return(globF,globG,globT,globM)

def search(f_low, f_hi, gamma_low, gamma_hi, datafile,
           tmplt_dur, outputfile, df=1.0, dg=0.1):
    """
    METHOD: Takes as input the upper and lower values of frequency 
    and gammas, constructs a bank of templates using this range of values, 
    and then computes the maximum of the match for each templates in the 
    bank and its corresponding time and returns that

    PARAMETERS:
    -----------
    f_low: Lower bound of the frequency grid
    f_hi: Upper bound of the frequency grid
    gamma_low: Lower bound of the gamma grid
    gamma_hi: Upper bound of the gamma grid
    datafile: The JSON file with the data time series
    tmplt_dur: The duration of the templates
    df: Step-size in frequency (default = 1.0)
    dg: Step-size in gamma (default = 0.1)
    outputfile: The txt file with the two dimensional search results
    """
    f = np.arange(f_low, f_hi+df, df)
    g = np.arange(gamma_low, gamma_hi + dg, dg)

    fs = []
    gs = []
    Ms = []
    Ts = []

    Obj = Crosscor(datafile)
    for i in f:
        for j in g:
            Obj.template(i, j, tmplt_dur)
            t, m = match(Obj.d, Obj.y, Obj.dt)
            M = m[np.argmax(m)] # Max match
            T = t[np.argmax(m)] # Time associated with max match
            fs.append(i)
            gs.append(j)
            Ts.append(T)
            Ms.append(M)

    output = np.vstack((fs,gs,Ts,Ms)).T
    outputfile = "results/{}".format(outputfile)
    np.savetxt(outputfile, output, fmt="%f\t%f\t%f\t%f")

    max_index = np.argmax(Ms)
    Maxm = Ms[max_index]
    Maxt = Ts[max_index]
    Maxg = gs[max_index]
    Maxf = fs[max_index]

    return(Maxm, Maxt, Maxg, Maxf)

#search(90,105,0,1,"newdatafile.json"

def plot(txtfile,plotfile):
    """
    METHOD: Takes as input the search function's outputfile, loads it, then
    isolates the frequency and match values to use as the x and y of the plot
    
    PARAMETERS:
    ----------
    txtfile: (txt) File that holds frequency, gamma, time, and match values
    plotfile: (png) Plot of the frequency and match values of provided file
    """
    results = np.loadtxt(txtfile)
    f_array = results[:,0]
    m_array = results[:,3]
    pl.rcParams.update({'font.size':18})
    pl.figure(figsize=(20,15))
    pl.plot(f_array,m_array, linewidth=2)
    pl.xlabel("Frequency")
    pl.ylabel("Match")
    pl.savefig("figures/{}".format(plotfile))
