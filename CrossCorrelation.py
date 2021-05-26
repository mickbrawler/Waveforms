import json
import numpy as np
import pylab as pl

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
        
    def match(self):
        
        ii = 0
        time_slides = []
        M = []
        
        while len(self.d[ii:]) >= len(self.y):
            time_slides.append(ii*self.dt)
            M.append(np.sum(self.d[ii: len(self.y) + ii] * self.y))
            ii += 1
        
        self.M = np.array(M)
        self.time_slides = np.array(time_slides)
        
        return(self.time_slides, self.M)

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
            print("f = {}\t gamma = {}".format(i, j))
            Obj.template(i, j, tmplt_dur)
            t, m = Obj.match()
            M = m[np.argmax(m)] # Max match
            T = t[np.argmax(m)] # Time associated with max match
            fs.append(i)
            gs.append(j)
            Ts.append(T)
            Ms.append(M)

    output = np.vstack((fs,gs,Ts,Ms)).T
    np.savetxt(outputfile, output, fmt="%f\t%f\t%f\t%f")

    max_index = np.argmax(Ms)
    Maxm = Ms[max_index]
    Maxt = Ts[max_index]
    Maxg = gs[max_index]
    Maxf = fs[max_index]

    return(Maxm, Maxt, Maxg, Maxf)

#search(90,105,0,1,"newdatafile.json"
