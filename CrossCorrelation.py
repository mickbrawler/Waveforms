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
        match = []
        
        while len(self.d[ii:]) >= len(self.y):
            time_slides.append(ii*self.dt)
            match.append(np.sum(self.d[ii: len(self.y) + ii] * self.y))
            ii += 1
        
        self.match = np.array(match)
        self.time_slides = np.array(time_slides)
        
        return(self.time_slides, self.match) 

def search(f_low, f_hi, gamma_low, gamma_hi,datafile):
    """
    METHOD: Takes as input the upper and lower values of frequency and gammas, constructs 
    a bank of templates using this range of values, and then computes the maximum of the 
    match for each templates in the bank and its corresponding time and returns that
    """
    f = np.arange(f_low, f_hi+1, 1)
    g = np.arange(gamma_low, gamma_hi + 0.1, 0.1)

    fs = []
    gs = []
    Ms = []
    Ts = []

    Obj = Crosscor(datafile)
    for i in f:
        for j in g:
            Obj.template(i,j,40)
            t,m = Obj.match()
            M = m[np.argmax(m)] # Max match
            T = t[np.argmax(m)] # Time associated with max match
            fs.append(i)
            gs.append(j)
            Ts.append(T)
            Ms.append(M)

            output = np.vstack((fs,gs,Ts,Ms)).T
            np.savetxt("firstattempt.txt", output, fmt="%f\t%f\t%f\t%f")

            f,g,t,m = np.loadtxt("firstattempt.txt",unpack=True)

    Maxm = m[np.argmax(m)]
    Maxt = t[np.argmax(m)]
    Maxg = g[np.argmax(m)]
    Maxf = f[np.argmax(m)]
    
    return(Maxm,Maxt,Maxg,Maxf)

search(90,105,0,1,"newdatafile.json")
