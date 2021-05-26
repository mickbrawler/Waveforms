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
    
import importlib 

f = np.arange(90,105+1,1)
g = np.arange(0,1+.1,.1)

fs = []
gs = []
Ms = []
Ts = []

for i in f:
    for j in g:
        
        Obj = Crosscor("newdatafile.json")
        Obj.template(i,j,40)
        
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
