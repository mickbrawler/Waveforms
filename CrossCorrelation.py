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
