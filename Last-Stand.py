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

def Rho(data, template, dt):

    ii = 0
    time_slides = []
    R = []
        
    while len(data[ii:]) >= len(template):
        time_slides.append(ii*dt)

        R.append(np.sum((data[ii: len(template) + ii] * template) / (1 + ((data[ii:len(template) + ii] - template) ** 2))))
        ii += 1
        
    return(time_slides, R)

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
        self.y = 10*np.sin(w*t)*np.exp(-gamma*t)
    
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
    outputfile = "{}.txt".format(outputfile)
    np.savetxt(outputfile, output, fmt="%f\t%f\t%f\t%f")
    
    max_index = np.argmax(Ms)
    Maxm = Ms[max_index]
    Maxt = Ts[max_index]
    Maxg = gs[max_index]
    Maxf = fs[max_index]

    return(Maxf, Maxg, Maxt, Maxm)
#search(90,105,0,1,"newdatafile.json"

def heatmap(txtfile,plotfile):
    """
    METHOD: Takes as input a txt with all the frequency, gamma, time,
    and match results. This txt
    PARAMETERS:
    -----------
    txtfile: (String) Name of txt file holding f, g, t, m values
    plotfile: (String) Name given to created heatmap
    """
    # Read in the files
    f, g, t, m = np.loadtxt("{}".format(txtfile), unpack=True)
    f_vals = np.unique(f)
    g_vals = np.unique(g)
    # Construct a 2D grid of the two coordinates
    F_GRID, G_GRID = np.meshgrid(f_vals, g_vals)
    # Reshape the match array in the same shape as the coordinate-grid
    match_grid = m.reshape(len(f_vals), len(g_vals))
    # Plot the heatmap
    pl.rcParams.update({'font.size': 30})
    pl.figure(figsize=(12,10))
    pl.pcolormesh(F_GRID, G_GRID, match_grid.T, shading='auto')
    pl.xlabel('Frequency')
    pl.ylabel('$\\gamma$')
    pl.title("Cross-Correlation: b = 0")
    pl.colorbar()
    pl.savefig("{}.png".format(plotfile))