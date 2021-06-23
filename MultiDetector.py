import waveform
import numpy as np
import pylab as pl

class Crosscor:
    def __init__(self, dt, T_full, d):
        self.dt = dt
        self.tfull = T_full
        self.d = d 

    def template(self, A, f, gamma, duration):
        
        t = np.arange(0, duration + self.dt, self.dt)
        self.t = t
        w = 2 * np.pi * f
        self.y = A * np.sin(w*t)*np.exp(-gamma*t)

# Produces rho at each "slide"
def match(data, template, dt):

    """
    METHOD
    ======
    Uses the data array, template array, and dt float. Performs 
    cross correlation between segments of the data and the whole template.
    Chi Square analysis is utilized.
    PARAMETERS
    ==========
    data: (Array) Time series of the waveform with an embedded signal
    template: (Array)  Time series of the template
    dt: (Float) Resolution of time array of waveform
    OUTPUT
    ======
    Returns array of time slides and rho outputs from the 
    template.
    """
        
    ii = 0
    time_slides = []
    M = []
        
    while len(data[ii:]) >= len(template):
        time_slides.append(ii*dt)
        
        #M.append(np.sum((data[ii: len(template) + ii] * template)))
        M.append(np.sum((data[ii: len(template) + ii] * template) / (1 + ((data[ii:len(template) + ii] - template) ** 2))))
        ii += 1
        
    return(time_slides, M)

def multidetector(f, A, b, t0, tend, A_low, A_hi, f_low, f_hi, gamma_low, 
                  gamma_hi, tmplt_dur, ndet, gamma=0.0, phi0=0.0, N=1000, 
                  d_end_t=None, verbose=False, df=1.0, dg=0.1, da=1.0):
    

    Ts, Ms, As, Fs, Gs = [], [], [], [], [] 

    for n in range(ndet):

        # (needs arguments of waveform) (still produces json/png)
        dt, T_full, d = waveform.waveform(f,A,b,t0,tend,d_end_t,seed_number=n)
        # (needs arguments of matcharrays) 
        a = np.arange(A_low, A_hi+da, da)
        f = np.arange(f_low, f_hi+df, df)
        g = np.arange(gamma_low, gamma_hi+dg, dg)

        A = []
        F = []
        G = []
        M = []
        T = []

        Obj = Crosscor(dt, T_full, d)
        for h in a:
            for i in f:
                for j in g:
                    Obj.template(h, i, j, tmplt_dur)
                    t, m = match(Obj.d, Obj.y, Obj.dt)
                    T.append(t)
                    M.append(m)
                    A.append(h)
                    F.append(i)
                    G.append(j)
        
        # Here is all the results per waveform 
        Ts.append(T)
        Ms.append(M)
        As.append(A)
        Fs.append(F)
        Gs.append(G)
    
    A0 = As[0]
    f0 = Fs[0]
    g0 = Gs[0]
    t0 = np.array(Ts[0])
    
    M_array = np.array(Ms)
    combinedM = np.sum(M_array ** 2, axis = 0) ** .5

    M_M = []
    M_T = []
    step = 0

    for i in combinedM:
        x = np.argmax(i)
        M_T.append(t0[step,x])
        M_M.append(combinedM[step,x])
        step += 1
    
    max_Match = np.argmax(np.array(M_M))

    globA = A0[max_Match]
    globF = f0[max_Match]
    globG = g0[max_Match]
    globT = M_T[max_Match]
    globM = M_M[max_Match]

    #output = np.vstack((A,f1,g,T,M)).T
    #outputfile = "results/{}.txt".format(outputfile)
    #np.savetxt(outputfile, output, fmt="%f\t%f\t%f\t%f\t%f")

    return(globA, globF, globG, globT, globM)        
