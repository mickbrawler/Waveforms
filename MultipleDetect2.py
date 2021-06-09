import json
import numpy as np
import pylab as pl
from numba import jit

def multidetect(file1,file2):
    
    """
    PARAMETERS:
    -----------
    file1: (json) First waveform json file
    file2: (json) Second waveform json file
    """

    with open(file1, "r") as f:
        data = json.load(f)
    f1 = data["f"]
    g = data["g"]
    t = data["t"]
    M1 = data["m"]

    with open(file2, "r") as f:
        data = json.load(f)
    M2 = data["m"]
    
    # Method of finding common peaks
    # Currently gives no leeway to index of common peaks across detectors, only leeway to peaks' match values
    peakIndice0 = np.array([])
    peakIndice = np.array([])
    peakValue = np.array([])

    for m in M1:
        i = 0
        j = 0
        potentialIndice0 = np.array([])
        potentialIndice = np.array([])
        potentialValue = np.array([])
        
        while len(m[i:]) >= 3: # Don't let "window" go beyond match array
            if m[i] < m[i+1] & m[i+1] > m[i+2]:
                if (M2[j,i+1] > m[i+1] - 10 and M2[j,i+1] < m[i+1] + 10) & (M2[j,i] < M2[j,i+1] & M2[j,i+1] > M2[j,i+2]):
                #  Gives leeway to compatible peak values. See if peak is present in 2nd detector in same index as first
                    potentialIndice0 = np.append(potentialIndice0,j)
                    potentialIndice = np.append(potentialIndice,i+1)
                    potentialValue = np.append(potentialValue,m[i+1])
            i += 1
        
        if len(potentialValue) > 0:
            maxIndexm = np.argmax(potentialValue)
            peakIndice0 = np.append(peakIndice0,potentialIndice0[maxIndexm])
            peakIndice = np.append(peakIndice,potentialIndice[maxIndexm])
            peakValue = np.append(peakValue,potentialValue[maxIndexm])
            
        j += 1

    # Then get maximum of peak Value and its indexes, indexes can be used to find the time of the peak maximum. 
    # Index0 should function as index for frequency and gamma.

    if len(peakValue) > 0:
        maxIndexM = np.argmax(peakValue)
        MaxpeakIndice0 = peakIndice0[maxIndexM]
        MaxpeakIndice = peakIndice[maxIndexM]
        MaxpeakValue =  peakValue[maxIndexM] 
    
    return(f[MaxpeakIndice0],g[MaxpeakIndice0],t[MaxpeakIndice0,MaxpeakIndice],MaxpeakValue)
