import json
import numpy as np
import pylab as pl
from numba import jit

# Original bare-bones approach to Multiple Detector Problem
def multidetect(file1, file2, threshold, forgive):
    
    """
    METHOD: Takes files, a match threshold, and a forgiveness value to
    find common peaks accross each "detector's" match arrays

    PARAMETERS:
    -----------
    file1: (json) First waveform json file
    file2: (json) Second waveform json file

    OUTPUT: Returns the global maximum values of frequency, gamma, time,
    and the corresponding match value.
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
    
    M1 = np.array(M1)
    M2 = np.array(M2)

    M1 = (M1 > threshold) * M1
    M2 = (M2 > threshold) * M2

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
            if m[i] < m[i+1] and m[i+1] > m[i+2]:
                if M2[j,i+1] > m[i+1] - forgive and M2[j,i+1] < m[i+1] + forgive and M2[j,i] < M2[j,i+1] and M2[j,i+1] > M2[j,i+2]:
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

def window(file1,file2):
    """
    METHOD
    ======
    Window of set size will slide through data and append the maximum 
    match in it to a list. For the corresponding windows in each detector,
    their maximums will be squared, summed, and square rooted. Index of 
    max value must be utilized to obtain global maximums: frequency, 
    gamma, and time.

    PARAMETERS
    ==========
    file1: (json) First waveform json file
    file2: (json) Second waveform json file

    OUTPUT
    ======
    Returns global maximum values of frequency, gamma, time, and match
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
    
    M1 = np.array(M1)
    M2 = np.array(M2)
    
    i = 0
    j = 0
    windMatches = np.array([])
    MaxPeak1 = np.array([])
    MaxPeak2 = np.array([])
    MatchListMaxs = np.array([])
    jIndices = np.array([])
    iIndices = np.array([])
    for m in M1:
        
        while True:
            if len(m[i:i+10]) > len(m[i:]):
                # Alteration of window size
                M1Index = np.where(m == np.amax(m[i:i+len(m[i:])]))
                M2Index = np.where(M2[j] == np.amax(M2[j,i:i+len(M2[j,i:])]))
                
                iIndices = np.append(iIndices,i)
                jIndices = np.append(jIndices,j)

                MaxPeak1 = np.append(MaxPeak1,m[M1Index])
                MaxPeak2 = np.append(MaxPeak2,M2[j,M2Index])

                operation = ((m[M1Index] ** 2) + (M2[j,M2Index] ** 2)) ** .5
                windMatches = np.append(windMatches, operation)
                break
            
            else:

                M1Index = np.where(m == np.amax(m[i:i+10]))
                M2Index = np.where(M2[j] == np.amax(M2[j,i:i+10]))

                iIndices = np.append(iIndices,i)
                jIndices = np.append(jIndices,j)

                MaxPeak1 = np.append(MaxPeak1,m[M1Index])
                MaxPeak2 = np.append(MaxPeak2,M2[j,M2Index])

                operation = ((m[M1Index] ** 2) + (M2[j,M2Index] ** 2)) ** .5
                windMatches = np.append(windMatches, operation)

                i += 10
        
        MaxIndex = np.argmax(windMatches) # index of max oper result gives the best window, which gives best peaks
        if MaxPeak1[MaxIndex] > MaxPeak2[MaxIndex]:
            MatchListMaxs = np.append(MatchListMaxs,MaxPeak1[MaxIndex])
        else:
            MatchListMaxs = np.append(MatchListMaxs,MaxPeak2[MaxIndex])
        j += 1

    return(np.max(MatchListMaxs))
