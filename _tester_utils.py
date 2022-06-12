import numpy as np
from ecgdetectors import Detectors
from datetime import datetime

"""
From the detected R peaks the function works itself backwards to
calculate the median delay the detector introduces. This is used
for benchmarking to compensate for different delays the detctors
introduce.
"""
def calcMedianDelay(detected_peaks, unfiltered_ecg, search_samples):

    r_peaks = []
    window = int(search_samples)

    for i in detected_peaks:
        i = int(i)
        if i<window:
            section = unfiltered_ecg[:i]
            r_peaks.append(np.argmax(section))
        else:
            section = unfiltered_ecg[i-window:i]
            r_peaks.append(window - np.argmax(section))

    m = int(np.median(r_peaks))
    #print("Delay = ",m)
    return m


def evaluate_detector(test, annotation, delay, tol=0):

    test = np.unique(test)
    reference = np.unique(annotation)
    
    TP = 0

    for anno_value in test:
        test_range = np.arange(anno_value-tol-delay, anno_value+1+tol-delay)
        in1d = np.in1d(test_range, reference)
        if np.any(in1d):
            TP = TP + 1
    
    FP = len(test)-TP
    FN = len(reference)-TP 

    return TP, FP, FN


def get_time():

        time = str(datetime.now().time())
        time = time[:5]
        time = time.replace(':', '.')
        
        return time
