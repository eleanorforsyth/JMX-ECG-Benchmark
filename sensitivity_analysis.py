#!/usr/bin/python3
"""
Trad sensitivity analysis with true positives, false positives, and false negatives.
"""
import numpy as np

"""
From the detected R peaks the function works itself backwards to
calculate the median delay the detector introduces. This is used
for benchmarking to compensate for different delays the detctors
introduce.
"""
def calcMedianDelay(detected_peaks, anno):

    r_peaks = []

    for i in detected_peaks:
        d = np.abs(i-anno)
        ii = np.argmin(d)
        # print(ii,d[ii],d)
        r_peaks.append(d[ii])

    m = int(np.median(r_peaks))
    return m


"""
The central function evaluating true positive, false positive and false negative.
"""
def evaluate(detected_peaks, annotation, tol):

    delay = calcMedianDelay(detected_peaks, annotation)

    detected_peaks = np.unique(detected_peaks)
    annotation = np.unique(annotation)
    
    TP = 0

    for anno_value in annotation:
        test_range = np.arange(anno_value-tol+delay, anno_value+1+tol+delay)
        in1d = np.in1d(test_range, detected_peaks)
        if np.any(in1d):
            TP = TP + 1

    FP = len(detected_peaks)-TP
    FN = len(annotation)-TP

    return TP, FP, FN
