#!/usr/bin/python3
"""
Trad sensitivity analysis with true positives, false positives, and false negatives.
"""
import numpy as np
import util

"""
The central function evaluating true positive, false positive and false negative.
"""
def evaluate(detected_peaks, annotation, tol):

    delay = util.calcMedianDelay(detected_peaks, annotation)

    detected_peaks = np.unique(detected_peaks)
    annotation = np.unique(annotation)
    
    tp = 0

    for anno_value in annotation:
        test_range = np.arange(anno_value-tol+delay, anno_value+1+tol+delay)
        in1d = np.in1d(test_range, detected_peaks)
        if np.any(in1d):
            tp = tp + 1

    fp = len(detected_peaks)-tp
    fn = len(annotation)-tp

    sensitivity = False
    
    if (tp + fn) > 0:
        sensitivity = tp/(tp+fn)*100.0

    return (sensitivity, tp, fp, fn)
