"""
JMX analysis
============
Analyses a detector for Jitter, Missed beats and eXtra/spurious detections.
At the heart is the jitter detection algorithm on annotation/detection
pairs. A Jitter of 0ms gives a 100% score, A jitter of 4ms 80% score
and drops to 0% for 40ms jitter.
Missed beats and extra/spurios detections are then used to calculate Accuracy
as a normalised measure independent on the number of beats.
The overall score is then: JMX = Jitter/% * Accuracy/%.
"""
import numpy as np
import util
from scipy import stats

# Used to determine how many beats could have been at max heartrate.
# This is needed to calculate the true negatives (TN).
maxHR = 220

# The standard jitter time at the moment for standard detectors.
norm_jitter = 12E-3 # sec

a = 10 # number of annotated beats to trim from start
b = -5 # number of annotated beats to trim from end

# keys for the jmx dict:
key_jitter = "jitter" # temproal jitter in s
key_accuracy = "accuracy" # statistical accuracy
key_jmx = "jmx" # jmx score
key_tp = "TP" # True positives
key_tn = "TN" # True negatives
key_fp = "FP" # False positives
key_fn = "FN" # False negatives

def mapping_jitter(x):
    # Normalises and maps to benchmark value using 'poly' 3rd order polynomial

    return 1/(x+1)


def nearest_diff(annotation, nearest_match):
    # Calculates the nearest difference between values in two arrays and saves
    # index and sample position of nearest
    
    len_annotation=len(annotation)

    # Array which will contain matches for all annotations no matter if the actual detection is missing
    used_indices=[] 
    for i in range(len_annotation): # scan through 'source' peaks 
        diff = nearest_match - annotation[i] # subtract ith source array value from ALL nearest match values
        index = np.abs(diff).argmin() # return the index of the smallest difference value
        used_indices.append((index,np.abs(nearest_match[index]-annotation[i]))) # save as tuple in used_indices

    # Array which will contain unique instances by always choosing the one which has the shortest time
    # difference.
    # 
    unique_diffs=[]
    index_used=[]
    for j in used_indices:
        if not (j[0] in index_used):
            uni = []
            for k in used_indices:
                if k[0] == j[0]:
                    uni.append(k)
            i = np.argmin(uni,0)[1]
            unique_diffs.append(uni[i][1])
            index_used.append(j[0])
        
    return unique_diffs


def score(jitter,accuracy):
    """
    Calculates the JMX score by multiplying the normalised jitter
    with the accuracy which in turn is based on missing and extra beats.
    """
    jitter_score = mapping_jitter(jitter / norm_jitter) # normalised jitter 0..1
    return accuracy * jitter_score


def evaluate(det_posn, anno_R, fs, nSamples, trim=True):
    """
    JMX analysis of interval variation, missed beat and extra detection positions
    det_posn: the timestamps of the detector in sample positions
    anno_R: the ground truth in samples
    fs: sampling rate of the ECG file
    nSamples: number of samples in the ECG file
    returns:
    jmx[key_jitter]   : jitter in s
    jmx[key_tp]       : true positive beats
    jmx[key_tn]       : true negative beats
    jmx[key_fp]       : false positive beats
    jmx[key_fn]       : false negative beats
    jmx[key_accuracy] : accuracy
    jmx[key_jmx]      : JMX Score
    """

    # Median delay of the detection against the annotations
    delay_correction = util.calcMedianDelay(det_posn, anno_R)

    # Correction for detector delay
    det_posn = np.array(det_posn)-int(delay_correction) 

    # Trims 1st and last detections
    if trim==True:
        det_posn, anno_R = util.trim_after_detection(det_posn, anno_R, a, b)

    # Do we have enough detections?
    if len(det_posn)<=10:
        warning='WARNING: Less than ten detections'
        print(warning)

    # Number of annotated R peaks
    len_anno_R = len(anno_R)

    # Number of detected R peaks
    len_det_posn = len(det_posn)

    # return anno / detector pairs
    anno_det_pairs = nearest_diff(anno_R, det_posn) 
    
    differences_for_jitter=[]
        
    for d in anno_det_pairs: 

        difference = np.abs(d / fs)
        differences_for_jitter.append(difference)

    jmx = {}

    jmx[key_jitter] = stats.median_absolute_deviation(differences_for_jitter)
    fp = len_det_posn - len(differences_for_jitter) # all detections - true positive = false positive
    fn = len_anno_R - len(differences_for_jitter) # all detections
    tp = len(differences_for_jitter)
    maxBeats = nSamples / fs * maxHR / 60
    tn = maxBeats - (tp + fn + fp) # remaining samples
    jmx[key_tp] = tp
    jmx[key_tn] = tn
    jmx[key_fp] = fp
    jmx[key_fn] = fn
    if (tp + tn + fp + fn) > 0:
        accuracy = (tp + tn)/(tp + tn + fp + fn)
        jmx[key_accuracy] = accuracy
        jmx[key_jmx] = score(jmx[key_jitter],accuracy)
    else:
        jmx[key_accuracy] = False
        jmx[key_jmx] = False
    print(jmx)
    return jmx
