"""
JMX analysis
============
Analyses a detector for Jitter, Missed beats and eXtra/spurious detections.
Missed beats and extra/spurios detections are then used to calculate precision
as a normalised measure independent on the number of beats.
"""
import numpy as np
import util
from scipy import stats

# Used to determine how many beats could have been at max heartrate.
# This is needed to calculate the true negatives (TN).
maxHR = 220

# Cassirame et al. 2019 requires a jitter of less than 4ms for
# a detector which can analyse HRV. With the mapping curve below
# that gives then an 80% rating for the jitter. Of course jitter
# of zero is best and gives 100%.
norm_jitter = 4E-3 # sec

a = 2 # number of annotated beats to trim from start
b = -2 # number of annotated beats to trim from end

# keys for the jmx dict:
key_jitter = "jitter" # temproal jitter in s
key_accuracy = "accuracy" # statistical accuracy
key_jmx = "jmx" # jmx score
key_tp = "TP" # True positives
key_tn = "TN" # True negatives
key_fp = "FP" # False positives
key_fn = "FN" # False negatives

def mapping_curve():
    ratio_is_1 = 0.8 # gives 80% rating at norm_jitter
    x = np.array([0.0, 1.0, 6.0, 10.0]) # 'source' input points for piecewise mapping
    y = np.array([1.0, ratio_is_1, 0.2, 0.0]) # 'destination' output points for piecewise mapping
    z = np.polyfit(x, y, 3) # z holds the polynomial coefficients of 3rd order curve fit
    # i.e. z[0] = coeff of x^3, z[1] = coeff of x^2, z[2] = coeff of x, z[3] = constant
    # If piecewise mapping points are changed, check that polynomial approximation
    # curve is smooth and does not dip below zero - move 3rd input point (default 6.0)
    # if required.
    # To test plot 3rd order curve fit (uncomment to plot):
    return z


def mapping_jitter(x):
    # Normalises and maps to benchmark value using 'poly' 3rd order polynomial
    poly=mapping_curve()
    
    if x<=10.0: # Is normalised param less than 10x the global reference?
        j_mapped=(poly[0]*x*x*x)+(poly[1]*x*x)+(poly[2]*x)+poly[3]
    else:
        j_mapped=0.0 # if the input parameter is greater than 10x the global
        # reference, a returned value of 0.0 will signify failure as a detector,
        # and when multiplied, the overall benchmark will also be 0

    return j_mapped


def nearest_diff(source_array, nearest_match):
    # Calculates the nearest difference between values in two arrays and saves
    # index and sample position of nearest
    diff=[] # Temporary working array for all differences
    used_indices=[] # Array which will contain only one instance of nearest matches
    
    len_source_array=len(source_array)
    nearest=np.zeros((len_source_array,3)) # store nearest index and position
    
    for i in range (0,len_source_array): # scan through 'source' peaks 
        diff = nearest_match - source_array[i] # subtract ith source array value from ALL nearest match values
        index=np.abs(diff).argmin() # return the index of the smallest difference value
        nearest[i, 0] = index # store the value of that smallest difference in array 'nearest' at position 'i'
        nearest[i, 1] = nearest_match[index]
        nearest[i, 2] = source_array[i] # store actual 'source' position in nearest at position i
        if i==0 or index> nearest[i-1, 0]: # Eliminate any multiple matches
            used_indices.append((index, nearest_match[index])) # save as tuple in used_indices
        
    return nearest, used_indices


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

    # return nearest anno index and position
    nearest_anno, used_anno = nearest_diff(det_posn, anno_R)

    # return nearest anno index and position
    nearest_det, used_det = nearest_diff(anno_R, det_posn) 
    
    """ MISSED BEAT IDENTIFICATION """
    # Generate array of unpaired annotated beats = missed beats
    # First fill up the array with all annotation positions
    unused_anno = np.arange(len_anno_R)
    unused_anno = unused_anno.astype(int)
    for x in reversed(used_anno): # work backwards, deleting used annotations
        index1 = x[0] # index of used anno is first in tuple
        unused_anno = np.delete(unused_anno, index1)
    
    """ EXTRA BEAT - IDENTIFICATION """
    # Find 'unpaired' detections which are NOT the nearest match to annotated positions
    # starting point - now remove any paired detection positions
    extra_det_posn = det_posn 

    # work backwards through extra_det_posn deleting already 'paired' detection
    for i in range(len_det_posn-1, -1, -1): 
        # print(i)
        for x in used_det:
            if det_posn[i] == x[1]: # if detected position already 'used'
                extra_det_posn = np.delete(extra_det_posn, i, 0)
          
    """ TEMPORAL ANALYSIS SECTION """
    differences_for_jitter=[]
        
    for i in range(len(used_anno)): 

        # ground truth
        valid_anno = int(used_anno[i][1])

        # detection position
        valid_det = nearest_det[(used_anno[i][0])][1]
        
        difference = np.abs((valid_det - valid_anno) / fs)
        differences_for_jitter.append(difference)

    missed_beats = unused_anno #for clarity
    extra_beats = extra_det_posn #for clarity

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
    return jmx
