"""
JMX analysis
============
Analyses a detector for interval
analysis, where Jitter, Missed beats and eXtra/spurious detections are
identified. Jitter is taken as the difference (in samples) between the
annotated interval and the detected interval, and is not truly HRV as it is
calculated not just at rest.
"""

import numpy as np

def nearest_diff(source_array, nearest_match):
    # Calculates the nearest difference between values in two arrays and saves
    # index and sample position of nearest
    diff=[] # Temporary working array for all differences
    used_indices=[] # Arraywhich will contain only one instance of nearest matches
    
    len_source_array=len(source_array)
    nearest=np.zeros((len_source_array,3)) # store nearest index and position
    
    for i in range (0,len_source_array): # scan through 'source' peaks 
        diff = nearest_match - source_array[i] # subtract ith source array value from ALL nearest match values
        index=np.abs(diff).argmin() # return the index of the smallest difference value
        nearest[i, 0] = index # store the value of that smallest difference in array 'nearest' at position 'i'
        nearest[i, 1] = nearest_match[index]
        nearest[i,2] = source_array[i] # store actual 'source' position in nearest at position i
        if i==0 or index> nearest[i-1, 0]: # Eliminate any multiple matches
            used_indices.append((index, nearest_match[index])) # save as tuple in used_indices
        
    return nearest, used_indices

def evaluate(det_posn, anno_R):
    """
    JMX analysis of interval variation, missed beat and extra detection positions
    det_posn: the timestamps of the detector in sample positions
    anno_R: the ground truth in samples
    """
    len_anno_R = len(anno_R)

    all_anno = np.zeros((len_anno_R,2)) # store for nearest anno index and position
    
    all_anno[:, 0] = np.linspace(0, len_anno_R-1, len_anno_R)
    all_anno[:,1] = anno_R
    all_anno=all_anno.astype(int)
        
    len_det_posn = len(det_posn)
        
    nearest_anno, used_anno = nearest_diff(det_posn, anno_R) # return nearest anno index and position
    nearest_det, used_det = nearest_diff(anno_R, det_posn) # return nearest anno index and position
    
    """ MISSED BEAT IDENTIFICATION """
    # Generate array of unpaired annotated beats = missed beats
    unused_anno=all_anno # make 'working copy'
    for x in reversed(used_anno): # work backwards, deleting used annotations
        index1 = x[0] # index of used anno is first in tuple
        unused_anno = np.delete(unused_anno, index1, 0)
    
    """ EXTRA BEAT - IDENTIFICATION """
    # Find 'unpaired' detections which are NOT the nearest match to annotated positions
    extra_det_posn = det_posn # starting point - now remove any paired detection positions
    
    for i in range(len_det_posn-1, -1, -1): # work backwards through extra_det_posn deleting already 'paired' detections
        # print(i)
        for x in used_det:
            if det_posn[i] == x[1]: # if detected position already 'used'
                extra_det_posn = np.delete(extra_det_posn, i, 0)
          
    """ TEMPORAL ANALYSIS SECTION """
    interval_differences_for_jitter=[]
        
    for i in range(1, len(used_anno)): 
       
        valid_interval_anno=int(used_anno[i][1] - used_anno[i-1][1])    
        valid_interval_det = (nearest_det[(used_anno[i][0])][1]) - (nearest_det[(used_anno[i-1][0])][1])
        
        difference=valid_interval_det-valid_interval_anno
        interval_differences_for_jitter.append(difference)

    missed_beats = unused_anno #for clarity
    extra_beats = extra_det_posn #for clarity
    
    return interval_differences_for_jitter, missed_beats, extra_beats

