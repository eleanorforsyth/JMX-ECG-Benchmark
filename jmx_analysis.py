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

a = 2 # number of annotated beats to trim from start
b = -2 # number of annotated beats to trim from end

norm_jmx = [2.991280672344543, 2.761072988147224, 4.480973175296319]

def det_delay(det_posn, anno_R):
    # Calculates the nearest difference and detect extra and missed detections
    diff={}
    len_det_posn=len(det_posn)
    delay=np.zeros(len_det_posn) # We have only one difference for each detected value
    
    for i in range (0,len_det_posn): # scan through detected peaks
        diff=det_posn[i] - anno_R # subtract ALL annotated values from the ith detection value
        index=np.abs(diff).argmin() # return the index of the smallest difference value
        delay[i]=diff[index] # store the value of that smallest difference in array 'val' at position 'i'
    delay=np.median(delay)
    print(delay)
    return delay


def trim_after_detection(detections, annotations, start_index, end_index):

    # start_index = annotated index to start at after trimming
    # end_index = annotated index to end at after trimming
    det_start_posn=int((annotations[start_index]+annotations[start_index-1])/2) # allow for detection half interval before start point annotation
    det_end_posn=int((annotations[end_index]+annotations[end_index+1])/2) # allow for detection half interval after end point annotation
    
    annotations_trimmed=annotations[start_index:(end_index+1)] # trim annotations to match trimmed detections
    detections_trimmed = detections[ (detections >= det_start_posn) & (detections <= det_end_posn) ] # remove detections with positions outwith range
    
    return detections_trimmed, annotations_trimmed




def mapping_curve():
    # equate mean point to cube root of 0.5 so that if all three parameters are average, when multiplied together we get 50% as an overall result
    for_x_is_1 = 0.5 ** (1. / 3) # i.e cube root of 0.5
    # To use individual parameters for comparison, use parameter^3, again to make 0.5=mean
    
    x = np.array([0.0, 1.0, 6.0, 10.0]) # 'source' input points for piecewise mapping
    # x = np.array([0.0, 1.0, 8.0, 15.0]) # alternative mapping for more detail in lower values.
    y = np.array([1.0, for_x_is_1, 0.2, 0.0]) # 'destination' output points for piecewise mapping
    z = np.polyfit(x, y, 3) # z holds the polynomial coefficients of 3rd order curve fit
    # i.e. z[0] = coeff of x^3, z[1] = coeff of x^2, z[2] = coeff of x, z[3] = constant
    # If piecewise mapping points are changed, check that polynomial approximation
    # curve is smooth and does not dip below zero - move 3rd input point (default 6.0)
    # if required.
    # To test plot 3rd order curve fit (uncomment to plot):
    return z


def normalise_and_map(param,norm):
    # Normalises and maps to benchmark value using 'poly' 3rd order polynomial
    poly=mapping_curve()
    param_norm = param / norm;
    if param_norm<=10.0: # Is normalised param less than 10x the global reference?
        x = param_norm    
        param_mapped=(poly[0]*x*x*x)+(poly[1]*x*x)+(poly[2]*x)+poly[3]
    else:
        param_mapped=0.0 # if the input parameter is greater than 10x the global
        # reference, a returned value of 0.0 will signify failure as a detector,
        # and when multiplied, the overall benchmark will also be 0
    
    return param_mapped


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


def evaluate(det_posn, anno_R, trim=True):
    """
    JMX analysis of interval variation, missed beat and extra detection positions
    det_posn: the timestamps of the detector in sample positions
    anno_R: the ground truth in samples
    """
    delay_correction = det_delay(det_posn, anno_R)
    det_posn = np.array(det_posn)-int(delay_correction) # Correction for detector delay
                    
    if trim==True:
        det_posn, anno_R = trim_after_detection(det_posn, anno_R, a, b)
        
    if len(det_posn)<=10:
        warning='WARNING: Less than ten detections'
        print(warning)
    
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
    
    return interval_differences_for_jitter, len(missed_beats), len(extra_beats)



def individual_score(jmx):
    """
    Takes a jmx and calculates the relative jmx benchmark of them between 0 and 1.
    """
    for i in range(len(jmx)):
        jmx = normalise_and_map(name_val, jmx_norm[i])

        
def total_score(jmx):
    jmx =  normalise_jmx(jmx)
    j = 1.0
    for i in jmx:
        j = j * i
    return j
