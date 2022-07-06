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
        r_peaks.append(d[ii])

    m = int(np.median(r_peaks))
    return m


def trim_after_detection(detections, annotations, start_index, end_index):

    # start_index = annotated index to start at after trimming
    # end_index = annotated index to end at after trimming
    det_start_posn=int((annotations[start_index]+annotations[start_index-1])/2) # allow for detection half interval before start point annotation
    det_end_posn=int((annotations[end_index]+annotations[end_index+1])/2) # allow for detection half interval after end point annotation
    
    annotations_trimmed=annotations[start_index:(end_index+1)] # trim annotations to match trimmed detections
    detections_trimmed = detections[ (detections >= det_start_posn) & (detections <= det_end_posn) ] # remove detections with positions outwith range
    
    return detections_trimmed, annotations_trimmed
