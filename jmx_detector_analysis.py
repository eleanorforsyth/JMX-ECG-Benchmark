# -*- coding: utf-8 -*-
"""

This code will run all subjects, all experiments, all leads recordings through
all detectors or a single detector as required.
For each recording (for which there are annotations) passed through a detector
the detection locations will be saved, and then these passed for interval
analysis, where jitter, missed beats and extra/spurious detections are
identified. Jitter is taken as the difference (in samples) between the
annotated interval and the detected interval, and is not truly HRV as it is
calculated not just at rest.
For each recording (as above) passed through a detector, the jitter, missed
beat sample locations and extra/spurious detection locations are all saved as
seperate csv files. This means that all 'raw' interval analysis data is
available for subsequent benchmarking, plotting or analysis by lead type,
experiment, etc as desired and has not been combined in a way which results in
loss of information.
The matched filter detector can use a default template or user generated
averaged PQRST shapes for each subject.
The positional data for missed and extra beats can be used to plot markers
onto source recordings to indicate where and how a detector is behaving with
actual ECG recordings.
This code can be used to generate 'global' average standard deviatrion and
mean missed beat and extra detection values which are then saved as a csv file
and are used as reference values for the new overall benchmarking method

Create a folder named 'saved_csv' in the current directory to save csv files to.

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ecg_gudb_database import GUDb
from ecgdetectors import Detectors
import pathlib # For local file use


""" FUNCTIONS """

def det_delay(det_posn, anno_R):
    # Calculates the nearest difference and detect extra and missed detections
    diff={}
    len_det_posn=len(det_posn)
    delay=np.zeros(len_det_posn) # We have only one difference for each detected value
    
    for i in range (0,len_det_posn): # scan through detected peaks
        diff=anno_R - det_posn[i] # subtract ith detection value from ALL annotated values
        index=np.abs(diff).argmin() # return the index of the smallest difference value
        delay[i]=diff[index] # store the value of that smallest difference in array 'val' at position 'i'
    mean_delay=np.mean(delay)
    return mean_delay


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

def trim_after_detection(detections, annotations, start_index, end_index):

    # start_index = annotated index to start at after trimming
    # end_index = annotated index to end at after trimming
    det_start_posn=int((annotations[start_index]+annotations[start_index-1])/2) # allow for detection half interval before start point annotation
    det_end_posn=int((annotations[end_index]+annotations[end_index+1])/2) # allow for detection half interval after end point annotation
    
    annotations_trimmed=annotations[start_index:(end_index+1)] # trim annotations to match trimmed detections
    detections_trimmed = detections[ (detections >= det_start_posn) & (detections <= det_end_posn) ] # remove detections with positions outwith range
    
    return detections_trimmed, annotations_trimmed


def interval_analysis(det_posn, anno_R):
    # Main analysis of interval variation and missed beat/extra detection positions
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

#%%
"""
***************************
    START OF MAIN CODE   
***************************
"""

fs = 250 #sampling rate
detectors = Detectors(fs)# Initialise detectors for 250Hz sample rate (GUDB)

current_dir = pathlib.Path(__file__).resolve()

#%% Initialise parameters for analysis

save_global_results = True # when 'True' saves global jitter, missed, extra values as csv and prints

trim=True # CHANGE TO FALSE IF YOU DONT WANT TO TRIM
a = 10 # number of annotated beats to trim from start
b = -5 # number of annotated beats to trim from end
# * Values chosen by observation of detector settling intervals required *

#initialise for plots (*if used*)
plt.rc('xtick',labelsize=10)
plt.rc('ytick',labelsize=10)

analysed=0 # overall count of analysed subjects

data_det_lead = {} # initialise for data to be saved by lead and detector
det_lead_stats = {} # initialise for data to be saved by lead and detector

# initialise for global average jitter standard deviation, and global mean missed and extra beats
jitter_std_dev_total = 0.0 # running total for standard deviations for jitter for each detector

global_jitter=[]
global_missed=[]
global_extra=[]

# Detectors, recording leads and experiments can be added/removed from lists as required
all_detectors=['two_average_detector', 'swt_detector', 'engzee_detector', 'christov_detector', 'hamilton_detector', 'pan_tompkins_detector', 'matched_filter_detector']
all_recording_leads=["einthoven_ii", "chest_strap_V2_V1"] # can be expanded if required
all_experiments = ["sitting","maths","walking","hand_bike","jogging"]

for detector in all_detectors: # initialise empty arrays:
    
    name_jitter_sub = 'jitter_accum_sub_' + detector
    exec(name_jitter_sub+' = []') # initialse for all detector names, jitter arrays
    name_missed_sub = 'missed_accum_sub_' + detector
    exec(name_missed_sub + ' = []') # initialse for all detector names, missed beat arrays
    name_extra_sub = 'extra_accum_sub_' + detector
    exec(name_extra_sub + ' = []') # initialse for all detector names, missed detection arrays

    for record_lead in all_recording_leads: # loop for all chosen leads
        
        name_jitter = 'jitter_accum_' + record_lead +'_' + detector
        exec(name_jitter+' = []') # initialse for all rec leads, det names (jitter arrays)
        name_missed = 'missed_accum_' + record_lead +'_' + detector
        exec(name_missed + ' = []') # initialse for all rec leads, det names (missed beat arrays)
        name_extra = 'extra_accum_' + record_lead +'_' + detector
        exec(name_extra + ' = []') # initialse for all rec leads, det names (missed detection arrays)
        
        jitter_all_exp=[]
        missed_all_exp=[]
        extra_all_exp=[]
        
        for experiment in all_experiments: # loop for all chosen experiments
            
            jitter=[]
            missed=[]
            extra=[]
            
            for subject_number in range(0, 25): # loop for all subjects
                
                # print('')
                print("Analysing subject %d, %s, %s" %(subject_number, experiment, record_lead))
    
                # creating class which loads the experiment
        
                # For online GUDB access
                ecg_class = GUDb(subject_number, experiment) 
            
                # For local GUDB file access:
                # from ecg_gla_database import Ecg # For local file use
                # data_path = str(pathlib.Path(__file__).resolve().parent.parent/'experiment_data')
                # ecg_class = Ecg(data_path, subject_number, experiment)
                
                # getting the raw ECG data numpy arrays from class
                chest_strap_V2_V1 = ecg_class.cs_V2_V1
                einthoven_i = ecg_class.einthoven_I
                einthoven_ii = ecg_class.einthoven_II
                einthoven_iii = ecg_class.einthoven_III
        
                # getting filtered ECG data numpy arrays from class
                ecg_class.filter_data()
                chest_strap_V2_V1_filt = ecg_class.cs_V2_V1_filt
                einthoven_i_filt = ecg_class.einthoven_I_filt
                einthoven_ii_filt = ecg_class.einthoven_II_filt
                einthoven_iii_filt = ecg_class.einthoven_III_filt
            
                data=eval(record_lead) # set data array (i.e. recording to be processed)
               
                if 'chest' in record_lead: 
                    if ecg_class.anno_cs_exists:
                        data_anno = ecg_class.anno_cs
                        exist=True
                        analysed=analysed+1
                    else:
                        exist=False
                        print('')
                        print("No chest strap annotations exist for subject %d, %s exercise" %(subject_number, experiment))
                        print('')
                else:
                    if ecg_class.anno_cables_exists:
                        data_anno = ecg_class.anno_cables
                        exist=True
                        analysed=analysed+1
                    else:
                        exist=False
                        print('')
                        print("No cables annotations exist for subject %d, %s exercise" %(subject_number, experiment))
                        print('')
                        
#%% Detection
        
        ### Applying detector to each subject ECG data set then correct for mean detector
        # delay as referenced to annotated R peak position
        # Note: the correction factor for each detector doesn't need to be exact,
        # but centres the detection point for finding the nearest annotated match
        # It may/will be different for different subjects and experiments
        
        # If trim is True, start and end will be trimmed to avoid detector
        # artefacts such as missed beats and extra detections while settling
            
                if exist==True: # only proceed if an annotation exists
                    # form of detector call: 'detected_peaks = detectors.two_average_detector(data)'
                    if detector == 'matched_filter_detector':
                        template_name = ('templates/sitting_Sub'+str(subject_number)+'_'+record_lead+'_PQRSTave_44.csv')
                        detected_peaks = detectors.matched_filter_detector(data, template_name)
                    else:
                        detected_peaks = eval('detectors.'+detector+'(data)') # call detector class for current detector
                    delay_correction=det_delay(detected_peaks, data_anno) # fetch mean delay
                    detected_peaks_corr=np.array(detected_peaks)+int(delay_correction) # Correction for detector delay
                    
                    if trim==True:
                        detected_peaks_corr_trim, data_anno_trim = trim_after_detection(detected_peaks_corr, data_anno, a, b)
                    else: # keep naming consistent even if trim not applied
                        data_anno_trim = data_anno
                        detected_peaks_corr_trim = detected_peaks_corr
                    if len(detected_peaks_corr_trim)<=10:
                        print()
                        warning='WARNING: Less than ten detections while using '+detector+' to analyse subject no. '+str(subject_number)+', '+record_lead+', '+experiment+'.'
                        print(warning)
                        print()
                    interval_results = interval_analysis(detected_peaks_corr_trim, data_anno_trim) # perform interval based analysis
                    
                    jitter=np.concatenate((jitter, interval_results[0])) # jitter results
                    missed.append(len(interval_results[1])) # missed beat results
                    extra.append(len(interval_results[2])) # extra detection results
                    
                    print(detector + ' done')
                    
            # ^ LOOP AROUND FOR NEXT SUBJECT
            
            
            # Save to csv files, all subjects concatenated/appended
            # https://stackoverflow.com/questions/27126511/add-columns-different-length-pandas/33404243
            jitter_list=jitter.tolist()
            missed_beats=missed
            extra_detections=extra
            jitter_df = pd.DataFrame({"jitter":jitter_list}) # since length different
            missed_extra_df = pd.DataFrame({"missed_beats":missed_beats, "extra_detections":extra_detections})
            categories_df = pd.concat([jitter_df, missed_extra_df], axis=1)
            file_name='saved_csv/'+detector+'_'+record_lead+'_'+experiment+'.csv'
            categories_df.to_csv(file_name, index=True)
            #print(new.head())
            
            # Concatenate all experiment results, for all subjects
            exec(name_jitter + ' = np.concatenate((' + name_jitter + ', jitter_list))') # jitter results
            exec(name_missed + ' = np.concatenate((' + name_missed + ', missed))') # jitter results
            exec(name_extra + ' = np.concatenate((' + name_extra + ', extra))') # jitter results
                        
        # ^ LOOP AROUND FOR NEXT EXPERIMENT
        
        
        # Add data for analysis by lead to (full array) 'data_det_lead' dictionary
        
        exec('data_det_lead["'+name_jitter+'"] = '+name_jitter)
        exec('data_det_lead["'+name_missed+'"] = '+name_missed)
        exec('data_det_lead["'+name_extra+'"] = '+name_extra)
        
        # Add stats for analysis by lead to dictionaries (SD, missed and extra mean values)
        # MAD is 'np.mean(abs(x - np.mean(x)))', but for a detected
        # interval=anno interval we can used: 'np.mean(abs(x))'
        # since matched interval difference=0
        exec(name_jitter+'_mad = np.mean(abs('+name_jitter+'))')
        exec(name_missed+'_mean = np.mean(np.asarray('+name_missed+'))')
        exec(name_extra+'_mean = np.mean('+name_extra+')')
        # Add stats for analysis to 'det_lead_stats' dictionary (SD, missed and extra mean dictionaries)
        exec('det_lead_stats["'+name_jitter+'_mad"] = '+name_jitter+'_mad')
        exec('det_lead_stats["'+name_missed+'_mean"] = '+name_missed+'_mean')
        exec('det_lead_stats["'+name_extra+'_mean"] = '+name_extra+'_mean')
        
        
        # det_lead_stats
        
        # for global calculations:
        exec(name_jitter_sub+ '= np.concatenate(('+name_jitter_sub+', ' + name_jitter +'))')
        exec(name_missed_sub+ '= np.concatenate(('+name_missed_sub+', ' + name_missed +'))')
        exec(name_extra_sub+ '= np.concatenate(('+name_extra_sub+', ' + name_extra +'))')        
        
    # ^ LOOP AROUND FOR NEXT LEAD
    
    exec('global_jitter=np.concatenate((global_jitter, ' + name_jitter_sub +'))')
    exec('global_missed=np.concatenate((global_missed, ' + name_missed_sub +'))')
    exec('global_extra=np.concatenate((global_extra, ' + name_extra_sub +'))')
    
# END DETECTOR LOOP

# Save data as csv files:
# https://www.edureka.co/community/65139/valueerror-arrays-same-length-valueerror-arrays-same-length
det_lead_df = pd.DataFrame.from_dict(data_det_lead, orient='index')
det_lead_df=det_lead_df.transpose() # transpose colums and rows
file_name='saved_csv/data_det_lead.csv'
det_lead_df.to_csv(file_name)

df_stats = pd.DataFrame(det_lead_stats, index=[0])
file_name='saved_csv/det_lead_stats.csv'
df_stats.to_csv(file_name)

if save_global_results==True:   
    # global_jitter_mad = np.mean(abs(global_jitter - np.mean(global_jitter))) (by definition)
    global_jitter_mad = np.mean(abs(global_jitter)) # since matched interval difference =0
    print('')
    print('global_jitter_mad = ', global_jitter_mad)
    
    global_missed_mean=np.mean(global_missed)
    print('')
    print('global_missed_mean = ', global_missed_mean)
    
    global_extra_mean=np.mean(global_extra)
    print('')
    print('global_extra_mean = ', global_extra_mean)
    
    # Save Global mean benchmark reference values 
    global_data={"Global_jitter_mad":global_jitter_mad, "Global_missed_mean":global_missed_mean, "Global_extra_mean":global_extra_mean}
    global_df = pd.DataFrame(global_data, index=[0])
    global_df.to_csv('saved_csv/global_results.csv')







