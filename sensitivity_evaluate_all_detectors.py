#!/usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
import json
from ecg_gudb_database import GUDb
from ecgdetectors import Detectors
import pathlib # For local file use
from multiprocessing import Process

# The JMX analysis for a detector
import sensitivity_analysis

# directory where the results are stored
resultsdir = "results"

try:
    os.mkdir(resultsdir)
except OSError as error:
    pass

fs = 250 #sampling rate

detectors = Detectors(fs) # Initialise detectors for 250Hz sample rate (GUDB)

current_dir = pathlib.Path(__file__).resolve()

# Detectors, recording leads and experiments can be added/removed from lists as required
all_recording_leads=["einthoven_ii", "chest_strap_V2_V1"] # can be expanded if required
all_experiments = ["sitting","maths","walking","hand_bike","jogging"]

def evaluate_detector(detector):

    detectorname = detector[1].__name__
    detectorfunc = detector[1]
    
    print("Processing:",detector[0])

    analysed=0 # overall count of analysed subjects

    sens_leads = {} # initialise for data to be saved by lead and detector

    for record_lead in all_recording_leads: # loop for all chosen leads
        
        sens_experiments = {}
        
        for experiment in all_experiments: # loop for all chosen experiments
            
            sens_subjects=[]
            
            for subject_number in range(0, 25): # loop for all subjects
                
                print("Analysing subject {}, {}, {}, {}".format(subject_number, experiment, record_lead, detector[0]))
    
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
                        print("No chest strap annotations exist for subject %d, %s exercise" %(subject_number, experiment))
                else:
                    if ecg_class.anno_cables_exists:
                        data_anno = ecg_class.anno_cables
                        exist=True
                        analysed=analysed+1
                    else:
                        exist=False
                        print("No cables annotations exist for subject %d, %s exercise" %(subject_number, experiment))
                        
                #%% Detection
        
                ### Applying detector to each subject ECG data set then correct for mean detector
                # delay as referenced to annotated R peak position
                # Note: the correction factor for each detector doesn't need to be exact,
                # but centres the detection point for finding the nearest annotated match
                # It may/will be different for different subjects and experiments

                if exist==True: # only proceed if an annotation exists
                    detected_peaks = detectorfunc(data) # call detector class for current detector
                    interval_results = sensitivity_analysis.evaluate(detected_peaks, data_anno, fs/10) # perform interval based analysis
                    sens_subjects.append(interval_results)
                    
            # ^ LOOP AROUND FOR NEXT SUBJECT

            sens_experiments[experiment] = sens_subjects
                        
        # ^ LOOP AROUND FOR NEXT EXPERIMENT
        
        # Add data for analysis by lead to (full array) 'data_det_lead' dictionary
        
        sens_leads[record_lead] = sens_experiments
        
    # ^ LOOP AROUND FOR NEXT LEAD
    serialized_data = json.dumps(sens_leads,indent="\t")
    f = open(resultsdir+"/sens_"+detectorname+".json","w")
    f.write(serialized_data)
    f.close


if (len(sys.argv)>1):
    evaluate_detector(detectors.detector_list[int(sys.argv[1])])
else:
    for detector in detectors.detector_list:
        pEvalDet = Process(target=evaluate_detector, args=(detector,))
        pEvalDet.start()
