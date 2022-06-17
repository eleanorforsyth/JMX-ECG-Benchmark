#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Calculates the individual jitters, accuracy and JMX values.
Outputs the averages and stores the inidividual ones in a file.
"""

import sys
import os
import numpy as np
import json
from ecg_gudb_database import GUDb
from ecgdetectors import Detectors
import pathlib # For local file use
from multiprocessing import Process

# The JMX analysis for a detector
import jmx_analysis

# directory where the results are stored
resultsdir = "results"

try:
    os.mkdir(resultsdir)
except OSError as error:
    pass

fs = 250 #sampling rate

detectors = Detectors(fs) # Initialise detectors for 250Hz sample rate (GUDB)

current_dir = pathlib.Path(__file__).resolve()

recording_leads = "einthoven_ii"
experiment = "sitting"

jmx_results = np.empty((0,2))

f = open("norm_calc.tsv","w")

for detector in detectors.detector_list:

    detectorname = detector[1].__name__
    detectorfunc = detector[1]
    
    print("Processing:",detector[0])

    for subject_number in range(0, 25): # loop for all subjects

        print("Analysing subject {}, {}, {}, {}".format(subject_number, experiment, recording_leads, detector[0]))

        # creating class which loads the experiment

        # For online GUDB access
        ecg_class = GUDb(subject_number, experiment) 

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

        data = einthoven_ii

        if ecg_class.anno_cables_exists:
            data_anno = ecg_class.anno_cables
            exist=True

        #%% Detection

        ### Applying detector to each subject ECG data set then correct for mean detector
        # delay as referenced to annotated R peak position
        # Note: the correction factor for each detector doesn't need to be exact,
        # but centres the detection point for finding the nearest annotated match
        # It may/will be different for different subjects and experiments

        if exist==True: # only proceed if an annotation exists
            detected_peaks = detectorfunc(data) # call detector class for current detector
            interval_results = jmx_analysis.evaluate(detected_peaks, data_anno, fs, len(data)) # perform interval based analysis
            jmx = np.array([interval_results[jmx_analysis.key_jitter],
                            interval_results[jmx_analysis.key_accuracy],
                            
            ])
            jmx_results = np.vstack( (jmx_results,jmx) )
            jmx_avg = np.average(jmx_results,axis=0)
            s = jmx_analysis.score(jmx_avg[0],jmx_avg[1])
            print("J = {:1.4f} sec, A = {:1.4f}, JMX = {:1.4f}".format(jmx[0],jmx[1],s))
            f.write("{}\t{}\t{}\t{}\n".format(jmx[0],jmx[1])
            f.flush()
print("FINAL: J = {:1.4f} sec, A = {:1.4f}".format(jmx_avg[0],jmx_avg[1]))
f.close()
