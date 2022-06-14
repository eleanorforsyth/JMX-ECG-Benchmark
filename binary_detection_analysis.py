#!/usr/bin/python3
"""
This script benchmarks all detectors with the GU database.
It runs it for both the chest strap and Einthoven II with loose cables.
"""
import pathlib
from multiprocessing import Process
import numpy as np
import pandas as pd
import os
from datetime import datetime
from ecgdetectors import Detectors
from ecg_gudb_database import GUDb

fs = 250 #sampling rate
resultsdir = 'results'

"""
From the detected R peaks the function works itself backwards to
calculate the median delay the detector introduces. This is used
for benchmarking to compensate for different delays the detctors
introduce.
"""
def calcMedianDelay(detected_peaks, anno, search_samples):

    r_peaks = []
    window = int(search_samples)

    for i in detected_peaks:
        d = np.abs(i-anno)
        ii = np.argmin(d)
        # print(ii,d[ii],d)
        r_peaks.append(d[ii])

    m = int(np.median(r_peaks))
    return m


def evaluate_detector(test, annotation, delay, tol=0):

    test = np.unique(test)
    reference = np.unique(annotation)
    
    TP = 0

    for anno_value in test:
        test_range = np.arange(anno_value-tol-delay, anno_value+1+tol-delay)
        in1d = np.in1d(test_range, reference)
        if np.any(in1d):
            TP = TP + 1
    
    FP = len(test)-TP
    FN = len(reference)-TP 

    return TP, FP, FN


def get_time():

        time = str(datetime.now().time())
        time = time[:5]
        time = time.replace(':', '.')
        
        return time


class Binary_test:
    """
    This class benchmarks detectors against the GU database.
    You need to install the API for the GUDB: https://github.com/berndporr/ECG-GUDB
    """
    
    def single_classifier_test(self, detector, tolerance, config="chest_strap"):

        max_delay_in_samples = fs / 3

        total_subjects = GUDb.total_subjects

        results = np.zeros((total_subjects, (4*len(GUDb.experiments))+1), dtype=int)

        for subject_number in range(0, total_subjects):
            progress = int(subject_number/float(total_subjects)*100.0)
            print("GUDB "+config+" progress: %i%%" % progress)

            results[subject_number, 0] = subject_number
            exp_counter = 1
            for experiment in GUDb.experiments:
                
                ecg_class = GUDb(subject_number, experiment)

                anno_exists = False
                if config=="chest_strap" and ecg_class.anno_cs_exists:
                    unfiltered_ecg = ecg_class.cs_V2_V1                   
                    anno = ecg_class.anno_cs
                    anno_exists = True
                elif config=="loose_cables" and ecg_class.anno_cables_exists:
                    unfiltered_ecg = ecg_class.einthoven_II 
                    anno = ecg_class.anno_cables
                    anno_exists = True
                elif config!="chest_strap" and config!="loose_cables":
                    raise RuntimeError("Config argument must be chest_strap or loose_cables!")
                    return results

                if anno_exists:                  

                    r_peaks = detector(unfiltered_ecg)

                    delay = calcMedianDelay(r_peaks, anno, max_delay_in_samples)
                    print("delay = ",delay)

                    # there must be a delay in all cases so anything below is a bad sign
                    if delay > 1:

                        TP, FP, FN = evaluate_detector(r_peaks, anno, delay, tol=tolerance)
                        TN = len(unfiltered_ecg)-(TP+FP+FN)

                        results[subject_number, exp_counter] = TP
                        results[subject_number, exp_counter+1] = FP
                        results[subject_number, exp_counter+2] = FN
                        results[subject_number, exp_counter+3] = TN

                exp_counter = exp_counter+4

        return results


    def classifer_test_all(self, tolerance, config="chest_strap"):

        try:
            os.mkdir(resultsdir)
        except OSError as error:
            pass

        output_names = ['TP', 'FP', 'FN', 'TN']

        detectors = Detectors(fs)

        total_results = np.zeros((GUDb.total_subjects, 4*len(GUDb.experiments)*len(detectors.detector_list)), dtype=int)

        counter = 0
        for det in detectors.detector_list:

            print('\n'+config+" "+det[0]+":")

            result = self.single_classifier_test(det[1], tolerance=tolerance, config=config)
            result = result[:, 1:]

            total_results[:, counter:counter+(4*len(GUDb.experiments))] = result

            counter = counter+(4*len(GUDb.experiments))        

        index_labels = np.arange(GUDb.total_subjects)
        col_labels = []

        for det in detectors.detector_list:
            for experiment_name in GUDb.experiments:
                for output_name in output_names:
                    label = det[1].__name__+" "+experiment_name+" "+output_name
                    col_labels.append(label)

        total_results_pd = pd.DataFrame(total_results, index_labels, col_labels, dtype=int)            
        total_results_pd.to_csv('results/binary'+config+'.csv', sep=',')

        return total_results_pd



def run_tests(leads):
    gu_test = Binary_test()
    gu_test.classifer_test_all(tolerance=(fs/10), config=leads)

pgustrap = Process(target=run_tests, args=('chest_strap',))
pgustrap.start()
pgucables = Process(target=run_tests, args=('loose_cables',))
pgucables.start()
    
pgustrap.join()
pgucables.join()
