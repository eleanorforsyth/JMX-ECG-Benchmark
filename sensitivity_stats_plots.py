#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from ecgdetectors import Detectors
import json

experiment_names = ['sitting','maths','walking','hand_bike','jogging']

einth = 'einthoven_ii'
cs = 'chest_strap_V2_V1'

detectors = Detectors()
det_names = [i[1].__name__ for i in detectors.get_detector_list()]
plot_names = [i[0] for i in detectors.get_detector_list()]

resultsdir = "results"

alpha = 0.05

min_sens = 90 # %

def get_sensitivities(detector_name, leads, experiment):
    f = open(resultsdir+"/sens_"+detector_name+".json","r")
    js = f.read()
    data = json.loads(js)
    s = [i[0] for i in data[leads][experiment]]
    return np.array(s)


def get_result(det, leads, experiment):
    
    m = []
    s = []
    for det in det_names:
        print(det,experiment,get_sensitivities(det, leads, experiment))
        m.append(np.mean(get_sensitivities(det, leads, experiment)))
        s.append(np.std(get_sensitivities(det, leads, experiment)))

    return np.array(m),np.array(s)

def print_stat(p):
    if p == None:
        print('--- & ',end='')
        return
    s = ""
    if p < alpha:
        s = "*"
    print('{:03.2f}{} & '.format(p,s),end='')

    
def calc_stats(leads, experiment):
    print("Stats:",leads, experiment)
    print("      & ",end='')
    for det1 in det_names:
        print(det1," & ",end='')
    print("\\\\")
    for det1 in det_names:
        r1 = get_sensitivities(det1, leads, experiment)
        t,p = stats.ttest_1samp(r1,min_sens,alternative='greater')
        print_stat(p)
    print()

    


def double_plot(data1, std1, data2, std2, y_label, legend1, legend2, title=None):
    fig, ax = plt.subplots()
    x_pos = np.arange(len(plot_names))

    fig.set_size_inches(10, 7)
    width = 0.4
    rects1 = ax.bar(x_pos, data1, width, yerr=std1, alpha=0.5, ecolor='black', capsize=10)
    rects2 = ax.bar(x_pos+width, data2, width, yerr=std2, alpha=0.5, ecolor='black', capsize=10)
    ax.set_ylim([0,150])
    ax.set_ylabel(y_label)
    ax.set_xlabel('Detector')
    ax.set_xticks(x_pos + width / 2)
    ax.set_xticklabels(plot_names)
    ax.legend((rects1[0], rects2[0]), (legend1, legend2))

    if title!=None:
        ax.set_title(title)

    plt.tight_layout()

    return rects1, rects2

def print_result(title,data,std,legend):
    print("Sensitivities:",title)
    for i in zip(legend,data,std):
        print("{}: {:1.1f}+/-{:1.1f}".format(i[0],i[1],i[2]))
    print()

cs_sitting_avg,cs_sitting_std = get_result(det_names, cs, 'sitting')
einthoven_sitting_avg,einthoven_sitting_std = get_result(det_names, einth, 'sitting')

cs_jogging_avg,cs_jogging_std = get_result(det_names, cs, 'jogging')
einthoven_jogging_avg,einthoven_jogging_std = get_result(det_names, einth, 'jogging')

print()

print_result('sitting Einthoven',einthoven_sitting_avg,einthoven_sitting_std,det_names)
print_result('jogging Einthoven',einthoven_jogging_avg,einthoven_jogging_std,det_names)

print_result('sitting chest strap',cs_sitting_avg,cs_sitting_std,det_names)
print_result('jogging chest strap',cs_jogging_avg,cs_jogging_std,det_names)

double_plot(einthoven_sitting_avg, einthoven_sitting_std,
            einthoven_jogging_avg,einthoven_jogging_std,
            'Sensitivity (%)', 'Sitting', 'Jogging', 'Einthoven')


double_plot(cs_sitting_avg, cs_sitting_std,
            cs_jogging_avg, cs_jogging_std,
            'Sensitivity (%)', 'Sitting', 'Jogging', 'Chest strap')


calc_stats(einth,"sitting")
calc_stats(einth,"jogging")

print()

calc_stats(cs,"sitting")
calc_stats(cs,"jogging")




plt.show()
