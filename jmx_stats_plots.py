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

def get_jmx(detector_name, leads, experiment):
    f = open(resultsdir+"/jmx_"+detector_name+".json","r")
    js = f.read()
    data = json.loads(js)
    s = []
    for i in data[leads][experiment]:
        if i["jmx"]:
            s.append(i["jmx"]*100)
    return np.array(s)


def get_result(det, leads, experiment):
    
    m = []
    s = []
    for det in det_names:
        print(det,experiment,get_jmx(det, leads, experiment))
        m.append(np.mean(get_jmx(det, leads, experiment)))
        s.append(np.std(get_jmx(det, leads, experiment)))

    return np.array(m),np.array(s)


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


cs_sitting_avg,gudb_cs_sitting_std = get_result(det_names, cs, 'sitting')
cable_sitting_avg,gudb_cable_sitting_std = get_result(det_names, einth, 'sitting')

cs_jogging_avg,gudb_cs_jogging_std = get_result(det_names, cs, 'jogging')
cable_jogging_avg,gudb_cable_jogging_std = get_result(det_names, einth, 'jogging')


double_plot(cs_sitting_avg, gudb_cs_sitting_std,
            cable_sitting_avg, gudb_cable_sitting_std,
            'JMX (%)', 'Chest Strap', 'Einthoven', 'Sitting')

double_plot(cs_jogging_avg, gudb_cs_jogging_std,
            cable_jogging_avg, gudb_cable_jogging_std,
            'JMX (%)', 'Chest Strap', 'Einthoven', 'Jogging')



plt.show()
