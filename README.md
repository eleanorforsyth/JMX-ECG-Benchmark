# JMX Benchmark Analysis for heartbeat detectors

A benchmarking analysis method that generates an overall benchmark for
ECG detector algorithms, based on measurements of temporal jitter (J),
missed beats (M) and extra detections (X), as these are the three ways
errors show up independent of the source of the error. The Glasgow
University GUDB ECG recordings database ([Howell and Porr,
2018](http://dx.doi.org/10.5525/gla.researchdata.716)) is used for
testing as it has annotated R-peaks for reference. The detector
algorithms are tested using recordings for all subjects, all
exercises, and Einthoven II and chest strap leads (Einthoven I and
Einthoven III can additionally be used if desired).

The benchmark gives a score between 0-100 where 100 is defined as the
ideal detector. The ideal detector has no extra beats, no missed
beats, and a mean absolute deviation (MAD) of zero for temporal jitter
when using our interval analysis method. The benchmark is independent
of application-specific performance to truly represent an overall
score that can be compared across applications encompassing all the
ways errors show up in ECG heartbeat detection.

## Installation

Install py- ecg-detectors ([Howell and Porr, 2019](https://doi.org/10.5281/zenodo.3353396)).

Linux / Mac:
```bash
pip3 install py-ecg-detectors [--user]
```
Windows:
```bash
pip install py-ecg-detectors [--user]
```
From source:
```bash
python3 setup.py install [--user]
```
*Use the option --user if you don't have system-wise write permission.*

---
Install ecg_gudb_database ([Howell and Porr, 2018](https://pypi.org/project/ecg-gudb-database/)) via the package manager [pip](https://pip.pypa.io/en/stable/) or pip3 which is a version of the pip installer for Python3.

```bash
pip install ecg_gudb_database
pip3 install ecg_gudb_database
```

---
Download the templates file and the folder is placed within the working directory. As well as individualised templates it includes the default templates for 250Hz and 360Hz.

## Usage

### jmx_evaluate_all_detectors.py

The code runs the detectors with all subjects, all leads, and all
experiments. The heartbeat detection locations are saved and the
scenarios that have annotated R-peaks are used for benchmarking the
detectors. After discounting the
extra/missed beats, the remaining detected beats could be considered
"true" detections

### jmx_stats_plots.py

The overall JMX Benchmark values for Einthoven
II and chest strap results are shown together on a bar graph for
comparison for sitting and jogging.
