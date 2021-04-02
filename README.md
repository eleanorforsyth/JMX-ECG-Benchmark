# JMX Benchmark Analysis for heartbeat detectors

A benchmarking analysis method that generates an overall benchmark for ECG detector algorithms, based on measurements of temporal jitter (J), missed beats (M) and extra detections (X), as these are the three ways errors show up independent of the source of the error. The Glasgow University GUDB ECG recordings database ([Porr and Luis, 2018](http://dx.doi.org/10.5525/gla.researchdata.716)) is used for testing as it has annotated R-peaks for reference. The detector algorithms are tested using recordings for all subjects, all exercises, and Einthoven II and chest strap leads (Einthoven I and Einthoven III can additionally be used if desired). Initially, the code should be run using a wide range of detector algorithms (seven in this case) to generate global mean values for the three error types, and these are subsequently used as the reference for the JMX Benchmark values generated for individual detector algorithms. Default values can be used to compare a single detector algorithm to the results obtained for the seven detector algorithms contained in ‘py-ecg-detectors’ ([Porr and Luis, 2019](https://doi.org/10.5281/zenodo.3353396)).

The benchmark gives a score between 0-100 where 100 is defined as the ideal detector. The ideal detector has no extra beats, no missed beats, and a mean absolute deviation (MAD) of zero for temporal jitter when using our interval analysis method. The benchmark is independent of application-specific performance to truly represent an overall score that can be compared across applications encompassing all the ways errors show up in ECG heartbeat detection.

## Installation

Install py- ecg-detectors ([Porr and Luis, 2019](https://doi.org/10.5281/zenodo.3353396)).

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

Install ecg_gudb_database ([Porr and Luis, 2018](https://pypi.org/project/ecg-gudb-database/)) via the package manager [pip](https://pip.pypa.io/en/stable/) or pip3 which is a version of the pip installer for Python3.

```bash
pip install ecg_gudb_database
pip3 install ecg_gudb_database
```


## Usage

### jmx_detector_analysis.py

The Python code jmx_detector_analysis.py uses the seven heartbeat detectors specified in py-ecg-detectors ([Porr and Luis, 2019](https://doi.org/10.5281/zenodo.3353396)). When run, it applies the detectors to the GUDB ECG database ([Porr and Luis, 2018](http://dx.doi.org/10.5525/gla.researchdata.716)). Initially, please create a folder named 'saved_csv' in the current directory to save the output CSV files. Additionally, it can be noted that the matched filter detector can use a default template or user-generated averaged PQRST shape for each subject.

The code runs the detectors with all subjects, all leads, and all experiments. The heartbeat detection locations are saved and the scenarios that have annotated R-peaks are used for benchmarking the detectors. The extra and missed beats are identified and counted for each of the detectors. For the temporal jitter, an interval analysis is used. For this, the extra beats and missed beats are excluded as they have already been accounted for. After discounting the extra/missed beats, the remaining detected beats could be considered "true" detections - but do not necessarily fall on a consistent location of the PQRST shape. The interval between the detected heartbeat points is found, and then the difference in samples is taken between that and the corresponding annotated interval. For each ECG recording with annotations passed through the detectors, the extra beats, the missed beats and the temporal jitter are all saved as separate CSV files. This allows for all the 'raw' interval analysis data to be available for subsequent optional analysis. 


The positional marker data for extra and missed beats can be used to plot markers onto source recordings to indicate where and how a detector is behaving with individual ECG recordings.

The code generates a 'global average' for the count of extra beats, missed beats, and the mean absolute deviation (MAD) of the temporal jitter from the interval analysis, which is saved as a CSV file called global_results.csv to be used as reference values for the overall benchmarking method. In the new overall benchmarking method the average of these three components from all seven detectors is later mapped to an equivalent benchmark score of 50. 

The breakdown scores of each of the individual detectors before being mapped against the global average are saved in det_lead_stats.csv. 

***Benchmarking your personal detector***

To benchmark a detector of your choice that is not one of the seven available from py-ecg-detectors, firstly, create a folder named 'saved_csv' in the working directory to store the generated CSV files. Add your detector algorithm code in the form of a function to the Class ‘Detectors’ in ecgdetectors.py. Within jmx_detector_analysis.py the new detector can either be added to the list ‘all_detectors’ to compare with others, or, replace all others in the list for individual benchmarking.

When testing a new detector ‘save_global_results’ should be set to ‘False’, as only a representative sample of trusted detectors should be used to generate the global mean values against which detectors are benchmarked, and adding an unknown or single detector would not generate a valid global value.

### jmx_stats_analysis.py

This code is run after jmx_detector_analysis.py. In the case that you chose to run that code with the seven available detectors from the GUDB database then it displays results from the detector data CSV files saved.

The pre-mapped scores are taken from the generated ‘det_lead_stats.csv’ file. Further values are taken from the ‘data_det_lead.csv’ file, which contains complete arrays of data in the same categories as ‘det_lead_stats.csv’.
The results for each detectors extra beats, missed beats and temporal jitter are normalised and then mapped using a curve to allow all three to be multiplied together. The mapping curve gives 1.0 for ideal performance, 0.5 for average performance, and 0 if the normalised performance value is greater than ten times the global mean value. The overall JMX Benchmark values for Einthoven II and chest strap results are shown together on a bar graph for comparison. The breakdown benchmark scores of results by category (missed, extra, temporal jitter) are plotted separately.

When jmx_stats_analysis.py is used to analyse the results of your personal detector, the same list of ‘all_detectors’ should be used. If no ‘global_results.csv’ file exists, saved default values generated from the GUDB and seven detectors in ecgdetectors.py will be used as references for the benchmark calculation.
