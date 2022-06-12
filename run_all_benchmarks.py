#!/usr/bin/python3
"""
This script benchmarks all detectors with both the MIT database and the GU database.
For the GU database it runs it for both the chest strap and Einthoven II with loose cables.
You need to download both the MIT arrhythmia database from: https://alpha.physionet.org/content/mitdb/1.0.0/
and the GU database from: http://researchdata.gla.ac.uk/716/
Both need to be placed below this directory: "../mit-bih-arrhythmia-database-1.0.0/" for the MITDB and
and "../dataset_716" for the GU database.
"""
import numpy as np
import pathlib
from tester_GUDB import GUDB_test
from ecgdetectors import Detectors
from multiprocessing import Process

def run_tests(leads):
    # GUDB database testing
    gu_test = GUDB_test()
    gu_detectors = Detectors(250)
    gu_test.classifer_test_all(tolerance=0, config=leads)

pgustrap = Process(target=run_tests, args=('chest_strap',))
pgustrap.start()
pgucables = Process(target=run_tests, args=('loose_cables',))
pgucables.start()
    
pgustrap.join()
pgucables.join()
