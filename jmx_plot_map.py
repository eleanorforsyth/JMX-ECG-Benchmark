#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt

# Plots the map which maps jitter in seconds to %

import jmx_analysis

jitter = np.linspace(0,100E-3,100)

jmx = []

for j in jitter:
    jmx.append(jmx_analysis.mapping_jitter(j / jmx_analysis.norm_jitter))

plt.plot(jitter,jmx)
plt.xlabel("jitter / s")
plt.ylabel("score")
plt.show()
