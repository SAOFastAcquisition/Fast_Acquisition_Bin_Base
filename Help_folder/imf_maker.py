import numpy as np
import os
import sys
import pandas as pd
import pickle
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.signal import argrelextrema as agl
from Help_folder.paths_via_class import DataPaths


if __name__ == '__main__':
    current_primary_file = '2023-02-10_01+20'
    current_primary_dir = '2023_02_10sun'
    main_dir = '2023'
    adr1 = DataPaths(current_primary_file, current_primary_dir, main_dir)
    converted_data_file_path = adr1.converted_data_file_path
    data_treatment_file_path = adr1.treatment_data_file_path

    path_npy = Path(str(converted_data_file_path) + '_spectrum.npy')
    spectrum = np.load(path_npy, allow_pickle=True)
    data = np.log10(spectrum[4, 4:199])
    l = len(data)
    arg = np.array([i for i in range(l)])
    idx_minimas = agl(data, np.less)[0]
    idx_maximas = agl(data, np.greater)[0]

    cs_min = CubicSpline(idx_minimas, data[idx_minimas])
    cs_max = CubicSpline(idx_maximas, data[idx_maximas])
    xs = np.arange(0, l, 1)
    data1 = data[xs] - (cs_max(xs) + cs_min(xs)) / 2

    fig = plt.figure(figsize=(12, 8))
    axes = fig.add_subplot()
    axes.grid(b=True, which='major', color='#666666', linestyle='-')
    # Show the minor grid lines with very faint and almost transparent grey lines
    axes.minorticks_on()
    axes.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)

    axes.set_title('Antenna Temperature', fontsize=18)
    # axes.set_yscale('log')
    # axes.set_xscale('log')
    axes.set_xlabel('Frequency, MHz', fontsize=18)
    axes.plot(arg, data)
    # axes.plot(idx_minimas, data[idx_minimas])
    # axes.plot(idx_maximas, data[idx_maximas])
    # axes.plot(xs, cs_min(xs))
    # axes.plot(xs, cs_max(xs))
    axes.plot(arg, data1)
    plt.show()
    pass
