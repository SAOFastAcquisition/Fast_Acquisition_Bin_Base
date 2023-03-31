import numpy as np
import os
import sys
import pickle
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.signal import argrelextrema as agl
from Help_folder.paths_via_class import DataPaths


def imf_gen(_data):
    from scipy.interpolate import make_interp_spline as mis
    _l = len(_data)
    cx = np.array([i for i in range(_l)])
    _idx_minimas = agl(_data, np.less)[0]
    _idx_maximas = agl(_data, np.greater)[0]

    # _cs_min = CubicSpline(_idx_minimas, _data[_idx_minimas])
    l, r = [(1, 0.0)], [(1, 0.0)]
    _cs_min = mis(_idx_minimas, _data[_idx_minimas], k=3, bc_type=(l, r))
    # _cs_max = CubicSpline(_idx_maximas, _data[_idx_maximas])
    _cs_max = mis(_idx_maximas, _data[_idx_maximas], k=3, bc_type=(l, r))
    _cs = (_cs_max(cx) + _cs_min(cx)) / 2

    _fig = plt.figure(figsize=(12, 8))
    _axes = _fig.add_subplot()
    _axes.grid(b=True, which='major', color='#666666', linestyle='-')
    # Show the minor grid lines with very faint and almost transparent grey lines
    _axes.minorticks_on()
    _axes.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)

    _axes.set_title('Antenna Temperature', fontsize=18)
    _axes.set_xlabel('Frequency, MHz', fontsize=18)
    _axes.plot(cx, _data)
    # axes.plot(idx_minimas, data[idx_minimas])
    # axes.plot(idx_maximas, data[idx_maximas])
    _axes.plot(cx, _cs)
    # axes.plot(xs, cs_max(xs))
    _axes.plot(cx, _data-_cs)
    plt.show()
    return _cs


def null_cross_sec(_data):
    from math import copysign
    only_ones_list = [copysign(1, element) for element in _data]

    summ = 0
    for i in range(len(only_ones_list) - 1):
        summ += abs(only_ones_list[i] + only_ones_list[i + 1])

    _sign_change_count = ((len(only_ones_list) - 1) * 2 - summ) / 2
    return _sign_change_count


if __name__ == '__main__':
    current_primary_file = '2023-02-10_01+20'
    current_primary_dir = '2023_02_10sun'
    main_dir = '2023'
    adr1 = DataPaths(current_primary_file, current_primary_dir, main_dir)
    converted_data_file_path = adr1.converted_data_file_path
    data_treatment_file_path = adr1.treatment_data_file_path

    path_npy = Path(str(converted_data_file_path) + '_spectrum.npy')
    spectrum = np.load(path_npy, allow_pickle=True)
    data = np.log10(spectrum[4, 5:197])
    l = len(data)
    arg = np.arange(0, l, 1)
    cs = imf_gen(data)
    data1 = data[arg] - cs
    sign_change_count1 = null_cross_sec(data1)
    n_min1 = len(agl(data1, np.less)[0])
    n_max1 = len(agl(data1, np.greater)[0])

    cs1 = imf_gen(data1)
    data2 = data1[arg] - cs1
    sign_change_count2 = null_cross_sec(data2)
    n_min2 = len(agl(data2, np.less)[0])
    n_max2 = len(agl(data2, np.greater)[0])

    cs2 = imf_gen(data2)
    data3 = data2[arg] - cs2
    sign_change_count3 = null_cross_sec(data3)
    n_min3 = len(agl(data3, np.less)[0])
    n_max3 = len(agl(data3, np.greater)[0])

    cs3 = imf_gen(data3)
    data4 = data3[arg] - cs3
    sign_change_count4 = null_cross_sec(data4)
    n_min4 = len(agl(data4, np.less)[0])
    n_max4 = len(agl(data4, np.greater)[0])

    cs4 = imf_gen(data4)
    data5 = data4[arg] - cs4
    sign_change_count5 = null_cross_sec(data5)
    n_min5 = len(agl(data5, np.less)[0])
    n_max5 = len(agl(data5, np.greater)[0])

    cs5 = imf_gen(data5)
    data6 = data5[arg] - cs5
    sign_change_count6 = null_cross_sec(data6)
    n_min6 = len(agl(data5, np.less)[0])
    n_max6 = len(agl(data5, np.greater)[0])

    cs0 = cs + cs1 + cs2 + cs3 + cs4 + cs5

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
    axes.plot(arg, cs0)
    # axes.plot(xs, cs_max(xs))
    axes.plot(arg, data - cs0)
    plt.show()
    pass
