import numpy as np
import os
import sys
import pickle
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.signal import argrelextrema as agl
from Help_folder.paths_via_class import DataPaths


def plot_imf(_x, _y, _y1, _y2):
    _fig = plt.figure(figsize=(12, 8))
    _axes = _fig.add_subplot()
    _axes.grid(b=True, which='major', color='#666666', linestyle='-')
    # Show the minor grid lines with very faint and almost transparent grey lines
    _axes.minorticks_on()
    _axes.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)

    _axes.set_title('Antenna Temperature', fontsize=18)
    _axes.set_xlabel('Frequency, MHz', fontsize=18)
    # plt.ylim([0, 1000000])
    _axes.plot(_x, _y)
    _axes.plot(_x, _y1)
    _axes.plot(_x, _y2)
    plt.show()


def imf_gen(_data):
    from scipy.interpolate import make_interp_spline as mis
    _l = len(_data)
    cx = np.array([i for i in range(_l)])
    _idx_minimas = agl(_data, np.less)[0]
    _idx_maximas = agl(_data, np.greater)[0]
    if len(_idx_maximas) <= 1 and len(_idx_minimas) <= 1:
        return print(f"IMF is over: min bin = {_idx_minimas}, max bin = {_idx_maximas}")
    # _cs_min = CubicSpline(_idx_minimas, _data[_idx_minimas])
    l, r = [(2, 0)], [(2, 0)]
    _cs_min = mis(_idx_minimas, _data[_idx_minimas], k=3, bc_type=(l, r))
    # _cs_max = CubicSpline(_idx_maximas, _data[_idx_maximas])
    _cs_max = mis(_idx_maximas, _data[_idx_maximas], k=3, bc_type=(l, r))
    _cs = (_cs_max(cx) + _cs_min(cx)) / 2

    # plot_imf(cx, _data, _cs, _data - _cs)
    return _cs


def imf_proc(_data):
    flag = 0
    sign_change_count = 2
    n_min = 1
    n_max = 1
    data0 = _data.copy()
    cs0 = []
    while flag == 0:
        try:
            cs = imf_gen(data0)
            if len(cs0) == 0:
                cs0 = cs
            else:
                cs0 = np.vstack((cs0, cs))
            data1 = data0 - cs
        except TypeError:
            print('IMFs is over')
            cs0 = np.zeros((1, l))
            break

        sign_change_count1 = null_cross_sec(data1)
        n_min1 = len(agl(data1, np.less)[0])
        n_max1 = len(agl(data1, np.greater)[0])
        if (sign_change_count1 - sign_change_count) == 0 and (n_min1 - n_min) == 0 and (n_max1 - n_max) == 0:
            flag = 1
        else:
            sign_change_count = sign_change_count1
            n_min = n_min1
            n_max = n_max1
            data0 = data1
    _cs01 = np.nansum(cs0, axis=0)
    return _cs01


def imf_decomp(_data):
    _data_cur = _data.copy()
    _imf = []
    while True:
        cs_cur = imf_proc(_data_cur)
        imf_cur = _data_cur - cs_cur
        if not len(_imf):
            _imf = imf_cur
        else:
            _imf = np.vstack((_imf, imf_cur))
        if not np.sum(cs_cur):
            break
        _data_cur = cs_cur
    return _imf


def null_cross_sec(_data):
    from math import copysign
    only_ones_list = [copysign(1, element) for element in _data]

    summa = 0
    for i in range(len(only_ones_list) - 1):
        summa += abs(only_ones_list[i] + only_ones_list[i + 1])

    _sign_change_count = ((len(only_ones_list) - 1) * 2 - summa) / 2
    return int(_sign_change_count)


if __name__ == '__main__':
    current_primary_file = '2023-02-10_01+20'
    current_primary_dir = '2023_02_10sun'
    main_dir = '2023'
    adr1 = DataPaths(current_primary_file, current_primary_dir, main_dir)
    converted_data_file_path = adr1.converted_data_file_path
    data_treatment_file_path = adr1.treatment_data_file_path
    e = 2.718281828459045
    dec = np.log(10)
    path_npy = Path(str(converted_data_file_path) + '_spectrum_time.npy')
    # path_npy = Path(str(converted_data_file_path) + '_scan_freq.npy')
    spectrum = np.load(path_npy, allow_pickle=True)
    spectrum_log = np.log10(spectrum)
    data = spectrum_log[:, 5:197]
    l = np.shape(data)[1]
    arg = np.arange(0, l, 1)

    imf1 = imf_decomp(data[4, :])
    # imf2 = imf_decomp(spectrum_log[4, 5:500])
    # imf3 = imf_decomp(spectrum_log[5, 5:500])

    plot_imf(arg[:], np.exp(dec * (np.sum(imf1, axis=0))),
                            np.exp(dec * imf1[3, :]) * (np.exp(dec * imf1[2, :]) - 1), np.exp(dec * imf1[0, :]))    #

    # plot_imf(arg, imf1[1, :], imf2[1, :], imf3[1, :])
    # plot_imf(arg, imf1[-1, :], imf2[-1, :], imf3[-1, :])

    pass
