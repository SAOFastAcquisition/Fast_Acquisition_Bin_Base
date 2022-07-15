import os
import numpy as np
import sys
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from pathlib import Path
from Supporting_func.stocks_coefficients import path_to_data


def del_random_mod(_s, _s0):
    _l = len(_s)
    _s = np.array(_s)
    _s[np.isnan(_s)] = 100
    for _i in range(1, _l - 2):
        if abs(_s[_i] - _s[_i - 1]) > 2 * abs(_s[_i + 1] - _s[_i - 1]):
            _s[_i] = (_s[_i - 1] + _s[_i + 1]) / 2
        if _s[_i] <= 0.01:
            _s[_i] = 10
    _s[0] = _s0
    _s[_l - 2:] = _s0
    return _s


if __name__ == '__main__':

    """ Расчет выравнивающих коэффициентов АЧХ приемника по шумовому сигналу от согласованной нагрузки на входе
    с учетом собственных шумов приемника."""

    current_data_dir = '2022'
    primary_data_dir = 'Primary_data'  # Каталог исходных данных (за определенный период, здесь - год)
    converted_data_dir = 'Converted_data'  # Каталог для записи результатов конвертации данных и заголовков
    data_treatment_dir = 'Data_treatment'  # Каталог для записи результатов обработки, рисунков

    current_primary_dir = '2022_06_28test'
    current_converted_dir = current_primary_dir + '_conv'
    current_converted_path = Path(converted_data_dir, current_converted_dir)
    current_treatment_dir = current_primary_dir + '_treat'
    current_treatment_path = Path(data_treatment_dir, current_treatment_dir)

    ngi_temperature_file_name = 'ngi_temperature.npy'  # Файл усредненной по базе шумовой температуры для ГШ
    receiver_temperature_file_name = 'receiver_temperature.npy'
    current_primary_file1 = '2022-06-28_01test'  # Файл с согласованной нагрузкой на обоих входах приемника
    current_primary_file2 = '2022-06-28_02test'  # Файл с согласованной нагрузкой и КЗ на входах приемника
    current_primary_file2 = '2022-06-28_03test'  # Файл с КЗ и согласованной нагрузкой на входах приемника

    converted_data_file_path, head_path = path_to_data(current_data_dir, current_converted_path)
    data_treatment_file_path, head_path = path_to_data(current_data_dir, current_treatment_path)
    ngi_temperature_path = Path(head_path, 'Alignment', ngi_temperature_file_name)
    receiver_temperature_path = Path(head_path, 'Alignment', receiver_temperature_file_name)

    date = current_primary_file1[0:10]

    delta_t = 8.3886e-3
    _delta_f = 7.8125
    aver_param_noise = 8
    temp0 = 300
    time_mask = [11, 19, 21, 29]

    # Загружаем результаты измерений при черном теле на входе рупора
    dir_horn_measure = '2022_06_27sun'
    file_horn_measure = '2022-06-27_00ant-04'
    horn_converted_dir = dir_horn_measure + '_conv'
    horn_converted_path = Path(converted_data_dir, horn_converted_dir)
    dir_horn_treatment = dir_horn_measure + '_treat'
    dir_horn_treat_path = Path(data_treatment_dir, dir_horn_treatment)
    file_horn_measure_path, head_path_horn = path_to_data(current_data_dir, horn_converted_path)
    path_horn = Path(file_horn_measure_path, file_horn_measure + '_spectrum.npy')

    time_mask_horn = [32.5, 47]
    num_mask_horn = [int(s / delta_t) for s in time_mask_horn]

    spectrum2 = np.load(path_horn, allow_pickle=True)
    with open(Path(file_horn_measure_path, file_horn_measure + '_head.bin'), 'rb') as inp:
        head = pickle.load(inp)
    att1, att2, att3 = head['att1'], head['att2'], head['att3']

    spectrum_BB = np.array([s[num_mask_horn[0]:num_mask_horn[1]] for s in spectrum2])
    a = spectrum_BB[0]
    a1 = a[:, 221]
    b = a1[a1 > 100]

    s_av_BB = [[np.mean(s[s > 100]) for s in s1.transpose()] for s1 in spectrum_BB]
    s_av_BB = [del_random_mod(s, 100) for s in s_av_BB]

    fig, ax = plt.subplots(1, figsize=(12, 6))
    ax.plot(s_av_BB[0])
    ax.plot(s_av_BB[1])
    ax.plot(s_av_BB[2])
    ax.plot(s_av_BB[3])
    plt.grid()
    plt.show()
    pass