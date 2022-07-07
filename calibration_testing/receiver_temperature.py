import os
import numpy as np
import sys
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from pathlib import Path
from Supporting_func.stocks_coefficients import path_to_data

if __name__ == '__main__':
    current_data_dir = '2022'
    primary_data_dir = 'Primary_data'  # Каталог исходных данных (за определенный период, здесь - год)
    converted_data_dir = 'Converted_data'  # Каталог для записи результатов конвертации данных и заголовков
    data_treatment_dir = 'Data_treatment'  # Каталог для записи результатов обработки, рисунков

    current_primary_dir = '2022_06_27test'
    current_converted_dir = current_primary_dir + '_conv'
    current_converted_path = Path(converted_data_dir, current_converted_dir)
    current_treatment_dir = current_primary_dir + '_treat'
    current_treatment_path = Path(data_treatment_dir, current_treatment_dir)

    ngi_temperature_file_name = 'ngi_temperature.npy'  # Файл усредненной по базе шумовой температуры для ГШ
    current_primary_file1 = '2022-06-27_01'  # Файл с согласованной нагрузкой на входе приемника
    current_primary_file2 = '2022-06-27_02'  # Файл с КЗ на входе приемника

    converted_data_file_path, head_path = path_to_data(current_data_dir, current_converted_path)
    data_treatment_file_path, head_path = path_to_data(current_data_dir, current_treatment_path)
    ngi_temperature_path = Path(head_path, 'Alignment', ngi_temperature_file_name)
    folder_align_path = Path(head_path, 'Alignment')
    date = current_primary_file1[0:10]

    delta_t = 8.3886e-3
    _delta_f = 7.8125
    aver_param_noise = 8
    t0 = 300
    time_mask = [11, 19, 11, 19]


    path1 = Path(converted_data_file_path, current_primary_file1 + '_spectrum.npy')
    spectrum1 = np.load(path1, allow_pickle=True)
    path2 = Path(converted_data_file_path, current_primary_file2 + '_spectrum.npy')
    spectrum2 = np.load(path2, allow_pickle=True)

    t0_matched, t1_matched = int(time_mask[0] / delta_t), int(time_mask[1] / delta_t)
    t0_short, t1_short = int(time_mask[2] / delta_t), int(time_mask[3] / delta_t)

    sp_matched0, sp_matched1, sp_matched2, sp_matched3 = spectrum1[0][t0_matched:t1_matched, :], \
                                                         spectrum1[1][t0_matched:t1_matched, :], \
                                                         spectrum1[2][t0_matched:t1_matched, :], \
                                                         spectrum1[3][t0_matched:t1_matched, :]
    sp_short0, sp_short1, sp_short2, sp_short3 = spectrum2[0][t0_short:t1_short, :], \
                                                 spectrum2[1][t0_short:t1_short, :], \
                                                 spectrum2[2][t0_short:t1_short, :], \
                                                 spectrum2[3][t0_short:t1_short, :]

    spm0_av, spm1_av, spm2_av,  spm3_av = np.mean(sp_matched0, axis=0), np.mean(sp_matched1, axis=0), \
                                          np.mean(sp_matched2, axis=0), np.mean(sp_matched3, axis=0)
    sps0_av, sps1_av, sps2_av, sps3_av = np.mean(sp_short0, axis=0), np.mean(sp_short1, axis=0), \
                                         np.mean(sp_short2, axis=0), np.mean(sp_short3, axis=0)

    # Исключение зон действия режекторных фильтров



    # sps1_av[k1:k2] = 2
    # sps1_av[k3:k4] = 2


    # sps0_av[k5:k6] = 2
    # spm0_av[k5:k6] = 2
    spm0_av[0:] = spm0_av[-1::-1]  # Правильный порядок отсчетов
    sps0_av[0:] = sps0_av[-1::-1]  # Правильный порядок отсчетов
    # для второй зоны Найквиста

    spm_left, spm_right = np.hstack((spm0_av, spm1_av)), np.hstack((spm0_av, spm1_av))
    sps_left, sps_right = np.hstack((sps0_av, sps1_av)), np.hstack((sps0_av, sps1_av))

    t_receiver = sps_left / (2 * spm_left - sps_left) * t0
    k1 = int((1090 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))
    k2 = int((1220 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))
    k3 = int((1520 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))
    k4 = int((1700 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))
    k6 = 1024
    k5 = int((770 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))
    t_receiver[0:20] = 100
    t_receiver[1024:1058] = 100
    t_receiver[-34:] = 100
    t_receiver[k1:k2] = 100
    t_receiver[k3:k4] = 100
    t_receiver[k5:k6] = 100

    fig, ax = plt.subplots(1, figsize=(12, 6))
    ax.plot(t_receiver)
    # ax.plot(spm_right)
    plt.grid()
    plt.show()

    pass
