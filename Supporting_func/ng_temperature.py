import os
import numpy as np
import sys
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from pathlib import Path
from Supporting_func.stocks_coefficients import path_to_data


def temperature_ngi(_spectrum, polarization, time_boders):
    _delta_t = 8.3886e-3
    _delta_f = 7.8125
    aver_param_noise = 8
    df0 = _delta_f / aver_param_noise
    time0, time1, time2, time3 = time_boders
    n0 = int((time0 + 1) / _delta_t)
    n1 = int((time1 - 1) / _delta_t)
    n2 = int((time1 + 1) / _delta_t)
    n3 = int((time2 - 1) / _delta_t)
    n4 = int((time2 + 1) / _delta_t)
    n5 = int((time3 - 1) / _delta_t)

    if polarization == 'left':
        s0, s1 = _spectrum[0], _spectrum[1]
    elif polarization == 'right':
        s0, s1 = _spectrum[2], _spectrum[3]

    av_gne0 = np.mean(s0[n0:n1, :], axis=0)
    av_gne1 = np.mean(s1[n0:n1, :], axis=0)
    av_ml0 = np.mean(s0[n2:n3, :], axis=0)
    av_ml1 = np.mean(s1[n2:n3, :], axis=0)
    av_gni0 = np.mean(s0[n4:n5, :], axis=0)
    av_gni1 = np.mean(s1[n4:n5, :], axis=0)
    temperature_scale0 = (av_gni0 - av_ml0) / av_gne0
    temperature_scale1 = (av_gni1 - av_ml1) / av_gne1


    k1 = int((90 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))
    k2 = int((220 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))
    k3 = int((540 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))
    k4 = int((700 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))
    temperature_scale1[k1:k2] = 0
    temperature_scale1[k3:k4] = 0

    k5 = 0
    k6 = int((230 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))
    temperature_scale0[k5:k6] = 0
    temperature_scale0[0:] = temperature_scale0[-1::-1]
    f0 = np.array([1000 + df0 / 2 + i * df0 for i in range(1024)])
    f1 = np.array([2000 + df0 / 2 + i * df0 for i in range(1024)])

    temp_nge0 = [temperature_nge(f) for f in f0]
    temp_nge1 = [temperature_nge(f) for f in f1]

    temp_ngi0 = temp_nge0 * temperature_scale0
    temp_ngi1 = temp_nge1 * temperature_scale1
    plt.plot(f0, temp_ngi0)
    plt.plot(f1, temp_ngi1)
    plt.show()
    pass
    return temp_ngi0, temp_ngi1


def temperature_nge(_f, t0=300):
    """ Функция возвращает шумовую температуру внешнего поверочного ГШ в зависимости от
    частоты в МГц, t0 - температура окружающей среды ГШ (по умолчанию 300К) На отрезках
    между известными значениями ИОШТ - линейная аппроксимация"""
    n_1000 = 32.06  # ИОШТ для 1000 МГц
    n_2000 = 31.44  # ИОШТ для 2000 МГц
    n_3000 = 29.89  # ИОШТ для 3000 МГц
    k1 = n_1000 - n_2000
    k2 = n_2000 - n_3000

    if 1000 <= _f <= 2000:
        t_nge = (n_1000 - k1 * (_f / 1000 - 1)) * t0
    elif 2000 < _f <= 3000:
        t_nge = (n_2000 - k2 * (_f / 1000 - 2)) * t0

    return t_nge


if __name__ == '__main__':
    ''' Создаем базу данных по измеренным шумовым температурам от внутреннего ГШ на входе радиометра
    для каждого экземпляра ГШ ***nge_temperature_base_name***, среднее значение шумовой температуры
    от этого ГШ по всем измерениям записывается в файл ***nge_temperature_file_name***'''

    current_data_file = '2022-06-18_07test'  # Имя файла с исходными текущими данными без расширения
    current_primary_dir = '2022_06_18test'
    current_data_dir = current_primary_dir + '_conv'  # Папка с текущими данными
    current_catalog = r'2022/Converted_data'  # Текущий каталог (за определенный период, здесь - год)

    ngi_temperature_base_name = 'ngi_temperature_base.npy'  # Базы данных по температуре ГШ на входе МШУ
    ngi_temperature_file_name = 'ngi_temperature.npy'  # Файл усредненной по базе шумовой температуры для ГШ
    file_path_data, head_path = path_to_data(current_catalog, current_data_dir)
    #                           ***********************
    # Создание или загрузка файла с шумовой температурой встроенного генератора шума
    ngi_temperature_base_path = Path(head_path, 'Alignment', ngi_temperature_base_name)
    ngi_temperature_path = Path(head_path, 'Alignment', ngi_temperature_file_name)
    columns_names_base = ['nge_id', 'date', 'att1', 'att2', 'att3',
                          'polar', 'attempt_num', 'low_band', 'upper_band']
    columns_names = ['nge_id', 'low_band', 'upper_band']
    if not os.path.isfile(ngi_temperature_base_path):
        ngi_temperature_base = pd.DataFrame(columns=columns_names_base)
    else:
        with open(ngi_temperature_base_path, 'rb') as inp:
            ngi_temperature_base = pickle.load(inp)

    if not os.path.isfile(ngi_temperature_path):
        ngi_temperature = pd.DataFrame(columns=columns_names)
    else:
        with open(ngi_temperature_path, 'rb') as inp:
            ngi_temperature = pickle.load(inp)
    #                           ************************

    spectrum = np.load(Path(file_path_data, current_data_file + '_spectrum.npy'), allow_pickle=True)
    with open(Path(file_path_data, current_data_file + '_head.bin'), 'rb') as inp:
        head = pickle.load(inp)
    nge_id = '01'
    attempt_num = current_data_file[11:13]
    idx = ngi_temperature_base.loc[(ngi_temperature_base.nge_id == nge_id)
                                   & (ngi_temperature_base.date == head['date'])
                                   & (ngi_temperature_base.att1 == head['att1'])
                                   & (ngi_temperature_base.att2 == head['att2'])
                                   & (ngi_temperature_base.att3 == head['att3'])
                                   & (ngi_temperature_base.polar == head['polar'])
                                   & (ngi_temperature_base.attempt_num == attempt_num)]

    if not len(idx):
        r = temperature_ngi(spectrum, head['polar'], [0, 20, 30, 50])
        temperature_row = {'nge_id': nge_id, 'date': current_data_file[:11],
                           'att1': head['att1'], 'att2': head['att2'], 'att3': head['att3'],
                           'polar': head['polar'], 'attempt_num': attempt_num,
                           'low_band': r[0], 'upper_band': r[1]
                           }
        ngi_temperature_base = ngi_temperature_base.append(temperature_row, ignore_index=True)

    else:
        pass

    #       *** Запись новых данных или обновление старых ***
    with open(ngi_temperature_base_path, 'wb') as out:
        pickle.dump(ngi_temperature_base, out)
    with open(ngi_temperature_path, 'wb') as out:
        pickle.dump(ngi_temperature, out)
    pass
