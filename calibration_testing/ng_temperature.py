import os
import numpy as np
import sys
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from pathlib import Path
from Supporting_func.stocks_coefficients import path_to_data


def temperature_ngi(_spectrum, polarization, time_borders):
    _delta_t = 8.3886e-3
    _delta_f = 7.8125
    aver_param_noise = 8
    df0 = _delta_f / aver_param_noise
    time0, time1, time2, time3 = time_borders
    n0 = int((time0 + 1) / _delta_t)  # На входе радиометра сигнал от эталонного ГШ
    n1 = int((time1 - 1) / _delta_t)
    n2 = int((time1 + 1) / _delta_t)  # На входе радиометра сигнал от согласованной нагрузки
    n3 = int((time2 - 1) / _delta_t)
    n4 = int((time2 + 1) / _delta_t)  # На входе радиометра сигнал от внутреннего ГШ
    n5 = int((time3 - 1) / _delta_t)

    if polarization == 'left':
        s0, s1 = _spectrum[0], _spectrum[1]
    elif polarization == 'right':
        s0, s1 = _spectrum[2], _spectrum[3]
    elif polarization == 'both':
        s0, s1 = _spectrum[0], _spectrum[1]

    # Вычисление коэффициентов привязки температуры шумового сигнала от внутреннего
    # ГШ к шумовой температуре эталонного ГШ
    av_gne0 = np.mean(s0[n0:n1, :], axis=0)
    av_gne1 = np.mean(s1[n0:n1, :], axis=0)
    av_ml0 = np.mean(s0[n2:n3, :], axis=0)
    av_ml1 = np.mean(s1[n2:n3, :], axis=0)
    av_gni0 = np.mean(s0[n4:n5, :], axis=0)
    av_gni1 = np.mean(s1[n4:n5, :], axis=0)
    temperature_scale0 = (av_gni0 - av_ml0) / av_gne0
    temperature_scale1 = (av_gni1 - av_ml1) / av_gne1

    # Исключение зон действия режекторных фильтров
    k1 = int((90 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))
    k2 = int((220 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))
    k3 = int((540 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))
    k4 = int((700 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))
    temperature_scale1[k1:k2] = 0
    temperature_scale1[k3:k4] = 0

    k5 = 0
    k6 = int((230 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))
    temperature_scale0[k5:k6] = 0
    temperature_scale0[0:] = temperature_scale0[-1::-1]  # Правильный порядок отсчетов
    # для второй зоны Найквиста

    #                   *************************************
    #   Формирование массивов значений шумовой температуры эталонного источника от частоты
    f0 = np.array([1000 + df0 / 2 + i * df0 for i in range(1024)])
    f1 = np.array([2000 + df0 / 2 + i * df0 for i in range(1024)])
    temp_nge0 = [temperature_nge(f) for f in f0]
    temp_nge1 = [temperature_nge(f) for f in f1]

    temp_ngi0 = temp_nge0 * temperature_scale0
    temp_ngi1 = temp_nge1 * temperature_scale1

    plt.plot(f0, temp_ngi0)
    plt.plot(f1, temp_ngi1)
    plt.show()

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


def average_ngi(_data):
    """ Усредняет по результатам серии измерений"""
    x, y = _data['low_band'], _data['upper_band']
    l = len(x)
    index_x, index_y = x.index, y.index
    array_x, array_y = x[index_x[0]], y[index_y[0]]
    for i in index_x[1:]:
        array_x = np.vstack([array_x, x[i]])
        array_y = np.vstack([array_y, y[i]])
    av_data_x = np.mean(array_x, axis=0)
    # Обнуление недостоверных значений из-за действия фильтров в начале и конце зон Найквиста
    av_data_x[0:20] = 0
    av_data_y = np.mean(array_y, axis=0)
    av_data_y[0:34] = 0
    av_data_y[-34:] = 0

    return av_data_x, av_data_y


def del_random_mod(_s, _s0):
    _l = len(_s)
    for _i in range(1, _l - 2):
        if abs(_s[_i] - _s[_i - 1]) > 2 * abs(_s[_i + 1] - _s[_i - 1]):
            _s[_i] = (_s[_i - 1] + _s[_i + 1]) / 2
        if _s[_i] <= 0.01:
            _s[_i] = _s0
    _s[0] = _s0
    _s[_l - 2:] = _s0
    return _s


def ngi_temperature_update(_ngi_selected, _ngi_id):
    if not os.path.isfile(ngi_temperature_path):
        ngi_temperature = pd.DataFrame(columns=columns_names)
    else:
        with open(ngi_temperature_path, 'rb') as inp:
            ngi_temperature = pickle.load(inp)

    ngi_temperature_low, ngi_temperature_upper = average_ngi(_ngi_selected)
    ngi_average_row = {'ngi_id': ngi_id, 'polar': head['polar'],
                       'low_band': ngi_temperature_low, 'upper_band': ngi_temperature_upper}

    idx_r = ngi_temperature.loc[(ngi_temperature.ngi_id == '01')
                                & (ngi_temperature.polar == head['polar'])].index

    if len(idx_r):
        ngi_temperature = ngi_temperature.drop(idx_r).append(ngi_average_row, ignore_index=True)
    else:
        ngi_temperature = ngi_temperature.append(ngi_average_row, ignore_index=True)
    with open(ngi_temperature_path, 'wb') as _out:
        pickle.dump(ngi_temperature, _out)
    pass


def plot_ngi(_data):
    # Arguments0, s1 = _spectrum[0], _spectrum[1]
    _delta_f = 7.8125
    aver_param_noise = 8
    df0 = _delta_f / aver_param_noise
    f0 = np.array([1000 + df0 / 2 + i * df0 for i in range(1024)])
    f1 = np.array([2000 + df0 / 2 + i * df0 for i in range(1024)])

    x, y = _data['low_band'][_data['polar'] == 'left'], _data['upper_band'][_data['polar'] == 'left']
    index_x, index_y = x.index, y.index
    array_x, array_y = x[index_x[0]], y[index_y[0]]
    f01 = f0
    f11 = f1
    for i in index_x[1:]:
        f01 = np.vstack([f01, f0])
        f11 = np.vstack([f11, f1])
        array_x = np.vstack([array_x, x[i]])
        array_y = np.vstack([array_y, y[i]])

    x1, y1 = _data['low_band'][_data['polar'] == 'right'], _data['upper_band'][_data['polar'] == 'right']
    index_x1, index_y1 = x1.index, y1.index
    array_x1, array_y1 = x1[index_x1[0]], y1[index_y1[0]]
    f02 = f0
    f12 = f1
    for i in index_x1[1:]:
        f02 = np.vstack([f02, f0])
        f12 = np.vstack([f12, f1])
        array_x1 = np.vstack([array_x1, x1[i]])
        array_y1 = np.vstack([array_y1, y1[i]])

    with open(ngi_temperature_path, 'rb') as _in:
        ngi_temperature = pickle.load(_in)

    ngi_left0 = ngi_temperature['low_band'][ngi_temperature['polar'] == 'left'][1]
    ngi_left1 = ngi_temperature['upper_band'][ngi_temperature['polar'] == 'left'][1]
    ngi_right0 = ngi_temperature['low_band'][ngi_temperature['polar'] == 'right'][0]
    ngi_right1 = ngi_temperature['upper_band'][ngi_temperature['polar'] == 'right'][0]

    fig, ax = plt.subplots(1, figsize=(12, 6))
    ax.plot(f0, ngi_left0)
    ax.plot(f1, ngi_left1)
    ax.plot(f0, ngi_right0)
    ax.plot(f1, ngi_right1)
    # ax.scatter(f01, array_x)
    # ax.scatter(f11, array_y)
    # ax.scatter(f02, array_x1)
    # ax.scatter(f12, array_y1)
    # line_minor = plt.plot(inp_scale / 1000)
    # plt.plot(inp_scale, inp_data, 'r--*', inp_scale, 'g:+')  # В кавычках компактно заданы цвет, стиль линии и маркеры
    # plt.setp(fig_main, linestyle=':', color='r')
    # plt.setp(fig, linestyle=':', color=(0, 1, 0, 0.9), marker='.', markerfacecolor='b')  # Четвертое число
    # задает прозрачность линии
    # plt.setp(line_minor, linestyle='-.')
    plt.grid()
    plt.show()

    pass


if __name__ == '__main__':
    ''' Создаем базу данных по измеренным шумовым температурам от внутреннего ГШ на входе радиометра
    для каждого экземпляра ГШ ***nge_temperature_base_name***, среднее значение шумовой температуры
    от этого ГШ по всем измерениям записывается в файл ***nge_temperature_file_name***'''

    current_data_file = '2022-11-18_01'  # Имя файла с исходными текущими данными без расширения
    current_primary_dir = '2022_11_18test'
    current_data_dir = current_primary_dir + '_conv'  # Папка с текущими данными
    current_catalog = r'2022/Test_and_calibration/Converted_data'  # Текущий каталог (за определенный период, здесь - год)

    ngi_temperature_base_name = 'ngi_temperature_base.npy'  # Базы данных по температуре ГШ на входе МШУ
    ngi_temperature_file_name = 'ngi_temperature.npy'  # Файл усредненной по базе шумовой температуры для ГШ
    file_path_data, head_path = path_to_data(current_catalog, current_data_dir)
    #                           ***********************
    # Создание или загрузка файла с шумовой температурой встроенного генератора шума
    ngi_id = '01'  # Идентификатор встроенного генератора шума
    timing = [0, 20, 30, 50]
    ngi_temperature_base_path = Path(head_path, 'Alignment', ngi_temperature_base_name)
    ngi_temperature_path = Path(head_path, 'Alignment', ngi_temperature_file_name)
    columns_names_base = ['nge_id', 'date', 'att1', 'att2', 'att3',
                          'polar', 'attempt_num', 'low_band', 'upper_band']
    columns_names = ['ngi_id', 'polar', 'low_band', 'upper_band']
    if not os.path.isfile(ngi_temperature_base_path):
        ngi_temperature_base = pd.DataFrame(columns=columns_names_base)
    else:
        with open(ngi_temperature_base_path, 'rb') as inp:
            ngi_temperature_base = pickle.load(inp)

    #                           ************************

    spectrum = np.load(Path(file_path_data, current_data_file + '_spectrum.npy'), allow_pickle=True)
    with open(Path(file_path_data, current_data_file + '_head.bin'), 'rb') as inp:
        head = pickle.load(inp)
    r = temperature_ngi(spectrum, head['polar'], timing)
    attempt_num = current_data_file[11:13]
    idx = ngi_temperature_base.loc[(ngi_temperature_base.nge_id == ngi_id)
                                   & (ngi_temperature_base.date == head['date'])
                                   & (ngi_temperature_base.att1 == head['att1'])
                                   & (ngi_temperature_base.att2 == head['att2'])
                                   & (ngi_temperature_base.att3 == head['att3'])
                                   & (ngi_temperature_base.polar == head['polar'])
                                   & (ngi_temperature_base.attempt_num == attempt_num)].index

    if not len(idx):
        r = temperature_ngi(spectrum, head['polar'], timing)
        temperature_row = {'nge_id': ngi_id, 'date': head['date'],
                           'att1': head['att1'], 'att2': head['att2'], 'att3': head['att3'],
                           'polar': head['polar'], 'attempt_num': attempt_num,
                           'low_band': r[0], 'upper_band': r[1]
                           }
        ngi_temperature_base = ngi_temperature_base.append(temperature_row, ignore_index=True)
        with open(ngi_temperature_base_path, 'wb') as out:
            pickle.dump(ngi_temperature_base, out)
    else:
        print('Данные этой записи уже внесены в базу')

    #               *** Усреднение по базе данных ***
    # шумовой температуры встроенного ГШ на входе Радиометра
    # раздельно для каналов левой и правой поляризаций.
    pw = input('Обновить усредненную шумовую температуру ГШ на входе РМ?\t')
    idx_drop = ngi_temperature_base.loc[ngi_temperature_base.att3 == 15].index
    if len(idx_drop):
        for n in idx_drop:
            ngi_temperature = ngi_temperature_base.drop([n])
    if pw == 'y':
        ngi_selected1 = ngi_temperature_base[ngi_temperature_base['nge_id'].isin([ngi_id])]
        head['polar'] = 'right'
        ngi_selected1l = ngi_selected1[ngi_selected1['polar'].isin([head['polar']])]
        ngi_temperature_update(ngi_selected1l, ngi_id)
        head['polar'] = 'left'
        ngi_selected1l = ngi_selected1[ngi_selected1['polar'].isin([head['polar']])]
        ngi_temperature_update(ngi_selected1l, ngi_id)

    #               *** Picture ***
    ngi_selected1 = ngi_temperature_base[ngi_temperature_base['nge_id'].isin([ngi_id])]
    plot_ngi(ngi_selected1)
    # plt.plot(ngi_selected1l['low_band'][1])
    # plt.plot(ngi_selected1l['upper_band'][1])
    # plt.plot(ngi_temperature_low_r)
    # plt.plot(ngi_temperature_upper_r)
    # plt.show()
    # *******************************************
    #       *** Запись новых данных или обновление старых ***

    with open(ngi_temperature_base_path, 'rb') as inp:
        ngi_temperature_base_wr = pickle.load(inp)
    pass
