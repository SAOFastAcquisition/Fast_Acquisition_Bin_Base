import os
import numpy as np
import sys
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from pathlib import Path
from Supporting_func.stokes_coefficients import path_to_data


def del_random_mod(_s, _s0):
    _l = len(_s)
    _s = np.array(_s)
    _s[np.isnan(_s)] = _s0
    for _i in range(1, _l - 2):
        if abs(_s[_i] - _s[_i - 1]) > 2 * abs(_s[_i + 1] - _s[_i - 1]):
            _s[_i] = (_s[_i - 1] + _s[_i + 1]) / 2
        if _s[_i] <= 0.01:
            _s[_i] = _s0
    _s[0] = _s0
    _s[_l - 2:] = _s0
    return _s


def func_path(data_dir, year='2022'):
    """ Функция принимает имя каталога с исходными данными. Возвращает пути к директориям,
    в которые будут записываться результаты обработки исходных данных"""
    _current_dir = year
    _primary_data_dir = 'Primary_data'  # Каталог исходных данных (за определенный период, здесь - год)
    _converted_data_dir = 'Converted_data'  # Каталог для записи результатов конвертации данных и заголовков
    _data_treatment_dir = 'Data_treatment'  # Каталог для записи результатов обработки, рисунков

    _current_primary_dir = data_dir
    _current_primary_path = Path(_primary_data_dir, _current_primary_dir)
    _current_converted_dir = _current_primary_dir + '_conv'
    _current_converted_path = Path(_converted_data_dir, _current_converted_dir)
    _current_treatment_dir = _current_primary_dir + '_treat'
    _current_treatment_path = Path(_data_treatment_dir, _current_treatment_dir)

    _primary_data_dir_path, _head_path = path_to_data(_current_primary_dir, _current_primary_path)
    _converted_data_dir_path, _head_path = path_to_data(_current_dir, _current_converted_path)
    _data_treatment_dir_path, _head_path = path_to_data(_current_dir, _current_treatment_path)
    pass
    return _primary_data_dir_path, _converted_data_dir_path, _data_treatment_dir_path, _head_path


def zone_deletion(_s, _flag):
    # Исключение зон действия режекторных фильтров
    k1 = int((1090 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))
    k2 = int((1220 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))
    k3 = int((1520 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))
    k4 = int((1700 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))
    k6 = 1024
    k5 = int((770 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))

    if _flag == 'matched':
        _s[0:20] = 100
        _s[1024:1058] = 100
        _s[-46:] = 100
        _s[k1:k2] = 100
        _s[k3:k4] = 100
        _s[k5:k6] = 100
    else:
        _s[0:20] = 50
        _s[1024:1058] = 50
        _s[-38:] = 50
        _s[k1:k2] = 50
        _s[k3:k4] = 50
        _s[k5:k6] = 50

    return _s


def horn_loss(_date):
    """ Работаем с файлами, содержащими КЗ на раскрыве рупора и согласованную нагрузку
        на входе приемника. Такая нагрузка эквивалентна черному телу в раскрыве рупора при
        равенстве физических температур черного тела и входного устройства, в котором
        локализованы омические потери антенны"""
    # _date = '2022-06-28'
    #                   *** Пути к папкам ***
    _primary_data_dir_path, _converted_data_dir_path, _data_treatment_dir_path, _head_path = \
        func_path(cal_dir)
    #                   *** Пути к файлам ***
    path_horn_short = Path(_converted_data_dir_path, horn_short_file + '_spectrum.npy')
    path_horn_matched = Path(_converted_data_dir_path, matched_load_file + '_spectrum.npy')
    _receiver_temperature_file_name = 'receiver_temperature.npy'
    _receiver_temperature_path = Path(_head_path, 'Alignment', _receiver_temperature_file_name)

    spectrum_short_p = np.load(path_horn_short, allow_pickle=True)
    spectrum_matched_p = np.load(path_horn_matched, allow_pickle=True)

    #       *** Поляризации сигналов КЗ на рупоре и согласованной нагрузки на входе ***
    with open(Path(_converted_data_dir_path, horn_short_file + '_head.bin'), 'rb') as inp:
        _head = pickle.load(inp)
    polar = _head['polar']
    with open(Path(_converted_data_dir_path, matched_load_file + '_head.bin'), 'rb') as inp:
        _head_matched = pickle.load(inp)
    polar_matched = _head_matched['polar']

    #   Берем из базы данных шумовую температуру приемника
    with open(_receiver_temperature_path, 'rb') as _inp:
        _data = pickle.load(_inp)
    try:
        receiver_temperature = _data['temperature'][_data.polar == _head['polar']][_data._date == _date].iloc[0]
    except IndexError:
        receiver_temperature = _data['temperature'][_data.polar == 'left'][_data._date == _date].iloc[0]
        receiver_temperature_right = _data['temperature'][_data.polar == 'right'][_data._date == _date].iloc[0]

    #       Если сигнал с одной поляризацией (левой или правой)
    if polar == 'left':
        spectrum_short_p = spectrum_short_p[0:2]
        spectrum_matched_p = spectrum_matched_p[0:2]
    elif polar == 'right':
        spectrum_short_p = spectrum_short_p[2:]
        spectrum_matched_p = spectrum_matched_p[2:]
    else:
        # Учесть, что происходит, если поляризации сигналов при КЗ и при согласованной нагрузке не совпадают
        # При КЗ сигнал в двух поляризациях polar == 'both', а согласованной нагрузке - в одной
        # polar_matched == 'left' or polar_matched == 'right'.
        if polar_matched == 'left':
            spectrum_short_p = spectrum_short_p[0:2]
            spectrum_matched_p = spectrum_matched_p[0:2]
        elif polar_matched == 'right':
            spectrum_short_p = spectrum_short_p[2:]
            spectrum_matched_p = spectrum_matched_p[2:]
        pass

        pass

    num_mask_short = [int(s / delta_t) for s in time_mask_short]
    num_mask_matched = [int(s / delta_t) for s in time_mask_matched]

    spectrum_short = np.array([s[num_mask_short[0]:num_mask_short[1]] for s in spectrum_short_p])
    spectrum_matched = np.array([s[num_mask_matched[0]:num_mask_matched[1]] for s in spectrum_matched_p])

    s_av_short = [[np.mean(s[s > 100]) for s in s1.transpose()] for s1 in spectrum_short]
    s_av_short = np.array([del_random_mod(s, 100) for s in s_av_short])
    s_av_matched = [[np.mean(s[s > 100]) for s in s1.transpose()] for s1 in spectrum_matched]
    s_av_matched = np.array([del_random_mod(s, 100) for s in s_av_matched])
    s_av_short[0] = s_av_short[0][-1::-1]
    s_av_matched[0] = s_av_matched[0][-1::-1]
    s_av_sh = np.hstack((s_av_short[0], s_av_short[1]))
    s_av_sh = zone_deletion(s_av_sh, 'short')
    s_av_mh = np.hstack((s_av_matched[0], s_av_matched[1]))
    s_av_mh = zone_deletion(s_av_mh, 'matched')

    #   Активные потери во входном устройстве (рупор + "гибрид")
    loss = np.sqrt((300 - receiver_temperature) / (300 + receiver_temperature) * s_av_mh / (s_av_mh - s_av_sh))

    if polar == 'both' and polar_matched == 'both':
        s_av_short[2] = s_av_short[2][-1::-1]
        s_av_matched[2] = s_av_matched[2][-1::-1]
        s_av_sh_r = np.hstack((s_av_short[2], s_av_short[3]))
        s_av_sh_r = zone_deletion(s_av_sh_r, 'short')
        s_av_mh_r = np.hstack((s_av_matched[2], s_av_matched[3]))
        s_av_mh_r = zone_deletion(s_av_mh_r, 'matched')
        #   Активные потери во входном устройстве (рупор + "гибрид")
        loss_r = np.sqrt((300 - receiver_temperature_right) / (300 + receiver_temperature_right) *
                         s_av_mh_r / (s_av_mh_r - s_av_sh_r))

    _fig, _ax = plt.subplots(1, figsize=(12, 6))
    # _ax.plot(s_av_sh)
    # _ax.plot(s_av_mh)
    _ax.plot(loss)
    # _ax.plot(loss_r)
    plt.grid()
    plt.show()
    pass


if __name__ == '__main__':

    """ Расчет выравнивающих коэффициентов АЧХ приемника по шумовому сигналу от согласованной нагрузки на входе
    с учетом собственных шумов приемника. Расчет активных потерь "рупор + гибрид" по измерениям реакции радиометра 
    на КЗ на апертуре рупора и согласованную нагрузку на входе приемника"""

    delta_t = 8.3886e-3
    _delta_f = 7.8125
    aver_param_noise = 8

    current_data_dir = '2022'
    current_primary_dir = '2022_06_28test'
    #       **************************************************************************
    #                *** КЗ на рупоре и согласованная нагрузка на входе ***
    #                  *** Считаем активные потери в "Рупор + гибрид" ***
    horn_short_file = '2022-06-22_03calibr'
    matched_load_file = '2022-06-22_08test'
    cal_dir = '2022_06_22calibr'
    time_mask_short = [21, 39]
    time_mask_matched = [21, 29]

    ngi_temperature_file_name = 'ngi_temperature.npy'  # Файл усредненной по базе шумовой температуры для ГШ
    receiver_temperature_file_name = 'receiver_temperature.npy'
    date = '2022-06-27'     # Дата измерения исходных данных для вычисления шумовой температуры приемника
                            # (2022-06-27 или 2022-06-28)

    horn_loss(date)
    #       ***************************************************************************

    current_primary_file1 = '2022-06-28_01test'  # Файл с согласованной нагрузкой на обоих входах приемника
    current_primary_file2 = '2022-06-28_02test'  # Файл с согласованной нагрузкой и КЗ на входах приемника
    current_primary_file2 = '2022-06-28_03test'  # Файл с КЗ и согласованной нагрузкой на входах приемника

    primary_data_dir_path, converted_data_dir_path, data_treatment_dir_path, head_path = \
        func_path(current_primary_dir)
    ngi_temperature_path = Path(head_path, 'Alignment', ngi_temperature_file_name)
    receiver_temperature_path = Path(head_path, 'Alignment', receiver_temperature_file_name)

    date = current_primary_file1[0:10]

    aver_param_noise = 8
    temp0 = 300
    time_mask = [11, 19, 21, 29]

    # Загружаем результаты измерений с черным телом на входе рупора
    dir_horn_measure = '2022_06_27sun'
    file_horn_measure = '2022-06-27_00ant-04'
    primary_data_dir_path, converted_data_dir_path, data_treatment_dir_path, head_path = \
        func_path(dir_horn_measure)
    path_horn = Path(converted_data_dir_path, file_horn_measure + '_spectrum.npy')

    time_mask_horn = [32.5, 47]
    num_mask_horn = [int(s / delta_t) for s in time_mask_horn]

    spectrum2 = np.load(path_horn, allow_pickle=True)
    with open(Path(converted_data_dir_path, file_horn_measure + '_head.bin'), 'rb') as inp:
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
