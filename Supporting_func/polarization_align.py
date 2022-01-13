import os
import numpy as np
import sys
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from pathlib import Path
from Supporting_func.stocks_coefficients import path_to_data


""" Скрипт принимает результат измерения неполяризованного шумового сигнала в двух поляризациях
    и возвращает коэффициенты, уравнивающие отклик радиометра на такой сигнал в левой и правой 
    круговых поляризациях. Опорной принимается левая поляризация. Коэффициенты записываются в файл 
    'I:\Fast_Acquisition\Alignment\Align_polar.bin'. Таким сигналом может быть излучение черного тела, 
    закрывающего рупор. Однако в данном случае это проходит не на всех частотах."""


def noise_kp(_spectrum_noise, _time_limits):
    """ Функция принимает файл измерения спектра мощности сигнала небо/черное тело
        на входе рупора радиометра и возвращает разность усредненных по времени
        значений спектров мощности сигнала черного тела и неба от частоты. Отсчеты
        от 0 до n0 - небо, от n1 до n2 - черное тело"""
    spectrum_noise = np.array(_spectrum_noise, dtype=np.int64)
    n_row, n_col = np.shape(spectrum_noise)
    s = np.zeros(n_col, dtype=np.int64)
    s1 = np.zeros(n_col)
    n0, n1, n2 = _time_limits
    # Усреднение по времени
    for i in range(n_col):
        s1[i] = np.sum(spectrum_noise[:n0, i]) / n0
        s[i] = np.sum(spectrum_noise[n1:n2, i]) / (n2 - n1)
        # s[i] -= s1[i]

    # fig, ax = plt.subplots(1, figsize=(12, 6))
    # ax.plot(s)
    # ax.plot(s1)
    # plt.show()
    return s


def align_func(_s_loc: object, aver_param_noise: object, _n_nyq: object, _time_limits) -> object:
    """ Функция возвращает спектр, в котором значения спектра в зоне действия
        режекторных фильтров приравнены к 1
    """
    # Исходные данные
    delta_f = 7.8125

    s_diff = noise_kp(_s_loc, _time_limits)

    # Исключение корректировки коэффициента усиления в зоне действия режекторных фильтров
    if _n_nyq == 3:
        n1 = int((90 - delta_f / 2 / aver_param_noise) // (delta_f / aver_param_noise))
        n2 = int((220 - delta_f / 2 / aver_param_noise) // (delta_f / aver_param_noise))
        n3 = int((540 - delta_f / 2 / aver_param_noise) // (delta_f / aver_param_noise))
        n4 = int((700 - delta_f / 2 / aver_param_noise) // (delta_f / aver_param_noise))
        s_diff[n3:n4] = 1
    else:
        n1 = int((80 - delta_f / 2 / aver_param_noise) // (delta_f / aver_param_noise))
        n2 = int((230 - delta_f / 2 / aver_param_noise) // (delta_f / aver_param_noise))
    s_diff[n1:n2] = 1

    # band = int(aver_param_noise / aver_param)
    # kp_band = [np.sum(kp_norm[band * i:band * (i + 1)]) / band for i in range(int(n_col / band))]

    # if _n_nyq == 3:
    #     freq = np.linspace(1000 * (_n_nyq - 1) + 3.9063 / aver_param_noise, 1000 * _n_nyq - 3.9063 / aver_param_noise, n_col)
    # else:
    #     freq = np.linspace(1000 * _n_nyq - 3.9063 / aver_param_noise, 1000 * (_n_nyq - 1) + 3.9063 / aver_param_noise, n_col)

    return s_diff


if __name__ == '__main__':
    # ******************** Путь к исходным данным *********************
    current_data_file = '2021-12-26_07+16'  # Имя файла с исходными текущими данными без расширения
    current_data_dir = '2021_12_26sun'  # Папка с текущими данными
    align_file_name = 'Align_polar.bin'  # Имя файла с текущими коэффициентами выравнивания АЧХ
    current_catalog = r'2021\Results'  # Текущий каталог (за определенный период, здесь - год)

    time_limits = [3600, 5500, 9100]  # От отсчета 0 и до отсчета time_limits[0] рупор открыт - "небо",
    # от отсчета time_limits[1] и до отсчета time_limits[2] рупор закрыт поглощающей пластиной - "черное тело"

    file_path_data, head_path = path_to_data(current_catalog, current_data_dir)
    folder_align_path = Path(head_path, 'Alignment')

    # Если файла хранениия коэффициентов не существует, то создаем его, если существует - загружаем
    columns_names = ['date', 'azimuth', 'att1', 'att2', 'att3',
                     'polar', 'polar_align_low', 'polar_align_high',
                     'flag_align']
    if not os.path.isfile(Path(folder_align_path, align_file_name)):
        calibration_frame = pd.DataFrame(columns=columns_names)
    else:
        with open(Path(folder_align_path, align_file_name), 'rb') as inp:
            calibration_frame = pickle.load(inp)
    # ************************************************************************************

    # Загрузка исходных данных ('_spectrum.npy' - сами данные, _head.bin - вспомогательные)
    spectrum = np.load(Path(file_path_data, current_data_file + '_spectrum.npy'), allow_pickle=True)
    with open(Path(file_path_data, current_data_file + '_head.bin'), 'rb') as inp:
        head = pickle.load(inp)
    n_aver = head['n_aver']
    aver_param = 2 ** (6 - n_aver)
    shift = head['shift']
    spectrum = spectrum * (2 ** shift)
    if head['polar'] == 'left' or head['polar'] == 'right':
        head['polar'] = 'half'

    # Объявление списка для хранения разности спектров небо/черное тело
    s_diff = [np.nan, np.nan, np.nan, np.nan]
    polar_align = [np.nan, np.nan]
    i = 0

    # Расчет разности спектров черное тело/небо по отдельности для каждой поляризации
    # и полосы частот 1-2 или 2-3 ГГц
    flag_align = 0  # Выравнивания по всему диапазону еще нет
    for s in spectrum:
        if i % 2:
            n_nyq = 3
        else:
            n_nyq = 2
        if np.size(s) > 1:
            s_unipol = s  # [s > 100]
            s_diff[i] = align_func(s_unipol, aver_param, n_nyq, time_limits)
        i += 1

    # Расчет коэффициентов для выравнивания поляризации каналов (1-2 ГГц и 2-3 ГГц).
    # Левая поляризация - опорная
    polar_align[0] = s_diff[2] / s_diff[0]
    polar_align[1] = s_diff[3] / s_diff[1]
    for j in [0, 1]:
        i = 0
        for el in polar_align[j]:
            if el > 2 or el < 0.5:
                polar_align[j][i] = 1
            i += 1
    flag_align = 1  # Произошло выравнивание поляризаций по всему диапазону

    # Рассчитанные коэффициенты вместе с исходной информацией записываем в виде словаря  для формирования
    # объекта Series и включения в сводную  таблицу корректирующих коэффициентов
    calibrate_row = {'date': current_data_file[:11], 'azimuth': current_data_file[-3:],
                     'att1': head['att1'], 'att2': head['att2'], 'att3': head['att3'],
                     'polar': head['polar'], 'polar_align_low': polar_align[0], 'polar_align_high': polar_align[1],
                     'flag_align': flag_align}
    # 'flag_align' - признак выравнивания по всему диапазону : 1 - сделано
    calibrate_row_ser = pd.Series(calibrate_row)

    # Определяем, есть ли в сводной таблице данные с такими же исходными параметрами. Если есть, то будет их
    # проверка на то, содержат они все коэффициенты или нет. Если нет, то объект Series будет полностью
    # вставлен в таблицу (объект DataFrame)

    idx = calibration_frame.loc[(calibration_frame.date == head['date'])
                                & (calibration_frame.att1 == head['att1'])
                                & (calibration_frame.att2 == head['att2'])
                                & (calibration_frame.att3 == head['att3'])
                                & (calibration_frame.polar == head['polar'])
                                & (calibration_frame.azimuth == current_data_file[-3:])].index

    if len(idx):
        r = calibration_frame.iloc[idx[0]]
        pass
    else:
        calibration_frame = calibration_frame.append(calibrate_row_ser, ignore_index=True)
        with open(Path(folder_align_path, align_file_name), 'wb') as out:
            pickle.dump(calibration_frame, out)

    with open(Path(folder_align_path, align_file_name), 'rb') as inp:
        calibration_frame_inp = pickle.load(inp)
    fig, ax = plt.subplots(1, figsize=(12, 6))
    ax.plot(polar_align[0])
    ax.plot(polar_align[1])
    plt.show()
    pass
