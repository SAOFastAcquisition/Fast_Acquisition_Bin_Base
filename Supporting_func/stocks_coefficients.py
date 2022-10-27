from typing import Any

import numpy as np
import os
import sys
import pickle
import matplotlib.pyplot as plt
import pandas as pd
from Supporting_func.afc_alignment import align_spectrum
from pathlib import Path
from Supporting_func.Fig_plot import fig_multi_axes
from Supporting_func.dict_calibr_from_csv import start_stop_calibr, calibration_temp

current_dir = Path.cwd()
sys.path.insert(0, current_dir)


class CustomError(Exception):
    pass


def initial_scan_cut(data):
    """ Функция принимает данные с левой поляризацией и определяет индексы центров импульсов
    меандра с левой поляризацией. Центры импульсов с правой поляризацией принимаются
    как лежащие посередине между между центрами с левой поляризацией.
    Функция возвращает массивы индексов центров с левой и правой поляризациями"""

    data_shape = np.shape(data)
    j = 0
    time_row = np.array(data[:, j])
    time_ind = np.argwhere(time_row > 100)

    # Выбор скана (выбор хорошо заполненного отсчетами скана) для определения
    # положения центра импульса меандра с левой поляризацией в отсчетах скана
    while np.size(time_ind) < 0.3 * data_shape[0]:
        time_row = np.array(data[:, j])
        j += 1
        # Выбор индексов ненулевых значений скана
        time_ind = np.argwhere(time_row > 100)

    # Отсечение начального участка скана до первого контролируемого переключения на
    # левую поляризацию
    i = 0
    while time_ind[i + 1] - time_ind[i] < 25:
        i += 1
    # Задание первого элемента для определения скачка индекса ненулевых членов
    i_prev = int(time_ind[i + 1])
    time_ind_cut = time_ind[time_ind >= i_prev]

    # ******************************************************************************
    # Определение центров импульсов меандра с левой поляризацией в отсчетах скана **
    mean_frame_ind_left = np.array([])
    frame_ind0 = np.array([])
    for i in time_ind_cut:
        # Фрагмент скана с последовательно меняющимися индексами (допускается скачок индекса не более чем на 25)
        if i - i_prev < 25:
            frame_ind0 = np.append(frame_ind0, i)
        else:
            # Средний индекс отсчетов за полупериод усреднения
            try:
                mean_ind = int((np.max(frame_ind0) + np.min(frame_ind0)) / 2)
            except ValueError:
                pass

            # if mean_ind - mean_frame_ind_left[-1] < 50 or mean_ind - mean_frame_ind_left[-1] > 70:
            #     if mean_frame_ind_left[-1] + 60 < data_shape1[0]:
            #         mean_ind = mean_frame_ind_left[-1] + 60
            #     else:
            #             mean_ind = data_shape1[0] - 2

            mean_frame_ind_left = np.append(mean_frame_ind_left, mean_ind)
            frame_ind0 = np.array([])
        i_prev = i
    len_left = np.size(mean_frame_ind_left)
    mean_frame_ind_right = np.zeros(len_left - 1)

    for k in range(len_left - 1):
        mean_frame_ind_right[k] = int((mean_frame_ind_left[k] + mean_frame_ind_left[k + 1]) / 2)
    # Отбрасываем последний элемент вектора центров импульсов левой поляризации, чтобы выровнять
    # по размеру с вектором центров импульсов правой поляризации
    mean_frame_ind_left0 = mean_frame_ind_left[: -1]
    return mean_frame_ind_left0, mean_frame_ind_right


def pol_intensity(data, mean_time_ind):
    """ Функция принимает исходную матрицу спектров левой или правой поляризаций и вектор индексов середин
    полупериодов соответствующих поляризации. Значения усредняются по полупериоду. В усреднении принимают
    принимают участие значения больше 100. Функция отдает матрицу спектров, усредненных по полупериоду
    переключения поляризаций. Размерность матрицы по первому индексу уменьшается примерно в 60
    раз для периода переключения поляризаций 0.5 сек."""

    data_shape = data.shape
    scan_len = np.size(mean_time_ind)
    pol_spectrum = np.ones((scan_len, data_shape[1])) * 10
    k2 = 0
    for j in range(data_shape[1]):
        time_row = np.array(data[:, j])
        k1 = 0
        for i in mean_time_ind:
            frame = time_row[int(i) - 15:int(i) + 14]
            pol_spectrum[k1, k2] = np.mean(frame[frame > 100])
            k1 += 1
        k2 += 1
        # plt.grid()
        # plt.plot(mean_time_ind, pol_spectrum[:, k2 - 1])
        # plt.show()
        pass
    return pol_spectrum


def path_to_data(current_catalog_in, current_data_dir_in):
    """ Функция принимает текущий каталог данных (за год или период) и папку текущих данных (за выбранный день).
    Определяет путь к папке текущих данных на конкретной машине и к корню каталога. """
    head_path1 = Path(r'H:\Fast_Acquisition')  # Путь к каталогу данных для домашнего ноута
    head_path1a = Path(r'E:\Fast_Acquisition')  # Путь к каталогу данных для домашнего ноута
    head_path2 = Path(r'/media/anatoly/Samsung_T5/Fast_Acquisition')  # Путь к каталогу данных для рабочего компа
    head_path3 = Path(r'D:\Fast_acquisition')  # Путь к каталогу данных для ноута ВМ
    head_path4 = Path(r'J:\Fast_Acquisition')  # Путь к каталогу данных для notebook 'Khristina'

    if head_path1.is_dir():
        head_path_out = head_path1
    elif head_path1a.is_dir():
        head_path_out = head_path1a
    elif head_path2.is_dir():
        head_path_out = head_path2
    elif head_path3.is_dir():
        head_path_out = head_path3
    elif head_path4.is_dir():
        head_path_out = head_path4
    else:
        raise CustomError('Path to data is not found!')

    file_path_data_out = Path(head_path_out, current_catalog_in, current_data_dir_in)
    return file_path_data_out, head_path_out


def freq_mask(_i):
    _n1 = 1
    _n2 = 6
    _freq_mask = [
        [1350],                                                               # [0]
        [2060, 2300, 2500, 2750, 2830, 2920],               # [1]
        [1020, 1260, 1340, 1430, 1540, 1670, 1750, 1930],                           # [2]
        [1000 * _n1 + 101 * _n2 + 20 * i for i in range(7)],                 # [3]
        [1050, 1465, 1535, 1600, 1700, 2265, 2550, 2700, 2800, 2920],         # [4]
        [1230, 1560, 2300, 2910],                                                               # [5]
        [1140, 1420, 1480, 2460, 2500, 2780],   # for Crab '2021-06-28_03+14' # [6]
        [1220, 1540, 1980, 2060, 2500, 2780],   # for Crab '2021-06-28_04+12' # [7]
        [1171, 1380, 1465, 1600, 1700, 2265, 2530, 2720, 2800, 2920]    # [8]
    ]
    return _freq_mask[_i]


def two_fig_plot():

    fig, axes = plt.subplots(2, 1, figsize=(12, 12))
    _i = 0
    max_y2 = np.nanmax(s0[:, num_mask])
    _i_max = len(num_mask)

    for _j in num_mask:
        f1 = freq_mask0[_i]
        text1 = 'f = ' + str(f1) + ' MHz'
        axes[0].plot(mean_frame_ind_pol, s0[:, _j], label=text1)
        axes[1].plot(mean_frame_ind_pol, s3[:, _j])
        axes[0].set_ylabel('Stocks_I')
        axes[1].set_ylabel('Stocks_V', color='darkred')
        axes[0].minorticks_on()
        axes[1].minorticks_on()
        axes[0].legend(loc='upper right')
        _i += 1
    axes[0].grid()
    axes[1].grid()
    axes[0].grid(which='minor',
                 axis='x',
                 color='k',
                 linestyle=':')
    plt.show()
    pass


def twin_fig_plot():

    for j in num_mask:
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()  # Создание второй оси ординат (справа). Ось абсцисс общая
        ax1.plot(mean_frame_ind_pol, s0[:, j], label='x(t)')
        ax2.plot(mean_frame_ind_pol, s3[:, j], label='y(t)', color='darkred')
        # ax1.plot([i for i in range(m)], s0[:, j], label='x(t)')
        # ax2.plot([i for i in range(m)], s3[:, j], label='y(t)', color='darkred')
        ax1.set_ylabel('Stocks_I')
        ax2.set_ylabel('Stocks_V', color='darkred')
        ax1.minorticks_on()
        f1 = int(1000 + 1000 / (n + 1) + j * 2000 / 1025)
        max_y2 = np.nanmax(s3[:, j])
        text1 = 'f = ' + str(f1) + ' MHz'
        plt.text(0, max_y2 / 2, text1, fontsize=12)  # Разрешение по времени
        ax1.grid()
        ax2.grid()
        ax1.grid(which='minor',
                 axis='x',
                 color='k',
                 linestyle=':')
        plt.show()


if __name__ == '__main__':
    align = 'n'
    channel_align = 'n'

    current_catalog = r'2022/Converted_data'  # Текущий каталог (за определенный период, здесь - год)
    current_primary_dir = '2022_06_18sun'
    current_data_dir = current_primary_dir + '_conv'  # Папка с текущими данными
    current_data_file = '2022-06-18_01+28'  # Имя файла с исходными текущими данными без расширения
    align_file_name: Any = 'Align_coeff.bin'  # Имя файла с текущими коэффициентами выравнивания АЧХ
    file_path_data, head_path = path_to_data(current_catalog, current_data_dir)
    path_to_stocks = Path(file_path_data, current_data_file + '_stocks.npy')
    path_to_stocks_left_txt = Path(file_path_data, current_data_file + '_left.txt')
    path_to_stocks_right_txt = Path(file_path_data, current_data_file + '_right.txt')
    freq_mask_list = freq_mask(8)
    freq_mask0 = np.array(freq_mask(8))

    if not (os.path.isfile(path_to_stocks)):

        #               **********************************************
        # ************** Загрузка матрицы спектров и установок (head) *************
        with open(Path(file_path_data, current_data_file + '_head.bin'), 'rb') as inp:
            head = pickle.load(inp)
        spectrum = np.load(Path(file_path_data, current_data_file + '_spectrum.npy'), allow_pickle=True)
        #               **********************************************
        channel_align = 'n'
        if align == 'y':
            path = Path(head_path, 'Alignment', align_file_name)
            spectrum1 = align_spectrum(spectrum[0], spectrum[1], spectrum[2], spectrum[3], head,
                                       path, 0)

        input_data_upper = {'left': spectrum[1],
                            'right': spectrum[3]}
        input_data_lower = {'left': spectrum[0],
                            'right': spectrum[2]}

        a1 = input_data_upper['left']
        b1 = input_data_upper['right']

        a = input_data_lower['left']
        # a = a[-1::-1][:]
        b = input_data_lower['right']
        # b = b[-1::-1][:]
        mean_frame_ind_left, mean_frame_ind_right = initial_scan_cut(a)

        c = pol_intensity(a, mean_frame_ind_left)
        d = pol_intensity(b, mean_frame_ind_right)
        c1 = pol_intensity(a1, mean_frame_ind_left)
        d1 = pol_intensity(b1, mean_frame_ind_right)
        c = np.hstack((c, c1))
        d = np.hstack((d, d1))
        #                               ****************
        # Вычисление выравнивающий коэффициентов по калибровочному сигналу - калибровочный сигнал д.б. неполяризованным
        if channel_align == 'y':
            dict_calibr_file_name = 'dict_calibr.csv'  # Имя файла с текущими коэффициентами выравнивания АЧХ
            path_to_csv = Path(file_path_data, dict_calibr_file_name)
            s = start_stop_calibr(current_data_file, path_to_csv)
            av_c_cal = (np.mean(c[s[0]:s[1], :], axis=0) + np.mean(c[s[2]:s[3], :], axis=0)) / 2
            av_d_cal = (np.mean(d[s[0]:s[1], :], axis=0) + np.mean(d[s[2]:s[3], :], axis=0)) / 2
            noise_coeff = av_c_cal / av_d_cal
            m, n = np.shape(c)
            for i in range(m):
                for j in range(n):
                    d[i, j] = d[i, j] * noise_coeff[j]
            #                               *****************

        np.savetxt(path_to_stocks_left_txt, c)
        np.savetxt(path_to_stocks_right_txt, d)
        # Параметры Стокса
        s0 = c + d
        s3 = c - d
        mean_frame_ind_pol = (mean_frame_ind_right + mean_frame_ind_left) // 2 * 0.008336
        stocks_coeff = pd.Series([s0, s3, mean_frame_ind_pol])
        np.save(path_to_stocks, stocks_coeff)

    else:
        stocks_coeff = np.load(path_to_stocks, allow_pickle=True)
        [s0, s3, mean_frame_ind_pol] = stocks_coeff

    m, n = np.shape(s0)
    freq_res = n   # Число отсчетов спектра шириной 2 ГГц по частоте
    num_mask = [int((s - 1000 * (1 + 1 / freq_res)) * freq_res / 2000) for s in freq_mask0]
    # calibration_temperature = [calibration_temp(f) for f in freq_mask0]
    # dict_calibr_file_name = 'dict_calibr.csv'  # Имя файла с текущими коэффициентами выравнивания АЧХ
    # path_to_csv = Path(file_path_data, dict_calibr_file_name)
    # s = start_stop_calibr(current_data_file, path_to_csv)
    # c = (s3 + s0) / 2   # левая поляризация
    #
    # for j in range(np.size(freq_mask0)):
    #     temp_mass = c[s[0]:s[1], num_mask[j]]
    #     av_c_cal = (np.mean(c[s[0]:s[1], num_mask[j]]) + np.mean(c[s[2]:s[3], num_mask[j]]))
    #     temp_coeff = calibration_temperature[j] / av_c_cal
    #     s0[:, num_mask[j]] = s0[:, num_mask[j]] * temp_coeff
    #     s3[:, num_mask[j]] = s3[:, num_mask[j]] * temp_coeff

    two_fig_plot()
    # twin_fig_plot()

    with open(Path(file_path_data, current_data_file + '_head.bin'), 'rb') as inp:
        head = pickle.load(inp)
    inform = [('time resol = ' + str(60 * 8.3886e-3) + 'sec'),
                ('freq resol = ' + str(int(2000 // (n + 1))) + 'MHz'),
                ('polarisation ' + head['polar']), 'align: ' + 'yes', 'kurtosis quality = ' + str(head['good_bound'])]
    current_treatment_dir = current_primary_dir + '_treat'
    current_treatment_path = Path('Data_treatment', current_treatment_dir)
    s0_selected = s3[:, num_mask]
    # fig_multi_axes(np.transpose(s0_selected), mean_frame_ind_pol, inform,  Path(current_treatment_path,
    #                                                                             current_data_file), freq_mask0, head)
    pass
