from typing import Any

import numpy as np
import os
import sys
import pickle
import matplotlib.pyplot as plt
import pandas as pd
from Supporting_func.afc_alignment import align_spectrum
from pathlib import Path

current_dir = Path.cwd()
sys.path.insert(0, current_dir)


class CustomError(Exception):
    pass


def initial_scan_cut(data):
    """ Функция принимает данные с левой поляризацией и определяет индексы центров импульсов
    меандра с левой поляризацией. Центры импульсов с правой поляризацией принимаются
    как лежащие посередине между между центрами с левой поляризацией.
    Функция отдает массивы индексов центров с левой и правой поляризациями"""

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
    head_path1 = Path(r'F:\Fast_Acquisition')  # Путь к каталогу данных для домашнего ноута
    head_path2 = Path(r'/media/anatoly/Samsung_T5/Fast_Acquisition')  # Путь к каталогу данных для рабочего компа
    head_path3 = Path(r'C:\SCIENCE\PYTHON 3\Fast_Acquisition')  # Путь к каталогу данных для ноута ВМ
    head_path4 = Path(r'I:\Fast_Acquisition')  # Путь к каталогу данных для notebook 'Khristina'
    if head_path1.is_dir():
        head_path_out = head_path1
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


if __name__ == '__main__':
    align = 'y'

    current_catalog = r'2021\Results'  # Текущий каталог (за определенный период, здесь - год)
    current_data_dir = '2021_03_27sun'  # Папка с текущими данными
    current_data_file = '2021-03-27_06+12'  # Имя файла с исходными текущими данными без расширения
    align_file_name: Any = 'Align_coeff.bin'  # Имя файла с текущими коэффициентами выравнивания АЧХ
    file_path_data, head_path = path_to_data(current_catalog, current_data_dir)


    #               **********************************************
    # ************** Загрузка матрицы спектров и установок (head) *************
    with open(Path(file_path_data, current_data_file + '_head.bin'), 'rb') as inp:
        head = pickle.load(inp)
    spectrum = np.load(Path(file_path_data, current_data_file + '_spectrum.npy'), allow_pickle=True)
    #               **********************************************

    if align == 'y':
        path = Path(head_path, 'Alignment', align_file_name)
        spectrum1 = align_spectrum(spectrum[0], spectrum[1], spectrum[2], spectrum[3], head,
                                   path)

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
    # Параметры Стокса
    s0 = c + d
    s3 = c - d
    mean_frame_ind_pol = (mean_frame_ind_right + mean_frame_ind_left) // 2 * 0.008125
    stocks_coeff = pd.Series([s0, s3, mean_frame_ind_pol])
    np.save(Path(file_path_data, current_data_file + '_stocks'), stocks_coeff)

    for j in range(0, 512, 20):
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        ax1.plot(mean_frame_ind_pol, s0[:, j], label='x(t)')
        ax2.plot(mean_frame_ind_pol, s3[:, j], label='y(t)', color='darkred')
        ax1.set_ylabel('Stocks_I')
        ax2.set_ylabel('Stocks_V', color='darkred')
        ax1.minorticks_on()
        ax1.grid()
        ax1.grid(which='minor',
                 axis='x',
                 color='k',
                 linestyle=':')
        plt.show()

    pass
