import numpy as np
import os
import sys
import pickle
import matplotlib.pyplot as plt
import pandas as pd
from Supporting_func.afc_alignment import align_spectrum

sys.path.insert(0, r'E:/rep1/Supporting_func')


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
    mean_frame_ind_right = np.zeros(len_left-1)

    for k in range(len_left-1):
        mean_frame_ind_right[k] = int((mean_frame_ind_left[k] + mean_frame_ind_left[k + 1])/2)
    # np.delete(mean_frame_ind_left, mean_frame_ind_left[-1], 0)
    mean_frame_ind_left0 = mean_frame_ind_left[: -1]
    return mean_frame_ind_left0, mean_frame_ind_right


def pol_intensity1(data):
    data_shape = data.shape
    pol_spectrum = np.array([])
    for j in range(data_shape[1]):
        time_row = np.array(data[:, j])
        # Выбор индексов ненулевых значений скана
        time_ind = np.argwhere(time_row > 100)
        if np.size(time_ind) > 0.35 * data_shape[0]:
            # j += 1
            # time_row = np.array(data[:, j])
            # time_ind = np.argwhere(time_row > 100)
            # Перевод массива в одномерный
            time_ind = np.ravel(time_ind)
            # Задание первого элемента для определения скачка индекса ненулевых членов
            i_prev = time_ind[0]
            # Фрагмент скана с последовательно меняющимися индексами (допускается скачок индекса не более чем на 5)
            frame_ind = np.array([], 'int32')
            frame_ind0 = np.array([], 'int32')

            mean_frame_ind = np.array([], 'int32')
            frame_mean = np.array([])

            # Обработка скана при фиксированной частоте
            for i in time_ind:
                # Фрагмент скана с последовательно меняющимися индексами (допускается скачок индекса не более чем на 20)
                if i - i_prev < 25:
                    frame_ind0 = np.append(frame_ind0, i)
                else:
                    # Средний индекс отсчетов за полупериод усреднения
                    if np.size(frame_ind0) == 0:
                        mean_ind = mean_frame_ind[-1] + 60
                        frame_ind0 = mean_ind
                    else:
                        mean_ind = int((np.max(frame_ind0) + np.min(frame_ind0)) / 2)

                    if i > 160 and ((mean_ind - mean_frame_ind[-1]) < 50 or (mean_ind - mean_frame_ind[-1]) > 70):
                        if mean_frame_ind[-1] + 60 < data_shape[0]:
                            mean_ind = mean_frame_ind[-1] + 60
                        else:
                            mean_ind = data_shape[0] - 2
                        frame_ind01 = time_ind[mean_ind - 15 < time_ind]
                        frame_ind0 = frame_ind01[frame_ind01 < mean_ind + 15]
                    try:
                        b = time_row[frame_ind0]
                    except IndexError:
                        print('stop')
                        pass
                    # Усреднение за полупериод интенсивности поляризации (левой или правой)
                    mean = np.mean(b[b > 100])
                    if np.size(b[b > 100]) == 0:
                            mean = 10

                    # try:
                    #     b[np.abs(time_row[frame_ind] - mean) > 3 * var] = 0
                    # except RuntimeWarning:
                    #     pass

                    frame_mean = np.append(frame_mean, mean)
                    mean_frame_ind = np.append(mean_frame_ind, mean_ind)
                    frame_ind = np.append(frame_ind, frame_ind0)
                    frame_ind0 = np.array([], 'int32')
                pass
                i_prev = i
            frame_mean.shape = -1, 1
            if j == 0:
                pol_spectrum = np.array(frame_mean)
            else:
                try:
                    pol_spectrum = np.append(pol_spectrum, frame_mean, 1)
                except ValueError:
                    pass

            pass
            # plt.grid()
            # plt.plot(mean_frame_ind, frame_mean)
            # plt.show()

            frame_mean = np.array([])
            # Заполнение вырезанных куртозисом частотных составляющих значением 10, чтобы сохранить
            # размер матрицы скан-спектр
        else:
            frame_mean = np.ones(np.size(pol_spectrum[:, j - 1])) * 10
            frame_mean.shape = -1, 1
            pol_spectrum = np.append(pol_spectrum, frame_mean, 1)
            pass
    return pol_spectrum


def pol_intensity(data, mean_time_ind):
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


if __name__ == '__main__':
    # head_path = 'E:\\Measure_res'
    align = 'y'
    head_path = r'F:\Fast_Acquisition\2020'
    file_path = r'F:\Fast_Acquisition\2021\Results' + r'\2021_03_27sun\2021-03-27_06+12'
    folder_align_path = r'F:\Fast_Acquisition\Alignment'
    align_file_name = r'\Align_coeff.bin'
    with open(file_path + '_head.bin', 'rb') as inp:
        head = pickle.load(inp)
    spectrum = np.load(file_path + '_spectrum.npy', allow_pickle=True)
    if align == 'y':
        path_output = folder_align_path + align_file_name
        # spectr_extr_left1, spectr_extr_left2, spectr_extr_right1, spectr_extr_right2 = \
        spectrum1 = align_spectrum(spectrum[0], spectrum[1], spectrum[2], spectrum[3], head,
                           path_output)
    # file_name0 = head_path + '\\Measure\\Fast_Acquisition\\2020_12_09test\\20201209_2chan_ML_-00_3'
    # file_name0_left = head_path + r'\2020_12_22sun\20201222_pol2_06_+04_3_left'
    # file_name0_right = head_path + r'\2020_12_22sun\20201222_pol2_06_+04_3_right'

    input_data_upper = {'left': spectrum[1],
                  'right': spectrum[3]}
    input_data_lower = {'left': spectrum[0],
                        'right': spectrum[2]}

    # a = input_data_upper['left']
    # b = input_data_upper['right']
    a = input_data_lower['left']
    # a = a[-1::-1][:]
    b = input_data_lower['right']
    # b = b[-1::-1][:]
    mean_frame_ind_left, mean_frame_ind_right = initial_scan_cut(a)
    c = pol_intensity(a, mean_frame_ind_left)
    d = pol_intensity(b, mean_frame_ind_right)

    # Параметры Стокса
    s0 = c + d
    s1 = c - d
    mean_frame_ind_pol = (mean_frame_ind_right + mean_frame_ind_left) // 2
    stocks_coeff = pd.Series([s0, s1, mean_frame_ind_pol])
    np.save(file_path + '_stocks', stocks_coeff)

    for j in range(0, 256, 10):
        plt.grid()
        plt.plot(mean_frame_ind_pol, s0[:, j])
        plt.plot(mean_frame_ind_pol, s1[:, j])
        plt.show()

    pass
