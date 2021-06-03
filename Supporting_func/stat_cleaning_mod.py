import numpy as np
import pandas as pd
import sys
import os
import pickle
import matplotlib.pyplot as plt
from Supporting_func.stocks_coefficients import initial_scan_cut
from pathlib import Path
from Supporting_func.stocks_coefficients import path_to_data

current_dir = Path.cwd()
sys.path.insert(0, current_dir)


def clean_func(data, frame_centers_loc, length_frame):
    """ Функция принимает исходные данные data и разбиение индексов данных на отрезки длиной length_frame
     с центрами frame_centers_loc. Для отрезка данных вычисляется стандартное отклонение от среднего и,
     на первом этапе, если отсчет отличается от среднего больше чем на три стандартных отклонения, он
     заменяется  на 0 или 10. В дальнейшем в статистике эти отсчеты не участвуют. Далее, снова вычисляется
     среднее для рассматриваемого отрезка данных и, если для отсчета отношение его отклонения от среднего
     к среднему на отрезке больше трех радиометрических выиграшей, то такой отсчет заменяется на 0 или 10
     и, в конечном итоге, отбрасывается из рассмотрения. """

    tau = 8.1925e-3  # Разрешение по времени (время накопления для одного отсчета)
    radiometric_advantage1 = np.sqrt(2 / tau / freq_resolution)  # радиометрический выигрыш для одного отсчета
    radiometric_advantage2 = np.sqrt(2 / tau / freq_resolution / length_frame)  # радиометрический выигрыш для
    # фрагмента длиной length_frame
    delta_ind = int(length_frame / 2)
    np.delete(frame_centers_loc, -1)

    #                               **************************
    # Сканы для двух частот (используются при отладке) работают вместе с отображением
    # их на рисунках  до статистической обработки и после нее
    if debugging == 'y':
        n = 213
        m = 214
        data01 = np.array(data[:, n], dtype=float)
        data02 = np.array(data[:, m], dtype=float)
        for i in range(np.size(data01)):
            if data01[i] < 100:
                data01[i] = np.nan
            if data02[i] < 100:
                data02[i] = np.nan
    #                               **************************

    for ind_frame in frame_centers_loc:
        frame = data[ind_frame - delta_ind:ind_frame + delta_ind, :]
        form = np.shape(frame)
        for j in range(form[1]):
        # for j in [n, m]:          # Используются при отладке
            frame_small = frame[:, j]
            mean_frame = np.mean(frame_small[frame_small > 100])
            sigma_mean_frame = np.sqrt(np.var(frame_small[frame_small > 100]))
            for i in range(form[0]):
                if np.abs(frame[i, j] - mean_frame) > sigma_mean_frame * 3:
                    frame[i, j] = 10
            frame_small = frame[:, j]
            mean_frame = np.mean(frame_small[frame_small > 100])
            for i in range(form[0]):
                c = frame[i, j]
                d = c / mean_frame - 1
                if np.abs(d) > (3 * radiometric_advantage1):
                    frame[i, j] = 0

        data[ind_frame - delta_ind:ind_frame + delta_ind] = frame
    pass

    #                              **************************
    # Отображения сканов для двух частот (используются при отладке)
    if debugging == 'y':
        fig, ax = plt.subplots(1, figsize=(12, 6))
        data1 = np.array(data[:, n], dtype=float)
        data2 = np.array(data[:, m], dtype=float)
        for i in range(np.size(data1)):
            if data1[i] < 100:
                data1[i] = np.nan
            if data2[i] < 100:
                data2[i] = np.nan
        ax.plot(data1)
        ax.plot(data2)
        ax.plot(data01 * 1.4)
        ax.plot(data02 * 1.4)
        ax.grid()
        plt.show()
    #                              **************************

    return data


def func_frame_unipolar(data, length_frame=63):
    if not length_frame % 2:
        raise ValueError('Число отсчетов должно быть нечетным')
    form_data = np.shape(data)
    n = form_data[0] // length_frame
    frame_centers_local = [length_frame // 2 + length_frame * k for k in range(n)]
    return frame_centers_local


def func_frame_centers():
    """ Если наблюдение велось в двух поляризациях, то функция определяет центры полупериодов с левой
    и правой поляризациями frame_centers_left1, frame_centers_right1. С индексом 1 - для
    полосы 1-2 ГГц, с индексом 2 - 2-3 ГГц.
    Для наблюдения в одной поляризации длина отрезка разбиения (период) по умолчению равна 63. Для незадействованной
    поляризации и полосы частот соответствующие frame_centers являются []"""

    #           **** Для диапазона частот 1-2 ГГц ****
    if np.size(spectrum_input[0]) > 1 and np.size(spectrum_input[2]) > 1:
        frame_centers_left1, frame_centers_right1 = initial_scan_cut(spectrum_input[0])
    elif np.size(spectrum_input[0]) > 1:
        frame_centers_left1 = func_frame_unipolar(spectrum_input[0])
        frame_centers_right1 = []
    elif np.size(spectrum_input[2]) > 1:
        frame_centers_right1 = func_frame_unipolar(spectrum_input[2])
        frame_centers_left1 = []
    else:
        frame_centers_left1 = []
        frame_centers_right1 = []

    #           **** Для диапазона частот 2-3 ГГц ****
    if np.size(spectrum_input[1]) > 1 and np.size(spectrum_input[3]) > 1:
        frame_centers_left2, frame_centers_right2 = initial_scan_cut(spectrum_input[1])
    elif np.size(spectrum_input[1]) > 1:
        frame_centers_left2 = func_frame_unipolar(spectrum_input[1])
        frame_centers_right2 = []
    elif np.size(spectrum_input[3]) > 1:
        frame_centers_right2 = func_frame_unipolar(spectrum_input[3])
        frame_centers_left2 = []
    else:
        frame_centers_left2 = []
        frame_centers_right2 = []

    frame_centers_local = [np.array(frame_centers_left1, dtype=np.int32),
                           np.array(frame_centers_left2, dtype=np.int32),
                           np.array(frame_centers_right1, dtype=np.int32),
                           np.array(frame_centers_right2, dtype=np.int32)]
    return frame_centers_local


if __name__ == '__main__':

    current_data_file = '2021-03-27_06+12'  # Имя файла с исходными текущими данными без расширения
    current_data_dir = '2021_03_27sun'  # Папка с текущими данными
    current_catalog = r'2021\Results'  # Текущий каталог (за определенный период, здесь - год)

    file_path_data, head_path = path_to_data(current_catalog, current_data_dir)

    spectrum_input = np.load(Path(file_path_data, current_data_file + '_spectrum.npy'), allow_pickle=True)
    with open(Path(file_path_data, current_data_file + '_head.bin'), 'rb') as inp:
        head = pickle.load(inp)
    n_aver = int(head['n_aver'])
    freq_resolution = 7.8125e6 / 2 ** (6 - n_aver)
    spectrum_out = [[], [], [], []]
    debugging = 'n'

    # Если запись содержит две поляризации, то длина отрезка разбиения (полупериод переключения
    # поляризаций) равна 34, если поляризация одна то длина отрезка - 63
    if ((np.size(spectrum_input[0]) > 1 and np.size(spectrum_input[2]) > 1) or
            (np.size(spectrum_input[1]) > 1 and np.size(spectrum_input[3]) > 1)):
        length_frame = 34
    else:
        length_frame = 63

    frame_centers = func_frame_centers()

    for i in range(4):
        if np.size(spectrum_input[i]) > 1:
            spectrum = clean_func(spectrum_input[i], frame_centers[i], length_frame)
            spectrum = np.array(spectrum, dtype=np.int32)
            spectrum_out[i] = spectrum
        print('Data ', i, 'calculated')
    print('calc over')
    path_to_data0 = str(Path(file_path_data, current_data_file))
    if not debugging == 'y':
        if not Path(file_path_data, current_data_file + '_spectrum_prime.npy').is_file():
            os.rename(path_to_data0 + '_spectrum.npy', path_to_data0 + '_spectrum_prime.npy')
            np.save(path_to_data0 + '_spectrum', spectrum_out)
        else:
            raise Error('File _spectrum_prime.npy is existing')
