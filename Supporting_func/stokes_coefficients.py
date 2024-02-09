from typing import Any
import numpy as np
import os
import gzip
import sys
import pickle
import matplotlib.pyplot as plt
import pandas as pd
from tkinter import *
from tkinter import messagebox as mb
from Supporting_func.afc_alignment import align_spectrum
from pathlib import Path
from scipy.fftpack import fft, ifft
from scipy.signal import lfilter, filtfilt
from Supporting_func.Fig_plot import fig_multi_axes, graph_contour_2d
from Supporting_func.dict_calibr_from_csv import start_stop_calibr, calibration_temp
from Help_folder.paths_via_class import DataPaths
from Supporting_func.sun_az_spead import sun_az_speed

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
    while np.size(time_ind) < 0.47 * data_shape[0]:
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
    data_shape1 = np.shape(time_ind_cut)
    # ******************************************************************************
    # Определение центров импульсов меандра с левой поляризацией в отсчетах скана **
    _mean_frame_ind_left = np.array([])
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

            if len(_mean_frame_ind_left) > 1:
                if mean_ind - _mean_frame_ind_left[-1] < 50 or mean_ind - _mean_frame_ind_left[-1] > 70:
                    if _mean_frame_ind_left[-1] + 60 < time_ind_cut[-1]:
                        mean_ind = _mean_frame_ind_left[-1] + 60
                    else:
                        mean_ind = time_ind_cut[-2]

            _mean_frame_ind_left = np.append(_mean_frame_ind_left, mean_ind)
            frame_ind0 = np.array([])
        i_prev = i
    len_left = np.size(_mean_frame_ind_left)
    _mean_frame_ind_right = np.zeros(len_left - 1)

    for k in range(len_left - 1):
        _mean_frame_ind_right[k] = int((_mean_frame_ind_left[k] + _mean_frame_ind_left[k + 1]) / 2)
    # Отбрасываем последний элемент вектора центров импульсов левой поляризации, чтобы выровнять
    # по размеру с вектором центров импульсов правой поляризации
    _mean_frame_ind_left0 = _mean_frame_ind_left[: -1]
    return _mean_frame_ind_left0, _mean_frame_ind_right


def pol_intensity(data, mean_time_ind):
    """ Функция принимает исходную матрицу спектров левой или правой поляризаций и вектор индексов середин
    полупериодов соответствующих поляризации. Значения усредняются по полупериоду. В усреднении принимают
    принимают участие значения больше 20. Функция отдает матрицу спектров, усредненных по полупериоду
    переключения поляризаций. Размерность матрицы по первому индексу уменьшается примерно в 30
    раз для периода переключения поляризаций 0.5 сек. В полупериоде, в котором значения данной поляризации
    не регистрировались, среднее значение интенсивности поляризации линейно аппроксимируется по двум соседним
    значениям этой поляризации"""

    data_shape = data.shape
    scan_len = np.size(mean_time_ind)
    pol_spectrum = np.ones((scan_len * 2 - 1, data_shape[1])) * 10
    k2 = 0
    for j in range(data_shape[1]):
        time_row = np.array(data[:, j])
        k1 = 0
        for i in mean_time_ind:
            frame = time_row[int(i) - 15:int(i) + 14]
            pol_spectrum[2 * k1, k2] = np.nanmean(frame[frame > 10])
            k1 += 1
        for k in range(scan_len - 1):
            pol_spectrum[2 * k + 1, k2] = (pol_spectrum[2 * k, k2] + pol_spectrum[2 * k + 2, k2]) / 2
        k2 += 1

    return pol_spectrum


def stokes_v_deviation(_s, _m):
    """
    Функция Вычисляет отклонение параметра Стокса V от тенденции. Тенденция получается в результате
    фильтрации с постоянной времени порядка 10 сек. Отдает отклонение
    :param _s:
    :param _m:
    :return:
    """
    shape = np.shape(_s)
    _s1 = np.ones(shape)
    for _i in range(shape[1]):
        _s1[:, _i] = maf_fir(_s[:, _i], _m)
        _s1[:, _i] = _s[:, _i] - _s1[:, _i]
    return _s1


def frame_ind_extend():
    """
    Функция объединяет последовательности отсчетов, соответствующих центрам полупериодов переключения левой и
    правой поляризаций в одну последовательность
    :return: int
    последовательность с чередующимися номерами отсчетов, соответствующих центрам полупериодов с левой
    и правой поляризациями
    """
    _mean_frame_ind = [0] * np.shape(c)[0]
    _mean_frame_ind[::2] = mean_frame_ind_right[:-1]
    _mean_frame_ind[1::2] = mean_frame_ind_left[1:]
    return _mean_frame_ind


def path_to_data(current_catalog_in, current_data_dir_in):
    """
    Функция принимает текущий каталог данных (за год или период) и папку текущих данных (за выбранный день).
    Определяет путь к папке текущих данных на конкретной машине и к корню каталога.
    """
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


def maf_fir(_s, _m=2):
    """
        Calculate moving average filter as FIR

        Parameters
        ----------
        _m : int
            moving average step
        _s : int
            data
        """

    return filtfilt(np.ones(_m - 1), 1, _s) / (_m - 1) / (_m - 1)


def freq_mask(_i):
    _n1 = 1
    _n2 = 5
    _freq_mask = [
        [1736],  # [0]
        [2060, 2300, 2500, 2750, 2830, 2920],  # [1]
        [1020, 1100, 1200, 1300, 1350, 1400, 1450, 1600],  # [2]
        [1000 * _n1 + 100 * _n2 + 90 + 5 * i for i in range(10)],  # [3]
        [1050, 1465, 1535, 1600, 1700, 2265, 2550, 2700, 2800, 2920],  # [4]
        [1230, 1560, 2300, 2910],  # [5]
        [1140, 1420, 1480, 2460, 2500, 2780],  # for Crab '2021-06-28_03+14' # [6]
        [1220, 1540, 1980, 2060, 2500, 2780],  # for Crab '2021-06-28_04+12' # [7]
        [1200, 1300, 1465, 1600, 1700, 2265, 2510, 2720, 2800, 2920]  # [8]
    ]
    return _freq_mask[_i]


def two_fig_plot(_path_to_fig_folder):
    _pic_name = pic_name(_path_to_fig_folder, 0, 'svg')
    _path_to_pic = Path(_path_to_fig_folder, _pic_name)
    fig, axes = plt.subplots(2, 1, figsize=(12, 12))
    _fig_folder = str(_path_to_fig_folder)
    title1, title2, title3 = title_func(_fig_folder, head)

    y_max = np.nanmax(s0[:, num_mask])
    y_min = np.nanmin(s0[:, num_mask])
    _i_max = len(num_mask)
    _i = 0
    for _j in num_mask:
        f1 = freq_mask0[_i]
        text1 = 'f = ' + str(f1) + ' MHz'
        axes[0].plot(mean_frame_ind_pol, s0[:, _j], label=text1)
        axes[1].plot(mean_frame_ind_pol, s3[:, _j])
        axes[0].legend(loc='upper right')
        _i += 1

    axes[0].set_title('Stokes Parameters ' + title1, fontsize=20)
    axes[0].set_ylabel('Stokes_I')
    axes[1].set_ylabel('Stokes_V', color='darkred')
    axes[0].minorticks_on()
    axes[1].minorticks_on()

    y1 = y_max - 2 * (y_max - y_min) / 10
    y2 = y_max - 3 * (y_max - y_min) / 10
    axes[0].text(0, y1, inform[0], fontsize=12)  # Разрешение по частоте
    axes[0].text(0, y2, inform[1], fontsize=12)  # Разрешение по времени

    axes[0].grid()
    axes[1].grid()
    axes[0].grid(which='minor',
                 axis='x',
                 color='k',
                 linestyle=':')
    axes[1].grid(which='minor',
                 axis='x',
                 color='k',
                 linestyle=':')
    plt.show()
    #                               ********************************
    #                        ************ Сохранение рисунка ****************
    fig.savefig(_path_to_pic)
    flag_save = save_question()
    if flag_save == 'no':
        if os.path.isfile(_path_to_pic):
            os.remove(_path_to_pic)
            print('Picture is not saved')
        else:
            print('File not found')
    else:
        print('Picture is saved')
    pass


def some_fig_plot(_path_to_fig_folder, _s_i, _s_v, _s_dv):
    """
    Функция принимает путь сохранения рисунка и три возможных последовательности для построения двух или трех
    графиков с общим аргументом mean_frame_ind_pol, который приходит от вызывающей функции. При этом наличие
    двух отображаемых на рис. последовательностей обязательно, последовательность _s_i присутствует всегда.
    :param _path_to_fig_folder: Путь к папке для сохранения рисунка
    :param _s_i:
    :param _s_v:
    :param _s_dv:
    :return:
    """
    _pic_name = pic_name(_path_to_fig_folder, 0, 'png')
    _path_to_pic = Path(_path_to_fig_folder, _pic_name)
    if _s_v is None or _s_dv is None:
        fig, axes = plt.subplots(2, 1, figsize=(12, 12))
    else:
        fig, axes = plt.subplots(3, 1, figsize=(12, 15))
    _fig_folder = str(_path_to_fig_folder)
    title1, title2, title3 = title_func(_fig_folder, head)

    y_max = np.nanmax(_s_i)
    y_min = np.nanmin(_s_i)
    _i_max = len(num_mask)
    _i = 0

    for _j in range(np.shape(_s_i)[1]):
        f1 = freq_mask0[_i]
        text1 = 'f = ' + str(f1) + ' MHz'
        axes[0].plot(mean_frame_ind_pol, _s_i[:, _j], label=text1)
        if _s_v is not None and _s_dv is not None:
            axes[1].plot(mean_frame_ind_pol, _s_v[:, _j])
            axes[2].plot(mean_frame_ind_pol, _s_dv[:, _j])
        elif _s_v is not None and _s_dv is None:
            axes[1].plot(mean_frame_ind_pol, _s_v[:, _j])
        else:
            axes[1].plot(mean_frame_ind_pol, _s_dv[:, _j])
        axes[0].legend(loc='upper right')
        _i += 1

    axes[0].set_title('Stokes Parameters ' + title1, fontsize=20)
    axes[0].set_ylabel('Stokes_I')
    if _s_v is not None:
        axes[1].set_ylabel('Stokes_V', color='darkred')
    else:
        axes[1].set_ylabel('Stokes_V Deviation', color='darkred')
    axes[0].minorticks_on()
    axes[1].minorticks_on()
    if _s_v is not None and _s_dv is not None:
        axes[2].set_ylabel('Stokes_V Deviation', color='darkred')
        axes[2].minorticks_on()
        axes[2].grid()
        axes[2].grid(which='minor',
                     axis='x',
                     color='k',
                     linestyle=':')
    y1 = y_max - 2 * (y_max - y_min) / 10
    y2 = y_max - 3 * (y_max - y_min) / 10
    axes[0].text(0, y1, inform[0], fontsize=12)  # Разрешение по частоте
    axes[0].text(0, y2, inform[1], fontsize=12)  # Разрешение по времени

    axes[0].grid()
    axes[1].grid()
    axes[0].grid(which='minor',
                 axis='x',
                 color='k',
                 linestyle=':')
    axes[1].grid(which='minor',
                 axis='x',
                 color='k',
                 linestyle=':')
    plt.show()
    #                               ********************************
    #                        ************ Сохранение рисунка ****************
    fig.savefig(_path_to_pic)
    flag_save = save_question()
    if flag_save == 'no':
        if os.path.isfile(_path_to_pic):
            os.remove(_path_to_pic)
            print('Picture is not saved')
        else:
            print('File not found')
    else:
        print('Picture is saved')
    pass


def simplest_fig(_x, _y, _z):
    fig, axes = plt.subplots(2, 1, figsize=(12, 12))
    axes[0].plot(_x, _y)
    axes[1].plot(_x, _z)
    axes[0].set_ylabel('Stokes_I')
    axes[1].set_ylabel('Stokes_V', color='darkred')
    axes[0].minorticks_on()
    axes[1].minorticks_on()
    axes[0].legend(loc='upper right')
    axes[0].grid()
    axes[1].grid()
    axes[0].grid(which='minor',
                 axis='x',
                 color='k',
                 linestyle=':')
    axes[1].grid(which='minor',
                 axis='x',
                 color='k',
                 linestyle=':')
    plt.show()


def twin_fig_plot():
    """
    Строит графики на одном поле, но с разными осями по 0у
    :return:
    """
    for j in num_mask:
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()  # Создание второй оси ординат (справа). Ось абсцисс общая
        ax1.plot(mean_frame_ind_pol, s0[:, j], label='x(t)')
        ax2.plot(mean_frame_ind_pol, s3[:, j], label='y(t)', color='darkred')
        # ax1.plot([i for i in range(m)], s0[:, j], label='x(t)')
        # ax2.plot([i for i in range(m)], s3[:, j], label='y(t)', color='darkred')
        ax1.set_ylabel('Stokes_I')
        ax2.set_ylabel('Stokes_V', color='darkred')
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


def pic_name(file_path, flag, format='png'):
    """
    Функция принимает папку в которую надо сохранить картинки с общим названием, различающиеся по номеру
    после названия, флаг, который определяет название, формат сохранения. Возвращает название файла с номером
    и расширением
    :param file_path: папка сохранения
    :param flag: название картинки (вид картинки: скан, спектр и т.п.)
    :param format: формат сохранения
    :return: имя файла с расширением, под которым будет сохранен рисунок
    """
    if flag == 1:
        add_pass0 = 'spectrum_00'
    elif flag == 2:
        add_pass0 = 'colour2D_00'
    elif flag == 3:
        add_pass0 = 'pic3D_00'
    else:
        add_pass0 = 'scan_stokes_00'

    l = len(add_pass0)
    add_pass1 = add_pass0 + '.' + format
    if not os.path.isfile(Path(file_path, add_pass1)):
        pass
    else:
        while os.path.isfile(Path(file_path, add_pass1)):
            num = int(add_pass0[l - 2:l]) + 1
            num_str = str(num)
            if num >= 10:
                add_pass0 = add_pass0[:l - 2] + num_str
            else:
                add_pass0 = add_pass0[:l - 2] + '0' + num_str
            add_pass1 = add_pass0 + '.' + format

    return add_pass1


def save_question():
    root = Tk()
    answer = mb.askquestion(
        title="Save control",
        message="Save picture?")
    root.mainloop()
    del root
    return answer


def title_func(file_name0, _head):
    az = file_name0[-3:]
    att1 = str(_head['att1'])
    att2 = str(_head['att2'])
    att3 = str(_head['att3'])
    date = _head['date']

    title1 = date + ', az = ' + az + ', Att = [' + att1 + ', ' + att2 + ', ' + att3 + ']'
    a = file_name0.find('sun', -50, -1)
    if not file_name0.find('sun', -50, -1) == -1:
        title2 = 'Sun intensity'
        title02 = 'Sun spectrum '
        if file_name0[-1:] == 'b':
            az = file_name0[-4:-1]
            title2 = 'Calibration'
            title02 = 'Calibration spectrum '
            title1 = date + ', az = ' + az + ', Black Body/Sky, Att = [' + att1 + ', ' + att2 + ', ' + att3 + ']'

    elif not file_name0.find('crab') == -1:
        title2 = 'Crab intensity'
        title02 = 'Crab spectrum '

    elif not file_name0.find('moon') == -1:
        title2 = 'Moon intensity'
        title02 = 'Moon spectrum '

    elif not file_name0.find('3C84') == -1:
        title2 = '3C84 intensity'
        title02 = '3C84 spectrum '

    elif not file_name0.find('calibration') == -1:
        title2 = 'Calibration'
        title02 = 'Calibration spectrum '
        title1 = date + ', az = ' + az + ', Att = [' + att1 + ', ' + att2 + ', ' + att3 + ']'

    elif not (file_name0.find('test') == -1):
        kind = file_name0[-2:]
        title2 = 'Test'
        title02 = 'Test spectrum'
        if kind == 'VG':
            power_vg = 0
            title1 = date + ', Vector Gen' + ', P = ' + str(power_vg) + 'dBm, ' \
                                                                        'Att = [' + att1 + ', ' + att2 + ', ' + att3 + ']'
        elif kind == 'NG':
            t_noise = 6300
            title1 = date + ', Noise Gen, ' + 'T = ' + str(
                t_noise) + 'K, Att = [' + att1 + ', ' + att2 + ', ' + att3 + ']'
        elif kind == 'ML':
            title1 = date + ', Matched Load' ', Att = [' + att1 + ', ' + att2 + ', ' + att3 + ']'
        elif kind == 'SC':
            title1 = date + ', Short Cut' ', Att = [' + att1 + ', ' + att2 + ', ' + att3 + ']'
        else:
            title1 = date + ', ' + ' , Att = [' + att1 + ', ' + att2 + ', ' + att3 + ']'

    else:
        title2 = 'Scan'
        title02 = 'Spectrum'

    return title1, title2, title02


def time_to_angle(_time, _data, _path=None, _az=0):
    if _path:
        _scale = sun_az_speed(_path, _az)
    else:
        _scale = 1900 / 180
    _time_sc = 200
    _angle = [-(t - _time_sc) * _scale for t in _time][-1::-1]
    _data = _data[-1::-1, :]
    return _angle, _data


if __name__ == '__main__':
    align = 'y'
    channel_align = 'y'
    noise_int_calibration = 'n'
    v_deviation = 'n'
    object_m = 'sun'
    freq_mask_list = freq_mask(8)
    freq_mask0 = np.array(freq_mask_list)

    current_data_file = '2024-01-05_13-24'  # Имя файла с исходными текущими данными без расширения
    main_dir = current_data_file[0:4]  # Каталог всех данных (первичных, вторичных) за год
    current_primary_dir = f'{main_dir}_{current_data_file[5:7]}_{current_data_file[8:10] + object_m}'

    align_file_name: Any = 'antenna_temperature_coefficients.npy'  # Имя файла с текущими коэффициентами
    # выравнивания АЧХ
    dict_calibr_file_name = 'dict_calibr.csv'  # Имя файла c таймингом калибровок по ГШ и по поляризации

    path_obj = DataPaths(current_data_file, current_primary_dir, main_dir)
    converted_dir_path = path_obj.converted_dir_path
    treatment_dir_path = path_obj.treatment_dir_path
    head_path = path_obj.head_path

    path_to_stokes = Path(converted_dir_path, current_data_file + '_stocks.npy')
    path_stokes_base = Path(converted_dir_path, 'stokes_base.npy')
    path_to_stokes_left_txt = Path(converted_dir_path, current_data_file + '_left.txt')
    path_to_stokes_right_txt = Path(converted_dir_path, current_data_file + '_right.txt')
    path_to_stocks_fig_folder = Path(treatment_dir_path, current_data_file)
    path_to_csv = Path(converted_dir_path, dict_calibr_file_name)
    s = start_stop_calibr(current_data_file, path_to_csv)

    if not (os.path.isfile(path_to_stokes)):
        #               **********************************************
        # ************** Загрузка матрицы спектров и установок (head) *************
        with open(Path(converted_dir_path, current_data_file + '_head.bin'), 'rb') as inp:
            head = pickle.load(inp)
        if os.path.exists(Path(converted_dir_path, current_data_file + '_spectrum.npy.gz')):
            filename_out = Path(converted_dir_path, current_data_file + '_spectrum.npy.gz')
            with gzip.open(filename_out, "rb") as fin:
                spectrum = np.load(fin, allow_pickle=True)
        else:
            spectrum = np.load(Path(converted_dir_path, current_data_file + '_spectrum.npy'), allow_pickle=True)
        #               **********************************************

        if align == 'y':
            path = Path(head_path, 'Alignment', align_file_name)
            spectrum1 = align_spectrum(spectrum[0], spectrum[1], spectrum[2], spectrum[3], head,
                                       path, 2)

        input_data_upper = {'left': spectrum1[1],
                            'right': spectrum1[3]}
        input_data_lower = {'left': spectrum1[0],
                            'right': spectrum1[2]}

        a1 = input_data_upper['left']
        b1 = input_data_upper['right']

        a = input_data_lower['left']
        a = a[:, -1::-1]
        b = input_data_lower['right']
        b = b[:, -1::-1]
        mean_frame_ind_left, mean_frame_ind_right = initial_scan_cut(a)

        c = pol_intensity(a, mean_frame_ind_left)[1:]
        d = pol_intensity(b, mean_frame_ind_right)[: -1]
        c1 = pol_intensity(a1, mean_frame_ind_left)[1:]
        d1 = pol_intensity(b1, mean_frame_ind_right)[: -1]
        c = np.hstack((c, c1))
        d = np.hstack((d, d1))
        mean_frame_ind = frame_ind_extend()
        ind_c1 = [s[4] <= el <= s[5] for el in mean_frame_ind]
        #                               ****************
        # Вычисление выравнивающий коэффициентов по калибровочному сигналу - калибровочный сигнал д.б. неполяризованным
        if channel_align == 'y':
            av_c_cal = np.nanmean(c[ind_c1, :], axis=0)
            av_d_cal = np.nanmean(d[ind_c1, :], axis=0)
            noise_coeff = av_c_cal / av_d_cal
            m, n = np.shape(c)
            for j in range(n):
                d[:, j] = d[:, j] * noise_coeff[j]

            #                               *****************
        # np.savetxt(path_to_stokes_left_txt, c)
        # np.savetxt(path_to_stokes_right_txt, d)
        # Параметры Стокса
        s0 = c + d
        s3 = c - d

        stokes_coeff = pd.Series([s0, s3, mean_frame_ind, noise_coeff])
        np.save(path_to_stokes, stokes_coeff)
        print('Stokes parameters are saved successfully')

    else:
        if os.path.exists(f'{str(path_to_stokes)}.gz'):
            filename_out = f'{str(path_to_stokes)}.gz'
            with gzip.open(filename_out, "rb") as fin:
                _spectrum = np.load(fin, allow_pickle=True)
        else:
            _spectrum = np.load(path_to_stokes, allow_pickle=True)

        stokes_coeff = np.load(path_to_stokes, allow_pickle=True)
        [s0, s3, mean_frame_ind, equalizing_factor] = stokes_coeff

    mean_frame_ind_pol = np.copy(mean_frame_ind)
    m, n = np.shape(s0)
    freq_res = n  # Число отсчетов спектра шириной 2 ГГц по частоте
    df = 3.90625
    freq = [1000 + df / 2 + df * i for i in range(n)]
    num_mask = [int((s - 1000 * (1 + 1 / freq_res)) * freq_res / 2000) for s in freq_mask0]

    calibration_temperature = [calibration_temp(f) for f in freq_mask0]
    c = (s3 + s0) / 2  # левая поляризация

    #                         *************************************
    #   *************** Калибровка антенной температуры по внутреннему ГШ *************
    #                         *************************************
    if noise_int_calibration == 'y':
        ind_c = [s[0] <= el <= s[1] or s[2] <= el <= s[3] for el in mean_frame_ind]
        for j in range(np.size(freq_mask0)):
            av_c_cal = np.nanmean(c[ind_c, num_mask[j]])
            temp_coeff = calibration_temperature[j] / av_c_cal
            s0[:, num_mask[j]] = s0[:, num_mask[j]] * temp_coeff
            s3[:, num_mask[j]] = s3[:, num_mask[j]] * temp_coeff
    #                                   ******************

    with open(Path(converted_dir_path, current_data_file + '_head.bin'), 'rb') as inp:
        head = pickle.load(inp)
    inform = [('time resol = ' + str(60 * 8.3886e-3) + 'sec'),
              ('freq resol = ' + str(int(2000 // (n + 1))) + 'MHz'),
              ('polarisation ' + head['polar']),
              'align: ' + 'yes',
              'kurtosis quality = ' + str(head['good_bound'])]

    if v_deviation == 'y':
        s3_dv = stokes_v_deviation(s3[:, num_mask], 31)
    else:
        s3_dv = None
    current_treatment_dir = current_primary_dir + '_treat'
    current_treatment_path = Path('Data_treatment', current_treatment_dir)
    s0_selected = s3[:, num_mask]
    mean_frame_ind_pol = mean_frame_ind_pol * 8.3886e-3
    mean_frame_ind_pol, s0_selected = time_to_angle(mean_frame_ind_pol, s0_selected,
                                                    current_data_file[:10],
                                                    int(current_data_file[-3::]))
    #                   **********************************************
    #                       ****** Графический вывод данных ******
    #                   **********************************************
    # twin_fig_plot()                               # График с двумя разномасштабными осями 0у (слева и справа)
    # two_fig_plot(path_to_stocks_fig_folder)       # Картинка с двумя графиками (I & V)
    # mean_frame_ind_pol = mean_frame_ind_pol[670:770]
    s3f = np.ones(np.shape(s3))
    for i in range(np.shape(s3)[1]):
        s3f[:, i] = maf_fir(s3[:, i], 9)
    if v_deviation == 'y':
        some_fig_plot(path_to_stocks_fig_folder, s0[-1::-1, num_mask], s3[-1::-1, num_mask], s3_dv[-1::-1, :])
    else:
        some_fig_plot(path_to_stocks_fig_folder, s0[-1::-1, num_mask], s3f[-1::-1, num_mask], s3_dv)
    # fig_multi_axes(np.transpose(s0_selected), mean_frame_ind_pol, inform,
    #                Path(current_treatment_path, current_data_file), freq_mask0, head)
    pass
    # s31 = np.ma.masked_array(s3, np.isnan(s3)) + 2000
    # graph_contour_2d(freq, mean_frame_ind_pol, s31, 0, 'ab',  treatment_dir_path, head)
