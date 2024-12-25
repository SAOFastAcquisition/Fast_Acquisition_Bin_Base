from typing import Any
import numpy as np
import os
import gzip
import sys
import pickle
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import pandas as pd
from tkinter import *
from tkinter import messagebox as mb
from pathlib import Path
from scipy.signal import filtfilt
from Supporting_func.Fig_plot import fig_multi_axes, graph_contour_2d
from Help_folder.paths_via_class import DataPaths
from Supporting_func.sun_az_spead import sun_az_speed

# from Data_treatment.data_visualisation_sfu import zone_deletion

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


def stokes_coefficients(_path_obj, _current_data_file):
    _path_to_stokes = Path(_path_obj.converted_dir_path, _current_data_file + '_stocks.npy')

    if not (os.path.isfile(_path_to_stokes)):
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

        input_data_upper = {'left': spectrum[1],
                            'right': spectrum[3]}
        input_data_lower = {'left': spectrum[0],
                            'right': spectrum[2]}

        a1 = input_data_upper['left']
        b1 = input_data_upper['right']

        a = input_data_lower['left']
        a = a[:, -1::-1]
        b = input_data_lower['right']
        b = b[:, -1::-1]
        _mean_frame_ind_left, _mean_frame_ind_right = initial_scan_cut(a)

        c = pol_intensity(a, _mean_frame_ind_left)[1:]
        d = pol_intensity(b, _mean_frame_ind_right)[: -1]
        c1 = pol_intensity(a1, _mean_frame_ind_left)[1:]
        d1 = pol_intensity(b1, _mean_frame_ind_right)[: -1]
        c = np.hstack((c, c1))
        d = np.hstack((d, d1))
        _mean_frame_ind = frame_ind_extend(_mean_frame_ind_left, _mean_frame_ind_right, c)

        #                               ****************
        # Вычисление выравнивающий коэффициентов по калибровочному сигналу - калибровочный сигнал д.б. неполяризованным

        flux_coeff_left = head['flux_coeff_left']
        flux_coeff_right = head['flux_coeff_right']
        _n = len(flux_coeff_left)
        for j in range(_n):
            d[:, j] = d[:, j] * flux_coeff_right[j]
            c[:, j] = c[:, j] * flux_coeff_left[j]
        noise_coeff = flux_coeff_left / flux_coeff_right

        #                               *****************
        #                               Параметры Стокса
        _s0 = c + d
        _s3 = c - d

        stokes_coeff = pd.Series([_s0, _s3, _mean_frame_ind, noise_coeff])
        np.save(_path_to_stokes, stokes_coeff)
        print('Stokes parameters are saved successfully')

    else:
        if os.path.exists(f'{str(_path_to_stokes)}.gz'):
            filename_out = f'{str(_path_to_stokes)}.gz'
            with gzip.open(filename_out, "rb") as fin:
                _spectrum = np.load(fin, allow_pickle=True)
        else:
            _spectrum = np.load(_path_to_stokes, allow_pickle=True)

        stokes_coeff = np.load(_path_to_stokes, allow_pickle=True)

    return stokes_coeff


def frame_ind_extend(_mean_frame_ind_right, _mean_frame_ind_left, _c):
    """
    Функция объединяет последовательности отсчетов, соответствующих центрам полупериодов переключения левой и
    правой поляризаций в одну последовательность
    :return: int
    последовательность с чередующимися номерами отсчетов, соответствующих центрам полупериодов с левой
    и правой поляризациями
    """
    _mean_frame_ind = [0] * np.shape(_c)[0]
    _mean_frame_ind[::2] = _mean_frame_ind_right[:-1]
    _mean_frame_ind[1::2] = _mean_frame_ind_left[1:]
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
    _n2 = 1
    _freq_mask = [
        [1736],  # [0]
        [2060, 2300, 2500, 2750, 2830, 2920],  # [1]
        [1020, 1100, 1200, 1300, 1350, 1400, 1450, 1600],  # [2]
        [1000 * _n1 + 100 * _n2 + 0 + 10 * i for i in range(10)],  # [3]
        [1050, 1465, 1535, 1600, 1700, 2265, 2550, 2700, 2800, 2920],  # [4]
        [1230, 1560, 2300, 2910],  # [5]
        [1140, 1420, 1480, 2460, 2500, 2780],  # for Crab '2021-06-28_03+14' # [6]
        [1220, 1540, 1980, 2060, 2500, 2780],  # for Crab '2021-06-28_04+12' # [7]
        [1200, 1300, 1465, 1600, 1700, 2265, 2510, 2720, 2800, 2920]  # [8]
    ]
    return _freq_mask[_i]


def two_fig_plot1(_x, _y1, _y2, _dict_pic, _head):
    _pic_name = pic_name(_dict_pic['path_fig_folder'], _dict_pic['flag'], _dict_pic['pic_format'])
    _path_to_pic = Path(_dict_pic['path_fig_folder'], _pic_name)
    fig, axes = plt.subplots(2, 1, figsize=(10, 14))
    _fig_folder = str(_dict_pic['path_fig_folder'])
    title1, title2, title3 = title_func(_fig_folder, _head)

    y_max = np.nanmax(_y1)
    y_min = np.nanmin(_y1)

    _i_max = len(num_mask)
    _i = 0
    for _j, _txt in enumerate(_dict_pic['line_labels']):
        _y = _y1[:, _j]
        axes[0].plot(_x, _y1[:, _j], label=_txt)
        _t = _y2[:, _j]
        axes[1].plot(_x, _y2[:, _j])
        axes[0].legend(loc=_dict_pic['legend_pos'])
        # _i += 1
    plt.show()
    axes[0].set_title(_dict_pic['title'] + title1, fontsize=20)
    axes[1].set_xlabel(_dict_pic['xlabel_ax1'], color='darkred')
    axes[0].set_ylabel(_dict_pic['ylabel_ax0'])
    axes[1].set_ylabel(_dict_pic['ylabel_ax1'], color='darkred')
    axes[0].minorticks_on()
    axes[1].minorticks_on()

    axes[0].legend(bbox_to_anchor=(1, 1), loc="upper left")

    y1 = y_max - 2 * (y_max - y_min) / 10
    y2 = y_max - 3 * (y_max - y_min) / 10
    y3 = y_max - 4 * (y_max - y_min) / 10
    axes[0].text(-2500, y1, inform[0], fontsize=12)  # Разрешение по частоте
    axes[0].text(-2500, y2, inform[1], fontsize=12)  # Разрешение по времени
    axes[0].text(-2500, y3, inform[3], fontsize=12)  # Калибровка
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


def legend_iter(_leg):
    yield from _leg


def mod_plt(_key, _obj, *_args, **_kwargs):
    dict_plt = {
        'linear': _obj.plot,
        'logy': _obj.semilogy
    }
    op_plot = dict_plt[_key]
    return op_plot(*_args, **_kwargs)


def two_fig_plot(*_args, _x1, _y3):
    """
    Функция выдает рисунок, состоящий из трех графиков: первые два - зависимости IV параметров Стокса или левой-правой
    поляризаций от частоты (спектры) или положения ДН телескопа на Солнце (сканы), третий - вспомогательный
    (информативный) показывает либо точки на скане к которым относится спектр, либо спектры в выбранных точках сканов
    :param _args:   _args[0] - аргумент для двух основных рисунков
                    _args[1] - ордината первого графика
                    _args[2] - ордината второго графика
                    _args[-2] - параметры графиков: заголовки, подписи к осям, единицы и т.д.
                    _args[-1] - параметры сбора и обработки данных
    :param _x1: аргумент третьего графика
    :param _y3: ордината третьего графика
    :return:
    """
    _x = _args[0]
    _y1 = _args[1]
    _y2 = _args[2]
    _dict_pic = _args[-2]
    _head = _args[-1]
    _pic_name = pic_name(_dict_pic['path_fig_folder'], _dict_pic['flag'], _dict_pic['pic_format'])
    _path_to_pic = Path(_dict_pic['path_fig_folder'], _pic_name)

    _fig = plt.figure(figsize=(10, 14))
    _fig_folder = str(_dict_pic['path_fig_folder'])
    title1, title2, title3 = title_func(_fig_folder, _head)
    gs = GridSpec(ncols=5, nrows=4, figure=_fig)

    _y_max1 = np.nanmax(_y1)
    _y_min1 = np.nanmin(_y1)
    _y_max2 = np.nanmax(_y2)
    _y_min2 = np.nanmin(_y2)

    _i_max = len(num_mask)
    _axes0 = plt.subplot(gs[0:2, 0:3])
    _axes1 = plt.subplot(gs[2:, 0:3])
    _label1 = legend_iter(_dict_pic['line_labels'])

    _key1 = _dict_pic['key1']
    _key2 = _dict_pic['key2']

    for s1, s2 in zip(_y1, _y2):
        _axes0.plot(_x, s1, label=next(_label1))
        # mod_plt(_key1, _axes0, _x, s1, label=next(_label1))
        _axes1.plot(_x, s2)
        _axes0.legend(loc=_dict_pic['legend_pos'])

    _pos = _dict_pic.get('pos_select')
    _pos_s = _dict_pic.get('pos_spectrum_select')
    if _pos:
        for s in _pos:
            _axes0.plot([s, s], [_y_min1, _y_max1],
                        '--', label=f'pos = {s} arcs')
            _axes1.plot([s, s], [_y_min2, _y_max2], '--')

    _axes0.set_title(_dict_pic['title'] + title1, fontsize=20)
    _axes1.set_xlabel(_dict_pic['xlabel_ax1'], color='darkred')
    _axes0.set_ylabel(_dict_pic['ylabel_ax0'])
    _axes1.set_ylabel(_dict_pic['ylabel_ax1'], color='darkred')
    _axes0.minorticks_on()
    _axes1.minorticks_on()

    _axes0.legend(bbox_to_anchor=(1, 1), loc="upper left")

    if len(_x1):
        _axes2 = plt.subplot(gs[2:, 3:])
        a = legend_iter(_dict_pic['line_labels_add'])
        for s in _y3:
            _axes2.plot(_x1, s, label=next(a))
            # mod_plt(_key2, _axes2, _x1, s, label=next(a))
        if _pos_s:
            _y_max1 = np.nanmax(_y3)
            for s in _pos_s:
                _axes2.plot([s, s], [_y_min1, _y_max1],
                            '--', label=f'pos = {s} arcs')
        _axes2.set_title(_dict_pic['title_add'], fontsize=14)
        _axes2.set_ylabel(_dict_pic['ylabel_ax2'])
        _axes2.set_xlabel(_dict_pic['xlabel_ax2'])
        _axes2.legend(bbox_to_anchor=(0.45, 1.1), loc="lower left")
        _axes2.minorticks_on()
        _axes2.grid()
        _axes2.grid(which='minor',
                    axis='x',
                    color='k',
                    linestyle=':')
        _axes2.tick_params(axis='both',  # Применяем параметры к обеим осям
                           which='major',  # Применяем параметры к основным делениям
                           direction='in',  # Рисуем деления внутри и снаружи графика
                           length=5,  # Длинна делений
                           width=1,  # Ширина делений
                           color='black',  # Цвет делений
                           pad=2,  # Расстояние между черточкой и ее подписью
                           # labelsize=f_size,  # Размер подписи
                           labelcolor='black',  # Цвет подписи
                           bottom=True,  # Рисуем метки снизу
                           top=False,  # сверху
                           left=True,  # слева
                           right=True,  # и справа
                           labelbottom=True,  # Рисуем подписи снизу
                           labeltop=False,  # сверху
                           labelleft=False,  # слева
                           labelright=True,  # и справа
                           labelrotation=0)  # Поворот подписей

    y1 = _y_max1 - 2 * (_y_max1 - _y_min1) / 10
    y2 = _y_max1 - 3 * (_y_max1 - _y_min1) / 10
    y3 = _y_max1 - 4 * (_y_max1 - _y_min1) / 10
    _axes0.text(-2500, y1, inform[0], fontsize=12)  # Разрешение по частоте
    _axes0.text(-2500, y2, inform[1], fontsize=12)  # Разрешение по времени
    _axes0.text(-2500, y3, inform[3], fontsize=12)  # Калибровка
    _axes0.grid()
    _axes1.grid()
    _axes0.grid(which='minor',
                axis='x',
                color='k',
                linestyle=':')
    _axes1.grid(which='minor',
                axis='x',
                color='k',
                linestyle=':')

    plt.show()
    #                               ********************************
    #                        ************ Сохранение рисунка ****************
    _fig.savefig(_path_to_pic)
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
    elif flag == 4:
        add_pass0 = 'scan_LR_polar_00'
    elif flag == 5:
        add_pass0 = 'spectrum_LR_polar_00'
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
    date = _head['_date']

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
        if _scale == 1:
            _scale = 1900 / 135
    else:
        _scale = 1900 / 135
    _time_sc = 200
    _angle = [-(t - _time_sc) * _scale for t in _time][-1::-1]
    _data = _data[-1::-1, :]
    return _angle, _data


def zone_deletion(_len):
    if True:
        # Исключение зон действия режекторных фильтров при правильном порядке отсчетов частоты во второй зоне Найквиста
        _df = 2000 / _len
        _f = np.array([1000 + _df / 2 + i * _df for i in range(_len)])
        _int = [s < 1025 or
                1770 < s < 2034 or
                2090 < s < 2230 or
                2525 < s < 2710 or
                s > 2954 for s in _f]

    return _int


if __name__ == '__main__':
    object_m = 'sun'
    current_data_file = '2024-03-11_01+24'  # Имя файла с исходными текущими данными без расширения

    freq_mask0 = np.array(freq_mask(8))  # Маска частот для основного рисунка со сканами
    pos_select = [-300, 250]           # Выбор позиций на Солнце для вспомогательного графика рисунка со сканами
    angle_mask = [-1000, -300, 1000]    # Выбор позиций на Солнце для основного рисунка со спектрами
    freq_select = [1500, 2300, 2850]    # Маска частот для вспомогательного графика рисунка со спектрами
    pic = 'LR'  # Выбор объектов для рисунка со сканами: 'LR' - левая и правая поляризации, 'IV' - параметры Стокса
    pic_spectrum = 'y'  #

    main_dir = current_data_file[0:4]  # Каталог всех данных (первичных, вторичных) за год
    current_primary_dir = f'{main_dir}_{current_data_file[5:7]}_{current_data_file[8:10] + object_m}'
    path_obj = DataPaths(current_data_file, current_primary_dir, main_dir)
    converted_dir_path = path_obj.converted_dir_path
    treatment_dir_path = path_obj.treatment_dir_path
    head_path = path_obj.head_path
    path_to_stocks_fig_folder = Path(treatment_dir_path, current_data_file)

    [s0, s3, mean_frame_ind, equalizing_factor] = stokes_coefficients(path_obj, current_data_file)
    mean_frame_ind_pol = np.copy(mean_frame_ind)

    m, n = np.shape(s0)  # n - Число отсчетов спектра шириной 2 ГГц по частоте
    df = 2000 / n  # разрешение по частоте
    freq = [1000 + df / 2 + df * i for i in range(n)]
    num_mask = [int((s - 1000 * (1 + 1 / n)) / df) for s in freq_mask0]

    # Исключение зон режекторных фильтров
    n_del = zone_deletion(len(freq))
    s0[:, n_del] = None
    s3[:, n_del] = None
    # Инверсия временных отсчетов для перехода к позиционным на Солнце
    s0 = s0[-1::-1, :]  # Stokes I
    s3 = s3[-1::-1, :]  # Stoked V

    c = (s3 + s0) / 2  # левая поляризация
    d = (s0 - s3) / 2  # правая поляризация
    #                                   ******************

    with open(Path(converted_dir_path, current_data_file + '_head.bin'), 'rb') as inp:
        head = pickle.load(inp)
    inform = [f'time resol = {60 * 8.3886e-3: .2f} sec',
              f'freq resol = {df: 1.2f} MHz',
              ('polarisation ' + head['polar']),
              'align: quiet Sun',
              'kurtosis quality = ' + str(head['good_bound'])]

    s0_selected = s3[:, num_mask]
    mean_frame_time = mean_frame_ind_pol * 8.3886e-3
    mean_frame_pos, s0_selected = time_to_angle(mean_frame_time, s0_selected,
                                                current_data_file[:10],
                                                int(current_data_file[-3::]))
    mean_f_pos = np.array(mean_frame_pos)
    pos_select_num = [np.where(mean_f_pos > s)[0][0] for s in pos_select]

    #                   **********************************************
    #                       ****** Графический вывод данных ******
    #                   **********************************************

    dict_pic_IV = {
        'pic_format': 'png',
        'flag': 0,  # Сохранение рисунка с подписью 'scan_stokes'
        'title': 'Stokes Parameters ',
        'title_add': 'Stokes IV spectrum, s.f.u. ',
        'key1': 'linear',
        'key2': 'logy',
        'ylabel_ax0': 'Stokes_I, s.f.u.',
        'ylabel_ax1': 'Stokes_V, s.f.u.',
        'ylabel_ax2': 'L&R polarizations, s.f.u.',
        'xlabel_ax1': 'Sun angle, arcs',
        'xlabel_ax2': 'Frequencies, MHz',
        'line_labels': [f'f = {f1} MHz' for f1 in freq_mask0],
        'legend_pos': 'upper right',
        'legend_add_pos': [0.6, 0.9],
        'pos_select': pos_select,
        'line_labels_add': [f'pos = {a} arcs - I' for a in pos_select] +
                           [f'pos = {a} arcs - V' for a in pos_select],
        'path_fig_folder': path_to_stocks_fig_folder
    }

    dict_pic_LR = {
        'pic_format': 'png',
        'flag': 4,  # Сохранение рисунка с подписью 'scan_LR_polar'
        'title': 'LR polarization ',
        'title_add': 'L&R spectrum, s.f.u. ',
        'key1': 'linear',
        'key2': 'logy',
        'ylabel_ax0': 'Left polarization, s.f.u.',
        'ylabel_ax1': 'Right polarization, s.f.u.',
        'ylabel_ax2': 'L&R polarizations, s.f.u.',
        'xlabel_ax1': 'Sun angle, arcs',
        'xlabel_ax2': 'Frequencies, MHz',
        'line_labels': [f'f = {f1} MHz' for f1 in freq_mask0],
        'legend_pos': 'upper right',
        'legend_add_pos': [0.6, 0.9],
        'pos_select': pos_select,
        'line_labels_add': [f'pos = {a} arcs - left' for a in pos_select] +
                           [f'pos = {a} arcs - right' for a in pos_select],
        'path_fig_folder': path_to_stocks_fig_folder
    }

    angle_num = [np.where(mean_f_pos > a)[0][0] for a in angle_mask]
    pos_sp_select = []

    freq_num_sel = [int((s - 1000 * (1 + 1 / n)) / df) for s in freq_select]

    dict_pic_spectrum = {
        'pic_format': 'png',
        'flag': 5,  # Сохранение рисунка с подписью 'spectrum_LR_polar'
        'title': 'Spectrum LR polarization ',
        'title_add': 'Stokes I scan, s.f.u. ',
        'key1': 'logy',
        'key2': 'linear',
        'ylabel_ax0': 'Left polarization, s.f.u.',
        'ylabel_ax1': 'Right polarization, s.f.u.',
        'ylabel_ax2': 'Stokes I, s.f.u.',
        'xlabel_ax1': 'Frequency, MHz',
        'xlabel_ax2': 'Sun position, arcs',
        'line_labels': [f'pos = {f1} arcs' for f1 in angle_mask],
        'legend_pos': 'upper right',
        'pos_spectrum_select': angle_mask,
        'line_labels_add': [f'freq = {a} MHz' for a in freq_select],
        'path_fig_folder': path_to_stocks_fig_folder
    }

    if pic == 'LR':
        args = [mean_frame_pos, c[:, num_mask].T, d[:, num_mask].T,
                dict_pic_LR, head]
        kwargs = {
            '_x1': freq,
            '_y3': np.vstack((c[pos_select_num, :], d[pos_select_num, :]))
        }
    elif pic == 'IV':
        args = [mean_frame_pos, s0[:, num_mask].T, s3[:, num_mask].T,
                dict_pic_IV, head]
        kwargs = {
            '_x1': freq,
            '_y3': np.vstack((s0[pos_select_num, :], s3[pos_select_num, :]))
        }

    two_fig_plot(*args, **kwargs)  # Картинка с двумя графиками (L & R)
    #                   *******************************************

    ord1 = c[angle_num, :]  #
    ord2 = d[angle_num, :]  #
    if pic_spectrum == 'y':
        args = [freq, ord1, ord2, dict_pic_spectrum, head]
        kwargs = {
            '_x1': mean_f_pos,
            '_y3': s0[:, freq_num_sel].T
        }
        two_fig_plot(*args, **kwargs)  # Картинка с двумя графиками (L & R)

    # fig_multi_axes(np.transpose(s0_selected), mean_frame_time, inform,
    #                Path(path_obj.treatment_data_file_path, current_data_file), freq_mask0, head)
    pass
    # s31 = np.ma.masked_array(s3, np.isnan(s3))
    # graph_contour_2d(freq, mean_frame_angle, s31[-1::-1, :], 0, 'ab',  treatment_dir_path, head)
