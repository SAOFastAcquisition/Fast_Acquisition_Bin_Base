import numpy as np
import os
import sys
import pandas as pd
import pickle
import gzip
import json as jsn
from datetime import datetime
from pathlib import Path
import matplotlib.pyplot as plt
from numpy.array_api import astype

from Supporting_func import Fig_plot as fp, align_spectrum, path_to_data
from Supporting_func.dict_calibr_from_csv import start_stop_calibr
from flux_dm import interpol_sfu
# from Supporting_func import align_spectrum, path_to_data
from Polyphase import low_freq_noise_spectrum, plot_low_freq_spec, plot_low_freq_spec_ab
from Polyphase.cic_filter import signal_filtering
from test_statistic import low_noise_spectra_base
from polys3d_demo import poly_graph3d
from Help_folder.paths_via_class import DataPaths

current_dir = Path.cwd()
home_dir = Path.home()

sys.path.insert(0, Path(current_dir, 'Supporting_func'))
sys.path.insert(0, Path(current_dir, 'Interface'))
sys.path.insert(0, Path(current_dir, 'Polyphase'))
start = datetime.now()

freq_spect_mask = [1171, 1380, 1465, 1535, 1600, 1700, 2265, 2550, 2700, 2800, 2920]
time_spect_mask = [47, 84.4, 104, 133, 133.05, 177.02, 177.38]


def parts_to_numpy(list_arr, len_list):
    """ Функция превращает список в массив numpy.
    Разбивает файл на меньшие части и обрабатывает их поотдельности. По ходу завершинея обработки частей
    происходит объединение обработанных частей."""
    n = int(len_list // 1e5)
    k = int(len_list % 1e5)
    numpy_arr = []
    for i in range(n + 1):
        if i == n:
            auxiliary = list_arr[int(i * 1e5):int(i * 1e5 + k)]
            auxiliary = np.array(auxiliary, dtype='int64')
        else:
            auxiliary = list_arr[int(i * 1e5):int((i + 1) * 1e5)]
            auxiliary = np.array(auxiliary, dtype='int64')
        l = np.size(numpy_arr)
        if l:
            numpy_arr = np.vstack([numpy_arr, auxiliary])
        else:
            numpy_arr = auxiliary

    return numpy_arr


def cut_spectrum(spectrum, n_aver):
    spectrum.pop(-1)
    n_frame_last = spectrum[-1][0]
    rest = (n_frame_last + 1) % 2 ** (6 - n_aver)
    if rest:
        for k in range(rest):
            spectrum.pop(-1)
    print(n_frame_last, spectrum[-1][0])
    return spectrum


def line_legend(freq_mask):
    N_col_leg = len(freq_mask)
    N_row_leg = len(time_spect_mask)
    legend_freq = [0] * N_col_leg
    legend_time = [0] * N_row_leg
    i1 = 0
    for i in freq_mask:
        legend_freq[i1] = str(i) + ' MHz'
        i1 += 1
    i1 = 0
    for i in time_spect_mask:
        legend_time[i1] = str(i) + ' sec'
        i1 += 1

    return legend_time, legend_freq


def form_spectr_sp1(spectr_extr, freq_spect_mask_in=freq_spect_mask, time_spect_mask_in=time_spect_mask):
    """ Возвращает s_freq - срезы частотного спектра в моменты времени time_spect_mask и s_time - сканы Солнца
    по времени на частотах freq_spect_mask с заданным разрешением по времени и частоте

    """
    ind_spec = []
    ind_time = []
    t_ng = 1
    N_col = np.shape(spectr_extr)[1]
    s_freq = np.zeros((len(time_spect_mask_in), N_col // kf))
    s_time = np.zeros((N_row // kt, len(freq_spect_mask_in)))
    j = 0
    for f in freq_spect_mask_in:
        if band_size_init == 'half':
            ind1 = (f - (N_Nyq - 1) * 1000 - delta_f / aver_param / 2) // (delta_f / aver_param)
        elif band_size_init == 'whole':
            ind1 = (f - 1000 - delta_f / aver_param / 2) // (delta_f / aver_param)
        ind = int(ind1)
        if ind > N_col - int(kf / 2) - 1:
            ind = N_col - int(kf / 2) - 1
        if ind < int(kf / 2):
            ind = int(kf / 2)
        i = 0
        while kt * (i + 1) < N_row:
            if kf == 1:
                s_time[i, j] = np.sum(
                    spectr_extr[i * kt:(i + 1) * kt, ind][spectr_extr[i * kt:(i + 1) * kt, ind] > 0.03])
                n_mesh = (spectr_extr[i * kt: (i + 1) * kt, ind] > 0.03).sum()
                if n_mesh == 0:
                    s_time[i, j] = 0.02
                else:
                    s_time[i, j] /= n_mesh
            else:
                s_time[i, j] = np.sum(spectr_extr[i * kt:(i + 1) * kt, ind - int(kf / 2):ind +
                                                                                         int(kf / 2)][
                                          spectr_extr[i * kt:(i + 1) * kt, ind - int(kf / 2):ind + int(kf / 2)] > 0.03])
                n_mesh = (spectr_extr[i * kt: (i + 1) * kt, ind - int(kf / 2):ind + int(kf / 2)] > 0.03).sum()
                if n_mesh == 0:
                    s_time[i, j] = 0.02
                else:
                    s_time[i, j] /= n_mesh
            i += 1
        ind_spec.append(ind)
        j += 1
    i = 0
    for t in time_spect_mask_in:
        ind = int(t // delta_t)
        if ind > N_row - kt / 2 - 1:
            ind = N_row - int(kt / 2) - 1
        if ind < (kt / 2):
            ind = int(kt / 2)
        j = 0
        while (j + 1) * kf < N_col:
            if kt == 1:
                s_freq[i, j] = np.sum(
                    spectr_extr[ind, j * kf:(j + 1) * kf][spectr_extr[ind, j * kf:(j + 1) * kf] > 0.03])
                n_mesh = (spectr_extr[ind, j * kf:(j + 1) * kf] > 0.03).sum()
                if n_mesh == 0:
                    s_freq[i, j] = 0.02
                else:
                    s_freq[i, j] /= n_mesh
            else:
                s_freq[i, j] = np.sum(spectr_extr[ind - int(kt / 2):ind + int(kt / 2), j * kf
                                                                                       :(j + 1) * kf][
                                          spectr_extr[ind - int(kt / 2):ind + int(kt / 2), j * kf:(j + 1) * kf] > 0.03])
                n_mesh = (spectr_extr[ind - int(kt / 2):ind + int(kt / 2), j * kf:(j + 1) * kf] > 0.03).sum()
                if n_mesh == 0:
                    s_freq[i, j] = 0.02
                else:
                    s_freq[i, j] /= n_mesh

            j += 1
        ind_time.append(ind)
        i += 1
    s_time = s_time.transpose()
    if head['att3'] == 0:
        a = 5.27e8
    elif head['att3'] == 5:
        a = 6.21e8
    else:
        a = 5.8e8
    a = 1
    return s_freq * (2 ** shift) * t_ng / a, s_time * (2 ** shift) * t_ng / a


def spectr_construction(Spectr, kf, kt):
    ''' Функция формирует спектр принятого сигнала с требуемым разрешением по частоте и времени. Максимальное
    разрешение отсчетов по времени 8192 мкс и по частоте 7,8125 МГц. Путем суммирования и усреднерия по kt*kf
    отсчетам разрешение по частоте и по времени в исходном спектре Spectr уменьшается в kf и kt раз,
    соответственно. Преобразованный спектр возвращается как S1. Для трехмерного изображения
    '''

    N_col1 = N_col // kf
    N_row1 = N_row // kt
    S1 = np.zeros((N_row1, N_col1))

    for i in range(N_row1):
        for j in range(N_col1):
            try:
                S1[i, j] = np.sum(Spectr[i * kt: (i + 1) * kt, j * kf: (j + 1) * kf])
                N_mesh = (Spectr[i * kt: (i + 1) * kt, j * kf: (j + 1) * kf] > 0.03).sum()
                if N_mesh == 0:
                    S1[i, j] = 2
                else:
                    S1[i, j] = S1[i, j] / N_mesh
                if S1[i, j] == 0:
                    S1[i, j] = 2
                # if (j > 3) & (S1[i, j] > 1.5 * np.sum(S1[i, j-3:j])//3):
                #     S1[i, j] = np.sum(S1[i, j-3:j])//3
                # if robust_filter == 'y':
                #     a = param_robust_filter
                #     if (i > 3) & (S1[i, j] < 1 / a * np.sum(S1[i - 3:i - 1, j]) // 2):
                #         S1[i, j] = np.sum(S1[i - 1, j])
                #     if (i > 3) & (S1[i, j] > a * np.sum(S1[i - 3:i - 1, j]) // 2):
                #         # print(S1[i - 3:i+1, j])
                #         S1[i, j] = np.sum(S1[i - 1, j])
                #         # print(S1[i, j])
                #         pass

            except IndexError as allert_message:
                print(allert_message, 'ind i = ', i, 'ind j = ', j)
                pass
            except ValueError as value_message:
                print(value_message, 'ind i = ', i, 'ind j = ', j)
                pass

    return S1  # // kt // kf


def path_to_fig(_path):
    """ Создает директорию для рисунков обрабатываемого наблюдения, если она до этого не была создана,
    название директории  совпадает с названием исходного файла данных наблюдения
    """
    if not os.path.isdir(_path):
        os.mkdir(_path)
    return


def preparing_data():
    """ Функция в зависимости от вида данных (полная полоса 1-3 ГГц, половинная полоса 1-2 или 2-3 ГГц,
    с двумя поляризациями или одной) выдает данные для построения графиков"""

    # Для полосы 1-3 ГГц и двух возможных поляризаций выдает по два спектра (1-2 и 2-3 ГГц) для каждой поляризации.
    # Если поляризация не задействована, то соответствующие спектры - пустые. Спектр 1-2 ГГц - в обратном порядке
    _path1 = Path(str(converted_data_file_path) + '_spectrum.npy')

    if os.path.exists(f'{str(_path1)}.gz'):
        filename_out = f'{str(_path1)}.gz'
        with gzip.open(filename_out, "rb") as fin:
            _spectrum = np.load(fin, allow_pickle=True)
    else:
        _spectrum = np.load(_path1, allow_pickle=True)

    with open(Path(str(converted_data_file_path) + '_head.bin'), 'rb') as inp:
        _head = pickle.load(inp)
    _n_aver = _head['n_aver']
    _band_size = _head['band_size']
    _polar = _head['polar']

    # Разделяем составляющие  записи в полной полосе и с возможными двумя поляризациями
    if num_of_polar == 2 and _band_size == 'whole':
        spectrum_left1 = np.array(_spectrum[0], dtype='int64')
        spectrum_left2 = np.array(_spectrum[1], dtype='float32')
        spectrum_right1 = np.array(_spectrum[2], dtype='int64')
        spectrum_right2 = np.array(_spectrum[3], dtype='float32')
    # Выдает спектры для левой и правой поляризаций шириной по 1 ГГц. При нумерации спектров учитывается
    # значение зоны Найквиста. С индексом "1" - 1-2 ГГц, с индексом "2" - 2-3 ГГц, как и для случая выше.
    # На выходе формально все 4 спектра, но для незадействованной полосы они пустые
    elif num_of_polar == 2 and _band_size == 'half':
        if N_Nyq == 2:
            spectrum_left1 = _spectrum[0]
            spectrum_right1 = _spectrum[1]
            spectrum_left2 = []
            spectrum_right2 = []
        else:
            spectrum_left2 = _spectrum[0]
            spectrum_right2 = _spectrum[1]
            spectrum_left1 = []
            spectrum_right1 = []
    pass
    return spectrum_left1, spectrum_left2, spectrum_right1, spectrum_right2, int(_n_aver), _band_size, _polar


def unite_spectrum(spec):
    spec1 = spec[0]
    spec2 = spec[1]
    spec3 = spec[2]
    spec4 = spec[3]
    ind = []
    if np.size(spec1) and np.size(spec2):
        n_row = np.min([np.shape(spec1)[0], np.shape(spec2)[0]])
        spec1 = spec1[:n_row]
        spec2 = spec2[:n_row]
        spec_left = np.hstack((spec1, spec2))
        ind.append('left_whole')
    elif np.size(spec1):
        spec_left = spec1
        ind.append('left_half1')
    elif np.size(spec2):
        spec_left = spec2
        ind.append('left_half2')
    else:
        spec_left = []
    if np.size(spec3) and np.size(spec4):
        n_row = np.min([np.shape(spec3)[0], np.shape(spec4)[0]])
        spec3 = spec3[:n_row]
        spec4 = spec4[:n_row]
        spec_right = np.hstack((spec3, spec4))
        ind.append('right_whole')
    elif np.size(spec3):
        spec_right = spec3
        ind.append('right_half1')
    elif np.size(spec4):
        spec_right = spec4
        ind.append('right_half2')
    else:
        spec_right = []
    if np.size(spec_left) and np.size(spec_right):
        shape_l = np.shape(spec_left)
        shape_r = np.shape(spec_right)
        n_row1 = np.max([shape_l[0], shape_r[0]])
        spec_left0 = np.full((n_row1, shape_l[1]), 2, dtype='float32')
        spec_right0 = np.full((n_row1, shape_r[1]), 2, dtype='float32')
        spec_left0[:shape_l[0]] = spec_left
        spec_right0[:shape_r[0]] = spec_right
        united_spec = pd.Series([spec_left0, spec_right0], ind)
    elif np.size(spec_left):
        united_spec = pd.Series([spec_left.astype('float32')], ind)
    else:
        united_spec = pd.Series([spec_right.astype('float32')], ind)
    print('Spectra are united')
    return united_spec


def del_random_mod(_s, _s0, _ds=2):
    """
    Функция принимает одномерный массив отсчетов и ограничивает скачки одного отсчета по отношению к соседним справа
    и слева. Если скачок i-того отсчета по отношению к отсчету i-1  больше в _ds раз, чем отсчета i+1 к отсчету i-1,
    то i-тый отсчет заменяется средним отсчетов i-1 и i+1.
    Кроме того, если отсчет принимает значение меньше порогового, то функция присваивает ему значение порога
    :param _s: исходный одномерный массив
    :param _s0: пороговое значение
    :param _ds: коэффициент ограничения скорости роста
    :return:
    """
    _l = len(_s)
    for _i in range(1, _l - 3):
        if abs(_s[_i] - _s[_i - 1]) > _ds * abs(_s[_i + 2] - _s[_i - 1]):
            _s[_i] = (_s[_i - 1] + _s[_i + 2]) / 2
        if _s[_i] <= _s0:
            _s[_i] = _s0
    _s[0] = _s0
    _s[_l - 1:] = _s0
    return _s


def treatment_null_mesh(_s, _n, _s0=0):
    """
    Функция заменяет элемент последовательности _s средним значением по промежутку 2*_n в окрестности
    заменяемого элемена из _s, когда элемент меньше порога _s0
    """
    _l = np.size(_s)
    for i in range(_l):
        if _s[i] < _s0:
            if _n <= i < _l - _n - 1:
                _s_loc = _s[i - _n: i + _n + 1][_s[i - _n: i + _n + 1] > 0]
                pass
            elif i < _n:
                _s_loc = _s[0: i + _n + 1][_s[0: i + _n + 1] > 0]
            elif i > _l - _n - 1:
                _s_loc = _s[i - _n - 2: -1][_s[i - _n - 2: -1] > 0]

            try:
                s_loc_av = np.mean(_s_loc)
                _s[i] = s_loc_av
            except ValueError:
                pass
                _s[i] = 1
    return _s


def freq_mask(_i):
    _n1 = 1
    _n2 = 0
    _freq_mask = [
        [2424],  # [0]
        [1245, 1375, 2500, 2820],  # [1] article to 'ab' Crab and 3C273
        [1080, 1140, 1360, 1420, 1620, 1780, 1980],  # [2]
        [1000 * _n1 + 100 * _n2 + 20 + 20 * i for i in range(10)],  # [3]
        [1050, 1465, 1535, 1600, 1700, 2265, 2550, 2700, 2800, 2920],  # [4]
        [1245, 1375, 2260, 2360, 2500, 2720, 2820, 2940],  # [5]
        [1140, 1420, 1480, 2460, 2500, 2780],  # for Crab '2021-06-28_03+14'    # [6]
        [1220, 1540, 1980, 2060, 2500, 2780],  # for Crab '2021-06-28_04+12'    # [7]
        [1200, 1380, 1465, 1600, 1700, 2265, 2490, 2720, 2800, 2860],  # [8]
        [1200, 1380, 2265, 2800]  # [9] article to 'ab' Sun
    ]
    return _freq_mask[_i]


def data_poly3d_prep(_spectrum_extr):
    _k, _m = np.shape(_spectrum_extr)
    freq_mask_poly = [1100, 1200, 1700, 1800, 2800, 2900]  # Маска для выбора спектров из трех диапазонов частот
    freq_mask_singles = [1100, 1200, 1300, 1700, 1780, 2080, 2300, 2500, 2720, 2800, 2920]  # Маска для одиночных частот
    mask_poly = 'y'
    _data_poly3d = []
    num_mask = []

    if mask_poly == 'y':
        i = 0
        _n = [0] * 6
        for s in freq_mask_poly:
            _n[i] = int((s - 1000) / 2000 * _m)
            if i % 2 == 1:
                if len(num_mask) == 0:
                    num_mask = [n for n in range(_n[i - 1], _n[i] + 1)]
                else:
                    num_mask = np.hstack((num_mask, [n for n in range(_n[i - 1], _n[i] + 1)]))
            i += 1
    else:
        num_mask = [int((s - 1000) / 2000 * _m) for s in freq_mask_singles]

    _data_poly3d = _spectrum_extr[:, num_mask]
    _freq_mask = freq[num_mask]

    return _data_poly3d, _freq_mask


def sun_calibration(_data, _file_name, _adr1):
    """

    :param _file_name: имя исходного первичного файла без расширения
    :param _adr1: объект адресов файловой системы данных
    :param _data: данные, которые надо откалибровать - спектры
    :return:
    """
    # Индексы отсчетов начала и конца калибровок. Инд. "4" и "5" соответствуют участку со спокойным Солнцем
    _path_to_csv = Path(_adr1.converted_dir_path, 'dict_calibr.csv')
    _path_to_head = Path(str(_adr1.converted_data_file_path) + '_head.bin')
    _s = start_stop_calibr(_file_name, _path_to_csv)
    with open(_path_to_head, 'rb') as _inp:
        _head = pickle.load(_inp)
    _l = 0
    for _sp in _data:
        print(_sp.dtype)
        _len_time, _len_freq = np.shape(_sp)
        _sfu = interpol_sfu(_len_freq)
        # Маска зон режекции аналогового тракта
        _del_zone = zone_deletion(_len_freq)
        # Индексы отсчетов начала и конца зоны калибровки по спокойному Солнцу
        ind_c1 = [_s[4] <= el <= _s[5] for el in [i for i in range(_len_time)]]
        # Калибровочный массив
        _spc = np.array(_sp[ind_c1, :], dtype='float32')
        # Среднее значение в попугаях спокойного Солнца от частоты
        _quiet_sun = np.array([np.nanmean(s[s > 100]) for s in _spc.T], dtype='float32')
        # Коэффициенты пересчета от попугаев к единицам потока спокойного Солнца
        _flux_coeff = np.array([_sfu[i] / _quiet_sun[i] for i in range(_len_freq)])
        # Пересчет
        _sp *= _flux_coeff
        if _l == 0:
            _head['flux_coeff_left'] = _flux_coeff
        else:
            _head['flux_coeff_right'] = _flux_coeff
        # Исключение из рассмотрения зон режекции и частот подавления куртозисом
        _sp[:, _del_zone] = 0.01
        _sp[_sp < 0.005] = 0.005

        _data[_l] = _sp
        _l += 1
    # if not ('flux_coeff_left' in _head) or not ('flux_coeff_right' in _head):
    with open(_path_to_head, 'wb') as out:
        pickle.dump(_head, out)
    return _data


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
    object = 'sun'
    current_primary_file = '2024-02-19_06+04'
    current_primary_dir = current_primary_file[0:4] + '_' + current_primary_file[5:7] + '_' + \
                          current_primary_file[8:10] + object
    main_dir = current_primary_file[0:4]  # Каталог всех данных (первичных, вторичных) за год
    # main_dir = r'2021/Results'           # Каталог (за определенный период, здесь - за 2021 год)
    date = current_primary_dir[0:10]
    adr1 = DataPaths(current_primary_file, current_primary_dir, main_dir)
    converted_data_file_path = adr1.converted_data_file_path
    data_treatment_file_path = adr1.treatment_data_file_path
    path_to_csv = Path(adr1.converted_dir_path, 'dict_calibr.csv')

    # !!!! ******************************************* !!!!
    # ****** Блок исходных параметров для обработки *******

    freq_res = 4  # Установка разрешения по частоте в МГц
    kt = 4  # Установка разрешения по времени в единицах минимального разрешения 8.3886e-3 сек
    delta_t = 8.3886e-3
    delta_f = 7.8125
    N_Nyq = 3
    freq_spect_mask = freq_mask(8)

    # att_val = [i * 0.5 for i in range(64)]
    # att_dict = {s: 10 ** (s / 10) for s in att_val}
    # *****************************************************

    band_size_init = 'whole'
    num_of_polar = 2
    # band_size = 'whole'   Параметр 'whole' означает работу в диапазоне 1-3 ГГц, 'half' - диапазон 1-2 или 2-3 ГГц
    # polar = 'both'        Принимает значения поляризаций: 'both', 'left', 'right'
    # *****************************************************
    output_picture_mode = 'y'
    save_data = 'n'  # Сохранение сканов в формате *.npy: 'y' / 'n'
    lf_filter = 'n'  # Применение НЧ фильтра для сглаживания сканов (скользящее среднее и др.): 'y' / 'n'
    low_noise_spectrum = 'n'  # Вывод графика НЧ спектра шумовой дорожки: 'y' / 'n'
    graph_3d_perm = 'n'
    contour_2d_perm = 'n'
    poly3d_perm = 'n'
    ab = 'n'  # Подготовка рисунков к публикации в АБ

    # *****************************************************
    # Чтение с диска, если спектры ранее извлекались,
    # или извлечение спектров из исходных записей
    spectr_extr_left1, spectr_extr_left2, spectr_extr_right1, spectr_extr_right2, n_aver, band_size, polar = \
        preparing_data()
    print('Data are prepared')
    aver_param = 2 ** (6 - n_aver)
    kf = int(freq_res / delta_f * aver_param)  # Установка разрешения по частоте в единицах максимального разрешения
    # для данного наблюдения delta_f/aver_param, где delta_f = 7.8125 МГц
    with open(Path(str(converted_data_file_path) + '_head.bin'), 'rb') as inp:
        head = pickle.load(inp)

    # Приведение порядка следования отсчетов по частоте к нормальному
    if np.size(spectr_extr_left1):
        N_row = np.shape(spectr_extr_left1)[0]
        for i in range(N_row):
            spectr_extr_left1[i][0:] = spectr_extr_left1[i][-1::-1]
    if np.size(spectr_extr_right1):
        N_row = np.shape(spectr_extr_right1)[0]
        for i in range(N_row):
            spectr_extr_right1[i][0:] = spectr_extr_right1[i][-1::-1]

    spectrum = pd.Series([spectr_extr_left1, spectr_extr_left2, spectr_extr_right1,
                          spectr_extr_right2])
    united_spectrum = unite_spectrum(spectrum)

    # Quiet sun calibration
    united_spectrum = sun_calibration(united_spectrum, current_primary_file, adr1)

    ser_ind = united_spectrum.index
    if len(ser_ind) == 2:
        spectrum_extr = united_spectrum[0] + united_spectrum[1]
    else:
        spectrum_extr = united_spectrum[0]

    # # ***********************************************
    # # ***        Графический вывод данных        ****
    # # ***********************************************

    # Динамическая маска (зависит от длины записи во времени)
    t_spect = N_row * delta_t
    time_spect_mask = [(lambda i: (t_spect * (i + 0.05)) // 7)(i) for i in range(7)]
    # time_spect_mask = [185.85, 186.1, 185.6, 185.35]  # az+20

    # Формирование спектров и сканов по маскам freq_spect_mask и time_spect_mask
    shift = head['shift']
    spectr_freq, spectr_time = form_spectr_sp1(spectrum_extr, freq_spect_mask, time_spect_mask)
    line_legend_time, line_legend_freq = line_legend(freq_spect_mask[:10])
    # for i in range(4):
    #     spectr_freq[i, :] = spectr_freq[2 * i + 1, :] - spectr_freq[2 * i, :]
    #     line_legend_time[i] = line_legend_time[2 * i + 1]
    # spectr_freq[0, :] = (spectr_freq[0, :] + spectr_freq[2, :]) / 2
    # spectr_freq[1, :] = spectr_freq[1, :] - spectr_freq[0, :]
    # spectr_freq[0, :] = spectr_freq[0, :] / 1000
    # line_legend_time = line_legend_time[0:2]
    # spectr_freq = spectr_freq[0:2, :]
    n_freq = len(time_spect_mask)
    n_time = len(freq_spect_mask)
    for i in range(n_freq):
        try:
            spectr_freq[i, :] = del_random_mod(spectr_freq[i, :], 0)
        except IndexError:
            pass
    for i in range(n_time):
        spectr_time[i, :] = del_random_mod(spectr_time[i, :], 0)
    print('spectr_freq, spectr_time are formed')
    # Формирование строк-аргументов по времени и частоте и легенды
    N_col = np.shape(spectrum_extr)[1]
    if band_size_init == 'half':
        freq = np.linspace(1000 * (N_Nyq - 1) + 3.9063 / aver_param * kf, 1000 * N_Nyq - 3.9063 / aver_param * kf,
                           N_col // kf)
    elif band_size_init == 'whole':
        freq = np.linspace(1000 + 3.9063 / aver_param * kf, 3000 - 3.9063 / aver_param * kf, N_col // kf)
    timeS = np.linspace(0, delta_t * N_row, N_row // kt)
    # timeS = [187.65, 188.15, 188.65, 189.15, 189.65]

    # ***************!! Вывод данных в текстовой форме !!*********************
    # path_txt = str(Path(converted_dir_path, current_data_file, '_scan.txt'))
    # path_npy1 = Path(str(converted_data_file_path) + '_spectrum_time.npy') # Спектры в фиксированные моменты времени
    path_npy2 = Path(str(converted_data_file_path) + '_scan_freq.npy')  # Сканы на фиксированных частотах
    path_npy_time = Path(str(converted_data_file_path) + '_time_count.npy')
    # print(path_txt)
    # np.savetxt(path_txt, )
    # np.save(path_npy2, spectr_time)                                        # Сканы на фиксированных частотах
    # np.save(path_npy_time, timeS)
    # np.save(path_npy1, spectr_freq)                                        # Спектры в фиксированные моменты времени
    # path_txt = str(Path(converted_dir_path, current_data_file, 'freq.txt'))
    # print(path_txt)
    # np.savetxt(path_txt, freq)
    # ***********************************************************************

    # ************************** !!! Вывод данных !!! ***********************
    # line_legend_time, line_legend_freq = line_legend(freq_spect_mask[:10])
    if not 'good_bound' in head:
        head['good_bound'] = 0.1
    a = head['good_bound']
    info_txt = [f'time resol = {delta_t * kt: .4f} sec',
                f'freq resol = {delta_f / aver_param * kf: 1.2f} MHz',
                f'polarisation {polar}', 'align: quiet Sun', f'kurtosis quality = {a}']
    path1 = data_treatment_file_path
    # ********************** Сохранение сканов в формате *.npy **************
    if save_data == 'y':
        np.save(path1, spectr_time)
        path_mesh = Path(data_treatment_file_path, 'Meshed_Spectrum', current_primary_file + '_meshed')
        path_mesh1 = Path(data_treatment_file_path, 'Meshed_Spectrum', current_primary_file + '_meshed' + '.npy')
        if os.path.isfile(path_mesh1):
            print('Meshed spectrum file exist')
        else:
            path_to_fig(Path(data_treatment_file_path, 'Meshed_Spectrum'))
            spectrum_mesh = spectr_construction(spectrum_extr, kf, kt)
            np.save(path_mesh, spectrum_mesh)
    # ***********************************************************************
    if lf_filter == 'y':
        spectr_time = signal_filtering(spectr_time, 1.0)

    # ***********************************************************************
    #               ****** Low noise spectra ******
    if low_noise_spectrum == 'y':
        spectrum_signal_av = low_freq_noise_spectrum(spectr_time, 32768 // 4)
        if kt == 1 & kf == 1:
            m, n = spectrum_signal_av.shape
            f_max = 1 / delta_t / 2
            f_min = f_max / n
            arg = np.linspace(f_min, f_max, n)
            low_noise_spectra_base(spectrum_signal_av, head, freq_spect_mask, arg, current_primary_file)
        # np.save(Path(path1, 'LN_spectrum'), spectrum_signal_av)
        if ab == 'y':
            plot_low_freq_spec_ab(spectrum_signal_av, delta_t * kt, path1, line_legend_freq)
        else:
            plot_low_freq_spec(spectrum_signal_av, delta_t * kt, path1, line_legend_freq)
    #                       *****************************

    if output_picture_mode == 'y':
        if ab == 'y':
            # fp.fig_plot_ab(spectr_freq, 0, freq, 1, info_txt, path1, head, line_legend_time)
            fp.fig_plot_ab(spectr_time, 0, timeS, 0, info_txt, path1, head, line_legend_freq)
        else:
            fp.fig_plot(spectr_freq, 0, freq, 1, info_txt, path1, head, line_legend_time)
            fp.fig_plot(spectr_time, 0, timeS, 0, info_txt, path1, head, line_legend_freq)
    # *********************************************************
    # ***            Многооконный вывод данных             ****
    # *********************************************************
    if output_picture_mode == 'no':
        t_start, t_stop = 50, 180
        n_start, n_stop = int(t_start / delta_t / kt), int(t_stop / delta_t / kt)
        if ab == 'y':
            fp.fig_multi_axes_ab(spectr_time[:10, n_start:n_stop], timeS[n_start:n_stop], info_txt,
                                 path1, freq_spect_mask, head)
        else:
            fp.fig_multi_axes(spectr_time[:10, n_start:n_stop], timeS[n_start:n_stop], info_txt,
                              path1, freq_spect_mask, head)

    # *********************************************************
    # ***        Вывод данных двумерный и трехмерный       ****
    # *********************************************************
    # Укрупнение  разрешения по частоте и времени для вывода в 2d и 3d
    if graph_3d_perm == 'y' or contour_2d_perm == 'y' or poly3d_perm == 'y':
        spectr_extr1 = spectr_construction(spectrum_extr, kf, kt)
    # Информация о временном и частотном резрешениях
    info_txt = [('time resol = ' + str(delta_t * kt) + 'sec'),
                ('freq resol = ' + str(delta_f / aver_param * kf) + 'MHz'),
                ('polarisation ' + polar)]

    if graph_3d_perm == 'y':
        t_start, t_stop = 175, 195
        n_start, n_stop = int(t_start / delta_t / kt), int(t_stop / delta_t / kt)
        f_start = 1000
        f_stop = 1400
        nf_start = int((f_start - 1000) / freq_res)
        nf_stop = int((f_stop - 1000) / freq_res)
        fp.graph_3d(freq[nf_start:nf_stop], timeS[n_start:n_stop], spectr_extr1[n_start:n_stop, nf_start:nf_stop],
                    3, path1, head)
    if contour_2d_perm == 'y':
        if ab == 'y':
            fp.graph_contour_2d_ab(freq, timeS, spectr_extr1, 0, info_txt, path1, head)
        else:
            fp.graph_contour_2d(freq, timeS, spectr_extr1, 0, info_txt, path1, head)

    if poly3d_perm == 'y':
        data_poly3d, freq_mask = data_poly3d_prep(spectr_extr1)
        poly_graph3d(timeS, data_poly3d, freq_mask)

    stop = datetime.now()
    print('\n Total time = ', stop - start)
