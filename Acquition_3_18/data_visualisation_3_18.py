import numpy as np
import os
import sys
import pandas as pd
import pickle
import json as jsn
from datetime import datetime
from pathlib import Path
import matplotlib.pyplot as plt
from Supporting_func import Fig_plot as fp, align_spectrum, path_to_data
# from Supporting_func import align_spectrum, path_to_data
from Polyphase import low_freq_noise_spectrum, plot_low_freq_spec, plot_low_freq_spec_ab
from Interface.window_handler import exec_app
from Polyphase.cic_filter import signal_filtering
# from test_statistic import low_noise_spectra_base
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


class DataPaths(object):

    def __init__(self, _data_file, _data_dir, _main_dir):
        if _data_dir.find('test') != -1 or _data_dir.find('calibration') != -1 or _data_dir.find('calibr') != -1:
            _main_dir = '2023'
        self.data_file_name = _data_file
        self.data_file_prime = _data_file + '.bin'
        self.data_dir = _data_dir
        self.main_dir = _main_dir
        self.head_path = path_to_data()
        self.primary_dir_path, self.primary_data_file_path = self.primary_paths()
        self.converted_dir_path, self.converted_data_file_path = self.converted_paths()
        self.treatment_dir_path, self.treatment_data_file_path = self.treat_paths()

    def primary_paths(self):
        _path = Path(self.head_path, self.main_dir, 'Primary_data_3_18', self.data_dir)
        create_folder(_path)
        if self.__check_paths():
            _primary_data_path = Path(_path, self.data_file_prime)
        else:
            raise CustomError('Head path not found!')
        return _path, _primary_data_path

    def converted_paths(self):
        _path = Path(self.head_path, self.main_dir, 'Converted_data_3_18', str(self.data_dir) + '_conv')
        create_folder(_path)
        if self.__check_paths():
            _convert_data_path = Path(_path, self.data_file_name)
        else:
            raise CustomError('Head path not found!')
        return _path, _convert_data_path

    def treat_paths(self):
        _path = Path(self.head_path, self.main_dir, 'Data_treatment_3_18', str(self.data_dir) + '_treat')
        create_folder(_path)
        if self.__check_paths():
            _treatment_data_path = Path(_path, self.data_file_name)
        else:
            raise CustomError('Head path not found!')
        return _path, _treatment_data_path

    def __check_paths(self):
        return not self.head_path == 'Err'


class CustomError(Exception):
    pass


def path_to_data():
    """
    Определяет путь на конкретной машине к корню каталога данных.
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
        return 'Err'
    return head_path_out


def create_folder(_path):
    if not os.path.isdir(_path):
        os.mkdir(_path)


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


def band_map():
    bm = [
        [1000, 2900],
        [2900, 5200],
        [5200, 6500],
        [6500, 8800],
        [8800, 11200],
        [11200, 13500],
        [13500, 15800],
        [15800, 18100]
    ]
    return bm


def form_spectrum(_sp, _freq_spect_mask_in=freq_spect_mask, _time_mask=time_spect_mask):
    """ Возвращает s_freq - срезы частотного спектра в моменты времени time_spect_mask и s_time - сканы Солнца
    по времени на частотах freq_spect_mask с заданным разрешением по времени и частоте

    """
    # Определение заполненных поддиапазонов
    _shape_sp = np.shape(_sp)
    _map_sp = []    # Номера поддиапазонов, в которых проводилось наблюдение

    for _i in range(_shape_sp[0]):
        for _j in range(_shape_sp[1]):
            _flag = 0
            if np.size(_sp[_i, _j]) > 2:
                _flag = 1
            if _flag and _i not in _map_sp:
                _map_sp.append(_i)

    # Определение наличия данных для сканов на частоах из freq_spect_mask
    _freq_map = {}                  # Частота анализа - ключ, номер поддиапазона - значение
    _band_map = band_map()
    for _i in _map_sp:
        for _f in _freq_spect_mask_in:
            if _band_map[_i][0] <= _f <= _band_map[_i][1]:
                _freq_map[_f] = _i
    _freq_mask = list(_freq_map.keys())   # Частоты, для которых доступны сканы во времени

    # Словарь, в котором ключ - номер поддиапазона, значение - количество частотных отсчетов в нем
    _N_col = {s: np.shape(_sp[s, 1])[1] for s in _map_sp}
    _N_raw = {s: np.shape(_sp[s, 1])[0] for s in _map_sp}
    _kf = int(freq_res / delta_f / aver_param)
    if not _kf:
        sys.exit(f'Frequency resolution freq_res={freq_res}MHz is too large. Decrease it, '
                 f'increase freq_res to {delta_f * aver_param}MHz at least')

    s_freq = [[[] for i in range(4)] for j in range(len(_freq_mask) * 2)]
    s_time = [[[] for i in range(4)] for j in range(len(_time_mask) * 2)]
    j = 0
    for _f in _freq_mask:

        ind = int((_f - _band_map[int(_freq_map[_f])][0] - delta_f * aver_param / 2) // (delta_f * aver_param))
        # - количество отсчетов по частоте в поддиапазоне с разрешением не лучше freq_res
        if ind > _N_col[int(_freq_map[_f])] - int(kf / 2) - 1:
            ind = _N_col [int(_freq_map[_f])]- int(kf / 2) - 1
        if ind < int(kf / 2):
            ind = int(kf / 2)
        i = 0
        for _i in range(2):
            _s_loc1 = []
            _time_count = []
            sp = _sp[_freq_map[_f], _i]
            if len(sp) < 2:
                s_time[j][2] = []
                s_time[j][0] = _f
                s_time[j][1] = _i
                j += 1
            else:
                while kt * (i + 1) < _N_raw[int(_freq_map[_f])]:
                    if kf == 1:
                        _s_loc = np.sum(sp[i * kt:(i + 1) * kt, ind][_sp[i * kt:(i + 1) * kt, ind] > 40])
                        n_mesh = (sp[i * kt: (i + 1) * kt, ind] > 40).sum()

                        if n_mesh == 0:
                            _s_loc = 2
                        else:
                            _s_loc /= n_mesh

                    else:
                        _s_loc = np.sum(sp[i * kt:(i + 1) * kt, ind - int(kf / 2):ind + int(kf / 2)][
                                                  sp[i * kt:(i + 1) * kt, ind - int(kf / 2):ind + int(kf / 2)] > 40])
                        n_mesh = (sp[i * kt: (i + 1) * kt, ind - int(kf / 2):ind + int(kf / 2)] > 40).sum()
                        if n_mesh == 0:
                            _s_loc = 2
                        else:
                            _s_loc /= n_mesh
                    i += 1
                    _s_loc1.append(_s_loc)
                    _time_count.append(kt * (i + 0.5) * delta_t)
                s_time[j][2] = np.array(_s_loc1)
                s_time[j][3] = np.array(_time_count)
                s_time[j][0] = _f
                s_time[j][1] = _i
                j += 1
    s_time = np.array(s_time)
    i = 0
    for t in _time_mask:
        ind = int(t // delta_t)
        if ind > N_row - kt / 2 - 1:
            ind = N_row - int(kt / 2) - 1
        if ind < (kt / 2):
            ind = int(kt / 2)
        j = 0
        for _i in range(2):
            _s_loc1 = []
            _freq_count = []
            sp = _sp[_freq_map[_f], _i]
            if len(sp) < 2:
                s_freq[j][2] = []
                s_freq[j][0] = t
                s_freq[j][1] = _i
                j += 1
            else:
                while (j + 1) * kf < N_col:
                    if kt == 1:
                        s_loc = np.sum(sp[ind, j * kf:(j + 1) * kf][sp[ind, j * kf:(j + 1) * kf] > 40])
                        n_mesh = (sp[ind, j * kf:(j + 1) * kf] > 40).sum()
                        if n_mesh == 0:
                            s_loc = 2
                        else:
                            s_loc /= n_mesh
                    else:
                        s_loc = np.sum(sp[ind - int(kt / 2):ind + int(kt / 2), j * kf :(j + 1) * kf][
                                                  sp[ind - int(kt / 2):ind + int(kt / 2), j * kf:(j + 1) * kf] > 40])
                        n_mesh = (sp[ind - int(kt / 2):ind + int(kt / 2), j * kf:(j + 1) * kf] > 40).sum()
                        if n_mesh == 0:
                            s_loc = 2
                        else:
                            s_loc /= n_mesh
                        _s_loc1.append(_s_loc)
                        _freq_count.append(kf * (i + 0.5) * delta_f)
                s_freq[j][2] = np.array(_s_loc1)
                s_freq[j][3] = np.array(_time_count)
                s_freq[j][0] = t
                s_freq[j][1] = _i
                j += 1
        ind_time.append(ind)
        i += 1
    s_freq = np.array(s_freq)

    return s_freq, s_time


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
                N_mesh = (Spectr[i * kt: (i + 1) * kt, j * kf: (j + 1) * kf] > 10).sum()
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

    # Для k-ой полосы и двух возможных поляризаций выдает по спектру для каждой поляризации.
    # Если поляризация не задействована, то соответствующие спектры - пустые. Спектры k = 3,4 - в обратном порядке
    _path1 = Path(str(converted_data_file_path) + '_spectrum.npy')
    _spectrum = np.load(_path1, allow_pickle=True)
    with open(Path(str(converted_data_file_path) + '_head.bin'), 'rb') as inp:
        _head = pickle.load(inp)
    _n_aver = _head['n_aver']
    _polar = _head['polar']

    # Разделяем составляющие  записи в полной полосе и с возможными двумя поляризациями

    # Выдает спектры для левой и правой поляризаций шириной по 1 ГГц. При нумерации спектров учитывается
    # значение зоны Найквиста. С индексом "1" - 1-2 ГГц, с индексом "2" - 2-3 ГГц, как и для случая выше.
    # На выходе формально все 4 спектра, но для незадействованной полосы они пустые

    return _spectrum, int(_n_aver), _polar


def band_separation(_s, _n_aver):
    """
    Функция принимает запись спектров для всей зоны Найквиста, отбрасывает спектры из защитной
    зоны частот, формируя матрицу сканов на частотах поддиапазона, например, 2.9 - 5.2 ГГц, шириной
    2.3 ГГц при том, что ширина зоны Найквиста - 2.7 ГГц
    :param _s: Матрица спектров по поддиапазонам
    :param _n_aver: количество усредняемых по частоте бинов для отсчета спектра
    :return: Матрицу спектров для поддиапазона без защитных зон
    """
    #
    _n_start = int(2432 / _n_aver)
    _n_stop = int(32768 / _n_aver) - _n_start
    _shape_s = np.shape(_s)
    for i in range(_shape_s[0]):
        for j in range(_shape_s[1]):
            if np.size(_s[i, j]):
                _N_row = np.shape(spectrum[i, j])[0]
                s1 = np.zeros((_N_row, _n_stop - _n_start))
                s = _s[i, j]
                for k in range(_N_row):
                    s1[k, :] = s[k, _n_start:_n_stop]
                _s[i, j] = s1

    return _s


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
        spec_left0 = np.full((n_row1, shape_l[1]), 2)
        spec_right0 = np.full((n_row1, shape_r[1]), 2)
        spec_left0[:shape_l[0]] = spec_left
        spec_right0[:shape_r[0]] = spec_right
        united_spec = pd.Series([spec_left0, spec_right0], ind)
    elif np.size(spec_left):
        united_spec = pd.Series([spec_left], ind)
    else:
        united_spec = pd.Series([spec_right], ind)
    print('Spectra are united')
    return united_spec


def ngi_temperature(_f, _case_id, _path):
    with open(_path, 'rb') as _inp:
        _data_saved = pickle.load(_inp)
    _arr_av_left = _data_saved['polyval_coeff'][data['case_id'] == _case_id][data['polar'] == 'left'].iloc[0]
    _arr_av_right = _data_saved['polyval_coeff'][data['case_id'] == _case_id][data['polar'] == 'right'].iloc[0]
    _temp_interp_left = np.polyval(_arr_av_left, _f)
    _temp_interp_right = np.polyval(_arr_av_right, _f)

    return _temp_interp_left, _temp_interp_right


def noise_self_calibration(_spectrum, _ngi_temperature_path):
    # Временные интервалы для калибровки по внутреннему ГШ
    t_cal0, t_cal1 = 0, 5.5  # Интервал Импульса ГШ, сек
    t_ground1, t_ground2 = 7, 12  # Интервал определения фона, сек
    t_cal0r, t_cal1r = 0, 5.5  # Интервал Импульса ГШ, сек
    t_ground1r, t_ground2r = 7, 12

    # Закрузка шумовой калибровочной температуры на входе приемника
    with open(_ngi_temperature_path, 'rb') as _inp:
        _ngi_temperature = pickle.load(_inp)

    temp_left = _ngi_temperature[_ngi_temperature['ngi_id'] == '01'][_ngi_temperature['polar'] == 'left']
    temp_right = _ngi_temperature[_ngi_temperature['ngi_id'] == '01'][_ngi_temperature['polar'] == 'right']

    temp_left0 = temp_left['low_band'][1]
    temp_left1 = temp_left['upper_band'][1]
    temp_right0 = temp_right['low_band'][0]
    temp_right1 = temp_right['upper_band'][0]
    _l = np.size(temp_left0)

    # Интервалы в отсчетах для калибровки по внутреннему ГШ
    n_cal0, n_cal1 = int(t_cal0 // delta_t), int(t_cal1 // delta_t)
    n_ground1, n_ground2 = int(t_ground1 // delta_t), int(t_ground2 // delta_t)

    # Пересчет данных из относительных единиц в температуру для левой поляризации
    if np.size(_spectrum[0]):
        s_left0 = _spectrum[0][n_cal0:n_cal1, :]
        s_left1 = _spectrum[1][n_cal0:n_cal1, :]
        s_ground0 = _spectrum[0][n_ground1:n_ground2, :]
        s_ground1 = _spectrum[1][n_ground1:n_ground2, :]

        _m, _n = np.shape(_spectrum[0])
        s_left0_av = np.array([np.mean(s_left0[:, i][s_left0[:, i] > 100], dtype=np.int64) for i in range(_n)])
        s_left0_av = treatment_null_mesh(s_left0_av, 10)
        s_left1_av = np.array([np.mean(s_left1[:, i][s_left1[:, i] > 100], dtype=np.int64) for i in range(_n)])
        s_left1_av = treatment_null_mesh(s_left1_av, 10)
        s_ground_av0 = np.array([np.mean(s_ground0[:, i][s_ground0[:, i] > 100], dtype=np.int64) for i in range(_n)])
        s_ground_av0 = treatment_null_mesh(s_ground_av0, 10)
        s_ground_av1 = np.array([np.mean(s_ground1[:, i][s_ground1[:, i] > 100], dtype=np.int64) for i in range(_n)])
        s_ground_av1 = treatment_null_mesh(s_ground_av1, 10)
        s_left0_av = s_left0_av - s_ground_av0
        s_left1_av = s_left1_av - s_ground_av1

        #                   **** Приведение разрешения по частоте ****
        #           калибровочной шумовой температуры к разрешению по частоте
        #                               преобразуемых данных
        _k = _l // _n
        temp_left0m = np.reshape(temp_left0, (_k, -1), order='F')
        temp_left1m = np.reshape(temp_left1, (_k, -1), order='F')
        temp_left0m_av = np.mean(temp_left0m, axis=0)
        temp_left1m_av = np.mean(temp_left1m, axis=0)
        coeff_left0 = temp_left0m_av / s_left0_av
        coeff_left1 = temp_left1m_av / s_left1_av

        # fig, ax = plt.subplots(1, figsize=(12, 6))
        #
        # plt.grid()
        # plt.show()

        for i in range(_m):
            _spectrum[0][i, :] = _spectrum[0][i, :] * coeff_left0

            _spectrum[1][i, :] = _spectrum[1][i, :] * coeff_left1

    # Пересчет данных из относительных единиц в температуру для правой поляризации
    if np.size(_spectrum[2]):
        n_cal0r, n_cal1r = int(t_cal0r // delta_t), int(t_cal1r // delta_t)
        n_ground1r, n_ground2r = int(t_ground1r // delta_t), int(t_ground2r // delta_t)
        s_right0 = _spectrum[2][n_cal0r:n_cal1r, :]
        s_right1 = _spectrum[3][n_cal0r:n_cal1r, :]
        s_ground2 = _spectrum[2][n_ground1r:n_ground2r, :]
        s_ground3 = _spectrum[3][n_ground1r:n_ground2r, :]
        _m1, _n1 = np.shape(_spectrum[2])

        s_right0_av = np.array([np.mean(s_right0[:, i][s_right0[:, i] > 100], dtype=np.int64) for i in range(_n1)])
        s_right0_av = treatment_null_mesh(s_right0_av, 10)
        s_right1_av = np.array([np.mean(s_right1[:, i][s_right1[:, i] > 100], dtype=np.int64) for i in range(_n1)])
        s_right1_av = treatment_null_mesh(s_right1_av, 10)
        s_ground_av2 = np.array([np.mean(s_ground2[:, i][s_ground2[:, i] > 100], dtype=np.int64) for i in range(_n1)])
        s_ground_av2 = treatment_null_mesh(s_ground_av2, 10)
        s_ground_av3 = np.array([np.mean(s_ground3[:, i][s_ground3[:, i] > 100], dtype=np.int64) for i in range(_n1)])
        s_ground_av3 = treatment_null_mesh(s_ground_av3, 10)
        s_right0_av = s_right0_av - s_ground_av2
        s_right1_av = s_right1_av - s_ground_av3

        _k1 = _l // _n1
        temp_right0m = np.reshape(temp_right0, (_k1, -1), order='F')
        temp_right1m = np.reshape(temp_right1, (_k1, -1), order='F')
        temp_right0m_av = np.mean(temp_right0m, axis=0)
        temp_right1m_av = np.mean(temp_right1m, axis=0)
        coeff_right0 = temp_right0m_av / s_right0_av
        coeff_right1 = temp_right1m_av / s_right1_av

        # fig, ax = plt.subplots(1, figsize=(12, 6))
        # ax.plot(coeff_left0)
        # ax.plot(coeff_left1)
        # ax.plot(coeff_right0)
        # ax.plot(coeff_right1)
        # plt.grid()
        # plt.show()

        for i in range(_m1):
            _spectrum[2][i, :] = _spectrum[2][i, :] * coeff_right0
            _spectrum[2][i, :][_spectrum[2][i, :] < 0] = 0
            _spectrum[3][i, :] = _spectrum[3][i, :] * coeff_right1
            _spectrum[3][i, :][_spectrum[3][i, :] < 0] = 0
        a = _spectrum[2]
        b = _spectrum[3]

    ant_coeff = np.loadtxt(ant_coeff_path)
    len_ant_coeff = np.size(ant_coeff)
    len1 = int(len_ant_coeff / 2)
    if np.size(spectrum[2]):
        k = _n1 / len1
    else:
        k = _n / len1

    # ant_coeff_low = ant_coeff[:len1]
    # ant1 = np.ravel(np.transpose(np.array([ant_coeff_low]*int(k))))
    # a = _spectrum[0]
    # for i in range(int(_m)):
    #     _spectrum[0][i, :] = _spectrum[0][i, :] / ant1
    # for i in range(int(_m1)):
    #     _spectrum[2][i, :] = _spectrum[2][i, :] / ant1
    # a1 = _spectrum[0]
    pass
    # fig, ax = plt.subplots(1, figsize=(12, 6))
    # ax.plot(a[50, :])
    # ax.plot(a1[50, :])
    # plt.grid()
    # plt.show()
    return _spectrum


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
    for _i in range(1, _l - 2):
        if abs(_s[_i] - _s[_i - 1]) > _ds * abs(_s[_i + 1] - _s[_i - 1]):
            _s[_i] = (_s[_i - 1] + _s[_i + 1]) / 2
        if _s[_i] <= _s0:
            _s[_i] = _s0
    _s[0] = _s0
    _s[_l - 2:] = _s0
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


def noise_black_body_calibration(_spectrum, _receiver_temperature_path):
    """
    Функция пересчитывает входные данные в градусы Кельвина по имеющейся в записи калибровке черным телом в
    промежуток времени t_cal0 ... t_cal1
    qweertr33
    :param _spectrum: входные данные
    :param _receiver_temperature_path: путь к файлу со значениями температуры собственных шумов приемника
    :return _spectrum: Входные данные, пересчитанные в градусы К по привязке к черному телу с вычтенной шумовой
    температурой приемника
    """

    n_av_loc = 10  # Половина длины последовательности отсчетов для усреднения функцией treatment_null_mesh(,_s0,)
    _s_threshold = 1e6  # Порог для значений исходного массива _spectrum для функции del_random_mod(,_n,)
    #                   *** Загрузка температуры собственных шумов приемника ***
    #              **** Через интерполяцию средней температуры по серии измерений ****
    #
    _s = _spectrum.iloc[0][0, :]
    if len(_s):
        _l1 = len(_s)
    else:
        _l1 = len(_spectrum.iloc[0][2, :])
    _l = int(2 * _l1)
    _receiver_temp = receiver_noise_temperature(receiver_polyval_coeff_path, _l, '03')
    temp_left0 = _receiver_temp[0: _l1]
    temp_left1 = _receiver_temp[_l1: _l]
    temp_right0 = _receiver_temp[0: _l1]
    temp_right1 = _receiver_temp[_l1: _l]
    #           *** Загрузка температуры собственныхшумов непосредственно ***
    #                  *** из результатов измерений ***
    # with open(_receiver_temperature_path, 'rb') as _inp:
    #     _receiver_temperature = pickle.load(_inp)
    # temp_left = _receiver_temperature['temperature'][_receiver_temperature['date'] ==
    #                                                  '2022-11-18'][_receiver_temperature['polar'] == 'left'].iloc[0]
    # temp_right = _receiver_temperature['temperature'][_receiver_temperature['date'] ==
    #                                                  '2022-11-18'][_receiver_temperature['polar'] == 'right'].iloc[0]
    # _l = np.size(temp_left)
    # temp_left0 = temp_left[0: _l1]
    # temp_left1 = temp_left[_l1: _l]
    # temp_right0 = temp_right[0: _l1]
    # temp_right1 = temp_right[_l1: _l]

    # Интервалы в отсчетах для калибровки по черному телу
    n_cal0, n_cal1 = int(t_cal0 // delta_t), int(t_cal1 // delta_t)

    # Пересчет данных из относительных единиц в температуру для левой поляризации
    if np.size(spectrum[0]):
        s_left0 = spectrum[0][n_cal0:n_cal1, :]
        s_left1 = spectrum[1][n_cal0:n_cal1, :]
        _m, _n = np.shape(_spectrum[0])
        s_left0_av = np.array([np.mean(s_left0[:, i][s_left0[:, i] > 100], dtype=np.int64) for i in range(_n)])
        # Исключение неоднозначности, если [np.mean(s_left0[:, i][s_left0[:, i] > 100] пуст
        s_left0_av = treatment_null_mesh(s_left0_av, n_av_loc)
        # Убираем одиночные выбросы отсчетов и приравниваем отсчеты ниже _s_threshold к _s_threshold
        s_left0_av = del_random_mod(s_left0_av, _s_threshold)
        s_left01_av = s_left0_av.copy()
        s_left0_av = del_random_mod(s_left0_av, _s_threshold)
        s_left1_av = np.array([np.mean(s_left1[:, i][s_left1[:, i] > 100], dtype=np.int64) for i in range(_n)])
        # Исключение неоднозначности, если [np.mean(s_left1[:, i][s_left1[:, i] > 100] пуст
        s_left1_av = treatment_null_mesh(s_left1_av, n_av_loc)
        s_left1_av = del_random_mod(s_left1_av, 1e6)

        #                   **** Приведение разрешения по частоте ****
        #           температуры собственных шумов приемника к разрешению по частоте
        #                               преобразуемых данных
        _k = _l1 // _n
        temp_left0m = np.reshape(temp_left0, (_k, -1), order='F')
        temp_left1m = np.reshape(temp_left1, (_k, -1), order='F')
        temp_left0m_av = np.mean(temp_left0m, axis=0)
        temp_left1m_av = np.mean(temp_left1m, axis=0)
        #           **** Определение коэффициентов пересчета произвольных единиц в градусы К ****
        coeff_left0 = (temp_left0m_av + 300) / s_left0_av
        coeff_left0 = del_random_mod(coeff_left0, 1e-6)
        coeff_left1 = (temp_left1m_av + 300) / s_left1_av
        coeff_left1 = del_random_mod(coeff_left1, 1e-6)

        # fig, ax = plt.subplots(1, figsize=(12, 6))
        # ax.plot(s_left0_av)
        # ax.plot(s_left01_av)
        # plt.grid()
        # plt.show()

        #               **** Пересчет в температуру с вычитанием собственных шумов приемника ****
        #     **** Температура собственных шумов приемника - табличная, ранее рассчитанная по измерениям ****
        for i in range(_m):
            a = _spectrum[0][i, :]
            b = _spectrum[1][i, :]
            #   Приводим вид сканов к "температурному" виду и Вычитаем температуру собственных шумов приемника
            _spectrum[0][i, :] = a * coeff_left0 - temp_left0m_av
            _spectrum[0][i, :][_spectrum[0][i, :] < 0] = 0
            _spectrum[1][i, :] = b * coeff_left1 - temp_left1m_av
            _spectrum[1][i, :][_spectrum[1][i, :] < 0] = 0
            pass

    # Пересчет данных из относительных единиц в температуру для правой поляризации
    if np.size(spectrum[2]):
        s_right0 = spectrum[2][n_cal0:n_cal1, :]
        s_right1 = spectrum[3][n_cal0:n_cal1, :]
        _m1, _n1 = np.shape(_spectrum[2])

        s_right0_av = np.array([np.mean(s_right0[:, i][s_right0[:, i] > 100], dtype=np.int64) for i in range(_n1)])
        s_right0_av = treatment_null_mesh(s_right0_av, n_av_loc)
        s_right0_av = del_random_mod(s_right0_av, _s_threshold)
        s_right1_av = np.array([np.mean(s_right1[:, i][s_right1[:, i] > 100], dtype=np.int64) for i in range(_n1)])
        s_right1_av = treatment_null_mesh(s_right1_av, n_av_loc)
        s_right1_av = del_random_mod(s_right1_av, _s_threshold)

        #                   **** Приведение разрешения по частоте ****
        #           температуры собственных шумов приемника к разрешению по частоте
        #                               преобразуемых данных
        _k1 = _l1 // _n1
        temp_right0m = np.reshape(temp_right0, (_k1, -1), order='F')
        temp_right1m = np.reshape(temp_right1, (_k1, -1), order='F')
        temp_right0m_av = np.mean(temp_right0m, axis=0)
        temp_right1m_av = np.mean(temp_right1m, axis=0)
        #           **** Определение коэффициентов пересчета произвольных единиц в градусы К ****
        coeff_right0 = (temp_right0m_av + 300) / s_right0_av
        coeff_right0 = del_random_mod(coeff_right0, 1e-6)
        coeff_right1 = (temp_right1m_av + 300) / s_right1_av
        coeff_right1 = del_random_mod(coeff_right1, 1e-6)

        #               **** Пересчет в температуру с вычитанием собственных шумов приемника ****
        #     **** Температура собственных шумов приемника - табличная, ранее рассчитанная по измерениям ****
        for i in range(_m1):
            a = _spectrum[2][i, :]
            b = _spectrum[3][i, :]
            _spectrum[2][i, :] = a * coeff_right0 - temp_right0m_av
            _spectrum[2][i, :][_spectrum[2][i, :] < 0] = 0
            _spectrum[3][i, :] = b * coeff_right1 - temp_right1m_av
            _spectrum[3][i, :][_spectrum[3][i, :] < 0] = 0

    pass
    # fig, ax = plt.subplots(1, figsize=(12, 6))
    # ax.plot(coeff_left0)
    # ax.plot(coeff_left1)
    # ax.plot(coeff_right0)
    # ax.plot(coeff_right1)
    # plt.grid()
    # plt.show()
    return _spectrum


def receiver_noise_temperature(_path, _n, _case):
    """
    Функция принимает путь к коэффициентам полинома аппроксимации температуры собственных шумов радиометра,
    длину вектора-аргумента частот и идентификатор серии измерений, по которой аппроксимируется температура.
    :param _path:   путь к коэффициентам полинома аппроксимации
    :param _n:      длина вектора-аргумента частот
    :param _case:   идентификатор серии измерений
    :return: расчетные значения температуры собственных шумов для вектора аргументов (по частоте)
    """
    with open(_path, 'rb') as _inp:
        _interpol_coefficients = pickle.load(_inp)
    _polinom_array = _interpol_coefficients.polyval_coeff[_interpol_coefficients.case_id == _case].iloc[0]
    _df = 2000 / _n
    _freq = [1000 + _df / 2 + _df * i for i in range(_n)]
    _receiver_temp = np.polyval(_polinom_array, _freq)
    return _receiver_temp


def freq_mask(_i):
    _n1 = 2
    _n2 = 4
    _freq_mask = [
        [3000, 3600, 4500, 5100, 6600, 7500, 8200, 8600],  # [0]

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


if __name__ == '__main__':

    align_file_name = 'antenna_temperature_coefficients.npy'  # Имя файла с текущими коэффициентами выравнивания АЧХ
    current_primary_file = '2023-06-25_03'
    current_primary_dir = '2023_06_25test'
    main_dir = '2023'
    date = current_primary_dir[0:10]
    adr1 = DataPaths(current_primary_file, current_primary_dir, main_dir)
    converted_data_file_path = adr1.converted_data_file_path
    data_treatment_file_path = adr1.treatment_data_file_path
    folder_align_path = Path(adr1.head_path, 'Alignment')

    ngi_temperature_file_name = 'ngi_temperature.npy'  # Файл усредненной по базе шумовой температуры для ГШ
    receiver_temperature_file = 'receiver_temperature.npy'  #
    receiver_polyval_file = 'receiver_temp_interpol.npy'
    ant_coeff_file = 'ant_afc.txt'

    ngi_temperature_path = Path(adr1.head_path, 'Alignment', ngi_temperature_file_name)
    receiver_temperature_path = Path(adr1.head_path, 'Alignment', receiver_temperature_file)
    receiver_polyval_coeff_path = Path(adr1.head_path, 'Alignment', receiver_polyval_file)
    ant_coeff_path = Path(adr1.head_path, 'Alignment', ant_coeff_file)

    # !!!! ******************************************* !!!!
    # ****** Блок исходных параметров для обработки *******

    freq_res = 22  # Установка разрешения по частоте в МГц
    kt = 64  # Установка разрешения по времени в единицах минимального разрешения 8.3886e-3 сек
    delta_t = 12.427567e-3  # sec
    delta_f = 0.082397461   # MHz

    att_val = [i * 0.5 for i in range(64)]
    att_dict = {s: 10 ** (s / 10) for s in att_val}
    freq_spect_mask = freq_mask(0)
    # *****************************************************
    polar = 'both'
    # polar = 'both'        Принимает значения поляризаций: 'both', 'left', 'right'
    # *****************************************************
    # output_picture_mode = 'y'
    # align = 'y'  # Выравнивание АЧХ усилительного тракта по калибровке от ГШ: 'y' / 'n'
    # noise_calibr = 'n'
    # black_body_calibr = 'n'
    # save_data = 'n'  # Сохранение сканов в формате *.npy: 'y' / 'n'
    # lf_filter = 'n'  # Применение НЧ фильтра для сглаживания сканов (скользящее среднее и др.): 'y' / 'n'
    # low_noise_spectrum = 'n'  # Вывод графика НЧ спектра шумовой дорожки: 'y' / 'n'
    # graph_3d_perm = 'n'
    # contour_2d_perm = 'n'
    # poly3d_perm = 'n'
    # ab = 'n'  # Подготовка рисунков к публикации в АБ
    a = 0
    # *****************************************************
    # Чтение с диска, если спектры ранее извлекались,
    # или извлечение спектров из исходных записей
    spectrum, n_aver, polar = preparing_data()
    print('Data are prepared')
    aver_param = 2 ** n_aver
    kf = int(freq_res / delta_f / aver_param)  # Установка разрешения по частоте в единицах максимального разрешения
    # для данного наблюдения delta_f*aver_param, где delta_f = 0.082397461 МГц
    with open(Path(str(converted_data_file_path) + '_head.bin'), 'rb') as inp:
        head = pickle.load(inp)

    # Выравнивание спектров по результатам шумовых измерений АЧХ
    # if align == 'y':
    #     if head['att3'] == 0:
    #         pos = 0
    #     elif head['att3'] == 5:
    #         pos = 1
    #     elif head['att3'] == 10:
    #         pos = 2
    #     elif head['att3'] == 15:
    #         pos = 3
    #     else:
    #         pos = 0
    #     path_output = Path(folder_align_path, align_file_name)
    #     spectr_extr_left1, spectr_extr_left2, spectr_extr_right1, spectr_extr_right2 = \
    #         align_spectrum(spectr_extr_left1, spectr_extr_left2, spectr_extr_right1, spectr_extr_right2,
    #                        head, path_output, pos)
    #     print('spectrum is aligned')
    a = 0
    # Приведение порядка следования отсчетов по частоте к нормальному,
    # условно, требуют смены порядка записи 1 и 4 - по совокупности признаков
    shape_sp = np.shape(spectrum)
    for i in range(shape_sp[0]):
        for j in range(shape_sp[1]):
            if np.size(spectrum[i, j]) and i in [1, 4]:
                N_row = np.shape(spectrum[i, j])[0]
                N_col = np.shape(spectrum[i, j])[1]
                s = spectrum[i, j]
                for k in range(N_row):
                    s[k, :] = s[k, -1::-1]
                spectrum[i, j] = s
                s1 = spectrum[i, j]
                aver_param = int(32768 / N_col)

    # Отбрасывание отсчетов защитной зоны поддиапазона наблюдения
    spectrum = band_separation(spectrum, aver_param)

    # if noise_calibr == 'y':
    #     spectrum = noise_self_calibration(spectrum, ngi_temperature_path)
    # if black_body_calibr == 'y':
    #     spectrum = noise_black_body_calibration(spectrum, receiver_temperature_path)
    #
    # united_spectrum = unite_spectrum(spectrum)
    # ser_ind = united_spectrum.index
    # if len(ser_ind) == 2:
    #     spectrum_extr = united_spectrum[0] + united_spectrum[1]
    # else:
    #     spectrum_extr = united_spectrum[0]

    # # ***********************************************
    # # ***        Графический вывод данных        ****
    # # ***********************************************
    a = 0
    # Динамическая маска (зависит от длины записи во времени)
    t_spect = N_row * delta_t
    time_spect_mask = [(lambda i: (t_spect * (i + 0.05)) // 7)(i) for i in range(7)]

    # Формирование спектров и сканов по маскам freq_spect_mask и time_spect_mask
    spectr_freq, spectr_time = form_spectrum(spectrum, freq_spect_mask, time_spect_mask)
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
            spectr_freq[i, :] = del_random_mod(spectr_freq[i, :], 40)
        except IndexError:
            pass
    for i in range(n_time):
        spectr_time[i, :] = del_random_mod(spectr_time[i, :], 40)
    print('spectr_freq, spectr_time are formed')
    # Формирование строк-аргументов по времени и частоте и легенды
    N_col = np.shape(spectrum_extr)[1]
    if band_size_init == 'half':
        freq = np.linspace(1000 * (N_Nyq - 1) + 3.9063 / aver_param * kf, 1000 * N_Nyq - 3.9063 / aver_param * kf,
                           N_col // kf)
    elif band_size_init == 'whole':
        freq = np.linspace(1000 + 3.9063 / aver_param * kf, 3000 - 3.9063 / aver_param * kf, N_col // kf)
    timeS = np.linspace(0, delta_t * N_row, N_row // kt)
    # timeS = []

    # ***************!! Вывод данных в текстовой форме !!*********************
    # path_txt = str(Path(converted_dir_path, current_data_file, '_scan.txt'))
    # path_npy1 = Path(str(converted_data_file_path) + '_spectrum_time.npy') # Спектры в фиксированные моменты времени
    # path_npy2 = Path(str(converted_data_file_path) + '_scan_freq.npy')     # Сканы на фиксированных частотах
    # print(path_txt)
    # np.savetxt(path_txt, )
    # np.save(path_npy2, spectr_time)                                        # Сканы на фиксированных частотах
    # np.save(path_npy1, spectr_freq)                                        # Спектры в фиксированные моменты времени
    # path_txt = str(Path(converted_dir_path, current_data_file, 'freq.txt'))
    # print(path_txt)
    # np.savetxt(path_txt, freq)
    # ***********************************************************************

    # ************************** !!! Вывод данных !!! ***********************
    # line_legend_time, line_legend_freq = line_legend(freq_spect_mask[:10])
    if not 'good_bound' in head:
        head['good_bound'] = 0.1
    info_txt = [('time resol = ' + str(delta_t * kt) + 'sec'),
                ('freq resol = ' + str(delta_f / aver_param * kf) + 'MHz'),
                ('polarisation ' + polar), 'align: ' + align, 'kurtosis quality = ' + str(head['good_bound'])]
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
                                 path1,
                                 freq_spect_mask, head)
        else:
            fp.fig_multi_axes(spectr_time[:10, n_start:n_stop], timeS[n_start:n_stop], info_txt,
                                 path1,
                                 freq_spect_mask, head)

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
            fp.graph_contour_2d(freq[12:115], timeS[277:311], spectr_extr1[277:311, 12:115], 0, info_txt, path1, head)

    if poly3d_perm == 'y':
        data_poly3d, freq_mask = data_poly3d_prep(spectr_extr1)
        poly_graph3d(timeS, data_poly3d, freq_mask)

    stop = datetime.now()
    print('\n Total time = ', stop - start)
