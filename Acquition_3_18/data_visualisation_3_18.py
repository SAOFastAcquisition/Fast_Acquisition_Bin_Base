import numpy as np
import os
import sys
import pandas as pd
import pickle
import json as jsn
from datetime import datetime
from pathlib import Path
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objs as go
import matplotlib.font_manager as font_manager
from Supporting_func import path_to_data
# from Supporting_func import align_spectrum, path_to_data
from Polyphase import low_freq_noise_spectrum, plot_low_freq_spec, plot_low_freq_spec_ab
from Polyphase.cic_filter import signal_filtering
# from test_statistic import low_noise_spectra_base
from polys3d_demo import poly_graph3d

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
            _main_dir = f'{_main_dir}'
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


def cut_spectrum(_spectrum, _n_aver):
    _spectrum.pop(-1)
    n_frame_last = _spectrum[-1][0]
    rest = (n_frame_last + 1) % (256 / _n_aver)
    if rest:
        for k in range(rest):
            _spectrum.pop(-1)
    print(n_frame_last, _spectrum[-1][0])
    return _spectrum


def band_map():
    bm = [
        [1000, 2900],
        [2900, 5200],
        [5200, 6500],
        [6500, 8800],
        [8800, 11100],
        [11100, 13400],
        [13400, 15700],
        [15700, 18000]
    ]
    return bm


def form_spectrum(_sp, _freq_spect_mask_in=freq_spect_mask, _time_mask=time_spect_mask):
    """ Возвращает s_freq - срезы частотного спектра в моменты времени time_spect_mask и s_time - сканы Солнца
    по времени на частотах freq_spect_mask с заданным разрешением по времени и частоте

    """
    # Определение заполненных поддиапазонов
    _shape_sp = np.shape(_sp)
    _map_sp = []  # Номера поддиапазонов, в которых проводилось наблюдение

    for _i in range(_shape_sp[0]):
        for _j in range(_shape_sp[1]):
            _flag = 0
            if np.size(_sp[_i, _j]) > 2:
                _flag = 1
            if _flag and _i not in _map_sp:
                _map_sp.append(_i)

    # Определение наличия данных для сканов на частоах из freq_spect_mask
    _freq_map = {}  # Частота анализа - ключ, номер поддиапазона - значение
    _band_map = band_map()
    for _i in _map_sp:
        for _f in _freq_spect_mask_in:
            if _band_map[_i][0] <= _f <= _band_map[_i][1]:
                _freq_map[_f] = _i
    _freq_mask = list(_freq_map.keys())  # Частоты, для которых доступны сканы во времени

    # Словарь, в котором ключ - номер поддиапазона, значение - количество частотных отсчетов в нем
    try:
        _N_col: dict = {int(_s): np.shape(_sp[_s, 0])[1] for _s in _map_sp}
    except IndexError:
        _N_col: dict = {int(_s): np.shape(_sp[_s, 1])[1] for _s in _map_sp}
    # Словарь, в котором ключ - номер поддиапазона, значение - количество временных отсчетов в нем
    try:
        _N_raw: dict = {int(_s): np.shape(_sp[_s, 0])[0] for _s in _map_sp}
    except IndexError:
        _N_raw: dict = {int(_s): np.shape(_sp[_s, 1])[0] for _s in _map_sp}
    _kf = int(freq_res / delta_f / aver_param)
    if not _kf:
        sys.exit(f'Frequency resolution freq_res={freq_res}MHz is too large. Decrease it, '
                 f'increase freq_res to {delta_f * aver_param}MHz at least')

    _scan = scan_former(_sp, _freq_map, _freq_mask, _band_map, _N_col, _N_raw)
    _spectrum = spectrum_former(_sp, _freq_map, _time_mask, _band_map, _N_col, _N_raw)

    return _spectrum, _scan


def scan_former(_sp, _freq_map, _freq_mask, _band_map, _n_col, _n_raw):
    """
    Формирование сканов с заданным разрешением по частоте и времени на фиксированных частотах из
    _freq_spect_mask_in.  В результате формируется список из списков, количество которых равно
    удвоенному числу частот в _freq_spect_mask_in (удвоение из-за двух поляризаций),
    первый элемент - частота скана,
    второй - поляризация (0 - левая, 1 - правая),
    третий - значения спектра скана,
    четвертый - моменты времени, в которые берутся значения спектров скана (шкала, аргумент)
    и так для каждой частоты в маске
    j   - индекс частоты в маске
    _i  - индекс по поляризации (0 или 1)
    i   - индекс по моментам времени, взятым с заданным разрешением
    :return: list, float64
    """
    #   Заготовка выдаваемого списка
    s_time = [
        [[] for i in range(4)] for j in range(len(_freq_mask) * 2)
    ]

    j = 0
    for _f in _freq_mask:
        ind = int((_f - _band_map[int(_freq_map[_f])][0] - delta_f * aver_param / 2) // (delta_f * aver_param))
        # - количество отсчетов по частоте в поддиапазоне с разрешением не лучше freq_res
        if ind > _n_col[int(_freq_map[_f])] - int(kf / 2) - 1:
            ind = _n_col[int(_freq_map[_f])] - int(kf / 2) - 1
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
                if _freq_map[_f] == 2:  # Для зоны 5100 - 6500 МГц разрешение по частоте и времени другое
                    _delta_t = 0.013981
                else:
                    _delta_t = 0.012428  # Стандартное разрешение по времени
                while kt * (i + 1) < _n_col[int(_freq_map[_f])]:
                    if kf == 1:
                        _s_loc = np.sum(sp[i * kt:(i + 1) * kt, ind][sp[i * kt:(i + 1) * kt, ind] > 40])
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
                    _time_count.append(kt * (i + 0.5) * _delta_t)
                    _s_loc1.append(_s_loc)
                    i += 1

                s_time[j][0] = _f  # scan frequency
                s_time[j][1] = _i  # scan polarization
                s_time[j][2] = np.array(_s_loc1)  # scan
                s_time[j][3] = np.array(_time_count)  # time sequence

                j += 1
    return s_time


def spectrum_former(_sp, _freq_map, _time_mask, _band_map, _n_col, _n_raw):
    """
    Формирование спектров с заданным разрешением по частоте и времени в фиксированные моменты времени  из
    _time_mask.  В результате формируется список из списков, количество которых равно
    удвоенному числу моментов времени в _time_mask (удвоение из-за двух поляризаций),
    первый элемент - время взятия спектра (маска по времени),
    второй - поляризация (0 - левая, 1 - правая),
    третий - значения спектра,
    четвертый - отсчеты частоты, на которых берутся значения спектров (частотная шкала, аргумент)
    и так для каждого значения из _time_mask и для двух поляризаций.
    _l  - индекс по поляризации (0 или 1)
    i   - индекс по моментам времени в маске
    :param _sp:
    :param _freq_map:
    :param _time_mask:
    :param _band_map: карта границ поддиапазонов
    :param _n_col:
    :param _n_raw:
    :return: list, float64
    """
    _s_freq = [
        [[] for i in range(8)] for j in range(len(_time_mask) * 2)
    ]
    i = 0
    for t in _time_mask:
        # _j - Номер частотного поддиапазона
        for _l in range(2):  # Цикл по поляризациям
            _s_loc1 = []
            _freq_count = []
            for _j in range(8):  # Цикл по поддиапазонам
                _delta_f = 0.082397
                if _j == 2:
                    _delta_f = 0.073242
                sp = _sp[_j, _l]
                if len(sp) > 2:
                    # Определение индекса середины зоны, в которой проводится усреднение по времени
                    ind = int(t // delta_t)
                    if ind > _n_raw[_j] - kt / 2 - 1:
                        ind = _n_raw[_j] - int(kt / 2) - 1
                    if ind < (kt / 2):
                        ind = int(kt / 2)
                    j = 0
                    while (j + 1) * kf < _n_col[_j]:  # Цикл по частоте - аргументу спектра
                        if kt == 1:
                            s_loc = np.sum(sp[ind, j * kf:(j + 1) * kf][sp[ind, j * kf:(j + 1) * kf] > 40])
                            n_mesh = (sp[ind, j * kf:(j + 1) * kf] > 40).sum()
                            if n_mesh == 0:
                                s_loc = 2
                            else:
                                s_loc /= n_mesh
                        else:
                            s_loc = np.sum(sp[ind - int(kt / 2):ind + int(kt / 2), j * kf:(j + 1) * kf][
                                               sp[ind - int(kt / 2):ind + int(kt / 2), j * kf:(j + 1) * kf] > 40])
                            n_mesh = (sp[ind - int(kt / 2):ind + int(kt / 2), j * kf:(j + 1) * kf] > 40).sum()
                            if n_mesh == 0:
                                s_loc = 2
                            else:
                                s_loc /= n_mesh
                            _s_loc1.append(s_loc)
                            _freq_count.append(_band_map[_j][0] + kf * (j + 0.5) * _delta_f * aver_param)
                        j += 1
                        if _j == 2:
                            _a = _freq_count[-1]
                            pass
            _s_freq[i][0] = t
            _s_freq[i][1] = _l
            _s_freq[i][2] = np.array(_s_loc1)
            _s_freq[i][3] = np.array(_freq_count)
            i += 1
    return _s_freq


def spectrum_construction(Spectr, kf, kt):
    """ Функция формирует спектр принятого сигнала с требуемым разрешением по частоте и времени. Максимальное
    разрешение отсчетов по времени порядка 12 мсек и по частоте 7,8125 МГц. Путем суммирования и усреднерия по kt*kf
    отсчетам разрешение по частоте и по времени в исходном спектре Spectr уменьшается в kf и kt раз,
    соответственно. Преобразованный спектр возвращается как S1. Для трехмерного изображения
    """

    N_col1 = n_col // kf
    N_row1 = n_row // kt
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

            except IndexError as allert_message:
                print(allert_message, 'ind i = ', i, 'ind j = ', j)
                pass
            except ValueError as value_message:
                print(value_message, 'ind i = ', i, 'ind j = ', j)
                pass

    return S1


def path_to_fig(_path):
    """ Создает директорию для рисунков обрабатываемого наблюдения, если она до этого не была создана,
    название директории  совпадает с названием исходного файла данных наблюдения
    """
    if not os.path.isdir(_path):
        os.mkdir(_path)
    return


def preparing_data():
    """ Функция выдает данные для построения графиков"""

    # Для k-ой полосы и двух возможных поляризаций выдает по спектру для каждой поляризации.
    # Если поляризация не задействована, то соответствующие спектры - пустые. Спектры k = 3,4 - в обратном порядке
    _path1 = Path(str(converted_data_file_path) + '_spectrum.npy')
    _spectrum = np.load(_path1, allow_pickle=True)
    with open(Path(str(converted_data_file_path) + '_head.bin'), 'rb') as inp:
        _head = pickle.load(inp)
    if '_n_aver' in _head:
        _n_aver = _head['_n_aver']
    else:
        _n_aver = _head['n_aver']
    _polar = _head['polar']

    # Разделяем составляющие  записи в полной полосе и с возможными двумя поляризациями
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
    _n_start0 = int(2432 / _n_aver)
    _n_stop0 = int(32768 / _n_aver) - _n_start0
    _n_start2 = int(4096 / _n_aver)
    _n_stop2 = int(21844 / _n_aver)
    _shape_s = np.shape(_s)
    for i in range(_shape_s[0]):
        _n_start = _n_start0
        _n_stop = _n_stop0
        if i == 2:
            _n_start = _n_start2
            _n_stop = _n_stop2
        for j in range(_shape_s[1]):
            if np.size(_s[i, j]):
                _N_row = np.shape(spectrum[i, j])[0]
                s1 = np.zeros((_N_row, _n_stop - _n_start))
                s = _s[i, j]
                for k in range(_N_row):
                    try:
                        s1[k, :] = s[k, _n_start:_n_stop]
                    except ValueError:
                        s1[k, :] = s1[k, :].fill(2)
                _s[i, j] = s1

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


def plotly_form(*args):
    _data, _legend, _info = args

    _fig = go.Figure()
    _m = len(_data)
    _l = 0
    for i in range(_m):
        if len(_data[i][2]) > 2:
            _fig.add_trace(go.Scatter(x=_data[i][3], y=_data[i][2], name=_legend[i]))
            _l += 1
    _fig.add_annotation(text=_info[0], xref="paper", yref="paper", x=1, y=0.95, showarrow=False)
    _fig.add_annotation(text=_info[1], xref="paper", yref="paper", x=1, y=0.9, showarrow=False)
    _fig.add_annotation(text=_info[4], xref="paper", yref="paper", x=1, y=0.85, showarrow=False)

    _fig.update_yaxes(type="log")
    _fig.show(renderer='browser')


def simple_plot(_inp_data):
    fig, ax = plt.subplots(1, figsize=(12, 6))
    line_color = ['green', 'blue', 'purple', 'lime', 'black', 'red', 'olivedrab', 'lawngreen', 'magenta', 'dodgerblue']

    _unit = [' sec', ' MHz']
    _t = 1
    if _inp_data[0, 0] < 1000:
        _t = 0
    legend = [str(_i) + _unit[_t] for _i in _inp_data[:, 0]]

    _i = 0
    _l = 0
    _m = np.shape(_inp_data)[0]
    for _i in range(_m):
        if len(_inp_data[_i][2]) > 2:
            ax.semilogy(_inp_data[_i][3], _inp_data[_i][2], color=line_color[_l], label=legend[_i])
            _l += 1

    font = font_manager.FontProperties(family='Comic Sans MS',
                                       weight='bold',
                                       style='normal', size=16)

    ax.legend(loc=10, prop=font, bbox_to_anchor=(1, 0.5))
    # line_minor = plt.plot(inp_data[1][3] / 1000)

    # plt.plot(inp_scale, inp_data, 'r--*', inp_scale, 'g:+')  # В кавычках компактно заданы цвет, стиль линии и маркеры
    # plt.setp(fig_main, linestyle=':', color='r')
    # plt.setp(fig_main, linestyle=':', color=(0, 1, 0, 0.9), marker='.', markerfacecolor='b')  # Четвертое число
    # задает прозрачность линии
    # plt.setp(line_minor, linestyle='-.')
    plt.grid()
    plt.show()


def freq_mask(_i):
    _n1 = 2
    _n2 = 4
    _freq_mask = [
        [3000, 3600, 4500, 5100, 6600, 7500, 8200, 8600],  # [0]
        [3700, 4850, 6350, 6500, 7000, 7100],  # [1]
        [5100, 6400, 7000, 8500, 9000, 1000, 11000, 13000]  # [2]
    ]
    return _freq_mask[_i]


if __name__ == '__main__':

    align_file_name = 'antenna_temperature_coefficients.npy'  # Имя файла с текущими коэффициентами выравнивания АЧХ
    current_primary_file = '2024-04-10_01'
    current_primary_dir = '2024_04_10test'
    main_dir = '2024'
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

    freq_res = 4  # Установка разрешения по частоте в МГц
    kt = 4  # Установка разрешения по времени в единицах минимального разрешения delta_t
    delta_t = 12.427567e-3  # sec
    delta_f = 0.082397461  # MHz

    att_val = [i * 0.5 for i in range(64)]
    att_dict = {s: 10 ** (s / 10) for s in att_val}
    pol_dict = {0: 'left', 1: 'right'}
    freq_spect_mask = freq_mask(2)
    # *****************************************************
    polar = 'both'
    # polar = 'both'        Принимает значения поляризаций: 'both', 'left', 'right'
    # *****************************************************

    a = 0
    # *****************************************************
    # Чтение с диска, если спектры ранее извлекались,
    # или извлечение спектров из исходных записей
    spectrum, n_aver, polar = preparing_data()
    print('Data are prepared')
    aver_param = n_aver  # Количество усредняемых частотных бинов для разрешения по частоте в записи
    kf = int(freq_res / delta_f / n_aver)  # Установка разрешения по частоте в единицах максимального разрешения
    # для данного наблюдения delta_f*aver_param, где delta_f = 0.082397461 МГц
    with open(Path(str(converted_data_file_path) + '_head.bin'), 'rb') as inp:
        head = pickle.load(inp)

    # Приведение порядка следования отсчетов по частоте к нормальному,
    # условно, требуют смены порядка записи 1 и 4 - по совокупности признаков
    shape_sp = np.shape(spectrum)
    for i in range(shape_sp[0]):
        for j in range(shape_sp[1]):
            if np.size(spectrum[i, j]) and i in [2]:
                n_row = np.shape(spectrum[i, j])[0]
                n_col = np.shape(spectrum[i, j])[1]
                s = spectrum[i, j]
                for k in range(n_row):
                    s[k, :] = s[k, -1::-1]
                spectrum[i, j] = s
                s1 = spectrum[i, j]

    # Отбрасывание отсчетов защитной зоны поддиапазона наблюдения
    spectrum = band_separation(spectrum, n_aver)

    # # ***********************************************
    # # ***        Графический вывод данных        ****
    # # ***********************************************

    # Динамическая маска (зависит от длины записи во времени)
    t_spect = n_row * delta_t
    time_spect_mask = [(lambda i: (t_spect * (i + 0.05)) // 7)(i) for i in range(7)]

    # Формирование спектров и сканов по маскам freq_spect_mask и time_spect_mask
    spectrum_f, scan = form_spectrum(spectrum, freq_spect_mask, time_spect_mask)

    n_freq = len(time_spect_mask)
    n_time = len(freq_spect_mask)

    print('spectr_freq, spectr_time are formed')

    # ************************** !!! Вывод данных !!! ***********************

    info_txt = [f'time resol = {delta_t * kt:.3f} sec',
                f'freq resol = {delta_f * aver_param * kf:.1f} MHz',
                f'polarisation {polar}', 'align: ' + '1', 'kurtosis quality = ' + str(head['good_bound'])]
    path1 = data_treatment_file_path

    # if output_picture_mode == '_y':
    # sp = spectrum_f[1:6]
    # sp1 = [s[2] for s in sp]
    # simple_plot(spectrum_f)
    # simple_plot(scan)
    unit = [' sec', ' MHz']
    _t = 1
    if spectrum_f[0][0] < 1000:
        _t = 0
    legend_sp = [str(s[0]) + unit[0] + ' ' + pol_dict[s[1]] for s in spectrum_f]    # for spectrum
    legend_sc = [str(s[0]) + unit[1] + ' ' + pol_dict[s[1]] for s in scan]          # for scans
    plotly_form(spectrum_f, legend_sp, info_txt)
    plotly_form(scan, legend_sc, info_txt)



    # *********************************************************
    # ***            Многооконный вывод данных             ****
    # *********************************************************

    # *********************************************************
    # ***        Вывод данных двумерный и трехмерный       ****
    # *********************************************************
    # Укрупнение  разрешения по частоте и времени для вывода в 2d и 3d

    stop = datetime.now()
    print('\n Total time = ', stop - start)
