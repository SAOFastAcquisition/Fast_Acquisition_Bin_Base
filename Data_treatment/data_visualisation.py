import numpy as np
import os
import sys
import pandas as pd
import pickle
import json as jsn
from datetime import datetime
from pathlib import Path
from Supporting_func import Fig_plot as fp, align_spectrum, path_to_data
# from Supporting_func import align_spectrum, path_to_data
from Interface import main
from Polyphase import low_freq_noise_spectrum, plot_low_freq_spec
from Interface.window_handler import exec_app
from Polyphase.cic_filter import signal_filtering
from test_statistic import low_noise_spectra_base

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
                s_time[i, j] = np.sum(spectr_extr[i * kt:(i + 1) * kt, ind])
                n_mesh = (spectr_extr[i * kt: (i + 1) * kt, ind] > 100).sum()
                if n_mesh == 0:
                    s_time[i, j] = 2
                else:
                    s_time[i, j] /= n_mesh
            else:
                s_time[i, j] = np.sum(spectr_extr[i * kt:(i + 1) * kt, ind - int(kf / 2):ind + int(kf / 2)])
                n_mesh = (spectr_extr[i * kt: (i + 1) * kt, ind - int(kf / 2):ind + int(kf / 2)] > 100).sum()
                if n_mesh == 0:
                    s_time[i, j] = 2
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
                s_freq[i, j] = np.sum(spectr_extr[ind, j * kf:(j + 1) * kf])
                n_mesh = (spectr_extr[ind, j * kf:(j + 1) * kf] > 10).sum()
                if n_mesh == 0:
                    s_freq[i, j] = 2
                else:
                    s_freq[i, j] /= n_mesh
            else:
                s_freq[i, j] = np.sum(spectr_extr[ind - int(kt / 2):ind + int(kt / 2), j * kf:(j + 1) * kf])
                n_mesh = (spectr_extr[ind - int(kt / 2):ind + int(kt / 2), j * kf:(j + 1) * kf] > 10).sum()
                if n_mesh == 0:
                    s_freq[i, j] = 2
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
                N_mesh = (Spectr[i * kt: (i + 1) * kt, j * kf: (j + 1) * kf] > 10).sum()
                if N_mesh == 0:
                    S1[i, j] = 2
                else:
                    S1[i, j] = S1[i, j] / N_mesh
                if S1[i, j] == 0:
                    S1[i, j] = 2
                # if (j > 3) & (S1[i, j] > 1.5 * np.sum(S1[i, j-3:j])//3):
                #     S1[i, j] = np.sum(S1[i, j-3:j])//3
                if robust_filter == 'y':
                    a = param_robust_filter
                    if (i > 3) & (S1[i, j] < 1 / a * np.sum(S1[i - 3:i - 1, j]) // 2):
                        S1[i, j] = np.sum(S1[i - 1, j])
                    if (i > 3) & (S1[i, j] > a * np.sum(S1[i - 3:i - 1, j]) // 2):
                        # print(S1[i - 3:i+1, j])
                        S1[i, j] = np.sum(S1[i - 1, j])
                        # print(S1[i, j])
                        pass

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
    _path1 = Path(converted_data_file_path, current_primary_file + '_spectrum.npy')
    spectrum = np.load(_path1, allow_pickle=True)
    with open(Path(converted_data_file_path, current_primary_file + '_head.bin'), 'rb') as inp:
        head = pickle.load(inp)
    n_aver = head['n_aver']
    band_size = head['band_size']
    polar = head['polar']
    # n_aver = head['n_aver']

    # Разделяем составляющие  записи в полной полосе и с возможными двумя поляризациями,
    # одновременно понижая разрядность данных, меняя их тип с int64 до int32 и уменьшая
    # занимаемую ими память
    if num_of_polar == 2 and band_size == 'whole':
        # if np.size(spectrum[0]) > 1:
        #     spectrum_left1 = (spectrum[0] / 1000).astype(np.int32)
        # else:
        spectrum_left1 = spectrum[0]
        # if np.size(spectrum[1]) > 1:
        #     spectrum_left2 = (spectrum[1] / 1000).astype(np.int32)
        # else:
        spectrum_left2 = spectrum[1]
        # if np.size(spectrum[2]) > 1:
        #     spectrum_right1 = (spectrum[2] / 1000).astype(np.int32)
        # else:
        spectrum_right1 = spectrum[2]
        # if np.size(spectrum[3]) > 1:
        #     spectrum_right2 = (spectrum[3] / 1000).astype(np.int32)
        # else:
        spectrum_right2 = spectrum[3]
    # Выдает спектры для левой и правой поляризаций шириной по 1 ГГц. При нумерации спектров учитывается
    # значение зоны Найквиста. С индексом "1" - 1-2 ГГц, с индексом "2" - 2-3 ГГц, как и для случая выше.
    # На выходе формально все 4 спектра, но для незадействованной полосы они пустые
    elif num_of_polar == 2 and band_size == 'half':
        if N_Nyq == 2:
            spectrum_left1 = spectrum[0]
            spectrum_right1 = spectrum[1]
            spectrum_left2 = []
            spectrum_right2 = []
        else:
            spectrum_left2 = spectrum[0]
            spectrum_right2 = spectrum[1]
            spectrum_left1 = []
            spectrum_right1 = []
    pass

    return spectrum_left1, spectrum_left2, spectrum_right1, spectrum_right2, int(n_aver), band_size, polar


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
        n_row1 = np.min([np.shape(spec_left)[0], np.shape(spec_right)[0]])
        spec_left = spec_left[:n_row1]
        spec_right = spec_right[:n_row1]
        united_spec = pd.Series([spec_left, spec_right], ind)
    elif np.size(spec_left):
        united_spec = pd.Series([spec_left], ind)
    else:
        united_spec = pd.Series([spec_right], ind)
    print('Spectra are united')
    return united_spec


def freq_mask(_i):
    _n1 = 1
    _n2 = 9
    _freq_mask = [
        [1900],                                                               # [0]
        [2060, 2220, 2300, 2500, 2560, 2700, 2800, 2880, 2980],               # [1]
        [1080, 1140, 1360, 1420, 1620, 1780, 1980],                           # [2]
        [1000 * _n1 + 100 * _n2 + 10 * i for i in range(10)],                 # [3]
        [1050, 1465, 1535, 1600, 1700, 2265, 2550, 2700, 2800, 2920],         # [4]
        [1230, 1560, 2300, 2910],                                                               # [5]
        [1140, 1420, 1480, 2460, 2500, 2780],   # for Crab '2021-06-28_03+14' # [6]
        [1220, 1540, 1980, 2060, 2500, 2780],   # for Crab '2021-06-28_04+12' # [7]
        [1171, 1380, 1465, 1600, 1700, 2265, 2530, 2720, 2800, 2920]    # [8]
    ]
    return _freq_mask[_i]


def data_poly3d_perm(_spectrum_extr):
    _k, _m = np.shape(_spectrum_extr)
    freq_mask_poly = [1100, 1200, 1700, 1800, 2800, 2900]
    i = 0
    _n = [0] * 6
    _data_poly3d = []
    for s in freq_mask_poly:
        _n[i] = int((s - 1000) / 2000 * _m)
        if i % 2 == 1:
            if len(_data_poly3d) ==0:
                _data_poly3d = _spectrum_extr[:, _n[i - 1]:_n[i]]
            else:
                _data_poly3d = np.hstack((_data_poly3d, _spectrum_extr[:, _n[i - 1]:_n[i]]))

        i += 1



    pass


if __name__ == '__main__':

    # parameters = main()
    # current_data_file = parameters['file_name']      # Имя файла с исходными текущими данными без расширения
    # current_data_dir = parameters['file_folder']          # Папка с текущими данными
    # freq_res = parameters['freq_res']  # Установка разрешения по частоте в МГц
    # kt = parameters['time_res'] // 8  # Установка разрешения по времени в единицах минимального разрешения
    # 8.1925e-3 сек
    # output_picture_mode = parameters['output_picture_mode'] == 'yes'
    align_file_name = 'Align_coeff.bin'         # Имя файла с текущими коэффициентами выравнивания АЧХ
    # current_catalog = r'2021/Results'           # Текущий каталог (за определенный период, здесь - год)

    current_data_dir = '2022'
    primary_data_dir = 'Primary_data'           # Каталог исходных данных (за определенный период, здесь - год)
    converted_data_dir = 'Converted_data'       # Каталог для записи результатов конвертации данных и заголовков
    data_treatment_dir = 'Data_treatment'       # Каталог для записи результатов обработки, рисунков

    current_primary_dir = '2022_03_26sun'
    current_converted_dir = current_primary_dir + '_conv'
    current_converted_path = Path(converted_data_dir, current_converted_dir)
    current_treatment_dir = current_primary_dir + '_treat'
    current_treatment_path = Path(data_treatment_dir, current_treatment_dir)

    current_primary_file = '2022-03-26_04+16'

    converted_data_file_path, head_path = path_to_data(current_data_dir, current_converted_path)
    data_treatment_file_path, head_path = path_to_data(current_data_dir, current_treatment_path)

    folder_align_path = Path(head_path, 'Alignment')
    date = current_primary_file[0:10]

    # !!!! ******************************************* !!!!
    # ****** Блок исходных параметров для обработки *******

    freq_res = 30  # Установка разрешения по частоте в МГц
    kt = 30  # Установка разрешения по времени в единицах минимального разрешения 8.1925e-3 сек
    delta_t = 8.3886e-3
    delta_f = 7.8125
    N_Nyq = 3
    att_val = [i * 0.5 for i in range(64)]
    att_dict = {s: 10 ** (s / 10) for s in att_val}
    freq_spect_mask = freq_mask(3)
    # *****************************************************

    band_size_init = 'whole'
    num_of_polar = 2
    # band_size = 'whole'   Параметр 'whole' означает работу в диапазоне 1-3 ГГц, 'half' - диапазон 1-2 или 2-3 ГГц
    # polar = 'both'        Принимает значения поляризаций: 'both', 'left', 'right'
    # *****************************************************
    output_picture_mode = 'b'
    align = 'y'  # Выравнивание АЧХ усилительного тракта по калибровке от ГШ: 'y' / 'n'
    noise_calibr = 'n'
    save_data = 'n'     # Сохранение сканов в формате *.npy: 'y' / 'n'
    lf_filter = 'n'     # Применение НЧ фильтра для сглаживания сканов (скользящее среднее и др.): 'y' / 'n'
    low_noise_spectrum = 'n'    # Вывод графика НЧ спектра шумовой дорожки: 'y' / 'n'
    robust_filter = 'n'
    graph_3d_perm = 'n'
    contour_2d_perm = 'n'
    poly3d_perm = 'y'

    # *****************************************************
    # Чтение с диска, если спектры ранее извлекались,
    # или извлечение спектров из исходных записей
    spectr_extr_left1, spectr_extr_left2, spectr_extr_right1, spectr_extr_right2, n_aver, band_size, polar = \
        preparing_data()
    print('Data are prepared')
    aver_param = 2 ** (6 - n_aver)
    kf = int(freq_res / delta_f * aver_param)   # Установка разрешения по частоте в единицах максимального разрешения
                                                # для данного наблюдения delta_f/aver_param, где delta_f = 7.8125 МГц
    with open(Path(converted_data_file_path, current_primary_file + '_head.bin'), 'rb') as inp:
        head = pickle.load(inp)

    # Выравнивание спектров по результатам шумовых измерений АЧХ
    if align == 'y':
        if head['att3'] == 5:
            pos = 0
        elif head['att3'] == 0:
            pos = 0
        else:
            pos = 1
        path_output = Path(folder_align_path, align_file_name)
        spectr_extr_left1, spectr_extr_left2, spectr_extr_right1, spectr_extr_right2 = \
            align_spectrum(spectr_extr_left1, spectr_extr_left2, spectr_extr_right1, spectr_extr_right2,
                           head, path_output, pos)
        print('spectrum is aligned')

    # Приведение порядка следования отсчетов по частоте к нормальному
    if np.size(spectr_extr_left1):
        N_row = np.shape(spectr_extr_left1)[0]
        for i in range(N_row):
            spectr_extr_left1[i][0:] = spectr_extr_left1[i][-1::-1]
    if np.size(spectr_extr_right1):
        N_row = np.shape(spectr_extr_right1)[0]
        for i in range(N_row):
            spectr_extr_right1[i][0:] = spectr_extr_right1[i][-1::-1]

    spectrum = pd.Series([spectr_extr_left1, spectr_extr_left2, spectr_extr_right1, spectr_extr_right2])

    united_spectrum = unite_spectrum(spectrum)
    ser_ind = united_spectrum.index
    if len(ser_ind) == 2:
        spectrum_extr = united_spectrum[0] + united_spectrum[1]
    else:
        spectrum_extr = united_spectrum[0]

    # if noise_calibr == 'y':
    #     spectr_time = calibration(t_cal, spectr_time)

    # # ***********************************************
    # # ***        Графический вывод данных        ****
    # # ***********************************************

    # Динамическая маска (зависит от длины записи во времени)
    t_spect = N_row * delta_t
    time_spect_mask = [(lambda i: (t_spect * (i + 0.05)) // 7)(i) for i in range(7)]
    # time_spect_mask = [5, 17, 31, 43, 58]
    # if band_size == 'whole':
    #   freq_spect_mask = []

    # Формирование спектров и сканов по маскам freq_spect_mask и time_spect_mask
    shift = head['shift']
    spectr_freq, spectr_time = form_spectr_sp1(spectrum_extr, freq_spect_mask, time_spect_mask)
    print('spectr_freq, spectr_time are formed')
    # Формирование строк-аргументов по времени и частоте и легенды
    N_col = np.shape(spectrum_extr)[1]
    if band_size_init == 'half':
        freq = np.linspace(1000 * (N_Nyq - 1) + 3.9063 / aver_param * kf, 1000 * N_Nyq - 3.9063 / aver_param * kf,
                           N_col // kf)
    elif band_size_init == 'whole':
        freq = np.linspace(1000 + 3.9063 / aver_param * kf, 3000 - 3.9063 / aver_param * kf, N_col // kf)
    timeS = np.linspace(0, delta_t * N_row, N_row // kt)

    # ***************!! Вывод данных в текстовой форме !!*********************
    # path_txt = str(Path(file_path_data, current_data_file, '_scan.txt'))
    # print(path_txt)
    # np.savetxt(path_txt, spectr_freq)
    # path_txt = str(Path(file_path_data, current_data_file, 'freq.txt'))
    # print(path_txt)
    # np.savetxt(path_txt, freq)
    # ***********************************************************************

    # ************************** !!! Вывод данных !!! ***********************
    line_legend_time, line_legend_freq = line_legend(freq_spect_mask[:10])
    info_txt = [('time resol = ' + str(delta_t * kt) + 'sec'),
                ('freq resol = ' + str(delta_f / aver_param * kf) + 'MHz'),
                ('polarisation ' + polar), 'align: ' + align]
    path1 = Path(data_treatment_file_path, current_primary_file)
    path_to_fig(data_treatment_file_path)
    path_to_fig(path1)
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
        spectr_time = signal_filtering(spectr_time, 0.003)
    if low_noise_spectrum == 'y':
        spectrum_signal_av = low_freq_noise_spectrum(spectr_time, 32768)
        if kt == 1 & kf == 1:
            m, n = spectrum_signal_av.shape
            f_max = 1 / delta_t / 2
            f_min = f_max / n
            arg = np.linspace(f_min, f_max, n)
            low_noise_spectra_base(spectrum_signal_av, head, freq_spect_mask, arg, current_primary_file)
        plot_low_freq_spec(spectrum_signal_av, delta_t * kt, path1, line_legend_freq)

    if output_picture_mode == 'y':
        fp.fig_plot(spectr_freq, 0, freq, 1, info_txt, path1, head, line_legend_time)
        fp.fig_plot(spectr_time, 0, timeS, 0, info_txt, path1, head, line_legend_freq)

    # *********************************************************
    # ***            Многооконный вывод данных             ****
    # *********************************************************
    if output_picture_mode == 'no':
        fp.fig_multi_axes(spectr_time, timeS, info_txt, Path(current_treatment_path, current_primary_file),
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
        fp.graph_3d(freq, timeS, spectr_extr1, 0, current_data_file, head)
    if contour_2d_perm == 'y':
        fp.graph_contour_2d(freq, timeS, spectr_extr1, 0)

    if poly3d_perm == 'y':
        data_poly3d = data_poly3d_perm(spectr_extr1)

    # if align == 'y':
    #     align_coeff1 = align_func1(spectr_freq[1, :], 'y', aver_param)
    #     spectr_extr = spectr_extr * align_coeff1

    # if graph_3d_perm == 'y':
    #     graph_3d(freq, timeS[n_start_flame:n_stop_flame], spectr_extr1[n_start_flame:n_stop_flame, :], 0)
    # fp.fig_multi_axes(spectr_time[:10, n_start_flame:n_stop_flame], timeS[n_start_flame:n_stop_flame],
    #                   info_txt, file_name0, freq_spect_mask[:10])


    stop = datetime.now()
    print('\n Total time = ', stop - start)
