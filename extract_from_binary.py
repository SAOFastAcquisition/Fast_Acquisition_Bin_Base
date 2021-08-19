import numpy as np
import os
import sys
import pandas as pd
import pickle
import json as jsn
# import matplotlib.pyplot as plt
from datetime import datetime
from pathlib import Path
from pandas._libs import json
from pathlib import Path
from Supporting_func import Fig_plot as fp
# from Supporting_func import stat_cleaning
from Supporting_func.afc_alignment import align_spectrum
from Supporting_func.stocks_coefficients import path_to_data

current_dir = Path.cwd()
home_dir = Path.home()

# from Supporting_func.afc_alignment1 import align_func1

sys.path.insert(0, Path(current_dir, 'Supporting_func'))
start = datetime.now()

current_data_file = '2021-06-28_14-20'      # Имя файла с исходными текущими данными без расширения
current_data_dir = '2021_06_28sun'          # Папка с текущими данными
align_file_name = 'Align_coeff.bin'         # Имя файла с текущими коэффициентами выравнивания АЧХ
current_catalog = r'2021\Results'           # Текущий каталог (за определенный период, здесь - год)

file_path_data, head_path = path_to_data(current_catalog, current_data_dir)
folder_align_path = Path(head_path, 'Alignment')
date = current_data_file[0:11]

# !!!! ******************************************* !!!!
# ****** Блок исходных параметров для обработки *******
kf = 4  # Установка разрешения по частоте
kt = 4  # Установка разрешения по времени
N_Nyq = 2   # Номер зоны Найквиста
shift = 10  # Усечение младших разрядов при обработке первичного бинарного файла данных
# *****************************************************

delta_t = 8.1925e-3
delta_f = 7.8125
num_of_polar = 2  # Параметр равен "1" для записей до 12.12.2020 и "2" для записей после 12.12.2020
if num_of_polar == 1:
    q = int(current_data_file[-1])
    N_Nyq = q
band_size_init = 'whole'
# band_size = 'whole'   Параметр 'whole' означает работу в диапазоне 1-3 ГГц, 'half' - диапазон 1-2 или 2-3 ГГц
# polar = 'both'        Принамает значения поляризаций: 'both', 'left', 'right'
robust_filter = 'n'
param_robust_filter = 1.1
align = 'y'  # Выравнивание АЧХ усилительного тракта по калибровке от ГШ 'y' / 'n'

noise_calibr = 'n'
graph_3d_perm = 'n'
contour_2d_perm = 'n'

if N_Nyq == 3:
    freq_spect_mask = [2120, 2300, 2700, 2820, 2900]  # 2060, 2750, 2760, 2770, 2780, 2790, 2800, 2810,
    # 2820, 2830, 2850, 2880, 2900, 2950# Временные сканы Солнца на этих частотах
elif band_size_init == 'whole':
    n1 = 1
    n2 = 1
    # freq_spect_mask = [2100, 2265, 2460, 2550, 2700, 2735, 2750, 2800, 2920]
    # freq_spect_mask = [1535,  2450, 2550, 2750,  2800, 2950]
    # freq_spect_mask = [1000 * n1 + 100 * n2 + 20 * i for i in range(10)]
    # freq_spect_mask = [1050, 1171, 1380, 1465, 1500, 1535, 1600, 1700, 1750, 1950]
    # freq_spect_mask = [1050, 1465, 1500, 1535, 1600, 1700, 1750, 1950]
    freq_spect_mask = [1050, 1465, 1535, 1600, 1700, 2265, 2550, 2700, 2800, 2920]
else:
    freq_spect_mask = [1050, 1171, 1380, 1465, 1535, 1600, 1700, 2265, 2550, 2700, 2800, 2920]

time_spect_mask = [47, 84.4, 104, 133, 133.05, 177.02, 177.38]  # Срез частотного спектра в эти моменты времени

att_val = [i * 0.5 for i in range(64)]
att_dict = {s: 10 ** (s / 10) for s in att_val}
pass





def extract():
    file_name = Path(file_path_data, current_data_file + '.bin')
    file_name_out = Path(file_path_data, current_data_file + '.txt')
    i = 0
    k = 0
    spectr = []
    frame = ' '

    try:
        if os.path.isfile(file_name) == 1:
            pass
        else:
            print('\n \t', file_name, ' not found!!!\n')

        f_in = open(file_name, 'rb')

        while frame:

            spectr_frame = []
            # Обработка кадра: выделение номера кадра, границ куртозиса, длины усреднения на ПЛИС
            # и 128-ми значений спектра в список spectr_frame на позиции [1:128]
            for k in range(129):
                frame = f_in.read(8)
                frame_int = int.from_bytes(frame, byteorder='little')
                if k == 0:
                    frame_num = frame_int & 0xFFFFFFF

                    # Выделение длины усреднения (количество усредняемых на ПЛИС отсчетов спектра = 2^n_aver)
                    # Выделение промежутка для значения куртозиса = [2 - bound_left/64, 2 + bound_right/64])
                    if i == 0:
                        n_aver = (frame_int & 0x3F00000000) >> 32
                        bound_left = (frame_int & 0x7FC000000000) >> (32 + 6)
                        bound_right = (frame_int & 0xFF800000000000) >> (32 + 6 + 9)
                    # Запись на первую позицию (с индексом 0) фрагмента спектра номера кадра frame_num
                    spectr_frame.append(frame_num)

                else:
                    spectr_val = (frame_int & 0x7FFFFFFFFFFFFF)
                    pp_good = (frame_int & 0xFF80000000000000) >> 55
                    spectr_frame.append(spectr_val)
                    pass

            spectr.append(spectr_frame)
            print(i, frame_num)
            i += 1

        pass

        spectr.pop(-1)
        N = len(spectr)
        n_frame_last = spectr[-1][0]
        rest = (n_frame_last + 1) % 2 ** (6 - n_aver)
        if rest:
            for k in range(rest):
                spectr.pop(-1)
        print(n_frame_last, spectr[-1][0])
    finally:
        f_in.close()
        pass

        spectr1 = convert_to_matrix(spectr, spectr[-1][0] + 1, n_aver)
    np.savetxt(file_name_out, spectr1, header=(str(n_aver) + '-n_aver ' + str(bound_left) + '-kurt'))

    return spectr1, n_aver


def extract_two_polar():
    file_name = file_name0 + '.bin'
    # file_name_out = file_name0 + '.txt'
    # *********** Если система работает с одной поляризацией ************
    if not file_name0.find('Ant2') == -1:
        antenna2_0 = 1
    else:
        antenna2_0 = 0
    # *******************************************************************
    i = 0
    k = 0
    spectr_left = []
    spectr_right = []
    attenuators = []
    frame = ' '

    try:
        if os.path.isfile(file_name) == 1:
            pass
        else:
            print('\n \t', file_name, ' not found!!!\n')

        f_in = open(file_name, 'rb')
        antenna = 0
        while frame:

            spectr_frame = []
            # Обработка кадра: выделение номера кадра, границ куртозиса, длины усреднения на ПЛИС
            # и 128-ми значений спектра в список spectr_frame на позиции [1:128]
            for k in range(130):
                frame = f_in.read(8)
                frame_int = int.from_bytes(frame, byteorder='little')
                if k == 0:
                    frame_num = frame_int & 0xFFFFFFF

                    # Выделение длины усреднения (количество усредняемых на ПЛИС отсчетов спектра = 2^n_aver)
                    # Выделение промежутка для значения куртозиса = [2 - bound_left/64, 2 + bound_right/64])
                    if i == 0:
                        n_aver = (frame_int & 0x3F00000000) >> 32
                        bound_left = (frame_int & 0x7FC000000000) >> (32 + 6)
                        bound_right = (frame_int & 0xFF800000000000) >> (32 + 6 + 9)
                    # Запись на первую позицию (с индексом 0) фрагмента спектра номера кадра frame_num
                    spectr_frame.append(frame_num)
                elif k == 1:
                    att_1 = frame_int & 0x3F
                    att_2 = (frame_int & 0xFC0) >> 6
                    att_3 = (frame_int & 0x3F000) >> 12
                    noise_gen_on = (frame_int & 0x40000) >> 18
                    antenna_before = antenna
                    antenna = (frame_int & 0x80000) >> 19
                    coupler = (frame_int & 0x100000) >> 20
                    band = (frame_int & 0x8000000000000000) >> 63
                    attenuators = [att_1, att_2, att_3]

                    pass

                else:
                    spectr_val = (frame_int & 0x7FFFFFFFFFFFFF)
                    pp_good = (frame_int & 0xFF80000000000000) >> 55
                    if pp_good / 256 < 0.78125:
                        spectr_val = 2
                    spectr_frame.append(spectr_val)
                    pass

            if antenna == 0 and (antenna_before - antenna == 0):
                spectr_left.append(spectr_frame)
            elif len(spectr_left) > 1 and ((antenna_before - antenna) != 0):
                spectr_left.pop(-1)

            if antenna == 1 and (antenna_before - antenna) == 0:
                spectr_right.append(spectr_frame)
            elif len(spectr_right) > 1 and ((antenna_before - antenna) != 0):
                spectr_right.pop(-1)
            print(i, frame_num, band)
            i += 1

        pass
        n_right = len(spectr_right)
        n_left = len(spectr_left)

        # В случае, если при работе с одной поляризацией ('Ant1' или 'Ant2') в переменную
        # antenna не записывается с какого входа берется сигнал (в любом случае antenna = 0),
        # то необходима следующая процедура перестановки значений переменных
        if n_right == 0 and antenna2_0 == 1:
            spectr_right = spectr_left
            spectr_left = []
            n_right = len(spectr_right)
            n_left = len(spectr_left)
        # **********************************************************************************
        if n_right > 1:
            spectr_right.pop(-1)
            n_frame_last = spectr_right[-1][0]
            rest = (n_frame_last + 1) % 2 ** (6 - n_aver)
            if rest:
                for k in range(rest):
                    spectr_right.pop(-1)
            print(n_frame_last, spectr_right[-1][0])
        if n_left > 1:
            spectr_left.pop(-1)
            n_frame_last = spectr_left[-1][0]
            rest = (n_frame_last + 1) % 2 ** (6 - n_aver)
            if rest:
                for k in range(rest):
                    spectr_left.pop(-1)
            print(n_frame_last, spectr_left[-1][0])
    finally:
        f_in.close()
        pass

    if n_left > 1:
        spectr1 = convert_to_matrix(spectr_left, spectr_left[-1][0] + 1, n_aver)
    else:
        spectr1 = []
    if n_right > 1:
        spectr2 = convert_to_matrix(spectr_right, spectr_right[-1][0] + 1, n_aver)
    else:
        spectr2 = []
    spectrum_extr = pd.Series([spectr1, spectr2])
    np.save(file_name0 + '_spectrum', spectrum_extr)
    head = np.array([n_aver, 6, 8])
    np.savetxt(file_name0 + '_head.txt', head)
    return spectrum_extr, n_aver


def extract_whole_band():
    file_name = Path(file_path_data, current_data_file + '.bin')
    # *********** Если система работает с одной поляризацией ************
    if not str(file_name).find('Ant2') == -1:
        antenna2_0 = 1
    else:
        antenna2_0 = 0
    # *******************************************************************
    i = 0
    spectrum_right_1 = []
    spectrum_left_1 = []
    spectrum_left_2 = []
    spectrum_right_2 = []
    attenuators = []
    frame = ' '

    try:
        if os.path.isfile(file_name) == 1:
            pass
        else:
            print('\n \t', file_name, ' not found!!!\n')
            return

        f_in = open(file_name, 'rb')
        antenna = 0
        while frame:
            spectr_frame = []
            # Обработка кадра: выделение номера кадра, границ куртозиса, длины усреднения на ПЛИС
            # и 128-ми значений спектра в список spectr_frame на позиции [2:129]
            for k in range(130):
                frame = f_in.read(8)
                frame_int = int.from_bytes(frame, byteorder='little')
                if k == 0:
                    frame_num = frame_int & 0xFFFFFFF

                    # Выделение длины усреднения (количество усредняемых на ПЛИС отсчетов спектра = 2^n_aver)
                    # Выделение промежутка для значения куртозиса = [2 - bound_left/64, 2 + bound_right/64])
                    if i == 0:
                        n_aver = (frame_int & 0x3F00000000) >> 32
                        bound_left = (frame_int & 0x7FC000000000) >> (32 + 6)
                        bound_right = (frame_int & 0xFF800000000000) >> (32 + 6 + 9)
                    # Запись на первую позицию (с индексом 0) фрагмента спектра номера кадра frame_num
                    spectr_frame.append(frame_num)
                elif k == 1:
                    att_1 = frame_int & 0x3F
                    att_1 = int((63 - att_1) / 2)
                    att_2 = (frame_int & 0xFC0) >> 6
                    att_2 = int((63 - att_2) / 2)
                    att_3 = (frame_int & 0x3F000) >> 12
                    att_3 = int((63 - att_3) / 2)
                    noise_gen_on = (frame_int & 0x40000) >> 18
                    antenna_before = antenna
                    antenna = (frame_int & 0x80000) >> 19
                    if antenna == 1:
                        pass
                    coupler = (frame_int & 0x100000) >> 20
                    band = (frame_int & 0x8000000000000000) >> 63
                    attenuators = [att_1, att_2, att_3]
                    if i == 10:
                        att01 = att_1
                        att02 = att_2
                        att03 = att_3
                    pass

                else:
                    spectrum_val = (frame_int & 0x7FFFFFFFFFFFFF) >> shift

                    # Отбросили "shift" младших разрядов двоичного представления или 3 разряда десятичного
                    # при "shift=10"
                    if band:
                        spectrum_val = int((spectrum_val * att_dict[att_3] * att_dict[att_1]))
                    else:
                        spectrum_val = int((spectrum_val * att_dict[att_3] * att_dict[att_2]))
                    if spectrum_val > 1000000000:
                        spectrum_val = 1000000000
                    pp_good = (frame_int & 0xFF80000000000000) >> 55
                    if pp_good / 256 < 0.1:
                        spectrum_val = 2
                    spectr_frame.append(spectrum_val)
                    pass

            if antenna == 0 and (antenna_before - antenna == 0):
                if band:
                    spectrum_left_2.append(spectr_frame)
                else:
                    spectrum_left_1.append(spectr_frame)
            if len(spectrum_left_1) > 1 and ((antenna_before - antenna) != 0):
                spectrum_left_1.pop(-1)
            if len(spectrum_left_2) > 1 and ((antenna_before - antenna) != 0):
                spectrum_left_2.pop(-1)

            if antenna == 1 and (antenna_before - antenna) == 0:
                if band:
                    spectrum_right_2.append(spectr_frame)
                else:
                    spectrum_right_1.append(spectr_frame)
            if len(spectrum_right_1) > 1 and ((antenna_before - antenna) != 0):
                spectrum_right_1.pop(-1)
            if len(spectrum_right_2) > 1 and ((antenna_before - antenna) != 0):
                spectrum_right_2.pop(-1)

            print(i, frame_num, band, attenuators)
            i += 1
        pass
        n_right1 = len(spectrum_right_1)
        n_left1 = len(spectrum_left_1)
        n_right2 = len(spectrum_right_2)
        n_left2 = len(spectrum_left_2)

        # В случае, если при работе с одной поляризацией ('Ant1' или 'Ant2') в переменную
        # antenna не записывается с какого входа берется сигнал (в любом случае antenna = 0),
        # то необходима следующая процедура перестановки значений переменных
        n_right = np.max([len(spectrum_right_1), len(spectrum_right_2)])
        if n_right == 0 and antenna2_0 == 1:
            spectrum_right_1 = spectrum_left_1
            spectrum_left_1 = []
            spectrum_right_2 = spectrum_left_2
            spectrum_left_2 = []
            n_right1 = len(spectrum_right_1)
            n_left1 = len(spectrum_left_1)
            n_right2 = len(spectrum_right_2)
            n_left2 = len(spectrum_left_2)

        # Приведение длины записи к величине кратной количеству частот
        if n_right1 > 1:
            spectrum_right_1 = cut_spectrum(spectrum_right_1, n_aver)
            # spectrum_right_1 = np.array(spectrum_right_1, dtype=np.int32)
            spectrum_right_1 = parts_to_numpy(spectrum_right_1, n_right1)
        if n_left1 > 1:
            spectrum_left_1 = cut_spectrum(spectrum_left_1, n_aver)
            spectrum_left_1 = np.array(spectrum_left_1, dtype=np.int32)
        if n_right2 > 1:
            spectrum_right_2 = cut_spectrum(spectrum_right_2, n_aver)
            # spectrum_right_2 = np.array(spectrum_right_2, dtype=np.int32)
            spectrum_right_2 = parts_to_numpy(spectrum_right_2, n_right2)
        if n_left2 > 1:
            spectrum_left_2 = cut_spectrum(spectrum_left_2, n_aver)
            spectrum_left_1 = np.array(spectrum_left_1, dtype=np.int32)
    finally:
        f_in.close()
        pass
    spectrum_extr = pd.Series([spectrum_left_1, spectrum_left_2, spectrum_right_1, spectrum_right_2])
    # head = [n_aver, shift, bound_left, bound_right, att01, att02, att03]
    band_size, polar, measure_kind = status_func(n_left1, n_left2, n_right1, n_right2)

    head = {'date': date,
            'measure_kind': measure_kind,    # Вид измерений: наблюдение Солнца, Луны, калибровка АЧХ
            'band_size': band_size,  # Параметр 'whole' означает работу в диапазоне 1-3 ГГц,
            # 'half_low' - диапазон 1-2, 'half_upper' - 2-3 ГГц
            'polar': polar,  # Принимает значения поляризаций: 'both', 'left', 'right'
            'cleaned': 'no',
            'n_aver': n_aver,
            'shift': shift,
            'kurtosis': bound_left,
            'att1': att01,
            'att2': att02,
            'att3': att03,
            'align_file_path': r'F:\Fast_Acquisition\Alignment\Align_coeff.bin',
            'align_coeff_pos': 5}
    return save_spectrum(spectrum_extr, head)


def parts_to_numpy(list_arr, len_list):
    n = int(len_list // 1e5)
    k = int(len_list % 1e5)
    numpy_arr = []
    for i in range(n + 1):
        if i == n:
            auxiliary = list_arr[int(i * 1e5):int(i * 1e5 + k)]
            auxiliary = np.array(auxiliary, dtype='int32')
        else:
            auxiliary = list_arr[int(i * 1e5):int((i + 1) * 1e5)]
            auxiliary = np.array(auxiliary, dtype='int32')
        l = np.size(numpy_arr)
        if l:
            numpy_arr = np.vstack([numpy_arr, auxiliary])
        else:
            numpy_arr = auxiliary

    return numpy_arr


def status_func(n_left1, n_left2, n_right1, n_right2):

    # Параметр 'whole' означает работу в диапазоне 1-3 ГГц,
    # 'half_low' - диапазон 1-2, 'half_upper' - 2-3 ГГц
    if (n_left1 > 1 and n_left2 > 1) or (n_right1 > 1 and n_right2 > 1):
        band_size = 'whole'
    if (n_left1 > 1 or n_right1 > 1) and (n_left2 == 0 and n_right2 == 0):
        band_size = 'half_low'
    if (n_left2 > 1 or n_right2 > 1) and (n_left1 == 0 and n_right1 == 0):
        band_size = 'half_upper'

    # polar Принамает значения поляризаций: 'both', 'left', 'right'
    if (n_left1 > 1 or n_left2 > 1) and (n_right1 == 0 or n_right2 == 0):
        polar = 'left'
    if (n_left1 == 0 or n_left2 == 0) and (n_right1 > 1 or n_right2 > 1):
        polar = 'right'
    if (n_left1 > 1 or n_left2 > 1) and (n_right1 > 1 or n_right2 > 1):
        polar = 'both'

    # Определение вида измерений: наблюдение Солнца, Луны, калибровка АЧХ
    measure_kind = ''
    file_name0 = str(Path(file_path_data, current_data_file))
    if file_name0.find('test') != -1:
        # l = file_name0.find('test')
        measure_kind = 'test'
    if file_name0.find('sun') != -1:
        measure_kind = 'Sun'
    if file_name0.find('moon') != -1:
        measure_kind = 'Moon'
    if file_name0.find('calibration') != -1:
        measure_kind = 'calibration'

    return band_size, polar, measure_kind


def save_spectrum(spectrum_extr, head):
    spectrum1 = spectrum_extr[0]
    spectrum2 = spectrum_extr[1]
    spectrum3 = spectrum_extr[2]
    spectrum4 = spectrum_extr[3]
    n_aver = head['n_aver']
    band_size = head['band_size']
    polar = head['polar']
    measure_kind = head['measure_kind']
    print(f'len_spectrum1 = {len(spectrum1)}, len_spectrum2 ={len(spectrum2)}, len_spectrum3 ={len(spectrum3)}, '
          f'len_spectrum4 ={len(spectrum4)}')
    if len(spectrum1) > 1:
        spectrum1_low = convert_to_matrix(spectrum1, spectrum1[-1][0] + 1, n_aver)
        pass
    else:
        spectrum1_low = []
    if len(spectrum2) > 1:
        spectrum1_high = convert_to_matrix(spectrum2, spectrum2[-1][0] + 1, n_aver)
    else:
        spectrum1_high = []
    if len(spectrum3) > 1:
        spectrum2_low = convert_to_matrix(spectrum3, spectrum3[-1][0] + 1, n_aver)
    else:
        spectrum2_low = []
    if len(spectrum4) > 1:
        spectrum2_high = convert_to_matrix(spectrum4, spectrum4[-1][0] + 1, n_aver)
    else:
        spectrum2_high = []
    spectrum_whole = pd.Series([spectrum1_low, spectrum1_high, spectrum2_low, spectrum2_high])
    np.save(Path(file_path_data, current_data_file + '_spectrum'), spectrum_whole)
    with open(Path(file_path_data, current_data_file + '_head.bin'), 'wb') as out:
        pickle.dump(head, out)
    jsn.dump(head, open(Path(file_path_data, current_data_file + '_head.txt'), "w"))
    # np.savetxt(file_name0 + '_head.bin', head)

    return spectrum_whole, n_aver, measure_kind, band_size, polar


def cut_spectrum(spectrum, n_aver):
    spectrum.pop(-1)
    n_frame_last = spectrum[-1][0]
    rest = (n_frame_last + 1) % 2 ** (6 - n_aver)
    if rest:
        for k in range(rest):
            spectrum.pop(-1)
    print(n_frame_last, spectrum[-1][0])
    return spectrum


def convert_to_matrix(S_total, counter, n_aver):
    """Функция принимает список списков S, который формируется в extract(file_name0) и превращает его в матрицу,
    строки которой - спектры с разрешением 7.8125/(2**(6-n_aver)) МГц, а столбцы - зависимость значения
    спектра на фиксированной частоте от времени. Разрешение по времени - 8192 мкс. Вместо пропущенных пакетов
    вставляет значение 2"""
    len_time = np.shape(S_total)[0]
    S = [[int(2)] * 128 for i in range(counter)]
    # S = [['NaN'] * 128 for i in range(counter)]
    for s in S_total:
        S[s[0]] = s[1:]
    aver_param_loc = 2 ** (6 - n_aver)
    n = 128 * aver_param_loc
    print(len(S))
    k = int(len(S) // n)
    print(f' n = {n}, k = {k}')
    parts = int(k // 100)
    rest = int(k % 100)
    s_agreed = []

    for i in range(parts +1):
        if i == parts:
            s_auxiliary = np.array(S[int(i * 100 * n): int((i * 100 + rest) * n)])
            s_ar = np.reshape(s_auxiliary, (-1, n))
        else:
            s_auxiliary = np.array(S[int(i * 100 * n): int((i + 1) * 100 * n)])
            s_ar = np.reshape(s_auxiliary, (-1, n))
        l = np.size(s_agreed)
        if l:
            s_agreed = np.vstack([s_agreed, s_ar])
        else:
            s_agreed = s_ar

    # s_agreed = np.array(S[0: k * n])
    print(f'len_s_agreed = {len(s_agreed)}, shape = {np.shape(s_agreed)}')
    # s_ar = np.reshape(s_agreed, (-1, n))
    return s_agreed


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
    return s_freq * (2 ** (shift - 10)), s_time * (2 ** (shift - 10))


def pic_title():
    title0 = file_name0[-19:-2]
    title1 = '  ' + title0[0:4] + '.' + title0[4:6] + '.' + title0[6:8] + \
             ' time=' + title0[9:11] + ':' + title0[11:13] + ' azimuth=' + title0[14:17]
    if not file_name0.find('sun') == -1:
        title2 = 'Sun intensity'
    elif not file_name0.find('crab') == -1:
        title2 = 'Crab intensity'
    elif not file_name0.find('calibr') == -1:
        title2 = 'Calibration'
        title0 = file_name0[-23:-2]
        title1 = '  ' + title0[0:4] + '.' + title0[4:6] + '.' + title0[6:8] + \
                 ' chanell att=' + title0[14:17] + ' source att=' + title0[18:21]
    elif not file_name0.find('test') == -1:
        title0 = file_name0[-24:-2]
        title2 = 'Test interference'
        title1 = '  ' + title0[0:4] + '.' + title0[4:6] + '.' + title0[6:8] + \
                 ' chanell att=' + title0[15:18] + ' source att=' + title0[19:22]
        pass
    else:
        title2 = []
    return title1, title2


def path_to_fig():
    """ Создает директорию для рисунков обрабатываемого наблюдения, если она до этого не была создана,
    название директории  совпадает с названием исходного файла данных наблюдения
    """
    if not os.path.isdir(Path(file_path_data, current_data_file)):
        os.mkdir(Path(file_path_data, current_data_file))
    return


def preparing_data():
    """ Функция в зависимости от вида данных (полная полоса 1-3 ГГц, половинная полоса 1-2 или 2-3 ГГц,
    с двумя поляризациями или одной) выдает данные для построения графиков"""

    # Для полосы 1-3 ГГц и двух возможных поляризаций выдает по два спектра (1-2 и 2-3 ГГц) для каждой поляризации.
    # Если поляризация не задействована, то соответствующие спектры - пустые. Спектр 1-2 ГГц - в обратном порядке
    if not (os.path.isfile(Path(file_path_data, current_data_file + '_spectrum.npy')) or os.path.isfile(Path(file_path_data, current_data_file + '_left1.npy'))):
        if num_of_polar == 2 and band_size_init == 'whole':
            spectrum, n_aver, measure_kind, band_size, polar = extract_whole_band()
        if num_of_polar == 2 and band_size_init == 'half':
            spectrum, n_aver = extract_two_polar()
    else:
        spectrum = np.load(Path(file_path_data, current_data_file + '_spectrum.npy'), allow_pickle=True)
        with open(Path(file_path_data, current_data_file + '_head.bin'), 'rb') as inp:
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

    return united_spec


# Чтение с диска, если спектры ранее извлекались,
# или извлечение спектров из исходных записей
spectr_extr_left1, spectr_extr_left2, spectr_extr_right1, spectr_extr_right2, n_aver, band_size, polar = \
    preparing_data()
aver_param = 2 ** (6 - n_aver)
with open(Path(file_path_data, current_data_file + '_head.bin'), 'rb') as inp:
    head = pickle.load(inp)

# Выравнивание спектров по результатам шумовых измерений АЧХ
if align == 'y':
    path_output = Path(folder_align_path, align_file_name)
    spectr_extr_left1, spectr_extr_left2, spectr_extr_right1, spectr_extr_right2 = \
        align_spectrum(spectr_extr_left1, spectr_extr_left2, spectr_extr_right1, spectr_extr_right2, head, path_output)

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

if noise_calibr == 'y':
    spectr_time = calibration(t_cal, spectr_time)

# # ***********************************************
# # ***        Графический вывод данных        ****
# # ***********************************************

# Динамическая маска (зависит от длины записи во времени)
t_spect = N_row * delta_t
time_spect_mask = [(lambda i: (t_spect * (i + 0.05)) // 7)(i) for i in range(7)]

# if band_size == 'whole':
#   freq_spect_mask = []

# Формирование спектров и сканов по маскам freq_spect_mask и time_spect_mask
spectr_freq, spectr_time = form_spectr_sp1(spectrum_extr, freq_spect_mask, time_spect_mask)
# np.save(file_name0 + '_spectr', spectr_time)
# Формирование строк-аргументов по времени и частоте и легенды
N_col = np.shape(spectrum_extr)[1]
if band_size_init == 'half':
    freq = np.linspace(1000 * (N_Nyq - 1) + 3.9063 / aver_param * kf, 1000 * N_Nyq - 3.9063 / aver_param * kf,
                       N_col // kf)
elif band_size_init == 'whole':
    freq = np.linspace(1000 + 3.9063 / aver_param * kf, 3000 - 3.9063 / aver_param * kf, N_col // kf)
timeS = np.linspace(0, delta_t * N_row, N_row // kt)

line_legend_time, line_legend_freq = line_legend(freq_spect_mask[:10])
info_txt = [('time resol = ' + str(delta_t * kt) + 'sec'),
            ('freq resol = ' + str(delta_f / aver_param * kf) + 'MHz'),
            ('polarisation ' + polar), 'align' + align]
path_to_fig()
fp.fig_plot(spectr_freq, 0, freq, 1, info_txt, Path(file_path_data, current_data_file), head, line_legend_time)
fp.fig_plot(spectr_time, 0, timeS, 0, info_txt, Path(file_path_data, current_data_file), head, line_legend_freq)
n_start_flame = int(t_start_flame // (delta_t * kt))
n_stop_flame = int(t_stop_flame // (delta_t * kt))

# *********************************************************
# ***        Вывод данных двумерный и трехмерный       ****
# *********************************************************

# Укрупнение  разрешения по частоте и времени для вывода в 2d и 3d
if graph_3d_perm == 'y' or contour_2d_perm == 'y':
    spectr_extr1 = spectr_construction(spectr_extr, kf, kt)
# Информация о временном и частотном резрешениях
info_txt = [('time resol = ' + str(delta_t * kt) + 'sec'),
            ('freq resol = ' + str(delta_f / aver_param * kf) + 'MHz'),
            ('polarisation ' + polar)]
path_to_fig()

if graph_3d_perm == 'y':
    graph_3d(freq, timeS, spectr_extr1, 0)
if contour_2d_perm == 'y':
    graph_contour_2d(freq, timeS, spectr_extr1, 0)

# if align == 'y':
#     align_coeff1 = align_func1(spectr_freq[1, :], 'y', aver_param)
#     spectr_extr = spectr_extr * align_coeff1

# if graph_3d_perm == 'y':
#     graph_3d(freq, timeS[n_start_flame:n_stop_flame], spectr_extr1[n_start_flame:n_stop_flame, :], 0)
# fp.fig_multi_axes(spectr_time[:10, n_start_flame:n_stop_flame], timeS[n_start_flame:n_stop_flame],
#                   info_txt, file_name0, freq_spect_mask[:10])
stop = datetime.now()
print('\n Total time = ', stop - start)
