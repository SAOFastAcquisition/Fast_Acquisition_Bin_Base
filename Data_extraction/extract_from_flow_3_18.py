import numpy as np
import os
import sys
import struct
import pandas as pd
import pickle
import json as jsn
from datetime import datetime
from pathlib import Path
from Supporting_func import path_to_data

current_dir = Path.cwd()
home_dir = Path.home()

sys.path.insert(0, Path(current_dir, 'Supporting_func'))
sys.path.insert(0, Path(current_dir, 'Interface'))


def extract_whole_band():
    file_name = Path(primary_data_file_path, current_primary_file + '.bin')
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
    frame_num_before = 0

    try:
        if os.path.isfile(file_name) == 1:
            pass
        else:
            print('\n \t', file_name, ' not found!!!\n')
            return

        f_in = open(file_name, 'rb')
        antenna = 0
        frame_num = 0
        while frame:
            spectr_frame = []
            # Обработка кадра: выделение номера кадра, границ куртозиса, длины усреднения на ПЛИС
            # и 128-ми значений спектра в список spectr_frame на позиции [2:129]


            frame_int = int.from_bytes(f_in.read(8), byteorder='little')
            frame_num = frame_int & 0xFFFFFFF
            ind = 0

            # ********* Обработка сбоя приема "нулевого" байта с номером кадра *********
            # Из f_in байты будут по 8
            # считываться до тех пор, пока число, прочитанное на месте в байте номера кадра, не
            # будет отличаться от предыдущего на величину меньше заданной. Если на протяжении 10
            # кадров такой байт не найдется (1300*8 байт), то обработка прерывается и выводятся
            # прочитанные до момента сбоя данные
            while abs(frame_num - frame_num_before) > 1000:
                frame_int = int.from_bytes(f_in.read(8), byteorder='little')
                frame_num = frame_int & 0xFFFFFFF
                ind += 1
                if ind > 1300:
                    print('Прервывание обработки из-за сбоя определения номера кадра')
                    break
            # ******************************************************************************

            # Выделение длины усреднения (количество усредняемых на ПЛИС отсчетов спектра = 2^n_aver)
            # Выделение промежутка для значения куртозиса = [2 - bound_left/64, 2 + bound_right/64])

            n_aver = (frame_int & 0xFF00000000) >> 32
            bound_left = (frame_int & 0x1FF0000000000) >> (32 + 8)
            bound_right = (frame_int & 0xFF800000000000) >> (32 + 8 + 9)
            # Запись на первую позицию (с индексом 0) фрагмента спектра номера кадра frame_num
            spectr_frame.append(frame_num)

            frame = f_in.read(8)
            frame_int = int.from_bytes(frame, byteorder='little')
            att_1 = frame_int & 0x3F
            att_1 = int((63 - att_1) / 2)
            att_2 = (frame_int & 0xFC0) >> 6
            att_2 = int((63 - att_2) / 2)
            att_3 = (frame_int & 0x3F000) >> 12
            att_3 = int((63 - att_3) / 2)
            antenna_before = antenna
            antenna = (frame_int & 0x80000) >> 19
            if antenna == 1:
                pass
            noise_gen_on = (frame_int & 0x100000) >> 20
            band = (frame_int & 0xF00000000) >> 32
            attenuators = [att_1, att_2, att_3]
            if i == 10:
                att01 = att_1
                att02 = att_2
                att03 = att_3
            pass

            frame1 = f_in.read(8*128)
            x = struct.unpack('ff'*128, frame1)
            verify_data = np.asarray(x).reshape(128, -1)
            spectrum_val = list(verify_data[:, 0])
            pp_good = verify_data[:, 1]
            # if band:
            #     spectrum_val = (spectrum_val * att_dict[att_3] * att_dict[att_1])
            # else:
            #     spectrum_val = (spectrum_val * att_dict[att_3] * att_dict[att_2])
            # if spectrum_val > 1000000000:
            #     spectrum_val = 1000000000
            #  pp_good = (frame_int & 0xFF800000) >> 23
            # if pp_good < pp_good_bound:
            #     spectrum_val = 2
            spectr_frame.extend(list(spectrum_val))
            pass

            if abs(frame_num_before - frame_num) > 1000:
                print('Прервывание обработки из-за сбоя определения номера кадра')
                break
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

            frame_num_before = frame_num
            # if att_1 == 31 & att_2 == 31 & att_3 == 31:
            #     break
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
            spectrum_left_1 = np.array(spectrum_left_1, dtype=np.int64)
        if n_right2 > 1:
            spectrum_right_2 = cut_spectrum(spectrum_right_2, n_aver)
            # spectrum_right_2 = np.array(spectrum_right_2, dtype=np.int32)
            spectrum_right_2 = parts_to_numpy(spectrum_right_2, n_right2)
        if n_left2 > 1:
            spectrum_left_2 = cut_spectrum(spectrum_left_2, n_aver)
            spectrum_left_2 = np.array(spectrum_left_2, dtype=np.int64)
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
            'good_bound': pp_good_bound,
            'att1': att01,
            'att2': att02,
            'att3': att03,
            'align_file_path': r'F:\Fast_Acquisition\Alignment\Align_coeff.bin',
            'align_coeff_pos': 5}
    return save_spectrum(spectrum_extr, head)


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


def status_func(n_left1, n_left2, n_right1, n_right2):

    # Параметр 'whole' означает работу в диапазоне 1-3 ГГц,
    # 'half_low' - диапазон 1-2, 'half_upper' - 2-3 ГГц
    if (n_left1 > 1 and n_left2 > 1) or (n_right1 > 1 and n_right2 > 1):
        band_size = 'whole'
    if (n_left1 > 1 or n_right1 > 1) and (n_left2 <= 1 and n_right2 <= 1):
        band_size = 'half_low'
    if (n_left2 > 1 or n_right2 > 1) and (n_left1 <= 1 and n_right1 <= 1):
        band_size = 'half_upper'

    # polar Принамает значения поляризаций: 'both', 'left', 'right'
    if (n_left1 > 1 or n_left2 > 1) and (n_right1 <= 1 or n_right2 <= 1):
        polar = 'left'
    if (n_left1 <= 1 or n_left2 <= 1) and (n_right1 > 1 or n_right2 > 1):
        polar = 'right'
    if (n_left1 > 1 or n_left2 > 1) and (n_right1 > 1 or n_right2 > 1):
        polar = 'both'

    # Определение вида измерений: наблюдение Солнца, Луны, калибровка АЧХ
    measure_kind = ''
    # file_name0 = str(Path(converted_dir_path, current_data_file))
    if current_primary_dir.find('test') != -1:
        # l = file_name0.find('test')
        measure_kind = 'test'
    if current_primary_dir.find('sun') != -1:
        measure_kind = 'Sun'
    if current_primary_dir.find('moon') != -1:
        measure_kind = 'Moon'
    if current_primary_dir.find('calibration') != -1 or current_primary_dir.find('calibr') != -1:
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

    # **** Создание папки для хранения конвертированных данных, если ее нет ****
    if not os.path.isdir(Path(converted_data_file_path)):
        os.mkdir(Path(converted_data_file_path))
    # **************************************************************************

    np.save(Path(converted_data_file_path, current_primary_file + '_spectrum'), spectrum_whole)
    with open(Path(converted_data_file_path, current_primary_file + '_head.bin'), 'wb') as out:
        pickle.dump(head, out)
    jsn.dump(head, open(Path(converted_data_file_path, current_primary_file + '_head.txt'), "w"))

    return print(f'Data in {current_primary_file} converted to numpy file and saved successfully')


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


def preparing_data():
    """ Функция в зависимости от вида данных (полная полоса 1-3 ГГц, половинная полоса 1-2 или 2-3 ГГц,
    с двумя поляризациями или одной) выдает данные для построения графиков"""

    # Для полосы 1-3 ГГц и двух возможных поляризаций выдает по два спектра (1-2 и 2-3 ГГц) для каждой поляризации.
    # Если поляризация не задействована, то соответствующие спектры - пустые. Спектр 1-2 ГГц - в обратном порядке
    path1 = Path(converted_data_file_path, current_primary_file + '_spectrum.npy')
    path2 = Path(converted_data_file_path, current_primary_file + '_left1.npy')
    if not (os.path.isfile(path1) or
            os.path.isfile(path2)):
        if num_of_polar == 2 and band_size_init == 'whole':
            extract_whole_band()
    else:
        print(f"Data with path {str(Path(primary_data_file_path, current_primary_file))} is converted to numpy file")

    return print(f'Data in {current_primary_file} converted to numpy file successfully')


if __name__ == '__main__':

    start = datetime.now()

    current_data_dir = '2023'
    primary_data_dir = 'Primary_data_3_18'           # Каталог исходных данных (за определенный период, здесь - год)
    converted_data_dir = 'Converted_data_3_18'       # Каталог для записи результатов конвертации данных и заголовков
    data_treatment_dir = 'Data_treatment3_18'       # Каталог для записи результатов обработки, рисунков

    current_primary_dir = '2023_06_25test'
    current_primary_file = '2023-06-25_04'
    # Переопределение каталога всех данных при калибровочных и тестовых наблюдениях
    # if current_primary_dir.find('test') != -1 or current_primary_dir.find('calibration') != -1 \
    #         or current_primary_dir.find('calibr') != -1:
    #     current_data_dir = '2022/Test_and_calibration'
    current_primary_path = Path(primary_data_dir, current_primary_dir)
    current_converted_dir = current_primary_dir + '_conv'
    current_converted_path = Path(converted_data_dir, current_converted_dir)

    primary_data_file_path, head_path = path_to_data(current_data_dir, current_primary_path)
    converted_data_file_path, head_path = path_to_data(current_data_dir, current_converted_path)

    align_file_name = 'Align_coeff.bin'         # Имя файла с текущими коэффициентами выравнивания АЧХ
    folder_align_path = Path(head_path, 'Alignment')

    date = current_primary_file[0:10]

    # !!!! ******************************************* !!!!
    # ****** Блок исходных параметров для обработки *******
    band_size_init = 'whole'
    num_of_polar = 2
    pp_good_bound = 0.5
    shift = 0
    # band_size = 'whole'   Параметр 'whole' означает работу в диапазоне 1-3 ГГц, 'half' - диапазон 1-2 или 2-3 ГГц
    # polar = 'both'        Принимает значения поляризаций: 'both', 'left', 'right'

    att_val = [i * 0.5 for i in range(64)]
    att_dict = {s: 10 ** (s / 10) for s in att_val}

    preparing_data()

    stop = datetime.now()
    print(f'Process duration = {stop - start} sec')
