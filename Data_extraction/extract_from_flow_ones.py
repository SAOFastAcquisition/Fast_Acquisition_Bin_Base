import numpy as np
import os
import sys
import pandas as pd
import pickle
import gzip
import glob as gb
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
    # attenuators = []
    frame = ' '
    frame_num_before = 0
    noise_gen_on_before = 0
    flag_registration = 0
    flag_registration_before = 0
    ng_counter1 = 0
    ng_counter2 = 0
    file_ind = 0
    try:
        if os.path.isfile(file_name) == 1 or os.path.exists(f'{str(file_name)}.gz') == 1:
            pass
        else:
            print('\n \t', file_name, ' not found!!!\n')
            return

        if os.path.exists(f'{str(file_name)}.gz'):
            f_in = gzip.open(f'{str(file_name)}.gz', "rb")
        else:
            f_in = open(file_name, 'rb')

        antenna = 0
        frame_num = 0
        while frame:
            spectr_frame = []
            # Обработка кадра: выделение номера кадра, границ куртозиса, длины усреднения на ПЛИС
            # и 128-ми значений спектра в список spectr_frame на позиции [2:129]
            for k in range(130):
                frame = f_in.read(8)
                frame_int = int.from_bytes(frame, byteorder='little')
                if k == 0:
                    frame_num = frame_int & 0xFFFFFFF
                    ind = 0

                    # ********* Обработка сбоя приема "нулевого" байта с номером кадра *********
                    # Из f_in байты будут
                    # считываться до тех пор, пока число, прочитанное на месте в байте номера кадра, не
                    # будет отличаться от предыдущего на величину меньше заданной. Если на протяжении 10
                    # кадров такой байт не найдется (1300 байт), то обработка прерывается и выводятся
                    # прочитанные до момента сбоя данные
                    while abs(frame_num - frame_num_before) > 1000:
                        frame = f_in.read(8)
                        frame_int = int.from_bytes(frame, byteorder='little')
                        frame_num = frame_int & 0xFFFFFFF
                        ind += 1
                        # delta_frame = frame_num - frame_num_before
                        # print(ind, delta_frame)
                        if ind > 1300:
                            print('Прервывание обработки из-за сбоя определения номера кадра')
                            break
                    # ******************************************************************************

                    # Выделение длины усреднения (количество усредняемых на ПЛИС отсчетов спектра = 2^n_aver)
                    # Выделение промежутка для значения куртозиса = [2 - bound_left/64, 2 + bound_right/64])
                    if i == 0:
                        n_aver = (frame_int & 0x3F00000000) >> 32
                        bound_left = (frame_int & 0x7FC000000000) >> (32 + 6)
                        bound_right = (frame_int & 0xFF800000000000) >> (32 + 6 + 9)
                    # Запись на первую позицию (с индексом 0) фрагмента спектра номера кадра frame_num

                elif k == 1:
                    att_1 = frame_int & 0x3F
                    att_1 = int((63 - att_1) / 2)
                    att_2 = (frame_int & 0xFC0) >> 6
                    att_2 = int((63 - att_2) / 2)
                    att_3 = (frame_int & 0x3F000) >> 12
                    att_3 = int((63 - att_3) / 2)
                    antenna_before = antenna
                    antenna = (frame_int & 0x80000) >> 19
                    if flag_registration:
                        print('frame_num1 = ', frame_num, 'antenna = ', antenna)
                    if antenna == 1:
                        pass
                    noise_gen_on = (frame_int & 0x100000) >> 20
                    if noise_gen_on - noise_gen_on_before == 1:
                        ng_counter1 += 1
                        if (ng_counter1 % 2 == 0) & (ng_counter2 % 2 == 1):
                            flag_registration = 1
                        print('NG on, frame num: ', frame_num, ' flag_registration =', flag_registration)
                    if noise_gen_on - noise_gen_on_before == -1:
                        ng_counter2 += 1
                        if (ng_counter2 // 2 == 1) & (ng_counter1 // 2 == 1) & (ng_counter2 % 2 == 1) & \
                                (ng_counter1 % 2 == 1):
                            flag_registration = 0
                            break
                            pass
                        print('NG off, frame num: ', frame_num, ' flag_registration =', flag_registration)
                    # Запись на первую позицию (с индексом 0) фрагмента спектра номера кадра frame_num
                    if flag_registration == 1:
                        spectr_frame.append(frame_num)

                    band = (frame_int & 0x8000000000000000) >> 63
                    # attenuators = [att_1, att_2, att_3]
                    if i == 10:
                        att01 = att_1
                        att02 = att_2
                        att03 = att_3
                    pass

                else:
                    if flag_registration == 1:
                        spectrum_val = (frame_int & 0x7FFFFFFFFFFFFF) >> shift

                        # Отбросили "shift" младших разрядов двоичного представления или 3 разряда десятичного
                        # при "shift=10"
                        if band:
                            spectrum_val = int((spectrum_val * att_dict[att_3] * att_dict[att_1]))
                        else:
                            spectrum_val = int((spectrum_val * att_dict[att_3] * att_dict[att_2]))
                        # if spectrum_val > 1000000000:
                        #     spectrum_val = 1000000000
                        pp_good = (frame_int & 0xFF80000000000000) >> 55
                        if pp_good / 256 < pp_good_bound:
                            spectrum_val = 2
                        spectr_frame.append(spectrum_val)
                pass
                #  Завершилась обработка кадра (фрейма) длиной 130 байт
            if abs(frame_num_before - frame_num) > 1000:
                print('Прервывание обработки из-за сбоя определения номера кадра')
                break

            # Запись результатов обработки кадра в массивы данных по признаку поляризации и полосы
            if antenna == 0 and (antenna_before - antenna == 0) and flag_registration == 1:
                if band:
                    spectrum_left_2.append(spectr_frame)
                else:
                    spectrum_left_1.append(spectr_frame)

            # Удаление отсчетов в моменты переключения поляризации
            if len(spectrum_left_1) > 1 and ((antenna_before - antenna) != 0):
                spectrum_left_1 = pop_exception(spectrum_left_1)
            if len(spectrum_left_2) > 1 and ((antenna_before - antenna) != 0):
                spectrum_left_2 = pop_exception(spectrum_left_2)

            if antenna == 1 and (antenna_before - antenna) == 0 and flag_registration == 1:
                if band:
                    try:
                        spectrum_right_2.append(spectr_frame)
                    except AttributeError as _err:
                        print(_err)
                        pass
                else:
                    try:
                        spectrum_right_1.append(spectr_frame)
                    except AttributeError as _err:
                        print(_err)
                        pass

            # Удаление отсчетов в моменты переключения поляризации
            if len(spectrum_right_1) > 1 and ((antenna_before - antenna) != 0):
                spectrum_right_1 = pop_exception(spectrum_right_1)
            if len(spectrum_right_2) > 1 and ((antenna_before - antenna) != 0):
                spectrum_right_2 = pop_exception(spectrum_right_2)
            # print(i, frame_num, band, attenuators)
            i += 1

            frame_num_before = frame_num
            noise_gen_on_before = noise_gen_on
            if flag_registration_before - flag_registration == 1:
                spectrum_right_1, spectrum_left_1, spectrum_right_2, spectrum_left_2 = \
                    one_spectrum(spectrum_right_1, spectrum_left_1,
                                 spectrum_right_2, spectrum_left_2, antenna2_0, n_aver)
                n_right1 = len(spectrum_right_1)
                n_left1 = len(spectrum_left_1)
                n_right2 = len(spectrum_right_2)
                n_left2 = len(spectrum_left_2)
                print('Data in azimuth are formed')

                # Сдвиг начала отсчета выделенной записи в 0
                i_num = []
                if n_right1 > 1:
                    i_num_r1 = spectrum_right_1[0, 0]
                    i_num.append(i_num_r1)
                if n_right2 > 1:
                    i_num_r2 = spectrum_right_2[0, 0]
                    i_num.append(i_num_r2)
                if n_left1 > 1:
                    i_num_l1 = spectrum_left_1[0, 0]
                    i_num.append(i_num_l1)
                if n_left2 > 1:
                    i_num_l2 = spectrum_left_2[0, 0]
                    i_num.append(i_num_l2)

                frame_num_start = min(i_num)
                if n_right1 > 1:
                    spectrum_right_1[:, 0] = spectrum_right_1[:, 0] - frame_num_start
                if n_right2 > 1:
                    spectrum_right_2[:, 0] = spectrum_right_2[:, 0] - frame_num_start
                if n_left1 > 1:
                    spectrum_left_1[:, 0] = spectrum_left_1[:, 0] - frame_num_start
                if n_left2 > 1:
                    spectrum_left_2[:, 0] = spectrum_left_2[:, 0] - frame_num_start

                spectrum_extr = pd.Series([spectrum_left_1, spectrum_left_2, spectrum_right_1, spectrum_right_2])
                # head = [n_aver, shift, bound_left, bound_right, att01, att02, att03]
                band_size, polar, measure_kind = status_func(n_left1, n_left2, n_right1, n_right2)

                head = {'_date': date,
                        'measure_kind': measure_kind,  # Вид измерений: наблюдение Солнца, Луны, калибровка АЧХ
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
                # Сохраняем выделенную двумя импульсами ГШ запись
                save_spectrum(spectrum_extr, head, file_ind)
                print(' ')
            # Переходим к поиску и выделению следующей записи
            #     file_ind += 1
            #     spectrum_right_1 = []
            #     spectrum_left_1 = []
            #     spectrum_left_2 = []
            #     spectrum_right_2 = []
            #     pass
            flag_registration_before = flag_registration
        pass
    finally:
        f_in.close()
        pass
    return


def pop_exception(_list):
    try:
        _list.pop(-1)
    except AttributeError as _err:
        print(_err)
        pass
    return _list


def one_spectrum(_spectrum_right_1, _spectrum_left_1, _spectrum_right_2, _spectrum_left_2, antenna2_0, _n_aver):
    """
    Функция осуществляет первичое форматирование выделенных данных, при необходимости переставляя левую и правую
    поляризации, приводя их к одной длине.
    :param _spectrum_right_1:
    :param _spectrum_left_1:
    :param _spectrum_right_2:
    :param _spectrum_left_2:
    :param antenna2_0:
    :param _n_aver:
    :return:
    """
    n_right1 = len(_spectrum_right_1)
    n_left1 = len(_spectrum_left_1)
    n_right2 = len(_spectrum_right_2)
    n_left2 = len(_spectrum_left_2)

    # В случае, если при работе с одной поляризацией ('Ant1' или 'Ant2') в переменную
    # antenna не записывается с какого входа берется сигнал (в любом случае antenna = 0),
    # то необходима следующая процедура перестановки значений переменных
    n_right = np.max([len(_spectrum_right_1), len(_spectrum_right_2)])
    if n_right == 0 and antenna2_0 == 1:
        _spectrum_right_1 = _spectrum_left_1
        _spectrum_left_1 = []
        _spectrum_right_2 = _spectrum_left_2
        _spectrum_left_2 = []
        n_right1 = len(_spectrum_right_1)
        n_left1 = len(_spectrum_left_1)
        n_right2 = len(_spectrum_right_2)
        n_left2 = len(_spectrum_left_2)

    # Приведение длины записи к величине кратной количеству частот
    if n_right1 > 1:
        _spectrum_right_1 = cut_spectrum(_spectrum_right_1, _n_aver)
        # spectrum_right_1 = np.array(spectrum_right_1, dtype=np.int32)
        _spectrum_right_1 = parts_to_numpy(_spectrum_right_1, n_right1)
    if n_left1 > 1:
        _spectrum_left_1 = cut_spectrum(_spectrum_left_1, _n_aver)
        _spectrum_left_1 = np.array(_spectrum_left_1, dtype=np.int64)
    if n_right2 > 1:
        _spectrum_right_2 = cut_spectrum(_spectrum_right_2, _n_aver)
        # spectrum_right_2 = np.array(spectrum_right_2, dtype=np.int32)
        _spectrum_right_2 = parts_to_numpy(_spectrum_right_2, n_right2)
    if n_left2 > 1:
        _spectrum_left_2 = cut_spectrum(_spectrum_left_2, _n_aver)
        _spectrum_left_2 = np.array(_spectrum_left_2, dtype=np.int64)
    pass
    return _spectrum_right_1, _spectrum_left_1, _spectrum_right_2, _spectrum_left_2


def parts_to_numpy(list_arr, len_list):
    """ Функция превращает список в массив numpy.
    Разбивает файл на меньшие части и обрабатывает их поотдельности. По ходу завершения обработки частей
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
    if current_primary_dir.find('test') != -1:
        measure_kind = 'test'
    if current_primary_dir.find('sun') != -1:
        measure_kind = 'Sun'
    if current_primary_dir.find('moon') != -1:
        measure_kind = 'Moon'
    if current_primary_dir.find('calibration') != -1 or current_primary_dir.find('calibr') != -1:
        measure_kind = 'calibration'

    return band_size, polar, measure_kind


def save_spectrum(spectrum_extr, head, _k):
    spectrum1 = spectrum_extr[0]
    spectrum2 = spectrum_extr[1]
    spectrum3 = spectrum_extr[2]
    spectrum4 = spectrum_extr[3]
    n_aver = head['n_aver']
    # band_size = head['band_size']
    # polar = head['polar']
    # measure_kind = head['measure_kind']
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

    for i in range(parts + 1):
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
    path_zip = Path(converted_data_file_path, current_primary_file + '_spectrum.npy.gz')
    path2 = Path(converted_data_file_path, current_primary_file + '_left1.npy')
    if not (os.path.isfile(path1) or os.path.isfile(path_zip)):
        if num_of_polar == 2 and band_size_init == 'whole':
            extract_whole_band()
    else:
        print(f"Data with path {str(Path(primary_data_file_path, current_primary_file))} is converted to numpy file")

    return print(f'Data in {current_primary_file} converted to numpy file successfully')


if __name__ == '__main__':
    """
    Скрипт последовательно считывает первичную запись сеанса и по признаку включения калибровочных импульсов ГШ 
    выделяет записи в азимутах и сохраняет в каталоге 'Converted_data' в формате .numpy. Каждую запись 
    ограничивают два импульса ГШ - в начале и конце записи. 
    Для присвоения имени выделяемым файлам записей должен в current_primary_dir присутствовать файл 
    azimuth_file_name с именем вида '2022-06-22_01+28-28az.txt', содержащий колонку с заголовком "azimuth" с 
    перечислением наблюдаемых азимутов в виде csv-файла с разделителем ",". Азимут записывается в виде _хх+az, 
    где хх - порядковый номер наблюдения, az - значение азимута со знаком + или -. 
    """
    start = datetime.now()

    primary_data_dir = 'Primary_data'  # Каталог исходных данных (за определенный период, здесь - год)
    converted_data_dir = 'Converted_data'  # Каталог для записи результатов конвертации данных и заголовков
    data_treatment_dir = 'Data_treatment'  # Каталог для записи результатов обработки, рисунков

    az_dict = {'+24': '_01', '+20': '_02', '+16': '_03', '+12': '_04', '+08': '_05', '+04': '_06', '+00': '_07',
               '-04': '_08', '-08': '_09', '-12': '_10', '-16': '_11', '-20': '_12', '-24': '_13'}
    object = 'sun'
    current_primary_file = '2024-08-30_01+24'
    current_primary_dir = current_primary_file[0:4] + '_' + current_primary_file[5:7] + '_' + \
                          current_primary_file[8:10] + object
    current_data_dir = current_primary_file[0:4]  # Каталог всех данных (первичных, вторичных) за год
    azimuth_file_name = current_primary_file + 'az.txt'
    # Переопределение каталога всех данных при калибровочных и тестовых наблюдениях
    if current_primary_dir.find('test') != -1 or current_primary_dir.find('calibration') != -1 \
            or current_primary_dir.find('calibr') != -1:
        current_data_dir = '2022/Test_and_calibration'

    current_primary_path = Path(primary_data_dir, current_primary_dir)
    current_converted_dir = current_primary_dir + '_conv'
    current_converted_path = Path(converted_data_dir, current_converted_dir)

    primary_data_file_path, head_path = path_to_data(current_data_dir, current_primary_path)

    converted_data_file_path, head_path = path_to_data(current_data_dir, current_converted_path)
    path_to_azimuth = Path(primary_data_file_path, azimuth_file_name)

    # align_file_name = 'Align_coeff.bin'  # Имя файла с текущими коэффициентами выравнивания АЧХ
    # folder_align_path = Path(head_path, 'Alignment')

    date = current_primary_file[0:10]

    # Переименование исходных файлов по шаблону скрипта
    paths = gb.glob(str(Path(primary_data_file_path, "*.bin.gz")))  # Set Pattern in glob() function
    if paths:
        for s in paths:
            s_num = s[-10:-7]
            s_new = str(Path(primary_data_file_path, date + az_dict[s_num] + s[-10::]))
            os.rename(s, s_new)
    paths = gb.glob(str(Path(primary_data_file_path, "*.bin.gz")))

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
    # azimuth_list = pd.read_csv(path_to_azimuth, delimiter=',')          # Перечень номеров записей с азимутими
    # output_filename_list = [_date + s for s in azimuth_list['azimuth']]
    if paths:
        for s in paths:
            current_primary_file = s[-23:-7]
            preparing_data()
    else:
        paths = gb.glob(str(Path(primary_data_file_path, "*.bin")))
        for s in paths:
            current_primary_file = s[-20:-4]
            preparing_data()
    stop = datetime.now()
    print(f'Process duration = {stop - start} sec')
