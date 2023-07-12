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

    i = 0
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
        spectrum = pd.DataFrame(index=[0, 1, 2, 3], columns=['left', 'right'])
        for i in [0, 1, 2, 3]:
            for j in ['left', 'right']:
                spectrum[j].loc[i] = []
        while frame:
            # Обработка кадра: выделение номера кадра, границ куртозиса, длины усреднения на ПЛИС
            # и 128-ми значений спектра на позиции [1:129]

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
            # Выделение промежутка для значения куртозиса = [2 - bound_left/128, 2 + bound_right/128])

            n_aver = (frame_int & 0xFF00000000) >> 32
            bound_left = (frame_int & 0x1FF0000000000) >> (32 + 8)
            bound_right = (frame_int & 0xFF800000000000) >> (32 + 8 + 9)


            frame = f_in.read(8)
            frame_int = int.from_bytes(frame, byteorder='little')
            att_1 = frame_int & 0x3F
            att_1 = int((63 - att_1) / 2)
            att_2 = (frame_int & 0xFC0) >> 6
            att_2 = int((63 - att_2) / 2)
            att_3 = (frame_int & 0x3F000) >> 12
            att_3 = int((63 - att_3) / 2)
            antenna = (frame_int & 0x80000) >> 19
            noise_gen_on = (frame_int & 0x100000) >> 20
            band = (frame_int & 0xF00000000) >> 32
            attenuators = [att_1, att_2, att_3]
            if i == 10:
                att01 = att_1
                att02 = att_2
                att03 = att_3
            if antenna == 0:
                polarization = 'left'
            else:
                polarization = 'right'
            pass

            frame1 = f_in.read(8*128)
            if len(frame1) != 1024:
                break
            x = struct.unpack('ff'*128, frame1)
            verify_data = np.asarray(x).reshape(128, -1)
            spectrum_val = list(verify_data[:, 0])
            pp_good = verify_data[:, 1]
            # if band:
            #     spectrum_val = (spectrum_val * att_dict[att_3] * att_dict[att_1])
            # else:
            #     spectrum_val = (spectrum_val * att_dict[att_3] * att_dict[att_2])
            spectrum_val[pp_good < pp_good_bound] = 2

            # Запись на первую позицию фрагмента спектра номера кадра frame_num
            spectrum[polarization].loc[band].extend([frame_num])
            # Запись нумерованного фрагмента спектра
            spectrum[polarization].loc[band].extend(spectrum_val)
            pass

            if abs(frame_num_before - frame_num) > 1000:
                print('Прервывание обработки из-за сбоя определения номера кадра')
                break

            print(i, len(frame1), frame_num, band, attenuators)
            i += 1

            frame_num_before = frame_num
            # if att_1 == 31 & att_2 == 31 & att_3 == 31:
            #     break
        pass

        # Приведение длины записи к величине кратной количеству частот

        for i in [0, 1, 2, 3]:
            for j in ['left', 'right']:
                if len(spectrum[j].loc[i]) > 1:
                    spectrum[j].loc[i] = cut_spectrum(spectrum[j].loc[i], n_aver)
                    spectrum[j].loc[i] = np.array(spectrum[j].loc[i])

    finally:
        f_in.close()
    pass

    spectrum_len = pd.DataFrame(index=[0, 1, 2, 3], columns=['left', 'right'])
    for i in [0, 1, 2, 3]:
        for j in ['left', 'right']:
            spectrum_len[j].loc[i] = len(spectrum[j].loc[i])
    polar, measure_kind = status_func(spectrum_len)

    head = {'date': date,
            'measure_kind': measure_kind,    # Вид измерений: наблюдение Солнца, Луны, калибровка АЧХ
            'polar': polar,  # Принимает значения поляризаций: 'both', 'left', 'right'
            'cleaned': 'no',
            'n_aver': n_aver,
            'kurtosis': bound_left,
            'good_bound': pp_good_bound,
            'att1': att01,
            'att2': att02,
            'att3': att03,
            'align_file_path': r'F:\Fast_Acquisition\Alignment\Align_coeff.bin',
            'align_coeff_pos': 5}
    return save_spectrum(spectrum, head)


def status_func(_sp_len):

    # polar Принамает значения поляризаций: 'both', 'left', 'right'
    for j in [0, 1, 2, 3]:
        if _sp_len['left'].loc[j] > 1 and _sp_len['right'].loc[j] > 1:
            polar = 'both'
            break
        if _sp_len['left'].loc[j] <= 1 and _sp_len['right'].loc[j] > 1:
            polar = 'right'
            break
        if _sp_len['left'].loc[j] > 1 and _sp_len['right'].loc[j] <= 1:
            polar = 'left'
            break

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

    return polar, measure_kind


def save_spectrum(_spectrum, _head):

    n_aver = _head['n_aver']
    # print(f'len_spectrum1 = {len(spectrum1)}, len_spectrum2 ={len(spectrum2)}, len_spectrum3 ={len(spectrum3)}, '
    #       f'len_spectrum4 ={len(spectrum4)}')

    for j in [0, 1, 2, 3]:
        for i in ['left', 'right']:
            if len(_spectrum[j].loc[i]) > 1:
                _spectrum[j].loc[i] = convert_to_matrix(_spectrum[j].loc[i],
                                                        _spectrum[j].loc[i][-1][0] + 1, n_aver)

    # **** Создание папки для хранения конвертированных данных, если ее нет ****
    if not os.path.isdir(Path(converted_data_file_path)):
        os.mkdir(Path(converted_data_file_path))
    # **************************************************************************

    np.save(Path(converted_data_file_path, current_primary_file + '_spectrum'), _spectrum)
    with open(Path(converted_data_file_path, current_primary_file + '_head.bin'), 'wb') as out:
        pickle.dump(_head, out)
    jsn.dump(_head, open(Path(converted_data_file_path, current_primary_file + '_head.txt'), "w"))

    return print(f'Data in {current_primary_file} converted to numpy file and saved successfully')


def cut_spectrum(_spectrum, n_aver):
    _len = len(_spectrum)
    _row = _len // 129
    _sp = [[_spectrum[i + j * 129] for i in range(129)] for j in range(_row)]
    n_frame_last = _sp[-1][0]
    rest = (n_frame_last + 1) % 2 ** (6 - n_aver)
    if rest:
        for k in range(rest):
            _sp.pop(-1)
    print(n_frame_last, _sp[-1][0])
    return _sp


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
    if not os.path.isfile(path1):
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

    att_val = [i * 0.5 for i in range(64)]
    att_dict = {s: 10 ** (s / 10) for s in att_val}

    preparing_data()

    stop = datetime.now()
    print(f'Process duration = {stop - start} sec')
