import pandas as pd
import os
import pickle
import math
import numpy as np
from matplotlib import pyplot as plt
import sys
from pathlib import Path
from Supporting_func import path_to_data


def mean_and_std(_t_start, _t_stop, _data):
    """
    Функция принимает двумерный массиив данных, начальный и конечный момент их считывания по индексу "0"  и
    вычисляет средние значения и средние стстистические отклонения от среднего по индексу "0". При этом выбросы
    большие, чем 0.1 * среднее, заменяет на знечение среднего. Возвращает отношение среднего отклонения к среднему
    (обратная величина к радиометрическому выигрышу) отклонение от среднего и среднее. Время дискретизации по оси "0"
    задается как delta_t в главной части.
    :param _t_start: Начало фрагмента по индексу "0"
    :param _t_stop: Окончание фрагмента по индексу "0"
    :param _data: Двумерный массив данных
    :return: Отношение среднего отклонения к среднему (обратная величина к радиометрическому выигрышу),
    отклонение от среднего, среднее
    """
    n_start = int(_t_start // delta_t)
    n_stop = int(_t_stop // delta_t)
    _data_clipped = _data[n_start:n_stop, :]
    _mean = np.mean(_data_clipped, axis=0)
    _shape = np.shape(_data_clipped)
    for i in range(_shape[1]):
        _data_clipped[abs(_data_clipped[:, i] - _mean[i]) > 0.5 * _mean[i], i] = _mean[i]
    _mean = np.mean(_data_clipped, axis=0)
    _std = np.std(_data[n_start:n_stop, :], axis=0)
    _dev = _std / _mean
    return _dev, _std, _mean


def paths(_current_primary_file, _current_primary_dir):
    current_data_dir = '2022'
    current_converted_dir = _current_primary_dir + '_conv'
    current_treatment_dir = _current_primary_dir + '_treat'
    converted_data_dir = 'Converted_data'  # Каталог для записи результатов конвертации данных и заголовков
    data_treatment_dir = 'Data_treatment'  # Каталог для записи результатов обработки, рисунков
    current_converted_path = Path(converted_data_dir, current_converted_dir)
    current_treatment_path = Path(data_treatment_dir, current_treatment_dir)
    converted_data_file_path, head_path = path_to_data(current_data_dir, current_converted_path)
    data_treatment_file_path, head_path = path_to_data(current_data_dir, current_treatment_path)

    _path_mesh1 = Path(data_treatment_file_path, 'Meshed_Spectrum', _current_primary_file + '_meshed' + '.npy')
    _path_head = Path(converted_data_file_path, _current_primary_file + '_head.bin')

    return _path_mesh1, _path_head


def noise_gen_test(_spectrum):
    """
    Принимает матрицу спектров с фиксированным шагом по частоте (ось "1") и времени (ось "0") при включении
    на входе радиометра по 15 сек. внешнего генератора шума 6000К, согласованной нагрузки и согласованной нагрузки
    вместе с встроенныи ГШ
    :param _spectrum:
    :return: относительное стандартное отклонение от среднего, стандартное отклонение от среднего и среднее
            (все в зависимости от частоты)
    """
    # Среднее, среднестатистическое отклонение для радиометра с ГШ 6000К на входе
    _t_start_extNG = 5
    _t_stop_extNG = 10
    dev0, b0, a0 = mean_and_std(_t_start_extNG, _t_stop_extNG, _spectrum)

    # Среднее, среднестатистическое отклонение для радиометра с согласованной нагрузкой на входе
    _t_start_ML = 20
    _t_stop_ML = 25
    dev1, b1, a1 = mean_and_std(_t_start_ML, _t_stop_ML, _spectrum)

    # Среднее, среднестатистическое отклонение для радиометра с внутренним ГШ на входе
    _t_start_intNG = 35
    _t_stop_intNG = 40
    dev2, b2, a2 = mean_and_std(_t_start_intNG, _t_stop_intNG, _spectrum)

    _dev = [dev0, dev1, dev2]
    _mean = [a0, a1, a2]
    _std = [b0, b1, b2]

    one = [1 / radiometric_advantage] * m
    arg = [1000 + delta_f * i for i in range(m)]

    # ******************** Рисунок *************************
    fig = plt.figure(figsize=(12, 8))
    axes = fig.add_subplot()
    axes.grid(b=True, which='major', color='#666666', linestyle='-')
    # Show the minor grid lines with very faint and almost transparent grey lines
    axes.minorticks_on()
    axes.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)

    axes.set_title('STD divided by MEAN', fontsize=18)
    axes.set_yscale('log')
    # axes.set_xscale('log')
    axes.set_xlabel('Freq, MHz', fontsize=18)
    axes.set_ylabel('STD divided by MEAN', fontsize=20)
    axes.set_ylim([0.1 / radiometric_advantage, 1])
    line_legend = ['External NG', 'Matched Load', 'Intrinsic NG']
    i = 0
    for s in _dev:
        axes.plot(arg, s, label=line_legend[i])
        i += 1
    axes.plot(arg, one, label='1/radiometric_advantage')
    axes.legend()
    plt.show()

    return _dev, _std, _mean


def low_noise_spectra_base(_spectrum, _head, _freq_mask, _arg, _current_primary_file):
    # Если файла хранениия коэффициентов не существует, то создаем его, если существует - загружаем
    # Формируем названия столбцов фрейма для хранения низкочастолтных спектров
    columns_names = ['file_name', 'att1', 'att2', 'att3', 'kurtosis', 'polar', 'config', 'freq_mask',
                     'arg', 'LN_spectra']
    config_dict = {'1': 'ADC only',
                   '2': 'ADC + Micran_ampl',
                   '3': 'ADC + Micran_ampl + preampl2',
                   '4': 'ADC + Micran_ampl + preampl2 + preampl1'}
    # Путь к папке, где будет находиться база со спектрами "LN_spectra_base.bin"
    path_to_ln_spectra, head_path = path_to_data('2022', 'Data_treatment')
    ln_spectra_file_name = 'LN_spectra_base.bin'

    with open(Path(path_to_ln_spectra, ln_spectra_file_name), 'rb') as inp:
        low_noise_spectra = pickle.load(inp)

    m = np.size(_freq_mask)
    for i in range(m):
        low_noise_spectrum = {'file_name': _current_primary_file, 'att1': _head['att1'], 'att2': _head['att2'],
                              'att3': _head['att3'], 'kurtosis': _head['kurtosis'], 'polar': _head['polar'],
                              'config': '1', 'freq_mask': _freq_mask[i], 'arg': _arg, 'LN_spectra': _spectrum[i, :]
                              }
        low_noise_ser = pd.Series(low_noise_spectrum)
        idx = None
        path1 = Path(path_to_ln_spectra, ln_spectra_file_name)
        if not os.path.isfile(path1):
            low_noise_spectra = pd.DataFrame(columns=columns_names)
            low_noise_spectra = low_noise_spectra.append(low_noise_ser, ignore_index=True)
            # low_noise_spectra = low_noise_spectra.drop(1)
        pass
        if low_noise_spectra.shape[0] == 0:
            low_noise_spectra = low_noise_spectra.append(low_noise_ser, ignore_index=True)
        else:
            a11 = _freq_mask
            a12 = low_noise_spectra.freq_mask
            a2 = low_noise_spectra.freq_mask
            idx = low_noise_spectra.loc[(low_noise_spectra.freq_mask == _freq_mask[i])
                                        & (low_noise_spectra.file_name == _current_primary_file)].index
        if idx == None:
            len_idx = 0
        else:
            len_idx = len(idx)
        if len_idx:
            print('Such low noise spectra is exist')
        else:
            low_noise_spectra = low_noise_spectra.append(low_noise_ser, ignore_index=True)
    # low_noise_spectra = low_noise_spectra.drop([1,2,3,4,5,6,7])
    with open(Path(path_to_ln_spectra, ln_spectra_file_name), 'wb') as out:
        pickle.dump(low_noise_spectra, out)

    pass
    return


if __name__ == '__main__':
    current_primary_file = '2022-01-27_11'
    current_primary_dir = '2022_01_27calibr'
    path_mesh, path_head = paths(current_primary_file, current_primary_dir)

    delta_t = 8.3886e-3
    delta_f = 7.8125
    radiometric_advantage = math.sqrt(delta_t * delta_f * 1e6)

    spectrum_mesh = spectrum = np.load(path_mesh, allow_pickle=True)
    l, m = np.shape(spectrum_mesh)
    with open(path_head, 'rb') as inp:
        head = pickle.load(inp)
    freq_mask = [1250, 1550, 2350, 2800]
    arg = [1000 + delta_f * i for i in range(m)]
    low_noise_spectra_base(spectrum_mesh, head, freq_mask, arg, current_primary_file)

    dev, std, mean = noise_gen_test(spectrum_mesh)

pass
