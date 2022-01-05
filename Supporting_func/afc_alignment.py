import os
import numpy as np
import pandas as pd
import pickle
from path_to_Yandex_Disk import path_to_YaDisk


def noise_kp(file_name, n_nyq, diff='n'):
    """ Функция принимает файл шумового измерения мощностной АЧХ и возвращает
        нормированный к максимуму коэффициент передачи тракта по мощности"""
    if not file_name.find(r'E:\Measure') == -1:
        file_name0 = file_name
    else:
        file_name0 = file_name + '.txt'
    spectr_noise_out = np.loadtxt(file_name0)
    n_row, n_col = np.shape(spectr_noise_out)
    s = np.zeros(n_col)
    s1 = np.zeros(n_col)
    # Усреднение по времени
    for i in range(n_col):
        if diff == 'n':
            s[i] = np.sum(spectr_noise_out[:500, i]) / 500
        else:
            s[i] = np.sum(spectr_noise_out[:1600, i]) / 1600
            s1[i] = np.sum(spectr_noise_out[2000:3600, i]) / 1600
            s[i] -= s1[i]

    s_max = np.max(s)
    kp_norm = s / s_max
    return kp_norm, n_col, s


def align_func(calibr_file_name: object, diff: object = 'n', aver_param: object = 2) -> object:
    """ Функция возвращает коэффициенты, выравнивающие исходную АЧХ

    """
    # Исходные данные
    # N_Nyq = 3
    delta_f = 7.8125
    # aver_param = 2

    # Определение пути к файлу, где лежат результаты измерения АЧХ
    # по мощности с генератором шума на входе тракта, чтение параметра
    # усреднения n_aver_noise для этого измерения
    # head_path = path_to_YaDisk()
    head_path = 'E:\\BinData'
    # file_name0 = head_path + '\\Measure\\Fast_Acquisition\\Calibration\\' + calibr_file_name

    if not calibr_file_name.find('left') == -1:
        file_name0 = calibr_file_name
        n_nyq = int(calibr_file_name[-10])
        f_in1 = open(file_name0)
    elif not calibr_file_name.find('right') == -1:
        file_name0 = calibr_file_name
        n_nyq = int(calibr_file_name[-11])
        f_in1 = open(file_name0)
    else:
        file_name0 = head_path + '\\2020_12_16calibr\\' + calibr_file_name
        n_nyq = int(calibr_file_name[-1])
        f_in1 = open(file_name0 + '.txt')

    n_aver_noise = int((f_in1.readline())[2])
    aver_param_noise = 2 ** (6 - n_aver_noise)
    f_in1.close()

    kp_norm, n_col, spectrum_cal = noise_kp(file_name0, n_nyq, diff)

    # Исключение корректировки коэффициента усиления в зоне действия режекторных фильтров
    if n_nyq == 3:
        n1 = int((80 - delta_f / 2 / aver_param_noise) // (delta_f / aver_param_noise))
        n2 = int((230 - delta_f / 2 / aver_param_noise) // (delta_f / aver_param_noise))
    else:
        n1 = int((100 - delta_f / 2 / aver_param_noise) // (delta_f / aver_param_noise))
        n2 = int((212 - delta_f / 2 / aver_param_noise) // (delta_f / aver_param_noise))
    kp_norm[n1:n2] = 1

    band = int(aver_param_noise / aver_param)
    kp_band = [np.sum(kp_norm[band * i:band * (i + 1)]) / band for i in range(int(n_col / band))]
    align_coeff = [1 / a for a in kp_band]

    if n_nyq == 3:
        freq = np.linspace(1000 * (n_nyq - 1) + 3.9063 / aver_param, 1000 * n_nyq - 3.9063 / aver_param, n_col)
    else:
        freq = np.linspace(1000 * n_nyq - 3.9063 / aver_param, 1000 * (n_nyq - 1) + 3.9063 / aver_param, n_col)

    return align_coeff, freq * 1000000, spectrum_cal


def align_spectrum(spectrum1, spectrum2, spectrum3, spectrum4, head, path_calibration):
    """ Принимает спектры левой (spectrum1, spectrum2) и правой (spectrum3, spectrum4) поляризаций.
    При этом spectrum1 и spectrum3 относятся ко второй зоне Найквиста и имеют обратный порядок следования
    отсчетов по частоте. По пути path_calibration загружаем выравнивающие коэффициенты для спектров.
    Ant1 - левая поляризация, Ant2 - правая, 2 - вторая зана Найквиста и обратный порядок следования
    коэффициентоа по частоте, 3 - третья"""

    # file_name_calibr1: str = path_calibration + r'\Calibr_Ant1_2.txt'
    # align_coeff = np.loadtxt(path_calibration)
    with open(path_calibration, 'rb') as inp:
        calibration_frame_inp = pickle.load(inp)
    r = calibration_frame_inp.iloc[1]
    align_coeff = [r['spectrum_left1'], r['spectrum_left2'], r['spectrum_right1'], r['spectrum_right2']]
    # align_coeff2 =
    # align_coeff3 =
    # align_coeff4 =
    len_calibr = np.size(align_coeff[1])
    if np.size(spectrum1):
        l1 = np.shape(spectrum1)[1]
    else:
        l1 = 0
    if np.size(spectrum2):
        l2 = np.shape(spectrum2)[1]
    else:
        l2 = 0
    if np.size(spectrum3):
        l3 = np.shape(spectrum3)[1]
    else:
        l3 = 0
    if np.size(spectrum4):
        l4 = np.shape(spectrum4)[1]
    else:
        l4 = 0

    # Проверка совпадения разрешения по частоте принимаемых функцией спектров и калибровочных
    # коэффифиентов
    len_freq_spectrum = np.array([l1, l2, l3, l4]).max()
    align_coeff_matched = [[], [], [], []]
    if int(len_calibr) != int(len_freq_spectrum):
        s = int(len_calibr / len_freq_spectrum)
        print(f"Вам надо уменьшить разрешение по частоте в {s} раз")
        # consent = str(input('Продолжить выполнение без выравнивания АЧХ (y/n)?'))
        # if consent == 'y':
        #     return spectrum1, spectrum2, spectrum3, spectrum4
        # else:
        #     pass
        j = 0

        for obj in align_coeff:
            align_coeff_matched[j] = [np.sum(obj[i*s:i*s+s]) / s for i in range(int(len_calibr / s))]
            j += 1
    else:
        align_coeff_matched = align_coeff
        # file_name_calibr2: str = path_calibration + r'\Calibr_Ant1_3.txt'
        # file_name_calibr3: str = path_calibration + r'\Calibr_Ant2_2.txt'
        # file_name_calibr4: str = path_calibration + r'\Calibr_Ant2_3.txt'


    if l1:
        spectrum1 = spectrum1 * align_coeff_matched[0]
    if l2:
        spectrum2 = spectrum2 * align_coeff_matched[1]
    if l3:
        spectrum3 = spectrum3 * align_coeff_matched[2]
    if l4:
        spectrum4 = spectrum4 * align_coeff_matched[3]

    return spectrum1, spectrum2, spectrum3, spectrum4


if __name__ == '__main__':
    align_coeff = align_func(2)
    pass
