import os
import numpy as np
from path_to_Yandex_Disk import path_to_YaDisk

""" Приложение принимает файл со spectrum1, spectrum2 - левая поляризация, полосы 1-2 и 2-3 ГГц,
соответственно, spectrum3, spectrum4 - правая поляризация, полосы 1-2 и 2-3 ГГц, соответственно. 
возвращает файл со спектрами, выровненными по АЧХ приемника, полученной в результате измерений 
с калибровочным  шумовым генератором. Местоположение выравнивающих коэффициентов и исходных данных 
измерений с ГШ привязываются ко входным данным"""



def align_spectrum(spectrum1, spectrum2, spectrum3, spectrum4, path_calibration):
    """ Принимает спектры левой (spectrum1, spectrum2) и правой (spectrum3, spectrum4) поляризаций.
    При этом spectrum1 и spectrum3 относятся ко второй зоне Найквиста и имеют обратный порядок следования
    отсчетов по частоте. По пути path_calibration загружаем выравнивающие коэффициенты для спектров.
    Ant1 - левая поляризация, Ant2 - правая, 2 - вторая зана Найквиста и обратный порядок следования
    коэффициентоа по частоте, 3 - третья"""

    file_name_calibr1: str = path_calibration + r'\Calibr_Ant1_2.txt'
    align_coeff1 = np.loadtxt(file_name_calibr1)
    len_calibr = np.size(align_coeff1)
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
    if int(len_calibr) != int(len_freq_spectrum):
        s = len_freq_spectrum / len_calibr

        print(f"Вам надо уменьшить разрешение по частоте в {s} раз")
        consent = str(input('Продолжить выполнение без выравнивания АЧХ (y/n)?'))
        if consent == 'y':
            return spectrum1, spectrum2, spectrum3, spectrum4
        else:
            pass

    else:
        file_name_calibr2: str = path_calibration + r'\Calibr_Ant1_3.txt'
        file_name_calibr3: str = path_calibration + r'\Calibr_Ant2_2.txt'
        file_name_calibr4: str = path_calibration + r'\Calibr_Ant2_3.txt'
        align_coeff2 = np.loadtxt(file_name_calibr2)
        align_coeff3 = np.loadtxt(file_name_calibr3)
        align_coeff4 = np.loadtxt(file_name_calibr4)

    if l1:
        spectrum1 = spectrum1 * align_coeff1
    if l2:
        spectrum2 = spectrum2 * align_coeff2
    if l3:
        spectrum3 = spectrum3 * align_coeff3
    if l4:
        spectrum4 = spectrum4 * align_coeff4

    return spectrum1, spectrum2, spectrum3, spectrum4


if __name__ == '__main__':
    align_coeff = align_func(2)
    pass
