import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy import interpolate
from Supporting_func.afc_alignment import align_func as af
sys.path.insert(0, r'E:/rep1/Supporting_func')


def noise_gen_spectrum(x1):
    # ********************* Учет спектра генератора шума *************************
    # ***** Превращение текстовых файлов данных от анализатора R&S в спектр ******
    # ********* генератора шума в диапазоне 1 ... 3 ГГц noise_gen от freq ********

    #  Путь к файлам мощности шумов (генератор шума + собственные шумы R&S
    # анализатора спектра) и (собственные шумы анализатора R&S)
    file_name_noise: str = path_input + r'\2020_12_24_09.txt'
    file_name_noise1: str = path_input + r'\2020_12_24_07.txt'
    file_name_noise_int: str = path_input + r'\2020_12_24_08.txt'

    with open(file_name_noise, "r") as file:
        contents_noise_gen = file.readlines()
    file.close()

    with open(file_name_noise1, "r") as file:
        contents_noise_gen1 = file.readlines()
    file.close()

    with open(file_name_noise_int, "r") as file:
        contents_noise_int = file.readlines()
    file.close()

    freq = np.zeros(617)
    noise_gen = np.zeros(617)  #
    noise_gen1 = np.zeros(617)  #
    noise_int = np.zeros(617)  #
    data_gen = contents_noise_gen[30:647]      # Спектр ГШ с ослаблением сигнала ГШ 3 дБ
    data_gen1 = contents_noise_gen1[30:647]    # Спектр ГШ с ослаблением сигнала ГШ 10 дБ
    data_int = contents_noise_int[30:647]      # Спектр собственных шумов R&S
    i = 0
    for s in data_gen:
        data1 = s.split(';')
        freq[i] = float(data1[0])
        noise_gen[i] = float(data1[1])
        i += 1

    i = 0
    for s in data_gen1:
        data1 = s.split(';')
        noise_gen1[i] = float(data1[1])
        i += 1

    i = 0
    for s in data_int:
        data1 = s.split(';')
        noise_int[i] = float(data1[1])
        i += 1

    noise_int = np.power(10, noise_int / 10)
    noise_gen = np.power(10, noise_gen / 10)
    noise_gen1 = np.power(10, noise_gen1 / 10)
    noise_gen -= noise_int  # В милливатах
    noise_gen *= 2
    noise_gen1 -= noise_int  # В милливатах
    noise_gen1 *= 10

    # f = interpolate.interp1d(freq, noise_gen, kind='cubic')
    f = interpolate.UnivariateSpline(freq, noise_gen1, w=None, bbox=[None, None], k=5, s=1, ext=0, check_finite=False)
    f2 = interpolate.UnivariateSpline(freq, noise_gen1, w=None, bbox=[None, None], k=5, s=1, ext=0, check_finite=False)
    x = np.arange(500000000, 3500000000, 4000000)
    y1 = f(x)
    y2 = f2(x1)

    fig, axes = plt.subplots(1, 1)
    plot_data = axes.plot(freq, noise_gen)
    plot_data1 = axes.plot(freq, noise_gen1)
    plot_data2 = axes.plot(x, y1)
    plot_data3 = axes.plot(x1, y2)
    plt.show()

    return y2


# Путь к папке с данными шумовых калибровок
path_input = r'E:\Measure_res\2020_12_24Calibrate'
# Путь к папке с коэффициентами корректировки АЧХ тракта по результатам
# шумовых измерений
path_output = r'E:\Measure_res\Calibr_coeff_2020_12'
# Путь к файлу данных шумовой калибровки входов "Ant1" или "Ant2"
# и зоны Найквиста "2" или "3"

dict_input = {'Ant1_2': r'\20201224_Ant1_NGen_2_left.txt',
              'Ant2_2': r'\20201224_Ant2_NGen_2_right.txt',
              'Ant1_3': r'\20201224_Ant1_NGen_3_left.txt',
              'Ant2_3': r'\20201224_Ant2_NGen_3_right.txt'}
dict_coeff = {}
dict_spectrum = {}
dict_spectrum_norm = {}
keys = dict_input.keys()
for key in keys:
    file_name0: str = path_input + dict_input[key]

    freq = np.zeros(256)
    dict_coeff[key], freq, dict_spectrum[key] = af(file_name0)
    dict_spectrum[key] *= dict_coeff[key]
    pass
dict_spectrum['Ant2_2'] /= 2   # В этом измерении ослабление было на 3 дБ меньше, чем в других
spectrum_level = [np.mean(dict_spectrum[s][5:20]) for s in keys]
spectrum_max = np.max(spectrum_level)
spectrum_level = [spectrum_max / s for s in spectrum_level]
i = 0
for s in keys:
    coeff = np.array(dict_coeff[s])
    coeff *= spectrum_level[i]
    # Путь записи файла с коэффициентами корректировки АЧХ тракта для входов "Ant1"
    # или "Ant2" и зоны Найквиста "2" или "3"
    file_name_out: str = path_output + r'\Calibr_' + s + '.txt'
    # Расчет и запись коэффициентов корректировки АЧХ в файл
    # При этом для второй зоны Найквиста - инверсный порядок по частоте
    # np.savetxt(file_name_out, coeff)

ngs = np.zeros(256)
ngs = noise_gen_spectrum(freq)
coeff *= ngs
coeff_max = np.max(coeff)
coeff /= coeff_max

# Расчет и запись коэффициентов корректировки АЧХ в файл
# При этом для второй зоны Найквиста - инверсный порядок по частоте


pass