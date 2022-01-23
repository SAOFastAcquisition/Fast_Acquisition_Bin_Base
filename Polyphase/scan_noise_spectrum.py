import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft
from pathlib import Path
from Supporting_func import path_to_data


def low_freq_noise_spectrum(_spectrum, _m):
    # _spectrum - сигнал радиометра в виде скана мощности в полосе частот по времени
    # _m - длина выборки для построения спектра низкочастотных шумов сигнала радиометра
    _shape = _spectrum.shape
    n = int(_shape[1] // _m * _m)
    _spectrum_signal_av = np.zeros((_shape[0], _m))
    for i in range(_shape[0]):
        a = _spectrum[i, 0:n]
        signal_rand_sample = np.reshape(a, (_m, -1))  # Разбиваем на реализации длиной _m отсчетов
        signal_rand_sample_t = np.transpose(signal_rand_sample)
        spectrum_signal_rand = fft(signal_rand_sample_t, _m)
        _spectrum_signal_av[i, :] = np.average(np.abs(spectrum_signal_rand ** 2), 0)
    return _spectrum_signal_av


def plot_low_freq_spec(_spectrum, _delta_t):
    m, n = _spectrum.shape
    f_max = 1 / _delta_t / 2
    f_min = f_max / n
    arg = np.linspace(f_min, f_max, n)
    axes = plt.subplots()
    for i in range(m):
        plt.plot(arg, _spectrum[i, :])
    plt.show()


if __name__ == '__main__':

    # ********************************** Путь к файлу данных ****************************
    current_catalog = r'2022\Results'  # Текущий каталог (за определенный период, здесь - год)
    current_data_dir = '2022_01_20test'  # Папка с текущими даннымиH:\Fast_Acquisition\2021\Results\2021_09_22test\2021-09-22_01_14bit_pm20
    current_data_file = '2022-01-20_01test'  # Имя файла с исходными текущими данными без расширенияH:\Fast_Acquisition\2022\Results\2022_01_20test
    file_path_data, head_path = path_to_data(current_catalog, current_data_dir)

    #           ************** Загрузка матрицы спектров *************
    spectrum = np.load(Path(file_path_data, current_data_file + '_spectrum.npy'), allow_pickle=True)
    #               **********************************************
    spectrum0 = spectrum[0]
    spectrum1 = spectrum[1]
    spectrum2 = spectrum[2]
    spectrum3 = spectrum[3]
    delta_t = 8.1925e-3

    spectrum_signal_av = low_freq_noise_spectrum(spectrum1, 1024)
    plot_low_freq_spec(spectrum_signal_av, delta_t)

    pass