import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from numpy import fft
from scipy.signal.windows import kaiser, flattop
import filters as ftr
from pathlib import Path
from scipy.fftpack import fft, fftshift
from Supporting_func.stocks_coefficients import path_to_data


def low_freq_filter(x, h):
    n_input = len(x)
    n_filter = len(h)
    # n - Длина входной последовательности
    # h - отклик фильтра НЧ

    y = [0] * n_input  # Выходная последовательность после интерполяции и НЧ фильтрации
    for n in range(0, n_input):
        for k in range(0, n_filter):
            try:
                y[n] += x[n - k - 1] * h[k]
            except IndexError:
                y[n] += 0
                print('ind n = ', n, 'ind k = ', k)
                pass
    return y


def filter_coeff(length_fft, filters_order, band_pass):
    #   length_fft - длина БПФ
    #   filters_order - Порядок фильтра
    #   band_pass - полоса фильтра в отсчетах
    h, fft_h = ftr.synt_filter(length_fft, filters_order, band_pass)  # Отклик прямоугольного модельного фильтра
    n_filter = len(h)
    h_short = h[n_filter // 2 - filters_order // 2:n_filter // 2 + filters_order // 2]  # Выбор количества значений
    # отклика, соответствующего порядку фильтра
    # y_int_wind = interpolate_filter(4, x_wind, h_short)
    # w_inter = [0.54 - 0.46 * np.cos(2 * np.pi * i / (n - 1)) for i in range(0, n * l_interpol)] # Окно Хэминга
    # w_inter = flattop(n)  # from scipy максимально плоское окно (применяется при полифазной обработке)
    return h_short


def random_signal(n):
    primary_sig = np.random.sample(n) - 0.5
    return primary_sig


def model_signal():
    m = 2 ** 16
    signal_sin = [0.1 * (np.sin(2 * np.pi * 0.2 * i) + 1.2 * np.sin(2 * np.pi * 0.29 * i)) for i in range(m)]
    signal_rand = random_signal(m) + signal_sin
    sig_mean = np.mean(signal_rand)
    sig_var = np.var(signal_rand)
    return signal_rand


model = 'n'
current_data_file = '2021-06-28_19-28'      # Имя файла с исходными текущими данными без расширения
current_data_dir = '2021_06_28sun'          # Папка с текущими данными
current_catalog = r'2021\Results'           # Текущий каталог (за определенный период, здесь - год)

file_path_data, head_path = path_to_data(current_catalog, current_data_dir)
spectrum = np.load(Path(file_path_data, current_data_file + '_spectrum.npy'), allow_pickle=True)
spectrum_loc = spectrum[2]
spectrum_trace_control = spectrum_loc[10]
# spectrum_signal_av1 = np.average(spectrum[2], 0)
spectrum_signal_av1 = spectrum_trace_control

if model == 'y':
    signal_rand = model_signal()
    h = filter_coeff(1024, 64, 900)
    fft_h = abs(fft(h, 1024) ** 2)
    signal_rand = low_freq_filter(signal_rand, h)
    signal_rand_sample = np.reshape(signal_rand, (-1, 1024))  # Разбиваем на реализации длиной 1024 отсчетов
    spectrum_signal_rand = fft(signal_rand_sample, 1024, axis=1)
    spectrum_signal = np.abs(spectrum_signal_rand ** 2)
    spectrum_signal_av = np.average(spectrum_signal, 0)
else:
    for i in range(np.size(spectrum)):
        l = np.shape(spectrum[i])
        if np.size(l) == 2 and l[1] > 4:
            h = filter_coeff(l[1], 64, int(l[1]//64))
            fft_h = abs(fft(h, l[1]) ** 2)
            spectrum_loc = spectrum[i]
            for j in range(l[0]):
                spectrum_line = spectrum_loc[j]
                spectrum_line = np.abs(low_freq_filter(spectrum_line, h))
                spectrum_loc[j] = spectrum_line
                pass
            spectrum[i] = spectrum_loc
    spectrum_loc = spectrum[2]
    spectrum_trace = spectrum_loc[10]
    spectrum_signal_av = spectrum_trace
    # spectrum_signal_av = np.average(spectrum[2], 0)


axes = plt.subplots()
plt.plot(spectrum_signal_av)
plt.plot(spectrum_signal_av1)
plt.show()

signal1 = fft(spectrum_signal_av, 512)
signal2 = fft(spectrum_signal_av1, 512)

axes = plt.subplots()
plt.plot(signal1)
plt.plot(signal2)
plt.plot(fft_h)
plt.show()

# np.save(Path(file_path_data, current_data_file + '_1spectrum'), spectrum_signal)
