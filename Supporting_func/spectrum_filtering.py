import numpy as np
from matplotlib import pyplot as plt
from numpy import fft
from scipy.signal.windows import kaiser, flattop
from scipy import signal

from pathlib import Path
from scipy.fftpack import fft, fftshift
from Supporting_func.stocks_coefficients import path_to_data


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


def kernel_cell(n, m):
    """ Функция отдает двумерное сверточное ядро размером (m x n) на основании максимально плоского фильтра"""
    w_inter_freq = np.array(flattop(n))
    w_inter_time = np.array(flattop(m))
    kernel_in = np.ones((m, n), dtype=np.float)
    for i in range(m):
        kernel_in[i] = w_inter_freq
    for i in range(n):
        kernel_in[:, i] = kernel_in[:, i] * w_inter_time
    kernel_in = kernel_in / np.sum(kernel_in)
    return kernel_in


model = 'n'
n = 32  # Фильтрация по частоте (постоянная фильтра примерно 2*n/3 отсчетов)
m = 1   # Фильтрация по времени (постоянная фильтра примерно 2*m/3 отсчетов)
current_data_file = '2021-06-28_20-30'  # Имя файла с исходными текущими данными без расширения
current_data_dir = '2021_06_28sun'      # Папка с текущими данными
current_catalog = r'2021\Results'       # Текущий каталог (за определенный период, здесь - год)

file_path_data, head_path = path_to_data(current_catalog, current_data_dir)
spectrum = np.load(Path(file_path_data, current_data_file + '_spectrum.npy'), allow_pickle=True)
spectrum_loc = spectrum[2]
spectrum_trace_control = spectrum_loc[10]
spectrum_signal_av1 = spectrum_trace_control

kernel = kernel_cell(n, m)

if model == 'y':
    signal_rand = model_signal()
    h = filter_coeff(1024, 8, 900)
    fft_h = abs(fft(h, 1024) ** 2)
    signal_rand = low_freq_filter(signal_rand, h)
    signal_rand_sample = np.reshape(signal_rand, (-1, 1024))  # Разбиваем на реализации длиной 1024 отсчетов
    spectrum_signal_rand = fft(signal_rand_sample, 1024, axis=1)
    spectrum_signal = np.abs(spectrum_signal_rand ** 2)
    spectrum_signal_av = np.average(spectrum_signal, 0)
else:
    for i in range(np.size(spectrum)-1):
        l = np.shape(spectrum[i])
        if np.size(l) == 2 and l[1] > 4:
            blurred = signal.fftconvolve(spectrum[i], kernel, mode='same')
            spectrum[i] = blurred
    spectrum_loc = spectrum[2]
    spectrum_trace = spectrum_loc[10]
    spectrum_signal_av = spectrum_trace

# kernel = np.outer(signal.gaussian(70, 8), signal.gaussian(70, 8))
l = np.shape(spectrum[3])
fft_h = abs(fft(kernel, l[1]) ** 2)

axes = plt.subplots()
plt.plot(spectrum_signal_av)
plt.plot(spectrum_signal_av1)
plt.show()

signal1 = fft(spectrum_signal_av, 512)
signal2 = fft(spectrum_signal_av1, 512)

axes = plt.subplots()
plt.plot(signal1)
plt.plot(signal2)
plt.plot(fft_h[0] * 1e15)
plt.show()

np.save(Path(file_path_data, current_data_file + '_1spectrum'), spectrum)
