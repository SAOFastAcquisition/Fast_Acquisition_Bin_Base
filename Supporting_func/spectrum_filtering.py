import numpy as np
from matplotlib import pyplot as plt
from numpy import fft
from scipy.signal.windows import kaiser, flattop
import filters as ftr
from scipy.fftpack import fft, fftshift


def low_freq_filter(x, h):
    n_input = len(x)
    n_filter = len(h)
    # n - Длина входной последовательности
    # l_interpol - коэффициент интерполяции
    # h - отклик фильтра НЧ для удаления лишних изображений

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



m = 2 ** 16
signal_sin = [0.1 * (np.sin(2 * np.pi * 0.3 * i) + 0.75 * np.sin(2 * np.pi * 0.17 * i)) for i in range(m)]
signal_rand = random_signal(m) + signal_sin
sig_mean = np.mean(signal_rand)
sig_var = np.var(signal_rand)

h = filter_coeff(1024, 256, 512)
signal_rand = low_freq_filter(signal_rand, h)


signal_rand_sample = np.reshape(signal_rand, (1024, -1))    # Разбиваем на реализации длиной 1024 отсчетов
spectrum_signal_rand = fft(signal_rand_sample, 1024)
spectrum_signal_av = np.average(np.abs(spectrum_signal_rand ** 2), 0)

axes = plt.subplots()
plt.plot(spectrum_signal_av)
plt.show()
