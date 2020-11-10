import numpy as np
from matplotlib import pyplot as plt
from numpy import fft
from scipy.signal.windows import kaiser, flattop
import filters as ftr


def interpolate_filter(l_interpol, x, h):
    n_input = len(x)
    n_filter = len(h)
    # n - Длина входной последовательности
    # l_interpol - коэффициент интерполяции
    # h - отклик фильтра НЧ для удаления лишних изображений

    y = [0] * (n_input*l_interpol)         # Выходная последовательность после интерполяции и НЧ фильтрации
    len_filter_part = int(n_filter//l_interpol)
    for n in range(0, n_input):
        for i in range(0, l_interpol):
            for k in range(0, len_filter_part):
                try:
                    y[n*l_interpol+i] += x[n+len_filter_part-k-1]*h[l_interpol * k + i]
                except IndexError:
                    y[n * l_interpol + i] += 0
                    print('ind i = ', i, 'ind j = ', n*l_interpol+i)
                    pass
    return y


if __name__ == '__main__':

    n = 128                 # Длина исходной последовательности
    m = 32                  # Порядок фильтра должен быть кратным коэффициенту интерполяции
    l_interpol = 4
                            # Исходная последовательность
    x = [np.sin((2*np.pi*((i+1)*0.20)))+ np.sin((2*np.pi*((i+1)*0.24))) for i in range(0, n)]   #  + np.sin((2*np.pi*((i+1)*0.160)))
    w = [0.54 - 0.46 * np.cos(2 * np.pi * i / (n - 1)) for i in range(0, n)]    # Окно
    # x_wind = [x[i] * w[i] for i in range(0,n)]          # Взвешенная входная последовательность

    h_rect = [1 for i in range(0, 32)]          # Отклик фильтра с АЧХ типа sin(x)/x
    h, fft_h = ftr.synt_filter(512, 32, 128)    # Отклик прямоугольного модельного фильтра
    n_filter = len(h)
    h_short =h[n_filter // 2 - m // 2:n_filter // 2 + m // 2]   # Выбор количества значений отклика, соответствующего
                                                                # порядку фильтра
    y_int = interpolate_filter(4, x, h_short)
    y_int_rect = interpolate_filter(4, x, h_rect)
    # y_int_wind = interpolate_filter(4, x_wind, h_short)
    # w_inter = [0.54 - 0.46 * np.cos(2 * np.pi * i / (n - 1)) for i in range(0, n * l_interpol)] # Окно Хэминга
    w_inter = flattop(n * l_interpol)  # from scipy максимально плоское окно (применяется при полифазной обработке)
    y_int1 = [y_int[i] * w_inter[i] for i in range(0, n * l_interpol)]
    y_int_np = np.array(y_int1)
    y1 = y_int_np.reshape(n, l_interpol) # 
    y1 = y1.sum(1)
    y_int_wind = np.fft.fft(y1)         # Спектр исходного сигнала после интерполяции на 4 и полифазной обработки

    synt_filt = np.fft.fft(y_int)
    synt_filter = np.fft.fft(h)
    synt_filt_rect = np.fft.fft(y_int_rect)
    # synt_filt_wind = np.fft.fft(y_int_wind)

    xf = np.fft.fft(x)

    figure, axes = plt.subplots(2, 2)

    axes[0, 0].plot(np.abs(np.fft.fft(x)))
    axes[0, 0].plot(np.abs(y_int_wind))

    axes[1, 0].plot(np.abs(synt_filt))
    axes[1, 0].plot(np.abs(synt_filter) * 100)

    axes[0, 1].plot(y_int)
    axes[0, 1].plot(y_int1)
    axes[0, 1].plot(y_int_rect)

    axes[1, 1].plot(x)
    axes[1, 1].plot(y1)
    # axes[1].plot(y1)
    plt.show()



