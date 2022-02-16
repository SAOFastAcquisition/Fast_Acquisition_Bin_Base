import numpy as np
from matplotlib import pyplot as plt
from numpy import fft

def filter_AFC(n, k):
    """ Прямоугольный фильтр НЧ длиной k при общей длине выборки n, приведенный к половиной частоте дискретизации """

    filter_sample_freq = [0] * n

    for i in range(k):
        if k % 2 == 0:
            filter_sample_freq[n // 2 + i - k // 2] = 1
        else:
            filter_sample_freq[n // 2 - 1 + i - k // 2] = 1
    return filter_sample_freq


def filter_AFC1(n, k):

    """ Прямоугольный фильтр длиной k при общей длине выборки n, приведенный к нулевой частоте и частоте
    дискретизации"""
    filter_sample_freq = [0] * n
    k1 = int(k//2)

    if k % 2 == 0:
        for i in range(k1):
            filter_sample_freq[i] = 1
            filter_sample_freq[n-k1+i] = 1
    else:
        for i in range(k1):
            filter_sample_freq[i] = 1
        for i in range(k1+1):
            filter_sample_freq[n - k1 - 1 + i] = 1
    return filter_sample_freq


def synt_filter(n, m, k=8192):
    """
    :param n: Общая длина заны Найквиста для фильтра - разрешение по частоте
    :param m: Количество ненулевых значений усеченного отклика идеального синтезируемого фильтра
    :param k: Ширина в отсчетах в частотной области прямоугольного идеального фильтра
    :return:
    """
    a = filter_AFC1(n, k)
    _response = [0] * (n - 0)
    _response_filt = [0] * n
    FFT = np.fft.ifft(a)  # Прообраз АЧХ во временной области (отклик фильтра) с правым крылом правее отсчета 0 и
    # левым крылом левее последнего отсчета n-1
    FFT1 = np.abs(FFT)
    # Центрирование отклика относительно отсчета n//2
    _response[0:n//2] = FFT[n//2+0:n]

    _response[n//2:n-0] = FFT[0:n//2-0]
    _response_filt[n//2-m//2:n//2+m//2] = _response[n//2-m//2:n//2+m//2]
    _response_filt_out = _response[n//2-m//2:n//2+m//2]

    return _response_filt, FFT1


if __name__ == '__main__':
    n = 4096                       # Длина выборки
    m =1024                          # Длина фильтрующей выборки отклика
    k = 128
    response_filt, FFT1 = synt_filter(n, m, k)
    # response[7:n//2] = FFT[n//2+7:n]
    # response[n//2:n-5] = FFT[0:n//2-5]
    synt_filter = np.fft.fft(response_filt)

    figure, axes = plt.subplots(1, 1)
    # axes.plot(a)
    axes.plot(np.abs(synt_filter))
    plt.show()

    figure, ax = plt.subplots(1, 1)
    ax.plot(FFT1)
    # ax.plot(np.abs(response))
    plt.show()

