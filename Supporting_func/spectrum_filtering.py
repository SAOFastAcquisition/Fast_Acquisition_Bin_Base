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


def some_visualisation():
    spectrum_loc1 = spectrum[2]
    spectrum_trace = spectrum_loc1[10]
    l = np.shape(spectrum[2])
    # Частотная характеристика ядра свертки по оси частот в отсчетах
    fft_h = abs(fft(kernel, l[1]) ** 2)
    shape_kernel = np.shape(kernel)
    n_control = int(shape_kernel[0] // 2)
    axes = plt.subplots()
    plt.plot(spectrum_trace_control)
    plt.plot(spectrum_trace)
    plt.show()

    signal1 = np.abs(fft(spectrum_trace, l[1]) ** 2)
    signal2 = np.abs(fft(spectrum_trace_control, l[1]) ** 2)

    axes = plt.subplots()
    plt.plot(signal1)
    plt.plot(signal2)
    plt.plot(fft_h[n_control])
    plt.show()


if __name__ == '__main__':
    """ Программа принимает двумерную спектрограмму (спектры, развернутые во времени) (*.num) и фильтрует ее по частоте 
    и времени, сворачивая с задаваемым двумерным ядром. Результат записывается в виде файла *.num с '_kernel_spectrum' 
    в имени. 
    """
    model = 'y'
    visualization = 'y'
    n = 32  # Фильтрация по частоте (постоянная фильтра примерно 2*n/3 отсчетов)
    m = 8   # Фильтрация по времени (постоянная фильтра примерно 2*m/3 отсчетов)
    current_data_file = '2021-06-28_06+00'  # Имя файла с исходными текущими данными без расширения
    current_data_dir = '2021_06_28sun'      # Папка с текущими данными
    current_catalog = r'2021\Results'       # Текущий каталог (за определенный период, здесь - год)

    if model == 'y':
        signal_rand = model_signal()
        scan_loc = np.reshape(signal_rand, (-1, 1024))  # Разбиваем на реализации длиной 1024 отсчетов
        spectrum_loc = np.abs(fft(scan_loc, 1024, axis=1) ** 2)
        spectrum = np.array([np.nan, np.nan, spectrum_loc, np.nan])
    else:
        file_path_data, head_path = path_to_data(current_catalog, current_data_dir)
        spectrum = np.load(Path(file_path_data, current_data_file + '_spectrum.npy'), allow_pickle=True)
    spectrum_loc = spectrum[2]
    spectrum_trace_control = spectrum_loc[10]


    kernel = kernel_cell(n, m)
    # spectrum_signal_rand = fft(signal_rand_sample, 1024, axis=1)

    for i in range(np.size(spectrum)-1):
        l = np.shape(spectrum[i])
        if np.size(l) == 2 and l[1] > 4:
            """Разбиение файла данных на части по времени для избежания проблем с памятью 
            вызываемой функции свертки signal.fftconvolve()"""
            m = l[0] // 3000
            k = l[0] % 3000
            auxiliary = spectrum[i]
            for j in range(m+1):
                if j == m:
                    auxiliary_part = auxiliary[j * 3000:j * 3000 + k, :]
                    # Сворачиваем с ядром и превращеем в 'int32' для экономии памяти
                    blurred = signal.fftconvolve(auxiliary_part, kernel, mode='same').astype('int32')
                    auxiliary[j * 3000:j * 3000 + k, :] = blurred
                else:
                    auxiliary_part = auxiliary[j * 3000:(j+1)*3000, :]
                    blurred = signal.fftconvolve(auxiliary_part, kernel, mode='same').astype('int32')
                    auxiliary[j * 3000:(j+1)*3000, :] = blurred
                pass
            spectrum[i] = auxiliary

    # kernel = np.outer(signal.gaussian(70, 8), signal.gaussian(70, 8))

    if visualization == 'y':
        some_visualisation()

    # np.save(Path(file_path_data, current_data_file + '_kernel_spectrum'), spectrum)
