import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.fftpack import fft
from pathlib import Path
# import matplotlib.font_manager as font_manager
# from matplotlib.ticker import MaxNLocator, ScalarFormatter, FixedLocator
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
    return _spectrum_signal_av / _m


def plot_low_freq_spec(_spectrum, _delta_t, _path_to_picture_folder, _line_legend):
    _path_to_picture = path_to_pic(_path_to_picture_folder)
    m, n = _spectrum.shape
    f_max = 1 / _delta_t / 2
    f_min = f_max / n
    arg = np.linspace(f_min, f_max, n)

    fig = plt.figure(figsize=(12, 8))
    axes = fig.add_subplot()
    axes.grid(b=True, which='major', color='#666666', linestyle='-')
    # Show the minor grid lines with very faint and almost transparent grey lines
    axes.minorticks_on()
    axes.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)

    axes.set_title('Temperature Low Noise Spectral Density', fontsize=18)
    axes.set_yscale('log')
    axes.set_xscale('log')
    axes.set_xlabel('Freq, Hz', fontsize=18)
    axes.set_ylabel('Spectral Density, $K^2/Hz$', fontsize=20)
    for i in range(m):
        axes.plot(arg[0:n // 2], _spectrum[i, 0:n // 2], label=_line_legend[i])
    axes.legend()
    plt.show()

    fig.savefig(_path_to_picture)


def plot_low_freq_spec_ab(_spectrum, _delta_t, _path_to_picture_folder, _line_legend):
    _path_ad = Path(_path_to_picture_folder, 'LN_spectrum.npy')
    _spectrum_ad = np.load(_path_ad, allow_pickle=True)
    _path_to_picture = path_to_pic(_path_to_picture_folder, 'eps')
    m, n = _spectrum.shape
    f_max = 1 / _delta_t / 2
    f_min = f_max / n
    arg = np.linspace(f_min, f_max, n)

    matplotlib.rcParams['font.family'] = 'serif'
    matplotlib.rcParams['font.size'] = '10'

    fig = plt.figure(figsize=(6, 4))
    axes = fig.add_subplot()
    axes.tick_params(axis='both',  # Применяем параметры к обеим осям
                     which='major',  # Применяем параметры к основным делениям
                     direction='in',  # Рисуем деления внутри и снаружи графика
                     length=3,  # Длинна делений
                     width=1,  # Ширина делений
                     color='black',  # Цвет делений
                     pad=2,  # Расстояние между черточкой и ее подписью
                     # labelsize=f_size,  # Размер подписи
                     labelcolor='black',  # Цвет подписи
                     bottom=True,  # Рисуем метки снизу
                     top=True,  # сверху
                     left=True,  # слева
                     right=True,  # и справа
                     labelbottom=True,  # Рисуем подписи снизу
                     labeltop=False,  # сверху
                     labelleft=True,  # слева
                     labelright=False,  # и справа
                     labelrotation=0)  # Поворот подписей
    axes.tick_params(axis='both',  # Применяем параметры к обеим осям
                     which='minor',  # Применяем параметры к вспомогательным делениям
                     direction='in',  # Рисуем деления внутри и снаружи графика
                     length=2,  # Длинна делений
                     width=1,  # Ширина делений
                     color='black',  # Цвет делений
                     pad=10,  # Расстояние между черточкой и ее подписью
                     labelsize=10,  # Размер подписи
                     labelcolor='black',  # Цвет подписи
                     bottom=True,  # Рисуем метки снизу
                     top=True,  # сверху
                     left=True,  # слева
                     right=True)  # и справа
    # axes.grid(b=False, which='major', color='#666666', linestyle='-')
    # Show the minor grid lines with very faint and almost transparent grey lines
    # axes.minorticks_off()
    # axes.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)

    # axes.set_title('Temperature Low Noise Spectral Density')
    axes.set_ylim([1, 1000000])
    axes.set_yscale('log')
    axes.set_xscale('log')
    axes.set_xlabel('Freq, Hz')
    axes.set_ylabel('Spectral Density, $K^2/Hz$')
    plt.text(0.7, 15, '1', style='italic')  # Семейство кривых
    plt.text(0.7, 10000, '2', style='italic')  #
    for i in range(m):
        axes.plot(arg[0:n // 2], _spectrum[i, 0:n // 2], color='black', linewidth=0.5)
        axes.plot(arg[0:n // 2], _spectrum_ad[i, 0:n // 2], color='black', linewidth=0.5)
    # axes.legend()
    plt.show()

    fig.savefig(_path_to_picture)


def path_to_pic(_file_path, format='png'):
    add_pass0 = 'LF_spectrum_00'
    _l = len(add_pass0)
    add_pass1 = add_pass0 + '.' + format
    if not os.path.isfile(Path(_file_path, add_pass1)):
        if not os.path.isdir(_file_path):
            os.makedirs(_file_path)
    else:
        while os.path.isfile(Path(_file_path, add_pass1)):
            num = int(add_pass0[_l - 2:_l]) + 1
            num_str = str(num)
            if num >= 10:
                add_pass0 = add_pass0[:_l - 2] + num_str
            else:
                add_pass0 = add_pass0[:_l - 2] + '0' + num_str
            add_pass1 = add_pass0 + '.' + format
    _path_to_picture = Path(_file_path, add_pass1)
    return _path_to_picture


if __name__ == '__main__':
    # ********************************** Путь к файлу данных ****************************
    current_catalog = r'2022\Results'  # Текущий каталог (за определенный период, здесь - год)
    current_data_dir = r'2022_01_20test'  # Папка с текущими даннымиH:\Fast_Acquisition\2021\Results\2021_09_22test\2021-09-22_01_14bit_pm20
    current_data_file = r'2022-01-20_02test'  # Имя файла с исходными текущими данными без расширенияH:\Fast_Acquisition\2022\Results\2022_01_20test
    file_path_data, head_path = path_to_data(current_catalog, current_data_dir)

    #           ************** Загрузка матрицы спектров *************
    spectrum = np.load(Path(file_path_data, current_data_file + '_spectrum.npy'), allow_pickle=True)
    #               **********************************************
    spectrum0 = spectrum[0]
    spectrum1 = spectrum[1]
    spectrum2 = spectrum[2]
    spectrum3 = spectrum[3]
    delta_t = 8.1925e-3
    spectrum_current = np.transpose(spectrum1[:, 20:22])
    spectrum_signal_av = low_freq_noise_spectrum(spectrum_current, 16384)
    legend = ['freq20', 'freq21']
    plot_low_freq_spec(spectrum_signal_av, delta_t, Path(file_path_data, current_data_file), legend)

    pass
