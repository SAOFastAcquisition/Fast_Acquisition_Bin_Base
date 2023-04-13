import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema as agl
from Help_folder.paths_via_class import DataPaths


def plot_imf(_x, _y, _y1, _y2):
    _fig = plt.figure(figsize=(12, 8))
    _axes = _fig.add_subplot()
    _axes.grid(b=True, which='major', color='#666666', linestyle='-')
    # Show the minor grid lines with very faint and almost transparent grey lines
    _axes.minorticks_on()
    _axes.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)

    _axes.set_title('Antenna Temperature', fontsize=18)
    _axes.set_xlabel('Frequency, MHz', fontsize=18)
    # plt.ylim([0, 1000000])
    _axes.plot(_x, _y)
    _axes.plot(_x, _y1)
    _axes.plot(_x, _y2)
    plt.show()


def some_fig_plot(_data, _path_to_fig_folder=None):
    """
    Функция принимает путь сохранения рисунка и три возможных последовательности для построения двух или трех
    графиков с общим аргументом mean_frame_ind_pol, который приходит от вызывающей функции. При этом наличие
    двух отображаемых на рис. последовательностей обязательно, последовательность _s_i присутствует всегда.
    :param _path_to_fig_folder: Путь к папке для сохранения рисунка
    :param _data:
    :param _path_to_fig_folder:
    :param :
    :return:
    """
    # _pic_name = pic_name(_path_to_fig_folder, 0, 'png')
    # _path_to_pic = Path(_path_to_fig_folder, _pic_name)
    # _fig_folder = str(_path_to_fig_folder)
    # title1, title2, title3 = title_func(_fig_folder, head)
    _l, _m = np.shape(_data)
    _df = 2000 / _m
    _freq = [1000 + _df / 2 + _df * _i for _i in range(_m)]
    fig, _axes = plt.subplots(_l, 1, figsize=(12, _l * 3))
    for _j in range(_l):
        _axes[_j].plot(_freq, _data[_j, :])
        _axes[_j].grid()
        _axes[_j].grid(which='minor',
                       axis='x',
                       color='k',
                       linestyle=':')
    _sum = np.sum(_data, axis=0)
    _axes[_l - 1].plot(_freq, _sum)

    # axes[0].set_title('Stokes Parameters ' + title1, fontsize=20)
    # axes[0].set_ylabel('Stokes_I')
    # if _s_v is not None:
    #     axes[1].set_ylabel('Stokes_V', color='darkred')
    # else:
    #     axes[1].set_ylabel('Stokes_V Deviation', color='darkred')
    # axes[0].minorticks_on()
    # axes[1].minorticks_on()
    # if _s_v is not None and _s_dv is not None:
    #     axes[2].set_ylabel('Stokes_V Deviation', color='darkred')
    #     axes[2].minorticks_on()
    #     axes[2].grid()
    #     axes[2].grid(which='minor',
    #                  axis='x',
    #                  color='k',
    #                  linestyle=':')
    # y1 = y_max - 2 * (y_max - y_min) / 10
    # y2 = y_max - 3 * (y_max - y_min) / 10
    # axes[0].text(0, y1, inform[0], fontsize=12)  # Разрешение по частоте
    # axes[0].text(0, y2, inform[1], fontsize=12)  # Разрешение по времени

    plt.show()
    # #                               ********************************
    # #                        ************ Сохранение рисунка ****************
    # fig.savefig(_path_to_pic)
    # flag_save = save_question()
    # if flag_save == 'no':
    #     if os.path.isfile(_path_to_pic):
    #         os.remove(_path_to_pic)
    #         print('Picture is not saved')
    #     else:
    #         print('File not found')
    # else:
    #     print('Picture is saved')
    # pass


def zone_deletion(_len):
    if tag:
        # Исключение зон действия режекторных фильтров при правильном порядке отсчетов частоты во второй зоне Найквиста
        _delta_f = 2000 / _len
        k1 = int((25 - _delta_f / 2) // _delta_f)  #
        k2 = int((770 - _delta_f / 2) // _delta_f)  #
        k3 = int((1034 - _delta_f / 2) // _delta_f)
        k4 = int((1090 - _delta_f / 2) // _delta_f)  #
        k5 = int((1230 - _delta_f / 2) // _delta_f)
        k6 = int((1525 - _delta_f / 2) // _delta_f)
        k7 = int((1710 - _delta_f / 2) // _delta_f)
        k8 = int((1954 - _delta_f / 2) // _delta_f)
        k9 = int(2000 / _delta_f) - 1
        _k = [0, k1, k2, k3, k4, k5, k6, k7, k8, k9]
    else:
        k1 = int(6 / delta_t)
        k2 = _len - k1
        _k = [0, k1, k2, _len - 1]
    return _k


def fill_zone_del(_data):
    """
    Принимает исходные спектры с пропусками в месте действия режекторных фильтров и заполняет их с помощью
    функции fill_func() для избежания разрывов первого рода
    :param _data:
    :return _data: исходные данные с удаленными разрывами первого рода
    """
    _len_data = len(_data)
    _k = zone_deletion(_len_data)
    if tag:
        _df = 2000 / _len_data
        _x_init = np.array([1000 + _df / 2 + _df * _i for _i in _k])
        _y_init = _data[_k]
        _len_init = len(_k)
        _y_init[0] = _y_init[1] - 0.1
        _y_init[_len_init - 1] = _y_init[_len_init - 2] - 0.15
        _kr = _k[0:2]
        reg0 = fill_func(_x_init[[0, 1]], _y_init[[0, 1]], _kr, 2)
        _data[_k[0]:_k[1] + 1] = reg0
        reg1 = fill_func(_x_init[[2, 3]], _y_init[[2, 3]], _k[2:4], 8)
        _data[_k[2]:_k[3] + 1] = reg1
        reg2 = fill_func(_x_init[[4, 5]], _y_init[[4, 5]], _k[4:6], 5)
        _data[_k[4]:_k[5] + 1] = reg2
        reg3 = fill_func(_x_init[[6, 7]], _y_init[[6, 7]], _k[6:8], 8)
        _data[_k[6]:_k[7] + 1] = reg3
        reg4 = fill_func(_x_init[[8, 9]], _y_init[[8, 9]], _k[8:], 2)
        _data[_k[8]:_k[9] + 1] = reg4
    else:
        _x_init = np.array([delta_t * _i for _i in _k])
        _y_init = _data[_k]
        _y_init[0] = _y_init[1]
        _y_init[3] = _y_init[2]
        reg0 = fill_func(_x_init[[0, 1]], _y_init[[0, 1]], _k[0:2], 8)
        _data[_k[0]:_k[1] + 1] = reg0
        reg1 = fill_func(_x_init[[2, 3]], _y_init[[2, 3]], _k[2:], 8)
        _data[_k[2]:_k[3] + 1] = reg1
    plt.plot(_data)
    plt.show()
    pass
    return _data


def fill_func(_x_init, _y_init, _k_init, _order):
    """
    Функция заполняет зону действия режекторных фильтров. Принимает начальные условия слева и справа и соответствующие
    им индексы исходнго вектора данных _k_init. Возвращает вектор для заполнения зоны.
    :param _x_init:
    :param _y_init:
    :param _k_init:
    :param _order:
    :return:
    """
    _k = np.array([i for i in range(_k_init[0], _k_init[1] + 1, 1)])
    _kl = _k_init[1] - _k_init[0]
    if tag:
        _y = [_y_init[0] + (_y_init[1] - _y_init[0]) / _kl * (i - _k_init[0]) \
              + 0.03 * np.cos(2 * 3.14 / _kl * (i - _k_init[0])) \
              + 0.01 * np.sin(2 * 3.14 * _order / _kl * (i - _k_init[0])) for i in _k]
    else:
        _y = [_y_init[0] + (_y_init[1] - _y_init[0]) / _kl * (i - _k_init[0]) \
              + 0.001 * np.sin(2 * 3.14 * _order / _kl * (i - _k_init[0])) for i in _k]
    # plt.plot(_y)
    # plt.show()
    pass
    return _y


def imf_gen(_data):
    """
    Принимает исходные данные и выделяет в них самую быструю компоненту
    :param _data:
    :return _cs: Самая быстрая компонента
    """
    from scipy.interpolate import make_interp_spline as mis
    _l = len(_data)
    cx = np.array([i for i in range(_l)])
    _idx_minimas = agl(_data, np.less)[0]
    _idx_maximas = agl(_data, np.greater)[0]
    if len(_idx_maximas) <= 1 and len(_idx_minimas) <= 1:
        return print(f"IMF is over: min bin = {_idx_minimas}, max bin = {_idx_maximas}")
    # _cs_min = CubicSpline(_idx_minimas, _data[_idx_minimas])
    l, r = [(2, 0)], [(2, 0)]
    _cs_min = mis(_idx_minimas, _data[_idx_minimas], k=3, bc_type=(l, r))
    # _cs_max = CubicSpline(_idx_maximas, _data[_idx_maximas])
    _cs_max = mis(_idx_maximas, _data[_idx_maximas], k=3, bc_type=(l, r))
    _cs = (_cs_max(cx) + _cs_min(cx)) / 2

    # plot_imf(cx, _data, _cs, _data - _cs)
    return _cs


def imf_proc(_data):
    """
    Выделяет собственно эмпирическую функцию, добиваясь того, чтобы у самой быстрой компоненты количество пересечений
    нуля равнялось сумме минимумов и максимумов и эти три числа повторились для двух иттераций
    :param _data:
    :return _cs01: Собственная эмпирическая функция
    """
    flag = 0
    sign_change_count = 2
    n_min = 1
    n_max = 1
    data0 = _data.copy()
    cs0 = []
    while flag == 0:
        try:
            cs = imf_gen(data0)
            if len(cs0) == 0:
                cs0 = cs
            else:
                cs0 = np.vstack((cs0, cs))
            data1 = data0 - cs
        except TypeError:
            print('IMFs is over')
            cs0 = np.zeros((1, l))
            break

        sign_change_count1 = null_cross_sec(data1)
        n_min1 = len(agl(data1, np.less)[0])
        n_max1 = len(agl(data1, np.greater)[0])
        if (sign_change_count1 - sign_change_count) == 0 and (n_min1 - n_min) == 0 and (n_max1 - n_max) == 0:
            flag = 1
        else:
            sign_change_count = sign_change_count1
            n_min = n_min1
            n_max = n_max1
            data0 = data1
    _cs01 = np.nansum(cs0, axis=0)
    return _cs01


def imf_decomp(_data):
    """
    Разложение исходного вектора данных на собственные эмпирические компоненты
    :param _data:
    :return:
    """
    _data_cur = _data.copy()
    _imf = []
    while True:
        cs_cur = imf_proc(_data_cur)
        imf_cur = _data_cur - cs_cur
        if not len(_imf):
            _imf = imf_cur
        else:
            _imf = np.vstack((_imf, imf_cur))
        if not np.sum(cs_cur):
            break
        _data_cur = cs_cur
    return _imf


def null_cross_sec(_data):
    """
    Функция подсчета количества пересечений нуля последовательностью _data
    :param _data:
    :return:
    """
    from math import copysign
    only_ones_list = [copysign(1, element) for element in _data]

    summa = 0
    for i in range(len(only_ones_list) - 1):
        summa += abs(only_ones_list[i] + only_ones_list[i + 1])

    _sign_change_count = ((len(only_ones_list) - 1) * 2 - summa) / 2
    return int(_sign_change_count)


if __name__ == '__main__':
    current_primary_file = '2023-02-10_01+20'
    current_primary_dir = '2023_02_10sun'
    main_dir = '2023'
    adr1 = DataPaths(current_primary_file, current_primary_dir, main_dir)
    converted_data_file_path = adr1.converted_data_file_path
    data_treatment_file_path = adr1.treatment_data_file_path
    e = 2.718281828459045
    delta_t = 8.3886e-3 * 32
    dec = np.log(10)
    #                               **************************
    # Загрузка исходных данных в виде спектров в фиксированные моменты времени '_spectrum_time.npy'
    # или сканов на фиксированных частотах '_scan_freq.npy'
    #                               **************************
    # path_npy = Path(str(converted_data_file_path) + '_spectrum_time.npy')
    path_npy = Path(str(converted_data_file_path) + '_scan_freq.npy')
    tag = 0
    if 'spectrum_time' in str(path_npy):
        tag = 1
    spectrum = np.load(path_npy, allow_pickle=True)
    mask = spectrum[3, :] == 40.0
    spectrum_log0 = np.log10(spectrum)
    data = spectrum_log0
    l = np.shape(data)[1]
    arg = np.arange(0, l, 1)
    #                               **************************
    imf00 = imf_decomp((data[6, :]))  # Разложение на собственные функции опорного спектра
    imf01 = imf_decomp(fill_zone_del(data[4, :]))
    imf02 = imf_decomp(fill_zone_del(data[5, :]))
    imf00[:, mask] = 0
    imf01[:, mask] = 0
    imf02[:, mask] = 0
    some_fig_plot(imf00)
    #                               **************************
    # data_mod1 = data[3, :] - imf0[0, :]
    # data_mod2 = data[3, :] - imf0[0, :] - imf0[1, :]
    # plot_imf(arg[:], np.exp(dec * data[5, :]), np.exp(dec * data_mod1), np.exp(dec * data_mod2))
    plot_imf(arg[:], imf00[0, :], imf01[0, :], imf02[0, :])
    plot_imf(arg[:], imf00[1, :], imf01[1, :], imf02[1, :])
    plot_imf(arg[:], imf00[2, :], imf01[2, :], imf02[2, :])
    plot_imf(arg[:], imf00[-1, :], imf01[-1, :], imf02[-1, :])
    plot_imf(arg[:], np.exp(dec * imf00[-1, :]), np.exp(dec * (data[5, :] - imf00[1, :])),
             np.exp(dec * data[5, :]))
    # imf1 = imf_decomp(data_mod1) - imf00[0, :]
    # imf2 = imf_decomp(data_mod2)
    # imf3 = imf_decomp(data[3, :])
    # plot_imf(arg[:], imf1[-1, :], imf2[-1, :], imf3[-1, :])

    pass
