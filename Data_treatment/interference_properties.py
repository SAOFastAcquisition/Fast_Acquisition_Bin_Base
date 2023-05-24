import numpy as np
import pandas as pd
import pickle
import os
from pathlib import Path
import matplotlib.pyplot as plt
from Help_folder.paths_via_class import DataPaths
from tkinter import *
from tkinter import messagebox as mb


def data_preparing(_path):
    spectrum0 = np.load(_path, allow_pickle=True)
    with open(Path(str(converted_data_file_path) + '_head.bin'), 'rb') as inp:
        head = pickle.load(inp)
    n1, n2 = int(6 / delta_t), int(16 / delta_t)
    _n = n2 - n1
    spectrum = np.array([s[n1:n2, :] for s in spectrum0])

    # Приведение порядка следования отсчетов по частоте к нормальному
    if np.size(spectrum[0]):
        _n_row = np.shape(spectrum[0])[0]
        for i in range(_n_row):
            spectrum[0][i][0:] = spectrum[0][i][-1::-1]
    if np.size(spectrum[2]):
        _n_row = np.shape(spectrum[2])[0]
        for i in range(_n_row):
            spectrum[2][i][0:] = spectrum[2][i][-1::-1]

    if len(spectrum[0]) > 1 and len(spectrum[2]) > 1:  # Обе поляризации
        dat = np.hstack((spectrum[0] + spectrum[2], spectrum[1] + spectrum[3]))
    if len(spectrum[0]) > 1 and not (len(spectrum[1])) > 1:  # Левая поляризация
        dat = np.hstack((spectrum[0], spectrum[2]))
    if not (len(spectrum[0])) > 1 and len(spectrum[1]) > 1:  # Правая поляризация
        dat = np.hstack((spectrum[1], spectrum[3]))

    _n_aver = head['n_aver']
    return dat, _n, _n_aver


def data_preparing1(_path, _t1, _t2):
    """
    Функция обращается к записи спектров наблюдения и выбирает промежуток времени _t1, _t2 для вычсления rms
    (полуширины шумовой дорожки) скана и возвращает спектры dat_left, dat_right левой и правой поляризаций в этом
    промежутке времени вместе с количеством отсчетов по времени (строк)
    :param _path:
    :param _t1:
    :param _t2:
    :return:
    """
    spectrum0 = np.load(_path, allow_pickle=True)
    with open(Path(str(converted_data_file_path) + '_head.bin'), 'rb') as inp:
        head = pickle.load(inp)
    n1, n2 = int(_t1 / delta_t), int(_t2 / delta_t)
    _n = n2 - n1
    # spectrum = np.array([s[n1:n2, :] for s in spectrum0])
    spectrum = pd.Series([[], [], [], []])
    i = 0
    for s in spectrum0:
        if len(s) > 1:
            spectrum[i] = s[n1:n2, :]
        else:
            spectrum[i] = []
        i += 1
    # Приведение порядка следования отсчетов по частоте к нормальному
    if np.size(spectrum[0]):
        _n_row = np.shape(spectrum[0])[0]
        for i in range(_n_row):
            spectrum[0][i][0:] = spectrum[0][i][-1::-1]
    if np.size(spectrum[2]):
        _n_row = np.shape(spectrum[2])[0]
        for i in range(_n_row):
            spectrum[2][i][0:] = spectrum[2][i][-1::-1]

    if len(spectrum[0]) > 1 and len(spectrum[2]) > 1:  # Обе поляризации
        dat_left = np.hstack((spectrum[0], spectrum[1]))
        dat_right = np.hstack((spectrum[2], spectrum[3]))
    if len(spectrum[0]) > 1 and not (len(spectrum[2])) > 1:  # Левая поляризация
        dat_left = np.hstack((spectrum[0], spectrum[1]))
        dat_right = []
    if not (len(spectrum[0])) > 1 and len(spectrum[1]) > 1:  # Правая поляризация
        dat_right = np.hstack((spectrum[2], spectrum[3]))
        dat_left = []
    _n_aver = head['n_aver']

    return dat_left, dat_right, _n, _n_aver


def relative_rms(_y, _n):
    _x = np.array([i for i in range(_n)]) * delta_t
    mx = (_n - 1) * (_n - 1) / 2 / _n * delta_t  # 1
    my = np.sum(_y, axis=0) / _n
    a2 = np.dot(_x.T, _x) / _n
    a11 = np.dot(_x.T, _y) / _n

    kk = (a11 - mx * my) / (a2 - mx ** 2)
    bb = my - kk * mx
    ff = np.array([kk * z * delta_t + bb for z in range(_n)])

    delta_ff = ff - _y
    res = np.multiply(delta_ff, delta_ff)
    _rms = res ** 0.5
    _relative_rms = np.sum(_rms, axis=0) / _n / my

    return _relative_rms


def relative_rms1(_y, _n):
    """
    Функция вычсляет относительное среднеквадратичное отклонение. Возможное изменение тренда (среднего) за время набора
    отсчетов аппроксимируется линейно по времени и rms вычисляется как разность исходных данных и аппроксимации
    среднего. Во внимание принимаются только значащие отсчеты выбранной поляризации. Все отсчеты со значением ниже
    порогового отбрасываются и при усреднении не учитываются.
    Возвращает относительное к среднему среднеквадратичное отклонение в зависимости от частоты.
    :param _y:
    :param _n:
    :return:
    """
    if len(_y) > 1:
        log_mask = _y > 100
    else:
        return []

    _shape = np.shape(_y)
    _x = np.array([i for i in range(_n)]) * delta_t
    _x_ext = np.zeros(_shape)

    mx = np.zeros(_shape[1])
    my = np.zeros(_shape[1])
    a2 = np.zeros(_shape[1])
    a11 = np.zeros(_shape[1])
    _n = np.ones(_shape[1])
    for i in range(_shape[1]):
        _x_ext[:, i] = _x * log_mask[:, i]
        _n[i] = int(np.sum(log_mask[:, i]))
        mx[i] = np.average(_x[log_mask[:, i]])
        my[i] = np.average(_y[log_mask[:, i], i])
        a2[i] = np.dot(_x_ext[:, i].T, _x_ext[:, i]) / _n[i]
        a11[i] = np.dot(_x_ext[:, i].T, _y[:, i]) / _n[i]

    kk = (a11 - mx * my) / (a2 - mx ** 2)
    bb = my - kk * mx
    ff = np.copy(_y)

    log = np.logical_not(log_mask)
    for i in range(_shape[1]):
        ff[:, i] = kk[i] * _x_ext[:, i] + bb[i] * log_mask[:, i] + 2 * log[:, i]

    delta_ff = ff - _y
    res = np.multiply(delta_ff, delta_ff)
    _rms = res ** 0.5
    _relative_rms = np.sum(_rms, axis=0) / _n / my

    return _relative_rms


def zone_deletion(_len):
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
    k9 = int(2000 / _delta_f)
    _k = [0, k1, k2, k3, k4, k5, k6, k7, k8, k9]
    _zone = [i for i in range(_k[0], _k[1])] + [i for i in range(_k[2], _k[3])] + [i for i in range(_k[4], _k[5])] + \
        [i for i in range(_k[6], _k[7])] + [i for i in range(_k[8], _k[9])]

    return _zone


def pic_name(file_path, flag, format='png'):
    """
    Функция принимает папку в которую надо сохранить картинки с общим названием, различающиеся по номеру
    после названия, флаг, который определяет название, формат сохранения. Возвращает название файла с номером
    и расширением
    :param file_path: папка сохранения
    :param flag: название картинки (вид картинки: скан, спектр и т.п.)
    :param format: формат сохранения
    :return: имя файла с расширением, под которым будет сохранен рисунок
    """
    if flag == 1:
        add_pass0 = 'spectrum_00'
    elif flag == 2:
        add_pass0 = 'colour2D_00'
    elif flag == 3:
        add_pass0 = 'pic3D_00'
    elif flag == 4:
        add_pass0 = 'rms_00'
    else:
        add_pass0 = 'scan_stokes_00'

    l = len(add_pass0)
    add_pass1 = add_pass0 + '.' + format
    if not os.path.isfile(Path(file_path, add_pass1)):
        pass
    else:
        while os.path.isfile(Path(file_path, add_pass1)):
            num = int(add_pass0[l - 2:l]) + 1
            num_str = str(num)
            if num >= 10:
                add_pass0 = add_pass0[:l - 2] + num_str
            else:
                add_pass0 = add_pass0[:l - 2] + '0' + num_str
            add_pass1 = add_pass0 + '.' + format

    return add_pass1


def some_fig_plot(_data1, _data2, _path_to_fig_folder=None):
    """
    Функция принимает путь сохранения рисунка и три возможных последовательности для построения двух или трех
    графиков с общим аргументом mean_frame_ind_pol, который приходит от вызывающей функции. При этом наличие
    двух отображаемых на рис. последовательностей обязательно, последовательность _s_i присутствует всегда.
    :param _path_to_fig_folder: Путь к папке для сохранения рисунка
    :param _data1:
    :param _data2:
    :param _path_to_fig_folder:
    :param :
    :return:
    """
    _pic_name = pic_name(str(_path_to_fig_folder), 4, 'png')
    _path_to_pic = Path(_path_to_fig_folder, _pic_name)
    _fig_folder = str(_path_to_fig_folder)
    # title1, title2, title3 = title_func(_fig_folder, head)
    if len(_data1) > 1:
        _m = int(np.shape(_data1)[0])
    else:
        _m = int(np.shape(_data2)[0])
    _df = 2000 / _m
    _freq = [1000 + _df / 2 + _df * _i for _i in range(_m)]
    _l = 1
    _r_level = np.sqrt(2 / _df / delta_t / 1e6)     # 1е6 переход к Гц от МГц
    fig, _axes = plt.subplots(_l, 1, figsize=(12, _l * 3))

    if len(_data1) > 1:
        _axes.semilogy(_freq, _data1)
    if len(_data2) > 1:
        _axes.semilogy(_freq, _data2)
    _axes.axhline(y=_r_level, xmin=0.05, xmax=0.95, color='r')
    _axes.grid()
    _axes.grid(which='minor',
               axis='x',
               color='k',
               linestyle=':')
    # _sum = np.sum(_data, axis=0)
    # _axes[_l - 1].plot(_freq, _sum)

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
    fig.savefig(_path_to_pic)
    flag_save = save_question()
    if flag_save == 'no':
        if os.path.isfile(_path_to_pic):
            os.remove(_path_to_pic)
            print('Picture is not saved')
        else:
            print('File not found')
    else:
        print('Picture is saved')
    pass


def save_question():
    root = Tk()
    answer = mb.askquestion(
        title="Save control",
        message="Save picture?")
    root.mainloop()
    del root
    return answer


def save_rms_base(_data_left, _data_right, _path_dir, _path_file):

    # **** Создание папки для хранения конвертированных данных, если ее нет ****
    if not os.path.isdir(Path(_path_dir)):
        os.mkdir(Path(_path_dir))
    # **************************************************************************

    _rms = pd.Series([_data_left, _data_right])
    np.save(str(_path_file) + '_rms', _rms)


if __name__ == '__main__':
    current_primary_file = '2022-12-23_04+14'
    current_primary_dir = '2022_12_23_3C273'
    main_dir = '2022'
    adr1 = DataPaths(current_primary_file, current_primary_dir, main_dir)
    converted_data_file_path = adr1.converted_data_file_path
    data_treatment_file_path = adr1.treatment_data_file_path
    converted_dir_path = adr1.converted_dir_path

    delta_t = 8.3886e-3
    delta_f = 7.8125

    path1 = Path(str(converted_data_file_path) + '_spectrum.npy')
    t1, t2 = 100, 110   # Интервал времени, на котором считаем rms скана

    if os.path.isfile(str(converted_data_file_path) + '_rms.npy'):
        rms_relative_left, rms_relative_right = np.load(str(converted_data_file_path) + '_rms.npy', allow_pickle=True)
    else:
        y_left, y_right, n, n_aver = data_preparing1(path1, t1, t2)
        rms_relative_left = relative_rms1(y_left, n)
        rms_relative_right = relative_rms1(y_right, n)
        save_rms_base(rms_relative_left, rms_relative_right, converted_dir_path, converted_data_file_path)

    #       *** Установка в зоне действия режекторных фильтров rms_relative = 1e-3 ***
    if len(rms_relative_left) > 1:
        dz = zone_deletion(len(rms_relative_left))
        rms_relative_left[dz] = 1e-3
    if len(rms_relative_right) > 1:
        dz = zone_deletion(len(rms_relative_right))
        rms_relative_right[dz] = 1e-3


    some_fig_plot(rms_relative_left, rms_relative_right, data_treatment_file_path)

    pass
