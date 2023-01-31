import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import os
from pathlib import Path
from Help_folder.paths_via_class import path_to_data
# from scipy.interpolate import splrep, BSpline


def plot_RTI(_data, _x):
    _d = _data['polyval_coeff']
    _date = _data['note']
    _len = len(_d)
    _coeff_array = _d.iloc[0]
    for _i in range(1, _len):
        _coeff_array = np.vstack([_coeff_array, _d.iloc[_i]])
    for _i in range(_len):
        _temp_interp = np.polyval(_coeff_array[_i, :], _x)
        _legend = _date.iloc[_i]
        plt.plot(x_new, _temp_interp, label=_legend)
    plt.grid()
    plt.legend()
    plt.show()
    pass


def plot_NGIT():
    """
    Строит Noise Generator Temperature Interpolation вместе с иходными данными значений средней шумовой
    температуры сигнала ГШ, впрыскиваемого в каналы левой и правой поляризаций
    :return:
    """
    plt.plot(x_new, temp_interp_left, label=f'Polynomial interpolation, rank = {polynomial_rank}, polar left')
    plt.plot(x_new, temp_interp_right, label=f'Polynomial interpolation, rank = {polynomial_rank}, polar right')
    plt.scatter(x_s, y_s_left, label='data polar = left')
    plt.scatter(x_s, y_s_right, label='data polar = right')
    plt.legend()
    plt.grid()
    plt.show()


def update_ng_temp(_series):
    if not os.path.isfile(ng_temperature_interpol_path):
        data_saved = pd.DataFrame(columns=column)
    else:
        with open(ng_temperature_interpol_path, 'rb') as inp:
            data_saved = pickle.load(inp)

    idx_r = data_saved.loc[((data_saved['case_id'] == _series['case_id'])
                            & (data_saved['polar'] == _series['polar']))].index
    if not len(idx_r):
        data_saved = data_saved.append(_series, ignore_index=True)
        with open(ng_temperature_interpol_path, 'wb') as _out:
            pickle.dump(data_saved, _out)


if __name__ == '__main__':
    ng_temperature_file_name = 'ngi_temperature1.npy'
    ng_temperature_interpol_name = 'ngi_temp_interpol.npy'
    head_path = path_to_data()
    ng_temperature_path = Path(head_path, 'Alignment', ng_temperature_file_name)
    ng_temperature_interpol_path = Path(head_path, 'Alignment', ng_temperature_interpol_name)

    case_id = '03'
    note = '2022_11_18-24'
    polynomial_rank = 80

    with open(ng_temperature_path, 'rb') as inp:
        data = pickle.load(inp)
    temp_left = data['temperature'][data['case_id'] == case_id][data['polar'] == 'left'].iloc[0]
    temp_right = data['temperature'][data['case_id'] == case_id][data['polar'] == 'right'].iloc[0]

    log_y = [not s == 0 for s in temp_left]     # Лог, по которому вырезаются отсчеты в зоне действия фильтров
    y_s_left = temp_left[log_y]                 # Из массива температур удалили отсчеты зоны действия фильтров
    y_s_right = temp_right[log_y]

    len_data = len(temp_left)
    df = 2000 / len_data
    x = np.array([1000 + df / 2 + i * df for i in range(len_data)])
    x_s = x[log_y]                              # Из значений частот удалили отсчеты зоны действия фильтров

    x_new = np.arange(1020, 2970, 1)            # Вектор-аргумент для расчета интерполированной кривой
    #                                       ******************
    #                   ***** Расчет коэффициентов интерполирующего полинома *****
    #                                       ******************
    arr_av_left = np.polyfit(x_s, y_s_left, polynomial_rank)
    arr_av_right = np.polyfit(x_s, y_s_right, polynomial_rank)

    temp_interp_left = np.polyval(arr_av_left, x_new)
    temp_interp_right = np.polyval(arr_av_right, x_new)
    #                                       ******************
    #                        ***** Запись коэффициентов полинома в файл *****
    #                                       ******************
    column = ['case_id', 'polyval_coeff', 'polar', 'note']
    series_to_data_left = pd.Series((case_id, arr_av_left, 'left', note), index=column)
    update_ng_temp(series_to_data_left)
    series_to_data_right = pd.Series((case_id, arr_av_right, 'right', note), index=column)
    update_ng_temp(series_to_data_right)
    #                                       ******************
    #                               ***** Визуализация результатов *****
    plot_NGIT()
    pass
