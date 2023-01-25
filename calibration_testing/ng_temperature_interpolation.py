import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import os
from pathlib import Path
from Help_folder.paths_via_class import path_to_data


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


if __name__ == '__main__':

    ng_temperature_file_name = 'ngi_temperature1.npy'
    ng_temperature_interpol_name = 'ngi_temp_interpol.npy'
    head_path = path_to_data()
    ng_temperature_path = Path(head_path, 'Alignment', ng_temperature_file_name)
    ng_temperature_interpol_path = Path(head_path, 'Alignment', ng_temperature_interpol_name)
    case_id = '03'
    case_id_add = '03'
    with open(ng_temperature_path, 'rb') as inp:
        data = pickle.load(inp)
    data_sample = data['temperature'][data['case_id'] == case_id]    # Массив температур [data['case_id'] == case_id]
    ind_data = data_sample.index
    ind_len = len(ind_data)
    temp_arr = data_sample.iloc[0]
    # ind_temp.pop(0)
    for i in range(1, ind_len):
        temp_arr = np.vstack([temp_arr, data_sample.iloc[i]])
    temp_aver = np.mean(temp_arr, axis=0)               # Средняя температура по результатам измерений
    log_y = [not s == 100 for s in temp_arr[0, :]]      # Лог, по которому вырезаются отсчеты в зоне действия фильтров
    y_s = temp_arr[:, log_y]                            # Из массива температур удалили отсчеты зоны действия фильтров
    temp_aver_s = temp_aver[log_y]
    len_data = len(temp_arr[0, :])
    df = 2000 / len_data
    x = np.array([1000 + df / 2 + i * df for i in range(len_data)])
    x_s = x[log_y]                                      # Из значений частот удалили отсчеты зоны действия фильтров

    x_new = np.arange(1000, 3000, 1)                    # Вектор-аргумент для расчета интерполированной кривой

    # arr0 = np.polyfit(x_s, y_s[0, :], 25)
    # y_calc0 = np.polyval(arr0, x_new)
    # arr1 = np.polyfit(x_s, y_s[1, :], 25)
    # y_calc1 = np.polyval(arr1, x_new)
    # arr2 = np.polyfit(x_s, y_s[2, :], 25)
    # y_calc2 = np.polyval(arr2, x_new)
    # arr3 = np.polyfit(x_s, y_s[3, :], 25)
    # y_calc3 = np.polyval(arr3, x_new)
    arr_av = np.polyfit(x_s, temp_aver_s, 11)
    temp_interp = np.polyval(arr_av, x_new)

    plt.plot(x_new, temp_interp, label='Polynomial interpolation, rank = 15')
    for i in range(ind_len):
        plt.scatter(x_s, y_s[i, :], label='data')
    plt.scatter(x_s, temp_aver_s, label='data')
    plt.legend()
    plt.grid()
    plt.show()
    column = ['case_id', 'polyval_coeff', 'note']
    series_to_data = pd.Series((case_id, arr_av, '2022_11_18-24'), index=column)
    if not os.path.isfile(ng_temperature_interpol_path):
        data_saved = pd.DataFrame(columns=column)
    else:
        with open(ng_temperature_interpol_path, 'rb') as inp:
            data_saved = pickle.load(inp)
    plot_RTI(data_saved, x_new)
    idx_r = data_saved.loc[(data_saved['case_id'] == series_to_data['case_id'])].index
    if not len(idx_r):
        data_saved = data_saved.append(series_to_data, ignore_index=True)
        with open(ng_temperature_interpol_path, 'wb') as _out:
            pickle.dump(data_saved, _out)
