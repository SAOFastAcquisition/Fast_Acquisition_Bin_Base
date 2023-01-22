import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import os
from pathlib import Path
from Help_folder.paths_via_class import path_to_data


if __name__ == '__main__':

    receiver_temperature_file_name = 'receiver_temperature1.npy'
    receiver_temperature_interpol_name = 'receiver_temp_interpol.npy'
    head_path = path_to_data()
    receiver_temperature_path = Path(head_path, 'Alignment', receiver_temperature_file_name)
    receiver_temperature_interpol_path = Path(head_path, 'Alignment', receiver_temperature_interpol_name)
    with open(receiver_temperature_path, 'rb') as inp:
        data = pickle.load(inp)

    len_series = len(data['temperature'])
    temp_arr = data['temperature'][0]                   # Массив температур
    for i in range(1, len_series):
        temp_arr = np.vstack([temp_arr, data['temperature'][i]])
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
    arr_av = np.polyfit(x_s, temp_aver_s, 25)
    temp_interp = np.polyval(arr_av, x_new)

    # plt.plot(x_new, y_calc0, label='Polynomial0')
    # plt.plot(x_new, y_calc1, label='Polynomial1')
    # plt.plot(x_new, y_calc2, label='Polynomial2')
    # plt.plot(x_new, y_calc3, label='Polynomial3')
    plt.plot(x_new, temp_interp, label='Polynomial_av')
    plt.scatter(x_s, y_s[0, :], label='data')
    plt.scatter(x_s, y_s[1, :], label='data')
    plt.scatter(x_s, y_s[2, :], label='data')
    plt.scatter(x_s, y_s[3, :], label='data')
    # plt.scatter(x_s, temp_aver_s, label='data')
    plt.legend()
    plt.show()
    column = ['case_id', 'polyval_coeff', 'note']
    _s = pd.Series(('01', arr_av, '2022_11_18-24'), index=column)
    if not os.path.isfile(receiver_temperature_interpol_path):
        _data = pd.DataFrame(columns=column)
    else:
        with open(receiver_temperature_interpol_path, 'rb') as inp:
            _data = pickle.load(inp)

    idx_r = _data.loc[(_data['case_id'] == _s['case_id'])].index
    if not len(idx_r):
        _data = _data.append(_s, ignore_index=True)
        with open(receiver_temperature_interpol_path, 'wb') as _out:
            pickle.dump(_data, _out)
