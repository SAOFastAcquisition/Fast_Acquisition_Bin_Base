import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle, matplotlib
import os
from pathlib import Path
from Help_folder.paths_via_class import path_to_data


def plot_RTI(_data, _x):
    """
    Отображает результаты Receiver Temperature Interpolation по наборам интерполяционных коэффициентов, полученных
    для шумовой температуры приемника по результатам серий измерений в разные дни.
    :param _data:   фрейм, включающий идентификатор case_id коэффициентов 'polyval_coeff', расчитанных np.polyfit( )
                    для различных наборов измеренных шумовых температур приемника
    :param _x:  вектор аргументов
    :return:    рисунок набора аппроксимационных кривых
    """
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


def plot_rnt():
    fig, ax = plt.subplots(1, figsize=(12, 6))
    ax.plot(x_new, temp_interp, label=f'Polynomial interpolation, rank = {polynomial_rank}')
    for _i in range(ind_len):
        ax.scatter(x_s, y_s[_i, :], label=f'data {_i}')
    ax.scatter(x_s, temp_aver_s, label='data averaged')
    ax.set_title('Receiver Noise Temperature', fontsize=20)
    ax.set_xlabel('Frequency, MHz', fontsize=18)
    ax.set_ylabel('Temperature, K', fontsize=18)
    ax.legend(fontsize=10)
    plt.grid()
    plt.show()


def plot_rnt_ab():
    matplotlib.rcParams['font.family'] = 'serif'
    matplotlib.rcParams['font.size'] = '10'
    fig, ax = plt.subplots(1, figsize=(6.4, 4))
    ax.plot(x_new, temp_interp, color='black', linewidth=1, linestyle='-')
    ax.tick_params(axis='both',  # Применяем параметры к обеим осям
                   which='major',  # Применяем параметры к основным делениям
                   direction='in',  # Рисуем деления внутри и снаружи графика
                   length=5,  # Длинна делений
                   width=1,  # Ширина делений
                   color='black',  # Цвет делений
                   pad=2,  # Расстояние между черточкой и ее подписью
                   # labelsize=_f_size,  # Размер подписи
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
    for _i in range(ind_len):
        ax.scatter(x_s, y_s[_i, :], 4, 'black', marker=matplotlib.markers.CARETDOWNBASE)
    # ax.scatter(x_s, temp_aver_s, label='data averaged')
    # ax.set_title('Receiver Noise Temperature', fontsize=10)
    ax.set_xlabel('Frequency, MHz')
    ax.set_ylabel('Temperature, K')
    # ax.legend(fontsize=10)
    # plt.grid()
    plt.tight_layout(pad=0.1, w_pad=0.1, h_pad=1.0)
    plt.show()


if __name__ == '__main__':

    receiver_temperature_file_name = 'receiver_temperature1.npy'
    receiver_temperature_interpol_name = 'receiver_temp_interpol.npy'
    head_path = path_to_data()
    receiver_temperature_path = Path(head_path, 'Alignment', receiver_temperature_file_name)
    receiver_temperature_interpol_path = Path(head_path, 'Alignment', receiver_temperature_interpol_name)
    case_id = '01'
    case_id_add = '02'
    with open(receiver_temperature_path, 'rb') as inp:
        data = pickle.load(inp)
    _ind1 = data[data.case_id == case_id].index
    if case_id_add:
        _ind1 = _ind1.append(data[data.case_id == case_id_add].index)
    data_sample = data.temperature.loc[_ind1]
    ind_len = len(_ind1)
    temp_arr = data_sample.iloc[0]
    for i in range(1, ind_len):
        temp_arr = np.vstack([temp_arr, data_sample.iloc[i]])
    temp_aver = np.mean(temp_arr, axis=0)  # Средняя температура по результатам измерений
    log_y = [not s == 100 for s in temp_arr[0, :]]  # Лог, по которому вырезаются отсчеты в зоне действия фильтров
    y_s = temp_arr[:, log_y]  # Из массива температур удалили отсчеты зоны действия фильтров
    temp_aver_s = temp_aver[log_y]
    len_data = len(temp_arr[0, :])
    df = 2000 / len_data
    x = np.array([1000 + df / 2 + i * df for i in range(len_data)])
    x_s = x[log_y]  # Из значений частот удалили отсчеты зоны действия фильтров

    x_new = np.arange(1000, 3000, 1)  # Вектор-аргумент для расчета интерполированной кривой

    polynomial_rank = 11
    arr_av = np.polyfit(x_s, temp_aver_s, polynomial_rank)
    temp_interp = np.polyval(arr_av, x_new)
    plot_rnt_ab()

    column = ['case_id', 'polyval_coeff', 'note']
    series_to_data = pd.Series((case_id, arr_av, '2022_06_27-28'), index=column)
    if not os.path.isfile(receiver_temperature_interpol_path):
        data_saved = pd.DataFrame(columns=column)
    else:
        with open(receiver_temperature_interpol_path, 'rb') as inp:
            data_saved = pickle.load(inp)
    # data_saved.drop(index=[3], axis=0, inplace=True)
    # plot_RTI(data_saved, x_new)
    idx_r = data_saved.loc[(data_saved['case_id'] == series_to_data['case_id'])].index
    if not len(idx_r):
        data_saved = data_saved.append(series_to_data, ignore_index=True)
        with open(receiver_temperature_interpol_path, 'wb') as _out:
            pickle.dump(data_saved, _out)
