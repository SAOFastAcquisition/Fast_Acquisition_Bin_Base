import numpy as np
import matplotlib.pyplot as plt
import pickle
from pathlib import Path
from Help_folder.paths_via_class import path_to_data


if __name__ == '__main__':

    receiver_temperature_file_name = 'receiver_temperature1.npy'
    head_path = path_to_data()
    receiver_temperature_path = Path(head_path, 'Alignment', receiver_temperature_file_name)
    with open(receiver_temperature_path, 'rb') as inp:
        data = pickle.load(inp)
    y0 = data['temperature'][0]
    y1 = data['temperature'][1]
    y2 = data['temperature'][2]
    y3 = data['temperature'][3]
    temp_arr = data['temperature'][0]       # Массив температур
    for i in range(1, 4):
        temp_arr = np.vstack([temp_arr, data['temperature'][i]])
    temp_aver = np.mean(temp_arr, axis=0)   # Средняя температура по результатам измерений
    log_y = [not s == 100 for s in y0]      # Лог, по которому вырезаются отсчеты в зоне действия фильтров
    y_s0 = y0[log_y]
    y_s1 = y1[log_y]
    y_s2 = y2[log_y]
    y_s3 = y3[log_y]
    temp_aver_s = temp_aver[log_y]
    len_data = len(y0)
    df = 2000 / len_data
    x = np.array([1000 + df / 2 + i * df for i in range(len_data)])
    x_s = x[log_y]

    x_new = np.arange(1000, 3000, 1)

    # arr0 = np.polyfit(x_s, y_s0, 25)
    # y_calc0 = np.polyval(arr0, x_new)
    # arr1 = np.polyfit(x_s, y_s1, 25)
    # y_calc1 = np.polyval(arr1, x_new)
    # arr2 = np.polyfit(x_s, y_s2, 25)
    # y_calc2 = np.polyval(arr2, x_new)
    # arr3 = np.polyfit(x_s, y_s3, 25)
    # y_calc3 = np.polyval(arr3, x_new)
    arr_av = np.polyfit(x_s, temp_aver_s, 25)
    temp_interp = np.polyval(arr_av, x_new)

    pass
    # plt.plot(x_new, y_calc0, label='Polynomial0')
    # plt.plot(x_new, y_calc1, label='Polynomial1')
    # plt.plot(x_new, y_calc2, label='Polynomial2')
    # plt.plot(x_new, y_calc3, label='Polynomial3')
    plt.plot(x_new, temp_interp, label='Polynomial_av')
    plt.scatter(x_s, y_s0, label='data')
    plt.scatter(x_s, y_s1, label='data')
    plt.scatter(x_s, y_s2, label='data')
    plt.scatter(x_s, y_s3, label='data')
    # plt.scatter(x_s, temp_aver_s, label='data')
    plt.legend()
    plt.show()
