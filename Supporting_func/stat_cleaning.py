import numpy as np
# import os
import matplotlib.pyplot as plt

path_to_data = 'E:\\Measure_res\\2020_12_20sun\\20201220_pol2_06_+20_2_left'

data_input = np.int64(np.loadtxt(path_to_data + '.txt'))
data_shape = data_input.shape
max_data = np.max(data_input)
kurt_frame = np.array([])
for j in range(data_shape[1]):
    time_row = np.array(data_input[:, j])
    # time_row -= 2
    # time_row[2] = 70000000
    # Выбор индексов ненулевых значений скана
    time_ind = np.argwhere(time_row > 100)
    # Перевод массива в одномерный
    time_ind = np.ravel(time_ind)
    # Задание первого элемента для определения скачка индекса ненулевых членов
    if np.size(time_ind) == 0:
        i_prev = 0
    else:
        i_prev = time_ind[0]
    # Фрагмент скана с последовательно меняющимися индексами (допускается скачок индекса не более чем на 5)
    frame_ind0 = np.array([], 'int32')
    mean_frame_ind = np.array([], 'int32')
    frame_var = np.array([])
    frame_mean = np.array([])
    for i in time_ind:
        if i - i_prev < 3:
            frame_ind0 = np.append(frame_ind0, i)
        else:
            b = time_row[frame_ind0]
            mean = np.mean(b[b > 100])
            mean_ind0 = int(np.mean(frame_ind0))
            var = np.sqrt(np.var(b[b > 100]))
            try:
                b[np.abs(time_row[frame_ind0] - mean) > 3 * var] = 0
            except RuntimeWarning:
                pass
            time_row[frame_ind0] = b
            frame_ind0 = np.array([], 'int32')
            # frame_ind = np.append(frame_ind, i)
            frame_var = np.append(frame_var, var)
            frame_mean = np.append(frame_mean, mean)
            mean_frame_ind = np.append(mean_frame_ind, mean_ind0)
        pass
        i_prev =+ i

    var_var = np.sqrt(np.var(frame_var))
    mean_mean = np.mean(frame_mean)
    kurt = var_var / mean_mean
    kurt_frame = np.append(kurt_frame, kurt)

    plt.grid()
    plt.plot(mean_frame_ind, frame_mean)
    plt.show()

    fig, ax = plt.subplots(1, figsize=(12, 6))
    y_max = np.max(data_input[:, j])
    y_min = np.min(data_input[:, j])
    x_min = 0
    x_max = data_shape[0]
    argument = np.arange(0, data_shape[0])
    ax.plot(argument, data_input[:, j])
    ax.plot(argument, time_row)
    if kurt < 0.25:
        data_input[:, j] = time_row
    else:
        data_input[:, j] = 10
    ax.plot(argument, data_input[:, j] * 0.9)
    # if j == 11:
    plt.show()
    frame_var = np.array([])
    frame_mean = np.array([])
    pass

max_data1 = np.max(data_input)
print(max_data, max_data1)
pass
# np.savetxt(path_to_data + '1.txt', data_input)
# np.savetxt(path_to_data + '_kurt.txt', kurt_frame)
