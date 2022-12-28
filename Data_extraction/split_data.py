import numpy as np
import pandas as pd
import pickle
from pathlib import Path
from Supporting_func.stocks_coefficients import path_to_data

current_data_dir = '2022'
converted_data_dir = 'Converted_data'  # Каталог для записи результатов конвертации данных и заголовков

current_primary_dir = '2022_12_22sun'
current_converted_dir = current_primary_dir + '_conv'
current_converted_path = Path(converted_data_dir, current_converted_dir)

current_primary_file = '2022-12-22_01+08-08'
csv_file = 'time_label'+current_primary_file[13:]
converted_data_file_path, head_path = path_to_data(current_data_dir, current_converted_path)

path1 = Path(converted_data_file_path, current_primary_file + '_spectrum.npy')
path_to_csv = Path(converted_data_file_path, csv_file + '.txt')

delta_t = 8.3886e-3

#           **** Считывание данных для разделения большого файла ****
csv = pd.read_csv(path_to_csv, delimiter=',')   # Перечень времен начала записей и их идентификаторов
time_start = [s for s in csv['time']]
# time_stop = [s + 460 for s in time_start]
time_stop = [s for s in csv['time_stop']]
num_start = [int(s / delta_t) for s in time_start]
num_stop = [int(s / delta_t) for s in time_stop]
output_filename_list = [current_primary_file[:10] + s for s in csv['azimuth']]
len_list = len(num_stop)

spectrum = np.load(path1, allow_pickle=True)
with open(Path(converted_data_file_path, current_primary_file + '_head.bin'), 'rb') as inp:
    head = pickle.load(inp)

for i in range(len_list):
    k = 0
    s1 = pd.Series([]*4)
    for s in spectrum:
        if len(s) > 1:
            s1[k] = s[num_start[i]:num_stop[i], :]
        else:
            s1[k] = []
        k += 1
    path_current = Path(converted_data_file_path, output_filename_list[i] + '_spectrum.npy')
    np.save(path_current, s1)
    with open(Path(converted_data_file_path, output_filename_list[i] + '_head.bin'), 'wb') as out:
        pickle.dump(head, out)
    print('Load file:', output_filename_list[i])
    pass
