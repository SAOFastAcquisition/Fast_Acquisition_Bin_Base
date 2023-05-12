import numpy as np
import pickle
from pathlib import Path
import matplotlib.pyplot as plt
from Help_folder.paths_via_class import DataPaths


if __name__ == '__main__':
    current_primary_file = '2023-02-10_01+20'
    current_primary_dir = '2023_02_10sun'
    main_dir = '2023'
    adr1 = DataPaths(current_primary_file, current_primary_dir, main_dir)
    converted_data_file_path = adr1.converted_data_file_path
    data_treatment_file_path = adr1.treatment_data_file_path

    delta_t = 8.3886e-3
    delta_f = 7.8125

    path1 = Path(str(converted_data_file_path) + '_spectrum.npy')
    spectrum = np.load(path1, allow_pickle=True)
    with open(Path(str(converted_data_file_path) + '_head.bin'), 'rb') as inp:
        head = pickle.load(inp)
    n1, n2 = int(6 / delta_t), int(16 / delta_t)
    n = n2 - n1

    # Приведение порядка следования отсчетов по частоте к нормальному
    if np.size(spectrum[0]):
        N_row = np.shape(spectrum[0])[0]
        for i in range(N_row):
            spectrum[0][i][0:] = spectrum[0][i][-1::-1]
    if np.size(spectrum[2]):
        N_row = np.shape(spectrum[2])[0]
        for i in range(N_row):
            spectrum[2][i][0:] = spectrum[2][i][-1::-1]

    if len(spectrum[0]) > 1 and len(spectrum[2]) > 1:       # Обе поляризации
        dat = np.hstack((spectrum[0] + spectrum[1], spectrum[2] + spectrum[3]))
    if len(spectrum[0]) > 1 and not(len(spectrum[1])) > 1:  # Левая поляризация
        dat = np.hstack((spectrum[0], spectrum[2]))
    if not(len(spectrum[0])) > 1 and len(spectrum[1]) > 1:  # Правая поляризация
        dat = np.hstack((spectrum[1], spectrum[3]))

    y = dat[n1:n2, :]
    n_aver = head['n_aver']

    x = np.array([i for i in range(n)]) * delta_t
    mx = (n-1) * (n-1) / 2 / n * delta_t    # 1
    my = np.sum(y, axis=0) / n
    a2 = np.dot(x.T, x) / n
    a11 = np.dot(x.T, y) / n

    kk = (a11 - mx * my) / (a2 - mx ** 2)
    bb = my - kk * mx
    pass