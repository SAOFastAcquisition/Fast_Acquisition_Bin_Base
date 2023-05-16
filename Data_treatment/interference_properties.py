import numpy as np
import pickle
from pathlib import Path
import matplotlib.pyplot as plt
from Help_folder.paths_via_class import DataPaths


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


def data_preparing1(_path):
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
        dat_left = np.hstack((spectrum[0], spectrum[1]))
        dat_rihgt = np.hstack((spectrum[2], spectrum[3]))
    if len(spectrum[0]) > 1 and not (len(spectrum[2])) > 1:  # Левая поляризация
        dat_left = np.hstack((spectrum[0], spectrum[1]))
        dat_rihgt = []
    if not (len(spectrum[0])) > 1 and len(spectrum[1]) > 1:  # Правая поляризация
        dat_rihgt = np.hstack((spectrum[2], spectrum[3]))
        dat_left = []
    _n_aver = head['n_aver']

    return dat_left, dat_rihgt, _n, _n_aver


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
    if len(_y) > 1:
        log_mask = _y > 100
    _shape = np.shape(_y)
    _x = np.array([i for i in range(_n)]) * delta_t
    _x_ext = np.zeros(_shape)
    for i in range(_shape[1]):
        _x_ext[:, i] = _x * log_mask[:, i]
    # mx0 = (_n - 1) * (_n - 1) / 2 / _n * delta_t
    mx = np.zeros(_shape[1])
    my = np.zeros(_shape[1])
    _n = np.ones(_shape[1])
    for i in range(_shape[1]):
        _n[i] = int(np.sum(log_mask[:, i]))
        mx[i] = np.average(_x[log_mask[:, i]])
        my[i] = np.average(_y[log_mask[:, i], i])

    a2 = np.dot(_x_ext.T, _x_ext) / _n
    a11 = np.dot(_x_ext.T, _y) / _n

    kk = (a11 - mx * my) / (a2 - mx ** 2)
    bb = my - kk * mx
    ff = np.copy(_y)

    # ff = np.array([kk * z * delta_t + bb for z in range(_n)])
    for i in range(_shape[1]):
        ff[:, i] = np.array([kk[:, i] * x + bb[:, i] for x in _x_ext[:, i]])
    delta_ff = ff - _y

    res = np.multiply(delta_ff, delta_ff)
    _rms = res ** 0.5
    _relative_rms = np.sum(_rms, axis=0) / _n / my

    return _relative_rms


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
    y_left, y_right, n, n_aver = data_preparing1(path1)
    res_relative = relative_rms1(y_left, n)

    pass