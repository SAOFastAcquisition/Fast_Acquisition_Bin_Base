"""

"""
import pandas
from pathlib import Path
import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.signal import lfilter, filtfilt, butter, freqs, remez, freqz
from Supporting_func.stokes_coefficients import path_to_data

current_catalog = Path.cwd()
sys.path.insert(0, current_catalog)


def zone_deletion(_s, _freq_quantity):
    _df = 2000 / _freq_quantity
    # Исключение зон действия режекторных фильтров
    k1 = int((1090 - _df / 2) // _df)
    k2 = int((1220 - _df / 2) // _df)
    k3 = int((1520 - _df / 2) // _df)
    k4 = int((1700 - - _df / 2) // _df)
    k6 = int(_freq_quantity / 2)
    k5 = int((770 - _df / 2) // _df)

    _s[0:5] = 1
    _s[k6:k6 + 8] = 1
    _s[-11:] = 1
    _s[k1:k2] = 1
    _s[k3:k4] = 1
    _s[k5:k6] = 1

    return _s


def del_random_mod(_s, _s0):
    _l = len(_s)
    for _i in range(1, _l - 2):
        if abs(_s[_i] - _s[_i - 1]) > 2 * abs(_s[_i + 1] - _s[_i - 1]):
            _s[_i] = (_s[_i - 1] + _s[_i + 1]) / 2
        if _s[_i] <= 0.01:
            _s[_i] = _s0
    _s[0] = _s0
    _s[_l - 2:] = _s0
    return _s


def maf_fir(_s, m=2):
    """
        Calculate moving average filter as FIR

        Parameters
        ----------
        m : int
            moving average step
        _s : int
            data
        """

    return filtfilt(np.ones(m - 1), 1, _s) / (m - 1) / (m - 1)


if __name__ == '__main__':
    current_data_dir = '2022'
    file_path_data, head_path = path_to_data(current_catalog, current_data_dir)
    #                   ***********************************
    # Путь к файлу хранения коэффициентов (он один для всех калибровок АЧХ каналов)
    folder_align_path = Path(head_path, 'Alignment')
    polar_coeff_file_name = 'Align_polar_2022_06_18.csv'
    _path_to_csv = Path(folder_align_path, polar_coeff_file_name)

    csv = pandas.read_csv(_path_to_csv, header=None, delimiter=',')
    freq_quantity = csv.shape[1]
    df = 2000 / 512
    freq = [1000 + df / 2 + i * df for i in range(512)]
    a = csv.drop(labels=[12], axis=0)
    b = a.mean()
    b = del_random_mod(b, 1)
    b = maf_fir(b, 7)
    b1 = zone_deletion(b, freq_quantity)

    fig, _ax = plt.subplots(1, figsize=(12, 6))
    _ax.plot(freq, b1)
    # _ax.plot(_data['temperature'][1])
    # _ax.plot(_data['temperature'][2])
    # _ax.plot(_data['temperature'][3])
    plt.grid()
    plt.show()
    pass
    pass
