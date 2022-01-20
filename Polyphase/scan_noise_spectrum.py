import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftshift
from pathlib import Path
import pickle
from Supporting_func import Fig_plot as fp, align_spectrum, path_to_data


current_dir = Path.cwd()
home_dir = Path.home()
current_catalog = r'2022\Results'  # Текущий каталог (за определенный период, здесь - год)
current_data_dir = '2022_01_20test'  # Папка с текущими даннымиH:\Fast_Acquisition\2021\Results\2021_09_22test\2021-09-22_01_14bit_pm20
current_data_file = '2022-01-20_02test'  # Имя файла с исходными текущими данными без расширенияH:\Fast_Acquisition\2022\Results\2022_01_20test
file_path_data, head_path = path_to_data(current_catalog, current_data_dir)

# ************** Загрузка матрицы спектров и установок (head) *************2022-01-20_01test_head
with open(Path(file_path_data, current_data_file + '_head.bin'), 'rb') as inp:
    head = pickle.load(inp)
spectrum = np.load(Path(file_path_data, current_data_file + '_spectrum.npy'), allow_pickle=True)
#               **********************************************
spectrum0 = spectrum[0]
spectrum1 = spectrum[1]
spectrum2 = spectrum[2]
spectrum3 = spectrum[3]

signal_rand_sample = np.reshape(spectrum1, (8192, -1))    # Разбиваем на реализации длиной 1024 отсчетов
spectrum_signal_rand = fft(signal_rand_sample, 8192)
spectrum_signal_av = np.average(np.abs(spectrum_signal_rand ** 2), 0)
# spectrum_signal_rand = fft(spectrum1[:, 100], 1024)
# spectrum_signal_av = np.abs(spectrum_signal_rand ** 2)

axes = plt.subplots()
plt.plot(spectrum_signal_av)
plt.show()
pass