import pandas as pd
import os
import numpy as np
from matplotlib import pyplot as plt
import sys
from scipy import stats
from pathlib import Path
from Supporting_func import Fig_plot as fp, path_to_data


current_primary_file = '2022-01-27_11'
current_primary_dir = '2022_01_27calibr'

current_data_dir = '2022'
converted_data_dir = 'Converted_data'       # Каталог для записи результатов конвертации данных и заголовков
data_treatment_dir = 'Data_treatment'       # Каталог для записи результатов обработки, рисунков

current_converted_dir = current_primary_dir + '_conv'
current_treatment_dir = current_primary_dir + '_treat'
current_converted_path = Path(converted_data_dir, current_converted_dir)
current_treatment_path = Path(data_treatment_dir, current_treatment_dir)
converted_data_file_path, head_path = path_to_data(current_data_dir, current_converted_path)
data_treatment_file_path, head_path = path_to_data(current_data_dir, current_treatment_path)

path_mesh1 = Path(data_treatment_file_path, 'Meshed_Spectrum', current_primary_file + '_meshed' + '.npy')

spectrum_mesh = spectrum = np.load(path_mesh1, allow_pickle=True)
delta_t = 8.3886e-3
delta_f = 7.8125
l, m = np.shape(spectrum_mesh)
# Среднее, среднестатистическое отклонение для радиометра с ГШ 6000К на входе
t_start_extNG = 5
t_stop_extNG = 10
n_start_extNG = int(t_start_extNG // delta_t)
n_stop_extNG = int(t_stop_extNG // delta_t)
a0 = np.mean(spectrum_mesh[n_start_extNG:n_stop_extNG, :], axis=0)
b0 = np.std(spectrum_mesh[n_start_extNG:n_stop_extNG, :], axis=0)
dev0 = b0 / a0

# Среднее, среднестатистическое отклонение для радиометра с согласованной нагрузкой на входе
t_start_ML = 20
t_stop_ML = 25
n_start_ML = int(t_start_ML // delta_t)
n_stop_ML = int(t_stop_ML // delta_t)
a1 = np.mean(spectrum_mesh[n_start_ML:n_stop_ML, :], axis=0)
b1 = np.std(spectrum_mesh[n_start_ML:n_stop_ML, :], axis=0)
dev1 = b1 / a1

# Среднее, среднестатистическое отклонение для радиометра с внутренним ГШ на входе
t_start_intNG = 35
t_stop_intNG = 40
n_start_intNG = int(t_start_intNG // delta_t)
n_stop_intNG = int(t_stop_intNG // delta_t)
a2 = np.mean(spectrum_mesh[n_start_intNG:n_stop_intNG, :], axis=0)
b2 = np.std(spectrum_mesh[n_start_intNG:n_stop_intNG, :], axis=0)
dev2 = b2 / a2

arg = [1000 + delta_f * i for i in range(256)]
fig = plt.figure(figsize=(12, 8))
axes = fig.add_subplot()
axes.grid(b=True, which='major', color='#666666', linestyle='-')
# Show the minor grid lines with very faint and almost transparent grey lines
axes.minorticks_on()
axes.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)

axes.set_title('STD divided by MEAN', fontsize=18)
axes.set_yscale('log')
# axes.set_xscale('log')
axes.set_xlabel('Freq, MHz', fontsize=18)
axes.set_ylabel('STD divided by MEAN', fontsize=20)
line_legend = ['External NG', 'Matched Load', 'Intrinsic NG']
axes.plot(arg, dev0, label=line_legend[0])
axes.plot(arg, dev1, label=line_legend[0])
axes.plot(arg, dev2, label=line_legend[0])
axes.legend()
plt.show()
pass