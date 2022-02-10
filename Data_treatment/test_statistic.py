import pandas as pd
import os
import numpy as np
import sys
from scipy import stats
from pathlib import Path
from Supporting_func import Fig_plot as fp, path_to_data


current_primary_file = '2022-01-27_14'
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

pass