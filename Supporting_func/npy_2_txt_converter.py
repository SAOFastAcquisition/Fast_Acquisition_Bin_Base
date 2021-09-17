import numpy as np
import pandas as pd
import os
import sys
from pathlib import Path
from Supporting_func.stocks_coefficients import path_to_data

current_dir = Path.cwd()
home_dir = Path.home()

# from Supporting_func.afc_alignment1 import align_func1

sys.path.insert(0, Path(current_dir, 'Supporting_func'))

current_data_file = '2021-06-28_22-30'      # Имя файла с исходными текущими данными без расширения
current_data_dir = '2021_06_28sun'          # Папка с текущими данными
current_catalog = r'2021/Results'           # Текущий каталог (за определенный период, здесь - год)

file_path_data, head_path = path_to_data(current_catalog, current_data_dir)



