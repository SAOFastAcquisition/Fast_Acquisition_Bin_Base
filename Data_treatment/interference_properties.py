import numpy as np
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

    mx = x.sum() / N    # 1
    my = y.sum() / N
    a2 = np.dot(x.T, x) / N
    a11 = np.dot(x.T, y) / N

    kk = (a11 - mx * my) / (a2 - mx ** 2)
    bb = my - kk * mx
    pass