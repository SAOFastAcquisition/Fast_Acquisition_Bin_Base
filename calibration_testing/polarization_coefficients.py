import pandas
from pathlib import Path
import sys
from Supporting_func.stocks_coefficients import path_to_data


current_catalog = Path.cwd()
sys.path.insert(0, current_catalog)


if __name__ == '__main__':
    current_data_dir = '2022'
    file_path_data, head_path = path_to_data(current_catalog, current_data_dir)
    #                   ***********************************
    # Путь к файлу хранения коэффициентов (он один для всех калибровок АЧХ каналов)
    folder_align_path = Path(head_path, 'Alignment')
    polar_coeff_file_name = 'Align_polar_2022_06_18.csv'
    _path_to_csv = Path(folder_align_path, polar_coeff_file_name)

    csv = pandas.read_csv(_path_to_csv, header=None, delimiter=',')
    df = 2000 / 512
    freq = [1000 + df / 2 + i * df for i in range(512)]
    a = csv.drop(labels=[12], axis=0)
    b = a.mean()
    pass