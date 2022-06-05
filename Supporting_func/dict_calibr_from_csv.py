import pandas
from pathlib import Path
import sys

current_dir = Path.cwd()
sys.path.insert(0, current_dir)


def path_to_data(current_catalog_in, current_data_dir_in):
    """ Функция принимает текущий каталог данных (за год или период) и папку текущих данных (за выбранный день).
    Определяет путь к папке текущих данных на конкретной машине и к корню каталога. """
    head_path1 = Path(r'H:\Fast_Acquisition')  # Путь к каталогу данных для домашнего ноута
    head_path1a = Path(r'E:\Fast_Acquisition')  # Путь к каталогу данных для домашнего ноута
    head_path2 = Path(r'/media/anatoly/Samsung_T5/Fast_Acquisition')  # Путь к каталогу данных для рабочего компа
    head_path3 = Path(r'D:\Fast_acquisition')  # Путь к каталогу данных для ноута ВМ
    head_path4 = Path(r'J:\Fast_Acquisition')  # Путь к каталогу данных для notebook 'Khristina'

    if head_path1.is_dir():
        head_path_out = head_path1
    elif head_path1a.is_dir():
        head_path_out = head_path1a
    elif head_path2.is_dir():
        head_path_out = head_path2
    elif head_path3.is_dir():
        head_path_out = head_path3
    elif head_path4.is_dir():
        head_path_out = head_path4
    else:
        raise CustomError('Path to data is not found!')

    file_path_data_out = Path(head_path_out, current_catalog_in, current_data_dir_in)
    return file_path_data_out, head_path_out


def calibration_temp(f):
    p0, p1, p2, p3 = [-3551, 6.9, -2.95e-3, 3.7e-7]
    temp = p0 + p1 * f + p2 * f * f + p3 * f ** 3
    return temp


def start_stop_calibr(*args):
    _current_file_name = args[0]
    _recording_id = _current_file_name[11:16]
    _path_to_csv = args[1]
    csv = pandas.read_csv(_path_to_csv, delimiter=',')
    idx = csv.loc[(csv.recording_id == _recording_id)].index
    n = idx[0]
    _start_stop = [csv.begin_start[idx][n], csv.begin_stop[idx][n], csv.end_start[idx][n], csv.end_stop[idx][n]]
    # s = csv.begin_start.iloc[idx[0]]
    return _start_stop


if __name__ == '__main__':
    current_data_file = '2022-04-29_08+00'  # Имя файла с исходными текущими данными без расширения
    current_primary_dir = '2022_04_29sun'
    current_data_dir = current_primary_dir + '_conv'  # Папка с текущими данными
    current_catalog = r'2022/Converted_data'  # Текущий каталог (за определенный период, здесь - год)

    dict_calibr_file_name = 'dict_calibr.csv'  # Имя файла с текущими коэффициентами выравнивания АЧХ
    file_path_data, head_path = path_to_data(current_catalog, current_data_dir)
    path_to_csv = Path(file_path_data, dict_calibr_file_name)
    start_stop = start_stop_calibr(current_data_file, path_to_csv)
    pass
