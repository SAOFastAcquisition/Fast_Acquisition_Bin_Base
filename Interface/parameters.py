from pathlib import Path


def param_dict():
    param_dict_val = {'time_resolution': 8,
                      'freq_resolution': 1,
                      'freq_mask': [[1080, 1140, 1360, 1420, 1620, 1780, 1980],
                                    [2060, 2220, 2300, 2500, 2530, 2710, 2800, 2880, 2980],
                                    [1050, 1465, 1535, 1600, 1700, 2265, 2530, 2720, 2800, 2920]],
                      'time_mask': [[10, 20, 30, 40, 50, 60, 70],
                                    [20, 40, 60, 80, 100, 120, 140]],
                      'path_to_catalog': str(path_to_data())}

    return param_dict_val


class CustomError(Exception):
    pass


def path_to_data():
    """ Определяет путь к каталогу текущих данных на конкретной машине и к корню каталога. """
    head_path1 = Path(r'E:\Fast_Acquisition')  # Путь к каталогу данных для домашнего ноута
    head_path2 = Path(r'/media/anatoly/Samsung_T5/Fast_Acquisition')  # Путь к каталогу данных для рабочего компа
    head_path3 = Path(r'C:\SCIENCE\PYTHON 3\Fast_Acquisition')  # Путь к каталогу данных для ноута ВМ
    head_path4 = Path(r'J:\Fast_Acquisition')  # Путь к каталогу данных для notebook 'Khristina'
    head_path5 = Path(r'H:\Fast_Acquisition')  # Путь к каталогу данных для нового ноута
    head_path5a = Path(r'E:\Fast_Acquisition')  # Путь к каталогу данных для нового ноута
    if head_path1.is_dir():
        head_path_out = head_path1
    elif head_path2.is_dir():
        head_path_out = head_path2
    elif head_path3.is_dir():
        head_path_out = head_path3
    elif head_path4.is_dir():
        head_path_out = head_path4
    elif head_path5.is_dir():
        head_path_out = head_path5
    elif head_path5a.is_dir():
        head_path_out = head_path5a
    else:
        raise CustomError('Path to data is not found!')

    path_intermediate = '2021/Results'
    head_path_out = Path(head_path_out, path_intermediate)
    # file_path_data_out = Path(head_path_out, current_catalog_in, current_data_dir_in)
    return head_path_out

