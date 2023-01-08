from pathlib import Path
import sys
import os


current_dir = Path.cwd()
sys.path.insert(0, current_dir)


class CustomError(Exception):
    pass


def path_to_data():
    """
    Определяет путь на конкретной машине к корню каталога данных.
    """
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
        return 'Err'
    return head_path_out


def create_folder(_path):
    if not os.path.isdir(_path):
        os.mkdir(_path)


class DataPaths(object):

    def __init__(self, _data_file, _data_dir, _main_dir):
        if _data_dir.find('test') != -1 or _data_dir.find('calibration') != -1 or _data_dir.find('calibr') != -1:
            _main_dir = '2022/Test_and_calibration'
        self.data_file_name = _data_file
        self.data_file_prime = _data_file + '.bin'
        self.data_dir = _data_dir
        self.main_dir = _main_dir
        self.head_path = path_to_data()
        self.primary_dir_path, self.primary_data_file_path = self.primary_paths()
        self.converted_dir_path, self.converted_data_file_path = self.converted_paths()
        self.treatment_dir_path, self.treatment_data_file_path = self.treat_paths()

    def primary_paths(self):
        _path = Path(self.head_path, self.main_dir, 'Primary_data', self.data_dir)
        create_folder(_path)
        if self.__check_paths():
            _primary_data_path = Path(_path, self.data_file_prime)
        else:
            raise CustomError('Head path not found!')
        return _path, _primary_data_path

    def converted_paths(self):
        _path = Path(self.head_path, self.main_dir, 'Converted_data', str(self.data_dir) + '_conv')
        create_folder(_path)
        if self.__check_paths():
            _convert_data_path = Path(_path, self.data_file_name)
        else:
            raise CustomError('Head path not found!')
        return _path, _convert_data_path

    def treat_paths(self):
        _path = Path(self.head_path, self.main_dir, 'Data_treatment', str(self.data_dir) + '_treat')
        create_folder(_path)
        if self.__check_paths():
            _treatment_data_path = Path(_path, self.data_file_name)
        else:
            raise CustomError('Head path not found!')
        return _path, _treatment_data_path

    def __check_paths(self):
        return not self.head_path == 'Err'


if __name__ == '__main__':
    data_file_name = '2022-12-24_01+08bb'
    data_dir = '2022_12_24sun'
    main_dir = '2022'
    date = data_dir[0:10]
    adr1 = DataPaths(data_file_name, data_dir, main_dir)

    ngi_temperature_file_name = 'ngi_temperature.npy'  # Файл усредненной по базе шумовой температуры для ГШ
    receiver_temperature_file = 'receiver_temperature.npy'  #
    ant_coeff_file = 'ant_afc.txt'

    ngi_temperature_path = Path(adr1.head_path, 'Alignment', ngi_temperature_file_name)
    receiver_temperature_path = Path(adr1.head_path, 'Alignment', receiver_temperature_file)
    ant_coeff_path = Path(adr1.head_path, 'Alignment', ant_coeff_file)
    pass

