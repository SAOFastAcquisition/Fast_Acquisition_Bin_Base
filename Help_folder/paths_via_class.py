from pathlib import Path
import sys


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


class DataPaths(object):

    def __init__(self, _data_file, _data_dir, _main_dir):
        self.data_file_prime = _data_file + '.bin'
        self.data_dir = _data_dir
        self.main_dir = _main_dir
        self.head_path = path_to_data()
        self.primary_data_path = self.primary_paths()

    def primary_paths(self):
        if self.__check_paths():
            _primary_data_path = Path(self.head_path, self.main_dir, 'Primary', self.data_dir, self.data_file_prime)
        else:
            raise CustomError('Head path not found!')
        return _primary_data_path

    def convert_catalog(self):

        pass

    def treat_catalog(self):

        pass

    def __check_paths(self):
        return not self.head_path == 'Err'


if __name__ == '__main__':
    data_file_name = 'a'
    data_dir = '2022_12_24sun'
    main_dir = '2022'
    adr1 = DataPaths(data_file_name, data_dir, main_dir)
    print(adr1.head_path)
    print(adr1.primary_data_path)

# current_primary_dir = '2022_12_22sun'
# current_primary_file = '2022-12-22_01+08bb'
#
# current_data_dir = '2022'
# # Переопределение каталога всех данных при калибровочных и тестовых наблюдениях
# if current_primary_dir.find('test') != -1 or current_primary_dir.find('calibration') != -1 \
#         or current_primary_dir.find('calibr') != -1:
#     current_data_dir = '2022/Test_and_calibration'
# primary_data_dir = 'Primary_data'  # Каталог исходных данных (за определенный период, здесь - год)
# converted_data_dir = 'Converted_data'  # Каталог для записи результатов конвертации данных и заголовков
# data_treatment_dir = 'Data_treatment'  # Каталог для записи результатов обработки, рисунков
#
# current_converted_dir = current_primary_dir + '_conv'
# current_converted_path = Path(converted_data_dir, current_converted_dir)
# current_treatment_dir = current_primary_dir + '_treat'
# current_treatment_path = Path(data_treatment_dir, current_treatment_dir)
#
# ngi_temperature_file_name = 'ngi_temperature.npy'  # Файл усредненной по базе шумовой температуры для ГШ
# receiver_temperature_file = 'receiver_temperature.npy'  #
# ant_coeff_file = 'ant_afc.txt'
#
# converted_data_file_path, head_path = path_to_data(current_data_dir, current_converted_path)
# data_treatment_file_path, head_path = path_to_data(current_data_dir, current_treatment_path)
# ngi_temperature_path = Path(head_path, 'Alignment', ngi_temperature_file_name)
# receiver_temperature_path = Path(head_path, 'Alignment', receiver_temperature_file)
# folder_align_path = Path(head_path, 'Alignment')
# ant_coeff_path = Path(folder_align_path, ant_coeff_file)
# date = current_primary_file[0:10]
