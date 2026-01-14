import gzip
import os, sys, shutil
import numpy as np
import glob as gb
from datetime import datetime
from pathlib import Path
# from Supporting_func import path_to_data
from Help_folder.paths_via_class import DataPaths

# class DataPaths(object):
#
#     def __init__(self, _data_file, _data_dir, _main_dir):
#         if _data_dir.find('test') != -1 or _data_dir.find('calibration') != -1 or _data_dir.find('calibr') != -1:
#             _main_dir = f'{_main_dir}'
#         self.data_file_name = _data_file
#         self.data_file_prime = _data_file + '.bin'
#         self.data_dir = _data_dir
#         self.main_dir = _main_dir
#         self.head_path = path_to_data()
#         self.primary_dir_path, self.primary_data_file_path = self.primary_paths()
#         self.converted_dir_path, self._converted_data_file_path = self.converted_paths()
#         self.treatment_dir_path, self.treatment_data_file_path = self.treat_paths()
#
#     def primary_paths(self):
#         _path = Path(self.head_path, self.main_dir, 'Primary_data_3_18', self.data_dir)
#         create_folder(_path)
#         if self.__check_paths():
#             _primary_data_path = Path(_path, self.data_file_prime)
#         else:
#             raise CustomError('Head path not found!')
#         return _path, _primary_data_path
#
#     def converted_paths(self):
#         _path = Path(self.head_path, self.main_dir, 'Converted_data_3_18', str(self.data_dir) + '_conv')
#         create_folder(_path)
#         if self.__check_paths():
#             _convert_data_path = Path(_path, self.data_file_name)
#         else:
#             raise CustomError('Head path not found!')
#         return _path, _convert_data_path
#
#     def treat_paths(self):
#         _path = Path(self.head_path, self.main_dir, 'Data_treatment_3_18', str(self.data_dir) + '_treat')
#         create_folder(_path)
#         if self.__check_paths():
#             _treatment_data_path = Path(_path, self.data_file_name)
#         else:
#             raise CustomError('Head path not found!')
#         return _path, _treatment_data_path
#
#     def __check_paths(self):
#         return not self.head_path == 'Err'


class CustomError(Exception):
    pass


def create_folder(_path):
    if not os.path.isdir(_path):
        os.mkdir(_path)


# def path_to_data():
#     """
#     Определяет путь на конкретной машине к корню каталога данных.
#     """
#     head_path1 = Path(r'G:\Fast_Acquisition')           # Путь к каталогу данных для домашнего ноута
#     head_path1a = Path(r'E:\Fast_Acquisition')          # Путь к каталогу данных для домашнего ноута
#     head_path1b = Path(r'G:\Fast_Acquisition_3_18')     # Путь к каталогу данных для домашнего ноута
#     head_path2 = Path(r'/media/anatoly/Samsung_T5/Fast_Acquisition')    # Путь к каталогу данных для рабочего компа
#     head_path2a = Path(r'/media/anatoly/T7/Fast_Acquisition')           # Путь к каталогу данных для рабочего компа
#     head_path3 = Path(r'D:\Fast_acquisition')  # Путь к каталогу данных для ноута ВМ
#     head_path4 = Path(r'J:\Fast_Acquisition')  # Путь к каталогу данных для notebook 'Khristina'
#
#     if head_path1.is_dir():
#         head_path_out = head_path1
#     elif head_path1a.is_dir():
#         head_path_out = head_path1a
#     elif head_path1b.is_dir():
#         head_path_out = head_path1b
#     elif head_path2.is_dir():
#         head_path_out = head_path2
#     elif head_path2a.is_dir():
#         head_path_out = head_path2a
#     elif head_path3.is_dir():
#         head_path_out = head_path3
#     elif head_path4.is_dir():
#         head_path_out = head_path4
#     else:
#         return 'Err'
#     return head_path_out


if __name__ == '__main__':
    start = datetime.now()

    date = '2025-12-17'
    spect = 'sun'        # Вид наблюдений: 'sun', 'test', 'calibr'
    main_dir = date[0:4]
    data_dir = f'{date[0:4]}_{date[5:7]}_{date[8:]}{spect}'
    path_obj = DataPaths(date, data_dir, main_dir)

    az_dict = {'+24': '_01', '+20': '_02', '+16': '_03', '+12': '_04', '+08': '_05', '+04': '_06', '+00': '_07',
               '-04': '_08', '-08': '_09', '-12': '_10', '-16': '_11', '-20': '_12', '-24': '_13'}

    paths = gb.glob(str(Path(path_obj.converted_dir_path, "*.npy")))  # Set Pattern in glob() function
    paths1 = gb.glob(str(Path(path_obj.primary_dir_path, "*.bin")))
    if paths:
        for s in paths:
            filename_in = s
            filename_out = f'{s}.gz'
            if not os.path.exists(filename_out):
                with open(Path(s), 'rb') as inp:
                    current_file = np.load(inp, allow_pickle=True)

                # обратите внимание как открывается выходной файл `gzip.open()`
                with open(filename_in, "rb") as fin, gzip.open(filename_out, "wb") as f_out:
                    # Читает файл по частям, экономя оперативную память
                    shutil.copyfileobj(fin, f_out)

                print(f"Несжатый размер: {str(s)} - {os.stat(filename_in).st_size}")
                print(f"Сжатый размер: {os.stat(filename_out).st_size}")

            if os.path.exists(s) & os.path.exists(f'{s}.gz'):
                os.remove(s)
                print(f"{s} удален.")
            # with gzip.open(filename_out, "rb") as fin:
            #     data = np.load(fin, allow_pickle=True)
            #     # data = fin.read()
            #     print(f"Несжатый размер: {sys.getsizeof(data)}")

    if paths1:
        for s in paths1:
            filename_in = s
            filename_out = f'{s}.gz'
            if not os.path.exists(filename_out):
                # обратите внимание как открывается выходной файл `gzip.open()`
                with open(filename_in, "rb") as fin, gzip.open(filename_out, "wb") as f_out:
                    # Читает файл по частям, экономя оперативную память
                    shutil.copyfileobj(fin, f_out)

                print(f"Несжатый размер: {str(s)} - {os.stat(filename_in).st_size}")
                print(f"Сжатый размер: {os.stat(filename_out).st_size}")

            if os.path.exists(s) & os.path.exists(f'{s}.gz'):
                os.remove(s)
                print(f"{s} удален.")

    stop = datetime.now()
    print(f'Process duration = {stop - start} hour:minute:sec')
