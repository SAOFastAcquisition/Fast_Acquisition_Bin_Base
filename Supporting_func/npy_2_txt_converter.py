import numpy as np
import os
import sys
from pathlib import Path
from Supporting_func.stocks_coefficients import path_to_data

"""Скрипт обращается к записанному файлу в формате *.npy и записывает его в формате *.txt, создавая при этом для 
него папку, если таковой не существовало. Исходный файл является массивом *.npy и должен формально содержать все 
4 спектра - левая поляризация в двух частотных полосах и правая так же. В текстовом формате спектры записываются 
поотдельности"""

current_dir = Path.cwd()
sys.path.insert(0, Path(current_dir, 'Supporting_func'))

if __name__ == '__main__':

    current_data_file = '2021-06-27_19-28'      # Имя файла с исходными текущими данными без расширения
    current_data_dir = '2021_06_27sun'          # Папка с текущими данными
    current_catalog = r'2021/Results'           # Текущий каталог (за определенный период, здесь - год)
    date = current_data_file[0:10]
    year = current_data_file[0:4]

    folder_data_path, head_path = path_to_data(current_catalog, current_data_dir)
    file_data_path = Path(folder_data_path, current_data_file + '_spectrum.npy')
    folder_data_txt_path = Path(head_path, year, date + 'sun_txt')
    file_suffixes = ['_left1.txt', '_left2.txt', '_right1.txt', '_right2.txt']

    #               ******** Считывание *********
    spectrum = np.load(file_data_path, allow_pickle=True)
    #               ********** Запись ***********
    if not os.path.isdir(folder_data_txt_path):
        os.makedirs(folder_data_txt_path)
    for i in range(4):
        np.savetxt(Path(folder_data_txt_path, current_data_file + file_suffixes[i]), spectrum[i])



