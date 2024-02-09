import numpy as np
import os
import sys
import pickle
from pathlib import Path
from Supporting_func.stokes_coefficients import path_to_data

"""Скрипт обращается к записанному файлу в формате *.bin и записывает его в формате *.txt, создавая при этом для 
него папку, если таковой не существовало. Исходный файл является массивом *.bin и должен формально содержать все 
4 набора данных (скрипт создается для конвертации в текстовый формат коэффициентов, выравнивающих ФЧХ аналогового 
тракта радиометра) - левая поляризация в двух частотных полосах и правая так же. В текстовом формате коэффициенты
записываются поотдельности"""

current_dir = Path.cwd()
sys.path.insert(0, Path(current_dir, 'Supporting_func'))

if __name__ == '__main__':

    # current_data_file = '2021-06-27_19-28'      # Имя файла с исходными текущими данными без расширения
    # current_data_dir = '2021_06_27sun'          # Папка с текущими данными
    # current_catalog = r'2021/Results'           # Текущий каталог (за определенный период, здесь - год)
    current_data_file = r'Align_coeff'  # Имя файла с исходными текущими данными без расширения
    current_data_dir = r''  # Папка с текущими данными
    current_catalog = r'Alignment'  # Текущий каталог (за определенный период, здесь - год)
    pos = 1

    # _date = current_data_file[0:10]
    # year = current_data_file[0:4]

    folder_data_path, head_path = path_to_data(current_catalog, current_data_dir)
    file_data_path = Path(folder_data_path, current_data_file + '.bin')
    folder_data_txt_path = Path(head_path, current_catalog, current_data_file + str(pos) + '_txt')
    file_suffixes = ['_left1.txt', '_left2.txt', '_right1.txt', '_right2.txt']

    #               ******** Считывание *********
    with open(file_data_path, 'rb') as inp:
        calibration_frame_inp = pickle.load(inp)
    r = calibration_frame_inp.iloc[pos]
    align_coeff = [r['spectrum_left1'], r['spectrum_left2'], r['spectrum_right1'], r['spectrum_right2']]
    #               ********** Запись ***********
    if not os.path.isdir(folder_data_txt_path):
        os.makedirs(folder_data_txt_path)
    for i in range(4):
        np.savetxt(Path(folder_data_txt_path, current_data_file + file_suffixes[i]), align_coeff[i])
