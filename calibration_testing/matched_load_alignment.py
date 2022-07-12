import os
import numpy as np
import sys
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from pathlib import Path
from Supporting_func.stocks_coefficients import path_to_data


if __name__ == '__main__':

    """ Расчет выравнивающих коэффициентов АЧХ приемника по шумовому сигналу от согласованной нагрузки на входе
    с учетом собственных шумов приемника."""

    current_data_dir = '2022'
    primary_data_dir = 'Primary_data'  # Каталог исходных данных (за определенный период, здесь - год)
    converted_data_dir = 'Converted_data'  # Каталог для записи результатов конвертации данных и заголовков
    data_treatment_dir = 'Data_treatment'  # Каталог для записи результатов обработки, рисунков

    current_primary_dir = '2022_06_28test'
    current_converted_dir = current_primary_dir + '_conv'
    current_converted_path = Path(converted_data_dir, current_converted_dir)
    current_treatment_dir = current_primary_dir + '_treat'
    current_treatment_path = Path(data_treatment_dir, current_treatment_dir)

    ngi_temperature_file_name = 'ngi_temperature.npy'  # Файл усредненной по базе шумовой температуры для ГШ
    receiver_temperature_file_name = 'receiver_temperature.npy'
    current_primary_file1 = '2022-06-28_01test'  # Файл с согласованной нагрузкой на обоих входах приемника
    current_primary_file2 = '2022-06-28_02test'  # Файл с согласованной нагрузкой и КЗ на входах приемника
    current_primary_file2 = '2022-06-28_03test'  # Файл с КЗ и согласованной нагрузкой на входах приемника

    converted_data_file_path, head_path = path_to_data(current_data_dir, current_converted_path)
    data_treatment_file_path, head_path = path_to_data(current_data_dir, current_treatment_path)
    ngi_temperature_path = Path(head_path, 'Alignment', ngi_temperature_file_name)
    receiver_temperature_path = Path(head_path, 'Alignment', receiver_temperature_file_name)

    date = current_primary_file1[0:10]

    delta_t = 8.3886e-3
    _delta_f = 7.8125
    aver_param_noise = 8
    temp0 = 300
    time_mask = [11, 19, 21, 29]