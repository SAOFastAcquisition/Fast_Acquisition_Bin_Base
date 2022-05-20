import os
import numpy as np
import sys
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from pathlib import Path
from Supporting_func.stocks_coefficients import path_to_data


def noise_kp(spectrum_noise_out, diff='n'):
    """ Функция принимает файл шумового измерения мощностной АЧХ и возвращает
        нормированный к максимуму коэффициент передачи тракта по мощности"""
    spectrum_noise_out = np.array(spectrum_noise_out, dtype=np.int64)
    n_row, n_col = np.shape(spectrum_noise_out)
    s = np.zeros(n_col, dtype=np.int64)
    s1 = np.zeros(n_col)
    # Усреднение по времени
    for i in range(n_col):
        if diff == 'n':
            s_loc = spectrum_noise_out[:, i]
            s_loc_mod = s_loc[s_loc > 100]
            # s[i] = np.sum(spectrum_noise_out[:500, i]) / 500
            s[i] = np.sum(s_loc_mod[1700:2700]) / 1000
        else:
            s[i] = np.sum(spectrum_noise_out[:1600, i]) / 1600
            s1[i] = np.sum(spectrum_noise_out[2000:3600, i]) / 1600
            s[i] -= s1[i]

    s_max = np.max(s)
    kp_norm = s / s_max
    # plt.plot(s)
    # plt.show()
    # plt.plot(kp_norm)
    # plt.show()
    return kp_norm, n_col, s_max


def align_func(s_loc: object, aver_param_noise: object, diff: object = 'n') -> object:
    """ Функция возвращает коэффициенты, выравнивающие исходную АЧХ

    """
    # Исходные данные
    delta_f = 7.8125

    kp_norm, n_col, s_max = noise_kp(s_loc, diff)

    # Исключение корректировки коэффициента усиления в зоне действия режекторных фильтров
    if n_nyq == 3:
        n1 = int((90 - delta_f / 2 / aver_param_noise) // (delta_f / aver_param_noise))
        n2 = int((220 - delta_f / 2 / aver_param_noise) // (delta_f / aver_param_noise))
        n3 = int((540 - delta_f / 2 / aver_param_noise) // (delta_f / aver_param_noise))
        n4 = int((700 - delta_f / 2 / aver_param_noise) // (delta_f / aver_param_noise))
        kp_norm[n3:n4] = 1
    else:
        n1 = int((80 - delta_f / 2 / aver_param_noise) // (delta_f / aver_param_noise))
        n2 = int((230 - delta_f / 2 / aver_param_noise) // (delta_f / aver_param_noise))
    kp_norm[n1:n2] = 1

    # band = int(aver_param_noise / aver_param)
    # kp_band = [np.sum(kp_norm[band * i:band * (i + 1)]) / band for i in range(int(n_col / band))]
    align_coeff = np.array([1 / a for a in kp_norm])

    if n_nyq == 3:
        freq = np.linspace(1000 * (n_nyq - 1) + 3.9063 / aver_param, 1000 * n_nyq - 3.9063 / aver_param, n_col)
    else:
        freq = np.linspace(1000 * n_nyq - 3.9063 / aver_param, 1000 * (n_nyq - 1) + 3.9063 / aver_param, n_col)

    return align_coeff, s_max


def align_visualization(coeff_set):
    n_col = 1024
    n_aver = 3
    aver_param = 2 ** (6 - n_aver)
    align_lower_left = np.array([calibration_frame.spectrum_left1.iloc[s] for s in coeff_set])
    align_upper_left = np.array([calibration_frame.spectrum_left2.iloc[s] for s in coeff_set])
    align_left = np.concatenate((align_lower_left, align_upper_left), axis=1)
    align_lower_right = np.array([calibration_frame.spectrum_right1.iloc[s] for s in coeff_set])
    align_upper_right = np.array([calibration_frame.spectrum_right2.iloc[s] for s in coeff_set])
    align_right = np.concatenate((align_lower_right, align_upper_right), axis=1)

    freq2 = np.linspace(1000 * (3 - 1) + 3.9063 / aver_param, 1000 * 3 - 3.9063 / aver_param, n_col)
    freq1 = np.linspace(1000 * 2 - 3.9063 / aver_param, 1000 * (2 - 1) + 3.9063 / aver_param, n_col)
    freq = np.concatenate((freq1, freq2))

    # fig, ax = plt.subplots(1, figsize=(12, 6))
    # ax.plot(freq, align_right[0])
    # ax.plot(freq, align_right[1])
    # ax.plot(freq, align_right[2])
    # ax.plot(align_coeff[1])
    # ax.plot(align_coeff[2])
    # ax.plot(align_coeff[3])
    # plt.show()

    pass


def convert_to_txt(path_to_align, _index):
    with open(path_to_align, 'rb') as _inp:
        _calibration_frame = pickle.load(_inp)

    align_left1 = _calibration_frame.iloc[_index].spectrum_left1
    align_left2 = _calibration_frame.iloc[_index].spectrum_left2
    align_right1 = _calibration_frame.iloc[_index].spectrum_right1
    align_right2 = _calibration_frame.iloc[_index].spectrum_right2
    align = [align_left1, align_left2, align_right1, align_right2]
    align_name = ['align_left1', 'align_left2', 'align_right1', 'align_right2']
    i = 0
    for _s in align:
        file_name_out: str = align_name[i] + '.txt'
        # Расчет и запись коэффициентов корректировки АЧХ в файл
        # При этом для второй зоны Найквиста - инверсный порядок по частоте
        path_to_align_txt = str(path_to_align)[0:-5] + '2_txt'
        path_align_txt = Path(path_to_align_txt, file_name_out) # r'H:\Fast_Acquisition\Alignment\Align_coeff2_txt'
        np.savetxt(path_align_txt, _s)
        i += 1



# ******************** Путь к исходным данным *********************
current_data_file = '2022-03-27_02'  # Имя файла с исходными текущими данными без расширения
current_data_dir = '2022_03_27calibr_conv'  # Папка с текущими данными
align_file_name = 'Align_coeff.bin'  # Имя файла с текущими коэффициентами выравнивания АЧХ
if current_data_file[0:4] == '2021':
    current_catalog = r'2021\Results'  # Текущий каталог (за определенный период, здесь - год)
if current_data_file[0:4] == '2022':
    current_catalog = r'2022/Converted_data'

file_path_data, head_path = path_to_data(current_catalog, current_data_dir)

#                   ***********************************
# Путь к файлу хранения коэффициентов (он один для всех калибровок АЧХ каналов)
folder_align_path = Path(head_path, 'Alignment')

visualization_set = [0, 1]

# Если файла хранениия коэффициентов не существует, то создаем его, если существует - загружаем
columns_names = ['date', 'att1', 'att2', 'att3',
                 'spectrum_left1', 'polar', 'spectrum_left2', 'spectrum_right1',
                 'spectrum_right2',
                 'max_left1', 'max_left2', 'max_right1', 'max_right2', 'flag_align']
if not os.path.isfile(Path(folder_align_path, align_file_name)):
    calibration_frame = pd.DataFrame(columns=columns_names)
else:
    with open(Path(folder_align_path, align_file_name), 'rb') as inp:
        calibration_frame = pickle.load(inp)
# ************************************************************************************

# align_visualization(visualization_set)

# Загрузка исходных данных ('_spectrum.npy' - сами данные, _head.bin - вспомогательные)
spectrum = np.load(Path(file_path_data, current_data_file + '_spectrum.npy'), allow_pickle=True)
with open(Path(file_path_data, current_data_file + '_head.bin'), 'rb') as inp:
    head = pickle.load(inp)
n_aver = head['n_aver']
aver_param = 2 ** (6 - n_aver)
if head['polar'] == 'left' or head['polar'] == 'right':
    head['polar'] = 'half'
align_coeff = [np.nan, np.nan, np.nan, np.nan]
s_max_band = np.empty(4)
s_max_band[:] = np.nan
i = 0

# Расчет выравнивающих коэффициентов по отдельности для каждой поляризации и полосы частот 1-2 или 2-3 ГГц
flag_align = 0  # Выравнивания по всему диапазону еще нет
for s in spectrum:
    if i % 2:
        n_nyq = 3
    else:
        n_nyq = 2
    if np.size(s) > 1:
        s_unipol = s  # [s > 100]
        align_coeff[i], s_max_band[i] = align_func(s_unipol, aver_param)
    i += 1

# Выравнивание коэффициентов по диапазону для каждого канала (левой и правой поляризаций)
if s_max_band[0] and s_max_band[1]:
    band_align1 = s_max_band[0] / s_max_band[1]
    align_coeff[1] = align_coeff[1] * band_align1
if s_max_band[2] and s_max_band[3]:
    band_align2 = s_max_band[2] / s_max_band[3]
    align_coeff[3] = align_coeff[3] * band_align2
# Выравнивание коэффициентов по каналам (левой и правой поляризаций)
if not np.isnan(s_max_band[0]) and not np.isnan(s_max_band[2]):
    channels_align = s_max_band[0] / s_max_band[2]
    align_coeff[2] = align_coeff[2] * channels_align
    align_coeff[3] = align_coeff[3] * channels_align
    flag_align = 1  # Произошло выравнивание по всему диапазону и для всех поляризаций

# print(align_coeff[0])
fig, ax = plt.subplots(1, figsize=(12, 6))
ax.plot(align_coeff[0])
ax.plot(align_coeff[1])
ax.plot(align_coeff[2])
ax.plot(align_coeff[3])
plt.show()

# Рассчитанные коэффициенты вместе с исходной информацией записываем в виде словаря  для формирования
# объекта Series и включения в сводную  таблицу корректирующих коэффициентов
calibrate_row = {'date': current_data_file[:10], 'att1': head['att1'], 'att2': head['att2'], 'att3': head['att3'],
                 'polar': head['polar'], 'spectrum_left1': align_coeff[0], 'spectrum_left2': align_coeff[1],
                 'spectrum_right1': align_coeff[2], 'spectrum_right2': align_coeff[3],
                 'max_left1': s_max_band[0], 'max_left2': s_max_band[1],
                 'max_right1': s_max_band[2], 'max_right2': s_max_band[3],
                 'flag_align': flag_align}  # 'flag_align' - признак выравнивания по всему диапазону : 1 - сделано
calibrate_row_ser = pd.Series(calibrate_row)

# Определяем, есть ли в сводной таблице данные с такими же исходными параметрами. Если есть, то будет их
# проверка на то, содержат они все коэффициенты или нет. Если нет, то объект Series будет полностью
# вставлен в таблицу (объект DataFrame)

if not 'polar' in calibration_frame.columns:
    calibration_frame.polar = ''
idx = calibration_frame.loc[(calibration_frame.date == head['date'])
                            & (calibration_frame.att1 == head['att1'])
                            & (calibration_frame.att2 == head['att2'])
                            & (calibration_frame.att3 == head['att3'])
                            ].index  # & (calibration_frame.polar == head['polar'])

if len(idx):
    r = calibration_frame.iloc[idx[0]]
    # Определяем, есть ли пустые поля в выделенной из таблицы строке. Если есть, то эти поля будут
    # заполнены из имеющегося объекта  Series
    ax_bool = r.isnull()
    r[ax_bool] = calibrate_row_ser[ax_bool]
    # Остались ли в r незаполненные поля
    ax_bool = r.isnull()

    if ('True' not in ax_bool) and (r['flag_align']) == 0:
        # Выравнивание коэффициентов по каналам (левой и правой поляризаций)
        channels_align = r['max_left1'] / r['max_right1']
        r['spectrum_right1'] = r['spectrum_right1'] * channels_align
        r['spectrum_right2'] = r['spectrum_right2'] * channels_align
        r['flag_align'] = 1
        # Если пустые поля заполнились ('True' not in ax_bool) и ранее эта строчка не была полностью сформирована
        # ('flag_align': 0), то строка с номером idx[0] будет удалена из объекта DataFrame и вставлена новая,
        # полностью заполненная r
        calibration_frame = calibration_frame.drop(idx).append(r, ignore_index=True)
else:
    calibration_frame = calibration_frame.append(calibrate_row_ser, ignore_index=True)
# calibration_frame = calibration_frame.drop(axis=0, index=[2, 3])
with open(Path(folder_align_path, align_file_name), 'wb') as out:
    pickle.dump(calibration_frame, out)

with open(Path(folder_align_path, align_file_name), 'rb') as inp:
    calibration_frame_inp = pickle.load(inp)

pass
# convert_to_txt(Path(folder_align_path, align_file_name), 0)
