import os
import numpy as np
import sys
import pandas as pd
import pickle
import matplotlib.pyplot as plt


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
            s[i] = np.sum(s_loc_mod[:500]) / 500
        else:
            s[i] = np.sum(spectrum_noise_out[:1600, i]) / 1600
            s1[i] = np.sum(spectrum_noise_out[2000:3600, i]) / 1600
            s[i] -= s1[i]

    s_max = np.max(s)
    kp_norm = s / s_max
    # plt.plot(s)
    #     # plt.show()
    #     # plt.plot(kp_norm)
    #     # plt.show()
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
    # , freq * 1000000, spectrum_cal


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

    fig, ax = plt.subplots(1, figsize=(12, 6))
    ax.plot(freq, align_right[0])
    ax.plot(freq, align_right[1])
    ax.plot(freq, align_right[2])
    # ax.plot(align_coeff[1])
    # ax.plot(align_coeff[2])
    # ax.plot(align_coeff[3])
    plt.show()

    pass


# Путь к исходным данным
folder_path = r'F:\Fast_Acquisition\2021\Results\2021_04_15test'
file_name = r'\2021-04-15_14'
file_path = folder_path + file_name

#                   ***********************************
# Путь к файлу хранения коэффициентов (он один для всех калибровок АЧХ каналов)
folder_align_path = r'F:\Fast_Acquisition\Alignment'
align_file_name = r'\Align_coeff.bin'

visualization_set = [0, 1, 2]

# Если файла хранениия коэффициентов не существует, то создаем его, если существует - загружаем
columns_names = ['date', 'att1', 'att2', 'att3',
                 'spectrum_left1', 'polar', 'spectrum_left2', 'spectrum_right1',
                 'spectrum_right2',
                 'max_left1', 'max_left2', 'max_right1', 'max_right2', 'flag_align']
if not os.path.isfile(folder_align_path + align_file_name):
    calibration_frame = pd.DataFrame(columns=columns_names)
else:
    with open(folder_align_path + align_file_name, 'rb') as inp:
        calibration_frame = pickle.load(inp)
# ************************************************************************************

align_visualization(visualization_set)

# Загрузка исходных данных ('_spectrum.npy' - сами данные, _head.bin - вспомогательные)
spectrum = np.load(file_path + '_spectrum.npy', allow_pickle=True)
with open(file_path + '_head.bin', 'rb') as inp:
    head = pickle.load(inp)
n_aver = head['n_aver']
aver_param = 2 ** (6 - n_aver)
if head['polar'] == 'left' or head['polar'] == 'right':
    head['polar'] = 'half'
align_coeff = [np.nan, np.nan, np.nan, np.nan]
s_max_band = np.empty(4)
s_max_band[:] = np.nan
i = 0

# Расчет выравнивающих коэффициентов поотдельности для каждой поляризации и полосы частот 1-2 или 2-3 ГГц
flag_align = 0  # Выравнивания по всему диапазону еще нет
for s in spectrum:
    if i % 2:
        n_nyq = 2
    else:
        n_nyq = 3
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
calibrate_row = {'date': file_name[1:11], 'att1': head['att1'], 'att2': head['att2'], 'att3': head['att3'],
                 'polar': head['polar'], 'spectrum_left1': align_coeff[0], 'spectrum_left2': align_coeff[1],
                 'spectrum_right1': align_coeff[2], 'spectrum_right2': align_coeff[3],
                 'max_left1': s_max_band[0], 'max_left2': s_max_band[1],
                 'max_right1': s_max_band[2], 'max_right2': s_max_band[3],
                 'flag_align': flag_align}  # 'flag_align' - признак выравнивания по всему диапазону : 1 - сделано
calibrate_row_ser = pd.Series(calibrate_row)

# Определяем, есть ли в сводной таблице данные с такими же исходными параметрами. Если есть, то будет их
# проверка на то, содержат они все коэффициенты или нет. Если нет, то объект Series будет полностью
# вставлен в таблицу (объект DataFrame)
idx = calibration_frame.loc[(calibration_frame.date == head['date'])
                            & (calibration_frame.att1 == head['att1'])
                            & (calibration_frame.att2 == head['att2'])
                            & (calibration_frame.att3 == head['att3'])
                            & calibration_frame.polar == head['polar']].index

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

with open(folder_align_path + align_file_name, 'wb') as out:
    pickle.dump(calibration_frame, out)

with open(folder_align_path + align_file_name, 'rb') as inp:
    calibration_frame_inp = pickle.load(inp)

pass
