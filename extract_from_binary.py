import numpy as np
import struct
import os
import scipy.io
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from datetime import datetime
import Fig_plot as fp
from path_to_Yandex_Disk import path_to_YaDisk
from afc_alignment import align_func
from afc_alignment1 import align_func1

start = datetime.now()

# head_path = path_to_YaDisk()
head_path = 'E:\\BinData'       # E:\BinData\2020_06_30test
# file_name0 = head_path + '\\Measure\\Fast_Acquisition\\2020_08_04test\\20200804_load1_NS_-06_2'
file_name0 = head_path + '\\2020_08_04test\\20200804_load2_NS_-00_2'
calibration_file_name = 'Noise_afc_-02_-33-3'
q = int(file_name0[-1])

# D:\YandexDisk\Measure\Fast_Acquisition\06022020calibr
# !!!! ******************************************* !!!!
# ****** Блок исходных параметров для обработки *******
kf = 1  # Установка разрешения по частоте
kt = 16  # Установка разрешения по времени
N_Nyq = q  # Номер зоны Найквиста
# *****************************************************

delta_t = 8.1925e-3
delta_f = 7.8125
robust_filter = 'n'
param_robust_filter = 1.1
align = 'n'
noise_calibr = 'n'
graph_3d_perm = 'n'
contour_2d_perm = 'n'

t_start_flame = 103.6
t_stop_flame = 105

if N_Nyq == 3:
    freq_spect_mask = [2250, 2315, 2370, 2500, 2600, 2670, 2690, 2730, 2735, 2740]  # 2060, 2750, 2760, 2770, 2780, 2790, 2800, 2810,
                       # 2820, 2830, 2850, 2880, 2900, 2950# Временные сканы Солнца на этих частотах
else:
    freq_spect_mask = [1050, 1171, 1380, 1465, 1500, 1535, 1600, 1700, 1950, 2000]

# Динамическая маска (зависит от длины записи во времени)
# time_spect_mask = [(lambda i: (t_spect * (i + 0.05) / 7) // 10 * 10)(i) for i in range(7)]
time_spect_mask = [47, 84.4, 104, 133, 133.05, 177.02, 177.38]  # Срез частотного спектра в эти моменты времени
# 173, 173.6, 173.8, 174.38
# t_cal = [0, 14, 17, 31]         # Для скана "20200318-1353_-24-3"
# t_cal = [0, 13, 17, 35]
def extract(file_name0):
    file_name = file_name0 + '.bin'
    file_name_out = file_name0 + '.txt'
    i = 0
    k = 0
    spectr = []
    frame = ' '

    try:
        if os.path.isfile(file_name) == 1:
            pass
        else:
            print('\n \t', file_name, ' not found!!!\n')

        f_in = open(file_name, 'rb')

        while frame:

            spectr_frame = []
            # Обработка кадра: выделение номера кадра, границ куртозиса, длины усреднения на ПЛИС
            # и 128-ми значений спектра в список spectr_frame на позиции [1:128]
            for k in range(129):
                frame = f_in.read(8)
                frame_int = int.from_bytes(frame, byteorder='little')
                if k == 0:
                    frame_num = frame_int & 0xFFFFFFF

                    # Выделение длины усреднения (количество усредняемых на ПЛИС отсчетов спектра = 2^n_aver)
                    # Выделение промежутка для значения куртозиса = [2 - bound_left/64, 2 + bound_right/64])
                    if i == 0:
                        n_aver = (frame_int & 0x3F00000000) >> 32
                        bound_left = (frame_int & 0x7FC000000000) >> (32 + 6)
                        bound_right = (frame_int & 0xFF800000000000) >> (32 + 6 + 9)
                    # Запись на первую позицию (с индексом 0) фрагмента спектра номера кадра frame_num
                    spectr_frame.append(frame_num)

                else:
                    spectr_val = (frame_int & 0x7FFFFFFFFFFFFF)
                    pp_good = (frame_int & 0xFF80000000000000) >> 55
                    spectr_frame.append(spectr_val)
                    pass

            spectr.append(spectr_frame)
            print(i, frame_num)
            i += 1

        pass

        spectr.pop(-1)
        N = len(spectr)
        n_frame_last = spectr[-1][0]
        rest = (n_frame_last + 1) % 2**(6 - n_aver)
        if rest:
            for k in range(rest):
                spectr.pop(-1)
        print(n_frame_last, spectr[-1][0])
    finally:
        f_in.close()
        pass

        spectr1 = convert_to_matrix(spectr, spectr[-1][0] + 1, n_aver)
    np.savetxt(file_name_out, spectr1, header=(str(n_aver) + '-n_aver ' + str(bound_left) + '-kurt'))

    return spectr1, n_aver


def convert_to_matrix(S_total, counter, n_aver):
    '''Функция принимает список списков S, который формируется в extract(file_name0) и превращает его в матрицу,
    строки которой - спектры с разрешением 7.8125/(2**(6-n_aver)) МГц, а столбцы - зависимость значения
    спектра на фиксированной частоте от времени. Разрешение по времени - 8192 мкс. Вместо пропущенных пакетов
    вставляет значение 2'''

    S = [[int(2)] * 128 for i in range(counter)]
    for s in S_total:
        S[s[0]] = s[1:]
    n = 128 * (2 ** (6 - n_aver))
    print(len(S))
    s_ar = np.reshape(S, (-1, n))
    return s_ar


def spectr_construction(Spectr, kf, kt):
    ''' Функция формирует спектр принятого сигнала с требуемым разрешением по частоте и времени. Максимальное
    разрешение отсчетов по времени 8192 мкс и по частоте 7,8125 МГц. Путем суммирования и усреднерия по kt*kf
    отсчетам разрешение по частоте и по времени в исходном спектре Spectr уменьшается в kf и kt раз,
    соответственно. Преобразованный спектр возвращается как S1.
    '''

    N_col1 = N_col // kf
    N_row1 = N_row // kt
    S1 = np.zeros((N_row1, N_col1))

    for i in range(N_row1):
        for j in range(N_col1):
            try:
                S1[i, j] = np.sum(Spectr[i * kt: (i + 1) * kt, j * kf: (j + 1) * kf])
                N_mesh = (Spectr[i * kt: (i + 1) * kt, j * kf: (j + 1) * kf] != 0).sum()
                if N_mesh == 0:
                    S1[i, j] = 2
                else:
                    S1[i, j] = S1[i, j] / N_mesh
                if S1[i, j] == 0:
                    S1[i, j] = 2
                # if (j > 3) & (S1[i, j] > 1.5 * np.sum(S1[i, j-3:j])//3):
                #     S1[i, j] = np.sum(S1[i, j-3:j])//3
                if robust_filter == 'y':
                    a = param_robust_filter
                    if (i > 3) & (S1[i, j] < 1 / a * np.sum(S1[i - 3:i - 1, j]) // 2):
                        S1[i, j] = np.sum(S1[i - 1, j])
                    if (i > 3) & (S1[i, j] > a * np.sum(S1[i - 3:i - 1, j]) // 2):
                        # print(S1[i - 3:i+1, j])
                        S1[i, j] = np.sum(S1[i - 1, j])
                        # print(S1[i, j])
                        pass

            except IndexError as allert_message:
                print(allert_message, 'ind i = ', i, 'ind j = ', j)
                pass
            except ValueError as value_message:
                print(value_message, 'ind i = ', i, 'ind j = ', j)
                pass

    return S1  # // kt // kf


def line_legend(freq_mask):
    N_col_leg = len(freq_mask)
    N_row_leg = len(time_spect_mask)
    legend_freq = [0] * N_col_leg
    legend_time = [0] * N_row_leg
    i1 = 0
    for i in freq_mask:
        legend_freq[i1] = str(i) + ' MHz'
        i1 += 1
    i1 = 0
    for i in time_spect_mask:
        legend_time[i1] = str(i) + ' sec'
        i1 += 1

    return legend_time, legend_freq


def form_spectr_sp():
    """ Возвращает s_freq - срезы частотного спектра в моменты времени time_spect_mask и s_time - сканы Солнца
    по времени на частотах freq_spect_mask

    """
    ind_spec = []
    ind_time = []
    s_freq = np.zeros((len(time_spect_mask), N_col // kf))
    s_time = np.zeros((N_row // kt, len(freq_spect_mask)))
    i = 0
    for f in freq_spect_mask:
        ind = int((f - (N_Nyq - 1) * 1000) // (delta_f / aver_param * kf))
        if ind > N_col // kf - 1:
            ind = -1
        s_time[:, i] = spectr_extr1[:, ind]
        ind_spec.append(ind)
        i += 1
    i = 0
    for t in time_spect_mask:
        ind = int(t // (delta_t * kt))
        if ind > N_row // kt - 1:
            ind = -1
        s_freq[i, :] = spectr_extr1[ind, :]
        ind_time.append(ind)
        i += 1
    s_time = s_time.transpose()
    return s_freq, s_time


def form_spectr_sp1(freq_spect_mask_in=freq_spect_mask, time_spect_mask_in=time_spect_mask):
    """ Возвращает s_freq - срезы частотного спектра в моменты времени time_spect_mask и s_time - сканы Солнца
    по времени на частотах freq_spect_mask с заданным разрешением по времени и частоте

    """
    ind_spec = []
    ind_time = []
    s_freq = np.zeros((len(time_spect_mask_in), N_col // kf))
    s_time = np.zeros((N_row // kt, len(freq_spect_mask_in)))
    j = 0
    for f in freq_spect_mask_in:
        ind1 = (f - (N_Nyq - 1) * 1000 - delta_f / aver_param / 2) // (delta_f / aver_param)
        ind = int(ind1)
        if ind > N_col - int(kf / 2) - 1:
            ind = N_col - int(kf / 2) - 1
        if ind < int(kf / 2):
            ind = int(kf / 2)
        i = 0
        while kt * (i + 1) < N_row:
            if kf == 1:
                s_time[i, j] = np.sum(spectr_extr[i * kt:(i + 1) * kt, ind])
            else:
                s_time[i, j] = np.sum(spectr_extr[i * kt:(i + 1) * kt, ind - int(kf / 2):ind + int(kf / 2)])
            i += 1
        ind_spec.append(ind)
        j += 1
    i = 0
    for t in time_spect_mask_in:
        ind = int(t // delta_t)
        if ind > N_row - kt / 2 - 1:
            ind = N_row - int(kt / 2) - 1
        if ind < (kt / 2):
            ind = int(kt / 2)
        j = 0
        while (j + 1) * kf < N_col:
            if kt == 1:
                s_freq[i, j] = np.sum(spectr_extr[ind, j * kf:(j + 1) * kf])
            else:
                s_freq[i, j] = np.sum(spectr_extr[ind - int(kt / 2):ind + int(kt / 2), j * kf:(j + 1) * kf])
            j += 1
        ind_time.append(ind)
        i += 1
    s_time = s_time.transpose()
    return s_freq // kt // kf, s_time // kt // kf


def pic_title():
    title0 = file_name0[-19:-2]
    title1 = '  ' + title0[0:4] + '.' + title0[4:6] + '.' + title0[6:8] + \
             ' time=' + title0[9:11] + ':' + title0[11:13] + ' azimuth=' + title0[14:17]
    if not file_name0.find('sun') == -1:
        title2 = 'Sun intensity'
    elif not file_name0.find('crab') == -1:
        title2 = 'Crab intensity'
    elif not file_name0.find('calibr') == -1:
        title2 = 'Calibration'
        title0 = file_name0[-23:-2]
        title1 = '  ' + title0[0:4] + '.' + title0[4:6] + '.' + title0[6:8] + \
                 ' chanell att=' + title0[14:17] + ' source att=' + title0[18:21]
    elif not file_name0.find('test') == -1:
        title0 = file_name0[-24:-2]
        title2 = 'Test interference'
        title1 = '  ' + title0[0:4] + '.' + title0[4:6] + '.' + title0[6:8] + \
                 ' chanell att=' + title0[15:18] + ' source att=' + title0[19:22]
        pass
    else:
        title2 = []
    return title1, title2


def graph_3d(*args):
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    xval, yval, z, s = args
    x, y = np.meshgrid(xval, yval)
    ax.zaxis._set_scale('log')  # Расставляет tiks логарифмически
    title1, title2 = pic_title()
    ax.set_title(title2 + ' ' + title1, fontsize=20)
    ax.text2D(0.05, 0.75, info_txt[0], transform=ax.transAxes, fontsize=16)
    ax.text2D(0.05, 0.65, info_txt[1], transform=ax.transAxes, fontsize=16)
    ax.set_xlabel('Frequency, MHz', fontsize=16)
    ax.set_ylabel('Time, s', fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=14)
    # cmap = plt.get_cmap('jet')
    if s:
        surf = ax.plot_surface(x, y, z, rstride=2, cstride=2, cmap=cm.plasma)
        plt.savefig(file_name0 + '_wK' + '.png', format='png', dpi=100)
        return
    surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.jet)
    add_path0 = fp.path_to_pic(file_name0 + '\\', 3)
    plt.savefig(file_name0 + '\\' + add_path0, format='png', dpi=100)
    plt.show()
    return


def graph_contour_2d(*args):
    import matplotlib.font_manager as font_manager
    xval, yval, z, s = args
    x, y = np.meshgrid(xval, yval)
    z = np.log10(z)

    levels = MaxNLocator(nbins=15).tick_values(z.min(), z.max())
    # pick the desired colormap, sensible levels, and define a normalization
    # instance which takes data values and translates those into levels.
    cmap = plt.get_cmap('jet')

    fig, ax1 = plt.subplots(1, figsize=(12, 6))

    cf = ax1.contourf(x, y, z, levels=levels, cmap=cmap)

    x_min = xval[1]
    y1 = yval[0] + (yval[-1] - yval[0]) * 0.05
    y2 = yval[0] + (yval[-1] - yval[0]) * 0.1
    fig.colorbar(cf, ax=ax1)
    title1, title2 = pic_title()
    ax1.set_title(title2 + ' ' + title1, fontsize=20)
    ax1.set_xlabel('Freq, MHz', fontsize=18)
    ax1.set_ylabel('Time, s', fontsize=18)

    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)
    plt.tick_params(axis='both', which='major', labelsize=16)

    plt.text(x_min, y1, info_txt[0], fontsize=16)
    plt.text(x_min, y2, info_txt[1], fontsize=16)

    # adjust spacing between subplots so `ax1` title and `ax0` tick labels
    # don't overlap
    fig.tight_layout()
    add_path0 = fp.path_to_pic(file_name0 + '\\', 2, 'png')
    fig.savefig(file_name0 + '\\' + add_path0)
    plt.show()
    return

    # Модуль проверки: формировалась ли ранее матрица спектра по времени и частоте
    # если - нет, то идем в extract(file_name0), если - да, то загружаем


def calibration(t_cal, spectrum):
    size_spectrum = spectrum.shape
    level_matrix = np.ones((size_spectrum[0], 2))
    cal_level = np.zeros(size_spectrum[0])
    n_cal = np.zeros(4)
    for i in range(4):
        n_cal[i] = int(t_cal[i] // (delta_t * kt))
    for sp_row in range(size_spectrum[0]):
        level_matrix[sp_row, 0] = np.average(spectrum[sp_row, int(n_cal[0]):int(n_cal[1])])
        level_matrix[sp_row, 1] = np.average(spectrum[sp_row, int(n_cal[2]):int(n_cal[3])])
        cal_level[sp_row] = np.abs(level_matrix[sp_row, 0] - level_matrix[sp_row, 1])
    max_cal_level = np.max(cal_level)
    cal_level = cal_level / max_cal_level

    for sp_row in range(size_spectrum[0]):
        spectrum[sp_row, :] = (spectrum[sp_row, :] - np.min(level_matrix[sp_row, :])) / cal_level[sp_row]
    return spectrum


def self_calibration1():
    """ Принимает среднее значение спектра за первые 5 - 10 сек наблюдения, пока не идет Солнце, при максимальном
    разрешении по частоте kf = 1

    :return:
    """
    file_calibr = 'self_calibr.txt'  # Файл, в который записываются спектры наблюдений в отсутствие Солнца
    if kt <= 500 and kf > 1:
        return

    ind_observation_id = file_name0.rfind('\\') + 1  # Индекс, с которого вычленяется идентификатор записи
    # наблюдения из ее полного пути
    observation_id = file_name0[ind_observation_id:]  # Идентификатор записи наблюдения, использованной
    # для учета АЧХ тракта
    file_observation_names = file_name0[0:ind_observation_id - 1] + '\\observ_name.txt'  # Файл, в который
    # записываются идентификаторы записей наблюдений, используемых для учета АЧХ тракта
    file_observ_name = open(file_observation_names, 'a+')
    file_observ_name.seek(0)
    observation_list = file_observ_name.read()
    if not observation_list.count(observation_id):
        # if not os.path.isfile(file_observation_names):
        if not os.path.isfile(file_calibr):
            np.savetxt(file_calibr, spectr_freq[0, :])
        else:
            self_calibr = np.loadtxt(file_calibr)
            try:
                self_calibr1 = np.vstack((self_calibr, spectr_freq[0, :]))
            except ValueError:
                print('Function self_calibration not append new data')
                file_observ_name.close()
                return
            np.savetxt(file_calibr, self_calibr1)
            file_observ_name.write(observation_id + '\n')
    file_observ_name.close()
    return


def path_to_fig():
    """ Создает директорию для рисунков обрабатываемого наблюдения, если она до этого не была создана,
    название директории  совпадает с названием исходного файла данных наблюдения
    """
    if not os.path.isdir(file_name0):
        os.mkdir(file_name0)
    return


if not os.path.isfile(file_name0 + '.txt'):
    spectr_extr, n_aver = extract(file_name0)
else:
    spectr_extr = np.loadtxt(file_name0 + '.txt')
    f_in1 = open(file_name0 + '.txt')
    n_aver = int((f_in1.readline())[2])
    f_in1.close()

aver_param = 2 ** (6 - n_aver)

if align == 'y':
    align_coeff = align_func(calibration_file_name, 'y', aver_param)
    spectr_extr = spectr_extr * align_coeff

print('spectr_extr.shape = ', spectr_extr.shape)

# Приведение номеров отсчетов спектра в соответствие с ростом частоты слева направо
N_row, N_col = np.shape(spectr_extr)
t_spect = N_row * delta_t
if N_Nyq % 2 == 0:
    for i in range(N_row):
        spectr_extr[i][0:] = spectr_extr[i][-1::-1]

# Формирование спектров и сканов по маскам freq_spect_mask и time_spect_mask
# Динамическая маска (зависит от длины записи во времени)
# time_spect_mask = [(lambda i: (t_spect * (i + 0.05) / 7) // 10 * 10)(i) for i in range(7)]
time_spect_mask = [(lambda i: (t_spect * (i + 0.05)) // 7)(i) for i in range(7)]
spectr_freq, spectr_time = form_spectr_sp1(freq_spect_mask, time_spect_mask)
# np.savetxt(file_name0+'_scan'+'.txt', spectr_time)
# np.savetxt(file_name0+'_spectr'+'.txt', spectr_extr1)
if noise_calibr == 'y':
    spectr_time = calibration(t_cal, spectr_time)

# Формирование строк-аргументов по времени и частоте и легенды
freq = np.linspace(1000 * (N_Nyq - 1) + 3.9063 / aver_param * kf, 1000 * N_Nyq - 3.9063 / aver_param * kf, N_col // kf)
timeS = np.linspace(0, delta_t * N_row, N_row // kt)

# ***********************************************
# ***        Графический вывод данных        ****
# ***********************************************

# Укрупнение  разрешения по частоте и времени для вывода в 2d и 3d
if graph_3d_perm == 'n' or contour_2d_perm == 'n':
    spectr_extr1 = spectr_construction(spectr_extr, kf, kt)
# Информация о временном и частотном резрешениях
info_txt = [('time resol = ' + str(delta_t * kt) + 'sec'),
            ('freq resol = ' + str(delta_f / aver_param * kf) + 'MHz')]
path_to_fig()

if graph_3d_perm == 'y':
    graph_3d(freq, timeS, spectr_extr1, 0)
if contour_2d_perm == 'y':
    graph_contour_2d(freq, timeS, spectr_extr1, 0)


# if align == 'y':
#     align_coeff1 = align_func1(spectr_freq[1, :], 'y', aver_param)
#     spectr_extr = spectr_extr * align_coeff1


line_legend_time, line_legend_freq = line_legend(freq_spect_mask[:10])

fp.fig_plot(spectr_freq, 0, freq, 1, info_txt, file_name0, line_legend_time)
fp.fig_plot(spectr_time, 0, timeS, 0, info_txt, file_name0, line_legend_freq)
n_start_flame = int(t_start_flame // (delta_t * kt))
n_stop_flame = int(t_stop_flame // (delta_t * kt))
# if graph_3d_perm == 'y':
#     graph_3d(freq, timeS[n_start_flame:n_stop_flame], spectr_extr1[n_start_flame:n_stop_flame, :], 0)
# fp.fig_multi_axes(spectr_time[:10, n_start_flame:n_stop_flame], timeS[n_start_flame:n_stop_flame],
#                   info_txt, file_name0, freq_spect_mask[:10])
stop = datetime.now()
print('\n Total time = ', stop - start)
