import numpy as np
from numpy.core._multiarray_umath import ndarray

from archive.path_to_Yandex_Disk import path_to_YaDisk
import Fig_plot as fp
from matplotlib.collections import PolyCollection
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from matplotlib import colors as mcolors
import os

# *******************!! Индексы выводимых на графики фантомов !!***********************
phan_map = [0, 1, 2, 3, 4, 5]         #
# ********************!! Индексы выводимых на графики частот !!************************
freq_range = 2          # '1' - нижние 10 частот, '2' - верхние 10 частот
# ********************=====================================************************
flare = 1               # '1' - Отображать неоднородность и ближайшую окрестность

def poly_graph(xs, y, zs):
    fig = plt.figure()
    ax = fig.gca(projection='3d')


    def cc(arg):
        return mcolors.to_rgba(arg, alpha=0.6)

    k, l = np.shape(y)
    verts = []
    for i in range(k):
        ys = y[i, :]
        verts.append(list(zip(xs, ys)))

    poly = PolyCollection(verts, facecolors=['r', 'g', 'c', 'y'])
    poly.set_alpha(0.7)
    ax.add_collection3d(poly, zs=zs, zdir='y')

    ax.set_xlabel('X')
    ax.set_xlim3d(132, 134)
    ax.set_ylabel('Y')
    ax.set_ylim3d(2000, 3000)
    ax.set_zlabel('Z')
    ax.set_zlim3d(0, 2)

    plt.show()


def lin_graph_3d(x, y, z):
    n = np.size(x)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    i = 0
    for f in y:
        ft = [f] * n
        ax.plot(x, ft, z[i, :] - 1, label='parametric curve')
        i += 1
    ax.set_xlabel('X')
    ax.set_xlim3d(132, 134)
    ax.set_ylabel('Y')
    ax.set_ylim3d(2750, 2900)
    ax.set_zlabel('Z')
    ax.set_zlim3d(0, 1)
    plt.show()


def fig_plot1(spectr1, flare, argument, flag, inform, file_name0, line_legend=[2] * 20):

    size_sp1 = spectr1.shape
    freq_line_sp1 = size_sp1[0]
    phantom_num_loc = size_sp1[1]

    title0 = file_name0[-19:-2]

    title1 = '  '+title0[0:4]+'.'+title0[4:6]+'.'+title0[6:8] +\
             ' time='+title0[9:11]+':'+title0[11:13]+' azimuth='+title0[14:17]

    fig, axes = plt.subplots(2, 3, figsize=(12, 6))

    line_colore = ['green', 'blue', 'purple', 'lime', 'black', 'red', 'olivedrab', 'lawngreen', 'magenta', 'dodgerblue']
    if not file_name0.find('sun') == -1:
        title2 = 'Sun intensity'
    elif not file_name0.find('crab') == -1:
        title2 = 'Crab intensity'
    elif not file_name0.find('calibr') == -1:
        title2 = 'Calibration'
        title0 = file_name0[-23:-2]
        title1 = '  ' + title0[0:4] + '.' + title0[4:6] + '.' + title0[6:8] + \
                 ' channel att=' + title0[14:17] + ' source att=' + title0[18:21]
        pass
    else:
        title2 = [ ]
    fig.suptitle(title2+' scan'+title1, y=1.0, fontsize=24)

    i_phantom = 0
    for num_phantom in phan_map:
        m_leg = 0
        for i_freq in range(freq_line_sp1):
            axes[i_phantom // 3, i_phantom % 3].\
                plot(argument[i_phantom, :int(n_finish[num_phantom] - n_start[num_phantom])],
                spectr1[i_freq, i_phantom, :int(n_finish[num_phantom] - n_start[num_phantom])] + 0.2 * m_leg,
                            color=line_colore[m_leg], label=line_legend[i_freq])

            axes[i_phantom // 3, i_phantom % 3].set_title(num_phantom)
            # Show the major grid lines with dark grey lines
            axes[i_phantom // 3, i_phantom % 3].grid(b=True, which='major', color='#666666', linestyle='-')
            # Show the minor grid lines with very faint and almost transparent grey lines
            axes[i_phantom // 3, i_phantom % 3].minorticks_on()
            axes[i_phantom // 3, i_phantom % 3].grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)
            # axes[i_phantom // 3, i_phantom % 3].set_yticks([0, 1])
            # axes[i_phantom // 3, i_phantom % 3].set_ylim(0, 4)
            # axes[i_phantom // 3, i_phantom % 3].set_xlabel('t, sec', fontsize=10)
            axes[i_phantom // 3, i_phantom % 3].xaxis.set_label_coords(1.05, -0.025)
            if flare == 1:
                axes[i_phantom // 3, i_phantom % 3].set_xlim(t_start_burn[num_phantom], t_finish_burn[num_phantom])
            m_leg += 1
            xticks = axes[i_phantom // 3, i_phantom % 3].get_xticks().tolist()
            xticks[-1] = ''
            axes[i_phantom // 3, i_phantom % 3].set_xticklabels(xticks)
            axes[i_phantom // 3, i_phantom % 3].annotate('t, sec', xy=(0.95, -0.04), ha='left', va='top',
                                                         xycoords='axes fraction', fontsize=10)

        i_phantom += 1
    axes[0, 2].legend(loc=10, bbox_to_anchor=(1.2, 0.5))

    plt.show()

    # Управление шрифтом легенды
    font = font_manager.FontProperties(family='Comic Sans MS',
                                       weight='bold',
                                       style='normal', size=16)
    add_pass1 = path_to_pic(file_name0 + '\\', 0)
    fig.savefig(file_name0 + '\\' + add_pass1)


def path_to_pic(file_path,  flag, format='png'):
    if flag == 1:
        add_pass0 = 'spectrum_00'
    elif flag == 2:
        add_pass0 = 'colour2D_00'
    elif flag == 3:
        add_pass0 = 'pic3D_00'
    else:
        add_pass0 = 'scan_00'

    l = len(add_pass0)
    add_pass1 = add_pass0 + '.' + format
    if not os.path.isfile(file_path + add_pass1):
        pass
    else:
        while os.path.isfile(file_path + add_pass1):
            num = int(add_pass0[l-2:l]) + 1
            num_str = str(num)
            if num >= 10:
                add_pass0 = add_pass0[:l - 2] + num_str
            else:
                add_pass0 = add_pass0[:l-2] + '0' + num_str
            add_pass1 = add_pass0 + '.' + format

    return add_pass1


def line_legend():
    N_col_leg = len(freq_spect_mask)

    legend_freq = [0] * N_col_leg

    i1 = 0
    for i in freq_spect_mask:
        legend_freq[i1] = str(i) + ' MHz'
        i1 += 1

    return legend_freq


# Начало и конец изучаемого отрезка записи с фантомом
t_start = np.array([46, 58, 83, 99, 102.8, 116, 128, 174, 174, 200, 200, 219, 225])
t_finish = np.array([50, 63, 87, 102, 106.8, 121, 138, 190, 190, 214, 214, 228, 232])
# Локализация фантома
t_start_burn = np.array([46.8, 58.8, 84, 99.8, 103.6, 117.2, 132.8, 176, 181, 202.5, 207, 221, 228])
t_finish_burn = np.array([47.6, 60, 85, 100.8, 104.6, 118.5, 133.6, 178.4, 185, 206.5, 210, 225, 230])
# Частотная маска
freq_spect_mask0 = [2060, 2250, 2315, 2370, 2440, 2500, 2525, 2580, 2695, 2750, 2760, 2770, 2780,
                   2790, 2800, 2810, 2820, 2830, 2850, 2880, 2915, 2950]    # Маска для скана в азимуте +2
freq_spect_mask01 = [2060, 2250, 2315, 2370, 2500,  2580, 2670, 2720, 2750, 2760, 2770,  2780,
                2790, 2800, 2810, 2820, 2830, 2850, 2880, 2915, 2950]       # Маска для скана в азимуте +6
if freq_range == 1:
    freq_spect_mask = freq_spect_mask0[2:12]
else:
    freq_spect_mask = freq_spect_mask0[12:]

delta_t = 8.1925e-3
kt = 1
delta_t_real = delta_t * kt
delta_f = 7.8125 / 4
kf = 1
aver_param = 1

poly_gr = 'n'
lin_gr = 'n'
multi_gr = 'y'

head_path = path_to_YaDisk()
file_name0 = head_path + '\\Measure\\Fast_Acquisition\\2020_03_19sun\\20200319-1209_+02-3'
scan_init = np.loadtxt(file_name0+'_scan'+'.txt')
l, m = np.shape(scan_init)
flares_count = len(t_start)     # Число фантомов (вспышек)

n_start: ndarray = np.zeros(flares_count)   # Массив номеров начальных отсчетов отрезков записи
n_finish = np.zeros(flares_count)           # Массив номеров конечных отсчетов отрезков записи
n_start_burn = np.zeros(flares_count)       # Массив номеров начальных отсчетов фантомов
n_finish_burn = np.zeros(flares_count)      # Массив номеров конечных отсчетов фантомов
for i in range(flares_count):
    n_start[i] = int(t_start[i] // delta_t_real)
    n_finish[i] = int(t_finish[i] // delta_t_real)
    n_start_burn[i] = int(t_start_burn[i] // delta_t_real)
    n_finish_burn[i] = int(t_finish_burn[i] // delta_t_real)

# Вычисление невозмущенного уровня отрезка записи
background = np.zeros((l, flares_count))
ground_samples = int(max(n_start_burn - n_start + n_finish - n_finish_burn))
scan_ground = np.zeros((l, flares_count, int(max(n_start_burn - n_start + n_finish - n_finish_burn))))
# Выделение значений из отрезков записи для каждого фантома и соответствующих отсчетов времени
scan_init_burn = np.zeros((l, flares_count, int(max(n_finish - n_start))))
# Выделение значений из отрезков записи для каждого фантома и соответствующих отсчетов времени
timeS = np.zeros((flares_count, int(max(n_finish - n_start))))
# Объединение левого и правого невозмущенных крыльев записи для вычисления невозмущенного уровня отрезка
for i in range(flares_count):
    len_gr = int(n_finish[i]-n_start[i]+n_start_burn[i]-n_finish_burn[i])
    scan_ground[:, i, 0:len_gr] = \
        np.hstack((scan_init[:, int(n_start[i]):int(n_start_burn[i])],
                   scan_init[:, int(n_finish_burn[i]):int(n_finish[i])]))
background = np.sum(scan_ground, axis=2) // (scan_ground != 0).sum(2)
# Нормировка отрезка записи к невозмущенному уровню и формирование массива таких записей с осями: частота анализа,
# номер записи/фантома, время
i2 = 0
for i in range(l):
    for i1 in range(flares_count):
        scan_init_burn[i, i1, 0:int(n_finish[i1] - n_start[i1])] = scan_init[i, int(n_start[i1]):int(n_finish[i1])]\
                                                                    / background[i, i1]
    i2 += 1
for i in range(flares_count):
    timeS[i, :int(n_finish[i] - n_start[i])] = np.linspace(t_start[i], t_finish[i], int(n_finish[i] - n_start[i]))
# Массив значений для выводимых на графики фантомов
scan_init_phantom = np.zeros((l, len(phan_map), int(max(n_finish - n_start))))
# Массив значений времени для выводимых на графики фантомов
timeS_phantom = np.zeros((len(phan_map), int(max(n_finish - n_start))))
i2 = 0
for s in phan_map:          # Формирование этих массивов
    scan_init_phantom[:, i2, :] = scan_init_burn[:, int(s), :]
    timeS_phantom[i2, :] = timeS[int(s), :]
    i2 += 1

info_txt = [('time resol = ' + str(delta_t * kt) + 'sec'),
            ('freq resol = ' + str(delta_f / aver_param * kf) + 'MHz')]
line_legend_freq = line_legend()
# fig_plot1(scan_init_flare[1:11, :, :], 1, timeS_flare, 0, info_txt, file_name0, line_legend_freq[1:11])
if freq_range == 1:
    fig_plot1(scan_init_phantom[2:12, :, :], flare, timeS_phantom, 0, info_txt, file_name0, line_legend_freq)
else:
    fig_plot1(scan_init_phantom[12:, :, :], flare, timeS_phantom, 0, info_txt, file_name0, line_legend_freq)
t_spec_start = 132
t_spec_stop = 134
n_spec_start = int(t_spec_start // delta_t_real)
n_spec_stop = int(t_spec_stop // delta_t_real)
timeS1 = np.linspace(t_spec_start, t_spec_stop, n_spec_stop - n_spec_start)
scan_init_burn1 = np.zeros((l, n_spec_stop - n_spec_start))
for i in range(l):
    scan_init_burn1[i, :] = scan_init[i, n_spec_start:n_spec_stop]/background[i]


if multi_gr == 'y':
    scan_init_burn2 = np.zeros((l, n_spec_stop - n_spec_start))
    k = 0
    for i in range(4, 14):
        scan_init_burn2[i, :] = scan_init[i, n_spec_start:n_spec_stop] / background[i] + 0.2 * k
        k += 1
    fp.fig_plot(scan_init_burn2[4:14, :], flare, timeS1, 0, info_txt, file_name0, line_legend_freq[4:14])
pass
if poly_gr == 'y':
    poly_graph(timeS, scan_init_burn, freq_spect_mask)
if lin_gr == 'y':
    lin_graph_3d(timeS1, freq_spect_mask, scan_init_burn1)
