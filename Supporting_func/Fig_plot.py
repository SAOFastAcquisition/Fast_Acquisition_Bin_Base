import matplotlib.pyplot as plt
import numpy as np
import matplotlib.font_manager as font_manager
from matplotlib.ticker import MaxNLocator, ScalarFormatter, FixedLocator
import pylab
import os
from pathlib import Path

# from IPython.display import set_matplotlib_formats
# set_matplotlib_formats('svg')

''' На входе программы несколько пар строк с данными и координатами, причем количество точек, в общем случае, разное для 
    разных пар. В программу передается легенда и название осей. 
    Программа создает сетку, надписывает оси, дает легенду.
'''


def fig_plot(spectr1, burn, argument, flag, inform, file_name0_path, head, line_legend=[2] * 20):
    file_name0 = str(file_name0_path)
    size_sp1 = spectr1.shape

    for i in range(size_sp1[0]):
        for j in range(size_sp1[1]):
            if spectr1[i, j] < 2:
                spectr1[i, j] = 'NaN'

    freq_line_sp1 = size_sp1[0]

    fig, ax = plt.subplots(1, figsize=(12, 6))

    line_color = ['green', 'blue', 'purple', 'lime', 'black', 'red', 'olivedrab', 'lawngreen', 'magenta', 'dodgerblue']

    # Show the major grid lines with dark grey lines
    plt.grid(b=True, which='major', color='#666666', linestyle='-')

    # Show the minor grid lines with very faint and almost transparent grey lines
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)

    title1, title2, title02 = title_func(file_name0, head)

    y_max = np.nanmax(spectr1)
    y_min = np.nanmin(spectr1)
    x_min = argument.min()
    x_max = argument.max()

    if flag:
        ax.set_xlabel('Freq, MHz', fontsize=18)
        ax.set_yscale('log')
        ax.set_ylabel('Power Spectrum', fontsize=20)
        ax.set_title(title02 + title1, fontsize=18)
        y1 = y_min * 2
        y2 = y_min
        y3 = y_max - (y_max - y_min) / 10
        y4 = y_min * 4
    else:
        # pylab.xlim(x_min, x_max + 100)
        plt.legend(loc='upper right')
        ax.set_xlabel('Time, sec', fontsize=18)
        if burn == 1:
            ax.set_yticks([0, 1])
            ax.set_ylim(0, 4)
        ax.set_ylabel('Intensity', fontsize=18)
        ax.set_title(title2 + ' scan ' + title1, fontsize=20)
        y1 = y_max - 3 * (y_max - y_min) / 10
        y2 = y_max - 4 * (y_max - y_min) / 10
        y3 = y_max - (y_max - y_min) / 10
        y4 = y_max - (y_max - y_min) / 5

    plt.text(x_min, y1, inform[0], fontsize=16)  # Разрешение по частоте
    plt.text(x_min, y2, inform[1], fontsize=16)  # Разрешение по времени
    plt.text(x_min, y3, inform[2], fontsize=16)  # Информация о поляризации
    plt.text(x_min, y4, inform[3], fontsize=16)  # Информация о статистической чистке сканов
    m = 0
    for i in range(freq_line_sp1):
        ax.plot(argument, spectr1[i, :], color=line_color[m], label=line_legend[i])
        m += 1

    # Управление шрифтом легенды
    font = font_manager.FontProperties(family='Comic Sans MS',
                                       weight='bold',
                                       style='normal', size=16)
    ax.legend(prop=font)
    ax.legend(loc=10, prop=font, bbox_to_anchor=(1, 0.5))

    plt.show()

    add_pass1 = path_to_pic(file_name0 + '\\', flag)
    fig.savefig(Path(file_name0, add_pass1))


def title_func(file_name0, head):
    az = file_name0[-3:]
    att1 = str(head['att1'])
    att2 = str(head['att2'])
    att3 = str(head['att3'])
    date = head['date'][:-1]

    title1 = date + ', az = ' + az + ', Att = [' + att1 + ', ' + att2 + ', ' + att3 + ']'
    if not file_name0.find('sun') == -1:
        title2 = 'Sun intensity'
        title02 = 'Sun spectrum '
        if file_name0[-1:] == 'b':
            az = file_name0[-4:-1]
            title2 = 'Calibration'
            title02 = 'Calibration spectrum '
            title1 = date + ', az = ' + az + ', Black Body/Sky, Att = [' + att1 + ', ' + att2 + ', ' + att3 + ']'

    elif not file_name0.find('crab') == -1:
        title2 = 'Crab intensity'
        title02 = 'Crab spectrum '

    elif not file_name0.find('calibration') == -1:
        title2 = 'Calibration'
        title02 = 'Calibration spectrum '
        title1 = date + ', az = ' + az + ', Att = [' + att1 + ', ' + att2 + ', ' + att3 + ']'

    elif not (file_name0.find('test') == -1):
        kind = file_name0[-2:]
        title2 = 'Test'
        title02 = 'Test spectrum'
        if kind == 'VG':
            power_vg = 0
            title1 = date + ', Vector Gen' + ', P = ' + str(power_vg) + 'dBm, ' \
                'Att = [' + att1 + ', ' + att2 + ', ' + att3 + ']'
        elif kind == 'NG':
            t_noise = 6300
            title1 = date + ', Noise Gen, ' + 'T = ' + str(t_noise) + 'K, Att = [' + att1 + ', ' + att2 + ', ' + att3 + ']'
        elif kind == 'ML':
            title1 = date + ', Matched Load' ', Att = [' + att1 + ', ' + att2 + ', ' + att3 + ']'
        elif kind == 'SC':
            title1 = date + ', Short Cut' ', Att = [' + att1 + ', ' + att2 + ', ' + att3 + ']'
        else:
            title1 = date + ', ' + ' , Att = [' + att1 + ', ' + att2 + ', ' + att3 + ']'

    else:
        title2 = 'Scan'
        title02 = 'Spectrum'

    return title1, title2, title02


def fig_multi_axes(spectr1, argument, inform, file_name0path, freq_mask, head, n_row_pic=2, n_col_pic=3):
    """ Функция строит многооконный рисунок, принимает numpy-массив в качестве значений (первый индекс массива -
    число графиков-окон), аргумент, информацию о свойствах графиков (разрешение по времени и частоте, поляризация), путь
    к файлу данных, характеристики файла данных head, список параметра freq_mask, каждому значению из которого
    соответствует свой график-окно. Максимальное количество окон задается n_row_pic, n_col_pic, и по умолчанию -
     2 * 3 = 6.
     Выдает рисунок как объект Figure matplotlib и сохраняет его в формате png в папку одноименную с file_name0path в
     той же директории, где находится файл данных"""

    file_name0 = str(file_name0path)
    size_sp1 = spectr1.shape
    freq_line_sp1 = size_sp1[0]

    title1, title2, title02 = title_func(file_name0, head)

    fig, axes = plt.subplots(n_row_pic, n_col_pic, figsize=(12, 6))

    # line_colore = ['green', 'blue', 'purple', 'lime', 'black', 'red', 'olivedrab', 'lawngreen',
    # 'magenta', 'dodgerblue']

    fig.suptitle(title2 + ' scan' + title1, y=1.0, fontsize=24)

    for i_freq in range(freq_line_sp1):
        axes[i_freq // n_col_pic, i_freq % n_col_pic].plot(argument, spectr1[i_freq, :])

        axes[i_freq // n_col_pic, i_freq % n_col_pic].set_title(str(freq_mask[i_freq]) + ' MHz')
        # Show the major grid lines with dark grey lines
        axes[i_freq // n_col_pic, i_freq % n_col_pic].grid(b=True, which='major', color='#666666', linestyle='-')
        # Show the minor grid lines with very faint and almost transparent grey lines
        axes[i_freq // n_col_pic, i_freq % n_col_pic].minorticks_on()
        axes[i_freq // n_col_pic, i_freq % n_col_pic].grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)
        # axes[i_phantom // 3, i_phantom % 3].set_yticks([0, 1])
        # axes[i_phantom // 3, i_phantom % 3].set_ylim(0, 4)
        # axes[i_phantom // 3, i_phantom % 3].set_xlabel('t, sec', fontsize=10)

        xticks = axes[i_freq // n_col_pic, i_freq % n_col_pic].get_xticks().tolist()
        xticks[-2:] = ''
        axes[i_freq // n_col_pic, i_freq % n_col_pic].xaxis.set_major_locator(FixedLocator(xticks))
        # Установление формата отображения меток по оси 0у - при порядке выше 4 он выносится на верх оси
        sf = ScalarFormatter()
        sf.set_powerlimits((-4, 4))
        axes[i_freq // n_col_pic, i_freq % n_col_pic].yaxis.set_major_formatter(sf)
        # *************************************************************************

        axes[i_freq // n_col_pic, i_freq % n_col_pic].xaxis.set_label_coords(1.05, -0.025)
        axes[i_freq // n_col_pic, i_freq % n_col_pic].annotate('t, sec', xy=(0.95, -0.05), ha='left', va='top',
                                               xycoords='axes fraction', fontsize=10)

    axes[-1, -1].axis('off')
    fig.text(0.68, 0.25, inform[2], fontsize=14)
    fig.text(0.68, 0.2, inform[0], fontsize=14)
    fig.text(0.68, 0.15, inform[1], fontsize=14)
    # plt.figtext(0.11, 0.25, 'time')
    # axes[0, 2].legend(loc=10, bbox_to_anchor=(1.2, 0.5))
    plt.subplots_adjust(hspace=0.4)
    plt.show()

    # Управление шрифтом легенды
    font = font_manager.FontProperties(family='Comic Sans MS',
                                       weight='bold',
                                       style='normal', size=16)
    add_pass1 = path_to_pic(file_name0, 0)
    fig.savefig(Path(file_name0, add_pass1))


def path_to_pic(file_path, flag, format='png'):
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
    if not os.path.isfile(Path(file_path, add_pass1)):
        pass
    else:
        while os.path.isfile(Path(file_path, add_pass1)):
            num = int(add_pass0[l - 2:l]) + 1
            num_str = str(num)
            if num >= 10:
                add_pass0 = add_pass0[:l - 2] + num_str
            else:
                add_pass0 = add_pass0[:l - 2] + '0' + num_str
            add_pass1 = add_pass0 + '.' + format

    return add_pass1


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


if __name__ == '__main__':
    fig_plot()
    pass
