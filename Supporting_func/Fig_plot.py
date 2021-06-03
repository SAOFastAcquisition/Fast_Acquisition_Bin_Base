import matplotlib.pyplot as plt
import numpy as np
import matplotlib.font_manager as font_manager
from matplotlib.ticker import MaxNLocator
import pylab
import os
from pathlib import Path

# from IPython.display import set_matplotlib_formats
# set_matplotlib_formats('svg')

''' На входе программы несколько пар строк с данными и координатами, причем количество точек, в общем случае, разное для 
    разных пар. В программу передается легенда и название осей. 
    Программа создает сетку, надписывает оси, дает легенду.
'''


def fig_plot(spectr1, burn, argument, flag, inform, file_name0_path, line_legend=[2] * 20):
    file_name0 = str(file_name0_path)
    size_sp1 = spectr1.shape
    for i in range(size_sp1[0]):
        for j in range(size_sp1[1]):
            if spectr1[i, j] < 100:
                spectr1[i, j] = 'NaN'

    freq_line_sp1 = size_sp1[0]

    # title0 = file_name0[-19:-2]
    # title1 = '  '+title0[0:4]+'.'+title0[4:6]+'.'+title0[6:8] +\
    #          ' time='+title0[9:11]+':'+title0[11:13]+' azimuth='+title0[14:17]
    title0 = file_name0[-22:-2]
    title1 = '  ' + title0[0:4] + '.' + title0[4:6] + '.' + title0[6:8] + \
             ' mode='+title0[9:13]+' attenuation='+title0[14:16]+' azimuth='+title0[17:20]
    fig, ax = plt.subplots(1, figsize=(12, 6))

    line_colore = ['green', 'blue', 'purple', 'lime', 'black', 'red', 'olivedrab', 'lawngreen', 'magenta', 'dodgerblue']

    # Show the major grid lines with dark grey lines
    plt.grid(b=True, which='major', color='#666666', linestyle='-')

    # Show the minor grid lines with very faint and almost transparent grey lines
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)

    y_max = np.nanmax(spectr1)
    y_min = np.nanmin(spectr1)
    x_min = argument.min()
    x_max = argument.max()

    if not file_name0.find('sun') == -1:
        title2 = 'Sun intensity'
    elif not file_name0.find('crab') == -1:
        title2 = 'Crab intensity'
    elif not file_name0.find('calibr') == -1:
        title2 = 'Calibration'
        title0 = file_name0[-23:-2]
        title1 = '  ' + title0[0:4] + '.' + title0[4:6] + '.' + title0[6:8] + \
                 ' channel att=' + title0[14:17] + ' source att=' + title0[18:21]

    elif not (file_name0.find('Calibrate') == -1):
        if not (file_name0.find('Ant1') == -1):
            title2 = 'Calibration Left_Pol'
        if not (file_name0.find('Ant2') == -1):
                title2 = 'Calibration Right_Pol'
        if not (file_name0.find('pol2') == -1):
            title2 = 'Calibration L&R Pol' # 20201224_Ant1_HotL_3
        title0 = file_name0[-20:-2]
        title1 = '  ' + title0[0:4] + '.' + title0[4:6] + '.' + title0[6:8] + \
                 ' mode=' + title0[9:13] + ' source=' + title0[14:18]

    elif not (file_name0.find('test') == -1):
        title2 = 'Interference Test'
        title0 = file_name0[-24:-2]
        title1 = '  ' + title0[0:4] + '.' + title0[4:6] + '.' + title0[6:8] + \
                 ' channel att=' + title0[15:18] + ' source dBm=' + title0[19:22]
        pass
    else:
        title2 = []

    if flag:
        ax.set_xlabel('Freq, MHz', fontsize=20)
        ax.set_yscale('log')
        ax.set_ylabel('Power Spectrum', fontsize=20)
        ax.set_title('Spectrum'+title1, fontsize=24)
        y1 = y_min + 0.15 * (y_max - y_min) / 10
        y2 = y_min
        y3 = y_max - (y_max - y_min) / 10
    else:
        # pylab.xlim(x_min, x_max + 100)
        plt.legend(loc='upper right')
        ax.set_xlabel('Time, sec', fontsize=20)
        if burn == 1:
            ax.set_yticks([0, 1])
            ax.set_ylim(0, 4)
        ax.set_ylabel('Intensity', fontsize=20)
        ax.set_title(title2+' scan'+title1, fontsize=24)
        y1 = y_min + (y_max - y_min) / 10
        y2 = y_min
        y3 = y_max - (y_max - y_min) / 10

    plt.text(x_min, y1, inform[0], fontsize=16)
    plt.text(x_min, y2, inform[1], fontsize=16)
    plt.text(x_min, y3, inform[2], fontsize=16)


    m = 0
    for i in range(freq_line_sp1):
        ax.plot(argument, spectr1[i, :], color=line_colore[m], label=line_legend[i])
        m += 1

    # Управление шрифтом легенды
    font = font_manager.FontProperties(family='Comic Sans MS',
                                       weight='bold',
                                       style='normal', size=16)
    # ax.legend(prop=font)
    ax.legend(loc=10, prop=font, bbox_to_anchor=(1, 0.5))

    plt.show()

    add_pass1 = path_to_pic(file_name0 + '\\', flag)
    fig.savefig(file_name0 + '\\' + add_pass1)


def fig_multi_axes(spectr1, argument, inform, file_name0path, freq_mask):
    file_name0 = str(file_name0path)
    size_sp1 = spectr1.shape
    freq_line_sp1 = size_sp1[0]

    title0 = file_name0[-19:-2]
    title1 = '  '+title0[0:4]+'.'+title0[4:6]+'.'+title0[6:8] +\
             ' time='+title0[9:11]+':'+title0[11:13]+' azimuth='+title0[14:17]

    fig, axes = plt.subplots(3, 4, figsize=(12, 6))

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

    for i_freq in range(freq_line_sp1):
        axes[i_freq // 4, i_freq % 4].plot(argument, spectr1[i_freq, :])

        axes[i_freq // 4, i_freq % 4].set_title(str(freq_mask[i_freq]) + ' MHz')
        # Show the major grid lines with dark grey lines
        axes[i_freq // 4, i_freq % 4].grid(b=True, which='major', color='#666666', linestyle='-')
        # Show the minor grid lines with very faint and almost transparent grey lines
        axes[i_freq // 4, i_freq % 4].minorticks_on()
        axes[i_freq // 4, i_freq % 4].grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)
        # axes[i_phantom // 3, i_phantom % 3].set_yticks([0, 1])
        # axes[i_phantom // 3, i_phantom % 3].set_ylim(0, 4)
        # axes[i_phantom // 3, i_phantom % 3].set_xlabel('t, sec', fontsize=10)
        axes[i_freq // 4, i_freq % 4].xaxis.set_label_coords(1.05, -0.025)

        xticks = axes[i_freq // 4, i_freq % 4].get_xticks().tolist()
        xticks[-2:] = ''
        axes[i_freq // 4, i_freq % 4].set_xticklabels(xticks)
        axes[i_freq // 4, i_freq % 4].annotate('t, sec', xy=(0.95, -0.05), ha='left', va='top',
                                               xycoords='axes fraction', fontsize=10)

    axes[0, 2].legend(loc=10, bbox_to_anchor=(1.2, 0.5))
    plt.subplots_adjust(hspace=0.4)
    plt.show()

    # Управление шрифтом легенды
    font = font_manager.FontProperties(family='Comic Sans MS',
                                       weight='bold',
                                       style='normal', size=16)
    add_pass1 = path_to_pic(file_name0, 0)
    fig.savefig(Path(file_name0, add_pass1))


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
    if not os.path.isfile(Path(file_path, add_pass1)):
        pass
    else:
        while os.path.isfile(Path(file_path, add_pass1)):
            num = int(add_pass0[l-2:l]) + 1
            num_str = str(num)
            if num >= 10:
                add_pass0 = add_pass0[:l - 2] + num_str
            else:
                add_pass0 = add_pass0[:l-2] + '0' + num_str
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