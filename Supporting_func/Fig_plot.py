import matplotlib.pyplot as plt
import numpy as np
import matplotlib.font_manager as font_manager
import pylab
import os
# from IPython.display import set_matplotlib_formats
# set_matplotlib_formats('svg')

''' На входе программы несколько пар строк с данными и координатами, причем количество точек, в общем случае, разное для 
    разных пар. В программу передается легенда и название осей. 
    Программа создает сетку, надписывает оси, дает легенду.
'''


def fig_plot(spectr1, burn, argument, flag, inform, file_name0, line_legend=[2] * 20):

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

    y_max = spectr1.max()
    y_min = spectr1.min()
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


def fig_multi_axes(spectr1, argument, inform, file_name0, freq_mask):

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
