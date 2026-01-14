import os
import glob as gb
import pandas as pd
import numpy as np
from pathlib import Path
import pickle
import gzip, shutil
import matplotlib.pyplot as plt
import plotly, plotly.offline as plo
from plotly.graph_objs import Scatter, Layout
import plotly.graph_objs as go

# from Data_extraction.c2023_12_15_03p16 import a
from stokes_coefficients import some_fig_plot, title_func
from sun_az_spead import sun_az_speed
from Help_folder.paths_via_class import DataPaths

dt = 8.3886e-3
df = 2000 / 512


def load_stokes(_path='2023-10-25_05-24_stocks.npy'):
    if _path[-2:] == 'gz':
        with gzip.open(Path(_path), "rb") as fin:
            _data = np.load(fin, allow_pickle=True)
    else:
        _data = np.load(Path(_path), allow_pickle=True)
    _stokes_I = _data[0]
    _stokes_V = _data[1]
    _time_counter = _data[2]
    # plt.plot(_y0[:, 150])
    # plt.grid('both')
    # plt.show()
    return _stokes_I, _stokes_V, _time_counter


def save_reformat_stokes(_path_obj, _path_stokes_base):
    """

    :param _path_obj:
    :param _path_stokes_base:
    :return:
    """
    _paths = gb.glob(str(Path(_path_obj.converted_dir_path, "*stocks.npy.gz")))  # Set Pattern in glob() function
    if not _paths:
        _paths = gb.glob(str(Path(_path_obj.converted_dir_path, "*stocks.npy")))
    if not _paths:
        print('save_reformat_stokes: No Stokes coefficient files in folder')
        exit(1)

    for _path_curr in _paths:
        _data_file = os.path.basename(_path_curr)
        _stokes_I, _stokes_V, _time_count = load_stokes(_path_curr)
        _dict_stokes = {
            '_date': _data_file[0:10],
            'azimuth': _data_file[13:16],
            'stokes_I': [_stokes_I],
            'stokes_V': [_stokes_V],
            'time_count': [_time_count]
        }
        _frame_stokes_curr = pd.DataFrame(_dict_stokes)

        if not os.path.isfile(_path_stokes_base):
            _stokes_base = _frame_stokes_curr
        else:
            if os.path.isfile(_path_stokes_base):
                with open(_path_stokes_base, 'rb') as _inp:
                    _stokes_base = pickle.load(_inp)
            elif os.path.isfile(f'{_path_stokes_base}.gz'):
                with gzip.open(f'{_path_stokes_base}.gz', 'rb') as _inp:
                    _stokes_base = pickle.load(_inp)

        _idx = _stokes_base.loc[(_stokes_base._date == _dict_stokes['_date'])
                                & (_stokes_base.azimuth == _dict_stokes['azimuth'])].index  #
        if not len(_idx):
            _stokes_base = pd.concat([_stokes_base, _frame_stokes_curr], axis=0, ignore_index=False)
        else:
            pass

        with open(_path_stokes_base, 'wb') as out:
            pickle.dump(_stokes_base, out)

    with open(_path_stokes_base, "rb") as fin, gzip.open(f'{_path_stokes_base}.gz', "wb") as f_out:
        # Читает файл по частям, экономя оперативную память
        shutil.copyfileobj(fin, f_out)


def plot_intensities(_arg, _x_l, _x_r, _param, _angle=0):
    _fig = plt.figure(figsize=[13, 8])
    #                           ****** Рисунок 1 левая поляризация ******
    _ax1 = plt.subplot(2, 1, 1)
    for _a, _b in zip(_x_l.T, _param):
        plt.plot(_arg, _a, label=f'left: freq = {_b} MHz')

    _ax1.set_ylim()
    plt.grid('both')
    plt.title(str(path_treatment)[-10:] + ' Intensities L & R')
    _ax1.set_ylabel('Intensity, K')
    plt.legend(loc=2)
    plt.text((_arg[0] + _arg[-1]) / 2, 1.16, f'Position angle on Sun {np.ceil(_angle)} arcs', size=12)
    plt.grid(which='major', color='#666666', linestyle='-')
    plt.grid(which='minor', color='#999999', linestyle='-', alpha=0.5)
    plt.minorticks_on()

    #                           ****** Рисунок 2 правая поляризация ******
    _ax2 = plt.subplot(2, 1, 2)
    for _a, _b in zip(_x_r.T, _param):
        plt.plot(_arg, _a, label=f'right: freq = {_b} MHz')
    #                                       ************
    plt.grid('on')
    _ax2.set_ylabel('Intensity, K')
    _ax2.set_ylim(ymax=1.2)
    plt.xlabel('Frequency, MHz')
    plt.legend(loc=2)
    plt.grid(which='major', color='#666666', linestyle='-')
    plt.grid(which='minor', color='#999999', linestyle='-', alpha=0.5)
    plt.minorticks_on()
    plt.show()

    return _fig, path_treatment, 7, 'png'


def plotly_pic(_x, _y):
    _y0 = _y[:, 50]
    _n = [20, 50, 55]
    _f = [1000 + df / 2 + i * df for i in range(512)]
    fig = go.Figure()
    i = 0
    for s in _y.T:
        fig.add_trace(go.Scatter(
            x=_x[:, 0], y=s,
            mode='lines',
            name=f'freq={_f[i]:4.0f}'
        ))
        i += 1
    # Set the _y-axis to log scale
    fig.update_layout(
        yaxis_type='log',
        xaxis_title='X',
        yaxis_title='Y',
        title='Log Scale Plot'
    )
    # fig.write_image("plot.png", format='png')
    fig.write_html("plot.html")
    # Show the plot
    fig.show()



def filter_freq(_i):
    _freq_0 = freq_mask(_i)
    df = 7.825 / 2
    _freq = np.array([1000. + df * (0.5 + _k) for _k in range(512)])
    _k = 0
    for _s in _freq_0:

        filt_cur = (_freq == _freq[np.where(_freq >= _s)[0][0]])
        if not _k:
            _filt_freq = filt_cur
        else:
            _filt_freq = filt_cur | _filt_freq
        _k += 1

    return _filt_freq


def filter_azimuth(_azimuth, _base):
    _k = 0
    for _s in _azimuth:
        filt_cur = _base['azimuth'].astype('int') == _s
        if not _k:
            _filt_az = filt_cur
        else:
            _filt_az = filt_cur | _filt_az
        _k += 1

    return _filt_az


def filter_position(_time_num, _angle0):
    _time_seq = _time_num * dt
    _angle_seq = time_to_angle(_time_seq)
    _k = 0
    for _s in _angle0:
        filt_cur = (_angle_seq == _angle_seq[np.where(_angle_seq >= _s)[0][0]])
        if not _k:
            _filt_pos = filt_cur
        else:
            _filt_pos = filt_cur | _filt_pos
        _k += 1

    return _filt_pos


def time_to_angle(_time_count, _date, _az, _time_culmination=200):
    if _date:
        _scale = sun_az_speed(_date, _az)
        if _scale[0] == 1:
            _scale = 1900 / 135
    else:
        _scale = 1950 / 135
    _time = np.array(_time_count) * dt
    _angle = [-(t - _time_culmination) * _scale for t in _time][-1::-1]
    _angle = np.array(_angle)
    return _angle


def freq_mask(_i):
    _n1 = 1
    _n2 = 5
    _freq_mask = [
        [1736],  # [0]
        [2060, 2300, 2500, 2750, 2830, 2920],  # [1]
        [1020, 1100, 1200, 1300, 1350, 1400, 1450, 1600],  # [2]
        [1000 * _n1 + 100 * _n2 + 90 + 5 * i for i in range(10)],  # [3]
        [1050, 1465, 1535, 1600, 1700, 2265, 2550, 2700, 2800, 2920],  # [4]
        [1230, 1560, 2300, 2910],  # [5]
        [1140, 1420, 1480, 2460, 2500, 2780],  # for Crab '2021-06-28_03+14' # [6]
        [1220, 1540, 1980, 2060, 2500, 2780],  # for Crab '2021-06-28_04+12' # [7]
        [1200, 1300, 1465, 1600, 1700, 2265, 2510, 2720, 2800, 2920]  # [8]
    ]
    return _freq_mask[_i]


if __name__ == '__main__':
    """
    Скрипт объединяет в один DataFrame все ранее вычисленные в данном азимуте параметры Стокса и номера 
    временных отсчетов, соответствующих середине пачек взятия отсчетов (примерно по 30) поляризаций поочередно
    """
    current_data_file = '2024-02-14_03+16'
    date = '2024-02-14'
    main_dir = date[0:4]
    data_dir = f'{date[0:4]}_{date[5:7]}_{date[8:]}sun'

    path_obj = DataPaths(date, data_dir, main_dir)
    path_stokes_base = Path(path_obj.converted_dir_path, 'stokes_base.npy')
    path_to_stokes_fig_folder = Path(path_obj.treatment_dir_path, current_data_file)
    path_treatment = path_obj.treatment_data_file_path

    with open(Path(path_obj.converted_dir_path, current_data_file + '_head.bin'), 'rb') as inp:
        head = pickle.load(inp)

    reformat_stokes = 'n'
    pic_demonstrate = '_y'
    az = [-24, -20, -16]

    if reformat_stokes == '_y':
        save_reformat_stokes(path_obj, path_stokes_base)

    if pic_demonstrate == '_y':
        if os.path.isfile(f'{path_stokes_base}.gz'):
            with gzip.open(f'{path_stokes_base}.gz', "rb") as inp:
                stokes_base = pickle.load(inp)
        elif os.path.isfile(path_stokes_base):
            with open(path_stokes_base, 'rb') as inp:
                stokes_base = pickle.load(inp)
        else:
            print('main: Base not found')
            quit('Base not found')

        filt_freq = filter_freq(8)
        filt_az = filter_azimuth(az, stokes_base)
        # В качестве индекса в отобранных данных используем значение азимута
        filtered_data = stokes_base[filt_az].set_index(stokes_base[filt_az].azimuth)
        stokes_I = filtered_data['stokes_I']
        stokes_V = filtered_data['stokes_V']

        i = 0
        for sI, sV in zip(stokes_I, stokes_V):
            sI_filt, sV_filt = sI[:, filt_freq], sV[:, filt_freq]
            s_left, s_right = (sI_filt + sV_filt) / 2, (sI_filt - sV_filt) / 2
            if not i:
                s_left_int = pd.Series([s_left], index=[str(az[i])])
                s_right_int = pd.Series([s_right], index=[str(az[i])])
            else:
                s_left_int[str(az[i])] = s_left
                s_right_int[str(az[i])] = s_right
            i += 1
    si = stokes_I[stokes_I.index.astype('int') == az[2]][0]
    sv = stokes_V[stokes_I.index.astype('int') == az[2]][0]
    with open(Path(path_obj.converted_dir_path, current_data_file + '_head.bin'), 'rb') as inp:
        head = pickle.load(inp)
    time_count = filtered_data[filtered_data['azimuth'].astype(int) == az[0]].time_count[0]
    arg = time_to_angle(time_count, date, az)
    # some_fig_plot(path_to_stokes_fig_folder, arg, si[-1::-1], sv[-1::-1])
    # plot_intensities(arg, si[-1::-1], sv[-1::-1], freq_mask(8))
    plotly_pic(arg, si[-1::-1])
    pass
