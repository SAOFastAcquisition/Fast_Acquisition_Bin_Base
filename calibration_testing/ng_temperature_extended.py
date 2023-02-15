import os
import numpy as np
import sys
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from pathlib import Path
from Supporting_func.stocks_coefficients import path_to_data


def temperature_ngi(_spectrum, polarization, time_borders, _path):
    '''

    :param _spectrum:    Калибровочные данные с участием эталонного ГШ
    :param polarization: Поляризация калибровочных данных _spectrum
    :param time_borders: Времена включения/переключения/выключения источников шума на входе приемника:
                         эталонный ГШ, согласованная нагрузка, внутренний ГШ (time0, time1, time2, time3)
    :param _path: Путь для сохранения калибровочных по эталонному ГШ коэффициентов перехода к
                  антенным температурам
    :return:
    '''
    _delta_t = 8.3886e-3
    _delta_f = 7.8125
    aver_param_noise = 8
    df0 = _delta_f / aver_param_noise
    time0, time1, time2, time3 = time_borders
    n0 = int((time0 + 1) / _delta_t)  # На входе радиометра сигнал от эталонного ГШ
    n1 = int((time1 - 1) / _delta_t)
    n2 = int((time1 + 1) / _delta_t)  # На входе радиометра сигнал от согласованной нагрузки
    n3 = int((time2 - 1) / _delta_t)
    n4 = int((time2 + 1) / _delta_t)  # На входе радиометра сигнал от внутреннего ГШ
    n5 = int((time3 - 1) / _delta_t)

    if polarization == 'left':
        s0, s1 = _spectrum[0], _spectrum[1]
    elif polarization == 'right':
        s0, s1 = _spectrum[2], _spectrum[3]
    elif polarization == 'both':
        s0, s1 = _spectrum[0], _spectrum[1]

    # Вычисление коэффициентов привязки температуры шумового сигнала от внутреннего
    # ГШ к шумовой температуре эталонного ГШ
    av_gne0 = np.mean(s0[n0:n1, :], axis=0)
    av_gne1 = np.mean(s1[n0:n1, :], axis=0)
    av_ml0 = np.mean(s0[n2:n3, :], axis=0)
    av_ml1 = np.mean(s1[n2:n3, :], axis=0)
    av_gni0 = np.mean(s0[n4:n5, :], axis=0)
    av_gni1 = np.mean(s1[n4:n5, :], axis=0)
    # Коэффициент перехода от температуры эталонного ГШ к температуре калибровочного ГШ
    temperature_scale0 = (av_gni0 - av_ml0) / av_gne0
    temperature_scale1 = (av_gni1 - av_ml1) / av_gne1
    # Коэффициент перехода от температуры эталонного ГШ к температуре (собственные шумы усилителя + согласованная нагр)
    temperature_scale_rt0 = av_ml0 / av_gne0
    temperature_scale_rt1 = av_ml1 / av_gne1

    # Исключение зон действия режекторных и антиэлайзинговых фильтров
    # Третья зона Найквиста
    k0 = int((90 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))
    k1 = int((220 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))
    k2 = int((540 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))
    k3 = int((700 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))
    # Вторая зона Найквиста
    k4 = 0
    k5 = int((230 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))
    k_zone = [k0, k1, k2, k3, k4, k5]

    temperature_scale0, temperature_scale1 = filter_zone(temperature_scale0, temperature_scale1, k_zone)
    temperature_scale_rt0, temperature_scale_rt1 = filter_zone(temperature_scale_rt0, temperature_scale_rt1, k_zone)

    #                   *************************************
    #   Формирование массивов значений шумовой температуры эталонного источника от частоты
    f0 = np.array([1000 + df0 / 2 + i * df0 for i in range(1024)])
    f1 = np.array([2000 + df0 / 2 + i * df0 for i in range(1024)])
    temp_nge0 = [temperature_nge(f) for f in f0]
    temp_nge1 = [temperature_nge(f) for f in f1]
    #   Шумовая температура внутреннего источника шума
    temp_ngi0 = temp_nge0 * temperature_scale0
    temp_ngi1 = temp_nge1 * temperature_scale1
    temp_ngi0 = del_random_mod(temp_ngi0, 400)
    temp_ngi1 = del_random_mod(temp_ngi1, 400)
    #   Собственная шумовая температура приемника
    temp_rec0 = temp_nge0 * temperature_scale_rt0 - 300
    temp_rec1 = temp_nge1 * temperature_scale_rt1 - 300
    temp_rec0 = del_random_mod(temp_rec0, 20)
    temp_rec1 = del_random_mod(temp_rec1, 20)

    #                   ***********************************************
    #          ******* Калибровка антенной температуры по эталонному ГШ *******
    #                   ***********************************************
    #   Определение коэффициентов пересчета произвольных единиц в антенную температуру
    c0, c1 = ant_calibr_coeff(av_gne0, av_gne1, temp_nge0, temp_nge1, k_zone)
    #   Запись коэффициентов пересчета в файл
    save_ant_calibr_coeff(c0, c1, polarization, _path)
    #                   ***********************************************

    plt.plot(f0, c0)
    plt.plot(f1, c1)
    # plt.plot(f0, temp_rec0)
    # plt.plot(f1, temp_rec1)

    plt.grid()
    plt.show()

    return temp_ngi0, temp_ngi1


def filter_zone(_scale0, _scale1, _k):
    """
    Функция обнуляет коэффициенты в зонах действия режекторных и антиэлайзинговых фильтров
    :param _scale0: Отсчеты частоты во второй зоне Найквиста
    :param _scale1: Отсчеты частоты в первой зоне Найквиста
    :param _k: Зоны действия режекторных фильтров
    :return:
    """
    _scale1[_k[0]:_k[1]] = 0
    _scale1[_k[2]:_k[3]] = 0
    _scale0[_k[4]:_k[5]] = 0
    _scale1[0:32] = 0  # Антиэлайзинговый фильтр
    _scale1[-46:] = 0  # Антиэлайзинговый фильтр
    _scale0[-20:] = 0  # Антиэлайзинговый фильтр

    _scale0[0:] = _scale0[-1::-1]  # Правильный порядок отсчетов для второй зоны Найквиста

    return _scale0, _scale1


def ant_calibr_coeff(_s0, _s1, _t_nge0, _t_nge1, _k):
    """
    Функция отдает коэффициенты перехода к антенным температурам. Отсчет частоты во второй зоне - обратный.
    Принимает усредненные по времени наблюдения значения спектра эталонного ГШ в произвольных единицах во второй и
    третьей зонах Найквиста и шумовую температуру ГШ. Отсчет частоты во второй зоне - обратный.
    :param _s0:     Усредненный спектр во второй зоне Найквиста
    :param _s1:     Усредненный спектр в третьей зоне Найквиста
    :param _t_nge0: Шумовая температура эталонного ГШ во второй зоне
    :param _t_nge1: Шумовая температура эталонного ГШ в первой зоне
    :param _k:      Границы действия режекторных фильтров в отсчетах по зонам Найквиста
    :return:        Коэффициенты перехода к антенным температурам
    """
    _c0 = [0] * 1024
    _c1 = [0] * 1024
    _t_nge0 = _t_nge0[-1::-1]
    _c0[_k[5]:-20] = _t_nge0[_k[5]:-20] / _s0[_k[5]:-20]
    _c1[32:_k[0]] = _t_nge1[32:_k[0]] / _s1[32:_k[0]]
    _c1[_k[1]:_k[2]] = _t_nge1[_k[1]:_k[2]] / _s1[_k[1]:_k[2]]
    _c1[_k[3]:-46] = _t_nge1[_k[3]:-46] / _s1[_k[3]:-46]
    # _c0 = _c0[-1::-1]
    pass
    return _c0, _c1


def save_ant_calibr_coeff(_c0, _c1, _polarisation, _path):
    """
    Функция сохраняет в файл atc_file_name = 'antenna_temperature_coefficients.npy' коэффициенты перехода к
    антенным температурам с учетом поляризации _c0, _c1
    :param _c0:             Коэффициенты для второй зоны Найквиста
    :param _c1:             Коэффициенты для первой зоны Найквиста
    :param _polarisation:   Поляризация коэффициентов
    :param _path:           Путь для сохранения
    :return:
    """

    _columns_names = ['date', 'att1', 'att2', 'att3', 'polar',
                      'spectrum_left1', 'spectrum_left2', 'spectrum_right1',
                      'spectrum_right2', 'flag_align']
    if not os.path.isfile(_path):
        calibration_frame = pd.DataFrame(columns=_columns_names)
    else:
        with open(_path, 'rb') as _inp:
            calibration_frame = pickle.load(_inp)

    align_coeff = [None, None, None, None]
    if _polarisation == 'left':
        align_coeff[0] = _c0
        align_coeff[1] = _c1
    else:
        align_coeff[2] = _c0
        align_coeff[3] = _c1
    flag_align = 0
    # Рассчитанные коэффициенты вместе с исходной информацией записываем в виде словаря  для формирования
    # объекта Series и включения в сводную  таблицу корректирующих коэффициентов
    calibrate_row = {'date': current_data_file[:10], 'att1': head['att1'], 'att2': head['att2'], 'att3': head['att3'],
                     'polar': head['polar'], 'spectrum_left1': align_coeff[0], 'spectrum_left2': align_coeff[1],
                     'spectrum_right1': align_coeff[2], 'spectrum_right2': align_coeff[3],
                     'flag_align': flag_align}  # 'flag_align' - признак выравнивания по всему диапазону : 1 - сделано
    calibrate_row_ser = pd.Series(calibrate_row)

    # Определяем, есть ли в сводной таблице данные с такими же исходными параметрами. Если есть, то будет их
    # проверка на то, содержат они все коэффициенты или нет. Если нет, то объект Series будет полностью
    # вставлен в таблицу (объект DataFrame)

    _idx = calibration_frame.loc[(calibration_frame.date == head['date'])
                                 & (calibration_frame.att1 == head['att1'])
                                 & (calibration_frame.att2 == head['att2'])
                                 & (calibration_frame.att3 == head['att3'])
                                 ].index  # & (calibration_frame.polar == head['polar'])

    if len(_idx):
        r = calibration_frame.iloc[_idx[0]]
        # Определяем, есть ли пустые поля в выделенной из таблицы строке. Если есть, то эти поля будут
        # заполнены из имеющегося объекта  Series
        ax_bool = r.isnull()
        r[ax_bool] = calibrate_row_ser[ax_bool]
        # Остались ли в r незаполненные поля
        ax_bool = r.isnull()

        if ('True' not in ax_bool) and (r['flag_align']) == 0:
            r['flag_align'] = 1
            r['polar'] = 'both'
            # Если пустые поля заполнились ('True' not in ax_bool) и ранее эта строчка не была полностью сформирована
            # ('flag_align': 0), то строка с номером idx[0] будет удалена из объекта DataFrame и вставлена новая,
            # полностью заполненная r
            calibration_frame = calibration_frame.drop(_idx).append(r, ignore_index=True)
    else:
        calibration_frame = calibration_frame.append(calibrate_row_ser, ignore_index=True)
    # calibration_frame = calibration_frame.drop(axis=0, index=[2, 3])
    with open(_path, 'wb') as _out:
        pickle.dump(calibration_frame, _out)
        print('Коэффициенты перехода к антенным температурам сохранены')


def temperature_nge(_f, t0=300):
    """ Функция возвращает шумовую температуру внешнего поверочного ГШ в зависимости от
    частоты в МГц, t0 - температура окружающей среды ГШ (по умолчанию 300К) На отрезках
    между известными значениями ИОШТ - линейная аппроксимация"""
    n_1000 = 32.06  # ИОШТ для 1000 МГц
    n_2000 = 31.44  # ИОШТ для 2000 МГц
    n_3000 = 29.89  # ИОШТ для 3000 МГц
    k1 = n_1000 - n_2000
    k2 = n_2000 - n_3000

    if 1000 <= _f <= 2000:
        t_nge = (n_1000 - k1 * (_f / 1000 - 1)) * t0
    elif 2000 < _f <= 3000:
        t_nge = (n_2000 - k2 * (_f / 1000 - 2)) * t0

    return t_nge


def average_ngi(_data):
    """ Усредняет по результатам серии измерений"""
    x, y = _data['low_band'], _data['upper_band']
    l = len(x)
    index_x, index_y = x.index, y.index
    array_x, array_y = x[index_x[0]], y[index_y[0]]
    for i in index_x[1:]:
        array_x = np.vstack([array_x, x[i]])
        array_y = np.vstack([array_y, y[i]])
    av_data_x = np.mean(array_x, axis=0)
    # Обнуление недостоверных значений из-за действия фильтров в начале и конце зон Найквиста
    av_data_x[0:20] = 0
    av_data_y = np.mean(array_y, axis=0)
    av_data_y[0:34] = 0
    av_data_y[-34:] = 0

    return np.hstack([av_data_x, av_data_y])


def del_random_mod(_s, _s0):
    _l = len(_s)
    for _i in range(1, _l - 2):
        if abs(_s[_i] - _s[_i - 1]) > 2 * abs(_s[_i + 1] - _s[_i - 1]):
            _s[_i] = (_s[_i - 1] + _s[_i + 1]) / 2
        if _s[_i] <= 0.01:
            _s[_i] = _s0
    _s[0] = _s0
    _s[_l - 2:] = _s0
    return _s


def ngi_temperature_update(_ngi_selected, _ngi_id, _case_id):
    if not os.path.isfile(ngi_temperature_path):
        _ngi_temperature = pd.DataFrame(columns=columns_names)
    else:
        with open(ngi_temperature_path, 'rb') as inp:
            _ngi_temperature = pickle.load(inp)

    _temperature = average_ngi(_ngi_selected)
    ngi_average_row = {'ngi_id': ngi_id, 'case_id': _case_id, 'polar': head['polar'],
                       'temperature': _temperature, 'note': note}

    idx_r = _ngi_temperature.loc[(_ngi_temperature.ngi_id == '01')
                                 & (_ngi_temperature.polar == head['polar'])
                                 & (_ngi_temperature.case_id == _case_id)].index

    if len(idx_r):
        _ngi_temperature = _ngi_temperature.drop(idx_r).append(ngi_average_row, ignore_index=True)
    else:
        _ngi_temperature = _ngi_temperature.append(ngi_average_row, ignore_index=True)
    with open(ngi_temperature_path, 'wb') as _out:
        pickle.dump(_ngi_temperature, _out)
    pass


def plot_ngi(_data):
    # Arguments0, s1 = _spectrum[0], _spectrum[1]
    _delta_f = 7.8125
    aver_param_noise = 8
    df0 = _delta_f / aver_param_noise
    f0 = np.array([1000 + df0 / 2 + i * df0 for i in range(1024)])
    f1 = np.array([2000 + df0 / 2 + i * df0 for i in range(1024)])
    f = np.hstack([f0, f1])
    x, y = _data['low_band'][_data['polar'] == 'left'], _data['upper_band'][_data['polar'] == 'left']
    index_x, index_y = x.index, y.index
    array_x, array_y = x[index_x[0]], y[index_y[0]]
    f01 = f0
    f11 = f1
    for i in index_x[1:]:
        f01 = np.vstack([f01, f0])
        f11 = np.vstack([f11, f1])
        array_x = np.vstack([array_x, x[i]])
        array_y = np.vstack([array_y, y[i]])

    x1, y1 = _data['low_band'][_data['polar'] == 'right'], _data['upper_band'][_data['polar'] == 'right']
    index_x1, index_y1 = x1.index, y1.index
    array_x1, array_y1 = x1[index_x1[0]], y1[index_y1[0]]
    f02 = f0
    f12 = f1
    for i in index_x1[1:]:
        f02 = np.vstack([f02, f0])
        f12 = np.vstack([f12, f1])
        array_x1 = np.vstack([array_x1, x1[i]])
        array_y1 = np.vstack([array_y1, y1[i]])

    with open(ngi_temperature_path, 'rb') as _in:
        ngi_temperature = pickle.load(_in)
    idx_case = ngi_temperature[ngi_temperature.case_id == case_id].index
    ngi_left = ngi_temperature['temperature'][ngi_temperature['polar'] == 'left'][idx_case[1]]
    ngi_right = ngi_temperature['temperature'][ngi_temperature['polar'] == 'right'][idx_case[0]]

    fig, ax = plt.subplots(1, figsize=(12, 6))
    ax.plot(f, ngi_left)
    ax.plot(f, ngi_right)
    ax.scatter(f01, array_x)
    ax.scatter(f11, array_y)
    ax.scatter(f02, array_x1)
    ax.scatter(f12, array_y1)
    # line_minor = plt.plot(inp_scale / 1000)
    # plt.plot(inp_scale, inp_data, 'r--*', inp_scale, 'g:+')  # В кавычках компактно заданы цвет, стиль линии и маркеры
    # plt.setp(fig_main, linestyle=':', color='r')
    # plt.setp(fig, linestyle=':', color=(0, 1, 0, 0.9), marker='.', markerfacecolor='b')  # Четвертое число
    # задает прозрачность линии
    # plt.setp(line_minor, linestyle='-.')
    plt.grid()
    plt.show()

    pass


if __name__ == '__main__':
    ''' 
    Создаем базу данных по измеренным шумовым температурам от внутреннего ГШ на входе радиометра
    для каждого экземпляра ГШ ***nge_temperature_base_name***, среднее значение шумовой температуры
    от этого ГШ по всем измерениям записывается в файл ***nge_temperature_file_name***
        Этот скрипт расширен расчетом температуры собственных шумов приемника по шумовой температуре на 
    нагруженном на согласованное сопротивление входе. Эталонным ГШ устанавливаем цену шкалы и из шумовой 
    температуры на входе вычитаем 300К
        Создана база коэффициентов пересчета спектров в произвольных единицах к антенным температурам с помощью
    эталонного ГШ, подключаемого ко входу приемника. Функция расчета: ant_calibr_coeff()
    функция создания файла и сохранения коэффициентов: save_ant_calibr_coeff()
        Все измерения выполняются с разрешением 1 МГц
    '''

    current_data_file = '2022-11-24_08'                        # Имя файла с исходными текущими данными без расширения
    current_primary_dir = '2022_11_24test'
    current_data_dir = current_primary_dir + '_conv'               # Папка с текущими данными
    current_catalog = r'2022/Test_and_calibration/Converted_data'  # Текущий каталог (за определенный период,
    # здесь - год)

    ngi_temperature_base_name = 'ngi_temperature_base1.npy' # Базы данных по температуре ГШ на входе МШУ
    ngi_temperature_file_name = 'ngi_temperature1.npy'      # Файл усредненной по базе шумовой температуры для ГШ
    atc_file_name = 'antenna_temperature_coefficients.npy'  # Файл коэффициентов пересчета к антенной температуре
    file_path_data, head_path = path_to_data(current_catalog, current_data_dir)

    #                           ***********************
    # Создание или загрузка файла с шумовой температурой встроенного генератора шума
    ngi_id = '01'                                           # Идентификатор встроенного генератора шума
    case_id = '03'      # Идентификатор набора данных, по которым усредняется температура собств. шумов приемника
    note = '2022_11_18-24'                                  # Дата или период измерений набора данных
    timing = [0, 20, 30, 50]
    ngi_temperature_base_path = Path(head_path, 'Alignment', ngi_temperature_base_name)
    ngi_temperature_path = Path(head_path, 'Alignment', ngi_temperature_file_name)
    atc_path = Path(head_path, 'Alignment', atc_file_name)  # atc - antenna temperature coefficients
    columns_names_base = ['nge_id', 'date', 'att1', 'att2', 'att3',
                          'polar', 'attempt_num', 'low_band', 'upper_band']
    columns_names = ['ngi_id', 'polar', 'temperature', 'case_id', 'note']
    if not os.path.isfile(ngi_temperature_base_path):
        ngi_temperature_base = pd.DataFrame(columns=columns_names_base)
    else:
        with open(ngi_temperature_base_path, 'rb') as inp:
            ngi_temperature_base = pickle.load(inp)

    #                           ************************

    spectrum = np.load(Path(file_path_data, current_data_file + '_spectrum.npy'), allow_pickle=True)
    with open(Path(file_path_data, current_data_file + '_head.bin'), 'rb') as inp:
        head = pickle.load(inp)
    r = temperature_ngi(spectrum, head['polar'], timing, atc_path)
    attempt_num = current_data_file[11:13]
    idx = ngi_temperature_base.loc[(ngi_temperature_base.nge_id == ngi_id)
                                   & (ngi_temperature_base.date == head['date'])
                                   & (ngi_temperature_base.att1 == head['att1'])
                                   & (ngi_temperature_base.att2 == head['att2'])
                                   & (ngi_temperature_base.att3 == head['att3'])
                                   & (ngi_temperature_base.polar == head['polar'])
                                   & (ngi_temperature_base.attempt_num == attempt_num)].index

    if not len(idx):
        r = temperature_ngi(spectrum, head['polar'], timing)
        temperature_row = {'nge_id': ngi_id, 'date': head['date'],
                           'att1': head['att1'], 'att2': head['att2'], 'att3': head['att3'],
                           'polar': head['polar'], 'attempt_num': attempt_num,
                           'low_band': r[0], 'upper_band': r[1]
                           }
        ngi_temperature_base = ngi_temperature_base.append(temperature_row, ignore_index=True)
        with open(ngi_temperature_base_path, 'wb') as out:
            pickle.dump(ngi_temperature_base, out)
    else:
        print('Данные этой записи уже внесены в базу')

    #               *** Усреднение по базе данных ***
    # шумовой температуры встроенного ГШ на входе Радиометра
    # раздельно для каналов левой и правой поляризаций.
    pw = input('Обновить усредненную шумовую температуру ГШ на входе РМ?\t')
    idx_drop = ngi_temperature_base.loc[ngi_temperature_base.att3 == 15].index
    if len(idx_drop):
        for n in idx_drop:
            ngi_temperature = ngi_temperature_base.drop([n])
    if pw == 'y':
        ngi_selected1 = ngi_temperature_base[
            ngi_temperature_base['nge_id'].isin([ngi_id])]  # [ngi_temperature_base['date']
        # == '2022-11-24']
        head['polar'] = 'right'
        ngi_selected1l = ngi_selected1[ngi_selected1['polar'].isin([head['polar']])]
        ngi_temperature_update(ngi_selected1l, ngi_id, case_id)
        head['polar'] = 'left'
        ngi_selected1l = ngi_selected1[ngi_selected1['polar'].isin([head['polar']])]
        ngi_temperature_update(ngi_selected1l, ngi_id, case_id)

    #               *** Picture ***
    ngi_selected1 = ngi_temperature_base[ngi_temperature_base['nge_id'].isin([ngi_id])]  # [ngi_temperature_base['date']
    # == '2022-11-24']
    plot_ngi(ngi_selected1)

    # *******************************************
    #       *** Запись новых данных или обновление старых ***

    with open(ngi_temperature_base_path, 'rb') as inp:
        ngi_temperature_base_wr = pickle.load(inp)
    pass
