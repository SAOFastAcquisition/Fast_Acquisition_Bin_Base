import os
import numpy as np
import sys
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from pathlib import Path
from Supporting_func.stocks_coefficients import path_to_data


def differ(_s, _num_mask):
    [_t0, _t1, _t2, _t3] = _num_mask
    sp_av1 = [np.mean(s[_t0:_t1, :], axis=0) for s in _s]
    sp_av2 = [np.mean(s[_t2:_t3, :], axis=0) for s in _s]
    a = sp_av1[0][500]
    b = sp_av2[2][500]
    if sp_av1[0][500] > 10 and sp_av1[0][500] > sp_av2[2][500]:
        _sp_matched = [sp_av1[0], sp_av1[1]]
        _sp_short = [sp_av2[2], sp_av2[3]]
        flag_matched = 'left'
    elif (sp_av1[0][500] < 10) and (sp_av1[2][500] > sp_av2[0][500]):
        _sp_matched = [sp_av1[2], sp_av1[3]]
        _sp_short = [sp_av2[0], sp_av2[1]]
        flag_matched = 'right'
    elif (sp_av1[0][500] < 10) and (sp_av1[2][500] < sp_av2[0][500]):
        _sp_matched = [sp_av2[0], sp_av2[1]]
        _sp_short = sp_av1[2], sp_av1[3]
        flag_matched = 'left'
    else:
        flag_matched = 'right'
        _sp_matched = [sp_av2[2], sp_av2[3]]
        _sp_short = [sp_av1[0], sp_av1[1]]

    if flag_matched == 'right':
        flag_short = 'left'
    else:
        flag_short = 'right'

    _sp_matched[0] = _sp_matched[0][-1::-1]  # Правильный порядок отсчетов
    _sp_short[0] = _sp_short[0][-1::-1]  # Правильный порядок отсчетов для второй зоны Найквиста

    _sp_m = np.hstack((_sp_matched[0], _sp_matched[1]))
    _sp_s = np.hstack((_sp_short[0], _sp_short[1]))

    _sp_m = zone_deletion(_sp_m, 'matched')
    _sp_s = zone_deletion(_sp_s, 'short')
    columns_names = ['date', 'spectrum', 'polar', 'load']

    _data = pd.DataFrame([[date, _sp_m, flag_matched, 'matched'],
                          [date, _sp_s, flag_short, 'short']],
                         columns=columns_names)

    # fig, ax = plt.subplots(1, figsize=(12, 6))
    # ax.plot(_data['spectrum'][0])
    # ax.plot(_data['spectrum'][1])
    # # ax.plot(_data['spectrum'][2])
    # # ax.plot(_data['spectrum'][3])
    # plt.grid()
    # plt.show()

    return _data


def receiver_temp_update(_s, _columns_names):
    """ Принимает объект Series, и, если в данных файла нет среза идентичного принимаемому,
    добавляет принимаемый объект в конец файла"""

    if not os.path.isfile(receiver_temperature_path):
        _data = pd.DataFrame(columns=_columns_names)
    else:
        with open(receiver_temperature_path, 'rb') as inp:
            _data = pickle.load(inp)

    idx_r = _data.loc[(_data['date'] == _s['date'])
                      & (_data['polar'] == _s['polar'])].index
    if not len(idx_r):
        _data = _data.append(_s, ignore_index=True)
        with open(receiver_temperature_path, 'wb') as _out:
            pickle.dump(_data, _out)

    # _fig, _ax = plt.subplots(1, figsize=(12, 6))
    # _ax.plot(_data['temperature'][0])
    # _ax.plot(_data['temperature'][1])
    # _ax.plot(_data['temperature'][2])
    # _ax.plot(_data['temperature'][3])
    # plt.grid()
    # plt.show()
    pass


def zone_deletion(_s, _flag):
    # Исключение зон действия режекторных фильтров
    k1 = int((1090 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))
    k2 = int((1220 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))
    k3 = int((1520 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))
    k4 = int((1700 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))
    k6 = 1024
    k5 = int((770 - _delta_f / 2 / aver_param_noise) // (_delta_f / aver_param_noise))

    if _flag == 'matched':
        _s[0:20] = 100
        _s[1024:1058] = 100
        _s[-46:] = 100
        _s[k1:k2] = 100
        _s[k3:k4] = 100
        _s[k5:k6] = 100
    else:
        _s[0:20] = 50
        _s[1024:1058] = 50
        _s[-38:] = 50
        _s[k1:k2] = 50
        _s[k3:k4] = 50
        _s[k5:k6] = 50

    return _s


def del_random(_s):
    _l = len(_s)
    for _i in range(1, _l - 2):
        if _s[_i] < 10 and _s[_i - 1] > 20 and _s[_i + 1] > 20:
            _s[_i] = (_s[_i - 1] + _s[_i + 1]) / 2
        if _s[_i] < 0:
            _s[_i] = 100
    return _s


def del_random_mod(_s, _s0):
    _l = len(_s)
    for _i in range(1, _l - 2):
        if abs(_s[_i] - _s[_i - 1]) > 2 * abs(_s[_i + 1] - _s[_i - 1]):
            _s[_i] = (_s[_i - 1] + _s[_i + 1]) / 2
        if _s[_i] <= 0.01:
            _s[_i] = 10
    _s[0] = _s0
    _s[_l - 2:] = _s0
    return _s


def receiver_temperature_calc(_data):
    a1 = np.array(_data['spectrum'][_data['polar'] == 'left'][_data['load'] == 'matched'].iloc[0])
    b1 = np.array(_data['spectrum'][_data['polar'] == 'left'][_data['load'] == 'short'].iloc[0])
    a2 = np.array(_data['spectrum'][_data['polar'] == 'right'][_data['load'] == 'matched'].iloc[0])
    b2 = np.array(_data['spectrum'][_data['polar'] == 'right'][_data['load'] == 'short'].iloc[0])

    # Расчет и сохранение в файл шумовой температуры приемника
    temp_left = b1 / (2 * a1 - b1) * temp0
    temp_left = del_random_mod(temp_left, 100)
    temp1 = pd.Series((date, temp_left, 'left'), index=['date', 'temperature', 'polar'])
    receiver_temp_update(temp1, ['date', 'temperature', 'polar'])
    temp_right = b2 / (2 * a2 - b2) * temp0
    temp_right = del_random_mod(temp_right, 100)
    temp2 = pd.Series((date, temp_right, 'right'), index=['date', 'temperature', 'polar'])
    receiver_temp_update(temp2, ['date', 'temperature', 'polar'])

    #       ***** Вызов функции расчета выравнивающих АЧХ коэффициентов *****
    align_coefficients_calc(a1, a2, b1, b2)

    fig, ax = plt.subplots(1, figsize=(12, 6))
    ax.plot(temp_left)
    ax.plot(temp_right)
    # # ax.plot(data['spectrum'][2])
    # # ax.plot(data['spectrum'][3])
    plt.grid()
    plt.show()
    return


def align_coefficients_calc(_a1, _a2, _b1, _b2):
    left_term_noise = _a1 - _b1 / 2
    right_term_noise = _a2 - _b2 / 2
    max_left = np.max(left_term_noise)
    max_right = np.max(right_term_noise)
    left_term_noise /= max_left
    right_term_noise /= max_right
    left_term_noise = del_random_mod(left_term_noise, 10)
    right_term_noise = del_random_mod(right_term_noise, 10)

    align_coeff_left = np.array([1 / a for a in left_term_noise])
    align_coeff_right = np.array([1 / a for a in right_term_noise])
    align_coefficient_update(align_coeff_left, align_coeff_right)
    fig, ax = plt.subplots(1, figsize=(12, 6))
    ax.plot(align_coeff_left)
    ax.plot(align_coeff_right)
    # # ax.plot(data['spectrum'][2])
    # # ax.plot(data['spectrum'][3])
    plt.grid()
    plt.show()

    return


def align_coefficient_update(_data1, _data2):
    with open(alignment_coeff_path, 'rb') as inp:
        _calibration_frame = pickle.load(inp)

    _columns_names = ['date', 'att1', 'att2', 'att3',
                      'spectrum_left1', 'polar', 'spectrum_left2', 'spectrum_right1',
                      'spectrum_right2',
                      'max_left1', 'max_left2', 'max_right1', 'max_right2', 'flag_align']

    _l = int(len(_data2))
    _l //= 2
    sp_left1 = _data1[:_l]
    sp_left1 = sp_left1[-1::-1]
    sp_left2 = _data1[_l: int(2 * _l)]
    sp_right1 = _data2[: _l]
    sp_right1 = sp_right1[-1::-1]
    sp_right2 = _data2[_l:int(2 * _l)]
    idx = _calibration_frame.loc[(_calibration_frame.date == date)
                                 & (_calibration_frame.att1 == att1)
                                 & (_calibration_frame.att2 == att2)
                                 & (_calibration_frame.att3 == att3)
                                 & (_calibration_frame.polar == 'both')].index
    if not len(idx):
        calibrate_row_ser = pd.Series((date, att1, att2, att3, sp_left1, 'both', sp_left2,
                                       sp_right1, sp_right2, 1, 1, 1, 1, 1), index=_columns_names)
        _calibration_frame = _calibration_frame.append(calibrate_row_ser, ignore_index=True)

        with open(alignment_coeff_path, 'wb') as out:
            pickle.dump(_calibration_frame, out)


if __name__ == '__main__':
    """ Для расчета температуры собственных шумов приемника используются два измерения отклика приемника:
    с согласованной нагрузкой на входе и КЗ на входе. На один вход подключается согласованная нагрузка, 
    на второй - КЗ. В промежуток времени time_mask[0] - time_mask[1] 
    включена одна поляризация (левая, как правило), в промежуток времени time_mask[2] - time_mask[3] - 
    другая (правая). Затем нагрузка на входах меняется местами и выполняется второе измерение. Т.о. для каждого 
    входа (левая и правая поляризации) получаем отклик приемника на согласованную нагрузку и КЗ"""

    current_data_dir = '2022'
    primary_data_dir = 'Primary_data'  # Каталог исходных данных (за определенный период, здесь - год)
    converted_data_dir = 'Converted_data'  # Каталог для записи результатов конвертации данных и заголовков
    data_treatment_dir = 'Data_treatment'  # Каталог для записи результатов обработки, рисунков

    current_primary_dir = '2022_06_27test'
    current_converted_dir = current_primary_dir + '_conv'
    current_converted_path = Path(converted_data_dir, current_converted_dir)
    current_treatment_dir = current_primary_dir + '_treat'
    current_treatment_path = Path(data_treatment_dir, current_treatment_dir)

    ngi_temperature_file_name = 'ngi_temperature.npy'  # Файл усредненной по базе шумовой температуры для ГШ
    receiver_temperature_file_name = 'receiver_temperature.npy'
    current_primary_file1 = '2022-06-27_02'  # Файл с согласованной нагрузкой и КЗ на входах приемника
    current_primary_file2 = '2022-06-27_03'  # Файл с КЗ и согласованной нагрузкой на входах приемника

    converted_data_file_path, head_path = path_to_data(current_data_dir, current_converted_path)
    data_treatment_file_path, head_path = path_to_data(current_data_dir, current_treatment_path)

    ngi_temperature_path = Path(head_path, 'Alignment', ngi_temperature_file_name)
    alignment_coeff_path = Path(head_path, 'Alignment', 'Align_coeff.bin')
    receiver_temperature_path = Path(head_path, 'Alignment', receiver_temperature_file_name)

    date = current_primary_file1[0:10]

    delta_t = 8.3886e-3
    _delta_f = 7.8125
    aver_param_noise = 8
    temp0 = 300
    time_mask = [11, 19, 21, 29]

    #           ********* Загружаем результаты двух измерений **********
    path1 = Path(converted_data_file_path, current_primary_file1 + '_spectrum.npy')
    path2 = Path(converted_data_file_path, current_primary_file2 + '_spectrum.npy')
    spectrum1 = np.load(path1, allow_pickle=True)
    spectrum2 = np.load(path2, allow_pickle=True)
    with open(Path(converted_data_file_path, current_primary_file1 + '_head.bin'), 'rb') as inp:
        head = pickle.load(inp)
    att1, att2, att3 = head['att1'], head['att2'], head['att3']

    #           Граничные отсчеты
    t0, t1 = int(time_mask[0] / delta_t), int(time_mask[1] / delta_t)
    t2, t3 = int(time_mask[2] / delta_t), int(time_mask[3] / delta_t)
    num_mask = [int(t / delta_t) for t in time_mask]

    #  Отклик приемника на КЗ и согласованную нагрузку каналов левой и правой поляризаций
    data = differ(spectrum1, num_mask)
    data = data.append(differ(spectrum2, num_mask), ignore_index=True)
    #                            *****************
    receiver_temperature_calc(data)

    pass
