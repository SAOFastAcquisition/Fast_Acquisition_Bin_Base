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

    _sp_matched[0] = _sp_matched[0][-1::-1]   # Правильный порядок отсчетов
    _sp_short[0] = _sp_short[0][-1::-1]       # Правильный порядок отсчетов для второй зоны Найквиста

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

    _fig, _ax = plt.subplots(1, figsize=(12, 6))
    _ax.plot(_data['temperature'][0])
    _ax.plot(_data['temperature'][1])
    _ax.plot(_data['temperature'][2])
    _ax.plot(_data['temperature'][3])
    plt.grid()
    plt.show()
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
    pass
    return _s


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
    receiver_temperature_path = Path(head_path, 'Alignment', receiver_temperature_file_name)

    date = current_primary_file1[0:10]

    delta_t = 8.3886e-3
    _delta_f = 7.8125
    aver_param_noise = 8
    temp0 = 300
    time_mask = [11, 19, 21, 29]

    #           ********* Загружаем результаты двух измерений **********
    path1 = Path(converted_data_file_path, current_primary_file1 + '_spectrum.npy')
    spectrum1 = np.load(path1, allow_pickle=True)
    path2 = Path(converted_data_file_path, current_primary_file2 + '_spectrum.npy')
    spectrum2 = np.load(path2, allow_pickle=True)
    #           Граничные отсчеты
    t0, t1 = int(time_mask[0] / delta_t), int(time_mask[1] / delta_t)
    t2, t3 = int(time_mask[2] / delta_t), int(time_mask[3] / delta_t)
    num_mask = [int(t / delta_t) for t in time_mask]

    #  Отклик приемника на КЗ и согласованную нагрузку каналов левой и правой поляризаций
    data = differ(spectrum1, num_mask)
    data = data.append(differ(spectrum2, num_mask), ignore_index=True)
    #                            *****************
    a1 = np.array(data['spectrum'][data['polar'] == 'left'][data['load'] == 'matched'].iloc[0])
    b1 = np.array(data['spectrum'][data['polar'] == 'left'][data['load'] == 'short'].iloc[0])
    a2 = np.array(data['spectrum'][data['polar'] == 'right'][data['load'] == 'matched'].iloc[0])
    b2 = np.array(data['spectrum'][data['polar'] == 'right'][data['load'] == 'short'].iloc[0])

    # Расчет и сохранение в файл шумовой температуры приемника
    temp_left = b1 / (2 * a1 - b1) * temp0
    temp1 = pd.Series((date, temp_left, 'left'), index=['date', 'temperature', 'polar'])
    receiver_temp_update(temp1, ['date', 'temperature', 'polar'])
    temp_right = b2 / (2 * a2 - b2) * temp0
    temp2 = pd.Series((date, temp_right, 'right'), index=['date', 'temperature', 'polar'])
    receiver_temp_update(temp2, ['date', 'temperature', 'polar'])
    #                       ****************************
    fig, ax = plt.subplots(1, figsize=(12, 6))
    ax.plot(temp_left)
    ax.plot(temp_right)
    # ax.plot(data['spectrum'][2])
    # ax.plot(data['spectrum'][3])
    plt.grid()
    plt.show()
    pass
