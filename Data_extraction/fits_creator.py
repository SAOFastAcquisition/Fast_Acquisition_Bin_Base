from astropy.io import fits
import numpy as np
from pathlib import Path
import os
import gzip
import ast
import pickle
import glob as gb
from Help_folder.paths_via_class import DataPaths
from Supporting_func.stokes_coefficients import time_to_angle

df = 2000 / 512
dt = 8.3886e-3


def head_creator(head_fits, head_data):
    #
    _fits_words = import_fits_word()
    _fits_words['OFF_TIME'] = _fits_words.pop('FEED_OFFSET_TIME')
    _a = _fits_words.keys()
    for _s in _a:
        head_fits[_s] = (_fits_words[_s])

    head_fits['DATE'] = (head_data['date'], 'Date of observation')
    head_fits['TELESCOP'] = ('RATAN_600', 'DM complex')
    # head_fits['MEASKIND'] = (head_data['measure_kind'], 'Object of observation')
    head_fits['UNITS'] = ('sfu', 'Data units')
    head_fits['BAND'] = band(head_data['band_size'])
    head_fits['POLAR'] = (head_data['polar'], 'Polarization')
    head_fits['CLEAN'] = ('no', 'Additional data cleaning')
    head_fits['ALIGNPA'] = (head_data['align_file_path'], 'Align coefficients')
    head_fits['ALIGNPOS'] = (head_data['align_coeff_pos'], 'Align coefficients position')

    head_fits['DTIME'] = (8.3886e-3, 'Time resolution, s')
    head_fits['DFREQ'] = (df, 'Frequency resolution, MHz')

    head_fits['KURTOSIS'] = (head_data['kurtosis'], 'Half wide of kurtosis interval')
    head_fits['ATT1'] = (head_data['att1'], 'Common attenuation')
    head_fits['ATT2'] = (head_data['att2'], 'Attenuation in left channel')
    head_fits['ATT3'] = (head_data['att3'], 'Attenuation in right channel')
    head_fits['COMMENT'] = ': NAXIS1 - IV-representation, I - ind 0, V - ind 1'

    return head_fits


def head_f_creator(_head_fits):
    _head_fits['COMMENT'] = ': Frequency (counts in MHz)'
    return _head_fits


def head_pos_creator(_head_fits):
    _head_fits['COMMENT'] = ': Position on Sun (counts in arcs)'
    return _head_fits


def import_fits_word():
    _path = gb.glob(str(Path(path_obj.primary_dir_path, f"*{data_file[-3:]}.desc")))[0]
    with open(_path) as file:
        d = file.read()
    ind = d.find('{')
    # Словарь с параметрами наблюдения
    res_dict = ast.literal_eval(d[ind:])

    return res_dict['fits_words']


def band(string):
    if string == 'whole':
        size = '1-3 GGz'
    elif string == 'half_low':
        size = '1-2 GGz'
    elif string == 'half_upper':
        size = '2-3 GGz'
    return size


def control(_path):
    # Проверка
    # Загрузили fits-file
    hdu_list = fits.open(_path)
    # Информация о составе
    hdu_list.info()
    # Первая (главная) составляющая - спектры
    hdu_r = hdu_list[0]
    _data = hdu_r.data
    head0 = hdu_r.header
    print(*repr(head0))
    # Отсчеты частоты
    hdu_f = hdu_list[1]
    _data1 = hdu_f.data
    head1 = hdu_f.header
    print(*repr(head1))
    # Отсчеты позиции на Солнце
    hdu_pos = hdu_list[2]
    _data2 = hdu_pos.data
    head2 = hdu_list[2].header
    print(*repr(head2))
    # plt.imshow(f_spectrum.data[:, :, 0], origin='lower')
    # plt.show()
    # plt.imshow(f_spectrum.data[:, :, 1], origin='lower')
    # plt.show()


if __name__ == '__main__':

    data_file = '2024-01-05_11-16'
    main_dir = data_file[0:4]
    data_dir = f'{data_file[0:4]}_{data_file[5:7]}_{data_file[8:10]}sun'

    path_obj = DataPaths(data_file, data_dir, main_dir)
    path_stokes = Path(str(path_obj.converted_data_file_path) + '_stocks.npy')
    path_head = Path(str(path_obj.converted_data_file_path) + '_head.bin')
    path_fits = Path(str(path_obj.converted_data_file_path) + '.fits')

    if os.path.exists(f'{str(path_stokes)}.gz'):
        filename_out = f'{str(path_stokes)}.gz'
        with gzip.open(filename_out, "rb") as fin:
            _spectrum = np.load(fin, allow_pickle=True)
    else:
        _spectrum = np.load(path_stokes, allow_pickle=True)

    with open(path_head, 'rb') as inp:
        head = pickle.load(inp)
    # for key, value in head.items():
    #     print(f'head key :: {key} \tvalue::\t{value}')
    spectrum0 = np.array(_spectrum[0][::-1, :])
    sp1 = spectrum0.reshape(spectrum0.shape[0], spectrum0.shape[1], -1)
    spectrum1 = np.array(_spectrum[1][::-1, :])
    sp2 = spectrum1.reshape(spectrum1.shape[0], spectrum0.shape[1], -1)
    sp = np.concatenate((sp1, sp2), axis=2, out=None, dtype=None, casting="same_kind")
    _time = np.array(_spectrum[2]) * dt
    pos_count, s = time_to_angle(_time, sp2, data_file[:10], int(data_file[-3::]))
    f_count = np.array([1000 + df / 2 + df * i for i in range(512)])
    # pos_count = [1, 2]

    f_spectrum = fits.PrimaryHDU(sp)
    f_spectrum.header = head_creator(f_spectrum.header, head)

    f_count_hdu = fits.ImageHDU(f_count)
    f_count_hdu.header = head_f_creator(f_count_hdu.header)

    pos_count_hdu = fits.ImageHDU(pos_count)
    pos_count_hdu.header = head_pos_creator(pos_count_hdu.header)

    hdu_list = fits.HDUList([f_spectrum, f_count_hdu, pos_count_hdu])
    hdu_list.writeto(path_fits, overwrite=True)

    control(path_fits)

