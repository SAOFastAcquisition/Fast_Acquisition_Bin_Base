from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import pickle
from Supporting_func.stokes_coefficients import path_to_data

def symple_plot(_data):
    _data_shape = np.shape(_data)

    _fig, _ax = plt.subplots(1, figsize=(12, 6))
    for _i in range(_data_shape[0]):
        _ax.plot(_data[_i, :])

    plt.show()


def fits_reader_andrew():
    _fits_file = Path(r"D:\YandexDisk\ИТ\Солнце\Нага\2024-08-01_121957_sun+00_out_1.fits")

    try:
        with fits.open(_fits_file) as _hdul:
            _primary_hdu = _hdul[0]
            _compressed_image_hdu = _hdul[1]
            _table_hdu = _hdul[2]

            _header = _primary_hdu.header
            _data_array = _compressed_image_hdu.data
            _table_data = _table_hdu.data
    except Exception as e:
        print(f"Error {e}")
    print(_hdul.info())
    _hdr = _hdul[0].header
    print(_hdr.cards)
    return _header, _data_array, _table_data


if __name__ == '__main__':
    # fits_image_filename = '20240212_122749_sun+0_out.fits'
    # fits_image_filename = r'D:\YandexDisk\ИТ\Солнце\Нага\2024-08-01_121957_sun+00_out.fits'
    # hdul = fits.open(fits_image_filename)
    # print(hdul.info())
    # hdr = hdul[0].header
    # print(f'hdr[0] = {hdr[5]}')
    # a = hdul[0].data
    # b = hdul[1].data
    # c = hdr.cards
    # a0 = a[:, 0, :]
    # a1 = a[:, 1, :]
    # print(c)
    # symple_plot(a0)
    # hdul.close()

    a, b, c = fits_reader_andrew()
    b0 = b[:, 0, :]
    b1 = b[:, 1, :]
    symple_plot(b0)
    pass
