from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from Supporting_func.stocks_coefficients import path_to_data


current_data_file = '2021-06-28_07-04'      # Имя файла с исходными текущими данными без расширения
current_data_dir = '2021_06_28sun'          # Папка с текущими данными
align_file_name = 'Align_coeff.bin'         # Имя файла с текущими коэффициентами выравнивания АЧХ
current_catalog = r'2021/Results'           # Текущий каталог (за определенный период, здесь - год)

file_path_data, head_path = path_to_data(current_catalog, current_data_dir)
# file_name = Path(file_path_data, current_data_file + '.npy')

spectrum = np.load(Path(file_path_data, current_data_file + '_spectrum.npy'), allow_pickle=True)
spectrum1 = np.array(spectrum[2])
hdulist = fits.PrimaryHDU(spectrum1)
pass
print(hdulist)
hdulist.info()
hdu = hdulist[2]
print(repr(hdu.header))
print(list(hdu.header.keys()))
a = hdu.header
b = hdu.data
print(hdu.data)
m = hdu.data.shape
print(m)
plt.imshow(hdu.data[:,:], origin='lower')
plt.show()
