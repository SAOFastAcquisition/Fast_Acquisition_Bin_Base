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
spectrum1 = np.array(spectrum[3])
f_spectrum = fits.PrimaryHDU(spectrum1)
pass
print(f_spectrum)
f_spectrum.update('NAXIS1', 512, "frequency axis")
# hdulist.info()

print(repr(f_spectrum.header))
print(list(f_spectrum.header.keys()))
a = f_spectrum.header
b = f_spectrum.data
print(f_spectrum.data)
m = f_spectrum.data.shape
print(f'shape of f_spectrum.data is: {m}')
plt.imshow(f_spectrum.data[:, :], origin='lower')
plt.show()

head = {'date': date,
        'measure_kind': measure_kind,    # Вид измерений: наблюдение Солнца, Луны, калибровка АЧХ
        'band_size': band_size,  # Параметр 'whole' означает работу в диапазоне 1-3 ГГц,
        # 'half_low' - диапазон 1-2, 'half_upper' - 2-3 ГГц
        'polar': polar,  # Принимает значения поляризаций: 'both', 'left', 'right'
        'cleaned': 'no',
        'n_aver': n_aver,
        'shift': shift,
        'kurtosis': bound_left,
        'att1': att01,
        'att2': att02,
        'att3': att03,
        'align_file_path': r'F:\Fast_Acquisition\Alignment\Align_coeff.bin',
        'align_coeff_pos': 5}