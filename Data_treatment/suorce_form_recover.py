import numpy as np
from pathlib import Path
from matplotlib import pyplot as plt
from scipy.fftpack import fft, ifft
from Help_folder.paths_via_class import DataPaths


def simplest_fig(_x, _y, _z):
    fig, axes = plt.subplots(2, 1, figsize=(12, 12))
    axes[0].plot(_x, _y)
    axes[1].plot(_x, _z)
    axes[0].set_ylabel('Stocks_I')
    axes[1].set_ylabel('Stocks_V', color='darkred')
    axes[0].minorticks_on()
    axes[1].minorticks_on()
    axes[0].legend(loc='upper right')
    axes[0].grid()
    axes[1].grid()
    axes[0].grid(which='minor',
                 axis='x',
                 color='k',
                 linestyle=':')
    axes[1].grid(which='minor',
                 axis='x',
                 color='k',
                 linestyle=':')
    plt.show()


if __name__ == "__main__":
    current_primary_file = '2022-12-23_01+28'
    current_primary_dir = '2022_12_23sun'
    main_dir = '2022'
    # main_dir = r'2021/Results'           # Каталог (за определенный период, здесь - за 2021 год)
    date = current_primary_dir[0:10]
    adr1 = DataPaths(current_primary_file, current_primary_dir, main_dir)
    converted_data_file_path = adr1.converted_data_file_path
    data_treatment_file_path = adr1.treatment_data_file_path
    data_saved_path = Path(data_treatment_file_path, 'Meshed_Spectrum', current_primary_file + '_meshed.npy')

    kt = 4
    delta_t0 = 8.3886e-3
    delta_t = delta_t0 * kt

    spectrum = np.load(data_saved_path)
    spectrum_one = spectrum[:, 24]
    shape_spectrum = np.shape(spectrum)

    n_freq = shape_spectrum[1]
    delta_f = 2000/n_freq   # unit = 'MHz'
    time_to_angle_coeff = 180 / 1920
    time = [i * delta_t for i in range(shape_spectrum[0])]
    a1 = [1.1 - np.sin((i * 2 / shape_spectrum[0]) * 3.14 / 2) for i in range(int(shape_spectrum[0] / 2))]
    a2 = [0.1 + np.sin((i * 2 / shape_spectrum[0]) * 3.14 / 2) for i in range(int(shape_spectrum[0] / 2))]
    a = np.array(a2 + a1)

    spectrum_ft = fft(spectrum_one)
    spectrum_ift = ifft(spectrum_ft / a)
    simplest_fig(time, spectrum_one, spectrum_ift)
    pass
