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


def gauss(_x, _s, _a=1, _x0=0):
    return _a*np.exp(-((_x-_x0)/_s)**2)


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
    spectrum_one = spectrum[:, 492]
    shape_spectrum = np.shape(spectrum)

    n_freq = shape_spectrum[1]
    delta_f = 2000/n_freq   # unit = 'MHz'
    time = [i * delta_t for i in range(shape_spectrum[0])]

    center = 215                                    # sec time
    sun_wide = 190                                  # sec time
    center_arc = 0
    sun_wide_arc = 1920                                             # arcsec
    time_to_angle_coeff = sun_wide_arc / sun_wide                   # arcsec/sec
    angle_per_sample = time_to_angle_coeff * delta_t / 3600 / 57.2  # rad
    n_angle = int(6.28 / angle_per_sample)
    n_angle_center = int(n_angle // 2)
    delta_angle = delta_t * time_to_angle_coeff
    n_time_center = int(center / delta_t - 1)
    n_wide = int(150 / delta_t)
    # sun_centered = spectrum_one[n_time_center - n_wide - 1:n_time_center + n_wide]
    sun_centered = [0] * n_angle
    #           *** Совмещение точки кульминации с центральным отсчетом ***
    #                размещение скана Солнца посередине зоны Найквиста
    sun_centered[n_angle_center - n_wide - 1:n_angle_center + n_wide] = \
        spectrum_one[n_time_center - n_wide - 1:n_time_center + n_wide]
    #           *** Заполнение "хвостов" зоны Найквиста со сканом  ***
    for i in range(n_angle_center - n_wide):
        sun_centered[i] = sun_centered[n_angle_center - n_wide + 1]
        sun_centered[n_angle_center + n_wide + i] = sun_centered[n_angle_center + n_wide - 1]

    angle = np.array([t * angle_per_sample for t in range(n_angle)])

    main_lobe1 = gauss(angle, 60 / 3600 / 57.2, 300, n_angle_center * angle_per_sample)
    main_lobe2 = gauss(angle, 20 / 3600 / 57.2, 300, n_angle_center * angle_per_sample)
    for i in range(len(main_lobe1)):
        if main_lobe1[i] < 1e-4:
            main_lobe1[i] = 1e-4
    for i in range(len(main_lobe2)):
        if main_lobe2[i] < 1e-4:
            main_lobe2[i] = 1e-4
    # a = np.array(a2 + a1)
    a_1 = fft(main_lobe1)   # Передаточная функция 1 телескопа (имеется)
    a_2 = fft(main_lobe2)   # Передаточная функция 2 телескопа (желательная)
    for i in range(len(main_lobe1)):
        if abs(a_1[i]) < 1e-4:
            a_1[i] = 1e-4
    for i in range(len(main_lobe2)):
        if abs(a_2[i]) < 1e-4:
            a_2[i] = 1e-4
    r = a_2 / a_1
    for i in range(n_angle):
        if abs(r[i]) > 10:
            r[i:-i] = 0
            r_inv = [0] * i
            r_inv = r[i:0:-1]
            r[n_angle - i - 1:-1] = r_inv
            break

    main_lobe_ift2 = ifft(a_2)
    spectrum_ft = fft(sun_centered)
    spectrum_ift = ifft(spectrum_ft * r)
    simplest_fig(angle, abs(spectrum_ift), sun_centered)
    pass
