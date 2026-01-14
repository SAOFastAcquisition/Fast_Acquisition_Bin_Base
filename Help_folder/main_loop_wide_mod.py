import numpy as np
import matplotlib.pyplot as plt
import pickle
import gzip
import os
from scipy.optimize import curve_fit
from pathlib import Path
from Supporting_func import path_to_data
from Polyphase.cic_filter import signal_filtering


def save_list(filename, data, method='json'):
    """Сохраняет список в файл разными методами"""
    if method == 'json':
        import json
        with open(filename, 'w', encoding='utf-8') as file:
            json.dump(data, file, indent=4)
    elif method == 'pickle':
        import pickle
        with open(filename, 'wb') as file:
            pickle.dump(data, file)
    elif method == 'text':
        with open(filename, 'w', encoding='utf-8') as file:
            for item in data:
                file.write(str(item) + '\n')


# ------------------------------------------------------------
# 1. Гауссова функция для аппроксимации
# ------------------------------------------------------------
def gauss(_x, _A, _x0, _sigma):
    return _A * np.exp(-(_x - _x0) ** 2 / (2 * _sigma ** 2))


# ------------------------------------------------------------
# 2. Функция для аппроксимации одного столбца (одной ДН)
# ------------------------------------------------------------
def fit_main_lobe_gauss(_y, _x=None, _threshold=0.1, _max_half_width=None):
    """
    _y: 1D array, нормированная мощность в разах
    _x: 1D array, индексы отсчётов (если None, то [0, 1, 2, ...])
    _threshold: порог по мощности для обрезания лепестка
    _max_half_width: максимальная полуширина области для фита (в отсчётах)

    Возвращает:
        A, x0, sigma (в отсчётах), y_fit, x_fit (область фита)
    """
    if _x is None:
        _x = np.arange(len(_y))

    # Находим максимум
    peak_idx = np.argmax(_y)
    peak_val = _y[peak_idx]

    # Определяем границы лепестка по порогу
    left = peak_idx
    right = peak_idx
    while left > 0 and _y[left] > _threshold:
        left -= 1
    while right < len(_y) - 1 and _y[right] > _threshold:
        right += 1

    # Если задана _max_half_width, сужаем область
    if _max_half_width is not None:
        left = max(left, peak_idx - _max_half_width)
        right = min(right, peak_idx + _max_half_width)

    # Область для фита
    fit_idx = np.arange(left, right + 1, dtype=int)
    x_fit = _x[fit_idx]
    y_fit_data = _y[fit_idx]

    # Начальные приближения
    A0 = peak_val
    x00 = _x[peak_idx]
    # Грубая оценка sigma: половина ширины по уровню exp(-0.5)≈0.606 от максимума
    sigma0 = (right - left) / 4.0  # эвристика

    # Аппроксимация
    try:
        popt, pcov = curve_fit(gauss, x_fit, y_fit_data, p0=[A0, x00, sigma0])
        A, x0, sigma = popt
        # Гаусс на всей области
        y_fit_full = gauss(_x, A, x0, sigma)
    except Exception as e:
        print(f"Ошибка curve_fit: {e}")
        A, x0, sigma = A0, x00, sigma0
        y_fit_full = gauss(_x, A, x0, sigma)

    return A, x0, sigma, y_fit_full, x_fit, y_fit_data


def get_initial_data(_converted_data_file_path, _current_primary_file):
    """ Функция в зависимости от вида данных (полная полоса 1-3 ГГц, половинная полоса 1-2 или 2-3 ГГц,
    с двумя поляризациями или одной) выдает данные для построения графиков"""

    # Для полосы 1-3 ГГц и двух возможных поляризаций выдает по два спектра (1-2 и 2-3 ГГц) для каждой поляризации.
    # Если поляризация не задействована, то соответствующие спектры - пустые. Спектр 1-2 ГГц - в обратном порядке
    _path1 = Path(_converted_data_file_path, _current_primary_file + '_spectrum.npy')
    # spectrum = np.load(_path1, allow_pickle=True)
    if os.path.exists(f'{str(_path1)}.gz'):
        _filename_out = f'{str(_path1)}.gz'
        with gzip.open(_filename_out, "rb") as _fin:
            spectrum = np.load(_fin, allow_pickle=True)
    else:
        _spectrum = np.load(_path1, allow_pickle=True)
    with open(Path(_converted_data_file_path, _current_primary_file + '_head.bin'), 'rb') as inp:
        _head = pickle.load(inp)
    _n_aver = _head['n_aver']
    _band_size = _head['band_size']
    _polar = _head['polar']
    # n_aver = head['n_aver']

    # Разделяем составляющие  записи в полной полосе и с возможными двумя поляризациями

    _spectrum_left1 = _spectrum[0]
    _spectrum_left2 = _spectrum[1]
    _spectrum_right1 = _spectrum[2]
    _spectrum_right2 = _spectrum[3]

    pass
    return _spectrum_left1, _spectrum_left2, _spectrum_right1, _spectrum_right2, int(_n_aver), _band_size, _polar


def data_preparing(_current_primary_dir, _current_primary_file):


    current_data_dir = '2025/Test_and_calibration'
    primary_data_dir = 'Primary_data'  # Каталог исходных данных (за определенный период, здесь - год)
    converted_data_dir = 'Converted_data'  # Каталог для записи результатов конвертации данных и заголовков
    data_treatment_dir = 'Data_treatment'  # Каталог для записи результатов обработки, рисунков



    current_converted_dir = _current_primary_dir + '_conv'
    current_converted_path = Path(converted_data_dir, current_converted_dir)
    current_treatment_dir = _current_primary_dir + '_treat'
    current_treatment_path = Path(data_treatment_dir, current_treatment_dir)

    converted_data_file_path, head_path = path_to_data(current_data_dir, current_converted_path)
    data_treatment_file_path, head_path = path_to_data(current_data_dir, current_treatment_path)
    spectrum_extr_left1, spectrum_extr_left2, spectrum_extr_right1, spectrum_extr_right2, n_aver, band_size, polar = \
        get_initial_data(converted_data_file_path, _current_primary_file)

    nc = int(44 / delta_t)
    nl = int(30 / delta_t)
    nr = int(58 / delta_t)

    if polar == 'left':
        spectrum1 = spectrum_extr_left1[nl:nr, :]
        spectrum2 = spectrum_extr_left2[nl:nr, :]
    else:
        spectrum1 = spectrum_extr_right1[nl:nr, :]
        spectrum2 = spectrum_extr_right2[nl:nr, :]

    spectrum1 = signal_filtering(spectrum1.T).T
    spectrum2 = signal_filtering(spectrum2.T).T

    spectrum1_max = np.nanmax(spectrum1, axis=0)
    spectrum2_max = np.nanmax(spectrum2, axis=0)
    indicis1 = np.where(spectrum1_max < 1e9)
    indicis2 = np.where(spectrum2_max < 3e8)
    spectrum1[:, indicis1] = 0
    spectrum2[:, indicis2] = 0

    spectrum1 = spectrum1 / spectrum1_max
    spectrum2 = spectrum2 / spectrum2_max

    return spectrum1, spectrum2


if __name__ == '__main__':

    current_primary_dir = '2025_12_21test'
    current_primary_file = '2025-12-21_01'

    delta_t = 8.3886e-3
    delta_f = 7.8125 / 2

# ------------------------------------------------------------
# 4. Аппроксимация всех столбцов
# ------------------------------------------------------------

    sp1, sp2 = data_preparing(current_primary_dir, current_primary_file)
    dat = (sp1, sp2)

    sp1_d = sp1
    row1, col1 = np.shape(sp1_d)

    x_global = np.arange(row1)

    freq1 = [2000 - delta_f * (i + 0.5) for i in range(col1)]
    freq2 = [2000 + delta_f * (i + 0.5) for i in range(col1)]
    # plt.figure(figsize=(10, 6))

    for j in range(2):
        params = []  # список параметров (A, x0, sigma) для каждой частоты
        obj = dat[j]
        for i in range(col1):
            y = obj[:, i]
            A, x0, sigma, y_fit, x_fit, y_fit_data = fit_main_lobe_gauss(
                y, _x=x_global, _threshold=0.1, _max_half_width=1500
            )
            b = gauss(x_fit, A, x0, sigma)
            c = np.std(b - y_fit_data)
            if c > 0.05:
                sigma = np.nan
            if sigma < 100:
                sigma = np.nan
            params.append((A, x0, sigma, c))

            # Графики
        #     if c < 0.02:
        #         plt.plot(x_global, _y, 'o', markersize=3, label=f'Данные freq{i+1}' if i == 0 else "")
        #         plt.plot(x_global, y_fit, '-', linewidth=2, label=f'Гаусс фит freq{i + 1}')
        #         plt.plot(x_fit, y_fit_data, 's', markersize=4, alpha=0.7)  # область фита
        #
        # plt.xlabel('Отсчёты')
        # plt.ylabel('Нормированная мощность')
        # plt.title('Аппроксимация главного лепестка гауссовой кривой')
        # plt.legend()
        # plt.grid(alpha=0.3)
        # plt.show()

        # ------------------------------------------------------------
        # 5. Вывод параметров
        # ------------------------------------------------------------
        print("Параметры гауссовой аппроксимации:")
        print("Частота \t A \t x0 (отсч.) \t sigma (отсч.)")
        for i, (A, x0, sigma, c) in enumerate(params):
            print(f"{i + 1} \t\t {A:.4f} \t {x0:.2f} \t\t {sigma:.2f} \t\t {c:.4f}")
        if j == 0:
            params1 = params
        else:
            params2 = params

    sigma_p1 = [a[2] for a in params1]
    sigma_p2 = [a[2] for a in params2]
    freq1.reverse()
    sigma_p1.reverse()

    sigma_p = np.array(sigma_p1 + sigma_p2)
    freq = np.array(freq1 + freq2)
    ind = np.isnan(sigma_p)
    sigma_p_mod = sigma_p[~ind]
    freq_mod = freq[~ind]

    # plt.plot(freq1, sigma_p1)
    # plt.plot(freq2, sigma_p2)
    # plt.show()
    coefficients = np.polyfit(np.log(freq_mod), np.log(sigma_p_mod), 1)
    appr_sigma = np.exp(coefficients[1]) * freq_mod ** coefficients[0]
    appr_sigma0 = np.exp(coefficients[1]) * 1000 ** coefficients[0] * 1000 / freq_mod

    hpbw_h = 200 / 1000**coefficients[0] * sigma_p_mod / np.exp(coefficients[1])

    delta_sigma = sigma_p_mod - appr_sigma

    fig = plt.figure(figsize=(7, 4))
    ax = fig.add_subplot()
    ax.plot(freq_mod, sigma_p_mod, 'o', markersize=2, label='расчетные значения сигма')
    ax.plot(freq_mod, hpbw_h, 'o', markersize=2, label='расчетные значения HPBW horizont')
    ax.plot(freq_mod, appr_sigma, '-', linewidth=2, label='аппроксимация степенной функцией')
    # ax.plot(freq_mod, appr_sigma0, '-', linewidth=2)
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.grid(alpha=0.3)
    plt.xlabel('Частота, МГц')
    plt.ylabel('Сигма, отсчеты')
    plt.title('Параметр сигма главного лепестка гауссовой кривой')
    ax.legend()
    ax.grid(True, which='major', linestyle='-', linewidth=0.9, color='black')
    ax.grid(True, which='minor', linestyle=':', linewidth=0.9, color='red')
    plt.show()

    plt.plot(freq_mod, delta_sigma)
    # plt.plot(freq2, sigma_p2)
    plt.show()
    sigma_file_path = Path(r"I:\Fast_Acquisition\2025\Test_and_calibration\Converted_data\2025_12_21test_conv",
                           f'{current_primary_file}_hpbw.bin')
    freq_file_path = Path(r"I:\Fast_Acquisition\2025\Test_and_calibration\Converted_data\2025_12_21test_conv",
                           f'{current_primary_file}_freq.bin')

    save_list(sigma_file_path, hpbw_h, method='pickle')
    save_list(freq_file_path, freq_mod, method='pickle')