import numpy as np
import matplotlib.pyplot as plt


def noise_kp(spectrum_frame, diff='n'):
    """ Функция принимает файл шумового измерения мощностной АЧХ и возвращает
        нормированный к максимуму коэффициент передачи тракта по мощности"""

    n_row, n_col = np.shape(spectrum_frame)
    s = np.zeros(n_col)
    s1 = np.zeros(n_col)
    # Усреднение по времени
    for i in range(n_col):
        if diff == 'n':
            s[i] = np.sum(spectrum_frame[:, i]) / n_row
        else:
            s[i] = np.sum(spectrum_frame[:, i]) / n_row
            s1[i] = np.sum(spectrum_frame[:, i]) / n_row
            s[i] -= s1[i]
    aver_level = np.sum(s[:])
    s /= aver_level
    fig, ax = plt.subplots(1, figsize=(12, 6))
    ax.plot(i in range(n_col), s)

    s_max = np.max(s)
    kp_norm = s / s_max
    return kp_norm, n_col


def align_func1(spectrum_frame, diff='n', aver_param=2):
    """ Функция возвращает коэффициенты, выравнивающие исходную АЧХ

    """
    # Исходные данные
    # N_Nyq = 3
    delta_f = 7.8125
    # aver_param = 2

    # Определение пути к файлу, где лежат результаты измерения АЧХ
    # по мощности с генератором шума на входе тракта, чтение параметра
    # усреднения n_aver_noise для этого измерения

    n_Nyq = 3
    # f_in1 = open(file_name0 + '.txt')
    # n_aver_noise = int((f_in1.readline())[2])
    n_aver_noise = 4
    aver_param_noise = 2 ** (6 - n_aver_noise)
    # f_in1.close()

    kp_norm, n_col = noise_kp(spectrum_frame, aver_param, diff)

    # Исключение корректировки коэффициента усиления в зоне действия режекторных фильтров
    if n_Nyq == 3:
        n1 = int((80 - delta_f / 2 / aver_param_noise) // (delta_f / aver_param_noise))
        n2 = int((230 - delta_f / 2 / aver_param_noise) // (delta_f / aver_param_noise))
    else:
        n1 = int((810 - delta_f / 2 / aver_param_noise) // (delta_f / aver_param_noise))
        n2 = int((890 - delta_f / 2 / aver_param_noise) // (delta_f / aver_param_noise))
    kp_norm[n1:n2] = 1

    band = int(aver_param_noise / aver_param)
    kp_band = [np.sum(kp_norm[band * i:band * (i + 1)]) / band for i in range(int(n_col / band))]

    align_coeff = [1 / a for a in kp_band]
    return align_coeff


if __name__ == '__main__':
    align_coeff = align_func(2)
    pass
