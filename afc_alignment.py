import os
import numpy as np
from path_to_Yandex_Disk import path_to_YaDisk


def noise_kp(file_name, diff='n'):
    """ Функция принимает файл шумового измерения мощностной АЧХ и возвращает
        нормированный к максимуму коэффициент передачи тракта по мощности"""
    spectr_noise_out = np.loadtxt(file_name + '.txt')
    n_row, n_col = np.shape(spectr_noise_out)
    s = np.zeros(n_col)
    s1 = np.zeros(n_col)
    # Усреднение по времени
    for i in range(n_col):
        if diff == 'n':
            if int(file_name[-1]) == 2:
                s[i] = np.sum(spectr_noise_out[1600:, i]) / (n_row - 1600)
            else:
                s[i] = np.sum(spectr_noise_out[:1600, i]) / 1600
        else:
            s[i] = np.sum(spectr_noise_out[:1600, i]) / 1600
            s1[i] = np.sum(spectr_noise_out[2000:3600, i]) / 1600
            s[i] -= s1[i]

    s_max = np.max(s)
    kp_norm = s / s_max
    return kp_norm, n_col


def align_func(calibr_file_name, diff='n', aver_param=2):
    """ Функция возвращает коэффициенты, выравнивающие исходную АЧХ

    """
    # Исходные данные
    # N_Nyq = 3
    delta_f = 7.8125
    # aver_param = 2

    # Определение пути к файлу, где лежат результаты измерения АЧХ
    # по мощности с генератором шума на входе тракта, чтение параметра
    # усреднения n_aver_noise для этого измерения
    head_path = path_to_YaDisk()
    file_name0 = head_path + '\\Measure\\Fast_Acquisition\\Calibration\\' + calibr_file_name
    n_Nyq = int(calibr_file_name[-1])
    f_in1 = open(file_name0 + '.txt')
    n_aver_noise = int((f_in1.readline())[2])
    aver_param_noise = 2 ** (6 - n_aver_noise)
    f_in1.close()

    kp_norm, n_col = noise_kp(file_name0, diff)

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
