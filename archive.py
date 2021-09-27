def afc_interpol(freq, Kp0, freq_i):
    i0 = 0
    Kp_i = 0
    Kp = [0] * len(freq_i)
    for k in range(len(freq_i)):
        for i in range(i0, len(freq)):
            if freq_i[k] >= freq[i] and freq_i[k] < freq[i + 1]:
                # print(k, freq_i[k], freq_i[k+1])
                Kp_i = Kp0[i] + (Kp0[i + 1] - Kp0[i]) / (freq[i + 1] - freq[i]) * (freq_i[k] - freq[i])
                Kp[k] = Kp_i
                i0 = i

        continue
    return Kp


def afc_correction(freq):
    try:
        mat1 = scipy.io.loadmat('E:\\YandexDisk-svit-commerc\\Piton_Progects\\Fast_Esquition\\Pic_16_09_19\\Kp0.mat')
    except FileNotFoundError:
        try:
            mat1 = scipy.io.loadmat(
                'C:\\Users\\PC\\YandexDisk\\Piton_Progects\\Fast_Esquition\\Pic_16_09_19\\Kp0.mat')
        except FileNotFoundError as ffe:
            print(ffe)
            try:
                mat1 = scipy.io.loadmat(
                    'D:\\YandexDisk\\Piton_Progects\\Fast_Esquition\\Pic_16_09_19\\K_chl0.mat')

            except FileNotFoundError as ffe:
                print(ffe)

    try:
        mat2 = scipy.io.loadmat('E:\\YandexDisk-svit-commerc\\Piton_Progects\\Fast_Esquition\\Pic_16_09_19\\freq.mat')
    except FileNotFoundError:
        try:
            mat2 = scipy.io.loadmat(
                'D:\\YandexDisk\\Piton_Progects\\Fast_Esquition\\Pic_16_09_19\\freq.mat')
        except FileNotFoundError as ffe:
            print(ffe)

    Kp0 = [s for s1 in np.array(mat1['K_ch_l0']) for s in s1[0:1]]
    freq_ch = [s / 1000000 for s1 in np.array(mat2['freq']) for s in s1]

    # Kp = [0] * len(freq)
    # for k in range(len(freq)):
    Kp = afc_interpol(freq_ch, Kp0, freq)

    Kp_max = max(Kp)
    Kp = np.asarray(Kp)
    Kp = 10 ** ((Kp_max-Kp)/10)
    # Kp =  [10 **((s - Kp_max)/10 for s in Kp)]
    return Kp


def fig_plot(x, y):
    # fig = plt.subplot()
    fig, ax = plt.subplots(1, figsize=(12, 6))
    # fig.suptitle('Graffics', fontsize=18)



    # Show the major grid lines with dark grey lines
    plt.grid(b=True, which='major', color='#666666', linestyle='-')

    # Show the minor grid lines with very faint and almost transparent grey lines
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)

    ax.set_xlabel('freq, MHz', fontsize=20)
    ax.set_ylabel('S', fontsize=20)
    # ax.set_yscale('log')
    ax.set_title('Graffics', fontsize=24)

    ax.plot(x, y, color='green', marker='o', markerfacecolor='red', label='Time = 0.2 sec')
    # ax.plot(x, y, 'ro-', label='Time = 0.2 sec') # Запись попроще, почти как в Матлаб
    ax.plot()
    ax.legend()

    plt.show()


def spectr_construction(Spectr, kf, kt):
    ''' Функция формирует спектр принятого сигнала с требуемым разрешением по частоте и времени. Максимальное
    разрешение отсчетов по времени 8192 мкс и по частоте 7,8125 МГц. Путем суммирования и усреднерия по kt*kf
    отсчетам разрешение по частоте и по времени в исходном спектре Spectr уменьшается в kf и kt раз,
    соответственно. Преобразованный спектр возвращается как S1.
    '''

    N_col1 = N_col // kf
    N_row1 = N_row // kt
    S1 = np.zeros((N_row1, N_col1))

    for i in range(N_row1):
        for j in range(N_col1):
            try:
                S1[i, j] = np.sum(Spectr[i * kt: (i + 1) * kt, j * kf: (j + 1) * kf])
                N_mesh = (Spectr[i * kt: (i + 1) * kt, j * kf: (j + 1) * kf] > 10).sum()
                if N_mesh == 0:
                    S1[i, j] = 2
                else:
                    S1[i, j] = S1[i, j] / N_mesh
                if S1[i, j] == 0:
                    S1[i, j] = 2
                # if (j > 3) & (S1[i, j] > 1.5 * np.sum(S1[i, j-3:j])//3):
                #     S1[i, j] = np.sum(S1[i, j-3:j])//3
                if robust_filter == 'y':
                    a = param_robust_filter
                    if (i > 3) & (S1[i, j] < 1 / a * np.sum(S1[i - 3:i - 1, j]) // 2):
                        S1[i, j] = np.sum(S1[i - 1, j])
                    if (i > 3) & (S1[i, j] > a * np.sum(S1[i - 3:i - 1, j]) // 2):
                        # print(S1[i - 3:i+1, j])
                        S1[i, j] = np.sum(S1[i - 1, j])
                        # print(S1[i, j])
                        pass

            except IndexError as allert_message:
                print(allert_message, 'ind i = ', i, 'ind j = ', j)
                pass
            except ValueError as value_message:
                print(value_message, 'ind i = ', i, 'ind j = ', j)
                pass

    return S1  # // kt // kf


def calibration(t_cal, spectrum):
    size_spectrum = spectrum.shape
    level_matrix = np.ones((size_spectrum[0], 2))
    cal_level = np.zeros(size_spectrum[0])
    n_cal = np.zeros(4)
    for i in range(4):
        n_cal[i] = int(t_cal[i] // (delta_t * kt))
    for sp_row in range(size_spectrum[0]):
        level_matrix[sp_row, 0] = np.average(spectrum[sp_row, int(n_cal[0]):int(n_cal[1])])
        level_matrix[sp_row, 1] = np.average(spectrum[sp_row, int(n_cal[2]):int(n_cal[3])])
        cal_level[sp_row] = np.abs(level_matrix[sp_row, 0] - level_matrix[sp_row, 1])
    max_cal_level = np.max(cal_level)
    cal_level = cal_level / max_cal_level

    for sp_row in range(size_spectrum[0]):
        spectrum[sp_row, :] = (spectrum[sp_row, :] - np.min(level_matrix[sp_row, :])) / cal_level[sp_row]
    return spectrum


def self_calibration1():
    """ Принимает среднее значение спектра за первые 5 - 10 сек наблюдения, пока не идет Солнце, при максимальном
    разрешении по частоте kf = 1

    :return:
    """
    file_calibr = 'self_calibr.txt'  # Файл, в который записываются спектры наблюдений в отсутствие Солнца
    if kt <= 500 and kf > 1:
        return

    ind_observation_id = file_name0.rfind('\\') + 1  # Индекс, с которого вычленяется идентификатор записи
    # наблюдения из ее полного пути
    observation_id = file_name0[ind_observation_id:]  # Идентификатор записи наблюдения, использованной
    # для учета АЧХ тракта
    file_observation_names = file_name0[0:ind_observation_id - 1] + '\\observ_name.txt'  # Файл, в который
    # записываются идентификаторы записей наблюдений, используемых для учета АЧХ тракта
    file_observ_name = open(file_observation_names, 'a+')
    file_observ_name.seek(0)
    observation_list = file_observ_name.read()
    if not observation_list.count(observation_id):
        # if not os.path.isfile(file_observation_names):
        if not os.path.isfile(file_calibr):
            np.savetxt(file_calibr, spectr_freq[0, :])
        else:
            self_calibr = np.loadtxt(file_calibr)
            try:
                self_calibr1 = np.vstack((self_calibr, spectr_freq[0, :]))
            except ValueError:
                print('Function self_calibration not append new data')
                file_observ_name.close()
                return
            np.savetxt(file_calibr, self_calibr1)
            file_observ_name.write(observation_id + '\n')
    file_observ_name.close()
    return


# 173, 173.6, 173.8, 174.38
# t_cal = [0, 14, 17, 31]         # Для скана "20200318-1353_-24-3"
# t_cal = [0, 13, 17, 35]


def image_filter(img_src):
    import cv2
    # read image
    # img_src = cv2.imread('sample.jpg')
    # prepare the 5x5 shaped filter
    # kernel = np.array([[1, 1, 1, 1, 1],
    #                    [1, 1, 1, 1, 1],
    #                    [1, 1, 1, 1, 1],
    #                    [1, 1, 1, 1, 1],
    #                    [1, 1, 1, 1, 1]])
    kernel = np.array([[1, 1, 1, 1],
                       [1, 1, 1, 1],
                       [1, 1, 1, 1],
                       [1, 1, 1, 1]])
    kernel = kernel / sum(kernel)  # filter the source image
    img_rst = cv2.filter2D(img_src, -1, kernel)
    # save result image
    # cv2.imwrite('result.jpg', img_rst)
    return img_rst
    # Источник: https://tonais.ru/library/filtratsiya-izobrazheniy-s-ispolzovaniem-svertki-opencv-v-python


def low_freq_filter(x, h):
    """ФНЧ с импульсной характеристикой h. Входная последовательность - x"""

    n_input = len(x)
    n_filter = len(h)
    # n - Длина входной последовательности
    # h - отклик фильтра НЧ

    y = [0] * n_input  # Выходная последовательность после интерполяции и НЧ фильтрации
    for n in range(0, n_input):
        for k in range(0, n_filter):
            try:
                y[n] += x[n - k - 1] * h[k]
            except IndexError:
                y[n] += 0
                print('ind n = ', n, 'ind k = ', k)
                pass
    return y


def filter_coeff(length_fft, filters_order, band_pass):
    """ Отдает укороченную импульсную характеристику h_short согласно заданному порядку КИХ фильтра"""
    import filters as ftr
    #   length_fft - длина БПФ
    #   filters_order - Порядок фильтра
    #   band_pass - полоса фильтра в отсчетах
    h, fft_h = ftr.synt_filter(length_fft, filters_order, band_pass)  # Отклик прямоугольного модельного фильтра
    n_filter = len(h)
    h_short = h[n_filter // 2 - filters_order // 2:n_filter // 2 + filters_order // 2]  # Выбор количества значений
    # отклика, соответствующего порядку фильтра
    # w_inter = [0.54 - 0.46 * np.cos(2 * np.pi * i / (n - 1)) for i in range(0, n * l_interpol)] # Окно Хэминга
    # w_inter = flattop(n)  # from scipy максимально плоское окно (применяется при полифазной обработке)
    return h_short


def pic_title():
    title0 = file_name0[-19:-2]
    title1 = '  ' + title0[0:4] + '.' + title0[4:6] + '.' + title0[6:8] + \
             ' time=' + title0[9:11] + ':' + title0[11:13] + ' azimuth=' + title0[14:17]
    if not file_name0.find('sun') == -1:
        title2 = 'Sun intensity'
    elif not file_name0.find('crab') == -1:
        title2 = 'Crab intensity'
    elif not file_name0.find('calibr') == -1:
        title2 = 'Calibration'
        title0 = file_name0[-23:-2]
        title1 = '  ' + title0[0:4] + '.' + title0[4:6] + '.' + title0[6:8] + \
                 ' chanell att=' + title0[14:17] + ' source att=' + title0[18:21]
    elif not file_name0.find('test') == -1:
        title0 = file_name0[-24:-2]
        title2 = 'Test interference'
        title1 = '  ' + title0[0:4] + '.' + title0[4:6] + '.' + title0[6:8] + \
                 ' chanell att=' + title0[15:18] + ' source att=' + title0[19:22]
        pass
    else:
        title2 = []
    return title1, title2

