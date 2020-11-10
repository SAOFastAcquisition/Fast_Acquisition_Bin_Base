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