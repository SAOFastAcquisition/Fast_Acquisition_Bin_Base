import scipy.io
import numpy as np
import matplotlib.pyplot as plt


def afc_int(freq, Kp0, freq_i):
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

N_Nyq = 2
kf = 1
N_col = 128

mat1 = scipy.io.loadmat('C:\\Users\\PC\\YandexDisk\\Piton_Progects\\Fast_Esquition\\Pic_16_09_19\\Kp0.mat')
mat2 = scipy.io.loadmat('C:\\Users\\PC\\YandexDisk\\Piton_Progects\\Fast_Esquition\\Pic_16_09_19\\freq.mat')
print(mat1)
# a = mat('Kp0')
# print(a)

Kp0 = [s for s1 in np.array(mat1['Kp0']) for s in s1[0:1]]
freq = [s / 1000000 for s1 in np.array(mat2['freq']) for s in s1]

freq_esq = np.linspace(1000 * (N_Nyq - 1) + 3.9063 * kf, 1000 * N_Nyq - 3.9063 * kf, N_col // kf)


Kp = afc_int(freq, Kp0, freq_esq)

Kp_max = max(Kp)
Kp = [s - Kp_max for s in Kp]

plt.plot(freq_esq, Kp)
plt.show()

