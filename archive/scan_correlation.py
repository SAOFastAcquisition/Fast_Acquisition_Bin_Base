import numpy as np
from archive.path_to_Yandex_Disk import path_to_YaDisk
import matplotlib.pyplot as plt

head_path = path_to_YaDisk()
file_name0 = head_path + '\\Measure\\Fast_Acquisition\\2020_03_18sun\\20200318-1353_-24-3'
scan_init = np.loadtxt(file_name0+'_scan'+'.txt')
size = np.shape(scan_init)


# ****** Блок исходных параметров для обработки *******
kf = 4  # Установка разрешения по частоте
kt = 1  # Установка разрешения по времени
  # Номер зоны Найквиста
# *****************************************************

delta_t = 8.1925e-3
delta_f = 7.8125
t_start_flame = 0
t_stop_flame = 270

n_start_flame = int(t_start_flame // (delta_t * kt))
n_stop_flame = int(t_stop_flame // (delta_t * kt))
x = scan_init[1, n_start_flame:n_stop_flame]
y = scan_init[3, n_start_flame:n_stop_flame]
z = scan_init[5, n_start_flame:n_stop_flame]
val = np.vstack((x, y, z))
#corr_func = np.corrcoef(scan_init[16:, n_start_flame:n_stop_flame], rowvar=True)
pass
delta_n = 10
corr_func = np.zeros((5, 5, int((n_stop_flame - n_start_flame - 60) // delta_n) + 1))
i = 0
delta_n = 10
while n_start_flame + i * delta_n + 60 <= n_stop_flame:
    corr_func[:, :, i] = np.corrcoef(scan_init[0:5, n_start_flame + i * delta_n:n_start_flame + i * delta_n + 60],
                                     rowvar=True)
    i += 1
pass
arg = [(t_start_flame + (int(i) * 10 + 30) * delta_t) for i in range((n_stop_flame - n_start_flame - 60) // delta_n + 1)]
fig, ax = plt.subplots(1, figsize=(12, 6))
for i1 in range(1, 5):
    ax.plot(arg, corr_func[0, int(i1), :])
plt.show()
# (t_start_flame + (int(i) * 10 + 30) * delta_t)