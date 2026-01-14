import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from path_to_Yandex_Disk import path_to_YaDisk


def fig_plot(x, y):
    fig, ax = plt.subplots(1, figsize=(12, 6))

    n = file_name0.rfind('az-')
    az = file_name0[n+2:n+5]

    # Show the major grid lines with dark grey lines
    plt.grid(b=True, which='major', color='#666666', linestyle='-')

    # Show the minor grid lines with very faint and almost transparent grey lines
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)

    ax.set_xlabel('freq, MHz', fontsize=20)
    ax.set_ylabel('Spectral power, dBm', fontsize=20)
    # ax.set_yscale('log')
    ax.set_title('Interference', fontsize=24)

    ax.plot(x, y, color='green', marker='o', markerfacecolor='red', label='Azimuth = '+az)
    font = font_manager.FontProperties(family='Comic Sans MS',
                                       weight='bold',
                                       style='normal', size=20)
    ax.legend(prop=font)
    # ax.plot(_x, _y, 'ro-', label='Time = 0.2 sec') # Запись попроще, почти как в Матлаб
    # ax.plot()
    # ax.legend()

    plt.show()
    return fig


head_path = path_to_YaDisk()
file_name0 = head_path + '\\Measure\\Fast_Esquition\\26122019interference\\az-14-2_05sun'
file = file_name0 + '.txt'  # C:\Users\PC\YandexDisk\Piton_Progects\Fast_Esquition
f=open(file,"r")    # E:\\YandexDisk-svit-commerc\\Measure\\Fast_Esquition\\26122019interference\\az-14-2_05sun'
lines=f.readlines()
result=[]
for x in lines[26:]:
    result.append(x.split(';')[0:2])
f.close()
result1 = [float(s) for s1 in result for s in s1]
result2 = np.reshape(result1,(len(result1)//2,2))
fig = fig_plot(result2[:,0], result2[:,1])
fig.savefig(file_name0)

pass







