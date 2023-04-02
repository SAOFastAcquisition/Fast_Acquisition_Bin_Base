import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

path = r'J:\Fast_Acquisition\2023\Interference\2023_03_18_02_200_3500MHz.csv'
a = pd.read_csv(path)
a['TraceA(Hz)'] = a['TraceA(Hz)'] / 1000000

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.size'] = '10'

ax = a.plot(x='TraceA(Hz)', y='TraceA(dBm)', color='k', linewidth=1, figsize=(6.4, 3.8),
            xlabel='Frequency, MHz', ylabel='Relative Power, dB', legend='')
ax.tick_params(axis='both',  # Применяем параметры к обеим осям
               which='major',  # Применяем параметры к основным делениям
               direction='in',  # Рисуем деления внутри и снаружи графика
               length=5,  # Длинна делений
               width=1,  # Ширина делений
               color='black',  # Цвет делений
               pad=2,  # Расстояние между черточкой и ее подписью
               # labelsize=f_size,  # Размер подписи
               labelcolor='black',  # Цвет подписи
               bottom=True,  # Рисуем метки снизу
               top=True,  # сверху
               left=True,  # слева
               right=True,  # и справа
               labelbottom=True,  # Рисуем подписи снизу
               labeltop=False,  # сверху
               labelleft=True,  # слева
               labelright=False,  # и справа
               labelrotation=0)  # Поворот подписей
plt.grid()
plt.tight_layout(pad=0.2, w_pad=0.2, h_pad=1.0)
plt.show()
pass