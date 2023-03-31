import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

path = r'H:\Fast_Acquisition\2023\Interference\2023_03_18_01_200_3500MHz.csv'
a = pd.read_csv(path)

a.plot(x='TraceA(Hz)', y='TraceA(dBm)')
plt.grid()
plt.show()
pass