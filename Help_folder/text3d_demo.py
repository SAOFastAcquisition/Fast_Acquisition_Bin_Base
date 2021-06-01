'''
======================
Text annotations in 3D
======================

Demonstrates the placement of text annotations on a 3D plot.

Functionality shown:
- Using the text function with three types of 'zdir' values: None,
  an axis name (ex. 'x'), or a direction tuple (ex. (1, 1, 0)).
- Using the text function with the color keyword.
- Using the text2D function to place text on a fixed position on the ax object.
'''

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import pickle
import pandas as pd

fig = plt.figure()
ax = fig.gca(projection='3d')

# Demo 1: zdir
zdirs = (None, 'x', 'y', 'z', (1, 1, 0), (1, 1, 1))
xs = (1, 4, 4, 9, 4, 1)
ys = (2, 5, 8, 10, 1, 2)
zs = (10, 3, 8, 9, 1, 8)

for zdir, x, y, z in zip(zdirs, xs, ys, zs):
    label = '(%d, %d, %d), dir=%s' % (x, y, z, zdir)
    ax.text(x, y, z, label, zdir)

# Demo 2: color
ax.text(9, 0, 0, "red", color='red')

# Demo 3: text2D
# Placement 0, 0 would be the bottom left, 1, 1 would be the top right.
ax.text2D(0.05, 0.95, "2D Text", transform=ax.transAxes)

# Tweaking display region and labels
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.set_zlim(0, 10)
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')

# plt.show()
att = [i * 0.5 for i in range(64)]
s = {s: 10 ** (s / 10) for s in att}

for f in att:
    f = 2
head = {'n_aver': 2,
        'kurtosis': 20,
        'att1': 3,
        'att2': 3,
        'att3': 6,
        'band': 'half1',  # 'half1' - 1-2 GHz, 'half2' - 2-3 GHz, 'whole' - 1-3 GHz
        'polarization': 'both',  # 'left', 'right', 'both'
        'align_file_path': r'E:\Measure_res\2021_03_28sun\2021-03-28_02+24_.bin',
        'align_coeff_pos': 5}

with open(head['align_file_path'], 'wb') as out:
    pickle.dump(head, out)

with open(head['align_file_path'], 'rb') as inp:
    a_in = pickle.load(inp)
print(a_in)

columns_names = ['Date', 'att1', 'att2', 'att3',
                 'spectrum_left1', 'spectrum_left2', 'spectrum_right1',
                 'spectrum_right2',
                 'max_left1', 'max_left2', 'max_right1', 'max_right2']
df = pd.DataFrame(columns=columns_names)
df['Date'] = '02-02-02'
# calibrate_row = {'Date': file_name[1:12], 'att1': head['att1'], 'att2': head['att2'], 'att3': head['att3'],
#                  'spectrum_left1': align_coeff[0], 'spectrum_left2': align_coeff[1], 'spectrum_right1': [],
#                  'spectrum_right2': [],
#                  'max_left1': [], 'max_left2': [], 'max_right1': [], 'max_right2': []}

sf = pd.Series(index=columns_names)
sf['att1'] = head['att1']
sf['att2'] = head['att2']
sf['att3'] = head['att3']
df = df.append(sf, ignore_index=True)

file_name0 = r'F:\Fast_Acquisition\2021\Results\2021_04_15test\2021-04-15_05'
data_ind = file_name0.rindex(r'\2021')
data = file_name0[data_ind+1:data_ind+11]

pass
