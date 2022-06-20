import numpy as np
import pandas as pd
import matplotlib.pylab as plt

import datashader as ds
import datashader.transfer_functions as tf
from datashader.utils import export_image

from functools import partial

background = "white"
export = partial(export_image, background=background, export_path="../export")

N = 100000
df = pd.DataFrame(np.random.random((N, 3)), columns=['x', 'y', 'z'])

f, ax = plt.subplots(2, 2)
ax_r = ax.ravel()

ax_r[0].scatter(df['x'], df['y'], df['z'].mean(), cmap=plt.get_cmap('jet'))
# ax_r[1].hist(df['x'])
# ax_r[2].hist(df['y'])
# ax_r[3].plot(df['z'])

cvs = ds.Canvas(plot_width=250, plot_height=300)
agg = cvs.points(df, 'x', 'y', ds.mean('z'))
# a = export(tf.shade(agg, cmap=['blue', 'red'], how='eq_hist'), 'test')
# cmap = plt.get_cmap('jet')
a = export(tf.shade(agg, cmap=plt.get_cmap('jet'), how='eq_hist'), 'test')
plt.show()
