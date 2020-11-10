"""
=============================================
Generate polygons to fill under 3D line graph
=============================================

Demonstrate how to create polygons which fill the space under a line
graph. In this example polygons are semi-transparent, creating a sort
of 'jagged stained glass' effect.
"""

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import numpy as np


def poly_graph(xs, y, zs):
    fig = plt.figure()
    ax = fig.gca(projection='3d')


    def cc(arg):
        return mcolors.to_rgba(arg, alpha=0.6)

    k, l = np.shape(y)
    verts = []
    for i in range(k):
        ys = y[i, :]
        verts.append(list(zip(xs, ys)))

    poly = PolyCollection(verts, facecolors=['r', 'g', 'c', 'y'])
    poly.set_alpha(0.7)
    ax.add_collection3d(poly, zs=zs, zdir='y')

    ax.set_xlabel('X')
    ax.set_xlim3d(132, 134)
    ax.set_ylabel('Y')
    ax.set_ylim3d(2000, 3000)
    ax.set_zlabel('Z')
    ax.set_zlim3d(0, 2)

    plt.show()


# xs = np.arange(0, 10, 0.4)
# y = [1]
# zs = [0.0, 1.0, 2.0, 3.0]
# i = 0
# for z in zs:
#     ys = np.random.rand(len(xs))
#     ys[0], ys[-1] = 0, 0
#     y[i, :] = ys

# poly_graph(xs, y, zs)
