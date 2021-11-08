"""
=============================================
Generate polygons to fill under 3D line graph
=============================================

Demonstrate how to create polygons which fill the space under a line
graph. In this example polygons are semi-transparent, creating a sort
of 'jagged stained glass' effect.
"""
import os
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import numpy as np
from tkinter import *
from tkinter import messagebox as mb


def save_fig(func):
    def wrapper(*args, **kwargs):
        figure: 'matplotlib.figure.Figure' = func(*args, **kwargs)
        figure.savefig('Poly3D.png')
        print('Poly3D is saved')
        del figure
        if save_question() == 'no':
            if os.path.isfile('Poly3D.png'):
                os.remove('Poly3D.png')
                print("success")
            else:
                print('File not found')
    return wrapper


def cc(arg):
    return mcolors.to_rgba(arg, alpha=0.6)


@save_fig
def poly_graph3d():
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    xs = np.arange(0, 10, 0.4)
    verts = []
    zs = [0.0, 1.0, 2.0, 3.0]
    for z in zs:
        ys = np.random.rand(len(xs))
        ys[0], ys[-1] = 0, 0
        verts.append(list(zip(xs, ys)))

    poly = PolyCollection(verts, facecolors=[cc('r'), cc('g'), cc('b'),
                                             cc('y')])
    poly.set_alpha(0.7)
    ax.add_collection3d(poly, zs=zs, zdir='y')

    ax.set_xlabel('X')
    ax.set_xlim3d(0, 10)
    ax.set_ylabel('Y')
    ax.set_ylim3d(-1, 4)
    ax.set_zlabel('Z')
    ax.set_zlim3d(0, 1)

    plt.show()
    print(f'type fig: {type(fig)}')
    # fig.savefig('Poly3D.png')
    return fig


def save_question():
    root = Tk()
    answer = mb.askquestion(
        title="Save control",
        message="Save picture?")
    root.mainloop()
    return answer


if __name__ == '__main__':
    poly_graph3d()
    #
    # print(save_question())

