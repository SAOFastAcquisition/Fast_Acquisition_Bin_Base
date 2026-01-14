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
def poly_graph3d_model():
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
                                             cc('_y')])
    poly.set_alpha(0.7)
    ax.add_collection3d(poly, zs=zs, zdir='_y')

    ax.set_xlabel('Time, sec')
    # ax.set_xlim3d(0, 10)
    ax.set_ylabel('Frequency, MHz')
    # ax.set_ylim3d(-1, 4)
    ax.set_zlabel('Spectrum, arb. un.')
    # ax.set_zlim3d(0, 1)

    plt.show()
    print(f'type fig: {type(fig)}')
    # fig.savefig('Poly3D.png')
    return fig


def poly_graph3d(*args):
    xs = args[0]    # Вектор значений аргумента (в нашем случае время)
    ys_arr = args[1]    # Массив значений ординаты при значениях параметра из вектора zs (скан по времени на фиксированной
                    # частоте)
    zs = args[2]    # Вектор значений параметра (в нашем случае частоты)
    _l = len(zs)
    verts = []
    for i in range(_l):
        ys = ys_arr[:, i]
        verts.append(list(zip(xs, ys)))

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    colours_map = [0]*_l
    for i in range(_l // 2):
        colours_map[2 * i] = cc('r')
        colours_map[2 * i + 1] = cc('b')
        if _l % 2:
            colours_map[_l - 1] = cc('r')
    # colours_map = [cc('r'), cc('g'), cc('b'), cc('_y'),
    #                cc('r'), cc('g'), cc('b'), cc('_y'), cc('g'),
    #                cc('r'), cc('g'), cc('b'),
    #                cc('_y')]

    poly = PolyCollection(verts, facecolors=colours_map)
    # poly = PolyCollection(verts)
    poly.set_alpha(0.6)
    ax.add_collection3d(poly, zs=zs, zdir='_y')

    ax.set_xlabel('Time, sec')
    ax.set_xlim3d(0, 400)
    ax.set_ylabel('Frequency, MHz')
    ax.set_ylim3d(1000, 3000)
    ax.set_zlabel('Spectrum, arb. un.')
    ax.set_zlim3d(0, 2e8)

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
    poly_graph3d_model()
    #
    # print(save_question())

