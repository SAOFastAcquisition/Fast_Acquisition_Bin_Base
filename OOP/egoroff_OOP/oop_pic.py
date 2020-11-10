
import matplotlib.pyplot as plt


class Pic:
    def __init__(self, x, y, legend='curve'):
        self.x = x
        self.y = y
        self.legend = legend

    def line(self):
        ax = plt.plot()
        ax.plot(self.x, self.y, label=self.legend)
        plt.show()
        return ax


arg1 = [i for i in range(100)]
ord1 = [0.1 * i for i in range(100)]
ax = Pic(arg1, ord1)
plt.show()