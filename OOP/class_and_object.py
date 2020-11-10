import matplotlib.pyplot as plt


class Gates:
    # arg = [i for i in range(100)]
    # ord = [0.1 * i for i in range(100)]
    def __init__(self, arg=[i for i in range(100)], ord=[i for i in range(100)]):
        if len(arg) == len(ord):
            self.arg = arg
            self.ord = ord
        else:
            raise ValueError('Размеры аргумента и ординаты должны быть равными')

    def graph(self):
        fig, ax = plt.subplots(1, figsize=(12, 6))
        if hasattr(self, 'legend'):
            ax.plot(self.arg, self.ord, label=self.legend)
            plt.legend(loc='upper right')
        else:
            ax.plot(self.arg, self.ord)
        plt.show()

    def set_legend(self, legend):
        self.legend = legend

arg = [i for i in range(100)]
ord = [0.1 * i for i in range(100)]
c = Gates(arg, ord)
legend = 'Graph'
c.set_legend(legend)
#  Создание нового аттрубута stock
Gates.stock = 'clock'

c.graph()
