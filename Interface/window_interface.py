from Window_D import Ui_MainWindow
import sys
from PyQt5 import QtWidgets
import os


class ExampleApp(QtWidgets.QMainWindow, Ui_MainWindow):

    def __init__(self):
        # Это здесь нужно для доступа к переменным, методам
        # и т.д. в файле design.py
        super().__init__()
        self.setupUi(self)  # Это нужно для инициализации нашего дизайна
        self.btn_set_parameters.clicked.connect(self.set_parameter_handler)  # Выполнить функцию browse_folder
        # при нажатии кнопки
        self.freq_res = 1
        self.time_res = 4

    def set_atr(self, name, value):
        self.__dict__[name] = value
        # print(self.__dict__[name])
        pass

    def set_parameter_handler(self):
        res_time = self.lne_time_resolution.text()
        res_frequency = self.lne_frequency_resolution.text()
        self.set_atr('freq_res', res_frequency)
        self.set_atr('time_res', res_time)


def main():
    app = QtWidgets.QApplication(sys.argv)  # Новый экземпляр QApplication
    window = ExampleApp()  # Создаём объект класса ExampleApp
    window.show()  # Показываем окно
    app.exec_()  # и запускаем приложение
    return int(window.time_res), int(window.freq_res)


if __name__ == '__main__':
    time_res, frequency_res = main()
    print(time_res, frequency_res)

# # объект приложения
# app = QApplication(sys.argv)
#
# # QMainWindow объект
# main_window = QMainWindow()
#
# # Это класс Ui_MainWindow, реализованный дизайнером qt
# ui_components = Ui_MainWindow()
# # Вызвать метод setupUi () для регистрации в объекте QMainWindow
# ui_components.setupUi(main_window)

# Шоу
# main_window.show()
#
# sys.exit(app.exec_())


# if __name__ == '__main__':
#     print(call_wind())


# self.listWidget.clear()  # На случай, если в списке уже есть элементы
# directory = QtWidgets.QFileDialog.getExistingDirectory(self, "ClickButton")
# # открыть диалог выбора директории и установить значение переменной
# # равной пути к выбранной директории
#
# if directory:  # не продолжать выполнение, если пользователь не выбрал директорию
#     for file_name in os.listdir(directory):  # для каждого файла в директории
#         self.listWidget.addItem(file_name)   # добавить файл в listWidget