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
        self.btnBrowse.clicked.connect(self.browse_folder)  # Выполнить функцию browse_folder
        # при нажатии кнопки

    def browse_folder(self):
        self.listWidget.clear()  # На случай, если в списке уже есть элементы
        directory = QtWidgets.QFileDialog.getExistingDirectory(self, "ClickButton")
        # открыть диалог выбора директории и установить значение переменной
        # равной пути к выбранной директории

        if directory:  # не продолжать выполнение, если пользователь не выбрал директорию
            for file_name in os.listdir(directory):  # для каждого файла в директории
                self.listWidget.addItem(file_name)   # добавить файл в listWidget


def main():
    app = QtWidgets.QApplication(sys.argv)  # Новый экземпляр QApplication
    a = sys.argv
    window = ExampleApp()  # Создаём объект класса ExampleApp
    window.show()  # Показываем окно
    app.exec_()  # и запускаем приложение


if __name__ == '__main__':

    main()
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