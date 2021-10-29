from Interface.Window_D import Ui_MainWindow
import sys
from pathlib import Path
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import *
from PyQt5.QtCore import Qt
from Interface.parameters import param_dict
import pickle
import os

__all__ = ['main', 'ExampleApp']


class PermitAPP(QtWidgets.QMainWindow):

    def __init__(self):
        super().__init__()

        self.window = QtWidgets.QMainWindow()
        self.window.btn1


class ExampleApp(QtWidgets.QMainWindow, Ui_MainWindow):

    def __init__(self, param_dict_str):
        # Это здесь нужно для доступа к переменным, методам
        # и т.д. в файле design.py
        super().__init__()
        self.setupUi(self)  # Это нужно для инициализации нашего дизайна

        # Начальные установки
        self.freq_res = param_dict_str['freq_res']
        self.time_res = param_dict_str['time_res']
        self.frequency_mask = []
        self.index_freq_mask = None
        self.time_mask = []
        self.index_time_mask = None
        self.file_name = ' '
        self.file_folder = ' '
        self.catalog = param_dict_str['path_to_catalog']
        self.lne_frequency_resolution.setText(self.freq_res)
        self.lne_time_resolution.setText(self.time_res)
        self.gB_pattrens.adjustSize()
        self.param_dict = param_dict_str
        print(f'path to catalog: {self.catalog}')

        # self.btn = QPushButton('Attention!', self)
        # Работа с масками по частоте и времени. Устанавливаем начальные маски по частоте, которые взяты из
        # словаря param_dict_str по ключу 'freq_mask'.
        self.tableWidget_freq_patterns.setColumnWidth(0, 450)
        freq_mask = param_dict_str['freq_mask']
        k = 0
        print(f'ind_freq_mask: {self.index_freq_mask}, freq_mask = {self.frequency_mask}')
        for unit in freq_mask:
            item_freq = self.tableWidget_freq_patterns.item(k, 0)
            item_freq.setText(unit)
            # newItem.setForeground(QBrush(QColor(255, 0, 0)))
            self.tableWidget_freq_patterns.setItem(k, 0, item_freq)
            k += 1

        self.tableWidget_time_patterns.setColumnWidth(0, 450)
        time_mask = param_dict_str['time_mask']
        k = 0
        for unit in time_mask:
            item_time = self.tableWidget_time_patterns.item(k, 0)
            item_time.setText(unit)
            # newItem.setForeground(QBrush(QColor(255, 0, 0)))
            self.tableWidget_time_patterns.setItem(k, 0, item_time)
            k += 1

        self.btn_load_setup.adjustSize()
        self.btn_load_setup.clicked.connect(self.set_initial_setup)
        # Кнопки поиска/выбора файла для обработки и передачи управления обработчику выбора параметров
        self.btn_find_file.clicked.connect(self.find_processing_file)
        self.cbx_save_current_parameters.adjustSize()

        self.btn_set_parameters.adjustSize()
        self.btn_set_parameters.clicked.connect(self.set_parameter_handler)  # Выполнить функцию set_parameter_handler
        # при нажатии кнопки

    def __set_attr(self, name, value):
        self.__dict__[name] = value
        # print(self.__dict__[name])
        pass

    def set_initial_setup(self):
        init_dict = self.param_dict
        with open('save_param.bin', 'rb') as inp:
            param_set = pickle.load(inp)[-1]
        print(f'param_set = {param_set}')
        self.__set_attr('index_freq_mask', param_set['ind_freq_checkbox'])
        self.__set_attr('index_time_mask', param_set['ind_time_checkbox'])
        self.__set_attr('frequency_mask', param_set['freq_mask'])
        self.__set_attr('time_mask', param_set['time_mask'])
        self.__set_attr('freq_res', param_set['freq_res'])
        self.__set_attr('time_res', param_set['time_res'])
        self.__set_attr('file_name', param_set['file_name'])
        self.__set_attr('file_folder', param_set['file_folder'])
        self.__set_attr('path_to_catalog', param_set['path_to_catalog'])
        self.lne_frequency_resolution.setText(self.freq_res)
        self.lne_time_resolution.setText(self.time_res)
        self.__set_attr('catalog', str(Path(param_set['path_to_catalog'], param_set['file_folder'])))

        self.lne_oserv_file.setText(self.lne_oserv_file.text() + self.file_name)
        self.lne_file_name.setText(self.lne_file_name.text() + self.file_name)

        a = param_set['file_folder']
        print(f'path to catalog: {self.catalog}, file folder: {a}')
        print(f'ind_freq_mask: {self.index_freq_mask}, freq_mask = {self.frequency_mask}')
        freq_mask = init_dict['freq_mask']
        if self.index_freq_mask is not None:
            i = self.index_freq_mask
            try:
                freq_mask[i] = self.frequency_mask
            except IndexError:
                freq_mask.append(self.frequency_mask)

            # self.tableWidget_freq_patterns                 # setCheckState(i, 0, 2)
        k = 0
        for unit in freq_mask:
            item_freq = self.tableWidget_freq_patterns.item(k, 0)
            item_freq.setText(unit)
            # newItem.setForeground(QBrush(QColor(255, 0, 0)))
            self.tableWidget_freq_patterns.setItem(k, 0, item_freq)
            k += 1

        time_mask = init_dict['time_mask']
        k = 0
        if self.index_time_mask is not None:
            i = self.index_time_mask
            try:
                time_mask[i] = self.time_mask
            except IndexError:
                time_mask.append(self.time_mask)
            # self.tableWidget_time_patterns.setCheckState(i, 0, 2)
        for unit in time_mask:
            item_time = self.tableWidget_time_patterns.item(k, 0)
            item_time.setText(unit)
            # newItem.setForeground(QBrush(QColor(255, 0, 0)))
            self.tableWidget_time_patterns.setItem(k, 0, item_time)
            k += 1
        pass

    def __save_latest_setup(self, state):
        if state == Qt.Checked:
            print(state, type(state))
            print(f'Save_latest - file_folder: {self.file_folder}, file_name: {self.file_name}')
            # file = open('custom_parameters.txt', 'a')
            parameter_dict = {'freq_res': self.freq_res,
                              'time_res': self.time_res,
                              'freq_mask': self.frequency_mask,
                              'time_mask': self.time_mask,
                              'file_name': self.file_name,
                              'file_folder': self.file_folder,
                              'path_to_catalog': self.catalog,
                              'ind_freq_checkbox': self.index_freq_mask,
                              'ind_time_checkbox': self.index_time_mask}
            self.__save_parameters(parameter_dict)
        pass

    def __save_parameters(self, data):
        if not (os.path.isfile('save_param.bin')):
            head = [None]
            with open('save_param.bin', 'wb') as out:
                pickle.dump(head, out)

        with open('save_param.bin', 'rb') as inp:
            head = pickle.load(inp)
            head.append(data)
            print(head)
        with open('save_param.bin', 'wb') as out:
            pickle.dump(head, out)

    def set_parameter_handler(self):

        # Обработка выбранных на закладке параметров разрешения
        res_time = self.lne_time_resolution.text()
        res_frequency = self.lne_frequency_resolution.text()
        if int(res_time) < 8:
            message_box_time = QMessageBox.information(self, 'Time resolution verification',
                                                       'Your resolution too small or negative! '
                                                       'Set the resolution greater then 8')
            message_box_time.move(self, 200, 400)
        else:
            self.__set_attr('time_res', res_time)
        if int(res_frequency) < 1:
            message_box_freq = QMessageBox.critical(self, 'Frequency resolution verification',
                                                    'Your resolution too small or negative! '
                                                    'Set the resolution greater then 1')
        else:
            self.__set_attr('freq_res', res_frequency)

        # Работа с QTableWidget для выбора или редактирования масок по частоте и времени.
        # Маска По частоте
        n = self.tableWidget_freq_patterns.rowCount()
        count_choose = 0
        for i in range(n):
            item_freq = self.tableWidget_freq_patterns.item(i, 0)
            a_freq = item_freq.checkState()
            print(a_freq, item_freq)

            if a_freq == 2:
                self.__set_attr('index_freq_mask', i)
                count_choose += 1
                if count_choose > 1:
                    print("Wrong choose frequency mask!!!")
                else:
                    freq_mask_choose_str = item_freq.text()
                    self.__set_attr('frequency_mask', freq_mask_choose_str)
        if count_choose == 0:
            print('Choose frequency pattern')
            message_box_freq = QMessageBox.critical(self, 'Frequency pattern verification',
                                                    'Choose frequency pattern! '
                                                    'Set the frequency pattern!')
        # Маска По времени
        n = self.tableWidget_time_patterns.rowCount()
        count_choose = 0
        for i in range(n):
            item_time = self.tableWidget_time_patterns.item(i, 0)
            a_time = item_time.checkState()
            print(a_time, item_time)

            if a_time == 2:
                self.__set_attr('index_time_mask', i)
                count_choose += 1
                if count_choose > 1:
                    print("Wrong choose frequency mask!!!")
                else:
                    time_mask_choose_str = item_time.text()
                    self.__set_attr('time_mask', time_mask_choose_str)
        if count_choose == 0:
            print('Choose time pattern')
            message_box_time = QMessageBox.critical(self, 'Time pattern verification',
                                                    'Choose time pattern! '
                                                    'Set the time pattern!')
            # Сохранение текущих настроек
        self.__save_latest_setup(self.cbx_save_current_parameters.checkState())
        pass

    def find_processing_file(self):
        """ Функция обработчик принимает клик кнопкой и возвращает в виде аттрибутов объекта класса название
        файла для обработки и папки в котором он лежит"""
        path_to_folder = self.catalog  # Исходная папка поиска файла для обработки
        # Результат поиска
        res = QFileDialog.getOpenFileNames(self, 'Open file', path_to_folder, 'BinaryFile (*.bin);;NumpyFile (*.npy)')
        res_part = res[0][0].split('/')

        file_name = res_part[-1][:-4]
        file_folder = res_part[-2]
        self.__set_attr('file_name', file_name)
        self.__set_attr('file_folder', file_folder)
        print(f'file_folder: {file_folder}, file_name: {file_name}')
        self.lne_oserv_file.setText(self.lne_oserv_file.text() + self.file_name)
        self.lne_file_name.setText(self.lne_file_name.text() + self.file_name)

        print(f'file_folder: {self.file_folder}, file_name: {self.file_name}')
        pass


def main():
    app = QtWidgets.QApplication(sys.argv)  # Новый экземпляр QApplication
    param_dict_str = param_dict_to_str(param_dict())
    window = ExampleApp(param_dict_str)  # Создаём объект класса ExampleApp
    # if window.btn_load_setup.clicked.connect()
    window.show()  # Показываем окно
    # save_parameters(window)
    app.exec_()  # и запускаем приложение
    freq_mask_num = list_str_to_num(window.frequency_mask)
    time_mask_num = list_str_to_num(window.time_mask)

    # print(window.frequency_mask, type(window.frequency_mask))
    parameters_dict = {'freq_res': int(window.freq_res),
                       'time_res': int(window.time_res),
                       'freq_mask': freq_mask_num,
                       'time_mask': time_mask_num,
                       'file_name': window.file_name,
                       'file_folder': window.file_folder}
    return parameters_dict


def num_to_str(b):
    """ Принимает двумерный массив типа list, возвращает одномерный массив строк типа list """
    return list(map(lambda x: str(x)[1:-1], b))


def param_dict_to_str(dict):
    """ Принимает словарь с разноформатными параметрами обработки данных, возвращает словарь с параметрами
    преобразованными в строковый формат"""
    freq_mask_str = num_to_str(dict['freq_mask'])
    time_mask_str = num_to_str(dict['time_mask'])
    dict_str = {'freq_res': str(dict['freq_resolution']),
                'time_res': str(dict['time_resolution']),
                'freq_mask': freq_mask_str,
                'time_mask': time_mask_str,
                # 'file_folder': dict['file_folder'],
                'path_to_catalog': dict['path_to_catalog']}
    return dict_str


def list_str_to_num(list_str):
    return list(map(lambda x: int(x), list_str.split(', ')))


# def save_parameters(data):
#     if not (os.path.isfile('save_object.bin')):
#         head = [None]
#         with open('save_object.bin', 'wb') as out:
#             pickle.dump(head, out)
#
#     with open('save_object.bin', 'rb') as inp:
#         head = pickle.load(inp)
#         head.append(data)
#         print(head)
#     with open('save_object.bin', 'wb') as out:
#         pickle.dump(head, out)


if __name__ == '__main__':
    parameter = main()
    time_res = parameter['time_res']
    freq_res = parameter['freq_res']
    file_name = parameter['file_name']
    file_folder = parameter['file_folder']
    freq_mask = parameter['freq_mask']
    time_mask = parameter['time_mask']
    print(f'time_res = {time_res}, freq_res = {freq_res}, file_name = {file_name}, file_folder= {file_folder},'
          f' freq_mask = {freq_mask}, time_mask = {time_mask}')

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
