# from Window_D import Ui_MainWindow
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import *
import pickle
from window_interface import main
import os
import sys


class Dialog(QMessageBox):
    def __init__(self):
        super().__init__()
        self.control_picture = None
        if self.question(self, 'Picture control', 'Save the picture?') == QMessageBox.Yes:
            val = 'yes'
        else:
            val = 'no'
        self.__setattr('control_picture', val)
        pass

    def __setattr(self, key, value):
        self.__dict__[key] = value



def exec_app():
    app = QApplication(sys.argv)
    dialog_win = Dialog()

    dialog_win.show()
    app.exec()
    # print(f'control Picture = {dialog_win.control_picture}')
    return dialog_win.control_picture



if not (os.path.isfile('c_param.bin')):
    with open('c_param.bin', 'wb') as out:
        pickle.dump([], out)

parameter_dict = {'freq_resolution': 1,
                  'freq_mask': [[1080, 1140, 1360, 1420, 1620, 1780, 1980],
                                [2060, 2220, 2300, 2500, 2560, 2700, 2800, 2880, 2980],
                                [1050, 1465, 1535, 1600, 1700, 2265, 2550, 2700, 2800, 2920]],
                  'time_mask': [],
                  'path_to_catalog': 'hernya'}

# with open('c_param.bin', 'rb') as inp:
#     head = pickle.load(inp)
#     head.append(parameter_dict)
#     print(head)
# with open('c_param.bin', 'wb') as out:
#     pickle.dump(head, out)


def param_dict_to_str(dict):
    return list(map(lambda x: str(x)[1:-1], dict['freq_mask']))


mask_str = param_dict_to_str(parameter_dict)
print(mask_str)
with open('save_param.bin', 'rb') as inp:
    head = pickle.load(inp)

print(head[-1])
flag = exec_app()
print(f'flag = {flag}')
