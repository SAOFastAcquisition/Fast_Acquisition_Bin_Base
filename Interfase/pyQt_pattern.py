import sys
from PyQt5.QtWidgets import *

# app = QApplication(sys.argv)
#
# dlgMain = QDialog()
# dlgMain.setWindowTitle('Dialog Window')
#
#
# dlgMain.show()
# sys.exit"(app.exec_())


class DlgMain(QDialog):
    def __init__(self):
        super(DlgMain, self).__init__()

        self.setWindowTitle('Dialog Window')
        self.resize(400, 300)

        self.label_window_change = QLineEdit('Standard Window', self)
        self.label_window_change.move(120, 50)

        self.btn_change_title_Window = QPushButton('Change Window Title', self)
        self.btn_change_title_Window.move(120, 70)
        self.btn_change_title_Window.clicked.connect(self.change_window_title)

    def change_window_title(self):
        new_title = self.label_window_change.text()
        self.setWindowTitle(new_title)


def application():
    app = QApplication(sys.argv)

    dlgMain = DlgMain()
    dlgMain.show()

    sys.exit(app.exec_())


if __name__ == '__main__':
    application()

