# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'Window_D.ui'
#
# Created by: PyQt5 UI code generator 5.15.5
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(783, 905)
        MainWindow.setStyleSheet("background-color: rgb(101, 255, 255);")
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.centralwidget)
        self.verticalLayout.setObjectName("verticalLayout")
        self.tabWidget_1 = QtWidgets.QTabWidget(self.centralwidget)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.tabWidget_1.setFont(font)
        self.tabWidget_1.setTabPosition(QtWidgets.QTabWidget.North)
        self.tabWidget_1.setMovable(True)
        self.tabWidget_1.setObjectName("tabWidget_1")
        self.tab_parameters = QtWidgets.QWidget()
        font = QtGui.QFont()
        font.setPointSize(12)
        self.tab_parameters.setFont(font)
        self.tab_parameters.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.tab_parameters.setObjectName("tab_parameters")
        self.verticalLayoutWidget = QtWidgets.QWidget(self.tab_parameters)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(0, 190, 761, 251))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.gB_resolution = QtWidgets.QGroupBox(self.tab_parameters)
        self.gB_resolution.setGeometry(QtCore.QRect(130, 10, 451, 131))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.gB_resolution.setFont(font)
        self.gB_resolution.setStyleSheet("background-color: rgb(170, 255, 255);")
        self.gB_resolution.setObjectName("gB_resolution")
        self.lbl_frequency_resolution = QtWidgets.QLabel(self.gB_resolution)
        self.lbl_frequency_resolution.setGeometry(QtCore.QRect(10, 60, 181, 25))
        self.lbl_frequency_resolution.setStyleSheet("background-color: rgb(233, 185, 110);\n"
"font: 14pt \"Ubuntu\";\n"
"border-color: rgb(0, 0, 0);")
        self.lbl_frequency_resolution.setObjectName("lbl_frequency_resolution")
        self.lbl_time_resolution = QtWidgets.QLabel(self.gB_resolution)
        self.lbl_time_resolution.setGeometry(QtCore.QRect(10, 90, 181, 25))
        self.lbl_time_resolution.setStyleSheet("background-color: rgb(233, 185, 110);\n"
"font: 14pt \"Ubuntu\";\n"
"border-color: rgb(0, 0, 0);")
        self.lbl_time_resolution.setObjectName("lbl_time_resolution")
        self.lne_time_resolution = QtWidgets.QLineEdit(self.gB_resolution)
        self.lne_time_resolution.setGeometry(QtCore.QRect(210, 90, 141, 25))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.lne_time_resolution.setFont(font)
        self.lne_time_resolution.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.lne_time_resolution.setObjectName("lne_time_resolution")
        self.lne_frequency_resolution = QtWidgets.QLineEdit(self.gB_resolution)
        self.lne_frequency_resolution.setGeometry(QtCore.QRect(210, 60, 141, 25))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(10)
        sizePolicy.setHeightForWidth(self.lne_frequency_resolution.sizePolicy().hasHeightForWidth())
        self.lne_frequency_resolution.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.lne_frequency_resolution.setFont(font)
        self.lne_frequency_resolution.setStyleSheet("background-color: rgb(255, 255, 255);\n"
"border-color: rgb(0, 0, 0);")
        self.lne_frequency_resolution.setObjectName("lne_frequency_resolution")
        self.lbl_time_unit = QtWidgets.QLabel(self.gB_resolution)
        self.lbl_time_unit.setGeometry(QtCore.QRect(370, 90, 67, 25))
        self.lbl_time_unit.setStyleSheet("background-color: rgb(233, 185, 110);")
        self.lbl_time_unit.setObjectName("lbl_time_unit")
        self.lbl_frequency_unit = QtWidgets.QLabel(self.gB_resolution)
        self.lbl_frequency_unit.setGeometry(QtCore.QRect(370, 60, 67, 25))
        self.lbl_frequency_unit.setStyleSheet("background-color: rgb(233, 185, 110);\n"
"font: 14pt \"Ubuntu\";")
        self.lbl_frequency_unit.setObjectName("lbl_frequency_unit")
        self.lne_oserv_file = QtWidgets.QLineEdit(self.gB_resolution)
        self.lne_oserv_file.setGeometry(QtCore.QRect(50, 20, 351, 25))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.lne_oserv_file.setFont(font)
        self.lne_oserv_file.setStyleSheet("background-color: rgb(255, 170, 255);\n"
"background-color: rgb(255, 250, 183);")
        self.lne_oserv_file.setObjectName("lne_oserv_file")
        self.btn_set_parameters = QtWidgets.QPushButton(self.tab_parameters)
        self.btn_set_parameters.setGeometry(QtCore.QRect(280, 750, 181, 31))
        self.btn_set_parameters.setStyleSheet("background-color: rgb(169, 191, 64);\n"
"font: 14pt \"Ubuntu\";\n"
"alternate-background-color: rgba(166, 64, 191, 28);")
        self.btn_set_parameters.setAutoDefault(False)
        self.btn_set_parameters.setDefault(False)
        self.btn_set_parameters.setFlat(False)
        self.btn_set_parameters.setObjectName("btn_set_parameters")
        self.gB_pattrens = QtWidgets.QGroupBox(self.tab_parameters)
        self.gB_pattrens.setGeometry(QtCore.QRect(10, 150, 731, 481))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.gB_pattrens.setFont(font)
        self.gB_pattrens.setStyleSheet("background-color: rgb(170, 255, 255);")
        self.gB_pattrens.setObjectName("gB_pattrens")
        self.lbl_time_mask = QtWidgets.QLabel(self.gB_pattrens)
        self.lbl_time_mask.setGeometry(QtCore.QRect(0, 260, 731, 31))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.lbl_time_mask.setFont(font)
        self.lbl_time_mask.setStyleSheet("background-color: rgb(170, 170, 0);")
        self.lbl_time_mask.setObjectName("lbl_time_mask")
        self.tableWidget_time_patterns = QtWidgets.QTableWidget(self.gB_pattrens)
        self.tableWidget_time_patterns.setGeometry(QtCore.QRect(90, 300, 561, 181))
        self.tableWidget_time_patterns.setMinimumSize(QtCore.QSize(561, 0))
        self.tableWidget_time_patterns.setStyleSheet("background-color: rgb(252, 255, 199);")
        self.tableWidget_time_patterns.setObjectName("tableWidget_time_patterns")
        self.tableWidget_time_patterns.setColumnCount(1)
        self.tableWidget_time_patterns.setRowCount(4)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget_time_patterns.setVerticalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget_time_patterns.setVerticalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget_time_patterns.setVerticalHeaderItem(2, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget_time_patterns.setVerticalHeaderItem(3, item)
        item = QtWidgets.QTableWidgetItem()
        item.setText("0")
        self.tableWidget_time_patterns.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        item.setCheckState(QtCore.Qt.Unchecked)
        self.tableWidget_time_patterns.setItem(0, 0, item)
        item = QtWidgets.QTableWidgetItem()
        item.setCheckState(QtCore.Qt.Unchecked)
        self.tableWidget_time_patterns.setItem(1, 0, item)
        item = QtWidgets.QTableWidgetItem()
        item.setCheckState(QtCore.Qt.Unchecked)
        self.tableWidget_time_patterns.setItem(2, 0, item)
        item = QtWidgets.QTableWidgetItem()
        item.setCheckState(QtCore.Qt.Unchecked)
        self.tableWidget_time_patterns.setItem(3, 0, item)
        self.tableWidget_freq_patterns = QtWidgets.QTableWidget(self.gB_pattrens)
        self.tableWidget_freq_patterns.setGeometry(QtCore.QRect(90, 60, 561, 192))
        self.tableWidget_freq_patterns.setMinimumSize(QtCore.QSize(561, 0))
        self.tableWidget_freq_patterns.setStyleSheet("background-color: rgb(252, 255, 199);")
        self.tableWidget_freq_patterns.setObjectName("tableWidget_freq_patterns")
        self.tableWidget_freq_patterns.setColumnCount(1)
        self.tableWidget_freq_patterns.setRowCount(4)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget_freq_patterns.setVerticalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget_freq_patterns.setVerticalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget_freq_patterns.setVerticalHeaderItem(2, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget_freq_patterns.setVerticalHeaderItem(3, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignVCenter)
        font = QtGui.QFont()
        font.setPointSize(10)
        item.setFont(font)
        self.tableWidget_freq_patterns.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        item.setCheckState(QtCore.Qt.Unchecked)
        self.tableWidget_freq_patterns.setItem(0, 0, item)
        item = QtWidgets.QTableWidgetItem()
        item.setCheckState(QtCore.Qt.Unchecked)
        self.tableWidget_freq_patterns.setItem(1, 0, item)
        item = QtWidgets.QTableWidgetItem()
        item.setCheckState(QtCore.Qt.Unchecked)
        self.tableWidget_freq_patterns.setItem(2, 0, item)
        item = QtWidgets.QTableWidgetItem()
        item.setCheckState(QtCore.Qt.Unchecked)
        self.tableWidget_freq_patterns.setItem(3, 0, item)
        self.lbl_freq_mask = QtWidgets.QLabel(self.gB_pattrens)
        self.lbl_freq_mask.setGeometry(QtCore.QRect(0, 20, 731, 31))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.lbl_freq_mask.setFont(font)
        self.lbl_freq_mask.setStyleSheet("background-color: rgb(170, 170, 0);")
        self.lbl_freq_mask.setObjectName("lbl_freq_mask")
        self.cbx_save_current_parameters = QtWidgets.QCheckBox(self.tab_parameters)
        self.cbx_save_current_parameters.setGeometry(QtCore.QRect(100, 680, 261, 31))
        self.cbx_save_current_parameters.setStyleSheet("font: 12pt \"MS Shell Dlg 2\";\n"
"background-color: rgb(170, 170, 0);")
        self.cbx_save_current_parameters.setObjectName("cbx_save_current_parameters")
        self.groupBox = QtWidgets.QGroupBox(self.tab_parameters)
        self.groupBox.setGeometry(QtCore.QRect(470, 640, 221, 121))
        self.groupBox.setStyleSheet("background-color: rgb(170, 255, 255);\n"
"font: 10pt \"MS Shell Dlg 2\";")
        self.groupBox.setObjectName("groupBox")
        self.lbl_choise_picture = QtWidgets.QLabel(self.groupBox)
        self.lbl_choise_picture.setGeometry(QtCore.QRect(10, 20, 201, 31))
        self.lbl_choise_picture.setStyleSheet("font: 12pt \"MS Shell Dlg 2\";\n"
"background-color: rgb(170, 170, 0);")
        self.lbl_choise_picture.setObjectName("lbl_choise_picture")
        self.rbn_one_window = QtWidgets.QRadioButton(self.groupBox)
        self.rbn_one_window.setGeometry(QtCore.QRect(10, 60, 201, 21))
        self.rbn_one_window.setStyleSheet("font: 10pt \"MS Shell Dlg 2\";\n"
"background-color: rgb(255, 170, 0);")
        self.rbn_one_window.setChecked(True)
        self.rbn_one_window.setObjectName("rbn_one_window")
        self.rbn_multi_window = QtWidgets.QRadioButton(self.groupBox)
        self.rbn_multi_window.setGeometry(QtCore.QRect(10, 90, 201, 21))
        self.rbn_multi_window.setStyleSheet("font: 10pt \"MS Shell Dlg 2\";\n"
"background-color: rgb(255, 170, 0);")
        self.rbn_multi_window.setObjectName("rbn_multi_window")
        self.tabWidget_1.addTab(self.tab_parameters, "")
        self.tab_choice_input = QtWidgets.QWidget()
        self.tab_choice_input.setObjectName("tab_choice_input")
        self.btn_find_file = QtWidgets.QPushButton(self.tab_choice_input)
        self.btn_find_file.setGeometry(QtCore.QRect(0, 90, 231, 31))
        self.btn_find_file.setStyleSheet("font: 14pt \"MS Shell Dlg 2\";\n"
"background-color: rgb(183, 255, 191);")
        self.btn_find_file.setObjectName("btn_find_file")
        self.lne_file_name = QtWidgets.QLineEdit(self.tab_choice_input)
        self.lne_file_name.setGeometry(QtCore.QRect(260, 90, 311, 31))
        self.lne_file_name.setStyleSheet("font: 12pt \"MS Shell Dlg 2\";\n"
"background-color: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:0, stop:0 rgba(255, 178, 102, 255), stop:0.55 rgba(235, 148, 61, 255), stop:0.98 rgba(0, 0, 0, 255), stop:1 rgba(0, 0, 0, 0));\n"
"background-color: rgb(225, 255, 225);")
        self.lne_file_name.setObjectName("lne_file_name")
        self.btn_load_setup = QtWidgets.QPushButton(self.tab_choice_input)
        self.btn_load_setup.setGeometry(QtCore.QRect(200, 20, 161, 41))
        self.btn_load_setup.setStyleSheet("font: 12pt \"MS Shell Dlg 2\";\n"
"background-color: rgb(170, 170, 0);")
        self.btn_load_setup.setObjectName("btn_load_setup")
        self.tabWidget_1.addTab(self.tab_choice_input, "")
        self.verticalLayout.addWidget(self.tabWidget_1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 783, 26))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        self.tabWidget_1.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.gB_resolution.setTitle(_translate("MainWindow", "Resolution Box"))
        self.lbl_frequency_resolution.setText(_translate("MainWindow", "Frequency resolution"))
        self.lbl_time_resolution.setText(_translate("MainWindow", "Time resolution"))
        self.lbl_time_unit.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-size:14pt;\">ms</span></p></body></html>"))
        self.lbl_frequency_unit.setText(_translate("MainWindow", "MHz"))
        self.lne_oserv_file.setText(_translate("MainWindow", "Observation File: "))
        self.btn_set_parameters.setText(_translate("MainWindow", "Set parameters"))
        self.gB_pattrens.setTitle(_translate("MainWindow", "Frequency and Time Patterns Box"))
        self.lbl_time_mask.setText(_translate("MainWindow", "<html><head/><body><p align=\"center\">Time Mask, sec</p></body></html>"))
        item = self.tableWidget_time_patterns.verticalHeaderItem(0)
        item.setText(_translate("MainWindow", "Pattern1"))
        item = self.tableWidget_time_patterns.verticalHeaderItem(1)
        item.setText(_translate("MainWindow", "Pattern2"))
        item = self.tableWidget_time_patterns.verticalHeaderItem(2)
        item.setText(_translate("MainWindow", "Pattern3"))
        item = self.tableWidget_time_patterns.verticalHeaderItem(3)
        item.setText(_translate("MainWindow", "Pattern4"))
        __sortingEnabled = self.tableWidget_time_patterns.isSortingEnabled()
        self.tableWidget_time_patterns.setSortingEnabled(False)
        self.tableWidget_time_patterns.setSortingEnabled(__sortingEnabled)
        item = self.tableWidget_freq_patterns.verticalHeaderItem(0)
        item.setText(_translate("MainWindow", "Pattern1"))
        item = self.tableWidget_freq_patterns.verticalHeaderItem(1)
        item.setText(_translate("MainWindow", "Pattern2"))
        item = self.tableWidget_freq_patterns.verticalHeaderItem(2)
        item.setText(_translate("MainWindow", "Pattern3"))
        item = self.tableWidget_freq_patterns.verticalHeaderItem(3)
        item.setText(_translate("MainWindow", "Pattern4"))
        item = self.tableWidget_freq_patterns.horizontalHeaderItem(0)
        item.setText(_translate("MainWindow", "Frequencies, MHz"))
        __sortingEnabled = self.tableWidget_freq_patterns.isSortingEnabled()
        self.tableWidget_freq_patterns.setSortingEnabled(False)
        self.tableWidget_freq_patterns.setSortingEnabled(__sortingEnabled)
        self.lbl_freq_mask.setText(_translate("MainWindow", "<html><head/><body><p align=\"center\">Frequency Mask, MHz</p></body></html>"))
        self.cbx_save_current_parameters.setText(_translate("MainWindow", "Save current parameters"))
        self.groupBox.setTitle(_translate("MainWindow", "Output picture Box"))
        self.lbl_choise_picture.setText(_translate("MainWindow", "Output picture mode"))
        self.rbn_one_window.setText(_translate("MainWindow", "One window picture"))
        self.rbn_multi_window.setText(_translate("MainWindow", "Multi window picture"))
        self.tabWidget_1.setTabText(self.tabWidget_1.indexOf(self.tab_parameters), _translate("MainWindow", "Processing parameters"))
        self.btn_find_file.setText(_translate("MainWindow", "Find file for processing"))
        self.lne_file_name.setText(_translate("MainWindow", "Result:"))
        self.btn_load_setup.setText(_translate("MainWindow", "Load the last set-up"))
        self.tabWidget_1.setTabText(self.tabWidget_1.indexOf(self.tab_choice_input), _translate("MainWindow", "Input data"))
