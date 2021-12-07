import sys
from PyQt5.QtWidgets import QApplication, QMainWindow
from PyQt5 import QtCore, uic, QtGui, QtWidgets
from PyQt5.QtGui import QIcon, QPixmap
import femmt as ft


class WindingControls(QtWidgets.QWidget):
    def __init__(self, index, *args, **kwargs):
        super(WindingControls, self).__init__(*args, **kwargs)
        self.windingDynVerticalLayout = QtWidgets.QVBoxLayout(self)
        self.windingDynVerticalLayout.setContentsMargins(0, 0, 0, 0)
        self.windingDynVerticalLayout.setObjectName("windingDynVerticalLayout" + str(index))
        self.windingDynHorizontalLayout = QtWidgets.QHBoxLayout()
        self.windingDynHorizontalLayout.setObjectName("windingDynHorizontalLayout")
        self.windingTypeLabel = QtWidgets.QLabel(self)
        font = QtGui.QFont()
        font.setFamily("Segoe UI Semibold")
        font.setPointSize(8)
        self.windingTypeLabel.setFont(font)
        self.windingTypeLabel.setObjectName("windingTypeLabel")
        self.windingDynHorizontalLayout.addWidget(self.windingTypeLabel)
        self.windingTypeLineEdit = QtWidgets.QComboBox(self)
        self.windingTypeLineEdit.setObjectName("windingTypeLineEdit")
        self.windingDynHorizontalLayout.addWidget(self.windingTypeLineEdit)
        self.windingMaterialLineEdit = QtWidgets.QComboBox(self)
        self.windingMaterialLineEdit.setObjectName("windingMaterialLineEdit")
        self.windingDynHorizontalLayout.addWidget(self.windingMaterialLineEdit)
        self.windingDynVerticalLayout.addLayout(self.windingDynHorizontalLayout)
        self.windingDynGroupBox = QtWidgets.QGroupBox(self)
        self.windingDynGroupBox.setTitle("")
        self.windingDynGroupBox.setFlat(False)
        self.windingDynGroupBox.setObjectName("windingDynGroupBox")
        self.windingDynGroupBoxHLayout = QtWidgets.QHBoxLayout(self.windingDynGroupBox)
        self.windingDynGroupBoxHLayout.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        self.windingDynGroupBoxHLayout.setObjectName("windingDynGroupBoxHLayout")
        self.strandLabel = QtWidgets.QLabel(self.windingDynGroupBox)
        self.strandLabel.setFont(font)
        self.strandLabel.setObjectName("strandLabel")
        self.windingDynGroupBoxHLayout.addWidget(self.strandLabel)
        self.strandLineEdit = QtWidgets.QLineEdit(self.windingDynGroupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.strandLineEdit.sizePolicy().hasHeightForWidth())
        self.strandLineEdit.setSizePolicy(sizePolicy)
        self.strandLineEdit.setFont(font)
        self.strandLineEdit.setMinimumSize(QtCore.QSize(30, 22))
        self.strandLineEdit.setObjectName("strandLineEdit")
        self.windingDynGroupBoxHLayout.addWidget(self.strandLineEdit)
        self.mmLabelOne = QtWidgets.QLabel(self.windingDynGroupBox)
        self.mmLabelOne.setFont(font)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.mmLabelOne.sizePolicy().hasHeightForWidth())
        self.mmLabelOne.setSizePolicy(sizePolicy)
        self.mmLabelOne.setObjectName("mmLabelOne")
        self.windingDynGroupBoxHLayout.addWidget(self.mmLabelOne)
        self.strdiameterLabel = QtWidgets.QLabel(self.windingDynGroupBox)
        self.strdiameterLabel.setFont(font)
        self.strdiameterLabel.setObjectName("strdiameterLabel")
        self.windingDynGroupBoxHLayout.addWidget(self.strdiameterLabel)
        self.strdiameterLineEdit = QtWidgets.QLineEdit(self.windingDynGroupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.strdiameterLineEdit.sizePolicy().hasHeightForWidth())
        self.strdiameterLineEdit.setSizePolicy(sizePolicy)
        self.strdiameterLineEdit.setFont(font)
        self.strdiameterLineEdit.setMinimumSize(QtCore.QSize(30, 22))
        self.strdiameterLineEdit.setObjectName("strdiameterLineEdit")
        self.windingDynGroupBoxHLayout.addWidget(self.strdiameterLineEdit)
        self.mmLabelTwo = QtWidgets.QLabel(self.windingDynGroupBox)
        self.mmLabelTwo.setFont(font)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.mmLabelTwo.sizePolicy().hasHeightForWidth())
        self.mmLabelTwo.setSizePolicy(sizePolicy)
        self.mmLabelTwo.setObjectName("mmLabelTwo")
        self.windingDynGroupBoxHLayout.addWidget(self.mmLabelTwo)
        self.diameterLabel = QtWidgets.QLabel(self.windingDynGroupBox)
        self.diameterLabel.setFont(font)
        self.diameterLabel.setObjectName("diameterLabel")
        self.windingDynGroupBoxHLayout.addWidget(self.diameterLabel)
        self.diameterLineEdit = QtWidgets.QLineEdit(self.windingDynGroupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.diameterLineEdit.sizePolicy().hasHeightForWidth())
        self.diameterLineEdit.setSizePolicy(sizePolicy)
        self.diameterLineEdit.setFont(font)
        self.diameterLineEdit.setMinimumSize(QtCore.QSize(30, 22))
        self.diameterLineEdit.setObjectName("diameterLineEdit")
        self.windingDynGroupBoxHLayout.addWidget(self.diameterLineEdit)
        self.mmLabelThree = QtWidgets.QLabel(self.windingDynGroupBox)
        self.mmLabelThree.setFont(font)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.mmLabelThree.sizePolicy().hasHeightForWidth())
        self.mmLabelThree.setSizePolicy(sizePolicy)
        self.mmLabelThree.setObjectName("mmLabelThree")
        self.windingDynGroupBoxHLayout.addWidget(self.mmLabelThree)
        self.fillfactorLabel = QtWidgets.QLabel(self.windingDynGroupBox)
        self.fillfactorLabel.setFont(font)
        self.fillfactorLabel.setObjectName("fillfactorLabel")
        self.windingDynGroupBoxHLayout.addWidget(self.fillfactorLabel)
        self.fillfactorLineEdit = QtWidgets.QLineEdit(self.windingDynGroupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.fillfactorLineEdit.sizePolicy().hasHeightForWidth())
        self.fillfactorLineEdit.setSizePolicy(sizePolicy)
        self.fillfactorLineEdit.setFont(font)
        self.fillfactorLineEdit.setMinimumSize(QtCore.QSize(25, 22))
        self.fillfactorLineEdit.setText("")
        self.fillfactorLineEdit.setObjectName("fillfactorLineEdit")
        self.windingDynGroupBoxHLayout.addWidget(self.fillfactorLineEdit)
        self.percentLabel = QtWidgets.QLabel(self.windingDynGroupBox)
        self.percentLabel.setFont(font)
        self.percentLabel.setObjectName("percentLabel")
        self.windingDynGroupBoxHLayout.addWidget(self.percentLabel)
        self.windingDynVerticalLayout.addWidget(self.windingDynGroupBox)
        self.index = index
        self.retranslateUi()

    def retranslateUi(self):
        _translate = QtCore.QCoreApplication.translate
        if self.index == 1:
            text_to_attach = '(P)'
        elif self.index == 2:
            text_to_attach = '(S)'
        elif self.index == 3:
            text_to_attach = '(T)'
        self.windingTypeLabel.setText(_translate("Form", f"<html><head/><body><p><span style=\" color:#ff0000;\"> {self.index}. Winding Type {text_to_attach}</span></p></body></html>"))
        self.strandLabel.setText(_translate("Form", "<html><head/><body><p><span style=\" color:#ff0000;\">Strand No.</span></p></body></html>"))
        self.mmLabelOne.setText(_translate("Form", "mm"))
        self.strdiameterLabel.setText(_translate("Form", "<html><head/><body><p><span style=\" color:#ff0000;\">Str.Diameter</span></p></body></html>"))
        self.mmLabelTwo.setText(_translate("Form", "mm"))
        self.diameterLabel.setText(_translate("Form", "<html><head/><body><p><span style=\" color:#ff0000;\">Diameter</span></p></body></html>"))
        self.mmLabelThree.setText(_translate("Form", "mm"))
        self.fillfactorLabel.setText(_translate("Form", "<html><head/><body><p><span style=\" color:#ff0000;\">Fill Factor</span></p></body></html>"))
        self.percentLabel.setText(_translate("Form", " %"))


class AirGapControls(QtWidgets.QWidget):
    def __init__(self, index, *args, **kwargs):
        super(AirGapControls, self).__init__(*args, **kwargs)
        self.airgapDynamicHLayout = QtWidgets.QHBoxLayout(self)
        self.airgapDynamicHLayout.setContentsMargins(0, 0, 0, 0)
        self.airgapDynamicHLayout.setSpacing(2)
        self.airgapDynamicHLayout.setObjectName("airgapDynamicHLayout")
        # spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        # self.airgapDynamicHLayout.addItem(spacerItem)
        self.countDynLabel = QtWidgets.QLabel(self)
        font = QtGui.QFont()
        font.setPointSize(8)
        font.setFamily("Segoe UI Semibold")
        font.setItalic(True)
        self.countDynLabel.setFont(font)
        self.countDynLabel.setObjectName("countDynLabel")
        self.airgapDynamicHLayout.addWidget(self.countDynLabel)
        self.p2pDynLabel = QtWidgets.QLabel(self)
        font.setItalic(False)
        font.setBold(True)
        font.setWeight(75)
        self.p2pDynLabel.setFont(font)
        self.p2pDynLabel.setObjectName("p2pDynLabel")
        self.airgapDynamicHLayout.addWidget(self.p2pDynLabel)
        self.p2pDynLineEdit = QtWidgets.QLineEdit(self)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.p2pDynLineEdit.sizePolicy().hasHeightForWidth())
        self.p2pDynLineEdit.setSizePolicy(sizePolicy)
        self.p2pDynLineEdit.setMinimumSize(QtCore.QSize(20, 24))
        font.setBold(False)
        self.p2pDynLineEdit.setFont(font)
        self.p2pDynLineEdit.setObjectName("p2pDynLineEdit")
        self.airgapDynamicHLayout.addWidget(self.p2pDynLineEdit)
        self.mmDnyLabelOne = QtWidgets.QLabel(self)
        self.mmDnyLabelOne.setFont(font)
        self.mmDnyLabelOne.setObjectName("mmDnyLabelOne")
        self.airgapDynamicHLayout.addWidget(self.mmDnyLabelOne)
        spacerItem = QtWidgets.QSpacerItem(10, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        self.airgapDynamicHLayout.addItem(spacerItem)
        self.sizeDynLabel = QtWidgets.QLabel(self)
        font.setBold(True)
        font.setWeight(75)
        self.sizeDynLabel.setFont(font)
        self.sizeDynLabel.setObjectName("sizeDynLabel")
        self.airgapDynamicHLayout.addWidget(self.sizeDynLabel)
        self.sizeDynLineEdit = QtWidgets.QLineEdit(self)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.sizeDynLineEdit.sizePolicy().hasHeightForWidth())
        self.sizeDynLineEdit.setSizePolicy(sizePolicy)
        self.sizeDynLineEdit.setMinimumSize(QtCore.QSize(20, 24))
        font.setBold(False)
        self.sizeDynLineEdit.setFont(font)
        self.sizeDynLineEdit.setObjectName("sizeDynLineEdit")
        self.airgapDynamicHLayout.addWidget(self.sizeDynLineEdit)
        self.mmDynLabelTwo = QtWidgets.QLabel(self)
        self.mmDynLabelTwo.setFont(font)
        self.mmDynLabelTwo.setLineWidth(0)
        self.mmDynLabelTwo.setWordWrap(False)
        self.mmDynLabelTwo.setObjectName("mmDynLabelTwo")
        self.airgapDynamicHLayout.addWidget(self.mmDynLabelTwo)
        self.index = index
        self.retranslateUi()

    def retranslateUi(self):
        _translate = QtCore.QCoreApplication.translate
        self.countDynLabel.setText(_translate("Form", f"{self.index}."))
        self.p2pDynLabel.setText(_translate("Form", "<html><head/><body><p><span style=\" color:#00e3e3;\">P2P</span></p></body></html>"))
        self.mmDnyLabelOne.setText(_translate("Form", "mm"))
        self.sizeDynLabel.setText(_translate("Form", "<html><head/><body><p><span style=\" color:#00e3e3;\">size</span></p></body></html>"))
        self.mmDynLabelTwo.setText(_translate("Form", "<html><head/><body><p><span style=\" color:#000000;\">mm</span></p></body></html>"))


class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        uic.loadUi('FemmtGuiModel.ui', self)
        _translate = QtCore.QCoreApplication.translate
        # self.setWindowIcon(QIcon('Images\\logo.png'))
        self.setWindowTitle(_translate("MainWindow", "Femmt ToolBox"))
        pixmap = QPixmap('ferriteCore.png')
        self.coreImageLabel.setPixmap(pixmap)
        self.imageBoxImageLabel.setPixmap(pixmap)
        self.initialize_controls()

        """slots and events"""
        self.layersNumComboBox.currentTextChanged.connect(self.add_combo_boxes)
        self.windingNumComboBox.currentTextChanged.connect(self.add_winding_type_widgets)
        self.gapsNumComboBox.currentTextChanged.connect(self.add_air_gap_widgets)
        # self.aboutToQuit.connect(self.closeEvent)

    def initialize_controls(self):
        layers_count = ['', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10']
        winding_num_count = ['', '1', '2', '3']
        gap_num_count = ['', '1', '2', '3', '4']
        for option in layers_count:
            self.layersNumComboBox.addItem(option)
        for option in winding_num_count:
            self.windingNumComboBox.addItem(option)
        for option in gap_num_count:
            self.gapsNumComboBox.addItem(option)

    def add_combo_boxes(self, selection_count: str):
        if self.horizontalLayout_3.count():
            clear_layout(self.horizontalLayout_3)
            self.adjustSize()
        if not selection_count == '':
            for index in range(int(selection_count)):
                windingDropDown = QtWidgets.QComboBox(self.windingInfoBox)
                sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
                sizePolicy.setHorizontalStretch(0)
                sizePolicy.setVerticalStretch(0)
                sizePolicy.setHeightForWidth(windingDropDown.sizePolicy().hasHeightForWidth())
                windingDropDown.setSizePolicy(sizePolicy)
                windingDropDown.setMinimumSize(QtCore.QSize(20, 30))
                windingDropDown.setMaximumSize(QtCore.QSize(16777215, 16777215))
                windingDropDown.setStyleSheet("border : 1px solid black; border-top-left-radius : 5px; border-top-right-radius : 5px; border-bottom-left-radius:5px; border-bottom-right-radius : 5px;")
                windingDropDown.setFrame(True)
                windingDropDown.setObjectName("windingDropDown"+f'_{index}')
                self.horizontalLayout_3.addWidget(windingDropDown)
        self.adjustSize()

    def add_air_gap_widgets(self, selection_count: str):
        if self.airgapDynamicLayout.count():
            clear_layout(self.airgapDynamicLayout)
            self.adjustSize()
        if not selection_count == '':
            widget_length = int(selection_count)
            spacer_column = 1 if widget_length > 1 else 0
            for index in range(widget_length):
                air_gap_widget = AirGapControls(index+1)
                widget_position = divmod(index, 2)
                self.airgapDynamicLayout.addWidget(air_gap_widget, widget_position[0], widget_position[1]+spacer_column)
            if spacer_column:  # if there are more than one item then adding a spacer in between position (2*3 matrix with one 2nd column for spacer)
                spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
                self.airgapDynamicLayout.addItem(spacerItem, 0, 1)
        self.adjustSize()

    def add_winding_type_widgets(self, selection_count: str):
        if self.windingDynamicLayout.count():
            clear_layout(self.windingDynamicLayout)
            self.adjustSize()
        if not selection_count == '':
            for index in range(int(selection_count)):
                winding_type_widget = WindingControls(index+1)
                self.windingDynamicLayout.addWidget(winding_type_widget)
        self.adjustSize()


def clear_layout(layout):
    while layout.count():
        child = layout.takeAt(0)
        if isinstance(child, QtWidgets.QSpacerItem):
            layout.removeItem(child)   # for spacer item
        elif child.widget() or child:   # child for comboBoxType and child.widget() for custom class types
            child.widget().deleteLater()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    mainWindow = MainWindow()
    mainWindow.show()
    try:
        sys.exit(app.exec_())
    except SystemExit:
        print('Closing Window....')