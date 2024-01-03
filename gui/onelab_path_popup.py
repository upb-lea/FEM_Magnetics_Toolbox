"""Class and methods for the separate GUI window to ask for the ONELAB filepath at first run."""
from PyQt5 import uic
from PyQt5.QtWidgets import *
import sys
import os

class OnelabPathDialog(QDialog):
    """Class to open a separate GUI window to ask for the ONELAB filepath."""

    def __init__(self):
        super(QDialog, self).__init__()
        ui_file_path = os.path.join(os.path.dirname(__file__), "onelab_path_popup.ui")
        uic.loadUi(ui_file_path, self)

        self.directory = os.path.join(os.path.dirname(__file__))
        self.onelab_path_box.setText(self.directory)

        self.browse_button.clicked.connect(self.clicked_browse_button)

    def clicked_browse_button(self):
        """Select a directory for the onlab filepath inside the GUI."""
        directory = str(QFileDialog.getExistingDirectory(self, "Select Directory", directory=self.directory))
        if os.path.isdir(directory):
            self.directory = directory
            self.onelab_path_box.setText(directory)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = OnelabPathDialog()
    window.open()
    app.exec_()
