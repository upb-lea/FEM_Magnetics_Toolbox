"""Graphical user interface (GUI) for FEMMT."""
import sys
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QVBoxLayout
from PyQt5 import QtCore, uic, QtWidgets
from PyQt5.QtGui import QPixmap, QDoubleValidator, QIntValidator
import femmt as fmt
import json
import os
from typing import List
import PIL
import webbrowser

import materialdatabase as mdb
import matplotlib.pyplot as plt

from gui.onelab_path_popup import OnelabPathDialog

from femmt.examples.inductor_optimization import AutomatedDesign
from femmt.examples.inductor_optimization import load_fem_simulation_results, filter_after_fem

import mplcursors

database = mdb.MaterialDatabase()

float_validator = QDoubleValidator()
int_validator = QIntValidator()


def comma_str_to_point_float(input_str: str) -> float:
    """
    Workaround to convert a comma (depends on the user input) in a point (needed for python).

    :param input_str: input string
    :type input_str: str
    :return: converts comma to point and returns a float
    :rtype: float
    """
    output = ""
    for letter in input_str:
        if letter != ",":
            output += letter
        else:
            output += "."

    if output == "":
        output = 0

    return float(output)


class MatplotlibWidget(QWidget):
    """MatplotlibWidget class which inherits from QWidget and is used to implement a Matplotlib figure inside a QWidget."""

    def __init__(self, parent=None):
        super(MatplotlibWidget, self).__init__(parent)
        self.figure = Figure()
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.axis = self.figure.add_subplot(111)
        self.layout = QVBoxLayout(self)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        self.layout.addWidget(self.canvas)
        self.layout.addWidget(self.toolbar)
        self.divider = make_axes_locatable(self.axis)
        self.axis_cm = self.divider.append_axes("right", size="3%", pad=0.03)
        self.sm = plt.cm.ScalarMappable(cmap=cm.inferno)
        self.figure.colorbar(mappable=self.sm, cax=self.axis_cm)

class MainWindow(QMainWindow):
    """Global variable declaration."""

    litz_conductor_r = 0
    FEM_data_matrix = 0
    winding_factor = 0
    mc1 = 0
    freq = 0
    i_max = 0
    litz_strand_n = []
    litz_strand_r = []
    param = []
    no_of_turns = 0
    ad = 0

    def __init__(self, parent=None):

        self.default_gui_working_directory = os.path.join(os.path.dirname(__file__), "GUI_working_directory")
        if not os.path.exists(self.default_gui_working_directory):
            os.mkdir(self.default_gui_working_directory)

        super(MainWindow, self).__init__(parent)
        self.md_simulation_type_comboBox = None
        self.aut_simulation_type_comboBox = None
        ui_file_path = os.path.join(os.path.dirname(__file__), "femmt_gui.ui")
        uic.loadUi(ui_file_path, self)
        _translate = QtCore.QCoreApplication.translate
        # self.setWindowIcon(QIcon('Images\\logo.png'))
        self.setWindowTitle(_translate("MainWindow", "FEM Magnetics Toolbox"))
        pixmap = QPixmap('ferriteCore.png')
        # self.coreImageLabel.setPixmap(pixmap)
        # self.imageBoxImageLabel.setPixmap(pixmap)
        self.translation_dict = {
            # key: Used in FEMMT code
            # value: Used in GUI
            "inductor": "inductor",
            "transformer": "transformer",
            "integrated transformer": "transformer with integrated stray path",
            "litz": "Litz Wire",
            "solid": "Solid Round Wire",
            "implicit_litz_radius": 'Litz Radius',
            "implicit_ff": "Fill Factor",
            "implicit_strands_number": "Strands count",
            "center": "Center Placement",
            "random": "Random Placement",
            "percent": "Percent",
            "manually": "Manual Placement",
            "hexa": "Hexadimensional",
            "square": "Square",
            "+-10": "+/- 10%",
            "+-20": "+/- 20%",
            "excel": "MS Excel"
        }
        "******** menu bar *********"

        # ## Actions for Help
        self.action_contribute.triggered.connect(self.webbrowser_contribute)
        self.action_documentation.triggered.connect(self.webbrowser_documentation)
        self.action_report_bug.triggered.connect(self.webbrowser_bugreport)

        "******* Manual Design *********"

        "Signals in Definition Tab"
        # predefine working directory
        self.md_working_directory_lineEdit.setText(self.default_gui_working_directory)
        self.FEM_sim_working_dir_LineEdit.setText(self.default_gui_working_directory)
        self.aut_load_design_directoryname_lineEdit.setText(self.default_gui_working_directory)

        # simulation
        self.md_simulation_type_comboBox.currentTextChanged.connect(self.md_change_simulation_type)
        # core
        self.md_core_geometry_comboBox.currentTextChanged.connect(self.md_set_core_geometry_from_database)
        # windings
        self.md_winding1_type_comboBox.currentTextChanged.connect(self.md_winding1_change_wire_type)
        self.md_winding2_type_comboBox.currentTextChanged.connect(self.md_winding2_change_wire_type)
        self.md_winding1_litz_material_comboBox.currentTextChanged.connect(
            self.md_winding1_set_litz_parameters_from_litz_database)
        self.md_winding2_litz_material_comboBox.currentTextChanged.connect(
            self.md_winding2_set_litz_parameters_from_litz_database)
        self.md_winding1_implicit_litz_comboBox.currentTextChanged.connect(self.md_winding1_change_litz_implicit)
        self.md_winding2_implicit_litz_comboBox.currentTextChanged.connect(self.md_winding2_change_litz_implicit)

        # air gaps
        self.md_air_gap_count_comboBox.currentTextChanged.connect(self.md_change_air_gap_count)
        self.md_air_gap_placement_method_comboBox.currentTextChanged.connect(self.md_change_air_gap_placement)

        # visualization
        self.md_gmsh_visualisation_QPushButton.clicked.connect(self.md_gmsh_pre_visualisation)

        # Results
        self.Inductance_value_pushButton.clicked.connect(self.inductancecalc)

        "checkbox for Insulation mode"
        # needed for enabeling the insulation
        self.enable_insulation_checkbox.stateChanged.connect(self.on_insulation_state_changed)
        # Initialize flag_insulation based on the checkbox's initial state
        self.flag_insulation = self.enable_insulation_checkbox.isChecked()

        "Signals in Excitation Tab"
        self.md_dc_checkBox.stateChanged.connect(self.md_dc_enable)
        self.md_fk1_checkBox.stateChanged.connect(self.md_change_frequencies_1)
        self.md_fk2_checkBox.stateChanged.connect(self.md_change_frequencies_2)
        self.md_fk3_checkBox.stateChanged.connect(self.md_change_frequencies_3)
        self.md_fk4_checkBox.stateChanged.connect(self.md_change_frequencies_4)
        self.md_fk5_checkBox.stateChanged.connect(self.md_change_frequencies_5)
        self.md_fk6_checkBox.stateChanged.connect(self.md_change_frequencies_6)
        self.md_fk7_checkBox.stateChanged.connect(self.md_change_frequencies_7)
        self.md_fk8_checkBox.stateChanged.connect(self.md_change_frequencies_8)
        self.md_excitation_update_graph_Button.clicked.connect(self.md_redraw_input_signals)

        "Signals in Simulation Tab"
        self.md_simulation_QPushButton.clicked.connect(self.md_action_run_simulation)

        # initialize checkboxes. Need to be set to True and then to False, so represent a change and the signal
        # triggers the action to uncheck all boxes
        self.md_dc_checkBox.setChecked(True)
        self.md_fk1_checkBox.setChecked(True)
        self.md_fk2_checkBox.setChecked(True)
        self.md_fk3_checkBox.setChecked(True)
        self.md_fk4_checkBox.setChecked(True)
        self.md_fk5_checkBox.setChecked(True)
        self.md_fk6_checkBox.setChecked(True)
        self.md_fk7_checkBox.setChecked(True)
        self.md_fk8_checkBox.setChecked(True)
        # self.md_fk1_checkBox.setChecked(False)
        self.md_dc_checkBox.setChecked(False)
        self.md_fk2_checkBox.setChecked(False)
        self.md_fk3_checkBox.setChecked(False)
        self.md_fk4_checkBox.setChecked(False)
        self.md_fk5_checkBox.setChecked(False)
        self.md_fk6_checkBox.setChecked(False)
        self.md_fk7_checkBox.setChecked(False)
        self.md_fk8_checkBox.setChecked(False)

        "Init controls"
        self.md_initialize_controls()

        "Set validators in Definition tab"
        self.md_core_width_lineEdit.setValidator(float_validator)
        self.md_window_height_lineEdit.setValidator(float_validator)
        self.md_window_width_lineEdit.setValidator(float_validator)
        self.md_winding1_radius_lineEdit.setValidator(float_validator)
        self.md_winding1_strands_lineEdit.setValidator(int_validator)
        self.md_winding1_fill_factor_lineEdit.setValidator(float_validator)
        self.md_winding1_strand_radius_lineEdit.setValidator(float_validator)
        self.md_winding2_radius_lineEdit.setValidator(float_validator)
        self.md_winding2_strands_lineEdit.setValidator(int_validator)
        self.md_winding2_fill_factor_lineEdit.setValidator(float_validator)
        self.md_winding2_strand_radius_lineEdit.setValidator(float_validator)
        self.md_winding1_turns_lineEdit.setValidator(int_validator)
        self.md_winding2_turns_lineEdit.setValidator(int_validator)
        self.md_air_gap_1_length_lineEdit.setValidator(float_validator)
        self.md_air_gap_2_length_lineEdit.setValidator(float_validator)
        self.md_air_gap_3_length_lineEdit.setValidator(float_validator)
        self.md_air_gap_4_length_lineEdit.setValidator(float_validator)
        self.md_air_gap_5_length_lineEdit.setValidator(float_validator)
        self.md_air_gap_1_position_lineEdit.setValidator(float_validator)
        self.md_air_gap_2_position_lineEdit.setValidator(float_validator)
        self.md_air_gap_3_position_lineEdit.setValidator(float_validator)
        self.md_air_gap_4_position_lineEdit.setValidator(float_validator)
        self.md_air_gap_5_position_lineEdit.setValidator(float_validator)
        self.md_isolation_p2p_lineEdit.setValidator(float_validator)
        self.md_isolation_p2s_lineEdit.setValidator(float_validator)
        self.md_isolation_s2s_lineEdit.setValidator(float_validator)
        # change 4: s2p
        self.md_isolation_s2p_lineEdit.setValidator(float_validator)
        self.md_isolation_core2cond_top_lineEdit.setValidator(float_validator)
        self.md_isolation_core2cond_bot_lineEdit.setValidator(float_validator)
        self.md_isolation_core2cond_inner_lineEdit.setValidator(float_validator)
        self.md_isolation_core2cond_outer_lineEdit.setValidator(float_validator)
        "Set Validators in Excitation Tab"
        self.md_winding1_idc_lineEdit.setValidator(float_validator)
        self.md_winding1_ik1_lineEdit.setValidator(float_validator)
        self.md_winding1_ik2_lineEdit.setValidator(float_validator)
        self.md_winding1_ik3_lineEdit.setValidator(float_validator)
        self.md_winding1_ik4_lineEdit.setValidator(float_validator)
        self.md_winding1_ik5_lineEdit.setValidator(float_validator)
        self.md_winding1_ik6_lineEdit.setValidator(float_validator)
        self.md_winding1_ik7_lineEdit.setValidator(float_validator)
        self.md_winding1_ik8_lineEdit.setValidator(float_validator)
        self.md_winding1_pk1_lineEdit.setValidator(float_validator)
        self.md_winding1_pk2_lineEdit.setValidator(float_validator)
        self.md_winding1_pk3_lineEdit.setValidator(float_validator)
        self.md_winding1_pk4_lineEdit.setValidator(float_validator)
        self.md_winding1_pk5_lineEdit.setValidator(float_validator)
        self.md_winding1_pk6_lineEdit.setValidator(float_validator)
        self.md_winding1_pk7_lineEdit.setValidator(float_validator)
        self.md_winding1_pk8_lineEdit.setValidator(float_validator)

        self.md_winding2_idc_lineEdit.setValidator(float_validator)
        self.md_winding2_ik1_lineEdit.setValidator(float_validator)
        self.md_winding2_ik2_lineEdit.setValidator(float_validator)
        self.md_winding2_ik3_lineEdit.setValidator(float_validator)
        self.md_winding2_ik4_lineEdit.setValidator(float_validator)
        self.md_winding2_ik5_lineEdit.setValidator(float_validator)
        self.md_winding2_ik6_lineEdit.setValidator(float_validator)
        self.md_winding2_ik7_lineEdit.setValidator(float_validator)
        self.md_winding2_ik8_lineEdit.setValidator(float_validator)
        self.md_winding2_pk1_lineEdit.setValidator(float_validator)
        self.md_winding2_pk2_lineEdit.setValidator(float_validator)
        self.md_winding2_pk3_lineEdit.setValidator(float_validator)
        self.md_winding2_pk4_lineEdit.setValidator(float_validator)
        self.md_winding2_pk5_lineEdit.setValidator(float_validator)
        self.md_winding2_pk6_lineEdit.setValidator(float_validator)
        self.md_winding2_pk7_lineEdit.setValidator(float_validator)
        self.md_winding2_pk8_lineEdit.setValidator(float_validator)

        self.md_base_frequency_lineEdit.setValidator(int_validator)

        "Set Tool Tips in Definition tab"
        self.md_core_geometry_comboBox.setToolTip(
            "Chose a core geometry from the database. Chose 'Manual' to insert any parameters.")
        self.md_winding1_litz_material_comboBox.setToolTip(
            "Chose a litz from the database. Chose 'Manual' to insert any parameters")
        self.md_winding1_implicit_litz_comboBox.setToolTip(
            "To describe a strand, 3 arguments are sufficient. Select here which of the arguments should not be entered.")
        self.md_winding2_implicit_litz_comboBox.setToolTip(
            "To describe a strand, 3 arguments are sufficient. Select here which of the arguments should not be entered.")
        self.md_winding2_litz_material_comboBox.setToolTip(
            "Chose a litz from the database. Chose 'Manual' to insert any parameters")

        "Set Tool Tips in exitation tab"
        self.md_winding1_idc_lineEdit.setToolTip("DC Current")
        self.md_winding1_ik1_lineEdit.setToolTip("Amplitude base frequency")
        self.md_winding1_ik2_lineEdit.setToolTip("Amplitude 2 * base frequency")
        self.md_winding1_ik3_lineEdit.setToolTip("Amplitude 3 * base frequency")
        self.md_winding1_ik4_lineEdit.setToolTip("Amplitude 4 * base frequency")
        self.md_winding1_ik5_lineEdit.setToolTip("Amplitude 5 * base frequency")
        self.md_winding1_ik6_lineEdit.setToolTip("Amplitude 6 * base frequency")
        self.md_winding1_ik7_lineEdit.setToolTip("Amplitude 7 * base frequency")
        self.md_winding1_ik8_lineEdit.setToolTip("Amplitude 8 * base frequency")
        self.md_winding1_pk1_lineEdit.setToolTip("Phase for base frequency")
        self.md_winding1_pk2_lineEdit.setToolTip("Phase for 2 * base frequency")
        self.md_winding1_pk3_lineEdit.setToolTip("Phase for 3 * base frequency")
        self.md_winding1_pk4_lineEdit.setToolTip("Phase for 4 * base frequency")
        self.md_winding1_pk5_lineEdit.setToolTip("Phase for 5 * base frequency")
        self.md_winding1_pk6_lineEdit.setToolTip("Phase for 6 * base frequency")
        self.md_winding1_pk7_lineEdit.setToolTip("Phase for 7 * base frequency")
        self.md_winding1_pk8_lineEdit.setToolTip("Phase for 8 * base frequency")

        self.md_winding2_idc_lineEdit.setToolTip("DC Current")
        self.md_winding2_ik1_lineEdit.setToolTip("Amplitude base frequency")
        self.md_winding2_ik2_lineEdit.setToolTip("Amplitude 2 * base frequency")
        self.md_winding2_ik3_lineEdit.setToolTip("Amplitude 3 * base frequency")
        self.md_winding2_ik4_lineEdit.setToolTip("Amplitude 4 * base frequency")
        self.md_winding2_ik5_lineEdit.setToolTip("Amplitude 5 * base frequency")
        self.md_winding2_ik6_lineEdit.setToolTip("Amplitude 6 * base frequency")
        self.md_winding2_ik7_lineEdit.setToolTip("Amplitude 7 * base frequency")
        self.md_winding2_ik8_lineEdit.setToolTip("Amplitude 8 * base frequency")
        self.md_winding2_pk1_lineEdit.setToolTip("Phase for base frequency")
        self.md_winding2_pk2_lineEdit.setToolTip("Phase for 2 * base frequency")
        self.md_winding2_pk3_lineEdit.setToolTip("Phase for 3 * base frequency")
        self.md_winding2_pk4_lineEdit.setToolTip("Phase for 4 * base frequency")
        self.md_winding2_pk5_lineEdit.setToolTip("Phase for 5 * base frequency")
        self.md_winding2_pk6_lineEdit.setToolTip("Phase for 6 * base frequency")
        self.md_winding2_pk7_lineEdit.setToolTip("Phase for 7 * base frequency")
        self.md_winding2_pk8_lineEdit.setToolTip("Phase for 8 * base frequency")

        self.md_dc_checkBox.setToolTip("Enable/Disable DC current")
        self.md_fk1_checkBox.setToolTip("Enable/Disable base frequency")
        self.md_fk2_checkBox.setToolTip("Enable/Disable 2 * base frequency")
        self.md_fk3_checkBox.setToolTip("Enable/Disable 3 * base frequency")
        self.md_fk4_checkBox.setToolTip("Enable/Disable 4 * base frequency")
        self.md_fk5_checkBox.setToolTip("Enable/Disable 5 * base frequency")
        self.md_fk6_checkBox.setToolTip("Enable/Disable 6 * base frequency")
        self.md_fk7_checkBox.setToolTip("Enable/Disable 7 * base frequency")
        self.md_fk8_checkBox.setToolTip("Enable/Disable 8 * base frequency")

        "Signals in Thermal simulation Tab"
        self.md_therm_simulation_QPushButton.clicked.connect(self.therm_simulation)

        "******* Automated Design *********"

        "Adding options"
        self.aut_initialize_controls()

        "Signals in Definition Tab"
        if self.aut_simulation_type_comboBox.currentText() == self.translation_dict['inductor']:
            self.aut_winding2_enable(False)
        self.aut_simulation_type_comboBox.currentTextChanged.connect(self.aut_change_simulation_type)
        self.aut_core_material_data_listWidget.addItem("N95")
        self.aut_core_material_data_listWidget.addItem("N87")
        self.aut_core_material_data_listWidget.addItem("N49")

        self.aut_litz_data_listWidget.addItem("1.5x105x0.1")
        self.aut_litz_data_listWidget.addItem("1.4x200x0.071")
        self.aut_litz_data_listWidget.addItem("2.0x405x0.071")
        self.aut_litz_data_listWidget.addItem("2.0x800x0.05")

        self.aut_litz2_data_listWidget.addItem("1.5x105x0.1")
        self.aut_litz2_data_listWidget.addItem("1.4x200x0.071")
        self.aut_litz2_data_listWidget.addItem("2.0x405x0.071")
        self.aut_litz2_data_listWidget.addItem("2.0x800x0.05")

        self.aut_airgap_type_listWidget.addItem("Edge distributed")
        self.aut_airgap_type_listWidget.addItem("Centre distributed")

        self.aut_winding1_type_comboBox.currentTextChanged.connect(self.aut_winding1_change_wire_type)

        self.aut_core_geo_add_pushButton.clicked.connect(self.oncgeoMultipleClicked)
        self.aut_core_geo_manual_add_pushButton.clicked.connect(self.oncgeomanualMultipleClicked)
        self.aut_select_all_core_geo_pushButton.clicked.connect(self.cgeoselectall)
        self.aut_core_geometry_listWidget.itemDoubleClicked.connect(self.oncgeoClicked)
        self.aut_core_geo_basket_clear_all_pushbutton.clicked.connect(self.oncgeoClearallClicked)
        self.aut_core_geo_basket_clear_pushbutton.clicked.connect(self.oncgeoClearClicked)

        self.aut_core_material_add_pushButton.clicked.connect(self.oncmatMultipleClicked)
        self.aut_select_all_core_mat_pushButton.clicked.connect(self.cmatselectall)
        self.aut_core_material_data_listWidget.itemDoubleClicked.connect(self.oncmatClicked)
        self.aut_core_mat_basket_clear_all_pushbutton.clicked.connect(self.oncmatClearallClicked)
        self.aut_core_mat_basket_clear_pushbutton.clicked.connect(self.oncmatClearClicked)

        self.aut_add_litz_pushButton.clicked.connect(self.onl1MultipleClicked)
        self.aut_select_all_litz_pushButton.clicked.connect(self.litz1selectall)
        self.aut_litz_data_listWidget.itemDoubleClicked.connect(self.onl1Clicked)
        self.aut_litz_basket_clear_all_pushbutton.clicked.connect(self.onl1ClearallClicked)
        self.aut_litz_basket_clear_pushbutton.clicked.connect(self.onl1ClearClicked)

        self.aut_add_litz2_pushButton.clicked.connect(self.onl2MultipleClicked)
        self.aut_select_all_litz2_pushButton.clicked.connect(self.litz2selectall)
        self.aut_litz2_data_listWidget.itemDoubleClicked.connect(self.onl2Clicked)
        self.aut_litz2_basket_clear_all_pushbutton.clicked.connect(self.onl2ClearallClicked)
        self.aut_litz2_basket_clear_pushbutton.clicked.connect(self.onl2ClearClicked)

        self.aut_add_air_gap_type_pushButton.clicked.connect(self.onairgaptypeMultipleClicked)
        self.aut_select_all_airgap_type_pushButton.clicked.connect(self.airgaptypeselectall)
        self.aut_airgap_type_listWidget.itemDoubleClicked.connect(self.onairgaptypeClicked)
        self.aut_air_gap_type_basket_clear_all_pushbutton.clicked.connect(self.onairgaptypeClearallClicked)
        self.aut_air_gap_type_basket_clear_pushbutton.clicked.connect(self.onairgaptypeClearClicked)

        self.aut_winding1_implicit_litz_comboBox.currentTextChanged.connect(self.aut_winding1_change_litz_implicit)
        self.aut_winding1_change_litz_implicit(self.aut_winding1_implicit_litz_comboBox.currentText())

        "Set Validators in Definition Tab"
        self.aut_min_core_width_lineEdit.setValidator(float_validator)
        self.aut_max_core_width_lineEdit.setValidator(float_validator)
        self.aut_step_core_width_lineEdit.setValidator(float_validator)

        self.aut_min_window_height_lineEdit.setValidator(float_validator)
        self.aut_max_window_height_lineEdit.setValidator(float_validator)
        self.aut_step_window_height_lineEdit.setValidator(float_validator)

        self.aut_min_window_width_lineEdit.setValidator(float_validator)
        self.aut_max_window_width_lineEdit.setValidator(float_validator)
        self.aut_step_window_width_lineEdit.setValidator(float_validator)

        self.aut_winding1_radius_lineEdit.setValidator(float_validator)
        self.aut_winding1_strands_lineEdit.setValidator(float_validator)
        self.aut_winding1_fill_factor_lineEdit.setValidator(float_validator)
        self.aut_winding1_strand_radius_lineEdit.setValidator(float_validator)

        self.aut_winding2_radius_lineEdit.setValidator(float_validator)
        self.aut_winding2_strands_lineEdit.setValidator(float_validator)
        self.aut_winding2_fill_factor_lineEdit.setValidator(float_validator)
        self.aut_winding2_strand_radius_lineEdit.setValidator(float_validator)

        self.aut_min_winding1_turns_lineEdit.setValidator(float_validator)
        self.aut_max_winding1_turns_lineEdit.setValidator(float_validator)

        self.aut_min_winding2_turns_lineEdit.setValidator(float_validator)
        self.aut_max_winding2_turns_lineEdit.setValidator(float_validator)
        self.aut_step_winding2_turns_lineEdit.setValidator(float_validator)

        self.aut_min_air_gap_count_lineEdit.setValidator(float_validator)
        self.aut_max_air_gap_count_lineEdit.setValidator(float_validator)

        self.aut_air_gap_length_min_lineEdit.setValidator(float_validator)
        self.aut_air_gap_length_max_lineEdit.setValidator(float_validator)
        self.aut_air_gap_length_step_lineEdit.setValidator(float_validator)
        self.aut_air_gap_position_min_lineEdit.setValidator(float_validator)
        self.aut_air_gap_position_step_lineEdit.setValidator(float_validator)
        self.aut_air_gap_position_step_lineEdit.setValidator(float_validator)

        self.aut_isolation_p2p_lineEdit.setValidator(float_validator)
        self.aut_isolation_p2s_lineEdit.setValidator(float_validator)
        self.aut_isolation_s2s_lineEdit.setValidator(float_validator)

        self.aut_isolation_core2cond_top_lineEdit.setValidator(float_validator)
        self.aut_isolation_core2cond_bot_lineEdit.setValidator(float_validator)
        self.aut_isolation_core2cond_inner_lineEdit.setValidator(float_validator)
        self.aut_isolation_core2cond_outer_lineEdit.setValidator(float_validator)

        self.aut_temp_lineEdit.setValidator(float_validator)

        "Set Validators in Reluctance Models tab"
        self.aut_goal_inductance_val_lineEdit.setValidator(float_validator)
        self.aut_maximum_current_lineEdit.setValidator(float_validator)
        self.aut_switching_freq_lineEdit.setValidator(float_validator)
        self.aut_b_sat_lineEdit.setValidator(float_validator)
        self.aut_hysterisis_loss_options_lineEdit.setValidator(float_validator)

        "Signals in Reluctance Models tab"

        self.aut_simulate_pushButton.clicked.connect(self.automated_design_func_config)
        self.matplotlib_widget_aut_tab2 = MatplotlibWidget()

        "Signals in FEM Simulations tab"

        self.aut_pos_mod_sim_pushButton.clicked.connect(self.automated_design_fem_sim_config)
        self.matplotlib_widget_aut_tab3 = MatplotlibWidget()

        "Signals in Load(Results) tab"
        self.aut_load_design_pushButton.clicked.connect(self.load_designs_config)
        self.matplotlib_widget_aut_tab4 = MatplotlibWidget()

        "******* Database Section *********"

        self.dat_update_preview_pushbutton1.clicked.connect(self.datupdateraph1_config)
        self.dat_update_preview_pushbutton2.clicked.connect(self.datupdateraph2_config)
        self.dat_update_preview_pushbutton3.clicked.connect(self.datupdateraph3_config)

        self.matplotlib_widget_datmm1 = MatplotlibWidget()
        self.matplotlib_widget_datmm2 = MatplotlibWidget()

        self.matplotlib_widget_datdd1 = MatplotlibWidget()
        self.matplotlib_widget_datdd2 = MatplotlibWidget()
        self.matplotlib_widget_datdd3 = MatplotlibWidget()
        self.matplotlib_widget_datdd4 = MatplotlibWidget()

        self.matplotlib_widget_datdm = MatplotlibWidget()

        self.dat_core_material1_comboBox.currentTextChanged.connect(self.tempfluxinput1)
        self.dat_core_material2_comboBox.currentTextChanged.connect(self.tempfluxinput2)
        self.dat_core_material3_comboBox.currentTextChanged.connect(self.tempfluxinput3)
        self.dat_core_material4_comboBox.currentTextChanged.connect(self.tempfluxinput4)
        self.dat_core_material5_comboBox.currentTextChanged.connect(self.tempfluxinput5)

        self.dat_core_material1_comboBox_2.currentTextChanged.connect(self.test_setup_name1)
        self.dat_core_material2_comboBox_2.currentTextChanged.connect(self.test_setup_name2)
        self.dat_core_material3_comboBox_2.currentTextChanged.connect(self.test_setup_name3)
        self.dat_core_material4_comboBox_2.currentTextChanged.connect(self.test_setup_name4)
        self.dat_core_material5_comboBox_2.currentTextChanged.connect(self.test_setup_name5)

        self.dat_core_material1_comboBox_2.currentTextChanged.connect(self.tempfreqinput1)
        self.dat_core_material2_comboBox_2.currentTextChanged.connect(self.tempfreqinput2)
        self.dat_core_material3_comboBox_2.currentTextChanged.connect(self.tempfreqinput3)
        self.dat_core_material4_comboBox_2.currentTextChanged.connect(self.tempfreqinput4)
        self.dat_core_material5_comboBox_2.currentTextChanged.connect(self.tempfreqinput5)

        self.test_name_1_comboBox.currentTextChanged.connect(self.tempfreqinput1)
        self.test_name_2_comboBox.currentTextChanged.connect(self.tempfreqinput2)
        self.test_name_3_comboBox.currentTextChanged.connect(self.tempfreqinput3)
        self.test_name_4_comboBox.currentTextChanged.connect(self.tempfreqinput4)
        self.test_name_5_comboBox.currentTextChanged.connect(self.tempfreqinput5)

        self.dat_core_material_comboBox.currentTextChanged.connect(self.temp_dat_input)
        self.dat_core_material_comboBox.currentTextChanged.connect(self.temp_meas_input)

    #  **************************** Menu bar ************************************************************  #
    # ## Help actions ###
    def string_to_material_enum(self, material_name):
        """
        Convert a material name string to its corresponding enum value.

        :param material_name: The name of the material to convert.
        :type material_name: str
        :returns: The enum value corresponding to the material name.
        :rtype: enum
        :raises ValueError: If the material name is not valid.
        """
        try:
            return fmt.Material[material_name]
        except KeyError as e:
            raise ValueError(f"{material_name} is not a valid material name") from e

    def webbrowser_contribute(self):
        """Open the web browser to the GitHub FEMMT repository contribution page."""
        webbrowser.open('https://github.com/upb-lea/FEM_Magnetics_Toolbox')

    def webbrowser_bugreport(self):
        """Open the web browser to the GitHub FEMMT repository issue page."""
        webbrowser.open('https://github.com/upb-lea/FEM_Magnetics_Toolbox/issues')

    def webbrowser_documentation(self):
        """Open the web brower to the FEMMT documentation."""
        webbrowser.open('https://upb-lea.github.io/FEM_Magnetics_Toolbox/main/intro.html')

    #  **************************** Automated design tab ************************************************************  #

    def plot_volume_loss(self, data_matrix, matplotlib_widget):
        """
        Plot estimated normalised volume vs loss graph from reluctance model results.

        :param data_matrix: Matrix containing the design parameters
        :type data_matrix: array
        """
        matplotlib_widget.axis.set(xlabel="Volume / m\u00b3", ylabel="Loss / W", title=" Volume vs Loss")
        lines = matplotlib_widget.axis.plot(data_matrix[:, 30],
                                            data_matrix[:, 28], 'o')
        mplcursors.cursor(lines)
        matplotlib_widget.figure.tight_layout()
        matplotlib_widget.axis.grid()

    def plot_2d(self, matplotlib_widget, x_value: list, y_value: list, x_label: str, y_label: str, title: str,
                plot_color: str,
                z_value: list = None,
                z_label: str = None, inductance_value: list = None, annotations: list = None):
        """
        Visualize data in 2d plot with popover next to mouse position.

        :param x_value: Data points for x-axis
        :type x_value: list
        :param y_value: Data points for y-axis
        :type y_value: list
        :param z_value: Data points for z-axis
        :type z_value: list
        :param x_label: x-axis label
        :type x_label: str
        :param y_label: y-axis label
        :type y_label: str
        :param z_label: z-axis label
        :type z_label: str
        :param title: Title of the graph
        :type title: str
        :param inductance_value: Data points for inductance value corresponding to the (x, y, z): (Optional)
        :type inductance_value: list
        :param annotations: Annotations corresponding to the 3D points
        :type annotations: list
        :param plot_color: Color of the plot (the colors are based on 'fmt.colors_femmt_default')
        :type annotations: str
        """
        if annotations is None:
            names = [str(x) for x in list(range(len(x_value)))]
        else:
            temp_var = [int(x) for x in annotations]
            names = [str(x) for x in temp_var]

        if inductance_value is not None:
            l_label = 'L / H'

        if z_value is not None:
            z_value_str = [str(round(z, 3)) for z in z_value]

        if inductance_value is not None:
            l_value_str = [str(round(inductance, 6)) for inductance in inductance_value]

        x_value_str = [str(round(x, 6)) for x in x_value]
        y_value_str = [str(round(y, 3)) for y in y_value]

        if z_value is None:
            sc = matplotlib_widget.axis.scatter(x_value, y_value,
                                                c='#%02x%02x%02x' % fmt.colors_femmt_default[plot_color])

        else:
            sc = matplotlib_widget.axis.scatter(x_value, y_value, c=z_value, cmap=plot_color)
            cbar = matplotlib_widget.figure.colorbar(sc)
            cbar.ax.get_yaxis().labelpad = 15
            cbar.ax.set_ylabel(z_label, rotation=270)

        mplcursors.cursor(sc)
        matplotlib_widget.axis.set(xlabel=x_label, ylabel=y_label, title=title)
        annot = matplotlib_widget.axis.annotate("", xy=(0, 0), xytext=(20, 20), textcoords="offset points",
                                                bbox=dict(boxstyle="round", fc="w"),
                                                arrowprops=dict(arrowstyle="->"))
        annot.set_visible(False)
        matplotlib_widget.figure.tight_layout()
        matplotlib_widget.axis.grid()

        def update_annot(ind):
            """Create popover annotations in 2d plot."""
            pos = sc.get_offsets()[ind["ind"][0]]
            annot.xy = pos
            text = ""
            if z_label is None and inductance_value is None:
                text = "case: {}\n{}: {}\n{}:{}". \
                    format(" ".join([names[n] for n in ind["ind"]]),
                           x_label, " ".join([x_value_str[n] for n in ind["ind"]]),
                           y_label, " ".join([y_value_str[n] for n in ind["ind"]]))
            elif z_label is not None and inductance_value is None:
                text = "case: {}\n{}: {}\n{}:{}\n{}:{}". \
                    format(" ".join([names[n] for n in ind["ind"]]),
                           x_label, " ".join([x_value_str[n] for n in ind["ind"]]),
                           y_label, " ".join([y_value_str[n] for n in ind["ind"]]),
                           z_label, " ".join([z_value_str[n] for n in ind["ind"]]))
            elif z_label is None and inductance_value is not None:
                text = "case: {}\n{}: {}\n{}:{}\n{}:{}". \
                    format(" ".join([names[n] for n in ind["ind"]]),
                           x_label, " ".join([x_value_str[n] for n in ind["ind"]]),
                           y_label, " ".join([y_value_str[n] for n in ind["ind"]]),
                           l_label, " ".join([l_value_str[n] for n in ind["ind"]]))
            else:
                text = "case: {}\n{}: {}\n{}:{}\n{}:{}\n{}:{}". \
                    format(" ".join([names[n] for n in ind["ind"]]),
                           x_label, " ".join([x_value_str[n] for n in ind["ind"]]),
                           y_label, " ".join([y_value_str[n] for n in ind["ind"]]),
                           z_label, " ".join([z_value_str[n] for n in ind["ind"]]),
                           l_label, " ".join([l_value_str[n] for n in ind["ind"]]))
            annot.set_text(text)
            annot.get_bbox_patch().set_alpha(0.4)

    def automated_design_func(self, matplotlib_widget):
        """Create to accept input parameters from the definitions tab, to create matrix with all input combinations.

        A call is being made to the reluctance model to filter out the cases for FEM simulation.

        :param matplotlib_widget: To plot volume vs loss in reluctance models tab
        """
        # ########################################   {DESIGN PARAMETERS}   #################################################
        goal_inductance = comma_str_to_point_float(
            self.aut_goal_inductance_val_lineEdit.text())  # Automated design-Reluctacne model-Goal Inductance
        self.trans_dict = {
            # key: Used in FEMMT code
            # value: Used in GUI
            "+/- 10%": "10",
            "Edge distributed": "1",
            "Centre distributed": "2"
        }
        L_tolerance_percent = int(self.trans_dict[self.aut_rel_tolerance_val_comboBox.currentText()])
        self.i_max = comma_str_to_point_float(self.aut_maximum_current_lineEdit.text())
        # Max current amplitude with assumption of sinusoidal current waveform
        percent_of_B_sat = int(self.aut_b_sat_lineEdit.text())  # Percent of B_sat allowed in the designed core

        percent_of_total_loss = int(self.aut_hysterisis_loss_options_lineEdit.text())
        # Percent of total_loss allowed in FEM simulation

        self.freq = comma_str_to_point_float(self.aut_switching_freq_lineEdit.text())
        mu_imag = 100
        Cu_sigma = 5.96 * 1e7  # copper conductivity (sigma) @ 20 degree celsius

        # Set core-geometry from core database or/and manual entry
        min_core_w = comma_str_to_point_float(self.aut_min_core_width_lineEdit.text())
        max_core_w = comma_str_to_point_float(self.aut_max_core_width_lineEdit.text())
        step_core_w = int(self.aut_step_core_width_lineEdit.text())
        min_window_h = comma_str_to_point_float(self.aut_min_window_height_lineEdit.text())
        max_window_h = comma_str_to_point_float(self.aut_max_window_height_lineEdit.text())
        step_window_h = int(self.aut_step_window_height_lineEdit.text())
        min_window_w = comma_str_to_point_float(self.aut_min_window_width_lineEdit.text())
        max_window_w = comma_str_to_point_float(self.aut_max_window_width_lineEdit.text())
        step_window_w = int(self.aut_step_window_width_lineEdit.text())

        manual_core_w = list(np.linspace(min_core_w, max_core_w, step_core_w))
        manual_window_h = list(np.linspace(min_window_h, max_window_h, step_window_h))
        manual_window_w = list(np.linspace(min_window_w, max_window_w, step_window_w))

        db_core_names = []  # "PQ 40/40", "PQ 40/30"
        for i in range(self.aut_core_geometry_basket_listWidget.count()):
            db_core_names.append(self.aut_core_geometry_basket_listWidget.item(i).text())

        core_db = fmt.core_database()
        db_core_w = [core_db[core_name]["core_inner_diameter"] for core_name in db_core_names]
        db_window_h = [core_db[core_name]["window_h"] for core_name in db_core_names]
        db_window_w = [core_db[core_name]["window_w"] for core_name in db_core_names]

        core_w_list = db_core_w + manual_core_w
        window_h_list = db_window_h + manual_window_h
        window_w_list = db_window_w + manual_window_w

        # Set winding settings (Solid and Litz winding type)
        solid_conductor_r = [comma_str_to_point_float(self.aut_winding1_radius_lineEdit.text())]

        litz_db = fmt.litz_database()
        # litz_names = ["1.5x105x0.1"]  # "1.5x105x0.1", "1.4x200x0.071"
        litz_names = []
        for i in range(self.aut_litz_basket_listWidget.count()):
            litz_names.append(self.aut_litz_basket_listWidget.item(i).text())
        self.litz_conductor_r = [litz_db[litz_name]["conductor_radii"] for litz_name in litz_names]
        self.litz_strand_r = [litz_db[litz_name]["strand_radii"] for litz_name in litz_names]
        self.litz_strand_n = [litz_db[litz_name]["strands_numbers"] for litz_name in litz_names]

        winding_scheme = self.aut_winding1_scheme_comboBox.currentText()

        min_conductor_r = min(self.litz_conductor_r + solid_conductor_r)

        # Set air-gap and core parameters7
        no_turns_min = int(self.aut_min_winding1_turns_lineEdit.text())
        no_turns_max = int(self.aut_max_winding1_turns_lineEdit.text())
        no_turns_step = (no_turns_max - no_turns_min) + 1
        no_airgaps_min = int(self.aut_min_air_gap_count_lineEdit.text())
        no_airgaps_max = int(self.aut_max_air_gap_count_lineEdit.text())
        airgap_h_min = comma_str_to_point_float(self.aut_air_gap_length_min_lineEdit.text())
        airgap_h_max = comma_str_to_point_float(self.aut_air_gap_length_max_lineEdit.text())
        airgap_h_step = int(self.aut_air_gap_length_step_lineEdit.text())
        airgap_pos_min = int(self.aut_air_gap_position_min_lineEdit.text())
        airgap_pos_max = int(self.aut_air_gap_position_max_lineEdit.text())
        airgap_pos_step = int(self.aut_air_gap_position_step_lineEdit.text())

        no_of_turns_float = list((np.linspace(no_turns_min, no_turns_max, no_turns_step)))
        no_of_turns = [int(item) for item in no_of_turns_float]
        n_air_gaps = [no_airgaps_min, no_airgaps_max]  # Set No. of air-gaps (n)
        # when the user chooses [i,j]; where i = j, the n_air_gaps will be just i or j
        if no_airgaps_min == no_airgaps_max:
            n_air_gaps = [no_airgaps_min]
        air_gap_height = list(np.linspace(airgap_h_min, airgap_h_max, airgap_h_step))  # Set air-gap length in metre (l)
        air_gap_position = list(np.linspace(airgap_pos_min, airgap_pos_max,
                                            airgap_pos_step))  # Set air-gap position in percent w.r.t. core window height

        # material_names = []
        # for i in range(self.aut_core_material_basket_listWidget.count()):
        #     material_names.append(self.aut_core_material_basket_listWidget.item(i).text())
        # # material_names = ["N95"]
        # #material_enum = fmt.Material(material_names)
        # mu_rel = [database.get_material_attribute(material_name=material_names, attribute="initial_permeability")
        #           for material_name in material_names]

        # selected materials as enums
        selected_materials = []
        # Iterate over the items in the material list widget
        for i in range(self.aut_core_material_data_listWidget.count()):
            item = self.aut_core_material_data_listWidget.item(i)
            # Check if the item is selected by the user
            if item.isSelected():
                # Convert the string of the selected item to its corresponding enum value
                material_enum = self.string_to_material_enum(item.text())
                selected_materials.append(material_enum)
        # empty list to hold the initial real permeability values
        mu_rels = []
        # Iterate over the selected materials as enums
        for material_name in selected_materials:
            try:
                # initial permeability attribute from the database
                mu_rel = database.get_material_attribute(material_name=material_name, attribute="initial_permeability")
                mu_rels.append(mu_rel)
            except KeyError as e:
                print(f"Failed to retrieve data for {material_name}: {e}")

        component = self.aut_simulation_type_comboBox.currentText()
        # Set two types of equally distributed air-gaps (used only for air-gaps more than 1):
        # Type 1: Equally distributed air-gaps including corner air-gaps (eg: air-gaps-position = [0, 50, 100])
        # Type 2: Equally distributed air-gaps excluding corner air-gaps (eg: air-gaps-position = [25, 50, 75])
        # 'Type1 = with corner air-gaps; 'Type2' = without corner air-gaps; 'Type0' = single air-gap
        mult_air_gap_type = [2]  # Type1-Edge, Type2: Centre

        matplotlib_widget.axis.clear()
        self.layout = QVBoxLayout(self.fem_vol_loss_plotwidget)
        self.layout.addWidget(matplotlib_widget)
        try:
            matplotlib_widget.axis_cm.remove()
        except:
            pass

        fem_directory = self.FEM_sim_working_dir_LineEdit.text()

        if self.aut_simulation_type_comboBox.currentText() == self.translation_dict['inductor']:
            self.ad = AutomatedDesign(working_directory=fem_directory,
                                      magnetic_component='inductor',
                                      target_inductance=goal_inductance,
                                      frequency=self.freq,
                                      target_inductance_percent_tolerance=L_tolerance_percent,
                                      winding_scheme=winding_scheme,
                                      peak_current=self.i_max,
                                      percent_of_flux_density_saturation=percent_of_B_sat,
                                      percent_of_total_loss=percent_of_total_loss,
                                      database_core_names=db_core_names,
                                      database_litz_names=litz_names,
                                      solid_conductor_r=solid_conductor_r,
                                      manual_core_inner_diameter=manual_core_w,
                                      manual_window_h=manual_window_h,
                                      manual_window_w=manual_window_w,
                                      no_of_turns=no_of_turns,
                                      n_air_gaps=n_air_gaps,
                                      air_gap_height=air_gap_height,
                                      air_gap_position=air_gap_position,
                                      core_material=selected_materials,
                                      mult_air_gap_type=['center_distributed'],
                                      top_core_insulation=comma_str_to_point_float(
                                          self.aut_isolation_core2cond_top_lineEdit.text()),
                                      bot_core_insulation=comma_str_to_point_float(
                                          self.aut_isolation_core2cond_bot_lineEdit.text()),
                                      left_core_insulation=comma_str_to_point_float(
                                          self.aut_isolation_core2cond_inner_lineEdit.text()),
                                      right_core_insulation=comma_str_to_point_float(
                                          self.aut_isolation_core2cond_outer_lineEdit.text()),
                                      inner_winding_insulation=comma_str_to_point_float(
                                          self.aut_isolation_p2p_lineEdit.text()),
                                      temperature=comma_str_to_point_float(self.aut_temp_lineEdit.text()),
                                      manual_litz_conductor_r=[],
                                      manual_litz_strand_r=[],
                                      manual_litz_strand_n=[],
                                      manual_litz_fill_factor=[])

        # Create csv file of data_matrix_fem which consist of all the fem simulation cases details
        self.ad.write_data_matrix_fem_to_csv()

        self.plot_volume_loss(self.ad.data_matrix_4, matplotlib_widget)

        n_cases_0 = len(self.ad.data_matrix_0)
        self.ncases0_label.setText(f"{n_cases_0}")

        n_cases_2 = len(self.ad.data_matrix_2)
        self.ncases2_label.setText(f"{n_cases_2}")

        n_cases_3 = len(self.ad.data_matrix_3)
        self.ncases3_label.setText(f"{n_cases_3}")

        n_cases_FEM = len(self.ad.data_matrix_fem)
        self.fem_cases_label.setText(f"{n_cases_FEM}")

    def automated_design_fem_sim(self, matplotlib_widget):
        """
        Run the fem simulations of the filtered cases from reluctance models tab. Plot the volume vs loss in FEM simulations tab.

        :param matplotlib_widget: To plot volume vs loss in FEM simulations tab
        """
        # ##########################################   {FEM_SIMULATION}   ##################################################

        # Run FEM simulation of "self.data_matrix_fem"
        self.ad.fem_simulation()

        # Save simulation settings in json file for later review
        self.ad.save_automated_design_settings()
        design_directory = self.aut_load_design_directoryname_lineEdit.text()
        real_inductance, total_loss, total_volume, total_cost, labels, automated_design_settings = load_fem_simulation_results(
            working_directory=design_directory)

        matplotlib_widget = MatplotlibWidget()
        matplotlib_widget.axis.clear()
        self.layout = QVBoxLayout(self.plotwidget_5)
        self.layout.addWidget(matplotlib_widget)
        try:
            matplotlib_widget.axis_cm.remove()
        except:
            pass

        plot_data = filter_after_fem(inductance=real_inductance, total_loss=total_loss, total_volume=total_volume,
                                     total_cost=total_cost,
                                     annotation_list=labels, goal_inductance=self.ad.goal_inductance,
                                     percent_tolerance=20)

        self.plot_2d(matplotlib_widget, x_value=plot_data[:, 1], y_value=plot_data[:, 2], z_value=plot_data[:, 3],
                     x_label='Volume / m\u00b3', y_label='Loss / W', z_label='Cost / \u20ac', title='Volume vs Loss',
                     annotations=plot_data[:, 4], plot_color='RdYlGn_r', inductance_value=plot_data[:, 0])

    def load_designs(self, matplotlib_widget):
        """
         Plot the volume vs loss from the already run files of FEM simulations from the directory path, in the Load(results) tab.

        :param matplotlib_widget: To plot volume vs loss in the Load(results) tab
        """
        matplotlib_widget = MatplotlibWidget()
        matplotlib_widget.axis.clear()
        self.layout = QVBoxLayout(self.plotwidget_9)
        self.layout.addWidget(matplotlib_widget)
        try:
            matplotlib_widget.axis_cm.remove()
        except:
            pass

        design_directory = self.aut_load_design_directoryname_lineEdit.text()
        real_inductance, total_loss, total_volume, total_cost, labels, automated_design_settings = load_fem_simulation_results(
            working_directory=design_directory)

        plot_data = filter_after_fem(inductance=real_inductance, total_loss=total_loss, total_volume=total_volume,
                                     total_cost=total_cost,
                                     annotation_list=labels, goal_inductance=0.00012,
                                     percent_tolerance=20)

        self.plot_2d(matplotlib_widget, x_value=plot_data[:, 1], y_value=plot_data[:, 2], z_value=plot_data[:, 3],
                     x_label='Volume / m\u00b3', y_label='Loss / W', z_label='Cost / \u20ac', title='Volume vs Loss',
                     annotations=plot_data[:, 4], plot_color='RdYlGn_r', inductance_value=plot_data[:, 0])

    def check_onelab_config(self, geo: fmt.MagneticComponent):
        """Check the onlab configuration to enter or read the onelab filepath."""
        # Ask for onelab path (if no config file exists)
        if not os.path.isfile(geo.file_data.config_path):
            onelab_path_dialog = OnelabPathDialog()
            valid = onelab_path_dialog.exec_()
            if valid is None:
                raise Exception("Something went wrong while entering onelab path")
            elif valid == 1:
                onelab_path = onelab_path_dialog.directory
                if os.path.isdir(onelab_path):
                    onelab_path_dict = {"onelab": onelab_path}
                    with open(geo.file_data.config_path, 'w', encoding='utf-8') as fd:
                        json.dump(onelab_path_dict, fd, indent=2, ensure_ascii=False)
                else:
                    raise Exception(f"{onelab_path} is not a valid folder path")
            elif valid == 0:
                sys.exit(-1)
            else:
                raise Exception(f"Unknown return type from OnelabPathDialog: {valid}")

    def automated_design_func_config(self):
        """Call automated_design_func, when simulate button is pressed in Reluctance models tab."""
        self.automated_design_func(self.matplotlib_widget_aut_tab2)

    def automated_design_fem_sim_config(self):
        """Call automated_design_fem_sim, when simulate button is pressed in FEM simulations tab."""
        self.automated_design_fem_sim(self.matplotlib_widget_aut_tab3)

    def load_designs_config(self):
        """Call load_designs, when Load design button is pressed in Load(results) tab."""
        self.load_designs(self.matplotlib_widget_aut_tab4)

    #  **************************** Automated design tab initializations ********************************************  #
    def oncgeoClearallClicked(self):
        """Clear all entries."""
        self.aut_core_geometry_basket_listWidget.clear()

    def oncgeoClearClicked(self):
        """Add the manually selected core choice to the basket."""
        List_item = self.aut_core_geometry_basket_listWidget.selectedItems()
        for item in List_item:
            self.aut_core_geometry_basket_listWidget.takeItem(self.aut_core_geometry_basket_listWidget.row(item))

    def oncgeoMultipleClicked(self):
        """Accept multiple choices."""
        itemsTextList = [str(self.aut_core_geometry_basket_listWidget.item(i).text()) for i in
                         range(self.aut_core_geometry_basket_listWidget.count())]
        checkitems = [item.text() for item in self.aut_core_geometry_listWidget.selectedItems()]
        reqlist = list(set(checkitems).difference(itemsTextList))
        for i in reqlist:
            self.aut_core_geometry_basket_listWidget.addItem(i)
        else:
            pass

    def oncgeomanualMultipleClicked(self):
        """Manual."""
        items = []

        if self.aut_min_core_width_lineEdit.text() and self.aut_max_core_width_lineEdit.text() and self.aut_step_core_width_lineEdit.text() \
                and self.aut_min_window_height_lineEdit.text() and self.aut_max_window_height_lineEdit.text() and self.aut_step_window_height_lineEdit.text() \
                and self.aut_min_window_width_lineEdit.text() and self.aut_max_window_width_lineEdit.text() and self.aut_step_window_width_lineEdit.text():
            items.append(
                f"{self.aut_min_core_width_lineEdit.text()}*"
                f"{self.aut_min_window_height_lineEdit.text()}*{self.aut_min_window_width_lineEdit.text()} (min value)")
        itemsTextList = [str(self.aut_core_geometry_manual_basket_listWidget.item(i).text()) for i in
                         range(self.aut_core_geometry_manual_basket_listWidget.count())]
        for i in items:
            if i not in itemsTextList:
                self.aut_core_geometry_manual_basket_listWidget.addItem(i)

    def cgeoselectall(self):
        """Select all the choices."""
        self.aut_core_geometry_listWidget.selectAll()

    def oncgeoClicked(self):
        """Add a choice by click to the basket if it is not already in the basket."""
        itemsTextList = [str(self.aut_core_geometry_basket_listWidget.item(i).text()) for i in
                         range(self.aut_core_geometry_basket_listWidget.count())]
        checkitem = self.aut_core_geometry_listWidget.currentItem().text()
        if checkitem not in itemsTextList:
            self.aut_core_geometry_basket_listWidget.addItem(self.aut_core_geometry_listWidget.currentItem().text())
        else:
            pass

    def onairgaptypeClearallClicked(self):
        """Clear all the entries."""
        self.aut_airgap_type_basket_listwidget.clear()

    def onairgaptypeClearClicked(self):
        """Clear the selected entry."""
        List_item = self.aut_airgap_type_basket_listwidget.selectedItems()
        for item in List_item:
            self.aut_airgap_type_basket_listwidget.takeItem(self.aut_airgap_type_basket_listwidget.row(item))

    def onairgaptypeMultipleClicked(self):
        """Accept multiple choices."""
        itemsTextList = [str(self.aut_airgap_type_basket_listwidget.item(i).text()) for i in
                         range(self.aut_airgap_type_basket_listwidget.count())]
        checkitems = [item.text() for item in self.aut_airgap_type_listWidget.selectedItems()]
        reqlist = list(set(checkitems).difference(itemsTextList))
        for i in reqlist:
            self.aut_airgap_type_basket_listwidget.addItem(i)
        else:
            pass

    def onairgaptypeClicked(self):
        """Add a choice by click to the basket if it is not already in the basket."""
        itemsTextList = [str(self.aut_airgap_type_basket_listwidget.item(i).text()) for i in
                         range(self.aut_airgap_type_basket_listwidget.count())]
        checkitem = self.aut_airgap_type_listWidget.currentItem().text()
        if checkitem not in itemsTextList:
            self.aut_airgap_type_basket_listwidget.addItem(self.aut_airgap_type_listWidget.currentItem().text())
        else:
            pass

    def airgaptypeselectall(self):
        """Select all the choices."""
        self.aut_airgap_type_listWidget.selectAll()

    def oncmatClearallClicked(self):
        """Clear all entries."""
        self.aut_core_material_basket_listWidget.clear()

    def oncmatClearClicked(self):
        """Clear the selected entry."""
        List_item = self.aut_core_material_basket_listWidget.selectedItems()
        for item in List_item:
            self.aut_core_material_basket_listWidget.takeItem(self.aut_core_material_basket_listWidget.row(item))

    def oncmatMultipleClicked(self):
        """Accept multiple choices."""
        itemsTextList = [str(self.aut_core_material_basket_listWidget.item(i).text()) for i in
                         range(self.aut_core_material_basket_listWidget.count())]
        checkitems = [item.text() for item in self.aut_core_material_data_listWidget.selectedItems()]
        reqlist = list(set(checkitems).difference(itemsTextList))
        for i in reqlist:
            self.aut_core_material_basket_listWidget.addItem(i)
        else:
            pass

    def cmatselectall(self):
        """Select all the choices."""
        self.aut_core_material_data_listWidget.selectAll()

    def oncmatClicked(self):
        """Add a choice by click to the basket if it is not already in the basket."""
        itemsTextList = [str(self.aut_core_material_basket_listWidget.item(i).text()) for i in
                         range(self.aut_core_material_basket_listWidget.count())]
        checkitem = self.aut_core_material_data_listWidget.currentItem().text()
        if checkitem not in itemsTextList:
            self.aut_core_material_basket_listWidget.addItem(
                self.aut_core_material_data_listWidget.currentItem().text())
        else:
            pass

    def onl1ClearallClicked(self):
        """Clear all entries."""
        self.aut_litz_basket_listWidget.clear()

    def onl1ClearClicked(self):
        """Clear the selected entry."""
        List_item = self.aut_litz_basket_listWidget.selectedItems()
        for item in List_item:
            self.aut_litz_basket_listWidget.takeItem(self.aut_litz_basket_listWidget.row(item))

    def onl1MultipleClicked(self):
        """Accept multiple choices."""
        itemsTextList = [str(self.aut_litz_basket_listWidget.item(i).text()) for i in
                         range(self.aut_litz_basket_listWidget.count())]
        checkitems = [item.text() for item in self.aut_litz_data_listWidget.selectedItems()]
        reqlist = list(set(checkitems).difference(itemsTextList))
        for i in reqlist:
            self.aut_litz_basket_listWidget.addItem(i)
        else:
            pass

    def onl1Clicked(self):
        """Add a choice by click to the basket if it is not already in the basket."""
        itemsTextList = [str(self.aut_litz_basket_listWidget.item(i).text()) for i in
                         range(self.aut_litz_basket_listWidget.count())]
        checkitem = self.aut_litz_data_listWidget.currentItem().text()
        if checkitem not in itemsTextList:
            self.aut_litz_basket_listWidget.addItem(self.aut_litz_data_listWidget.currentItem().text())
        else:
            pass

    def litz1selectall(self):
        """Select all the choices."""
        self.aut_litz_data_listWidget.selectAll()

    def onl2ClearallClicked(self):
        """Clear all entries."""
        self.aut_litz2_basket_listWidget.clear()

    def onl2ClearClicked(self):
        """Clear the selected entry."""
        List_item = self.aut_litz2_basket_listWidget.selectedItems()
        for item in List_item:
            self.aut_litz2_basket_listWidget.takeItem(self.aut_litz2_basket_listWidget.row(item))

    def onl2MultipleClicked(self):
        """Accept multiple choices."""
        itemsTextList = [str(self.aut_litz2_basket_listWidget.item(i).text()) for i in
                         range(self.aut_litz2_basket_listWidget.count())]
        checkitems = [item.text() for item in self.aut_litz2_data_listWidget.selectedItems()]
        reqlist = list(set(checkitems).difference(itemsTextList))
        for i in reqlist:
            self.aut_litz2_basket_listWidget.addItem(i)
        else:
            pass

    def onl2Clicked(self):
        """Add a choice by click to the basket if it is not already in the basket."""
        itemsTextList = [str(self.aut_litz2_basket_listWidget.item(i).text()) for i in
                         range(self.aut_litz2_basket_listWidget.count())]
        checkitem = self.aut_litz2_data_listWidget.currentItem().text()
        if checkitem not in itemsTextList:
            self.aut_litz2_basket_listWidget.addItem(self.aut_litz2_data_listWidget.currentItem().text())
        else:
            pass

    def litz2selectall(self):
        """Select all the choices."""
        self.aut_litz2_data_listWidget.selectAll()

    def aut_initialize_controls(self) -> None:
        """
        Initialize the comboboxes with pre-defined values.

        :return: None
        :rtype: None
        """
        aut_simulation_type_options = [self.translation_dict['inductor']]
        aut_winding_material_options = [key for key in fmt.wire_material_database()]
        aut_winding_type_options = [self.translation_dict['litz'], self.translation_dict['solid']]
        aut_implicit_litz_options = [self.translation_dict["implicit_litz_radius"],
                                     self.translation_dict["implicit_ff"],
                                     self.translation_dict['implicit_strands_number']]
        aut_winding_scheme_options = [self.translation_dict["square"], self.translation_dict["hexa"]]
        aut_tolerance_val_options = [self.translation_dict['+-10'], self.translation_dict['+-20']]
        aut_core_geometry_options = [core_geometry for core_geometry in fmt.core_database()]

        get_material_list = database.material_list_in_database()
        get_material_list.insert(0, None)
        dat_core_material_options = get_material_list

        for option in dat_core_material_options:
            self.dat_core_material1_comboBox.addItem(option)
        for option in dat_core_material_options:
            self.dat_core_material2_comboBox.addItem(option)
        for option in dat_core_material_options:
            self.dat_core_material3_comboBox.addItem(option)
        for option in dat_core_material_options:
            self.dat_core_material4_comboBox.addItem(option)
        for option in dat_core_material_options:
            self.dat_core_material5_comboBox.addItem(option)

        for option in dat_core_material_options:
            self.dat_core_material1_comboBox_2.addItem(option)
        for option in dat_core_material_options:
            self.dat_core_material2_comboBox_2.addItem(option)
        for option in dat_core_material_options:
            self.dat_core_material3_comboBox_2.addItem(option)
        for option in dat_core_material_options:
            self.dat_core_material4_comboBox_2.addItem(option)
        for option in dat_core_material_options:
            self.dat_core_material5_comboBox_2.addItem(option)

        for option in dat_core_material_options:
            self.dat_core_material_comboBox.addItem(option)

        for option in aut_core_geometry_options:
            self.aut_core_geometry_listWidget.addItem(option)
        for option in aut_simulation_type_options:
            self.aut_simulation_type_comboBox.addItem(option)
        for option in aut_winding_material_options:
            self.aut_winding1_material_comboBox.addItem(option)
            self.aut_winding2_material_comboBox.addItem(option)
        for option in aut_winding_type_options:
            self.aut_winding1_type_comboBox.addItem(option)
            self.aut_winding2_type_comboBox.addItem(option)
        for option in aut_implicit_litz_options:
            self.aut_winding1_implicit_litz_comboBox.addItem(option)
            self.aut_winding2_implicit_litz_comboBox.addItem(option)
        for option in aut_winding_scheme_options:
            self.aut_winding1_scheme_comboBox.addItem(option)
            self.aut_winding2_scheme_comboBox.addItem(option)
        for option in aut_tolerance_val_options:
            self.aut_rel_tolerance_val_comboBox.addItem(option)
        for option in aut_tolerance_val_options:
            self.aut_load_tolerance_val_comboBox.addItem(option)

        self.aut_min_core_width_lineEdit.setPlaceholderText("Minimum value")
        self.aut_max_core_width_lineEdit.setPlaceholderText("Maximum value")
        self.aut_step_core_width_lineEdit.setPlaceholderText("Step value")
        self.aut_min_window_height_lineEdit.setPlaceholderText("Minimum value")
        self.aut_max_window_height_lineEdit.setPlaceholderText("Maximum value")
        self.aut_step_window_height_lineEdit.setPlaceholderText("Step value")
        self.aut_min_window_width_lineEdit.setPlaceholderText("Minimum value")
        self.aut_max_window_width_lineEdit.setPlaceholderText("Maximum value")
        self.aut_step_window_width_lineEdit.setPlaceholderText("Step value")
        self.aut_min_winding1_turns_lineEdit.setPlaceholderText("Minimum value")
        self.aut_max_winding1_turns_lineEdit.setPlaceholderText("Maximum value")
        self.aut_min_winding2_turns_lineEdit.setPlaceholderText("Minimum value")
        self.aut_max_winding2_turns_lineEdit.setPlaceholderText("Maximum value")
        self.aut_step_winding2_turns_lineEdit.setPlaceholderText("Step value")
        self.aut_min_air_gap_count_lineEdit.setPlaceholderText("Minimum value")
        self.aut_max_air_gap_count_lineEdit.setPlaceholderText("Maximum value")
        self.aut_goal_inductance_val_lineEdit.setPlaceholderText(" Value in Henry")
        self.aut_air_gap_length_min_lineEdit.setPlaceholderText("Minimum value")
        self.aut_air_gap_length_max_lineEdit.setPlaceholderText("Maximum value")
        self.aut_air_gap_length_step_lineEdit.setPlaceholderText("Step value")
        self.aut_air_gap_position_min_lineEdit.setPlaceholderText("Minimum value")
        self.aut_air_gap_position_max_lineEdit.setPlaceholderText("Maximum value")
        self.aut_air_gap_position_step_lineEdit.setPlaceholderText("Step value")

    def aut_winding1_change_litz_implicit(self, implicit_type_from_combo_box: str) -> None:
        """
        Enable / Disable input parameter fields for different "implicit xyz" types in case of litz wire.

        :param implicit_type_from_combo_box: input type to implicit
        :type implicit_type_from_combo_box: str
        :return: None
        :rtype: None
        """
        if implicit_type_from_combo_box == self.translation_dict['implicit_litz_radius']:
            self.aut_winding1_strands_lineEdit.setEnabled(True)
            self.aut_winding1_fill_factor_lineEdit.setEnabled(True)
            self.aut_winding1_strand_radius_lineEdit.setEnabled(True)
            self.aut_winding1_radius_lineEdit.setEnabled(False)
        if implicit_type_from_combo_box == self.translation_dict['implicit_strands_number']:
            self.aut_winding1_strands_lineEdit.setEnabled(False)
            self.aut_winding1_fill_factor_lineEdit.setEnabled(True)
            self.aut_winding1_strand_radius_lineEdit.setEnabled(True)
            self.aut_winding1_radius_lineEdit.setEnabled(True)
        if implicit_type_from_combo_box == self.translation_dict['implicit_ff']:
            self.aut_winding1_strands_lineEdit.setEnabled(True)
            self.aut_winding1_fill_factor_lineEdit.setEnabled(False)
            self.aut_winding1_strand_radius_lineEdit.setEnabled(True)
            self.aut_winding1_radius_lineEdit.setEnabled(True)

    def aut_winding1_change_wire_type(self, wire_type_from_combot_box: str) -> None:
        """
        Enable / Disable input parameter for litz/solid wire.

        :param wire_type_from_combot_box: wire type
        :type wire_type_from_combot_box: str
        :return: None
        :rtype: None
        """
        self.aut_winding1_change_litz_implicit(self.aut_winding1_implicit_litz_comboBox.currentText())
        if wire_type_from_combot_box == self.translation_dict['litz']:
            self.aut_winding1_strands_lineEdit.setEnabled(True)
            self.aut_winding1_implicit_litz_comboBox.setEnabled(True)
            self.aut_winding1_fill_factor_lineEdit.setEnabled(True)
            self.aut_winding1_strand_radius_lineEdit.setEnabled(True)
            self.aut_winding1_radius_lineEdit.setEnabled(True)
            self.aut_litz_data_listWidget.setEnabled(True)
            self.aut_winding1_change_litz_implicit(self.aut_winding1_implicit_litz_comboBox.currentText())

        elif wire_type_from_combot_box == self.translation_dict['solid']:
            self.aut_winding1_strands_lineEdit.setEnabled(False)
            self.aut_winding1_implicit_litz_comboBox.setEnabled(False)
            self.aut_winding1_fill_factor_lineEdit.setEnabled(False)
            self.aut_winding1_strand_radius_lineEdit.setEnabled(False)
            self.aut_winding1_radius_lineEdit.setEnabled(True)
            self.aut_litz_data_listWidget.setEnabled(False)

    def aut_change_simulation_type(self, simulation_type_from_combo_box: str) -> None:
        """
        Perform action when signal of aut_simulation_type_comboBox text has changed.

        Action will be enabling / disabling user inputs for not-used windings.

        :param simulation_type_from_combo_box:
        :type simulation_type_from_combo_box: str
        :return: None
        :rtype: None
        """
        if simulation_type_from_combo_box == self.translation_dict['inductor']:
            self.aut_winding2_enable(False)

        elif simulation_type_from_combo_box == self.translation_dict['transformer']:
            # set winding definitions of winding 2 to editable
            self.aut_winding2_enable(True)

        elif simulation_type_from_combo_box == self.translation_dict['integrated transformer']:
            # set winding definitions of winding 2 to editable
            self.aut_winding2_enable(True)

    def aut_winding2_enable(self, status: bool) -> None:
        """
        Enable/disable all fields being in contact with winding 2.

        :param status: True / False
        :type status: bool
        :return: None
        :rtype: None
        """
        # set winding definitions of winding 2 (enable and visible)
        self.aut_winding2_material_comboBox.setEnabled(status)
        self.aut_winding2_type_comboBox.setEnabled(status)
        self.aut_min_winding2_turns_lineEdit.setEnabled(status)
        self.aut_max_winding2_turns_lineEdit.setEnabled(status)
        self.aut_step_winding2_turns_lineEdit.setEnabled(status)
        self.aut_winding2_strands_lineEdit.setEnabled(status)
        self.aut_winding2_radius_lineEdit.setEnabled(False)
        self.aut_winding2_fill_factor_lineEdit.setEnabled(status)
        self.aut_winding2_strand_radius_lineEdit.setEnabled(status)
        self.aut_winding2_groupBox.setVisible(status)

        # Set turns of winding 2 (enable and visible)
        self.aut_min_winding2_turns_lineEdit.setVisible(status)
        self.aut_max_winding2_turns_lineEdit.setVisible(status)
        self.aut_step_winding2_turns_lineEdit.setVisible(status)
        self.aut_winding2_scheme_comboBox.setVisible(status)
        self.aut_winding2_turns_label.setVisible(status)
        self.aut_winding2_scheme_label.setVisible(status)

        # set isolation of winding 2 (enable and visible)
        self.aut_isolation_s2s_lineEdit.setEnabled(status)
        self.aut_isolation_p2s_lineEdit.setEnabled(status)
        self.aut_isolation_s2s_lineEdit.setVisible(status)
        self.aut_isolation_p2s_lineEdit.setVisible(status)
        self.aut_isolation_s2s_label.setVisible(status)
        self.aut_isolation_p2s_label.setVisible(status)

        self.aut_litz2_basket_listWidget.setEnabled(status)
        self.aut_litz2_basket_listWidget.setVisible(status)
        self.aut_litz2_basket_clear_all_pushbutton.setVisible(status)
        self.aut_litz2_basket_clear_pushbutton.setVisible(status)
        self.aut_litzbasket2_label.setVisible(status)

    #  **************************** Database tab ********************************************************************  #

    def datupdateraph1(self, matplotlib_widget1, matplotlib_widget2, matplotlib_widget3, matplotlib_widget4):
        """
        Datasheet-datasheet plot.

        :param matplotlib_widget1: for the first plot of relative power loss vs B
        :param matplotlib_widget2: for the second plot of relative power loss vs temperature
        :param matplotlib_widget3: for the third plot of relative power loss vs frequency
        :param matplotlib_widget4: for the fourth plot of B vs H
        """
        matplotlib_widget1.axis.clear()
        self.layout = QVBoxLayout(self.plotwidget)
        self.layout.addWidget(matplotlib_widget1)
        try:
            matplotlib_widget1.axis_cm.remove()
        except:
            pass

        mat1_name = self.dat_core_material1_comboBox.currentText()
        mat2_name = self.dat_core_material2_comboBox.currentText()
        mat3_name = self.dat_core_material3_comboBox.currentText()
        mat4_name = self.dat_core_material4_comboBox.currentText()
        mat5_name = self.dat_core_material5_comboBox.currentText()

        mat1_temp = comma_str_to_point_float(self.aut_temp_m1_comboBox.currentText())
        mat2_temp = comma_str_to_point_float(self.aut_temp_m2_comboBox.currentText())
        mat3_temp = comma_str_to_point_float(self.aut_temp_m3_comboBox.currentText())
        mat4_temp = comma_str_to_point_float(self.aut_temp_m4_comboBox.currentText())
        mat5_temp = comma_str_to_point_float(self.aut_temp_m5_comboBox.currentText())

        mat1_flux = comma_str_to_point_float(self.aut_flux_m1_comboBox.currentText())
        mat2_flux = comma_str_to_point_float(self.aut_flux_m2_comboBox.currentText())
        mat3_flux = comma_str_to_point_float(self.aut_flux_m3_comboBox.currentText())
        mat4_flux = comma_str_to_point_float(self.aut_flux_m4_comboBox.currentText())
        mat5_flux = comma_str_to_point_float(self.aut_flux_m5_comboBox.currentText())

        # print(f"mat1_name: {mat1_name},{mat2_name},{mat3_name},{mat4_name},{mat5_name}")
        # print(f"mat1_temp: {mat1_temp},{mat2_temp},{mat3_temp},{mat4_temp},{mat5_temp}")
        # print(f"mat1_flux: {mat1_flux},{mat2_flux},{mat3_flux},{mat4_flux},{mat5_flux}")

        materials_used_list = []
        material_list = [mat1_name, mat2_name, mat3_name, mat4_name, mat5_name]
        for items in material_list:
            if items:
                materials_used_list.append(items)
        # print(materials_used_list)

        database.compare_core_loss_flux_density_data(matplotlib_widget1, material_list=materials_used_list,
                                                     temperature_list=[mat1_temp, mat2_temp, mat3_temp, mat4_temp,
                                                                       mat5_temp])
        # self.matplotlib_widget1.axis.legend(fontsize=13)
        matplotlib_widget1.axis.grid()
        matplotlib_widget1.figure.canvas.draw_idle()
        matplotlib_widget1.figure.tight_layout()

        ################################################################################################################

        matplotlib_widget2.axis.clear()
        self.layout = QVBoxLayout(self.plotwidget_2)
        self.layout.addWidget(matplotlib_widget2)
        try:
            matplotlib_widget2.axis_cm.remove()
        except:
            pass

        flux_list = [mat1_flux, mat2_flux, mat3_flux, mat4_flux, mat5_flux]
        # print(f"flux_list: {flux_list}")
        database.compare_core_loss_temperature(matplotlib_widget2, material_list=materials_used_list,
                                               flux_density_list=[mat1_flux, mat2_flux, mat3_flux, mat4_flux,
                                                                  mat5_flux])
        # self.matplotlib_widget2.axis.legend(fontsize=13)
        matplotlib_widget2.axis.grid()
        matplotlib_widget2.figure.canvas.draw_idle()
        matplotlib_widget2.figure.tight_layout()

        ################################################################################################################

        matplotlib_widget3.axis.clear()
        self.layout = QVBoxLayout(self.plotwidget_3)
        self.layout.addWidget(matplotlib_widget3)
        try:
            matplotlib_widget3.axis_cm.remove()
        except:
            pass

        database.compare_core_loss_frequency(matplotlib_widget3, material_list=materials_used_list,
                                             temperature_list=[mat1_temp, mat2_temp, mat3_temp, mat4_temp, mat5_temp],
                                             flux_density_list=[mat1_flux, mat2_flux, mat3_flux, mat4_flux, mat5_flux])
        # self.matplotlib_widget3.axis.legend(fontsize=13)
        matplotlib_widget3.axis.grid()
        matplotlib_widget3.figure.canvas.draw_idle()
        matplotlib_widget3.figure.tight_layout()

        ################################################################################################################

        matplotlib_widget4.axis.clear()
        self.layout = QVBoxLayout(self.plotwidget_4)
        self.layout.addWidget(matplotlib_widget4)
        try:
            matplotlib_widget4.axis_cm.remove()
        except:
            pass

        database.compare_b_h_curve(matplotlib_widget4, material_list=materials_used_list,
                                   temperature_list=[mat1_temp, mat2_temp, mat3_temp, mat4_temp, mat5_temp])
        # self.matplotlib_widget4.axis.legend(fontsize=13)
        matplotlib_widget4.axis.grid()
        matplotlib_widget4.figure.canvas.draw_idle()
        matplotlib_widget4.figure.tight_layout()

    def datupdateraph2(self, matplotlib_widget1, matplotlib_widget2):
        """
        Measurement-Measurement plot.

        :param matplotlib_widget1: Fot the first plot of uR/u0 vs B
        :param matplotlib_widget2: for the second plot of uR/u0 vs B
        """
        matplotlib_widget1.axis.clear()
        self.layout = QVBoxLayout(self.plotwidget_13)
        self.layout.addWidget(matplotlib_widget1)
        try:
            matplotlib_widget1.axis_cm.remove()
        except:
            pass

        mat1_name = self.dat_core_material1_comboBox_2.currentText()
        mat2_name = self.dat_core_material2_comboBox_2.currentText()
        mat3_name = self.dat_core_material3_comboBox_2.currentText()
        mat4_name = self.dat_core_material4_comboBox_2.currentText()
        mat5_name = self.dat_core_material5_comboBox_2.currentText()

        mat1_test_name = self.test_name_1_comboBox.currentText()
        mat2_test_name = self.test_name_2_comboBox.currentText()
        mat3_test_name = self.test_name_3_comboBox.currentText()
        mat4_test_name = self.test_name_4_comboBox.currentText()
        mat5_test_name = self.test_name_5_comboBox.currentText()

        mat1_temp = comma_str_to_point_float(self.aut_temp_m1_comboBox_2.currentText())
        mat2_temp = comma_str_to_point_float(self.aut_temp_m2_comboBox_2.currentText())
        mat3_temp = comma_str_to_point_float(self.aut_temp_m3_comboBox_2.currentText())
        mat4_temp = comma_str_to_point_float(self.aut_temp_m4_comboBox_2.currentText())
        mat5_temp = comma_str_to_point_float(self.aut_temp_m5_comboBox_2.currentText())

        mat1_freq = comma_str_to_point_float(self.aut_freq_m1_comboBox.currentText())
        mat2_freq = comma_str_to_point_float(self.aut_freq_m2_comboBox.currentText())
        mat3_freq = comma_str_to_point_float(self.aut_freq_m3_comboBox.currentText())
        mat4_freq = comma_str_to_point_float(self.aut_freq_m4_comboBox.currentText())
        mat5_freq = comma_str_to_point_float(self.aut_freq_m5_comboBox.currentText())

        # print(f"mat1_name: {mat1_name},{mat2_name},{mat3_name},{mat4_name},{mat5_name}")
        # print(f"mat1_temp: {mat1_temp},{mat2_temp},{mat3_temp},{mat4_temp},{mat5_temp}")
        # print(f"mat1_freq: {mat1_freq},{mat2_freq},{mat3_freq},{mat4_freq},{mat5_freq}")

        materials_used_list = []
        material_list = [mat1_name, mat2_name, mat3_name, mat4_name, mat5_name]
        for items in material_list:
            if items:
                materials_used_list.append(items)
        # print(materials_used_list)

        database.compare_permeability_measurement_data(matplotlib_widget1, material_list=materials_used_list,
                                                       measurement_name=[mat1_test_name, mat2_test_name, mat3_test_name,
                                                                         mat4_test_name, mat5_test_name],
                                                       frequency_list=[mat1_freq, mat2_freq, mat3_freq, mat4_freq,
                                                                       mat5_freq],
                                                       temperature_list=[mat1_temp, mat2_temp, mat3_temp, mat4_temp,
                                                                         mat5_temp],
                                                       plot_real_part=True)

        matplotlib_widget1.axis.grid()
        matplotlib_widget1.figure.canvas.draw_idle()
        matplotlib_widget1.figure.tight_layout()

        ################################################################################################################

        matplotlib_widget2.axis.clear()
        self.layout = QVBoxLayout(self.plotwidget_14)
        self.layout.addWidget(matplotlib_widget2)
        try:
            matplotlib_widget2.axis_cm.remove()
        except:
            pass
        database.compare_permeability_measurement_data(matplotlib_widget2, material_list=materials_used_list,
                                                       measurement_name=[mat1_test_name, mat2_test_name, mat3_test_name,
                                                                         mat4_test_name, mat5_test_name],
                                                       frequency_list=[mat1_freq, mat2_freq, mat3_freq, mat4_freq,
                                                                       mat5_freq],
                                                       temperature_list=[mat1_temp, mat2_temp, mat3_temp, mat4_temp,
                                                                         mat5_temp],
                                                       plot_real_part=False)
        matplotlib_widget2.axis.grid()
        matplotlib_widget2.figure.canvas.draw_idle()
        matplotlib_widget2.figure.tight_layout()

    def datupdateraph3(self, matplotlib_widget):
        """Datasheet-Measurement plot.

        :param matplotlib_widget: To plot relative power loss vs B
        """
        matplotlib_widget.axis.clear()
        self.layout = QVBoxLayout(self.plotwidget_15)
        self.layout.addWidget(matplotlib_widget)
        try:
            matplotlib_widget.axis_cm.remove()
        except:
            pass

        mat_dat_temp = comma_str_to_point_float(self.aut_temp_dat_comboBox.currentText())
        mat_meas_temp = comma_str_to_point_float(self.aut_temp_meas_comboBox.currentText())
        mat_name = self.dat_core_material_comboBox.currentText()

        database.compare_core_loss_flux_datasheet_measurement(matplotlib_widget, material=mat_name,
                                                              temperature_list=[mat_dat_temp, mat_meas_temp], measurement_name=fmt.MeasurementSetup.LEA_LK)

        matplotlib_widget.axis.grid()
        matplotlib_widget.figure.canvas.draw_idle()
        matplotlib_widget.figure.tight_layout()

    def datupdateraph1_config(self):
        """Call datupdateraph1, when Update preview button is pressed in Datasheet-Datasheet tab."""
        self.datupdateraph1(self.matplotlib_widget_datdd1, self.matplotlib_widget_datdd2,
                            self.matplotlib_widget_datdd3, self.matplotlib_widget_datdd4)

    def datupdateraph2_config(self):
        """Call datupdateraph2, when Update preview button is pressed in Measurement-Measurement tab."""
        self.datupdateraph2(self.matplotlib_widget_datmm1, self.matplotlib_widget_datmm2)

    def datupdateraph3_config(self):
        """Call datupdateraph3, when Update preview button is pressed in Datasheet-Measurement tab."""
        self.datupdateraph3(self.matplotlib_widget_datdm)

    def tempfluxinput1(self):
        """Get the flux and temperature of a particular material selected."""
        mat_text1 = self.dat_core_material1_comboBox.currentText()

        get_temp1_list = []
        get_flux1_list = []
        if mat_text1:
            get_temp1_list = database.drop_down_list(material_name=mat_text1, datatype=None, measurement_name=None,
                                                     comparison_type="dvd", temperature=True)
            get_flux1_list = database.drop_down_list(material_name=mat_text1, datatype=None, measurement_name=None,
                                                     comparison_type="dvd", flux_density=True)

        aut_temp_options1 = get_temp1_list
        aut_flux_options1 = get_flux1_list

        temp_str = [f'{item:.2f}' for item in aut_temp_options1]
        flux_str = [f'{item:.4f}' for item in aut_flux_options1]
        self.aut_temp_m1_comboBox.clear()
        self.aut_flux_m1_comboBox.clear()
        for option in temp_str:
            self.aut_temp_m1_comboBox.addItem(option)

        for option in flux_str:
            self.aut_flux_m1_comboBox.addItem(option)

    def tempfluxinput2(self):
        """Get the flux and temperature of a particular material selected."""
        mat_text2 = self.dat_core_material2_comboBox.currentText()
        get_temp2_list = []
        get_flux2_list = []
        if mat_text2:
            get_temp2_list = database.drop_down_list(material_name=mat_text2, datatype=None, measurement_name=None,
                                                     comparison_type="dvd", temperature=True)
            get_flux2_list = database.drop_down_list(material_name=mat_text2, datatype=None, measurement_name=None,
                                                     comparison_type="dvd", flux_density=True)
        aut_temp_options2 = get_temp2_list
        aut_flux_options2 = get_flux2_list

        temp_str = [f'{item:.2f}' for item in aut_temp_options2]
        flux_str = [f'{item:.4f}' for item in aut_flux_options2]
        self.aut_temp_m2_comboBox.clear()
        self.aut_flux_m2_comboBox.clear()
        for option in temp_str:
            self.aut_temp_m2_comboBox.addItem(option)

        for option in flux_str:
            self.aut_flux_m2_comboBox.addItem(option)

    def tempfluxinput3(self):
        """Get the flux and temperature of a particular material selected."""
        mat_text3 = self.dat_core_material3_comboBox.currentText()
        get_temp3_list = []
        get_flux3_list = []
        if mat_text3:
            get_temp3_list = database.drop_down_list(material_name=mat_text3, datatype=None, measurement_name=None,
                                                     comparison_type="dvd", temperature=True)
            get_flux3_list = database.drop_down_list(material_name=mat_text3, datatype=None, measurement_name=None,
                                                     comparison_type="dvd", flux_density=True)
        aut_temp_options3 = get_temp3_list
        aut_flux_options3 = get_flux3_list

        temp_str = [f'{item:.2f}' for item in aut_temp_options3]
        flux_str = [f'{item:.4f}' for item in aut_flux_options3]
        self.aut_temp_m3_comboBox.clear()
        self.aut_flux_m3_comboBox.clear()
        for option in temp_str:
            self.aut_temp_m3_comboBox.addItem(option)

        for option in flux_str:
            self.aut_flux_m3_comboBox.addItem(option)

    def tempfluxinput4(self):
        """Get the flux and temperature of a particular material selected."""
        mat_text4 = self.dat_core_material4_comboBox.currentText()
        get_temp4_list = []
        get_flux4_list = []
        if mat_text4:
            get_temp4_list = database.drop_down_list(material_name=mat_text4, datatype=None, measurement_name=None,
                                                     comparison_type="dvd", temperature=True)
            get_flux4_list = database.drop_down_list(material_name=mat_text4, datatype=None, measurement_name=None,
                                                     comparison_type="dvd", flux_density=True)
        aut_temp_options4 = get_temp4_list
        aut_flux_options4 = get_flux4_list

        temp_str = [f'{item:.2f}' for item in aut_temp_options4]
        flux_str = [f'{item:.4f}' for item in aut_flux_options4]
        self.aut_temp_m4_comboBox.clear()
        self.aut_flux_m4_comboBox.clear()
        for option in temp_str:
            self.aut_temp_m4_comboBox.addItem(option)

        for option in flux_str:
            self.aut_flux_m4_comboBox.addItem(option)

    def tempfluxinput5(self):
        """Get the flux and temperature of a particular material selected."""
        mat_text5 = self.dat_core_material5_comboBox.currentText()
        get_temp5_list = []
        get_flux5_list = []
        if mat_text5:
            get_temp5_list = database.drop_down_list(material_name=mat_text5, datatype=None, measurement_name=None,
                                                     comparison_type="dvd", temperature=True)
            get_flux5_list = database.drop_down_list(material_name=mat_text5, datatype=None, measurement_name=None,
                                                     comparison_type="dvd", flux_density=True)
        aut_temp_options5 = get_temp5_list
        aut_flux_options5 = get_flux5_list

        temp_str = [f'{item:.2f}' for item in aut_temp_options5]
        flux_str = [f'{item:.4f}' for item in aut_flux_options5]
        self.aut_temp_m5_comboBox.clear()
        self.aut_flux_m5_comboBox.clear()
        for option in temp_str:
            self.aut_temp_m5_comboBox.addItem(option)

        for option in flux_str:
            self.aut_flux_m5_comboBox.addItem(option)

    def test_setup_name1(self):
        """Get test setup names from database for particular material."""
        mat_text1 = self.dat_core_material1_comboBox_2.currentText()
        # Temporarily disconnect the signal to prevent triggering changes
        self.test_name_1_comboBox.currentTextChanged.disconnect(self.tempfreqinput1)
        # Clear existing data in the comboBox
        self.test_name_1_comboBox.clear()

        names_list = []

        if mat_text1:
            names_list = database.find_measurement_names(material_name=mat_text1, datatype="complex_permeability")

        for option in names_list:
            self.test_name_1_comboBox.addItem(option)
        # Reconnect the signal after the updates
        self.test_name_1_comboBox.currentTextChanged.connect(self.tempfreqinput1)

    def test_setup_name2(self):
        """Get test setup names from database for particular material."""
        mat_text2 = self.dat_core_material2_comboBox_2.currentText()
        # Temporarily disconnect the signal to prevent triggering changes
        self.test_name_2_comboBox.currentTextChanged.disconnect(self.tempfreqinput2)
        # Clear existing data in the comboBox
        self.test_name_2_comboBox.clear()

        names_list = []

        if mat_text2:
            names_list = database.find_measurement_names(material_name=mat_text2, datatype="complex_permeability")

        for option in names_list:
            self.test_name_2_comboBox.addItem(option)
        # Reconnect the signal after the updates
        self.test_name_2_comboBox.currentTextChanged.connect(self.tempfreqinput2)

    def test_setup_name3(self):
        """Get test setup names from database for particular material."""
        mat_text3 = self.dat_core_material3_comboBox_2.currentText()
        # Temporarily disconnect the signal to prevent triggering changes
        self.test_name_3_comboBox.currentTextChanged.disconnect(self.tempfreqinput3)
        # Clear existing data in the comboBox
        self.test_name_3_comboBox.clear()

        names_list = []

        if mat_text3:
            names_list = database.find_measurement_names(material_name=mat_text3, datatype="complex_permeability")

        for option in names_list:
            self.test_name_3_comboBox.addItem(option)
        # Reconnect the signal after the updates
        self.test_name_3_comboBox.currentTextChanged.connect(self.tempfreqinput3)

    def test_setup_name4(self):
        """Get test setup names from database for particular material."""
        mat_text4 = self.dat_core_material4_comboBox_2.currentText()
        # Temporarily disconnect the signal to prevent triggering changes
        self.test_name_4_comboBox.currentTextChanged.disconnect(self.tempfreqinput4)
        # Clear existing data in the comboBox
        self.test_name_4_comboBox.clear()

        names_list = []

        if mat_text4:
            names_list = database.find_measurement_names(material_name=mat_text4, datatype="complex_permeability")

        for option in names_list:
            self.test_name_4_comboBox.addItem(option)
        # Reconnect the signal after the updates
        self.test_name_4_comboBox.currentTextChanged.connect(self.tempfreqinput4)

    def test_setup_name5(self):
        """Get test setup names from database for particular material."""
        mat_text5 = self.dat_core_material5_comboBox_2.currentText()
        # Temporarily disconnect the signal to prevent triggering changes
        self.test_name_5_comboBox.currentTextChanged.disconnect(self.tempfreqinput5)
        # Clear existing data in the comboBox
        self.test_name_5_comboBox.clear()

        names_list = []

        if mat_text5:
            names_list = database.find_measurement_names(material_name=mat_text5, datatype="complex_permeability")

        for option in names_list:
            self.test_name_5_comboBox.addItem(option)
        # Reconnect the signal after the updates
        self.test_name_5_comboBox.currentTextChanged.connect(self.tempfreqinput5)

    def tempfreqinput1(self):
        """Get the frequency and temperature of a particular material selected."""
        mat_text1 = self.dat_core_material1_comboBox_2.currentText()
        test_setup_text1 = self.test_name_1_comboBox.currentText()

        get_temp1_list = []
        get_freq1_list = []
        if mat_text1:
            get_temp1_list = database.drop_down_list(material_name=mat_text1, comparison_type="mvm",
                                                     datatype="complex_permeability", measurement_name=test_setup_text1,
                                                     temperature=True)
            get_freq1_list = database.drop_down_list(material_name=mat_text1, comparison_type="mvm",
                                                     datatype="complex_permeability", measurement_name=test_setup_text1,
                                                     frequency=True)
        aut_temp_options1 = get_temp1_list
        aut_freq_options1 = get_freq1_list

        temp_str = [f'{item:.2f}' for item in aut_temp_options1]
        freq_str = [f'{item:.2f}' for item in aut_freq_options1]
        self.aut_temp_m1_comboBox_2.clear()
        self.aut_freq_m1_comboBox.clear()
        for option in temp_str:
            self.aut_temp_m1_comboBox_2.addItem(option)

        for option in freq_str:
            self.aut_freq_m1_comboBox.addItem(option)

    def tempfreqinput2(self):
        """Get the frequency and temperature of a particular material selected."""
        mat_text2 = self.dat_core_material2_comboBox_2.currentText()
        test_setup_text2 = self.test_name_2_comboBox.currentText()
        get_temp2_list = []
        get_freq2_list = []
        if mat_text2:
            get_temp2_list = database.drop_down_list(material_name=mat_text2, datatype="complex_permeability",
                                                     measurement_name=test_setup_text2, comparison_type="mvm",
                                                     temperature=True)
            get_freq2_list = database.drop_down_list(material_name=mat_text2, datatype="complex_permeability",
                                                     measurement_name=test_setup_text2, comparison_type="mvm",
                                                     frequency=True)
        aut_temp_options2 = get_temp2_list
        aut_freq_options2 = get_freq2_list

        temp_str = [f'{item:.2f}' for item in aut_temp_options2]
        freq_str = [f'{item:.2f}' for item in aut_freq_options2]
        self.aut_temp_m2_comboBox_2.clear()
        self.aut_freq_m2_comboBox.clear()
        for option in temp_str:
            self.aut_temp_m2_comboBox_2.addItem(option)

        for option in freq_str:
            self.aut_freq_m2_comboBox.addItem(option)

    def tempfreqinput3(self):
        """Get the frequency and temperature of a particular material selected."""
        mat_text3 = self.dat_core_material3_comboBox_2.currentText()
        test_setup_text3 = self.test_name_3_comboBox.currentText()
        get_temp3_list = []
        get_freq3_list = []
        if mat_text3:
            get_temp3_list = database.drop_down_list(material_name=mat_text3, datatype="complex_permeability",
                                                     measurement_name=test_setup_text3, comparison_type="mvm",
                                                     temperature=True)
            get_freq3_list = database.drop_down_list(material_name=mat_text3, datatype="complex_permeability",
                                                     measurement_name=test_setup_text3, comparison_type="mvm",
                                                     frequency=True)
        aut_temp_options3 = get_temp3_list
        aut_freq_options3 = get_freq3_list

        temp_str = [f'{item:.2f}' for item in aut_temp_options3]
        freq_str = [f'{item:.2f}' for item in aut_freq_options3]
        self.aut_temp_m3_comboBox_2.clear()
        self.aut_freq_m3_comboBox.clear()
        for option in temp_str:
            self.aut_temp_m3_comboBox_2.addItem(option)

        for option in freq_str:
            self.aut_freq_m3_comboBox.addItem(option)

    def tempfreqinput4(self):
        """Get the frequency and temperature of a particular material selected."""
        mat_text4 = self.dat_core_material4_comboBox_2.currentText()
        test_setup_text4 = self.test_name_4_comboBox.currentText()
        get_temp4_list = []
        get_freq4_list = []
        if mat_text4:
            get_temp4_list = database.drop_down_list(material_name=mat_text4, datatype="complex_permeability",
                                                     measurement_name=test_setup_text4, comparison_type="mvm",
                                                     temperature=True)
            get_freq4_list = database.drop_down_list(material_name=mat_text4, datatype="complex_permeability",
                                                     measurement_name=test_setup_text4, comparison_type="mvm",
                                                     frequency=True)
        aut_temp_options4 = get_temp4_list
        aut_freq_options4 = get_freq4_list

        temp_str = [f'{item:.2f}' for item in aut_temp_options4]
        freq_str = [f'{item:.2f}' for item in aut_freq_options4]
        self.aut_temp_m4_comboBox_2.clear()
        self.aut_freq_m4_comboBox.clear()
        for option in temp_str:
            self.aut_temp_m4_comboBox_2.addItem(option)

        for option in freq_str:
            self.aut_freq_m4_comboBox.addItem(option)

    def tempfreqinput5(self):
        """Get the frequency and temperature of a particular material selected."""
        mat_text5 = self.dat_core_material5_comboBox_2.currentText()
        test_setup_text5 = self.test_name_5_comboBox.currentText()
        get_temp5_list = []
        get_freq5_list = []
        if mat_text5:
            get_temp5_list = database.drop_down_list(material_name=mat_text5, datatype="complex_permeability",
                                                     measurement_name=test_setup_text5, comparison_type="mvm",
                                                     temperature=True)
            get_freq5_list = database.drop_down_list(material_name=mat_text5, datatype="complex_permeability",
                                                     measurement_name=test_setup_text5, comparison_type="mvm",
                                                     frequency=True)
        aut_temp_options5 = get_temp5_list
        aut_freq_options5 = get_freq5_list

        temp_str = [f'{item:.2f}' for item in aut_temp_options5]
        freq_str = [f'{item:.2f}' for item in aut_freq_options5]
        self.aut_temp_m5_comboBox_2.clear()
        self.aut_freq_m5_comboBox.clear()
        for option in temp_str:
            self.aut_temp_m5_comboBox_2.addItem(option)

        for option in freq_str:
            self.aut_freq_m5_comboBox.addItem(option)

        ########################################################################

    def temp_dat_input(self):
        """Get the database temperature of a particular material selected."""
        mat_text1 = self.dat_core_material_comboBox.currentText()

        get_temp1_list = []
        if mat_text1:
            get_temp1_list = database.drop_down_list(material_name=mat_text1, comparison_type="dvd", temperature=True,
                                                     datatype="datasheet", measurement_name=None)
        aut_temp_options1 = get_temp1_list

        temp_str = [f'{item:.2f}' for item in aut_temp_options1]
        self.aut_temp_dat_comboBox.clear()
        for option in temp_str:
            self.aut_temp_dat_comboBox.addItem(option)

    def temp_meas_input(self):
        """Get the measurement temperature of a particular material selected."""
        mat_text1 = self.dat_core_material_comboBox.currentText()

        get_temp1_list = []
        if mat_text1:
            get_temp1_list = database.drop_down_list(material_name=mat_text1, comparison_type="mvm", temperature=True,
                                                     datatype=fmt.MeasurementDataType.ComplexPermeability, measurement_name=fmt.MeasurementSetup.LEA_LK)
        aut_temp_options1 = get_temp1_list

        temp_str = [f'{item:.2f}' for item in aut_temp_options1]
        self.aut_temp_meas_comboBox.clear()
        for option in temp_str:
            self.aut_temp_meas_comboBox.addItem(option)

    #  **************************** Manual design tab initializations ***********************************************  #

    def md_initialize_controls(self) -> None:
        """
        Initialize the comboboxes with pre-defined values.

        :return: None
        :rtype: None
        """
        md_simulation_type_options = [self.translation_dict['inductor'], self.translation_dict['transformer']]
        md_core_material_options = ['N95', 'N49', 'N87']
        md_winding_material_options = [key for key in fmt.wire_material_database()]
        md_winding_type_options = [self.translation_dict['litz'], self.translation_dict['solid']]
        md_implicit_litz_options = [self.translation_dict["implicit_litz_radius"], self.translation_dict["implicit_ff"],
                                    self.translation_dict['implicit_strands_number']]
        md_air_gap_method_options = [self.translation_dict["percent"], self.translation_dict['manually']]
        md_air_gap_counts_options = ['0', '1', '2', '3', '4', '5']
        md_winding_scheme_options = [self.translation_dict["square"], self.translation_dict["hexa"]]
        md_core_geometry_options = [core_geometry for core_geometry in fmt.core_database()]
        md_core_geometry_options.insert(0, 'Manual')

        md_winding_litz_material_options = [key for key in fmt.litz_database()]
        md_winding_litz_material_options.insert(0, 'Manual')

        for option in md_simulation_type_options:
            self.md_simulation_type_comboBox.addItem(option)
        for option in md_core_material_options:
            self.md_core_material_comboBox.addItem(option)
        for option in md_winding_material_options:
            self.md_winding1_material_comboBox.addItem(option)
            self.md_winding2_material_comboBox.addItem(option)
        for option in md_winding_type_options:
            self.md_winding1_type_comboBox.addItem(option)
            self.md_winding2_type_comboBox.addItem(option)
        for option in md_implicit_litz_options:
            self.md_winding1_implicit_litz_comboBox.addItem(option)
            self.md_winding2_implicit_litz_comboBox.addItem(option)
        for option in md_air_gap_method_options:
            self.md_air_gap_placement_method_comboBox.addItem(option)
        for option in md_air_gap_counts_options:
            self.md_air_gap_count_comboBox.addItem(option)
        for option in md_winding_scheme_options:
            self.md_winding1_scheme_comboBox.addItem(option)
            self.md_winding2_scheme_comboBox.addItem(option)
        for option in md_winding_litz_material_options:
            self.md_winding1_litz_material_comboBox.addItem(option)
            self.md_winding2_litz_material_comboBox.addItem(option)
        for option in md_core_geometry_options:
            self.md_core_geometry_comboBox.addItem(option)

    # ----------------------------------------------------------
    # Definition tab
    # ----------------------------------------------------------
    def md_change_simulation_type(self, simulation_type_from_combo_box: str) -> None:
        """
        Change simulation type when signal of md_simulation_type_comboBox text has changed.

        Action will be enabling / disabling user inputs for not-used windings.

        :param simulation_type_from_combo_box:
        :type simulation_type_from_combo_box: str
        :return: None
        :rtype: None
        """
        if simulation_type_from_combo_box == self.translation_dict['inductor']:
            self.md_winding2_enable(False)
            # enable Inductance Value for an Inductor case
            self.groupBox_7.setVisible(True)

        elif simulation_type_from_combo_box == self.translation_dict['transformer']:
            # set winding definitions of winding 2 to editable
            self.md_winding2_enable(True)
            # disable Inductance Value for a transformer case
            self.groupBox_7.setVisible(False)

        elif simulation_type_from_combo_box == self.translation_dict['integrated transformer']:
            # set winding definitions of winding 2 to editable
            self.md_winding2_enable(True)

    def on_insulation_state_changed(self):
        """Update the flag_insulation based on the checkbox state."""
        self.flag_insulation = self.enable_insulation_checkbox.isChecked()

    def md_winding2_enable(self, status: bool) -> None:
        """
        Enable/disable all fields being in contact with winding 2.

        :param status: True / False
        :type status: bool
        :return: None
        :rtype: None
        """
        # set winding definitions of winding 2 (enable and visible)
        self.md_winding2_material_comboBox.setEnabled(status)
        self.md_winding2_type_comboBox.setEnabled(status)
        self.md_winding2_turns_lineEdit.setEnabled(status)
        self.md_winding2_implicit_litz_comboBox.setEnabled(status)
        self.md_winding2_strands_lineEdit.setEnabled(status)
        self.md_winding2_radius_lineEdit.setEnabled(False)
        self.md_winding2_fill_factor_lineEdit.setEnabled(status)
        self.md_winding2_strand_radius_lineEdit.setEnabled(status)
        self.md_winding2_litz_material_comboBox.setEnabled(status)
        self.md_winding2_groupBox.setVisible(status)

        # Set turns of winding 2 (enable and visible)
        self.md_winding2_turns_lineEdit.setVisible(status)
        self.md_winding2_scheme_comboBox.setVisible(status)
        self.md_winding2_turns_label.setVisible(status)
        self.md_winding2_scheme_label.setVisible(status)

        # set current shapes of winding 2 (enable and visible)
        self.md_winding2_current_groupBox.setVisible(status)

        # set isolation of winding 2 (enable and visible)
        self.md_isolation_s2s_lineEdit.setEnabled(status)
        self.md_isolation_p2s_lineEdit.setEnabled(status)
        self.md_isolation_s2s_lineEdit.setVisible(status)
        self.md_isolation_p2s_lineEdit.setVisible(status)
        # enable s2p only for winding 2
        self.md_isolation_s2p_lineEdit.setVisible(status)
        self.md_isolation_s2s_label.setVisible(status)
        self.md_isolation_p2s_label.setVisible(status)
        self.md_isolation_s2p_label.setVisible(status)

    def md_gmsh_pre_visualisation(self):
        """Pre-visualize when update preview button is pressed in the definitions tab."""
        geo = self.md_setup_geometry()
        print(f"geo:{geo}")

        # geo.create_model(freq=100000, visualize_before=False, do_meshing=False, save_png=True)
        geo.create_model(freq=comma_str_to_point_float(self.md_base_frequency_lineEdit.text()), pre_visualize_geometry=False,
                         save_png=True)
        print(f"geo.file_data.hybrid_color_visualize_file: {geo.file_data.hybrid_color_visualize_file}")
        image_pre_visualisation = PIL.Image.open(geo.file_data.hybrid_color_visualize_file)

        px = image_pre_visualisation.load()
        image_width, image_height = image_pre_visualisation.size

        for cut_x_left in range(0, image_width):
            if px[cut_x_left, image_height / 2] != (255, 255, 255):
                cut_x_left -= 1
                break
        for cut_x_right in reversed(range(0, image_width)):
            if px[cut_x_right, image_height / 2] != (255, 255, 255):
                cut_x_right += 1
                break

        for cut_y_bot in reversed(range(0, image_height)):
            if px[image_width / 2, cut_y_bot] != (255, 255, 255):
                cut_y_bot += 1
                break
        for cut_y_top in range(0, image_height):
            if px[image_width / 2, cut_y_top] != (255, 255, 255):
                cut_y_top -= 1
                break

        im_crop = image_pre_visualisation.crop((cut_x_left, cut_y_top, cut_x_right, cut_y_bot))
        im_crop.save(geo.file_data.hybrid_color_visualize_file, quality=95)

        pixmap = QPixmap(geo.file_data.hybrid_color_visualize_file)
        self.md_gmsh_visualisation_QLabel.setPixmap(pixmap)
        self.md_gmsh_visualisation_QLabel.setMask(pixmap.mask())
        self.md_gmsh_visualisation_QLabel.show()

    def md_set_core_geometry_from_database(self):
        """Set a core geometry for the material database."""
        core_dict = fmt.core_database()
        core_type = self.md_core_geometry_comboBox.currentText()

        if core_type != 'Manual':
            core = core_dict[core_type]

            self.md_core_width_lineEdit.setText(str(core["core_inner_diameter"]))
            self.md_window_height_lineEdit.setText(str(core["window_h"]))
            self.md_window_width_lineEdit.setText(str(core["window_w"]))
            self.md_core_width_lineEdit.setEnabled(False)
            self.md_window_height_lineEdit.setEnabled(False)
            self.md_window_width_lineEdit.setEnabled(False)
        else:
            self.md_core_width_lineEdit.setEnabled(True)
            self.md_window_height_lineEdit.setEnabled(True)
            self.md_window_width_lineEdit.setEnabled(True)

    def md_winding1_set_litz_parameters_from_litz_database(self):
        """Set litz parameters from material database for winding 1 in manual design."""
        litz_dict = fmt.litz_database()
        litz_type = self.md_winding1_litz_material_comboBox.currentText()

        if litz_type != 'Manual':
            litz = litz_dict[litz_type]

            self.md_winding1_strands_lineEdit.setText(str(litz["strands_numbers"]))
            self.md_winding1_strand_radius_lineEdit.setText(str(litz["strand_radii"]))
            self.md_winding1_radius_lineEdit.setText(str(litz["conductor_radii"]))
            self.md_winding1_fill_factor_lineEdit.setText(str(litz["ff"]))

            # for key, value in enumerate(["implicit_litz_radius", "implicit_ff", 'implicit_strands_number']):
            #     if value == litz["implicit"]:
            #         self.md_winding1_implicit_litz_comboBox.setCurrentIndex(key)

            self.md_winding1_radius_lineEdit.setEnabled(False)
            self.md_winding1_implicit_litz_comboBox.setEnabled(False)
            self.md_winding1_strands_lineEdit.setEnabled(False)
            self.md_winding1_fill_factor_lineEdit.setEnabled(False)
            self.md_winding1_strand_radius_lineEdit.setEnabled(False)
        else:
            self.md_winding1_change_litz_implicit(self.md_winding1_implicit_litz_comboBox.currentText())
            self.md_winding1_implicit_litz_comboBox.setEnabled(True)

    def md_winding1_change_litz_implicit(self, implicit_type_from_combo_box: str) -> None:
        """
        Enable / Disable input parameter fields for different "implicit xyz" types in case of litz wire.

        :param implicit_type_from_combo_box: input type to implicit
        :type implicit_type_from_combo_box: str
        :return: None
        :rtype: None
        """
        if implicit_type_from_combo_box == self.translation_dict['implicit_litz_radius']:
            self.md_winding1_strands_lineEdit.setEnabled(True)
            self.md_winding1_fill_factor_lineEdit.setEnabled(True)
            self.md_winding1_strand_radius_lineEdit.setEnabled(True)
            self.md_winding1_radius_lineEdit.setEnabled(False)
        if implicit_type_from_combo_box == self.translation_dict['implicit_strands_number']:
            self.md_winding1_strands_lineEdit.setEnabled(False)
            self.md_winding1_fill_factor_lineEdit.setEnabled(True)
            self.md_winding1_strand_radius_lineEdit.setEnabled(True)
            self.md_winding1_radius_lineEdit.setEnabled(True)
        if implicit_type_from_combo_box == self.translation_dict['implicit_ff']:
            self.md_winding1_strands_lineEdit.setEnabled(True)
            self.md_winding1_fill_factor_lineEdit.setEnabled(False)
            self.md_winding1_strand_radius_lineEdit.setEnabled(True)
            self.md_winding1_radius_lineEdit.setEnabled(True)

    def md_winding1_change_wire_type(self, wire_type_from_combot_box: str) -> None:
        """
        Enable / Disable input parameter for litz/solid wire.

        :param wire_type_from_combot_box: wire type
        :type wire_type_from_combot_box: str
        :return: None
        :rtype: None
        """
        self.md_winding1_change_litz_implicit(self.md_winding1_implicit_litz_comboBox.currentText())
        if wire_type_from_combot_box == self.translation_dict['litz']:
            self.md_winding1_strands_lineEdit.setEnabled(True)
            self.md_winding1_implicit_litz_comboBox.setEnabled(True)
            self.md_winding1_fill_factor_lineEdit.setEnabled(True)
            self.md_winding1_strand_radius_lineEdit.setEnabled(True)
            self.md_winding1_radius_lineEdit.setEnabled(True)
            self.md_winding1_litz_material_comboBox.setEnabled(True)
            self.md_winding1_change_litz_implicit(self.md_winding1_implicit_litz_comboBox.currentText())

        elif wire_type_from_combot_box == self.translation_dict['solid']:
            self.md_winding1_strands_lineEdit.setEnabled(False)
            self.md_winding1_implicit_litz_comboBox.setEnabled(False)
            self.md_winding1_fill_factor_lineEdit.setEnabled(False)
            self.md_winding1_strand_radius_lineEdit.setEnabled(False)
            self.md_winding1_radius_lineEdit.setEnabled(True)
            self.md_winding1_litz_material_comboBox.setEnabled(False)

    def md_winding2_set_litz_parameters_from_litz_database(self):
        """Set litz parameters from material database for winding 2 in manual design."""
        litz_dict = fmt.litz_database()
        litz_type = self.md_winding2_litz_material_comboBox.currentText()
        if litz_type != 'Manual':
            litz = litz_dict[litz_type]

            self.md_winding2_strands_lineEdit.setText(str(litz["strands_numbers"]))
            self.md_winding2_strand_radius_lineEdit.setText(str(litz["strand_radii"]))
            self.md_winding2_radius_lineEdit.setText(str(litz["conductor_radii"]))
            self.md_winding2_fill_factor_lineEdit.setText(str(litz["ff"]))

            for key, value in enumerate(["implicit_litz_radius", "implicit_ff", 'implicit_strands_number']):
                if value == litz["implicit"]:
                    self.md_winding2_implicit_litz_comboBox.setCurrentIndex(key)

            self.md_winding2_radius_lineEdit.setEnabled(False)
            self.md_winding2_implicit_litz_comboBox.setEnabled(False)
            self.md_winding2_strands_lineEdit.setEnabled(False)
            self.md_winding2_fill_factor_lineEdit.setEnabled(False)
            self.md_winding2_strand_radius_lineEdit.setEnabled(False)

        else:
            self.md_winding2_change_litz_implicit(self.md_winding2_implicit_litz_comboBox.currentText())
            self.md_winding2_implicit_litz_comboBox.setEnabled(True)

    def md_winding2_change_litz_implicit(self, implicit_type_from_combo_box: str) -> None:
        """
        Enable / Disable input parameter fields for different "implicit xyz" types in case of litz wire:.

        :param implicit_type_from_combo_box: input type to implicit
        :type implicit_type_from_combo_box: str
        :return: None
        :rtype: None
        """
        if implicit_type_from_combo_box == self.translation_dict['implicit_litz_radius']:
            self.md_winding2_strands_lineEdit.setEnabled(True)
            self.md_winding2_fill_factor_lineEdit.setEnabled(True)
            self.md_winding2_strand_radius_lineEdit.setEnabled(True)
            self.md_winding2_radius_lineEdit.setEnabled(False)
        if implicit_type_from_combo_box == self.translation_dict['implicit_strands_number']:
            self.md_winding2_strands_lineEdit.setEnabled(False)
            self.md_winding2_fill_factor_lineEdit.setEnabled(True)
            self.md_winding2_strand_radius_lineEdit.setEnabled(True)
            self.md_winding2_radius_lineEdit.setEnabled(True)
        if implicit_type_from_combo_box == self.translation_dict['implicit_ff']:
            self.md_winding2_strands_lineEdit.setEnabled(True)
            self.md_winding2_fill_factor_lineEdit.setEnabled(False)
            self.md_winding2_strand_radius_lineEdit.setEnabled(True)
            self.md_winding2_radius_lineEdit.setEnabled(True)

    def md_winding2_change_wire_type(self, wire_type_from_combot_box: str) -> None:
        """
        Enable / Disable input parameter for litz/solid wire.

        :param wire_type_from_combot_box: wire type
        :type wire_type_from_combot_box: str
        :return: None
        :rtype: None
        """
        self.md_winding2_change_litz_implicit(self.md_winding2_implicit_litz_comboBox.currentText())
        if wire_type_from_combot_box == self.translation_dict['litz']:
            self.md_winding2_strands_lineEdit.setEnabled(True)
            self.md_winding2_implicit_litz_comboBox.setEnabled(True)
            self.md_winding2_fill_factor_lineEdit.setEnabled(True)
            self.md_winding2_strand_radius_lineEdit.setEnabled(True)
            self.md_winding2_radius_lineEdit.setEnabled(False)
            self.md_winding2_litz_material_comboBox.setEnabled(True)
            self.md_winding2_change_litz_implicit(self.md_winding2_implicit_litz_comboBox.currentText())

        elif wire_type_from_combot_box == self.translation_dict['solid']:
            self.md_winding2_strands_lineEdit.setEnabled(False)
            self.md_winding2_implicit_litz_comboBox.setEnabled(False)
            self.md_winding2_fill_factor_lineEdit.setEnabled(False)
            self.md_winding2_strand_radius_lineEdit.setEnabled(False)
            self.md_winding2_radius_lineEdit.setEnabled(True)
            self.md_winding2_litz_material_comboBox.setEnabled(False)

    def md_change_air_gap_count(self, air_gap_count_from_combo_box: str) -> None:
        """
        Set the number of editable air gap fields in dependence of the air gap count combobox.

        :param air_gap_count_from_combo_box: Number of editable air gaps
        :type air_gap_count_from_combo_box: str
        :return: None
        :rtype: None
        """
        self.md_air_gap_placement_method_comboBox.setEnabled(False)
        self.md_air_gap_1_enable(False)
        self.md_air_gap_2_enable(False)
        self.md_air_gap_3_enable(False)
        self.md_air_gap_4_enable(False)
        self.md_air_gap_5_enable(False)
        if int(air_gap_count_from_combo_box) >= 1:
            self.md_air_gap_1_enable(True)
            self.md_air_gap_placement_method_comboBox.setEnabled(True)
        if int(air_gap_count_from_combo_box) >= 2:
            self.md_air_gap_2_enable(True)
        if int(air_gap_count_from_combo_box) >= 3:
            self.md_air_gap_3_enable(True)
        if int(air_gap_count_from_combo_box) >= 4:
            self.md_air_gap_4_enable(True)
        if int(air_gap_count_from_combo_box) >= 5:
            self.md_air_gap_5_enable(True)

    def md_change_air_gap_placement(self, air_gap_placement_from_combo_box: str) -> None:
        """
        Change the labels in case of different air gap placement methods.

        :param air_gap_placement_from_combo_box: Air gap placement method
        :type air_gap_placement_from_combo_box: str
        :return: None
        :rtype: None
        """
        if air_gap_placement_from_combo_box == self.translation_dict['manually']:
            md_air_gap_placement_text = 'Position'
        elif air_gap_placement_from_combo_box == self.translation_dict['percent']:
            md_air_gap_placement_text = 'Position in [%]'
        self.md_air_gap_1_position_label.setText(md_air_gap_placement_text)
        self.md_air_gap_2_position_label.setText(md_air_gap_placement_text)
        self.md_air_gap_3_position_label.setText(md_air_gap_placement_text)
        self.md_air_gap_4_position_label.setText(md_air_gap_placement_text)
        self.md_air_gap_5_position_label.setText(md_air_gap_placement_text)

    def md_air_gap_1_enable(self, status: bool) -> None:
        """
        Enable / Disable the input fields for air gap No. 1.

        :param status: True for enable fields, False for disable fields
        :type status: bool
        :return: None
        :rtype: None
        """
        self.md_air_gap_1_length_lineEdit.setEnabled(status)
        self.md_air_gap_1_position_lineEdit.setEnabled(status)
        self.md_air_gap_1_length_label.setVisible(status)
        self.md_air_gap_1_length_lineEdit.setVisible(status)
        self.md_air_gap_1_position_label.setVisible(status)
        self.md_air_gap_1_position_lineEdit.setVisible(status)

    def md_air_gap_2_enable(self, status: bool) -> None:
        """
        Enable / Disable the input fields for air gap No. 2.

        :param status: True for enable fields, False for disable fields
        :type status: bool
        :return: None
        :rtype: None
        """
        self.md_air_gap_2_length_lineEdit.setEnabled(status)
        self.md_air_gap_2_position_lineEdit.setEnabled(status)
        self.md_air_gap_2_length_label.setVisible(status)
        self.md_air_gap_2_length_lineEdit.setVisible(status)
        self.md_air_gap_2_position_label.setVisible(status)
        self.md_air_gap_2_position_lineEdit.setVisible(status)

    def md_air_gap_3_enable(self, status: bool) -> None:
        """
        Enable / Disable the input fields for air gap No. 3.

        :param status: True for enable fields, False for disable fields
        :type status: bool
        :return: None
        :rtype: None
        """
        self.md_air_gap_3_length_lineEdit.setEnabled(status)
        self.md_air_gap_3_position_lineEdit.setEnabled(status)
        self.md_air_gap_3_length_label.setVisible(status)
        self.md_air_gap_3_length_lineEdit.setVisible(status)
        self.md_air_gap_3_position_label.setVisible(status)
        self.md_air_gap_3_position_lineEdit.setVisible(status)

    def md_air_gap_4_enable(self, status: bool) -> None:
        """
        Enable / Disable the input fields for air gap No. 4.

        :param status: True for enable fields, False for disable fields
        :type status: bool
        :return: None
        :rtype: None
        """
        self.md_air_gap_4_length_lineEdit.setEnabled(status)
        self.md_air_gap_4_position_lineEdit.setEnabled(status)
        self.md_air_gap_4_length_label.setVisible(status)
        self.md_air_gap_4_length_lineEdit.setVisible(status)
        self.md_air_gap_4_position_label.setVisible(status)
        self.md_air_gap_4_position_lineEdit.setVisible(status)

    def md_air_gap_5_enable(self, status: bool) -> None:
        """
        Enable / Disable the input fields for air gap No. 5.

        :param status: True for enable fields, False for disable fields
        :type status: bool
        :return: None
        :rtype: None
        """
        self.md_air_gap_5_length_lineEdit.setEnabled(status)
        self.md_air_gap_5_position_lineEdit.setEnabled(status)
        self.md_air_gap_5_length_label.setVisible(status)
        self.md_air_gap_5_length_lineEdit.setVisible(status)
        self.md_air_gap_5_position_label.setVisible(status)
        self.md_air_gap_5_position_lineEdit.setVisible(status)

    # ----------------------------------------------------------
    # Excitation tab
    # ----------------------------------------------------------

    def md_dc_enable(self, status: bool) -> None:
        """
        Enable / Disable the input fields for dc current.

        :param status: True for enable fields, False for disable fields
        :type status: bool
        :return: None
        :rtype: None
        """
        self.md_winding1_idc_lineEdit.setEnabled(status)
        self.md_winding2_idc_lineEdit.setEnabled(status)
        if self.md_simulation_type_comboBox.currentText() == self.translation_dict['inductor'] and status:
            self.md_winding2_idc_lineEdit.setEnabled(False)
            self.md_winding2_idc_lineEdit.setEnabled(False)
        else:
            self.md_winding2_idc_lineEdit.setEnabled(status)
            self.md_winding2_idc_lineEdit.setEnabled(status)

    def md_f1_enable(self, status: bool) -> None:
        """
        Enable / Disable the input fields for frequency/phase No. 1.

        :param status: True for enable fields, False for disable fields
        :type status: bool
        :return: None
        :rtype: None
        """
        self.md_winding1_ik1_lineEdit.setEnabled(status)
        self.md_winding1_pk1_lineEdit.setEnabled(status)
        if self.md_simulation_type_comboBox.currentText() == self.translation_dict['inductor'] and status:
            self.md_winding2_ik1_lineEdit.setEnabled(False)
            self.md_winding2_pk1_lineEdit.setEnabled(False)
        else:
            self.md_winding2_ik1_lineEdit.setEnabled(status)
            self.md_winding2_pk1_lineEdit.setEnabled(status)

    def md_f2_enable(self, status: bool) -> None:
        """
        Enable / Disable the input fields for frequency/phase No. 2.

        :param status: True for enable fields, False for disable fields
        :type status: bool
        :return: None
        :rtype: None
        """
        self.md_winding1_ik2_lineEdit.setEnabled(status)
        self.md_winding1_pk2_lineEdit.setEnabled(status)
        if self.md_simulation_type_comboBox.currentText() == self.translation_dict['inductor'] and status:
            self.md_winding2_ik2_lineEdit.setEnabled(False)
            self.md_winding2_pk2_lineEdit.setEnabled(False)
        else:
            self.md_winding2_ik2_lineEdit.setEnabled(status)
            self.md_winding2_pk2_lineEdit.setEnabled(status)

    def md_f3_enable(self, status: bool) -> None:
        """
        Enable / Disable the input fields for frequency/phase No. 3.

        :param status: True for enable fields, False for disable fields
        :type status: bool
        :return: None
        :rtype: None
        """
        self.md_winding1_ik3_lineEdit.setEnabled(status)
        self.md_winding1_pk3_lineEdit.setEnabled(status)
        if self.md_simulation_type_comboBox.currentText() == self.translation_dict['inductor'] and status:
            self.md_winding2_ik3_lineEdit.setEnabled(False)
            self.md_winding2_pk3_lineEdit.setEnabled(False)
        else:
            self.md_winding2_ik3_lineEdit.setEnabled(status)
            self.md_winding2_pk3_lineEdit.setEnabled(status)

    def md_f4_enable(self, status: bool) -> None:
        """
        Enable / Disable the input fields for frequency/phase No. 4.

        :param status: True for enable fields, False for disable fields
        :type status: bool
        :return: None
        :rtype: None
        """
        self.md_winding1_ik4_lineEdit.setEnabled(status)
        self.md_winding1_pk4_lineEdit.setEnabled(status)
        if self.md_simulation_type_comboBox.currentText() == self.translation_dict['inductor'] and status:
            self.md_winding2_ik4_lineEdit.setEnabled(False)
            self.md_winding2_pk4_lineEdit.setEnabled(False)
        else:
            self.md_winding2_ik4_lineEdit.setEnabled(status)
            self.md_winding2_pk4_lineEdit.setEnabled(status)

    def md_f5_enable(self, status: bool) -> None:
        """
        Enable / Disable the input fields for frequency/phase No. 5.

        :param status: True for enable fields, False for disable fields
        :type status: bool
        :return: None
        :rtype: None
        """
        self.md_winding1_ik5_lineEdit.setEnabled(status)
        self.md_winding1_pk5_lineEdit.setEnabled(status)
        if self.md_simulation_type_comboBox.currentText() == self.translation_dict['inductor'] and status:
            self.md_winding2_ik5_lineEdit.setEnabled(False)
            self.md_winding2_pk5_lineEdit.setEnabled(False)
        else:
            self.md_winding2_ik5_lineEdit.setEnabled(status)
            self.md_winding2_pk5_lineEdit.setEnabled(status)

    def md_f6_enable(self, status: bool) -> None:
        """
        Enable / Disable the input fields for frequency/phase No. 6.

        :param status: True for enable fields, False for disable fields
        :type status: bool
        :return: None
        :rtype: None
        """
        self.md_winding1_ik6_lineEdit.setEnabled(status)
        self.md_winding1_pk6_lineEdit.setEnabled(status)
        if self.md_simulation_type_comboBox.currentText() == self.translation_dict['inductor'] and status:
            self.md_winding2_ik6_lineEdit.setEnabled(False)
            self.md_winding2_pk6_lineEdit.setEnabled(False)
        else:
            self.md_winding2_ik6_lineEdit.setEnabled(status)
            self.md_winding2_pk6_lineEdit.setEnabled(status)

    def md_f7_enable(self, status: bool) -> None:
        """
        Enable / Disable the input fields for frequency/phase No. 7.

        :param status: True for enable fields, False for disable fields
        :type status: bool
        :return: None
        :rtype: None
        """
        self.md_winding1_ik7_lineEdit.setEnabled(status)
        self.md_winding1_pk7_lineEdit.setEnabled(status)
        if self.md_simulation_type_comboBox.currentText() == self.translation_dict['inductor'] and status:
            self.md_winding2_ik7_lineEdit.setEnabled(False)
            self.md_winding2_pk7_lineEdit.setEnabled(False)
        else:
            self.md_winding2_ik7_lineEdit.setEnabled(status)
            self.md_winding2_pk7_lineEdit.setEnabled(status)

    def md_f8_enable(self, status: bool) -> None:
        """
        Enable / Disable the input fields for frequency/phase No. 8.

        :param status: True for enable fields, False for disable fields
        :type status: bool
        :return: None
        :rtype: None
        """
        self.md_winding1_ik8_lineEdit.setEnabled(status)
        self.md_winding1_pk8_lineEdit.setEnabled(status)
        if self.md_simulation_type_comboBox.currentText() == self.translation_dict['inductor'] and status:
            self.md_winding2_ik8_lineEdit.setEnabled(False)
            self.md_winding2_pk8_lineEdit.setEnabled(False)
        else:
            self.md_winding2_ik8_lineEdit.setEnabled(status)
            self.md_winding2_pk8_lineEdit.setEnabled(status)

    def md_change_frequencies_dc(self, status: int) -> None:
        """
        Change the frequency field in case of checking/unchecking the frequency-checkboxes.

        :param status: 0 for disabling, anything else for enabling freqency boxes
        :type status: int
        :return: None
        :rtype: None
        """
        self.md_dc_enable(False) if status == 0 else self.md_dc_enable(True)

    def md_change_frequencies_1(self, status: int) -> None:
        """
        Change the frequency field in case of checking/unchecking the frequency-checkboxes.

        :param status: 0 for disabling, anything else for enabling freqency boxes
        :type status: int
        :return: None
        :rtype: None
        """
        self.md_f1_enable(False) if status == 0 else self.md_f1_enable(True)

    def md_change_frequencies_2(self, status) -> None:
        """
        Change the frequency field in case of checking/unchecking the frequency-checkboxes.

        :param status: 0 for disabling, anything else for enabling freqency boxes
        :type status: int
        :return: None
        :rtype: None
        """
        self.md_f2_enable(False) if status == 0 else self.md_f2_enable(True)

    def md_change_frequencies_3(self, status) -> None:
        """
        Change the frequency field in case of checking/unchecking the frequency-checkboxes.

        :param status: 0 for disabling, anything else for enabling freqency boxes
        :type status: int
        :return: None
        :rtype: None
        """
        self.md_f3_enable(False) if status == 0 else self.md_f3_enable(True)

    def md_change_frequencies_4(self, status) -> None:
        """
        Change the frequency field in case of checking/unchecking the frequency-checkboxes.

        :param status: 0 for disabling, anything else for enabling freqency boxes
        :type status: int
        :return: None
        :rtype: None
        """
        self.md_f4_enable(False) if status == 0 else self.md_f4_enable(True)

    def md_change_frequencies_5(self, status) -> None:
        """
        Change the frequency field in case of checking/unchecking the frequency-checkboxes.

        :param status: 0 for disabling, anything else for enabling freqency boxes
        :type status: int
        :return: None
        :rtype: None
        """
        self.md_f5_enable(False) if status == 0 else self.md_f5_enable(True)

    def md_change_frequencies_6(self, status) -> None:
        """
        Change the frequency field in case of checking/unchecking the frequency-checkboxes.

        :param status: 0 for disabling, anything else for enabling freqency boxes
        :type status: int
        :return: None
        :rtype: None
        """
        self.md_f6_enable(False) if status == 0 else self.md_f6_enable(True)

    def md_change_frequencies_7(self, status) -> None:
        """
        Change the frequency field in case of checking/unchecking the frequency-checkboxes.

        :param status: 0 for disabling, anything else for enabling freqency boxes
        :type status: int
        :return: None
        :rtype: None
        """
        self.md_f7_enable(False) if status == 0 else self.md_f7_enable(True)

    def md_change_frequencies_8(self, status) -> None:
        """
        Change the frequency field in case of checking/unchecking the frequency-checkboxes.

        :param status: 0 for disabling, anything else for enabling freqency boxes
        :type status: int
        :return: None
        :rtype: None
        """
        self.md_f8_enable(False) if status == 0 else self.md_f8_enable(True)

    def md_redraw_input_signals(self) -> None:
        """
        Generate visual graphics for the input signals.

        Generates a graphic. This graphic is read and insertet to the gui.

        :return: None
        :rtype: None
        """
        winding1_frequency_list, winding1_amplitude_list, winding1_phi_rad_list, winding2_frequency_list, winding2_amplitude_list, \
            winding2_phi_rad_list = self.md_get_frequency_lists()

        fmt.plot_fourier_coefficients(winding1_frequency_list, winding1_amplitude_list, winding1_phi_rad_list,
                                      figure_directory="./md_winding_1.png")
        pixmap = QPixmap("./md_winding_1.png")
        self.md_graphic_winding_1.setPixmap(pixmap)
        self.md_graphic_winding_1.setMask(pixmap.mask())
        self.md_graphic_winding_1.show()
        if self.md_simulation_type_comboBox.currentText() != self.translation_dict['inductor']:
            fmt.plot_fourier_coefficients(winding2_frequency_list, winding2_amplitude_list, winding2_phi_rad_list,
                                          figure_directory="./md_winding_2.png")
            pixmap = QPixmap("./md_winding_2.png")
            self.md_graphic_winding_2.setPixmap(pixmap)
            self.md_graphic_winding_2.setMask(pixmap.mask())
            self.md_graphic_winding_2.show()

    # ----------------------------------------------------------
    # Simulation tab
    # ----------------------------------------------------------

    def md_get_frequency_lists(self) -> List:
        """
        Read frequency, amplitude and phase depending on the checked frequencies and return it as a list.

        :return: [winding1_frequency_list, winding1_amplitude_list, winding1_phi_rad_list, winding2_frequency_list,
             winding2_amplitude_list, winding2_phi_rad_list]
        :rtype: List
        """
        winding1_frequency_list = []
        winding1_amplitude_list = []
        winding1_phi_rad_list = []

        winding2_frequency_list = []
        winding2_amplitude_list = []
        winding2_phi_rad_list = []

        if self.md_fk1_checkBox.isChecked():
            winding1_frequency_list.append(1 * comma_str_to_point_float(self.md_base_frequency_lineEdit.text()))
            winding1_amplitude_list.append(comma_str_to_point_float(self.md_winding1_ik1_lineEdit.text()))
            winding1_phi_rad_list.append(comma_str_to_point_float(self.md_winding1_pk1_lineEdit.text()))
            if self.md_simulation_type_comboBox.currentText() != self.translation_dict['inductor']:
                winding2_frequency_list.append(1 * comma_str_to_point_float(self.md_base_frequency_lineEdit.text()))
                winding2_amplitude_list.append(comma_str_to_point_float(self.md_winding2_ik1_lineEdit.text()))
                winding2_phi_rad_list.append(comma_str_to_point_float(self.md_winding2_pk1_lineEdit.text()))
        if self.md_fk2_checkBox.isChecked():
            winding1_frequency_list.append(2 * comma_str_to_point_float(self.md_base_frequency_lineEdit.text()))
            winding1_amplitude_list.append(comma_str_to_point_float(self.md_winding1_ik2_lineEdit.text()))
            winding1_phi_rad_list.append(comma_str_to_point_float(self.md_winding1_pk2_lineEdit.text()))
            if self.md_simulation_type_comboBox.currentText() != self.translation_dict['inductor']:
                winding2_frequency_list.append(2 * comma_str_to_point_float(self.md_base_frequency_lineEdit.text()))
                winding2_amplitude_list.append(comma_str_to_point_float(self.md_winding2_ik2_lineEdit.text()))
                winding2_phi_rad_list.append(comma_str_to_point_float(self.md_winding2_pk2_lineEdit.text()))
        if self.md_fk3_checkBox.isChecked():
            winding1_frequency_list.append(3 * comma_str_to_point_float(self.md_base_frequency_lineEdit.text()))
            winding1_amplitude_list.append(comma_str_to_point_float(self.md_winding1_ik3_lineEdit.text()))
            winding1_phi_rad_list.append(comma_str_to_point_float(self.md_winding1_pk3_lineEdit.text()))
            if self.md_simulation_type_comboBox.currentText() != self.translation_dict['inductor']:
                winding2_frequency_list.append(3 * comma_str_to_point_float(self.md_base_frequency_lineEdit.text()))
                winding2_amplitude_list.append(comma_str_to_point_float(self.md_winding2_ik2_lineEdit.text()))
                winding2_phi_rad_list.append(comma_str_to_point_float(self.md_winding2_pk2_lineEdit.text()))
        if self.md_fk4_checkBox.isChecked():
            winding1_frequency_list.append(4 * comma_str_to_point_float(self.md_base_frequency_lineEdit.text()))
            winding1_amplitude_list.append(comma_str_to_point_float(self.md_winding1_ik4_lineEdit.text()))
            winding1_phi_rad_list.append(comma_str_to_point_float(self.md_winding1_pk4_lineEdit.text()))
            if self.md_simulation_type_comboBox.currentText() != self.translation_dict['inductor']:
                winding2_frequency_list.append(4 * comma_str_to_point_float(self.md_base_frequency_lineEdit.text()))
                winding2_amplitude_list.append(comma_str_to_point_float(self.md_winding2_ik4_lineEdit.text()))
                winding2_phi_rad_list.append(comma_str_to_point_float(self.md_winding2_pk4_lineEdit.text()))
        if self.md_fk5_checkBox.isChecked():
            winding1_frequency_list.append(5 * comma_str_to_point_float(self.md_base_frequency_lineEdit.text()))
            winding1_amplitude_list.append(comma_str_to_point_float(self.md_winding1_ik5_lineEdit.text()))
            winding1_phi_rad_list.append(comma_str_to_point_float(self.md_winding1_pk5_lineEdit.text()))
            if self.md_simulation_type_comboBox.currentText() != self.translation_dict['inductor']:
                winding2_frequency_list.append(5 * comma_str_to_point_float(self.md_base_frequency_lineEdit.text()))
                winding2_amplitude_list.append(comma_str_to_point_float(self.md_winding2_ik5_lineEdit.text()))
                winding2_phi_rad_list.append(comma_str_to_point_float(self.md_winding2_pk5_lineEdit.text()))
        if self.md_fk6_checkBox.isChecked():
            winding1_frequency_list.append(6 * comma_str_to_point_float(self.md_base_frequency_lineEdit.text()))
            winding1_amplitude_list.append(comma_str_to_point_float(self.md_winding1_ik6_lineEdit.text()))
            winding1_phi_rad_list.append(comma_str_to_point_float(self.md_winding1_pk6_lineEdit.text()))
            if self.md_simulation_type_comboBox.currentText() != self.translation_dict['inductor']:
                winding2_frequency_list.append(6 * comma_str_to_point_float(self.md_base_frequency_lineEdit.text()))
                winding2_amplitude_list.append(comma_str_to_point_float(self.md_winding2_ik6_lineEdit.text()))
                winding2_phi_rad_list.append(comma_str_to_point_float(self.md_winding2_pk6_lineEdit.text()))
        if self.md_fk7_checkBox.isChecked():
            winding1_frequency_list.append(7 * comma_str_to_point_float(self.md_base_frequency_lineEdit.text()))
            winding1_amplitude_list.append(comma_str_to_point_float(self.md_winding1_ik7_lineEdit.text()))
            winding1_phi_rad_list.append(comma_str_to_point_float(self.md_winding1_pk7_lineEdit.text()))
            if self.md_simulation_type_comboBox.currentText() != self.translation_dict['inductor']:
                winding2_frequency_list.append(7 * comma_str_to_point_float(self.md_base_frequency_lineEdit.text()))
                winding2_amplitude_list.append(comma_str_to_point_float(self.md_winding2_ik7_lineEdit.text()))
                winding2_phi_rad_list.append(comma_str_to_point_float(self.md_winding2_pk7_lineEdit.text()))
        if self.md_fk8_checkBox.isChecked():
            winding1_frequency_list.append(8 * comma_str_to_point_float(self.md_base_frequency_lineEdit.text()))
            winding1_amplitude_list.append(comma_str_to_point_float(self.md_winding1_ik8_lineEdit.text()))
            winding1_phi_rad_list.append(comma_str_to_point_float(self.md_winding1_pk8_lineEdit.text()))
            if self.md_simulation_type_comboBox.currentText() != self.translation_dict['inductor']:
                winding2_frequency_list.append(8 * comma_str_to_point_float(self.md_base_frequency_lineEdit.text()))
                winding2_amplitude_list.append(comma_str_to_point_float(self.md_winding2_ik8_lineEdit.text()))
                winding2_phi_rad_list.append(comma_str_to_point_float(self.md_winding2_pk8_lineEdit.text()))

        return winding1_frequency_list, winding1_amplitude_list, winding1_phi_rad_list, winding2_frequency_list, winding2_amplitude_list, winding2_phi_rad_list

    def wdg_scheme(self):
        """Choose the winding scheme."""
        if self.md_winding1_scheme_comboBox.currentText() == self.translation_dict["square"]:
            scheme = 'square'
        elif self.md_winding1_scheme_comboBox.currentText() == self.translation_dict["hexa"]:
            scheme = 'hexa'

    def md_setup_geometry(self):
        """
        Set up the core and conductor geometry depending on the GUI input parameters.

        returns: femmt MagneticComponent

        """
        # geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, is_gui=True)
        # self.check_onelab_config(geo)

        if self.md_simulation_type_comboBox.currentText() == self.translation_dict['inductor']:
            self.md_simulation_QLabel.setText('simulation startet...')

            # -----------------------------------------------
            # Core
            # -----------------------------------------------
            """
            geo.core.update(type="EI",
                            core_w=comma_str_to_point_float(self.md_core_width_lineEdit.text()),
                            window_h=comma_str_to_point_float(self.md_window_height_lineEdit.text()),
                            window_w=comma_str_to_point_float(self.md_window_width_lineEdit.text()))"""

            geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor,
                                        working_directory=self.md_working_directory_lineEdit.text(),
                                        verbosity=fmt.Verbosity.ToConsole)

            core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=comma_str_to_point_float(self.md_core_width_lineEdit.text()),
                                                            window_w=comma_str_to_point_float(self.md_window_width_lineEdit.text()),
                                                            window_h=comma_str_to_point_float(self.md_window_height_lineEdit.text()),
                                                            core_h=1)

            material_enum = fmt.Material(self.md_core_material_comboBox.currentText())

            core = fmt.Core(core_type=fmt.CoreType.Single,
                            core_dimensions=core_dimensions,
                            detailed_core_model=False,
                            material=material_enum, temperature=45, frequency=int(self.md_base_frequency_lineEdit.text()),
                            # permeability_datasource="manufacturer_datasheet",
                            permeability_datasource=fmt.MaterialDataSource.Measurement,
                            permeability_datatype=fmt.MeasurementDataType.ComplexPermeability,
                            permeability_measurement_setup=mdb.MeasurementSetup.LEA_LK,
                            permittivity_datasource=fmt.MaterialDataSource.Measurement,
                            permittivity_datatype=fmt.MeasurementDataType.ComplexPermittivity,
                            permittivity_measurement_setup=mdb.MeasurementSetup.LEA_LK, mdb_verbosity=fmt.Verbosity.Silent)

            geo.set_core(core)

            # -----------------------------------------------
            # Air Gaps
            # -----------------------------------------------

            air_gap_count = int(self.md_air_gap_count_comboBox.currentText())
            air_gap_heigth_array = []
            air_gap_position_array = []
            air_gap_position_tag_array = []

            if air_gap_count >= 1:
                md_air_gap_1_height = comma_str_to_point_float(self.md_air_gap_1_length_lineEdit.text())
                md_air_gap_1_position = comma_str_to_point_float(self.md_air_gap_1_position_lineEdit.text())

                air_gap_heigth_array.append(md_air_gap_1_height)
                air_gap_position_array.append(md_air_gap_1_position)
                air_gap_position_tag_array.append(0)

            if air_gap_count >= 2:
                md_air_gap_2_height = comma_str_to_point_float(self.md_air_gap_2_length_lineEdit.text())
                md_air_gap_2_position = comma_str_to_point_float(self.md_air_gap_2_position_lineEdit.text())

                air_gap_heigth_array.append(md_air_gap_2_height)
                air_gap_position_array.append(md_air_gap_2_position)
                air_gap_position_tag_array.append(0)

            if air_gap_count >= 3:
                md_air_gap_3_height = comma_str_to_point_float(self.md_air_gap_3_length_lineEdit.text())
                md_air_gap_3_position = comma_str_to_point_float(self.md_air_gap_3_position_lineEdit.text())

                air_gap_heigth_array.append(md_air_gap_3_height)
                air_gap_position_array.append(md_air_gap_3_position)
                air_gap_position_tag_array.append(0)

            if air_gap_count >= 4:
                md_air_gap_4_height = comma_str_to_point_float(self.md_air_gap_4_length_lineEdit.text())
                md_air_gap_4_position = comma_str_to_point_float(self.md_air_gap_4_position_lineEdit.text())

                air_gap_heigth_array.append(md_air_gap_4_height)
                air_gap_position_array.append(md_air_gap_4_position)
                air_gap_position_tag_array.append(0)

            if air_gap_count >= 5:
                md_air_gap_5_height = comma_str_to_point_float(self.md_air_gap_5_length_lineEdit.text())
                md_air_gap_5_position = comma_str_to_point_float(self.md_air_gap_5_position_lineEdit.text())

                air_gap_heigth_array.append(md_air_gap_5_height)
                air_gap_position_array.append(md_air_gap_5_position)
                air_gap_position_tag_array.append(0)

            if air_gap_count == 0:
                air_gaps = fmt.AirGaps(fmt.AirGapMethod.Center, core)
                # air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, None, 0.0005)
                geo.set_air_gaps(air_gaps)
                """
                geo.air_gaps.update(method="percent",
                                    n_air_gaps=0,
                                    air_gap_h=[],
                                    position_tag=[],
                                    air_gap_position=[])"""

            elif self.md_air_gap_placement_method_comboBox.currentText() == self.translation_dict["percent"] and air_gap_count >= 1:
                air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
                for i in range(1, air_gap_count + 1):
                    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, air_gap_heigth_array[i - 1],
                                         air_gap_position_array[i - 1])
                """
                air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
                air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, md_air_gap_1_position, md_air_gap_1_height)"""
                geo.set_air_gaps(air_gaps)
                """
                geo.air_gaps.update(method="percent",
                                    n_air_gaps=air_gap_count,
                                    air_gap_h=air_gap_heigth_array,
                                    position_tag=air_gap_position_tag_array,
                                    air_gap_position=air_gap_position_array)"""

            elif self.md_air_gap_placement_method_comboBox.currentText() == self.translation_dict["manually"] and air_gap_count >= 1:

                air_gaps = fmt.AirGaps(fmt.AirGapMethod.Manually, core)
                for i in range(1, air_gap_count + 1):
                    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, air_gap_heigth_array[i - 1],
                                         air_gap_position_array[i - 1])
                """
                air_gaps = fmt.AirGaps(fmt.AirGapMethod.Manually, core)
                air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, md_air_gap_1_position, md_air_gap_1_height)"""
                geo.set_air_gaps(air_gaps)
                """
                geo.air_gaps.update(method="manually",
                                    n_air_gaps=air_gap_count,
                                    air_gap_h=air_gap_heigth_array,
                                    position_tag=air_gap_position_tag_array,
                                    air_gap_position=air_gap_position_array)"""

            # ------------------------------------------------------
            # Set insulations
            # ------------------------------------------------------
            insulation = fmt.Insulation(flag_insulation=self.flag_insulation)
            insulation.add_core_insulations(comma_str_to_point_float(self.md_isolation_core2cond_top_lineEdit.text()),
                                            comma_str_to_point_float(self.md_isolation_core2cond_bot_lineEdit.text()),
                                            comma_str_to_point_float(self.md_isolation_core2cond_inner_lineEdit.text()),
                                            comma_str_to_point_float(self.md_isolation_core2cond_outer_lineEdit.text()))
            insulation.add_winding_insulations([[comma_str_to_point_float(self.md_isolation_p2p_lineEdit.text())]])
            geo.set_insulation(insulation)

            # ---------------------------------------------------------
            # Create winding window and virtual winding windows (vww)
            # ----------------------------------------------------------

            winding_window = fmt.WindingWindow(core, insulation)
            vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

            # ----------------------------------------------------------
            # Create conductor and set parameters: use solid wires
            # -----------------------------------------------------------

            winding_material_name = self.md_winding1_material_comboBox.currentText()
            if winding_material_name == 'Copper':
                winding_material_enum = fmt.Conductivity.Copper
            elif winding_material_name == 'Aluminium':
                winding_material_enum = fmt.Conductivity.Aluminium

            winding = fmt.Conductor(0, winding_material_enum)
            if self.md_winding1_type_comboBox.currentText() == self.translation_dict['solid']:
                winding.set_solid_round_conductor(
                    conductor_radius=comma_str_to_point_float(self.md_winding1_radius_lineEdit.text()),
                    conductor_arrangement=fmt.ConductorArrangement.Square)
            elif self.md_winding1_type_comboBox.currentText() == self.translation_dict['litz']:
                winding.set_litz_round_conductor(
                    conductor_radius=comma_str_to_point_float(self.md_winding1_radius_lineEdit.text()),
                    number_strands=int(self.md_winding1_strands_lineEdit.text()),
                    strand_radius=comma_str_to_point_float(self.md_winding1_strand_radius_lineEdit.text()),
                    fill_factor=None,
                    conductor_arrangement=fmt.ConductorArrangement.Square)

            # ----------------------------------------------------------------------
            # 7. add conductor to vww and add winding window to MagneticComponent
            # ----------------------------------------------------------------------
            # vww.set_winding(winding, comma_str_to_point_float(self.md_winding1_turns_lineEdit.text()), None)
            # change 1: Number of turns should be an integer, not as a float
            vww.set_winding(winding, int(self.md_winding1_turns_lineEdit.text()), None,
                            fmt.Align.ToEdges, placing_strategy=fmt.ConductorDistribution.HorizontalRightward_VerticalUpward, zigzag=True)
            geo.set_winding_windows([winding_window])

        elif self.md_simulation_type_comboBox.currentText() == 'transformer':
            self.md_simulation_QLabel.setText('simulation startet...')

            # 1. chose simulation type
            geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Transformer,
                                        working_directory=self.md_working_directory_lineEdit.text())

            # -----------------------------------------------
            # Core
            # -----------------------------------------------
            core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=comma_str_to_point_float(self.md_core_width_lineEdit.text()),
                                                            window_w=comma_str_to_point_float(self.md_window_width_lineEdit.text()),
                                                            window_h=comma_str_to_point_float(self.md_window_height_lineEdit.text()),
                                                            core_h=0.04)
            core = fmt.Core(core_dimensions=core_dimensions, mu_r_abs=3100, phi_mu_deg=12, sigma=1.2,
                            permeability_datasource=fmt.MaterialDataSource.Custom,
                            permittivity_datasource=fmt.MaterialDataSource.Custom,
                            detailed_core_model=True)
            geo.set_core(core)
            """
            geo.core.update(window_h = comma_str_to_point_float(self.md_window_height_lineEdit.text()),
                            window_w = comma_str_to_point_float(self.md_window_width_lineEdit.text()),
                            core_w = comma_str_to_point_float(self.md_core_width_lineEdit.text()),
                            mu_rel=3100, phi_mu_deg=12, sigma=0.6)"""

            # -----------------------------------------------
            # Air Gaps
            # -----------------------------------------------

            air_gap_count = int(self.md_air_gap_count_comboBox.currentText())
            air_gap_heigth_array = []
            air_gap_position_array = []
            air_gap_position_tag_array = []

            if air_gap_count >= 1:
                md_air_gap_1_height = comma_str_to_point_float(self.md_air_gap_1_length_lineEdit.text())
                md_air_gap_1_position = comma_str_to_point_float(self.md_air_gap_1_position_lineEdit.text())

                air_gap_heigth_array.append(md_air_gap_1_height)
                air_gap_position_array.append(md_air_gap_1_position)
                air_gap_position_tag_array.append(0)

            if air_gap_count >= 2:
                md_air_gap_2_height = comma_str_to_point_float(self.md_air_gap_2_length_lineEdit.text())
                md_air_gap_2_position = comma_str_to_point_float(self.md_air_gap_2_position_lineEdit.text())

                air_gap_heigth_array.append(md_air_gap_2_height)
                air_gap_position_array.append(md_air_gap_2_position)
                air_gap_position_tag_array.append(0)

            if air_gap_count >= 3:
                md_air_gap_3_height = comma_str_to_point_float(self.md_air_gap_3_length_lineEdit.text())
                md_air_gap_3_position = comma_str_to_point_float(self.md_air_gap_3_position_lineEdit.text())

                air_gap_heigth_array.append(md_air_gap_3_height)
                air_gap_position_array.append(md_air_gap_3_position)
                air_gap_position_tag_array.append(0)

            if air_gap_count >= 4:
                md_air_gap_4_height = comma_str_to_point_float(self.md_air_gap_4_length_lineEdit.text())
                md_air_gap_4_position = comma_str_to_point_float(self.md_air_gap_4_position_lineEdit.text())

                air_gap_heigth_array.append(md_air_gap_4_height)
                air_gap_position_array.append(md_air_gap_4_position)
                air_gap_position_tag_array.append(0)

            if air_gap_count >= 5:
                md_air_gap_5_height = comma_str_to_point_float(self.md_air_gap_5_length_lineEdit.text())
                md_air_gap_5_position = comma_str_to_point_float(self.md_air_gap_5_position_lineEdit.text())

                air_gap_heigth_array.append(md_air_gap_5_height)
                air_gap_position_array.append(md_air_gap_5_position)
                air_gap_position_tag_array.append(0)

            if air_gap_count == 0:

                air_gaps = fmt.AirGaps(fmt.AirGapMethod.Center, core)
                # air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, None, 0.0005)
                geo.set_air_gaps(air_gaps)
                """
                geo.air_gaps.update(method="percent",
                                    n_air_gaps=0,
                                    air_gap_h=[],
                                    position_tag=[],
                                    air_gap_position=[])"""

            elif self.md_air_gap_placement_method_comboBox.currentText() == self.translation_dict["percent"] and air_gap_count >= 1:

                air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
                for i in range(1, air_gap_count + 1):
                    # change 2
                    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, air_gap_heigth_array[i - 1],
                                         air_gap_position_array[i - 1])
                """
                air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
                air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, md_air_gap_1_position, md_air_gap_1_height)"""
                geo.set_air_gaps(air_gaps)
                """
                geo.air_gaps.update(method="percent",
                                    n_air_gaps=air_gap_count,
                                    air_gap_h=air_gap_heigth_array,
                                    position_tag=air_gap_position_tag_array,
                                    air_gap_position=air_gap_position_array)"""

            elif self.md_air_gap_placement_method_comboBox.currentText() == self.translation_dict["manually"] and air_gap_count >= 1:

                air_gaps = fmt.AirGaps(fmt.AirGapMethod.Manually, core)
                for i in range(1, air_gap_count + 1):
                    # change 2
                    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, air_gap_heigth_array[i - 1],
                                         air_gap_position_array[i - 1])
                """
                air_gaps = fmt.AirGaps(fmt.AirGapMethod.Manually, core)
                air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, md_air_gap_1_position, md_air_gap_1_height)"""
                geo.set_air_gaps(air_gaps)
                """
                geo.air_gaps.update(method="manually",
                                    n_air_gaps=air_gap_count,
                                    air_gap_h=air_gap_heigth_array,
                                    position_tag=air_gap_position_tag_array,
                                    air_gap_position=air_gap_position_array)"""

            # 4. set insulation
            insulation = fmt.Insulation(flag_insulation=self.flag_insulation)
            insulation.add_core_insulations(comma_str_to_point_float(self.md_isolation_core2cond_top_lineEdit.text()),
                                            comma_str_to_point_float(self.md_isolation_core2cond_bot_lineEdit.text()),
                                            comma_str_to_point_float(self.md_isolation_core2cond_inner_lineEdit.text()),
                                            comma_str_to_point_float(self.md_isolation_core2cond_outer_lineEdit.text()))
            insulation.add_winding_insulations([
                [comma_str_to_point_float(self.md_isolation_p2p_lineEdit.text()),
                 comma_str_to_point_float(self.md_isolation_p2s_lineEdit.text())],
                [comma_str_to_point_float(self.md_isolation_s2p_lineEdit.text()),
                 comma_str_to_point_float(self.md_isolation_s2s_lineEdit.text())]])
            geo.set_insulation(insulation)

            # 5. create winding window and virtual winding windows (vww)
            winding_window = fmt.WindingWindow(core, insulation)
            left, right = winding_window.split_window(fmt.WindingWindowSplit.HorizontalSplit)

            # 6. create conductors and set parameters
            winding1 = fmt.Conductor(0, fmt.Conductivity.Copper)
            winding1.set_solid_round_conductor(comma_str_to_point_float(self.md_winding1_radius_lineEdit.text()),
                                               fmt.ConductorArrangement.Square)

            winding2 = fmt.Conductor(1, fmt.Conductivity.Copper)
            winding2.set_solid_round_conductor(comma_str_to_point_float(self.md_winding2_radius_lineEdit.text()),
                                               fmt.ConductorArrangement.Square)

            # 7. add conductor to vww and add winding window to MagneticComponent
            left.set_winding(winding1, int(self.md_winding1_turns_lineEdit.text()), None, fmt.Align.ToEdges,
                             fmt.ConductorDistribution.VerticalUpward_HorizontalRightward, zigzag=False)
            right.set_winding(winding2, int(self.md_winding2_turns_lineEdit.text()), None, fmt.Align.ToEdges,
                              fmt.ConductorDistribution.VerticalUpward_HorizontalRightward, zigzag=False)
            geo.set_winding_windows([winding_window])

        elif self.md_simulation_type_comboBox.currentText() == 'integrated transformer':
            pass

        return geo

    def md_action_run_simulation(self) -> None:
        """
        Read all input parameters from the fields. Run the simulation in dependence of input fields.

        :return: None
        :rtype: None
        """
        geo = self.md_setup_geometry()
        # -----------------------------------------------
        # Simulation
        # -----------------------------------------------
        geo.create_model(freq=comma_str_to_point_float(self.md_base_frequency_lineEdit.text()), pre_visualize_geometry=False,
                         save_png=False)
        # geo.create_model(freq=comma_str_to_point_float(self.md_base_frequency_lineEdit.text()), visualize_before=False, do_meshing=True, save_png=False)

        winding1_frequency_list, winding1_amplitude_list, winding1_phi_rad_list, winding2_frequency_list, winding2_amplitude_list, \
            winding2_phi_rad_list = self.md_get_frequency_lists()
        print(winding1_frequency_list)
        print(winding1_amplitude_list)
        if len(winding1_frequency_list) == 1:
            if self.md_simulation_type_comboBox.currentText() == self.translation_dict['inductor']:
                geo.single_simulation(freq=winding1_frequency_list[0],
                                      current=[winding1_amplitude_list[0]],
                                      show_fem_simulation_results=True)
            elif self.md_simulation_type_comboBox.currentText() == self.translation_dict['transformer']:
                geo.single_simulation(freq=winding1_frequency_list[0],
                                      current=[winding1_amplitude_list[0], winding2_amplitude_list[0]],
                                      phi_deg=[winding1_phi_rad_list[0], winding2_phi_rad_list[0]])

        else:

            if self.md_simulation_type_comboBox.currentText() == self.translation_dict['inductor']:
                amplitude_list = []
                print(f"{winding1_amplitude_list = }")
                for amplitude_value in winding1_amplitude_list:
                    amplitude_list.append([amplitude_value])

                phase_rad_list = []
                for phase_value in winding1_phi_rad_list:
                    phase_rad_list.append([phase_value])
                geo.excitation_sweep(frequency_list=winding1_frequency_list, current_list_list=amplitude_list,
                                     phi_deg_list_list=phase_rad_list)

            elif self.md_simulation_type_comboBox.currentText() == self.translation_dict['transformer']:
                amplitude1_list = []
                for amplitude1_value, amplitude2_value in zip(winding1_amplitude_list, winding2_amplitude_list):
                    amplitude1_list.append([amplitude1_value, amplitude2_value])

                phase1_rad_list = []
                for phase1_value, phase2_value in zip(winding1_phi_rad_list, winding2_phi_rad_list):
                    phase1_rad_list.append([phase1_value, phase2_value])

                geo.excitation_sweep(frequency_list=winding1_frequency_list,
                                     current_list_list=amplitude1_list,
                                     phi_deg_list_list=phase1_rad_list)

            # geo.excitation_sweep(winding1_frequency_list, amplitude_list, phase_rad_list)

        # -----------------------------------------------
        # Read back results
        # -----------------------------------------------

        self.md_simulation_QLabel.setText('simulation fertig.')

        # loaded_results_dict = fmt.visualize_simulation_results(geo.file_data.femm_results_log_path, './results.png', show_plot=False)
        loaded_results_dict = fmt.visualize_simulation_results(geo.file_data.e_m_results_log_path,
                                                               geo.file_data.results_em_simulation, show_plot=False)
        # pixmap = QPixmap("./results.png")
        pixmap = QPixmap(geo.file_data.results_em_simulation)
        self.md_loss_plot_label.setPixmap(pixmap)
        self.md_loss_plot_label.setMask(pixmap.mask())
        self.md_loss_plot_label.show()

        # inductance = loaded_results_dict["single_sweeps"][0]["winding1"]["flux_over_current"][0]
        loss_core_eddy_current = loaded_results_dict["total_losses"]["eddy_core"]
        loss_core_hysteresis = loaded_results_dict["total_losses"]["hyst_core_fundamental_freq"]
        # loss_winding_1 = loaded_results_dict["total_losses"]["winding1"]["total"]
        self.md_loss_core_hysteresis_label.setText(f"Core Hysteresis loss: {loss_core_hysteresis} W")
        self.md_loss_core_eddy_current_label.setText(f"Core Eddy Current loss: {loss_core_eddy_current} W")
        # self.md_loss_winding1_label.setText(f"Winding 1 loss: {loss_winding_1} W")
        # self.md_inductance_label.setText(f"Primary Inductance: {inductance} H")
        inductances = []
        windings_loss = []
        for i in range(1, 3):  # for 3 windings
            winding_key = f"winding{i}"
            if winding_key in loaded_results_dict["single_sweeps"][0]:
                inductances.append(loaded_results_dict["single_sweeps"][0][winding_key]["flux_over_current"][0])
                windings_loss.append(loaded_results_dict["total_losses"][winding_key]["total"])

        # show them just for 2 windings in GUI
        if self.md_simulation_type_comboBox.currentText() == self.translation_dict['inductor']:
            self.md_loss_winding1_label.setText(f"Winding 1 loss: {windings_loss[0]} W")
            self.md_inductance1_label.setText(f"Primary Inductance: {inductances[0]} H")
        elif self.md_simulation_type_comboBox.currentText() == self.translation_dict['transformer']:
            self.md_loss_winding1_label.setText(f"Winding 1 loss: {windings_loss[0]} W")
            self.md_loss_winding2_label.setText(f"Winding 2 loss: {windings_loss[1]} W")
            self.md_inductance1_label.setText(f"Primary Inductance: {inductances[0]} H")
            self.md_inductance2_label.setText(f"Secondary Inductance: {inductances[1]} H")

        # log_path = geo.e_m_results_log_path
        # simulation_results = str(fmt.read_results_log(log_path))
        # print(simulation_results)
        # self.md_simulation_output_textBrowser.setText(simulation_results)

    def inductancecalc(self):
        """Calculate inductance from given geometries."""
        air_gap_count = int(self.md_air_gap_count_comboBox.currentText())
        air_gap_heigth_array = []
        air_gap_position_array = []
        air_gap_position_tag_array = []

        if air_gap_count >= 1:
            md_air_gap_1_height = comma_str_to_point_float(self.md_air_gap_1_length_lineEdit.text())
            md_air_gap_1_position = comma_str_to_point_float(self.md_air_gap_1_position_lineEdit.text())

            air_gap_heigth_array.append(md_air_gap_1_height)
            air_gap_position_array.append(md_air_gap_1_position)
            air_gap_position_tag_array.append(0)

        if air_gap_count >= 2:
            md_air_gap_2_height = comma_str_to_point_float(self.md_air_gap_2_length_lineEdit.text())
            md_air_gap_2_position = comma_str_to_point_float(self.md_air_gap_2_position_lineEdit.text())

            air_gap_heigth_array.append(md_air_gap_2_height)
            air_gap_position_array.append(md_air_gap_2_position)
            air_gap_position_tag_array.append(0)

        if air_gap_count >= 3:
            md_air_gap_3_height = comma_str_to_point_float(self.md_air_gap_3_length_lineEdit.text())
            md_air_gap_3_position = comma_str_to_point_float(self.md_air_gap_3_position_lineEdit.text())

            air_gap_heigth_array.append(md_air_gap_3_height)
            air_gap_position_array.append(md_air_gap_3_position)
            air_gap_position_tag_array.append(0)

        if air_gap_count >= 4:
            md_air_gap_4_height = comma_str_to_point_float(self.md_air_gap_4_length_lineEdit.text())
            md_air_gap_4_position = comma_str_to_point_float(self.md_air_gap_4_position_lineEdit.text())

            air_gap_heigth_array.append(md_air_gap_4_height)
            air_gap_position_array.append(md_air_gap_4_position)
            air_gap_position_tag_array.append(0)

        if air_gap_count >= 5:
            md_air_gap_5_height = comma_str_to_point_float(self.md_air_gap_5_length_lineEdit.text())
            md_air_gap_5_position = comma_str_to_point_float(self.md_air_gap_5_position_lineEdit.text())

            air_gap_heigth_array.append(md_air_gap_5_height)
            air_gap_position_array.append(md_air_gap_5_position)
            air_gap_position_tag_array.append(0)

        self.core_w = comma_str_to_point_float(self.md_core_width_lineEdit.text())
        self.window_w = comma_str_to_point_float(self.md_window_width_lineEdit.text())
        self.window_h = comma_str_to_point_float(self.md_window_height_lineEdit.text())
        n_turns = int(self.md_winding1_turns_lineEdit.text())
        method = (self.md_air_gap_placement_method_comboBox.currentText())
        # murel = database.get_initial_permeability(self.md_core_material_comboBox.currentText())
        # murel = database.get_initial_permeability("N95")
        air_gap_h = self.md_air_gap_1_length_lineEdit.text()
        air_gap_position = self.md_air_gap_1_position_lineEdit.text()

        # material_names = []
        # material_names.append(self.md_core_material_comboBox.currentText())
        # mu_rel_val = [database.get_material_attribute(material_name=material_name, attribute="initial_permeability")
        #               for material_name in material_names]
        # mu_rel = [int(item) for item in mu_rel_val]

        # Fetch the current material name from the combobox
        material_name_str = self.md_core_material_comboBox.currentText()
        # Convert the string to enum
        material_enum = self.string_to_material_enum(material_name_str)

        # initial permeability from the database using the enum material
        # No need for a loop as the user can choose only one material in manual design
        mu_rel_val = database.get_material_attribute(material_name=material_enum, attribute="initial_permeability")
        mu_rel = [int(mu_rel_val)]

        print(f"core_inner_diameter: {[self.core_w]}")
        print(f"window_h: {[self.window_h]}")
        print(f"window_w: {[self.window_w]}")
        print(f"no_of_turns: {[n_turns]}")
        print(f"n_air_gaps: {[air_gap_count]}")
        print(f"air_gap_h: {air_gap_heigth_array}")
        print(f"air_gap_position: {air_gap_position_array}")
        print(f"mu_rel: {mu_rel}")
        print(f"component_type: {self.md_simulation_type_comboBox.currentText()}")

        # mc1 = fmt.MagneticCircuit([self.core_w], [self.window_h], [self.window_w], [n_turns], [n_air_gaps],
        # [air_gap_h], [air_gap_position], [3000], [1]) #3000 - relative permeability of selected material

        if self.md_air_gap_placement_method_comboBox.currentText() == self.translation_dict['percent']:
            mc1 = fmt.MagneticCircuit(core_inner_diameter=[self.core_w], window_h=[self.window_h],
                                      window_w=[self.window_w], no_of_turns=[n_turns],
                                      n_air_gaps=[air_gap_count], air_gap_h=air_gap_heigth_array,
                                      air_gap_position=air_gap_position_array,
                                      mu_r_abs=mu_rel, mult_air_gap_type=[1, 2], air_gap_method='Percent',
                                      component_type=self.md_simulation_type_comboBox.currentText(),
                                      sim_type='single')  # 0.0149
        elif self.md_air_gap_placement_method_comboBox.currentText() == self.translation_dict['manually']:
            mc1 = fmt.MagneticCircuit(core_inner_diameter=[self.core_w], window_h=[self.window_h],
                                      window_w=[self.window_w],
                                      no_of_turns=[n_turns],
                                      n_air_gaps=[air_gap_count], air_gap_h=air_gap_heigth_array,
                                      air_gap_position=air_gap_position_array,
                                      mu_r_abs=mu_rel, mult_air_gap_type=[1, 2], air_gap_method='manually',
                                      component_type=self.md_simulation_type_comboBox.currentText(), sim_type='single')

        mc1.calculate_inductance()
        inductance = mc1.data_matrix[:, 9]

        self.Inductanceval_label.setText(f"{round(inductance[0], 10)} H")

    def therm_simulation(self):
        """Implement the thermal simulation."""
        # Thermal simulation:
        # The losses calculated by the magnetics simulation can be used to calculate the heat distribution of the given magnetic component
        # In order to use the thermal simulation, thermal conductivities for each material can be entered as well as a boundary temperature
        # which will be applied on the boundary of the simulation (dirichlet boundary condition).

        # The case parameter sets the thermal conductivity for a case which will be set around the core.
        # This could model some case in which the transformer is placed in together with a set potting material.
        # flag_insulation = True
        thermal_conductivity_dict = {
            "air": 0.0263,
            "case": {  # (epoxy resign) | transformer oil
                "top": 0.122,
                "top_right": 0.122,
                "right": 0.122,
                "bot_right": 0.122,
                "bot": 0.122
            },
            "core": 5,  # ferrite
            "winding": 400,  # copper
            "air_gaps": 180,  # aluminiumnitride
            "insulation": 0.42 if self.flag_insulation else None  # polyethylen
        }
        # Here the case size can be determined
        case_gap_top = 0.002
        case_gap_right = 0.0025
        case_gap_bot = 0.002
        # Here the boundary temperatures can be set, currently it is set to 20°C (around 293°K).
        # This does not change the results of the simulation (at least when every boundary is set equally) but will set the temperature offset.
        boundary_temperatures = {
            "value_boundary_top": 20,
            "value_boundary_top_right": 20,
            "value_boundary_right_top": 20,
            "value_boundary_right": 20,
            "value_boundary_right_bottom": 20,
            "value_boundary_bottom_right": 20,
            "value_boundary_bottom": 20
        }

        # In order to compare the femmt thermal simulation with a femm heat flow simulation the same boundary temperature should be applied.
        # Currently only one temperature can be applied which will be set on every boundary site.
        femm_boundary_temperature = 20

        # Here the boundary sides can be turned on (1) or off (0)
        # By turning off the flag a neumann boundary will be applied at this point with heat flux = 0
        boundary_flags = {
            "flag_boundary_top": 0,
            "flag_boundary_top_right": 0,
            "flag_boundary_right_top": 1,
            "flag_boundary_right": 1,
            "flag_boundary_right_bottom": 1,
            "flag_boundary_bottom_right": 1,
            "flag_boundary_bottom": 1
        }

        # In order for the thermal simulation to work an electro_magnetic simulation has to run before.
        # The em-simulation will create a file containing the losses.
        # When the losses file is already created and contains the losses for the current model, it is enough to run geo.create_model in
        # order for the thermal simulation to work (geo.single_simulation is not needed).
        # Obviously when the model is modified and the losses can be out of date and therefore the geo.single_simulation needs to run again.
        geo = self.md_setup_geometry()
        geo.create_model(freq=comma_str_to_point_float(self.md_base_frequency_lineEdit.text()), pre_visualize_geometry=False,
                         save_png=False)
        # electromagnetic simulation is needed for thermal simulation
        winding1_frequency_list, winding1_amplitude_list, winding1_phi_rad_list, winding2_frequency_list, winding2_amplitude_list, \
            winding2_phi_rad_list = self.md_get_frequency_lists()
        if len(winding1_frequency_list) == 1:
            if self.md_simulation_type_comboBox.currentText() == self.translation_dict['inductor']:
                geo.single_simulation(freq=winding1_frequency_list[0],
                                      current=[winding1_amplitude_list[0]],
                                      show_fem_simulation_results=True)
            elif self.md_simulation_type_comboBox.currentText() == self.translation_dict['transformer']:
                geo.single_simulation(freq=winding1_frequency_list[0],
                                      current=[winding1_amplitude_list[0], winding2_amplitude_list[0]],
                                      phi_deg=[winding1_phi_rad_list[0], winding2_phi_rad_list[0]])

        else:

            if self.md_simulation_type_comboBox.currentText() == self.translation_dict['inductor']:
                amplitude_list = []
                print(f"{winding1_amplitude_list = }")
                for amplitude_value in winding1_amplitude_list:
                    amplitude_list.append([amplitude_value])

                phase_rad_list = []
                for phase_value in winding1_phi_rad_list:
                    phase_rad_list.append([phase_value])
                geo.excitation_sweep(frequency_list=winding1_frequency_list, current_list_list=amplitude_list,
                                     phi_deg_list_list=phase_rad_list)

            elif self.md_simulation_type_comboBox.currentText() == self.translation_dict['transformer']:
                amplitude1_list = []
                for amplitude1_value, amplitude2_value in zip(winding1_amplitude_list, winding2_amplitude_list):
                    amplitude1_list.append([amplitude1_value, amplitude2_value])

                phase1_rad_list = []
                for phase1_value, phase2_value in zip(winding1_phi_rad_list, winding2_phi_rad_list):
                    phase1_rad_list.append([phase1_value, phase2_value])

                geo.excitation_sweep(frequency_list=winding1_frequency_list,
                                     current_list_list=amplitude1_list,
                                     phi_deg_list_list=phase1_rad_list)

        geo.thermal_simulation(thermal_conductivity_dict, boundary_temperatures, boundary_flags, case_gap_top,
                               case_gap_right, case_gap_bot, True,
                               color_scheme=fmt.colors_ba_jonas, colors_geometry=fmt.colors_geometry_ba_jonas, flag_insulation=self.flag_insulation)

        # Because the isolations inside of the winding window are not implemented in femm simulation.
        # The validation only works when the isolations for the FEMMT thermal simulation are turned off.
        # geo.femm_thermal_validation(thermal_conductivity_dict, femm_boundary_temperature, case_gap_top, case_gap_right, case_gap_bot)


def clear_layout(layout):
    """Clear the layout."""
    while layout.count():
        child = layout.takeAt(0)
        if isinstance(child, QtWidgets.QSpacerItem):
            layout.removeItem(child)  # for spacer item
        elif child.widget() or child:  # child for comboBoxType and child.widget() for custom class types
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