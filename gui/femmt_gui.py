import sys
import pandas as pd
import gmsh
import matplotlib.pyplot as plt
import numpy as np
from PyQt5.QtWidgets import QApplication, QMainWindow
from PyQt5 import QtCore, uic, QtGui, QtWidgets
from PyQt5.QtGui import QIcon, QPixmap, QDoubleValidator, QValidator, QIntValidator
import femmt as fmt
import json
from typing import List, Union, Optional
import os
import PIL

# import sys
# import matplotlib
#
# matplotlib.use('Qt5Agg')
#
# from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
# from matplotlib.figure import Figure

float_validator = QDoubleValidator()
int_validator = QIntValidator()


def comma_str_to_point_float(input_str: str) -> float:
    """
    Workaround to convert a comma (depends on the user input) in a point (needed for python)

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


class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.md_simulation_type_comboBox = None
        uic.loadUi('femmt_gui.ui', self)
        _translate = QtCore.QCoreApplication.translate
        #self.setWindowIcon(QIcon('Images\\logo.png'))
        self.setWindowTitle(_translate("MainWindow", "FEM Magnetics Toolbox"))
        pixmap = QPixmap('ferriteCore.png')
        #self.coreImageLabel.setPixmap(pixmap)
        #self.imageBoxImageLabel.setPixmap(pixmap)
        self.translation_dict =  {
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
            "+-5": "+/- 5%",
            "excel": "MS Excel"
        }

        "******* Manual Design *********"

        "Signals in Definition Tab"
        # simulation
        self.md_simulation_type_comboBox.currentTextChanged.connect(self.md_change_simulation_type)
        # core
        self.md_core_geometry_comboBox.currentTextChanged.connect(self.md_set_core_geometry_from_database)
        # windings
        self.md_winding1_type_comboBox.currentTextChanged.connect(self.md_winding1_change_wire_type)
        self.md_winding2_type_comboBox.currentTextChanged.connect(self.md_winding2_change_wire_type)
        self.md_winding1_litz_material_comboBox.currentTextChanged.connect(self.md_winding1_set_litz_parameters_from_litz_database)
        self.md_winding2_litz_material_comboBox.currentTextChanged.connect(self.md_winding2_set_litz_parameters_from_litz_database)
        self.md_winding1_implicit_litz_comboBox.currentTextChanged.connect(self.md_winding1_change_litz_implicit)
        self.md_winding2_implicit_litz_comboBox.currentTextChanged.connect(self.md_winding2_change_litz_implicit)

        # air gaps
        self.md_air_gap_count_comboBox.currentTextChanged.connect(self.md_change_air_gap_count)
        self.md_air_gap_placement_method_comboBox.currentTextChanged.connect(self.md_change_air_gap_placement)

        # visualization
        self.md_gmsh_visualisation_QPushButton.clicked.connect(self.md_gmsh_pre_visualisation)

        # Results
        self.Inductance_value_pushButton.clicked.connect(self.inductancecalc)


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
        self.md_isolation_s2s_lineEdit.setValidator(float_validator)
        self.md_isolation_p2s_lineEdit.setValidator(float_validator)
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
        self.md_core_geometry_comboBox.setToolTip("Chose a core geometry from the database. Chose 'Manual' to insert any parameters.")
        self.md_winding1_litz_material_comboBox.setToolTip("Chose a litz from the database. Chose 'Manual' to insert any parameters")
        self.md_winding1_implicit_litz_comboBox.setToolTip(
            "To describe a strand, 3 arguments are sufficient. Select here which of the arguments should not be entered.")
        self.md_winding2_implicit_litz_comboBox.setToolTip("To describe a strand, 3 arguments are sufficient. Select here which of the arguments should not be entered.")
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
        self.aut_simulation_type_comboBox.currentTextChanged.connect(self.aut_change_simulation_type)

    def aut_initialize_controls(self) -> None:
        """
        Initialize the comboboxes with pre-defined values.

        :return: None
        :rtype: None
        """
        aut_simulation_type_options = [self.translation_dict['inductor'], self.translation_dict['transformer']]
        aut_core_material_options = ['N95']
        aut_winding_material_options = [key for key in fmt.wire_material_database()]
        aut_winding_type_options = [self.translation_dict['litz'], self.translation_dict['solid']]
        aut_implicit_litz_options = [self.translation_dict["implicit_litz_radius"], self.translation_dict["implicit_ff"],
                                    self.translation_dict['implicit_strands_number']]
        aut_air_gap_method_options = [self.translation_dict["percent"]]
        aut_winding_scheme_options = [self.translation_dict["square"], self.translation_dict["hexa"]]
        aut_winding_litz_material_options = [key for key in fmt.litz_database()]
        aut_winding_litz_material_options.insert(0, 'Manual')
        aut_tolerance_val_options = [self.translation_dict['+-5']]

        for option in aut_simulation_type_options:
            self.aut_simulation_type_comboBox.addItem(option)
        for option in aut_core_material_options:
            self.aut_core_material_comboBox.addItem(option)
        for option in aut_winding_material_options:
            self.aut_winding1_material_comboBox.addItem(option)
            self.aut_winding2_material_comboBox.addItem(option)
        for option in aut_winding_type_options:
            self.aut_winding1_type_comboBox.addItem(option)
            self.aut_winding2_type_comboBox.addItem(option)
        for option in aut_implicit_litz_options:
            self.aut_winding1_implicit_litz_comboBox.addItem(option)
            self.aut_winding2_implicit_litz_comboBox.addItem(option)
        for option in aut_air_gap_method_options:
            self.aut_air_gap_placement_method_comboBox.addItem(option)
        for option in aut_winding_scheme_options:
            self.aut_winding1_scheme_comboBox.addItem(option)
            self.aut_winding2_scheme_comboBox.addItem(option)
        for option in aut_winding_litz_material_options:
            self.aut_winding1_litz_material_comboBox.addItem(option)
            self.aut_winding2_litz_material_comboBox.addItem(option)
        for option in aut_tolerance_val_options:
            self.aut_tolerance_val_comboBox.addItem(option)

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
        self.aut_step_winding1_turns_lineEdit.setPlaceholderText("Step value")
        self.aut_min_winding2_turns_lineEdit.setPlaceholderText("Minimum value")
        self.aut_max_winding2_turns_lineEdit.setPlaceholderText("Maximum value")
        self.aut_step_winding2_turns_lineEdit.setPlaceholderText("Step value")
        self.aut_min_air_gap_count_lineEdit.setPlaceholderText("Minimum value")
        self.aut_max_air_gap_count_lineEdit.setPlaceholderText("Maximum value")
        self.aut_goal_inductance_val_lineEdit.setPlaceholderText("Goal inductance")

        "Signals in FEM Simulations Tab"

        aut_download_options = [self.translation_dict['excel']]

        for option in aut_download_options:
            self.aut_download_comboBox.addItem(option)

        self.aut_pos_mod_sim_pushButton.clicked.connect(self.aut_reluctance_model_matrix)

        self.aut_pos_mod_download_pushButton.clicked.connect(self.aut_download_pos_model_data)

    def aut_reluctance_model_matrix(self):
        self.aut_rel_model_calc_timeline.setText("Simulating possible models..")
        matrix = [[0.0149, 0.0295, 0.01105, 0.0001],
                  [0.0149, 0.0295, 0.01105, 0.0002],
                  [0.0149, 0.0295, 0.01105, 0.0003],
                  [0.0149, 0.0295, 0.01105, 0.0004]]
        num_r = 4
        flag = 0
        m = []
        for i in range(num_r):
            m.append([float(x) for x in matrix[i]])
        for i in m:
            self.aut_action_run_simulation(i)
            flag += 1
        if flag == 4:
            self.aut_rel_model_calc_timeline.setText("Simulation ended.Ready to download!")


    def aut_action_run_simulation(self, sim_value):

        geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor)
        core = fmt.Core(core_w=sim_value[0], window_h=sim_value[1], window_w=sim_value[2],
                        mu_rel=3100, phi_mu_deg=12,
                        sigma=0.6)
        geo.set_core(core)

        # 3. set air gap parameters
        air_gaps = fmt.AirGaps(fmt.AirGapMethod.Manually, core)
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.001, sim_value[3])
        geo.set_air_gaps(air_gaps)

        # 4. set conductor parameters: use solid wires
        winding = fmt.Winding(8, 0, fmt.Conductivity.Copper, fmt.WindingType.Primary, fmt.WindingScheme.Square)
        winding.set_litz_conductor(None, 600, 35.5e-6, 0.6)
        # winding.set_solid_conductor(0.0015)
        geo.set_windings([winding])

        # 5. set isolations
        isolation = fmt.Isolation()
        isolation.add_core_isolations(0.001, 0.001, 0.002, 0.001)
        isolation.add_winding_isolations(0.0001)
        geo.set_isolation(isolation)

        # 5. create the model
        geo.create_model(freq=100000, visualize_before=False, save_png=False)

        # 6.a. start simulation
        geo.single_simulation(freq=100000, current=[3], show_results=True)

    def aut_download_pos_model_data(self):

        list1 = [0.0149, 0.0149, 0.0149]
        list2 = [0.0295, 0.0295, 0.0295]
        list3 = [0.01105, 0.01105, 0.01105]
        list4 = [0.0001, 0.0002, 0.0003]
        col1 = "core_w"
        col2 = "window_w"
        col3 = "window_h"
        col4 = "core_h"
        data = pd.DataFrame({col1: list1, col2: list2, col3: list3, col4: list4})
        data.to_excel('sample_data.xlsx', sheet_name='sheet1', index=False)
        self.aut_pos_model_download_status.setText("Downloaded!")



    def aut_change_simulation_type(self, simulation_type_from_combo_box: str) -> None:
        """
        Action performed when signal of aut_simulation_type_comboBox text has changed.
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
        self.aut_winding2_turns_lineEdit.setEnabled(status)
        self.aut_min_winding2_turns_lineEdit.setEnabled(status)
        self.aut_max_winding2_turns_lineEdit.setEnabled(status)
        self.aut_step_winding2_turns_lineEdit.setEnabled(status)
        self.aut_winding2_strands_lineEdit.setEnabled(status)
        self.aut_winding2_radius_lineEdit.setEnabled(False)
        self.aut_winding2_fill_factor_lineEdit.setEnabled(status)
        self.aut_winding2_strand_radius_lineEdit.setEnabled(status)
        self.aut_winding2_litz_material_comboBox.setEnabled(status)
        self.aut_winding2_groupBox.setVisible(status)

        # Set turns of winding 2 (enable and visible)
        self.aut_winding2_turns_lineEdit.setVisible(status)
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



    def md_initialize_controls(self) -> None:
        """
        Initialize the comboboxes with pre-defined values.

        :return: None
        :rtype: None
        """
        md_simulation_type_options = [self.translation_dict['inductor'], self.translation_dict['transformer']]
        md_core_material_options = ['N95']
        md_winding_material_options = [key for key in fmt.wire_material_database()]
        md_winding_type_options = [self.translation_dict['litz'], self.translation_dict['solid']]
        md_implicit_litz_options = [self.translation_dict["implicit_litz_radius"], self.translation_dict["implicit_ff"], self.translation_dict['implicit_strands_number']]
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
        Action performed when signal of md_simulation_type_comboBox text has changed.
        Action will be enabling / disabling user inputs for not-used windings.

        :param simulation_type_from_combo_box:
        :type simulation_type_from_combo_box: str
        :return: None
        :rtype: None
        """
        if simulation_type_from_combo_box == self.translation_dict['inductor']:
            self.md_winding2_enable(False)

        elif simulation_type_from_combo_box == self.translation_dict['transformer']:
            # set winding definitions of winding 2 to editable
            self.md_winding2_enable(True)

        elif simulation_type_from_combo_box == self.translation_dict['integrated transformer']:
            # set winding definitions of winding 2 to editable
            self.md_winding2_enable(True)



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
        self.md_isolation_s2s_label.setVisible(status)
        self.md_isolation_p2s_label.setVisible(status)

    def md_gmsh_pre_visualisation(self):
        geo = self.md_setup_geometry()
        #geo.create_model(freq=100000, visualize_before=False, do_meshing=False, save_png=True)
        geo.create_model(freq=comma_str_to_point_float(self.md_base_frequency_lineEdit.text()), visualize_before=False, save_png=True)
        image_pre_visualisation = PIL.Image.open(geo.hybrid_color_visualize_file)

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
        im_crop.save(geo.hybrid_color_visualize_file, quality=95)

        pixmap = QPixmap(geo.hybrid_color_visualize_file)
        self.md_gmsh_visualisation_QLabel.setPixmap(pixmap)
        self.md_gmsh_visualisation_QLabel.setMask(pixmap.mask())
        self.md_gmsh_visualisation_QLabel.show()

    def md_set_core_geometry_from_database(self):
        core_dict = fmt.core_database()
        core_type = self.md_core_geometry_comboBox.currentText()

        if core_type != 'Manual':
            core = core_dict[core_type]

            self.md_core_width_lineEdit.setText(str(core["core_w"]))
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
        litz_dict = fmt.litz_database()
        litz_type = self.md_winding1_litz_material_comboBox.currentText()

        if litz_type != 'Manual':
            litz = litz_dict[litz_type]

            self.md_winding1_strands_lineEdit.setText(str(litz["strands_numbers"]))
            self.md_winding1_strand_radius_lineEdit.setText(str(litz["strand_radii"]))
            self.md_winding1_radius_lineEdit.setText(str(litz["conductor_radii"]))
            self.md_winding1_fill_factor_lineEdit.setText(str(litz["ff"]))

            for key, value in enumerate(["implicit_litz_radius", "implicit_ff", 'implicit_strands_number']):
                if value == litz["implicit"]:
                    self.md_winding1_implicit_litz_comboBox.setCurrentIndex(key)

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
        Enables / Disables input parameter fields for different "implicit xyz" types in case of litz wire:
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

    def  md_winding1_change_wire_type(self, wire_type_from_combot_box: str) -> None:
        """
        Enables / Disables input parameter for litz/solid wire
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
        Enables / Disables input parameter fields for different "implicit xyz" types in case of litz wire:
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
        Enables / Disables input parameter for litz/solid wire
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
            self.md_winding2_radius_lineEdit.setEnabled(True)
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
        Sets the number of editable air gap fields in dependence of the air gap count combobox.
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
        Changes the labels in case of different air gap placement methods
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
        Enables / Disables the input fields for air gap No. 1

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
        Enables / Disables the input fields for air gap No. 2

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
        Enables / Disables the input fields for air gap No. 3

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
        Enables / Disables the input fields for air gap No. 4

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
        Enables / Disables the input fields for air gap No. 5

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
        Enables / Disables the input fields for dc current.

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
        Enables / Disables the input fields for frequency/phase No. 1

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
        Enables / Disables the input fields for frequency/phase No. 2

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
        Enables / Disables the input fields for frequency/phase No. 3

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
        Enables / Disables the input fields for frequency/phase No. 4

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
        Enables / Disables the input fields for frequency/phase No. 5

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
        Enables / Disables the input fields for frequency/phase No. 6

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
        Enables / Disables the input fields for frequency/phase No. 7

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
        Enables / Disables the input fields for frequency/phase No. 8

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
        Changes the frequency field in case of checking/unchecking the frequency-checkboxes

        :param status: 0 for disabling, anything else for enabling freqency boxes
        :type status: int
        :return: None
        :rtype: None
        """
        self.md_dc_enable(False) if status == 0 else self.md_dc_enable(True)

    def md_change_frequencies_1(self, status: int) -> None:
        """
        Changes the frequency field in case of checking/unchecking the frequency-checkboxes

        :param status: 0 for disabling, anything else for enabling freqency boxes
        :type status: int
        :return: None
        :rtype: None
        """
        self.md_f1_enable(False) if status == 0 else self.md_f1_enable(True)

    def md_change_frequencies_2(self, status) -> None:
        """
        Changes the frequency field in case of checking/unchecking the frequency-checkboxes

        :param status: 0 for disabling, anything else for enabling freqency boxes
        :type status: int
        :return: None
        :rtype: None
        """
        self.md_f2_enable(False) if status == 0 else self.md_f2_enable(True)

    def md_change_frequencies_3(self, status) -> None:
        """
        Changes the frequency field in case of checking/unchecking the frequency-checkboxes

        :param status: 0 for disabling, anything else for enabling freqency boxes
        :type status: int
        :return: None
        :rtype: None
        """
        self.md_f3_enable(False) if status == 0 else self.md_f3_enable(True)

    def md_change_frequencies_4(self, status) -> None:
        """
        Changes the frequency field in case of checking/unchecking the frequency-checkboxes

        :param status: 0 for disabling, anything else for enabling freqency boxes
        :type status: int
        :return: None
        :rtype: None
        """
        self.md_f4_enable(False) if status == 0 else self.md_f4_enable(True)

    def md_change_frequencies_5(self, status) -> None:
        """
        Changes the frequency field in case of checking/unchecking the frequency-checkboxes

        :param status: 0 for disabling, anything else for enabling freqency boxes
        :type status: int
        :return: None
        :rtype: None
        """
        self.md_f5_enable(False) if status == 0 else self.md_f5_enable(True)

    def md_change_frequencies_6(self, status) -> None:
        """
        Changes the frequency field in case of checking/unchecking the frequency-checkboxes

        :param status: 0 for disabling, anything else for enabling freqency boxes
        :type status: int
        :return: None
        :rtype: None
        """
        self.md_f6_enable(False) if status == 0 else self.md_f6_enable(True)

    def md_change_frequencies_7(self, status) -> None:
        """
        Changes the frequency field in case of checking/unchecking the frequency-checkboxes

        :param status: 0 for disabling, anything else for enabling freqency boxes
        :type status: int
        :return: None
        :rtype: None
        """
        self.md_f7_enable(False) if status == 0 else self.md_f7_enable(True)

    def md_change_frequencies_8(self, status) -> None:
        """
        Changes the frequency field in case of checking/unchecking the frequency-checkboxes

        :param status: 0 for disabling, anything else for enabling freqency boxes
        :type status: int
        :return: None
        :rtype: None
        """
        self.md_f8_enable(False) if status == 0 else self.md_f8_enable(True)

    def md_redraw_input_signals(self) -> None:
        """
        Generate visual graphics for the input signals
        Generates a graphic. This graphic is read and insertet to the gui.

        :return: None
        :rtype: None
        """

        winding1_frequency_list, winding1_amplitude_list, winding1_phi_rad_list, winding2_frequency_list, winding2_amplitude_list, winding2_phi_rad_list = self.md_get_frequency_lists()


        fmt.plot_fourier_coefficients(winding1_frequency_list, winding1_amplitude_list, winding1_phi_rad_list, figure_directory =  "./md_winding_1.png")
        pixmap = QPixmap("./md_winding_1.png")
        self.md_graphic_winding_1.setPixmap(pixmap)
        self.md_graphic_winding_1.setMask(pixmap.mask())
        self.md_graphic_winding_1.show()
        if self.md_simulation_type_comboBox.currentText() != self.translation_dict['inductor']:
            fmt.plot_fourier_coefficients(winding2_frequency_list, winding2_amplitude_list, winding2_phi_rad_list, figure_directory = "./md_winding_2.png")
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

        :return: [winding1_frequency_list, winding1_amplitude_list, winding1_phi_rad_list, winding2_frequency_list, winding2_amplitude_list, winding2_phi_rad_list]
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
        if self.md_winding1_scheme_comboBox.currentText() == self.translation_dict["square"]:
            scheme = 'square'
        elif self.md_winding1_scheme_comboBox.currentText() == self.translation_dict["hexa"]:
            scheme = 'hexa'

    def md_setup_geometry(self):
        """
        Sets up the core and conductor geometry depending on the GUI input parameters

        returns: femmt MagneticComponent

        """
        if self.md_simulation_type_comboBox.currentText() == self.translation_dict['inductor']:
            self.md_simulation_QLabel.setText('simulation startet...')

            #geo = fmt.MagneticComponent(component_type="inductor")
            geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor)

            # -----------------------------------------------
            # Core
            # -----------------------------------------------
            """
            geo.core.update(type="EI",
                            core_w=comma_str_to_point_float(self.md_core_width_lineEdit.text()),
                            window_h=comma_str_to_point_float(self.md_window_height_lineEdit.text()),
                            window_w=comma_str_to_point_float(self.md_window_width_lineEdit.text()))"""

            core = fmt.Core(core_w=comma_str_to_point_float(self.md_core_width_lineEdit.text()),
                            window_w=comma_str_to_point_float(self.md_window_width_lineEdit.text()),
                            window_h=comma_str_to_point_float(self.md_window_height_lineEdit.text()),
                            mu_rel=3100, phi_mu_deg=12,
                            sigma=0.6)
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

            if air_gap_count >=4:
                md_air_gap_4_height = comma_str_to_point_float(self.md_air_gap_4_length_lineEdit.text())
                md_air_gap_4_position = comma_str_to_point_float(self.md_air_gap_4_position_lineEdit.text())

                air_gap_heigth_array.append(md_air_gap_4_height)
                air_gap_position_array.append(md_air_gap_4_position)
                air_gap_position_tag_array.append(0)

            if air_gap_count >=5:
                md_air_gap_5_height = comma_str_to_point_float(self.md_air_gap_5_length_lineEdit.text())
                md_air_gap_5_position = comma_str_to_point_float(self.md_air_gap_5_position_lineEdit.text())

                air_gap_heigth_array.append(md_air_gap_5_height)
                air_gap_position_array.append(md_air_gap_5_position)
                air_gap_position_tag_array.append(0)


            if air_gap_count == 0:
                air_gaps = fmt.AirGaps(fmt.AirGapMethod.Center, core)
                #air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, None, 0.0005)
                geo.set_air_gaps(air_gaps)
                """
                geo.air_gaps.update(method="percent",
                                    n_air_gaps=0,
                                    air_gap_h=[],
                                    position_tag=[],
                                    air_gap_position=[])"""

            elif self.md_air_gap_placement_method_comboBox.currentText() == self.translation_dict["percent"] and air_gap_count >= 1:

                air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
                for i in range(1, air_gap_count+1):
                    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, air_gap_position_array[i-1], air_gap_heigth_array[i-1])
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
                for i in range(1, air_gap_count+1):
                    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, air_gap_position_array[i-1], air_gap_heigth_array[i-1])
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

            # -----------------------------------------------
            # Conductors
            # -----------------------------------------------


            if self.md_winding1_scheme_comboBox.currentText() == self.translation_dict["square"]:
                fmt.WindingScheme = fmt.WindingScheme.Square
            elif self.md_winding1_scheme_comboBox.currentText() == self.translation_dict["hexa"]:
                fmt.WindingScheme = fmt.WindingScheme.Hexagonal

            if self.md_winding1_type_comboBox.currentText() == self.translation_dict['solid']:
                self.md_simulation_QLabel.setText('setze conductors')
                winding = fmt.Winding(int(self.md_winding1_turns_lineEdit.text()), 0,
                                      fmt.Conductivity.Copper,
                                      fmt.WindingType.Primary,
                                      fmt.WindingScheme)
                cond_radii = comma_str_to_point_float(self.md_winding1_radius_lineEdit.text())
                winding.set_solid_conductor(cond_radii)
                geo.set_windings([winding])
                isolation = fmt.Isolation()
                isolation.add_core_isolations(comma_str_to_point_float(self.md_isolation_core2cond_top_lineEdit.text()),
                                              comma_str_to_point_float(self.md_isolation_core2cond_bot_lineEdit.text()),
                                              comma_str_to_point_float(self.md_isolation_core2cond_inner_lineEdit.text()),
                                              comma_str_to_point_float(self.md_isolation_core2cond_outer_lineEdit.text()))
                isolation.add_winding_isolations(comma_str_to_point_float(self.md_isolation_p2p_lineEdit.text()))
                geo.set_isolation(isolation)

                """
                geo.update_conductors(n_turns=[[int(self.md_winding1_turns_lineEdit.text())]],
                                      conductor_type=['solid'],
                                      conductor_radii=[comma_str_to_point_float(self.md_winding1_radius_lineEdit.text())],
                                      winding=["primary"],
                                      scheme=[scheme],
                                      core_cond_isolation=[comma_str_to_point_float(self.md_isolation_core2cond_top_lineEdit.text()),
                                                           comma_str_to_point_float(self.md_isolation_core2cond_bot_lineEdit.text()),
                                                           comma_str_to_point_float(self.md_isolation_core2cond_inner_lineEdit.text()),
                                                           comma_str_to_point_float(self.md_isolation_core2cond_outer_lineEdit.text())],
                                      cond_cond_isolation=[comma_str_to_point_float(self.self.md_isolation_p2p_lineEdit.text())],
                                      conductivity_sigma=[self.md_winding1_material_comboBox.currentText()])"""

            elif self.md_winding1_type_comboBox.currentText() == self.translation_dict['litz']:
                litz_para_type = ''
                if self.md_winding1_implicit_litz_comboBox.currentText() == self.translation_dict['implicit_litz_radius']:
                    litz_para_type = "implicit_litz_radius"
                elif self.md_winding1_implicit_litz_comboBox.currentText() == self.translation_dict[
                    'implicit_ff']:
                    litz_para_type = 'implicit_ff'
                elif self.md_winding1_implicit_litz_comboBox.currentText() == self.translation_dict[
                    'implicit_strands_number']:
                    litz_para_type = 'implicit_strands_number'
                winding = fmt.Winding(int(self.md_winding1_turns_lineEdit.text()), 0,
                                      fmt.Conductivity.Copper,
                                      fmt.WindingType.Primary,
                                      fmt.WindingScheme)
                cond_radii = comma_str_to_point_float(self.md_winding1_radius_lineEdit.text())
                winding.set_litz_conductor(None,
                                           comma_str_to_point_float(self.md_winding1_strands_lineEdit.text()),
                                           comma_str_to_point_float(self.md_winding1_strand_radius_lineEdit.text()),
                                           comma_str_to_point_float(self.md_winding1_fill_factor_lineEdit.text()))
                geo.set_windings([winding])
                isolation = fmt.Isolation()
                isolation.add_core_isolations(comma_str_to_point_float(self.md_isolation_core2cond_top_lineEdit.text()),
                                              comma_str_to_point_float(self.md_isolation_core2cond_bot_lineEdit.text()),
                                              comma_str_to_point_float(self.md_isolation_core2cond_inner_lineEdit.text()),
                                              comma_str_to_point_float(self.md_isolation_core2cond_outer_lineEdit.text()))
                isolation.add_winding_isolations(comma_str_to_point_float(self.md_isolation_p2p_lineEdit.text()))
                geo.set_isolation(isolation)


                """
                geo.update_conductors(n_turns=[[int(self.md_winding1_turns_lineEdit.text())]],
                                      conductor_type=['litz'],
                                      conductor_radii=[comma_str_to_point_float(self.md_winding1_radius_lineEdit.text())],
                                      winding=["primary"],
                                      scheme=[scheme],
                                      litz_para_type=[litz_para_type],
                                      ff=[comma_str_to_point_float(self.md_winding1_fill_factor_lineEdit.text())],
                                      strands_numbers=[comma_str_to_point_float(self.md_winding1_strands_lineEdit.text())],
                                      strand_radii=[comma_str_to_point_float(self.md_winding1_strand_radius_lineEdit.text())],
                                      core_cond_isolation=[comma_str_to_point_float(self.md_isolation_core2cond_top_lineEdit.text()),
                                                           comma_str_to_point_float(self.md_isolation_core2cond_bot_lineEdit.text()),
                                                           comma_str_to_point_float(self.md_isolation_core2cond_inner_lineEdit.text()),
                                                           comma_str_to_point_float(self.md_isolation_core2cond_outer_lineEdit.text())],
                                      cond_cond_isolation=[comma_str_to_point_float(self.md_isolation_p2p_lineEdit.text())],
                                      conductivity_sigma=[self.md_winding1_material_comboBox.currentText()])"""

        elif self.md_simulation_type_comboBox.currentText() == 'transformer':


            self.md_simulation_QLabel.setText('simulation startet...')

            geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Transformer)
            #geo = fmt.MagneticComponent(component_type="transformer")


            # -----------------------------------------------
            # Core
            # -----------------------------------------------
            core = fmt.Core(core_w=comma_str_to_point_float(self.md_core_width_lineEdit.text()),
                            window_w=comma_str_to_point_float(self.md_window_width_lineEdit.text()),
                            window_h=comma_str_to_point_float(self.md_window_height_lineEdit.text()),
                            mu_rel=3100, phi_mu_deg=12,
                            sigma=0.6)
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

            if air_gap_count >=4:
                md_air_gap_4_height = comma_str_to_point_float(self.md_air_gap_4_length_lineEdit.text())
                md_air_gap_4_position = comma_str_to_point_float(self.md_air_gap_4_position_lineEdit.text())

                air_gap_heigth_array.append(md_air_gap_4_height)
                air_gap_position_array.append(md_air_gap_4_position)
                air_gap_position_tag_array.append(0)

            if air_gap_count >=5:
                md_air_gap_5_height = comma_str_to_point_float(self.md_air_gap_5_length_lineEdit.text())
                md_air_gap_5_position = comma_str_to_point_float(self.md_air_gap_5_position_lineEdit.text())

                air_gap_heigth_array.append(md_air_gap_5_height)
                air_gap_position_array.append(md_air_gap_5_position)
                air_gap_position_tag_array.append(0)


            if air_gap_count == 0:
                air_gaps = fmt.AirGaps(fmt.AirGapMethod.Center, core)
                #air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, None, 0.0005)
                geo.set_air_gaps(air_gaps)
                """
                geo.air_gaps.update(method="percent",
                                    n_air_gaps=0,
                                    air_gap_h=[],
                                    position_tag=[],
                                    air_gap_position=[])"""

            elif self.md_air_gap_placement_method_comboBox.currentText() == self.translation_dict["percent"] and air_gap_count >= 1:

                air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
                for i in range(1, air_gap_count+1):
                    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, air_gap_position_array[i-1], air_gap_heigth_array[i-1])
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
                for i in range(1, air_gap_count+1):
                    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, air_gap_position_array[i-1], air_gap_heigth_array[i-1])
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

            # -----------------------------------------------
            # Conductors
            # -----------------------------------------------

            if self.md_winding1_scheme_comboBox.currentText() == self.translation_dict["square"]:
                fmt.WindingScheme = fmt.WindingScheme.Square
            elif self.md_winding1_scheme_comboBox.currentText() == self.translation_dict["hexa"]:
                fmt.WindingScheme = fmt.WindingScheme.Hexagonal
            if self.md_winding2_scheme_comboBox.currentText() == self.translation_dict["square"]:
                fmt.WindingScheme = fmt.WindingScheme.Square
            elif self.md_winding2_scheme_comboBox.currentText() == self.translation_dict["hexa"]:
                fmt.WindingScheme = fmt.WindingScheme.Hexagonal

            wdg1type = self.md_winding1_type_comboBox.currentText()
            wdg2type = self.md_winding2_type_comboBox.currentText()

            if wdg1type == self.translation_dict['solid']:
                if wdg2type == self.translation_dict['solid']:
                    self.md_simulation_QLabel.setText('setze conductors')

                    winding1 = fmt.Winding(int(self.md_winding1_turns_lineEdit.text()), 0, fmt.Conductivity.Copper, fmt.WindingType.Primary,
                                           fmt.WindingScheme)
                    winding1.set_solid_conductor(comma_str_to_point_float(self.md_winding1_radius_lineEdit.text()))

                    winding2 = fmt.Winding(0, int(self.md_winding2_turns_lineEdit.text()), fmt.Conductivity.Copper, fmt.WindingType.Secondary,
                                           fmt.WindingScheme)
                    winding2.set_solid_conductor(comma_str_to_point_float(self.md_winding2_radius_lineEdit.text()))

                    geo.set_windings([winding1, winding2])

                    isolation = fmt.Isolation()
                    isolation.add_core_isolations(comma_str_to_point_float(self.md_isolation_core2cond_top_lineEdit.text()),
                                                  comma_str_to_point_float(self.md_isolation_core2cond_bot_lineEdit.text()),
                                                  comma_str_to_point_float(self.md_isolation_core2cond_inner_lineEdit.text()),
                                                  comma_str_to_point_float(self.md_isolation_core2cond_outer_lineEdit.text()))
                    isolation.add_winding_isolations(comma_str_to_point_float(self.md_isolation_p2p_lineEdit.text()),
                                                     comma_str_to_point_float(self.md_isolation_s2s_lineEdit.text()),
                                                     comma_str_to_point_float(self.md_isolation_p2s_lineEdit.text()))
                    geo.set_isolation(isolation)
                    """
                    geo.update_conductors(n_turns=[[int(self.md_winding1_turns_lineEdit.text()), 0],
                                                   [0, int(self.md_winding2_turns_lineEdit.text())]],
                                      conductor_type=["solid", "solid"],
                                      litz_para_type=['implicit_litz_radius', 'implicit_litz_radius'],
                                      ff=[comma_str_to_point_float(self.md_winding1_fill_factor_lineEdit.text()),
                                          comma_str_to_point_float(self.md_winding2_fill_factor_lineEdit.text())],
                                      strands_numbers=[comma_str_to_point_float(self.md_winding1_strands_lineEdit.text()),
                                                       comma_str_to_point_float(self.md_winding2_strands_lineEdit.text())],
                                      strand_radii=[comma_str_to_point_float(self.md_winding1_strand_radius_lineEdit.text()),
                                                    comma_str_to_point_float(self.md_winding2_strand_radius_lineEdit.text())],
                                      conductor_radii=[comma_str_to_point_float(self.md_winding1_radius_lineEdit.text()),
                                                       comma_str_to_point_float(self.md_winding2_radius_lineEdit.text())],
                                      winding=["primary", "secondary"],
                                      scheme=[scheme1, scheme2],
                                      core_cond_isolation=[comma_str_to_point_float(self.md_isolation_core2cond_top_lineEdit.text()),
                                          comma_str_to_point_float(self.md_isolation_core2cond_bot_lineEdit.text()),
                                          comma_str_to_point_float(self.md_isolation_core2cond_inner_lineEdit.text()),
                                          comma_str_to_point_float(self.md_isolation_core2cond_outer_lineEdit.text())],
                                      cond_cond_isolation=[comma_str_to_point_float(self.md_isolation_p2p_lineEdit.text()),
                                                           comma_str_to_point_float(self.md_isolation_s2s_lineEdit.text()),
                                                           comma_str_to_point_float(self.md_isolation_p2s_lineEdit.text())],
                                      conductivity_sigma=[self.md_winding1_material_comboBox.currentText(),
                                                          self.md_winding2_material_comboBox.currentText()])"""

                elif wdg2type == self.translation_dict['litz']:
                    litz_para_type = ''
                    if self.md_winding2_implicit_litz_comboBox.currentText() == self.translation_dict[
                        'implicit_litz_radius']:
                        litz_para_type = "implicit_litz_radius"
                    elif self.md_winding2_implicit_litz_comboBox.currentText() == self.translation_dict['implicit_ff']:
                        litz_para_type = 'implicit_ff'
                    elif self.md_winding2_implicit_litz_comboBox.currentText() == self.translation_dict[
                        'implicit_strands_number']:
                        litz_para_type = 'implicit_strands_number'

                    winding1 = fmt.Winding(int(self.md_winding1_turns_lineEdit.text()), 0, fmt.Conductivity.Copper, fmt.WindingType.Primary,
                                           fmt.WindingScheme)
                    winding1.set_solid_conductor(comma_str_to_point_float(self.md_winding1_radius_lineEdit.text()))

                    winding2 = fmt.Winding(0, int(self.md_winding2_turns_lineEdit.text()), fmt.Conductivity.Copper, fmt.WindingType.Secondary,
                                           fmt.WindingScheme)
                    winding2.set_litz_conductor(None,
                                                comma_str_to_point_float(self.md_winding2_strands_lineEdit.text()),
                                                comma_str_to_point_float(self.md_winding2_strand_radius_lineEdit.text()),
                                                comma_str_to_point_float(self.md_winding2_fill_factor_lineEdit.text()))

                    geo.set_windings([winding1, winding2])

                    isolation = fmt.Isolation()
                    isolation.add_core_isolations(comma_str_to_point_float(self.md_isolation_core2cond_top_lineEdit.text()),
                                                  comma_str_to_point_float(self.md_isolation_core2cond_bot_lineEdit.text()),
                                                  comma_str_to_point_float(self.md_isolation_core2cond_inner_lineEdit.text()),
                                                  comma_str_to_point_float(self.md_isolation_core2cond_outer_lineEdit.text()))
                    isolation.add_winding_isolations(comma_str_to_point_float(self.md_isolation_p2p_lineEdit.text()),
                                                     comma_str_to_point_float(self.md_isolation_s2s_lineEdit.text()),
                                                     comma_str_to_point_float(self.md_isolation_p2s_lineEdit.text()))
                    geo.set_isolation(isolation)
                    """
                    geo.update_conductors(n_turns=[[int(self.md_winding1_turns_lineEdit.text()), 0],
                                                   [0, int(self.md_winding2_turns_lineEdit.text())]],
                                          conductor_type=["solid", "litz"],
                                          litz_para_type=['implicit_litz_radius', litz_para_type],
                                          ff=[comma_str_to_point_float(self.md_winding1_fill_factor_lineEdit.text()),
                                              comma_str_to_point_float(self.md_winding2_fill_factor_lineEdit.text())],
                                          strands_numbers=[
                                              comma_str_to_point_float(self.md_winding1_strands_lineEdit.text()),
                                              comma_str_to_point_float(self.md_winding2_strands_lineEdit.text())],
                                          strand_radii=[
                                              comma_str_to_point_float(self.md_winding1_strand_radius_lineEdit.text()),
                                              comma_str_to_point_float(self.md_winding2_strand_radius_lineEdit.text())],
                                          conductor_radii=[
                                              comma_str_to_point_float(self.md_winding1_radius_lineEdit.text()),
                                              comma_str_to_point_float(self.md_winding2_radius_lineEdit.text())],
                                          winding=["primary", "secondary"],
                                          scheme=[scheme1, scheme2],
                                          core_cond_isolation=[
                                              comma_str_to_point_float(self.md_isolation_core2cond_top_lineEdit.text()),
                                              comma_str_to_point_float(self.md_isolation_core2cond_bot_lineEdit.text()),
                                              comma_str_to_point_float(
                                                  self.md_isolation_core2cond_inner_lineEdit.text()),
                                              comma_str_to_point_float(
                                                  self.md_isolation_core2cond_outer_lineEdit.text())],
                                          cond_cond_isolation=[
                                              comma_str_to_point_float(self.md_isolation_p2p_lineEdit.text()),
                                              comma_str_to_point_float(self.md_isolation_s2s_lineEdit.text()),
                                              comma_str_to_point_float(self.md_isolation_p2s_lineEdit.text())],
                                          conductivity_sigma=[self.md_winding1_material_comboBox.currentText(),
                                                              self.md_winding2_material_comboBox.currentText()])"""

            elif wdg1type == self.translation_dict['litz']:
                if wdg2type == self.translation_dict['litz']:
                    litz_para_type = ''
                    if self.md_winding1_implicit_litz_comboBox.currentText() == self.translation_dict['implicit_litz_radius']:
                        litz_para_type = "implicit_litz_radius"
                    elif self.md_winding1_implicit_litz_comboBox.currentText() == self.translation_dict['implicit_ff']:
                        litz_para_type = 'implicit_ff'
                    elif self.md_winding1_implicit_litz_comboBox.currentText() == self.translation_dict['implicit_strands_number']:
                        litz_para_type = 'implicit_strands_number'

                    if self.md_winding2_implicit_litz_comboBox.currentText() == self.translation_dict['implicit_litz_radius']:
                        litz_para_type = "implicit_litz_radius"
                    elif self.md_winding2_implicit_litz_comboBox.currentText() == self.translation_dict['implicit_ff']:
                        litz_para_type = 'implicit_ff'
                    elif self.md_winding2_implicit_litz_comboBox.currentText() == self.translation_dict['implicit_strands_number']:
                        litz_para_type = 'implicit_strands_number'

                    winding1 = fmt.Winding(int(self.md_winding1_turns_lineEdit.text()), 0, fmt.Conductivity.Copper, fmt.WindingType.Primary,
                                           fmt.WindingScheme)
                    winding1.set_litz_conductor(None,
                                                comma_str_to_point_float(self.md_winding1_strands_lineEdit.text()),
                                                comma_str_to_point_float(self.md_winding1_strand_radius_lineEdit.text()),
                                                comma_str_to_point_float(self.md_winding1_fill_factor_lineEdit.text()))

                    winding2 = fmt.Winding(0, int(self.md_winding2_turns_lineEdit.text()), fmt.Conductivity.Copper, fmt.WindingType.Secondary,
                                           fmt.WindingScheme)
                    winding2.set_litz_conductor(None,
                                                comma_str_to_point_float(self.md_winding2_strands_lineEdit.text()),
                                                comma_str_to_point_float(self.md_winding2_strand_radius_lineEdit.text()),
                                                comma_str_to_point_float(self.md_winding2_fill_factor_lineEdit.text()))

                    geo.set_windings([winding1, winding2])

                    isolation = fmt.Isolation()
                    isolation.add_core_isolations(comma_str_to_point_float(self.md_isolation_core2cond_top_lineEdit.text()),
                                                  comma_str_to_point_float(self.md_isolation_core2cond_bot_lineEdit.text()),
                                                  comma_str_to_point_float(self.md_isolation_core2cond_inner_lineEdit.text()),
                                                  comma_str_to_point_float(self.md_isolation_core2cond_outer_lineEdit.text()))
                    isolation.add_winding_isolations(comma_str_to_point_float(self.md_isolation_p2p_lineEdit.text()),
                                                     comma_str_to_point_float(self.md_isolation_s2s_lineEdit.text()),
                                                     comma_str_to_point_float(self.md_isolation_p2s_lineEdit.text()))
                    geo.set_isolation(isolation)
                    """
                    geo.update_conductors(n_turns=[[int(self.md_winding1_turns_lineEdit.text()), 0],
                                                   [0, int(self.md_winding2_turns_lineEdit.text())]],
                                          conductor_type=["litz", "litz"],
                                          litz_para_type=[litz_para_type, litz_para_type],
                                          ff=[comma_str_to_point_float(self.md_winding1_fill_factor_lineEdit.text()),
                                              comma_str_to_point_float(self.md_winding2_fill_factor_lineEdit.text())],
                                          strands_numbers=[
                                              comma_str_to_point_float(self.md_winding1_strands_lineEdit.text()),
                                              comma_str_to_point_float(self.md_winding2_strands_lineEdit.text())],
                                          strand_radii=[
                                              comma_str_to_point_float(self.md_winding1_strand_radius_lineEdit.text()),
                                              comma_str_to_point_float(self.md_winding2_strand_radius_lineEdit.text())],
                                          conductor_radii=[
                                              comma_str_to_point_float(self.md_winding1_radius_lineEdit.text()),
                                              comma_str_to_point_float(self.md_winding2_radius_lineEdit.text())],
                                          winding=["primary", "secondary"],
                                          scheme=[scheme1, scheme2],
                                          core_cond_isolation=[
                                              comma_str_to_point_float(self.md_isolation_core2cond_top_lineEdit.text()),
                                              comma_str_to_point_float(self.md_isolation_core2cond_bot_lineEdit.text()),
                                              comma_str_to_point_float(
                                                  self.md_isolation_core2cond_inner_lineEdit.text()),
                                              comma_str_to_point_float(
                                                  self.md_isolation_core2cond_outer_lineEdit.text())],
                                          cond_cond_isolation=[
                                              comma_str_to_point_float(self.md_isolation_p2p_lineEdit.text()),
                                              comma_str_to_point_float(self.md_isolation_s2s_lineEdit.text()),
                                              comma_str_to_point_float(self.md_isolation_p2s_lineEdit.text())],
                                          conductivity_sigma=[self.md_winding1_material_comboBox.currentText(),
                                                              self.md_winding2_material_comboBox.currentText()])"""

                elif wdg2type == self.translation_dict['solid']:
                    litz_para_type = ''
                    if self.md_winding1_implicit_litz_comboBox.currentText() == self.translation_dict['implicit_litz_radius']:
                        litz_para_type = "implicit_litz_radius"
                    elif self.md_winding1_implicit_litz_comboBox.currentText() == self.translation_dict['implicit_ff']:
                        litz_para_type = 'implicit_ff'
                    elif self.md_winding1_implicit_litz_comboBox.currentText() == self.translation_dict['implicit_strands_number']:
                        litz_para_type = 'implicit_strands_number'

                    winding1 = fmt.Winding(int(self.md_winding1_turns_lineEdit.text()), 0, fmt.Conductivity.Copper, fmt.WindingType.Primary,
                                           fmt.WindingScheme)
                    winding1.set_litz_conductor(None,
                                                comma_str_to_point_float(self.md_winding1_strands_lineEdit.text()),
                                                comma_str_to_point_float(self.md_winding1_strand_radius_lineEdit.text()),
                                                comma_str_to_point_float(self.md_winding1_fill_factor_lineEdit.text()))

                    winding2 = fmt.Winding(0, int(self.md_winding2_turns_lineEdit.text()), fmt.Conductivity.Copper, fmt.WindingType.Secondary,
                                           fmt.WindingScheme)
                    winding2.set_solid_conductor(comma_str_to_point_float(self.md_winding2_radius_lineEdit.text()))

                    geo.set_windings([winding1, winding2])

                    isolation = fmt.Isolation()
                    isolation.add_core_isolations(comma_str_to_point_float(self.md_isolation_core2cond_top_lineEdit.text()),
                                                  comma_str_to_point_float(self.md_isolation_core2cond_bot_lineEdit.text()),
                                                  comma_str_to_point_float(self.md_isolation_core2cond_inner_lineEdit.text()),
                                                  comma_str_to_point_float(self.md_isolation_core2cond_outer_lineEdit.text()))
                    isolation.add_winding_isolations(comma_str_to_point_float(self.md_isolation_p2p_lineEdit.text()),
                                                     comma_str_to_point_float(self.md_isolation_s2s_lineEdit.text()),
                                                     comma_str_to_point_float(self.md_isolation_p2s_lineEdit.text()))
                    geo.set_isolation(isolation)

                    """
                    geo.update_conductors(n_turns=[[int(self.md_winding1_turns_lineEdit.text()), 0],
                                               [0, int(self.md_winding2_turns_lineEdit.text())]],
                                      conductor_type=["litz", "solid"],
                                      litz_para_type=[litz_para_type, 'implicit_litz_radius'],
                                      ff=[comma_str_to_point_float(self.md_winding1_fill_factor_lineEdit.text()),
                                          comma_str_to_point_float(self.md_winding2_fill_factor_lineEdit.text())],
                                      strands_numbers=[
                                          comma_str_to_point_float(self.md_winding1_strands_lineEdit.text()),
                                          comma_str_to_point_float(self.md_winding2_strands_lineEdit.text())],
                                      strand_radii=[
                                          comma_str_to_point_float(self.md_winding1_strand_radius_lineEdit.text()),
                                          comma_str_to_point_float(self.md_winding2_strand_radius_lineEdit.text())],
                                      conductor_radii=[
                                          comma_str_to_point_float(self.md_winding1_radius_lineEdit.text()),
                                          comma_str_to_point_float(self.md_winding2_radius_lineEdit.text())],
                                      winding=["primary", "secondary"],
                                      scheme=[scheme1, scheme2],
                                      core_cond_isolation=[
                                          comma_str_to_point_float(self.md_isolation_core2cond_top_lineEdit.text()),
                                          comma_str_to_point_float(self.md_isolation_core2cond_bot_lineEdit.text()),
                                          comma_str_to_point_float(self.md_isolation_core2cond_inner_lineEdit.text()),
                                          comma_str_to_point_float(self.md_isolation_core2cond_outer_lineEdit.text())],
                                      cond_cond_isolation=[
                                              comma_str_to_point_float(self.md_isolation_p2p_lineEdit.text()),
                                              comma_str_to_point_float(self.md_isolation_s2s_lineEdit.text()),
                                              comma_str_to_point_float(self.md_isolation_p2s_lineEdit.text())],
                                      conductivity_sigma=[self.md_winding1_material_comboBox.currentText(),
                                                          self.md_winding2_material_comboBox.currentText()])"""

        elif self.md_simulation_type_comboBox.currentText() == 'integrated transformer':
            pass

        return geo

    def md_action_run_simulation(self) -> None:
        """
        Read all input parameters from the fields.
        Run the simulation in dependence of input fields.

        :return: None
        :rtype: None
        """


        geo = self.md_setup_geometry()
        # -----------------------------------------------
        # Simulation
        # -----------------------------------------------
        geo.create_model(freq=comma_str_to_point_float(self.md_base_frequency_lineEdit.text()), visualize_before=False, save_png=False)
        #geo.create_model(freq=comma_str_to_point_float(self.md_base_frequency_lineEdit.text()), visualize_before=False, do_meshing=True, save_png=False)

        winding1_frequency_list, winding1_amplitude_list, winding1_phi_rad_list, winding2_frequency_list, winding2_amplitude_list, winding2_phi_rad_list = self.md_get_frequency_lists()


        if len(winding1_frequency_list) == 1:
            if self.md_simulation_type_comboBox.currentText() == self.translation_dict['inductor']:
                geo.single_simulation(freq=winding1_frequency_list[0],
                                      current=[winding1_amplitude_list[0]],
                                      show_results=True)
            elif self.md_simulation_type_comboBox.currentText() == self.translation_dict['transformer']:
                geo.single_simulation(freq=winding1_frequency_list[0],
                                      current=[winding1_amplitude_list[0], winding2_amplitude_list[0]],
                                      phi_deg=[winding1_phi_rad_list[0], winding2_phi_rad_list[0]])

                                      #phi_deg=[- 1.66257715 / np.pi * 180, 170])


            #geo.single_simulation(freq=winding1_frequency_list[0], current=winding1_amplitude_list)
            #geo.single_simulation(freq=250000, current=[4.14723021, 14.58960019],
                                 # phi_deg=[- 1.66257715 / np.pi * 180, 170])

        else:

            if self.md_simulation_type_comboBox.currentText() == self.translation_dict['inductor']:
                amplitude_list = []
                print(f"{winding1_amplitude_list = }")
                for amplitude_value in winding1_amplitude_list:
                    amplitude_list.append([amplitude_value])

                phase_rad_list = []
                for phase_value in winding1_phi_rad_list:
                    phase_rad_list.append([phase_value])
                geo.excitation_sweep(frequency_list=winding1_frequency_list, current_list_list=amplitude_list, phi_deg_list_list=phase_rad_list)

            elif self.md_simulation_type_comboBox.currentText() == self.translation_dict['transformer']:
                amplitude1_list = []
                for amplitude1_value, amplitude2_value in zip(winding1_amplitude_list, winding2_amplitude_list):
                    amplitude1_list.append([amplitude1_value, amplitude2_value])

                phase1_rad_list = []
                for phase1_value, phase2_value in zip(winding1_phi_rad_list, winding2_phi_rad_list):
                    phase1_rad_list.append([phase1_value, phase2_value])

                geo.excitation_sweep(frequency_list= winding1_frequency_list,
                                     current_list_list=amplitude1_list,
                                     phi_deg_list_list=phase1_rad_list)


            #geo.excitation_sweep(winding1_frequency_list, amplitude_list, phase_rad_list)

        # -----------------------------------------------
        # Read back results
        # -----------------------------------------------

        self.md_simulation_QLabel.setText('simulation fertig.')

        loaded_results_dict = fmt.visualize_simulation_results(geo.e_m_results_log_path, './results.png', show_plot=False)

        pixmap = QPixmap("./results.png")
        self.md_loss_plot_label.setPixmap(pixmap)
        self.md_loss_plot_label.setMask(pixmap.mask())
        self.md_loss_plot_label.show()

        inductance = loaded_results_dict["single_sweeps"][0]["winding1"]["self_inductivity"][0]
        loss_core_eddy_current = loaded_results_dict["total_losses"]["eddy_core"]
        loss_core_hysteresis = loaded_results_dict["total_losses"]["hyst_core_fundamental_freq"]
        loss_winding_1 = loaded_results_dict["total_losses"]["winding1"]["total"]

        self.md_loss_core_hysteresis_label.setText(f"Core Hysteresis loss: {loss_core_hysteresis} W")
        self.md_loss_core_eddy_current_label.setText(f"Core Eddy Current loss: {loss_core_eddy_current} W")
        self.md_loss_winding1_label.setText(f"Winding 1 loss: {loss_winding_1} W")
        self.md_inductance_label.setText(f"Inductance: {inductance} H")

        #log_path = geo.e_m_results_log_path
        # simulation_results = str(fmt.read_results_log(log_path))
        # print(simulation_results)
        # self.md_simulation_output_textBrowser.setText(simulation_results)


    def inductancecalc(self):
        self.core_w = comma_str_to_point_float(self.md_core_width_lineEdit.text())
        self.window_w=comma_str_to_point_float(self.md_window_width_lineEdit.text())
        self.window_h=comma_str_to_point_float(self.md_window_height_lineEdit.text())
        """
        n_turns=int(self.md_winding1_turns_lineEdit.text())
        method=self.md_air_gap_placement_method_comboBox.currentText()
        n_air_gaps=int(self.md_air_gap_count_comboBox.currentText())
        air_gap_h=4
        air_gap_position = 2
        """
        inductance = self.core_w+self.window_w+self.window_h
        self.Inductanceval_label.setText(f"{str(round(inductance,4))} H")


    def therm_simulation(self):
        # Thermal simulation:
        # The losses calculated by the magnetics simulation can be used to calculate the heat distribution of the given magnetic component
        # In order to use the thermal simulation, thermal conductivities for each material can be entered as well as a boundary temperature
        # which will be applied on the boundary of the simulation (dirichlet boundary condition).

        # The case parameter sets the thermal conductivity for a case which will be set around the core.
        # This could model some case in which the transformer is placed in together with a set potting material.
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
            "isolation": 0.42  # polyethylen
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
        geo.create_model(freq=comma_str_to_point_float(self.md_base_frequency_lineEdit.text()), visualize_before=False,
                         save_png=False)
        geo.thermal_simulation(thermal_conductivity_dict, boundary_temperatures, boundary_flags, case_gap_top,
                               case_gap_right, case_gap_bot, True, color_scheme=fmt.colors_ba_jonas,
                               colors_geometry=fmt.colors_geometry_ba_jonas)

        # Because the isolations inside of the winding window are not implemented in femm simulation.
        # The validation only works when the isolations for the FEMMT thermal simulation are turned off.
        #geo.femm_thermal_validation(thermal_conductivity_dict, femm_boundary_temperature, case_gap_top, case_gap_right, case_gap_bot)

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
