# Python standard libraries
import csv
import fileinput
import numpy as np
import os
import sys
import gmsh
import json
import warnings
import inspect
from matplotlib import pyplot as plt
from scipy.optimize import brentq
from scipy.integrate import quad
from typing import List, Dict
from datetime import datetime 

# Third parry libraries
from onelab import onelab

# Local libraries
import femmt.Functions as ff
from femmt.Mesh import Mesh
from femmt.Model import VirtualWindingWindow, WindingWindow, Core, Insulation, StrayPath, AirGaps, Conductor
from femmt.Enumerations import *
from femmt.Data import FileData, MeshData, AnalyticalCoreData
from femmt.TwoDaxiSymmetric import TwoDaxiSymmetric
from femmt.thermal import thermal_simulation

class MagneticComponent:
    """
    A MagneticComponent is the main object for all simulation purposes in femmt.

        - One or more "MagneticComponents" can be created
        - Each "MagneticComponent" owns its own instance variable values

    """
    # Initialization of all class variables
    # Common variables for all instances

    onelab_folder_path = None

    def __init__(self, component_type: ComponentType = ComponentType.Inductor, working_directory: str = None, silent: bool = False, is_gui = False):
        """
        :param component_type: Available options:
                               - "inductor"
                               - "transformer"
                               - "integrated_transformer" (Transformer with included stray-path)
        :type component_type: ComponentType
        :param working_directory: Sets the working directory
        :type working_directory: string
        """
        # Variable to set silent mode
        ff.set_silent_status(silent)

        ff.femmt_print(f"\n"
              f"Initialized a new Magnetic Component of type {component_type.name}\n"
              f"--- --- --- ---")

        # Get caller filepath when no working_directory was set
        if working_directory is None:
            caller_filename = inspect.stack()[1].filename 
            working_directory = os.path.join(os.path.dirname(caller_filename), "femmt")

        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # Create file paths class in order to handle all paths
        self.file_data = FileData(working_directory)

        # Breaking variable. This is set to False when the created model is not drawable (not valid). 
        self.valid = True

        # To make sure femm is only imported once
        self.femm_is_imported = False

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Component Geometry
        self.component_type = component_type  # "inductor", "transformer", "integrated_transformer" (or "three-phase-transformer")

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Components
        self.core = None                        # Contains all informations about the cores
        self.air_gaps = None                    # Contains every air gap
        self.windings = None                    # List of the different winding objects which the following structure: windings[0]: primary, windings[1]: secondary, windings[2]: tertiary ....
        self.insulation = None                   # Contains information about the needed insulations
        self.virtual_winding_windows = None     # Contains a list of every virtual_winding_window which was created
        self.stray_path = None                  # Contains information about the stray_path (only for integrated transformers)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Control Flags
        self.plot_fields = "standard"  # can be "standard" or False

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Excitation Parameters
        # Empty lists will be set when a winding window is added to the magnetic component
        self.imposed_reduced_frequency = None
        self.flag_excitation_type = None
        self.current = []                       # Defined for every conductor
        self.current_density = []               # Defined for every conductor
        self.voltage = []                       # Defined for every conductor
        self.frequency = None
        self.phase_deg = None                   # Default is zero, Defined for every conductor
        self.red_freq = None                    # [] * self.n_windings  # Defined for every conductor
        self.delta = None

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Materials
        self.mu0 = np.pi * 4e-7
        self.e0 = 8.8541878128e-12
        self.Ipeak = None
        # self.ki = None
        # self.alpha = None
        # self.beta = None
        self.t_rise = None
        self.t_fall = None
        self.f_switch = None

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # MeshData to store the mesh size for different points
        # Object is added in set_core
        self.mesh_data = None

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # -- Used for Litz Validation --
        self.sweep_frequencies = None

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Pre-Simulation
        self.reluctance_model = self.ReluctanceModel(self)

        # -- Results --
        self.L_11 = None
        self.L_22 = None
        self.M = None
        self.Pv = None
        # - Primary Concentrated
        self.n_conc = None
        self.L_s_conc = None
        self.L_h_conc = None

        # -- FEMM variables --
        self.tot_loss_femm = None

        self.onelab_setup(is_gui)
        self.onelab_client = onelab.client(__file__)

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   -  -  -  -  -  -  -  -  -  -  -
    # Thermal simulation
    def thermal_simulation(self, thermal_conductivity_dict: Dict, boundary_temperatures_dict: Dict, boundary_flags_dict: Dict, case_gap_top: float,
                           case_gap_right: float, case_gap_bot: float, show_results: bool = True, visualize_before: bool = False, color_scheme: Dict = ff.colors_femmt_default,
                           colors_geometry: Dict = ff.colors_geometry_femmt_default):
        """
        Starts the thermal simulation using thermal_simulation.py

        :param thermal_conductivity_dict: Contains the thermal conductivities for every region
        :type thermal_conductivity_dict: Dict
        :param boundary_temperatures_dict: Contains the tmperatures at each boundary line
        :type boundary_temperatures_dict: Dict
        :param boundary_flags_dict: Sets the boundary type (dirichlet or von neumann) for each boundary line
        :type boundary_flags_dict: Dict
        :param case_gap_top: Size of the top case
        :type case_gap_top: float
        :param case_gap_right: Size of the right case
        :type case_gap_right: float
        :param case_gap_bot: Size of the bot case
        :type case_gap_bot: float
        :param show_results: Shows thermal results in gmsh, defaults to True
        :type show_results: bool, optional
        :param visualize_before: Shows the thermal model before simulation, defaults to False
        :type visualize_before: bool, optional
        :param color_scheme: Color scheme for visualization, defaults to ff.colors_femmt_default
        :type color_scheme: Dict, optional
        :param colors_geometry: Color geometry for visualization, defaults to ff.colors_geometry_femmt_default
        :type colors_geometry: Dict, optional
        """
        # Create necessary folders
        self.file_data.create_folders(self.file_data.thermal_results_folder_path)

        self.mesh.generate_thermal_mesh(case_gap_top, case_gap_right, case_gap_bot, color_scheme, colors_geometry, visualize_before)

        if not os.path.exists(self.file_data.e_m_results_log_path):
            # Simulation results file not created
            raise Exception("Cannot run thermal simulation -> Magnetic simulation needs to run first (no results_log.json found")

        # Check if the results log path simulation settings fit the current simulation settings
        current_settings = MagneticComponent.encode_settings(self)
        del current_settings["working_directory"]
        del current_settings["date"]

        log_settings = None
        with open(self.file_data.e_m_results_log_path, "r") as fd:
            content = json.load(fd)
            log_settings = content["simulation_settings"]
        del log_settings["working_directory"]
        del log_settings["date"]
        
        if current_settings != log_settings:
            raise Exception(f"The settings from the log file {self.file_data.e_m_results_log_path} do not match the current simulation settings. \
                                Please re-run the magnetic simulation.")

        tags = {
            "core_tag": self.mesh.ps_core,
            "background_tag": self.mesh.ps_air,
            "winding_tags": self.mesh.ps_cond,
            "air_gaps_tag": self.mesh.ps_air_gaps if self.air_gaps.number > 0 else None,
            "boundary_regions": self.mesh.thermal_boundary_region_tags,
            "insulations_tag": self.mesh.ps_insulation if len(self.insulation.core_cond) == 4 else None
        }

        # Core area -> Is needed to estimate the heat flux
        # Power density for volumes W/m^3
        core_area = self.calculate_core_volume()

        # Set wire radii
        wire_radii = [winding.conductor_radius for winding in self.windings]

        thermal_parameters = {
            "file_data": self.file_data,
            "tags_dict": tags,
            "thermal_conductivity_dict": thermal_conductivity_dict,
            "boundary_temperatures": boundary_temperatures_dict,
            "boundary_flags": boundary_flags_dict,
            "boundary_physical_groups": {
                "top": self.mesh.thermal_boundary_ps_groups[0],
                "top_right": self.mesh.thermal_boundary_ps_groups[1],
                "right": self.mesh.thermal_boundary_ps_groups[2],
                "bot_right": self.mesh.thermal_boundary_ps_groups[3],
                "bot": self.mesh.thermal_boundary_ps_groups[4]
            },
            "core_area": core_area,
            "conductor_radii": wire_radii,
            "wire_distances": self.get_wire_distances(),
            "case_volume": self.core.r_outer * case_gap_top + self.core.core_h * case_gap_right + self.core.r_outer * case_gap_bot,
            "show_results": show_results,
            "print_sensor_values": False,
            "silent": ff.silent
        }

        thermal_simulation.run_thermal(**thermal_parameters)

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   -  -  -  -  -  -  -  -  -  -  -
    # Setup
    def onelab_setup(self, is_gui):
        """
        Either reads ONELAB parent folder path from config.json or asks the user to provide the ONELAB path it.
        Creates a config.json inside the site-packages folder at first run.
        """
        # check if config.json is available and not empty
        if os.path.isfile(self.file_data.config_path) and os.stat(self.file_data.config_path).st_size != 0:
            onelab_path = ""
            with open(self.file_data.config_path, "r") as fd:
                loaded_dict = json.loads(fd.read())
                onelab_path = loaded_dict['onelab']

            if os.path.exists(onelab_path) and os.path.isfile(os.path.join(onelab_path, "onelab.py")):
                # Path found
                self.file_data.onelab_folder_path = onelab_path
                return

        # Let the user enter the onelab_path:
        # Find out the onelab_path of installed module, or in case of running directly from git, find the onelab_path of git repository
        # loop until path is correct
        onelab_path_wrong = True
        path_wrong = True

        # This is needed because in the gui the input() command cannot be called (it would result in an infinite loop).
        # If the config file was not found just return out of the function. The config file will be added later by the gui handler.
        if is_gui:
            return

        while onelab_path_wrong:
            onelab_path = os.path.normpath(input("Enter the path of onelabs parent folder (path to folder which contains getdp, onelab executables): "))

            if os.path.exists(onelab_path):
                onelab_path_wrong = False
                break
            else:
                ff.femmt_print('onelab not found! Tool searches for onelab.py in the folder. Please re-enter path!')
        self.file_data.onelab_folder_path = onelab_path

        # Write the path to the config.json
        onelab_path_dict = {"onelab": onelab_path}
        with open(os.path.join(self.file_data.config_path), 'w', encoding='utf-8') as fd:
            json.dump(onelab_path_dict, fd, indent=2, ensure_ascii=False)

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   -  -  -  -  -  -  -  -  -  -  -
    # Geometry Parts
    def high_level_geo_gen(self, frequency: float = None, skin_mesh_factor: float = None):
        """ Updates the mesh data and creates the model and mesh objects

        :param frequency: Frequency used in the mesh denisty, defaults to None
        :type frequency: float, optional
        :param skin_mesh_factor: Used in the mesh density, defaults to None
        :type skin_mesh_factor: float, optional
        """
        # Always reset to valid
        self.valid = True

        # Update mesh data
        self.mesh_data.update_data(frequency, skin_mesh_factor)

        # Create model
        self.two_d_axi = TwoDaxiSymmetric(self.core, self.mesh_data, self.air_gaps, self.virtual_winding_windows, self.stray_path,
                                            self.insulation, self.component_type, len(self.windings))
        self.two_d_axi.draw_model()

        # Create mesh
        self.mesh = Mesh(self.two_d_axi, self.windings, self.core.correct_outer_leg, self.file_data, None, ff.silent)

    def mesh(self, frequency: float = None, skin_mesh_factor: float = None):
        """Generates model and mesh.

        :param frequency: Frequency used in the mesh denisty, defaults to None
        :type frequency: float, optional
        :param skin_mesh_factor: Used in the mesh density, defaults to None
        :type skin_mesh_factor: float, optional
        """
        self.high_level_geo_gen(frequency=frequency, skin_mesh_factor=skin_mesh_factor)
        if self.valid:
            self.mesh.generate_hybrid_mesh()
            self.mesh.generate_electro_magnetic_mesh()

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   -  -  -  -  -  -  -  -  -  -  -
    # Create Model
    def set_insulation(self, insulation: Insulation):
        """Adds the insulation to the model

        :param insulation: insulation object
        :type insulation: Insulation
        """
        if insulation.inner_winding_insulations is None or not insulation.inner_winding_insulations:
            raise Exception("insulations between the conductors must be set")

        if insulation.core_cond is None or not insulation.core_cond:
            raise Exception("insulations between the core and the conductors must be set")

        self.insulation = insulation

    def set_stray_path(self, stray_path: StrayPath):
        """Adds the stray path to the model

        :param stray_path: StrayPath object
        :type stray_path: StrayPath
        """
        self.stray_path = stray_path

    def set_air_gaps(self, air_gaps: AirGaps):
        """Adds the air_gaps to the model

        :param air_gaps: AirGaps object
        :type air_gaps: AirGaps
        """
        # Sorting air gaps from lower to upper
        air_gaps.midpoints.sort(key=lambda x: x[1])

        self.air_gaps = air_gaps

    def set_winding_window(self, winding_window: WindingWindow):
        """Adds the virtual winding windows to the model. Creates the windings list, which contains the conductors
        from the virtual winding windows but sorted by the winding_number (ascending).
        Sets empty lists for excitation parameters

        :param winding_window: WindingWindow object
        :type winding_window: WindingWindow
        """
        self.virtual_winding_windows = winding_window.virtual_winding_windows

        windings = []

        for vww in winding_window.virtual_winding_windows:
            if not vww.winding_is_set:
                raise Exception("Each virtual winding window needs to have a winding")
            for winding in vww.windings:
                if winding not in windings:
                    windings.append(winding)

        self.windings = sorted(windings, key = lambda x: x.winding_number)

        # Set excitation parameter lists
        self.current = [None] * len(windings)
        self.current_density = [None] * len(windings)
        self.voltage = [None] * len(windings)
        self.phase_deg = np.zeros(len(windings))

        # Default values for global_accuracy and padding
        self.mesh_data = MeshData(0.5, 1.5, self.mu0, self.core.core_inner_diameter, self.core.window_w, self.windings)

    def set_core(self, core: Core):
        """Adds the core to the model

        :param core: Core object
        :type core: Core
        """
        self.core = core

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   -  -  -  -  -  -  -  -  -  -  -
    # Pre-Processing
    class ReluctanceModel:
        """
        Depending on Core-Configurations, given number of turns and Inductance-goals, calculate air gap lengths
        """

        def __init__(self, component):
            """
            f.e.:   n_reluctances = 1 without stray path
                    n_reluctances = 3 with stray path

                    n_air_gaps = air_gap_division * n_reluctances or sth. like that for distributed air gaps

            """
            self.component = component

            # Control
            self.singularity = False
            self.visualize_waveforms = False
            self.visualize_airgaps = False
            self.visualize_max = False
            self.visualize_nom = False
            self.stray_path_parametrization = None  # "max_flux", "given_flux", "mean_flux"
            self.n_singularities = 0

            # Goals
            self.n_theo = None
            self.L_goal = None

            # Geometric/Material Parameters
            self.N = None
            self.max_length = None
            self.n_ag_per_rel = None
            self.air_gap_types = None
            self.air_gap_lengths = None
            self.material = None
            self.b_peaks = None
            self.b_max = None
            self.b_stray = None
            self.b_stray_rel_overshoot = 1
            self.A_stray = None
            self.A_core = None
            self.real_core_width = None

            # Excitation
            self.f_1st = None
            self.max_current = None
            self.nom_current = None
            self.max_phi = None
            self.nom_phi = None
            self.nom_current_1st = None
            self.nom_phase_1st = None

            # Results
            self.p_hyst_nom_1st = None
            self.p_hyst_nom = None

            self.component.file_data.create_folders(self.component.file_data.reluctance_model_folder_path)

        def calculate_air_gap_lengths_idealized(self, reluctances: List, types: str) -> List:
            """
            :param reluctances: List with reluctances
            :type reluctances: List
            :param types:
                - "round-round":
                - "round-inf":
                - "cyl-cyl":
            :type types: str
            :return: air gap length idealized
            :rtype: List
            """
            air_gap_lengths = [None] * len(reluctances)

            # Go through reluctances
            for n_reluctance, R_0 in enumerate(reluctances):
                if self.component.dimensionality == "2D":

                    if types[n_reluctance] == ("round-round" or "round-inf"):
                        A_core = (self.component.core.core_w / 2) ** 2 * np.pi
                        length = R_0 * self.component.mu0 * A_core

                    if types[n_reluctance] == "cyl-cyl":
                        # return R_0 * self.component.mu0 * w * np.pi * (r_o + r_i)
                        length = None  # TODO: Find explicit formula of idealized air gap length
                air_gap_lengths[n_reluctance] = length

            return air_gap_lengths

        def calculate_air_gap_lengths_with_sct(self, reluctances: List):
            """
            Method calculates air gap lengths according to the given reluctances.
            Method uses several instance variables of the Reluctance Model.

            .. TODO:: List with lists of air gap lengths important to always keep the order of elements [or use a dictionary] future use case integrated_transformer:
                [[l_top_1, l_top_2, ...], [l_bot_1, l_bot_2, ...], [l_stray_1, l_stray_2, ...]]
                use case now to be implemented:  [[l_top_1], [l_bot_1], [l_stray_1]]

            :param reluctances:
            :type reluctances: List
            :return: Dictionary with air gap names and the associated lengths
            """
            air_gap_lengths = {}
            b_peaks = {}

            # Go through reluctances
            for n_reluctance, R_0 in enumerate(reluctances):

                # Define the Reluctance function to be solved for the air gap length
                if self.component.dimensionality == "2D":

                    if n_reluctance == 2:
                        self.A_stray = self.component.stray_path.width * self.component.core.core_w * np.pi
                        A_part = self.A_stray
                        air_gap_name = "R_stray"

                    else:
                        self.A_core = (self.component.core.core_w / 2) ** 2 * np.pi
                        A_part = self.A_core

                        if n_reluctance == 0:
                            air_gap_name = "R_top"

                        if n_reluctance == 1:
                            air_gap_name = "R_bot"

                        # Define tablet height for SCT
                        h_basis = self.max_length[n_reluctance] - self.component.stray_path.width / 2

                    # Check for saturation
                    b_peak = self.max_phi[n_reluctance] / A_part
                    b_peaks[air_gap_name + "_b_peak"] = b_peak

                    # ff.femmt_print(f"{b_peak=}")
                    # ff.femmt_print(f"{A_part=}")

                    if b_peak > self.b_max:
                        air_gap_lengths[air_gap_name] = "saturated"
                        # ff.femmt_print("saturated")
                        # break the outer for loop
                        break

                    # Go through the distributed air gaps (of each reluctance)
                    for n_distributed in range(0, self.n_ag_per_rel[n_reluctance]):

                        if self.air_gap_types[n_reluctance][n_distributed] == "round-inf":
                            # ff.femmt_print(self.component.core.window_h,
                            #       self.component.stray_path.midpoint,
                            #       self.component.stray_path.width)

                            def r_sct(length):
                                return ff.r_round_inf(l=length,
                                                   r=self.component.core.core_w / 2,
                                                   sigma=ff.sigma(l=length,
                                                               w=self.component.core.core_w / 2,
                                                               R_equivalent=ff.r_basis(l=length,
                                                                                    w=self.component.core.core_w,
                                                                                    h=h_basis))
                                                   ) - R_0

                            if self.visualize_airgaps:
                                def r_sct_ideal(length):
                                    return ff.r_round_inf(l=length,
                                                       r=self.component.core.core_w / 2,
                                                       sigma=1
                                                       ) - R_0

                        if self.air_gap_types[n_reluctance][n_distributed] == "round-round":
                            # def r_sct(length):
                            #    return r_round_round(length)
                            pass

                        if self.air_gap_types[n_reluctance][n_distributed] == "cyl-cyl":
                            def r_sct(length):
                                return ff.r_cyl_cyl(l=length,
                                                 sigma=ff.sigma(l=length,
                                                             w=self.component.stray_path.width,
                                                             R_equivalent=ff.r_basis(l=length,
                                                                                  w=self.component.stray_path.width,
                                                                                  h=self.component.core.window_w
                                                                                    - length) / 2),
                                                 w=self.component.stray_path.width,
                                                 r_o=self.component.core.window_w + self.component.core.core_w / 2
                                                 ) - R_0

                            def r_sct_real(length):
                                return ff.r_cyl_cyl_real(l=length,
                                                      sigma=ff.sigma(l=length,
                                                                  w=self.component.stray_path.width,
                                                                  R_equivalent=ff.r_basis(l=length,
                                                                                       w=self.component.stray_path.width,
                                                                                       h=self.component.core.window_w
                                                                                         - length) / 2),
                                                      w=self.component.stray_path.width,
                                                      r_o=self.component.core.window_w + self.component.core.core_w / 2,
                                                      h_real_core=self.real_core_width
                                                      ) - R_0

                        # ff.femmt_print(f"\n  {air_gap_name}")
                        # ff.femmt_print(f"n_reluctance {n_reluctance}")
                        # ff.femmt_print(f"self.component.stray_path.width {self.component.stray_path.width}")
                        # ff.femmt_print(f"max_length[n_reluctance] {self.max_length[n_reluctance]}")
                        # ff.femmt_print(f"R_0 {R_0}")
                        # ff.femmt_print(f"r_sct(a) {r_sct(1e-6)}")
                        # ff.femmt_print(f"r_sct(b) {r_sct(self.max_length[n_reluctance])}")

                        # Check for different signs (zero crossing)
                        if r_sct(1e-6) * r_sct(self.max_length[n_reluctance]) > 0:
                            air_gap_lengths[air_gap_name] = "out of bounds"
                            # ff.femmt_print("out of bounds")

                        else:
                            if self.visualize_airgaps:
                                length_vec = np.linspace(1e-6, self.max_length[n_reluctance] / 4)
                                R_res = np.array([r_sct(length) + R_0 for length in length_vec])
                                R_res_ideal = np.array([r_sct_ideal(length) + R_0 for length in length_vec])
                                R_0_vec = np.array([R_0 for length in length_vec])
                                length_vec = length_vec * 1000

                                plt.figure(figsize=(5, 3))
                                plt.plot(length_vec, R_res * 1e-6, label=r"$R_{\mathrm{fringing}}$")
                                plt.plot(length_vec, R_res_ideal * 1e-6, label=r"$R_{\mathrm{ideal}}$")
                                plt.plot(length_vec, R_0_vec * 1e-6, label=r"$R_{\mathrm{goal}}$")
                                plt.xlabel(r"$l / \mathrm{mm}$")
                                plt.ylabel(r"$R / \mathrm{(\mu H)^{-1}}$")
                                plt.legend()
                                plt.grid()
                                plt.savefig(f"C:/Users/tillp/sciebo/Exchange Till/04_Documentation/Inkscape/Reluctance_Model/sct_resistance.pdf",
                                            bbox_inches="tight")

                                plt.show()

                            air_gap_lengths[air_gap_name] = brentq(r_sct, 1e-6, self.max_length[n_reluctance])
                            if air_gap_name == "R_stray":
                                air_gap_lengths[air_gap_name + "_real"] = brentq(r_sct_real, 1e-9, self.max_length[n_reluctance])

            # ff.femmt_print(f"{air_gap_lengths=}")
            return air_gap_lengths, b_peaks

        def get_core_loss(self):
            """
            Calculates the hysteresis loss corresponding to complex core parameters.
            """
            # ff.femmt_print(f"Calculate Analytical Core Losses:")

            # time domain
            Phi_top_nom, Phi_bot_nom, Phi_stray_nom = \
                self.phi_from_time_currents(self.nom_current[0], self.nom_current[1], self.visualize_nom)

            # Peaks of time domain
            Phi_top_nom_peak, Phi_bot_nom_peak, Phi_stray_nom_peak = \
                self.max_phi_from_time_phi(Phi_top_nom, Phi_bot_nom, Phi_stray_nom)

            # Phases of time domain
            phase_top, phase_bot, phase_stray = self.phases_from_time_phi(Phi_top_nom, Phi_bot_nom, Phi_stray_nom)

            # fundamental
            Phi_top_nom_1st, Phi_bot_nom_1st, Phi_stray_nom_1st = self.phi_fundamental()

            # Amplitudes of fundamental
            Phi_top_nom_1st_peak, Phi_bot_nom_1st_peak, Phi_stray_nom_1st_peak = \
                self.max_phi_from_time_phi(Phi_top_nom_1st, Phi_bot_nom_1st, Phi_stray_nom_1st)

            p_top_1st, p_bot_1st, p_stray_1st = self.hysteresis_loss(Phi_top_nom_1st_peak, Phi_bot_nom_1st_peak, Phi_stray_nom_1st_peak)
            # ff.femmt_print(p_top_1st, p_bot_1st, p_stray_1st)
            self.p_hyst_nom_1st = sum([p_top_1st, p_bot_1st, p_stray_1st])

            p_top, p_bot, p_stray = self.hysteresis_loss(Phi_top_nom_peak, Phi_bot_nom_peak, Phi_stray_nom_peak)
            # ff.femmt_print(p_top, p_bot, p_stray)
            self.p_hyst_nom = sum([p_top, p_bot, p_stray])

            self.visualize_phi_core_loss(Phi_top_nom, Phi_bot_nom, Phi_stray_nom,
                                         phase_top, phase_bot, phase_stray)

            # ff.femmt_print(f"Analytical Core Losses = \n\n")

        def visualize_phi_core_loss(self, Phi_top_init, Phi_bot_init, Phi_stray_init,
                                    phase_top, phase_bot, phase_stray):
            """
            Visualization of the fluxes used for the core loss calculation.
            """
            Phi_top_nom_peak, Phi_bot_nom_peak, Phi_stray_nom_peak = \
                self.max_phi_from_time_phi(Phi_top_init, Phi_bot_init, Phi_stray_init)

            # time
            t = np.linspace(min(self.nom_current[0][0]), max(self.nom_current[0][0]), 500) / self.f_1st / 2 / np.pi
            # phase = 0
            Phi_top = Phi_top_nom_peak * np.cos(t * self.f_1st * 2 * np.pi - phase_top)
            Phi_bot = Phi_bot_nom_peak * np.cos(t * self.f_1st * 2 * np.pi - phase_bot)
            Phi_stray = Phi_stray_nom_peak * np.cos(t * self.f_1st * 2 * np.pi - phase_stray)

            if self.visualize_nom:
                # figure, axis = plt.subplots(3, figsize=(4, 6))

                # t = np.array(t) * 10 ** 6 / 200000 / 2 / np.pi
                # for i, phi in enumerate([Phi_top, Phi_bot, Phi_stray]):
                #     axis[i].plot(t, 1000 * np.array(phi), label=r"$\mathcal{\phi}_{\mathrm{top}}$")
                #     axis[i].set_ylabel("Magnetic fluxes / mWb")
                #     axis[i].set_xlabel(r"$t$ / \mu s")
                #     axis[i].legend()
                #     axis[i].grid()

                plt.figure(figsize=(6, 3))
                plt.plot(t, 1000 * np.array(Phi_top), color="tab:blue", label=r"$\mathcal{\phi}_{\mathrm{top, fd}}$")
                plt.plot(t, 1000 * np.array(Phi_bot), color="tab:orange", label=r"$\mathcal{\phi}_{\mathrm{bot, fd}}$")
                plt.plot(t, 1000 * np.array(Phi_stray), color="tab:green", label=r"$\mathcal{\phi}_{\mathrm{stray, fd}}$")

                time = np.array(self.nom_current[0][0]) / 200000 / 2 / np.pi

                plt.plot(time, 1000 * np.array(Phi_top_init), ":", color="tab:blue", label=r"$\mathcal{\phi}_{\mathrm{top, td}}$")
                plt.plot(time, 1000 * np.array(Phi_bot_init), ":", color="tab:orange", label=r"$\mathcal{\phi}_{\mathrm{bot, td}}$")
                plt.plot(time, 1000 * np.array(Phi_stray_init), ":", color="tab:green", label=r"$\mathcal{\phi}_{\mathrm{stray, td}}$")

                plt.ylabel("Magnetic fluxes / mWb")
                plt.xlabel(r"$t$ / s")
                plt.yticks([-0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03])
                plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
                plt.grid()
                # ToDo: general path
                plt.savefig(f"C:/Users/tillp/sciebo/Exchange Till/04_Documentation/Reluctance_Model_Current_Shapes/core_loss.pdf",
                            bbox_inches="tight")

                plt.show()

        def phi_fundamental(self):
            """TODO Doc
            """
            # time
            t = np.linspace(min(self.nom_current[0][0]), max(self.nom_current[0][0]), 500) / self.f_1st / 2 / np.pi

            # time functions
            # f_I1 = interp1d(self.nom_current[0][0], self.nom_current[0][1])  # current_pairs_nom[0]
            # f_I2 = interp1d(self.nom_current[1][0], self.nom_current[1][1])  # current_pairs_nom[1]

            # 1st harmonic

            # Imax_in = max(time_nom_in[1])
            # I1 = Imax_in * np.cos(t * f + phase_pairs_nom[0][0])
            I1 = self.nom_current_1st[0] * np.cos(t * self.f_1st * 2 * np.pi + self.nom_phase_1st[0])

            # Imax_out = max(time_nom_out[1])
            # I2 = Imax_out * np.cos(t * f + phase_pairs_nom[0][1]+np.pi)
            I2 = self.nom_current_1st[1] * np.cos(t * self.f_1st * 2 * np.pi + self.nom_phase_1st[1] + np.pi)

            # Calc fluxes
            Phi_top = []
            Phi_bot = []
            Phi_stray = []

            for i in range(0, len(I1)):
                CurrentVector = [I1[i], I2[i]]

                #     [Phi_top, Phi_bot] = np.matmul(np.matmul(np.linalg.inv(R), N), self.current)
                [Phi_top_s, Phi_bot_s] = np.matmul(np.matmul(np.linalg.inv(np.transpose(self.N)), self.L_goal),
                                                   np.transpose(CurrentVector))
                Phi_stray_s = Phi_bot_s - Phi_top_s

                # Append to lists (for all time-steps)
                Phi_top.append(Phi_top_s)
                Phi_bot.append(Phi_bot_s)
                Phi_stray.append(Phi_stray_s)

            if self.visualize_nom:
                self.visualize_current_and_flux(t, Phi_top, Phi_bot, Phi_stray, I1, I2)

            return Phi_top, Phi_bot, Phi_stray

        def hysteresis_loss(self, Phi_top, Phi_bot, Phi_stray):
            """
            Returns the hysteresis losses. Assumptions: homogeneous flux, sinusoidal excitation
            """
            length_corner = self.A_core / self.component.core.core_w / np.pi
            ri = self.component.core.core_w / 2
            ro = ri + self.component.core.window_w

            # Top Part
            b_top = Phi_top / self.A_core
            Vol_top = self.A_core * (100 - self.component.stray_path.midpoint) / 100 * self.component.core.window_h + \
                      self.A_core * ((100 - self.component.stray_path.midpoint) / 100 * self.component.core.window_h - \
                                     self.air_gap_lengths["R_top"])

            p_top = 0.5 * self.component.mu0 * AnalyticalCoreData.f_N95_mu_imag(self.f_1st, b_top) * Vol_top * 2 * np.pi * self.f_1st * \
                    (b_top / self.component.core.mu_rel / self.component.mu0) ** 2 + \
                    self.p_loss_cyl(Phi_top, length_corner, ri, ro)[0]

            # Bot Part
            b_bot = Phi_bot / self.A_core
            Vol_bot = self.A_core * self.component.stray_path.midpoint / 100 * self.component.core.window_h + \
                      self.A_core * (self.component.stray_path.midpoint / 100 * self.component.core.window_h + \
                                     self.air_gap_lengths["R_bot"])
            p_bot = 0.5 * self.component.mu0 * AnalyticalCoreData.f_N95_mu_imag(self.f_1st, b_bot) * Vol_bot * 2 * np.pi * self.f_1st * \
                    (b_bot / self.component.core.mu_rel / self.component.mu0) ** 2 + \
                    self.p_loss_cyl(Phi_bot, length_corner, ri, ro)[0]

            # Stray Path
            p_stray = self.p_loss_cyl(Phi_stray, self.component.stray_path.width, ri,
                                      ro - self.air_gap_lengths["R_stray"])[0]

            return p_top, p_bot, p_stray

        def p_loss_cyl(self, Phi, w, ri, ro):
            """TODO Doc
            """

            def b(r, Phi, w):
                return Phi / (2 * np.pi * r * w)

            def p_loss_density(r, Phi, w):
                return 2 * np.pi * r * w * \
                       np.pi * self.f_1st * \
                       self.component.mu0 * AnalyticalCoreData.f_N95_mu_imag(self.f_1st, b(r, Phi, w)) * \
                       (b(r, Phi, w) / self.component.core.mu_rel / self.component.mu0) ** 2

            return quad(p_loss_density, ri, ro, args=(Phi, w), epsabs=1e-4)

        @staticmethod
        def max_phi_from_time_phi(Phi_top, Phi_bot, Phi_stray):
            """TODO Doc
            :param Phi_stray:
            :param Phi_bot:
            :param Phi_top:
            """
            Phi_top_peak = max([abs(ele) for ele in Phi_top])
            Phi_bot_peak = max([abs(ele) for ele in Phi_bot])
            Phi_stray_peak = max([abs(ele) for ele in Phi_stray])

            # ff.femmt_print(f"{Phi_top_peak=}\n"
            #       f"{Phi_bot_peak=}\n"
            #       f"{Phi_stray_peak=}\n")

            return Phi_top_peak, Phi_bot_peak, Phi_stray_peak

        def phases_from_time_phi(self, Phi_top, Phi_bot, Phi_stray):
            """TODO DOc
            Returns the phases of the peaks
            :param Phi_stray:
            :param Phi_bot:
            :param Phi_top:
            """
            # ff.femmt_print(np.array(Phi_top))
            # ff.femmt_print(np.array(Phi_top).argmax(axis=0))
            # ff.femmt_print(self.nom_current[0][0])

            phase_top = self.nom_current[0][0][np.array(Phi_top).argmax(axis=0)]
            phase_bot = self.nom_current[0][0][np.array(Phi_bot).argmax(axis=0)]
            phase_stray = self.nom_current[0][0][np.array(Phi_stray).argmax(axis=0)]

            # ff.femmt_print(np.array(Phi_top).argmax(axis=0))
            # ff.femmt_print(np.array(Phi_bot).argmax(axis=0))
            # ff.femmt_print(np.array(Phi_stray).argmax(axis=0))

            return phase_top, phase_bot, phase_stray

        def phi_from_time_currents(self, time_current_1, time_current_2, visualize=False):
            """TODO Doc
            :param visualize:
            :param time_current_2:
            :param time_current_1:
            """
            I1 = time_current_1[1]
            I2 = time_current_2[1]

            t = None
            if time_current_1[0] == time_current_2[0]:
                t = time_current_1[0]
            else:
                warnings.warn("Time values do not match equally for both currents!")

            Phi_top = []
            Phi_bot = []
            Phi_stray = []

            for i in range(0, len(I1)):
                # Negative sign is placed here
                CurrentVector = [I1[i], -I2[i]]

                #     [Phi_top, Phi_bot] = np.matmul(np.matmul(np.linalg.inv(R), N), self.max_current)
                [Phi_top_s, Phi_bot_s] = np.matmul(np.matmul(np.linalg.inv(np.transpose(self.N)), self.L_goal),
                                                   np.transpose(CurrentVector))
                Phi_stray_s = Phi_bot_s - Phi_top_s

                # Append to lists (for all time-steps)
                Phi_top.append(Phi_top_s)
                Phi_bot.append(Phi_bot_s)
                Phi_stray.append(Phi_stray_s)

            # Visualize
            if visualize:
                self.visualize_current_and_flux(t, Phi_top, Phi_bot, Phi_stray, I1, I2)

            return Phi_top, Phi_bot, Phi_stray

        @staticmethod
        def visualize_current_and_flux(t, Phi_top, Phi_bot, Phi_stray, I1, I2):
            figure, axis = plt.subplots(2, figsize=(4, 4))

            t = np.array(t) * 10 ** 6 / 200000 / 2 / np.pi

            axis[0].plot(t, 1000 * np.array(Phi_top), label=r"$\mathcal{\phi}_{\mathrm{top}}$")
            axis[0].plot(t, 1000 * np.array(Phi_bot), label=r"$\mathcal{\phi}_{\mathrm{bot}}$")
            axis[0].plot(t, 1000 * np.array(Phi_stray), label=r"$\mathcal{\phi}_{\mathrm{stray}}$")

            axis[0].set_yticks([-0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03])

            axis[1].plot(t, I1, label=r"$I_{\mathrm{in}}$")
            axis[1].plot(t, -np.array(I2), label=r"$I_{\mathrm{out}}$")
            # axis[1].plot(t, np.array(I2)/3.2 - np.array(I1), label=f"Im")

            axis[0].set_ylabel("Magnetic fluxes / mWb")
            axis[0].set_xlabel(r"$t$ / \mu s")
            axis[1].set_ylabel("Currents / A")
            axis[1].set_xlabel(r"$t$ / $\mathrm{\mu s}$")

            axis[0].legend()
            axis[1].legend()

            axis[0].grid()
            axis[1].grid()
            #ToDo: general path
            plt.savefig(f"C:/Users/tillp/sciebo/Exchange Till/04_Documentation/Reluctance_Model_Current_Shapes/{I1[0]}.pdf",
                        bbox_inches="tight")

            # plt.show()

        def stray_path_parametrization_two_d_axi(self):
            """ TODO Doc
            Method defines instance variables.
            """
            # Calculate stray path width

            if np.linalg.cond(np.transpose(self.N)) < 1 / sys.float_info.epsilon:

                # Get fluxes for max current
                Phi_top_max, Phi_bot_max, Phi_stray_max = self.phi_from_time_currents(self.max_current[0],
                                                                                      self.max_current[1],
                                                                                      visualize=self.visualize_max)

                Phi_top_max_peak, Phi_bot_max_peak, Phi_stray_max_peak = \
                    self.max_phi_from_time_phi(Phi_top_max, Phi_bot_max, Phi_stray_max)

                # Store max peak fluxes in instance variable -> used for saturation check
                self.max_phi = [Phi_top_max_peak, Phi_bot_max_peak, Phi_stray_max_peak]
                # ff.femmt_print(self.max_phi)

                # self.component.stray_path.width = (Phi_stray+np.abs(Phi_top + Phi_bot)) / self.b_stray /
                # (np.pi * self.component.core.core_w)

                # Calculate stray path width corresponding to ...
                # ... externally set magnetic flux density or ...
                if self.stray_path_parametrization == "given_flux":
                    self.component.stray_path.width = Phi_stray_max_peak / self.b_stray / (np.pi * self.component.core.core_w)
                # ... to mean flux density of top and bottom flux density
                if self.stray_path_parametrization == "mean":
                    b_stray_mean = (np.abs(Phi_top_max_peak) + np.abs(Phi_bot_max_peak)) / 2 / np.pi / (self.component.core.core_w / 2) ** 2
                    self.component.stray_path.width = Phi_stray_max_peak / b_stray_mean / (np.pi * self.component.core.core_w)
                if self.stray_path_parametrization == "max_flux":
                    phi_abs = np.array([np.abs(Phi_top_max_peak), np.abs(Phi_bot_max_peak)])
                    b_stray_max = phi_abs.max(0) / np.pi / (self.component.core.core_w / 2) ** 2
                    self.component.stray_path.width = Phi_stray_max_peak / b_stray_max / self.b_stray_rel_overshoot \
                                                      / (np.pi * self.component.core.core_w)

                # Max allowed lengths in each leg
                # TODO: Fit values a bit
                self.max_length = [self.component.core.window_h * (100 - self.component.stray_path.midpoint) / 100,
                                   self.component.core.window_h * self.component.stray_path.midpoint / 100,
                                   self.component.core.window_w / 2]

                # Air gap types: "round-round", "round-inf", "cyl-cyl"
                self.air_gap_types = [["round-inf"],  # top
                                      ["round-inf"],  # bot
                                      ["cyl-cyl"]]  # stray

                # TODO: R_top and R_bot can be realized with distributed air gaps
                self.n_ag_per_rel = [1, 1, 1]

            else:
                # Singular Matrices cannot be treated
                self.n_singularities += 1
                self.singularity = True
                # ff.femmt_print(N)
                [Phi_top_max_peak, Phi_bot_max_peak, Phi_stray_max_peak] = [0, 0, 0]  # TODO:Case treatment

        def get_air_gaps_from_winding_matrix(self):
            """
            Calculates air gap lengths according to a given winding matrix.
            Uses several instance variables.

            :return: Dictionary with air gap results or None
            """
            results = None
            # Core and Component type decide about air gap characteristics
            if self.component.dimensionality == "2D":

                # Check winding matrix

                # Inductor
                if self.N.shape == (1,) and self.component.component_type == ComponentType.Inductor:
                    raise NotImplemented

                # Transformer
                if self.N.shape == (2,) and self.component.component_type == ComponentType.Transformer:
                    raise NotImplemented

                # Dedicated stray path
                if self.N.shape == (2, 2) and self.component.component_type == ComponentType.IntegratedTransformer:

                    # Calculate goal reluctance matrix
                    R_matrix = ff.calculate_reluctances(N=self.N, L=self.L_goal)

                    # R_goal = [R_top, R_bot, R_stray]
                    R_goal = [R_matrix[0, 0] + R_matrix[0, 1], R_matrix[1, 1] + R_matrix[0, 1], -R_matrix[0, 1]]

                    # Stray path specific parameters
                    # ff.femmt_print(R_matrix)
                    # ff.femmt_print(R_goal)
                    # ff.femmt_print(self.L_goal)
                    self.stray_path_parametrization_two_d_axi()

                    # Check for negative Reluctances
                    if all(R >= 0 for R in R_goal) and self.singularity == False:

                        # Calculate the air gap lengths with the help of SCT
                        # air_gap_lengths is a dictionary with air gap names and the associated length
                        self.air_gap_lengths, self.b_peaks = self.calculate_air_gap_lengths_with_sct(reluctances=R_goal)

                        # ff.femmt_print(air_gap_lengths.values())

                        # Check for invalid data
                        if self.air_gap_lengths.values() in ['saturated', 'out of bounds']:

                            results = None
                            # ff.femmt_print("Invalid Data\n\n")

                        else:
                            # Width of the stray path is added to the result data
                            stray_path_width = {"stray_path_width": self.component.stray_path.width}

                            # Put together the single dictionaries
                            lengths_and_peaks = dict(self.b_peaks, **self.air_gap_lengths)
                            results = dict(stray_path_width, **lengths_and_peaks)

                    else:
                        results = None

                # Save resulting parameter set
                # air_gap_lengths_valid = np.asarray(air_gap_lengths_valid)
                # np.save('Reluctance_Model/air_gap_lengths_valid.npy', air_gap_lengths_valid)

                # ff.femmt_print(f"{results=}")
                return results

        def air_gap_design(self, L_goal, parameters_init, stray_path_parametrization=None, f_1st=None,
                           max_current=None, nom_current=None, b_max=None, b_stray=None, material=None,
                           visualize_waveforms=False, nom_current_1st=None, nom_phase_1st=None):
            """
            Performs calculation of air gap lengths according to given data.

            :param f_1st:
            :param nom_phase_1st:
            :param nom_current_1st:
            :param max_current:
            :param nom_current:
            :param material:
            :param visualize_waveforms:
            :param stray_path_parametrization: the width of the stray path may be parametrized by "max_flux" or
                "mean_flux" of the other two legs (top and bot) or by a "given_flux" which can be defined with param:b_stray
            :param b_max:
            :param parameters_init:
            :param b_stray:
            :param L_goal: list of inductance goals [inductor: single value L;
                                                     transformer: Inductance Matrix [[L_11, M], [M, L_22]]

            :return: Dictionary with resulting air gap (and stray path) parameters. Saturated and other invalid
                        parameter combinations are set to None
            """
            # Excitation
            self.f_1st = f_1st
            self.max_current = max_current
            self.nom_current = nom_current
            self.nom_current_1st = nom_current_1st
            self.nom_phase_1st = nom_phase_1st


            # Others
            self.stray_path_parametrization = stray_path_parametrization
            self.b_stray = b_stray
            self.b_max = b_max
            self.material = material

            # Visualization
            self.visualize_waveforms = visualize_waveforms
            if self.visualize_waveforms == "max":
                self.visualize_max = True
                self.visualize_nom = False
            if self.visualize_waveforms == "nom":
                self.visualize_max = False
                self.visualize_nom = True
            if self.visualize_waveforms == "all":
                self.visualize_max = True
                self.visualize_nom = True

            if not os.path.exists(self.reluctance_model_folder):
                os.mkdir(self.reluctance_model_folder)

            # Save initial parameter set
            # parameters_init = np.asarray(parameters_init)
            # np.save('Reluctance_Model/parameters_init.npy', parameters_init)

            # Save goal inductance values
            self.L_goal = np.asarray(L_goal)
            np.save(self.component.path + '/Reluctance_Model/goals.npy', self.L_goal)

            # Initialize result list
            parameter_results = []

            # Put to Terminal
            ff.femmt_print(f"\n"
                  f"--- ---\n"
                  f"Perform reluctance calculations\n\n"
                  # f"Goal inductance values                   : {L_goal}\n\n"
                  f"Number of initial reluctance parameters  : {len(parameters_init)}\n")

            # Update the core to use its internal core parameter calculation functionality
            # Set attributes of core with given keywords
            for index, parameters in enumerate(parameters_init):
                for key, value in parameters.items():
                    # ff.femmt_print(key, value)
                    if hasattr(self, key):
                        setattr(self, key, value)

                    if hasattr(self.component.core, key):
                        setattr(self.component.core, key, value)

                    if hasattr(self.component.stray_path, key):
                        setattr(self.component.stray_path, key, value)

                    if key == "N":
                        self.N = value

                self.singularity = False
                air_gap_results = self.get_air_gaps_from_winding_matrix()

                # Check if the result is valid
                if air_gap_results is None:
                    all_parameters = None
                else:
                    # Add frequency to results
                    wp_frequency = {"frequency": self.f_1st}
                    model_results = dict(air_gap_results, **wp_frequency)

                    # Calculate analytical Hysteresis Loss
                    self.get_core_loss()
                    losses = {"p_hyst_nom_1st": self.p_hyst_nom_1st, "p_hyst_nom": self.p_hyst_nom}

                    model_results = dict(model_results, **losses)

                    all_parameters = dict(parameters, **model_results)

                parameter_results.append(all_parameters)

            ff.femmt_print(f"Number of singularities: {self.n_singularities}\n")

            # Save Results including invalid parameters
            # ff.femmt_print(f"{parameter_results=}")

            np.save(self.component.path + '/Reluctance_Model/parameter_results.npy', parameter_results)

            # Filter all entries, that are None
            # Elements of list reluctance_parameters are either a dictionary or None
            valid_parameter_results = [x for x in parameter_results if x is not None]

            # Save resulting valid parameters
            np.save(self.component.path + '/Reluctance_Model/valid_parameter_results.npy', valid_parameter_results)

            ff.femmt_print(f"Number of valid parameters: {len(valid_parameter_results)}\n\n"
                  f"Ready with reluctance calculations\n"
                  f"--- ---\n")

            return valid_parameter_results

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    # Standard Simulations
    def create_model(self, freq: float, skin_mesh_factor: float = 0.5, visualize_before: bool = False,
                     save_png: bool = False, color_scheme: Dict = ff.colors_femmt_default,
                     colors_geometry: Dict = ff.colors_geometry_femmt_default):
        """
        Create a model from the abstract geometry description inside onelab including optional mesh generation

        :param freq: Frequency [Hz]
        :type freq: float
        :param skin_mesh_factor: [default to 0.5]
        :type skin_mesh_factor: float
        :param visualize_before: True for a pre-visualisation (e.g. check your geometry) and after this a simulation runs, False for a direct simulation
        :type visualize_before: bool
        :param: do_meshing: [default to True], internal use only (e.g. for GUI)
        :type do_meshing: bool
        :param save_png: True to save a png-figure, false for no figure
        :type save_png: bool
        :param color_scheme: colorfile (definition for red, green, blue, ...)
        :type color_scheme: Dict
        :param colors_geometry: definition for e.g. core is grey, winding is orange, ...
        :type colors_geometry: Dict
        """
        if self.core is None:
            raise Exception("A core class needs to be added to the magnetic component")
        if self.air_gaps is None:
            self.air_gaps = AirGaps(None, None)
            ff.femmt_print("No air gaps are added")
        if self.insulation is None:
            raise Exception("An insulation class need to be added to the magnetic component")
        if self.virtual_winding_windows is None:
            raise Exception("Virtual winding windows are not set properly. Please check the winding creation")

        self.high_level_geo_gen(frequency=freq, skin_mesh_factor=skin_mesh_factor)
        if self.valid:
            self.mesh.generate_hybrid_mesh(visualize_before=visualize_before, save_png=save_png, color_scheme=color_scheme, colors_geometry=colors_geometry)
        else:
            raise Exception("The model is not valid. The simulation won't start.")

    def pre_simulation(self):
        """
        - Complete "pre-simulation" call
        """
        self.high_level_geo_gen()
        self.excitation(frequency=100000, amplitude_list=1)  # arbitrary values: frequency and current
        self.file_communication()
        self.pre_simulate()

    def single_simulation(self, freq: float, current: List[float], phi_deg: List[float] = None, show_results = True):
        """

        Start a _single_ electromagnetic ONELAB simulation.
        :param NL_core:
        :param freq: frequency to simulate
        :type freq: float
        :param current: current to simulate
        :param skin_mesh_factor:
        :type skin_mesh_factor: float
        :param phi_deg: phase angle in degree
        :type phi_deg: List[float]
        """
        phi_deg = phi_deg or []

        self.mesh.generate_electro_magnetic_mesh()
        self.excitation(frequency=freq, amplitude_list=current, phase_deg_list=phi_deg)  # frequency and current
        self.file_communication()
        self.pre_simulate()
        self.simulate()
        self.calculate_and_write_log()
        if show_results:
            self.visualize()
        # results =
        # return results

    def excitation_sweep(self, frequency_list: List, current_list_list: List, phi_deg_list_list: List,
                         show_last: bool = False, return_results: bool = False, 
                         excitation_meshing_type: ExcitationMeshingType = None, skin_mesh_factor: float = 0.5, visualize_before: bool = False, save_png: bool = False,
                         color_scheme: Dict = ff.colors_femmt_default, colors_geometry: Dict = ff.colors_geometry_femmt_default) -> Dict:
        """
        Performs a sweep simulation for frequency-current pairs. Both values can
        be passed in lists of the same length. The mesh is only created ones (fast sweep)!

        :Example Code for Inductor:

        >>> import femmt as fmt
        >>> fs_list = [0, 10000, 30000, 60000, 100000, 150000]
        >>> amplitue_list_list = [[10], [2], [1], [0.5], [0.2], [0.1]]
        >>> phase_list_list = [[0], [10], [20], [30], [40], [50]]
        >>> geo.excitation_sweep(frequency_list=fs_list, current_list_list=amplitue_list_list, phi_deg_list_list=phase_list_list)

        :Example Code for Transformer with 2 windings:

        >>> import femmt as fmt
        >>> fs_list = [0, 10000, 30000, 60000, 100000, 150000]
        >>> amplitue_list_list = [[10, 2], [2, 1], [1, 0.5], [0.5, 0.25], [0.2, 0.1], [0.1, 0.05]]
        >>> phase_list_list = [[0, 170], [10, 180], [20, 190], [30, 200], [40, 210], [50, 220]]
        >>> geo.excitation_sweep(frequency_list=fs_list, current_list_list=amplitue_list_list, phi_deg_list_list=phase_list_list)

        :param frequency_list: Frequency in a list
        :type frequency_list: List
        :param current_list_list: current amplitude, must be a list in a list, see example!
        :type current_list_list: List
        :param phi_deg_list_list: phase in degree, must be a list in a list, see example!
        :type phi_deg_list_list: List
        :param show_last: shows last simulation in gmsh if set to True
        :type show_last: bool
        :param return_results: returns results in a dictionary
        :type return_results: bool
        :param visualize_before: show genarated mesh before the simulation is run
        :type visualize_before: bool
        :param meshing:
        :type meshing: bool
        :param color_scheme: colorfile (definition for red, green, blue, ...)
        :type color_scheme: Dict
        :param colors_geometry: definition for e.g. core is grey, winding is orange, ...
        :type colors_geometry: Dict


        :return: Results in a dictionary
        :rtype: Dict
        """
        # frequencies = frequencies or []
        # currents = currents or []
        # phi = phi or []
        if show_last:
            self.plot_fields = "standard"
        else:
            self.plot_fields = False

        # If one conductor is solid and no meshing type is given then change the meshing type to MeshEachFrequency
        # In case of litz wire, only the lowest frequency is meshed (frequency indepent due to litz-approximation)
        if excitation_meshing_type is None:
            for winding in self.windings:
                if winding.conductor_type == ConductorType.RoundSolid:
                    excitation_meshing_type = ExcitationMeshingType.MeshEachFrequency
                    break
                if winding.conductor_type == ConductorType.RoundLitz:
                    excitation_meshing_type = ExcitationMeshingType.MeshOnlyLowestFrequency

        if excitation_meshing_type == ExcitationMeshingType.MeshEachFrequency:
            for i in range(0, len(frequency_list)):
                self.high_level_geo_gen(frequency=frequency_list[i], skin_mesh_factor=skin_mesh_factor)
                if self.valid:
                    self.mesh.generate_hybrid_mesh(color_scheme, colors_geometry, visualize_before=visualize_before, save_png=save_png)
                    self.mesh.generate_electro_magnetic_mesh()
                
                self.excitation(frequency=frequency_list[i], amplitude_list=current_list_list[i],
                                    phase_deg_list=phi_deg_list_list[i])  # frequency and current
                self.file_communication()
                self.pre_simulate()
                self.simulate()
        else:
            if excitation_meshing_type == ExcitationMeshingType.MeshOnlyHighestFrequency:
                self.high_level_geo_gen(frequency=max(frequency_list), skin_mesh_factor=skin_mesh_factor)
            elif excitation_meshing_type == ExcitationMeshingType.MeshOnlyLowestFrequency:
                self.high_level_geo_gen(frequency=min(frequency_list), skin_mesh_factor=skin_mesh_factor)
            else:
                raise Exception(f"Unknown excitation meshing type {excitation_meshing_type}")
            if self.valid:
                self.mesh.generate_hybrid_mesh(color_scheme, colors_geometry, visualize_before=visualize_before, save_png=save_png)
                self.mesh.generate_electro_magnetic_mesh()

            for i in range(0, len(frequency_list)):
                self.excitation(frequency=frequency_list[i], amplitude_list=current_list_list[i],
                                phase_deg_list=phi_deg_list_list[i])  # frequency and current
                self.file_communication()
                self.pre_simulate()
                self.simulate()
                # self.visualize()

        if self.valid:
            self.calculate_and_write_log(sweep_number=len(frequency_list), currents=current_list_list, frequencies=frequency_list)

            if show_last:
                self.visualize()

            if return_results:
                # Return a Dictionary with the results
                results = {"j2F": self.load_result("j2F", len(frequency_list), "real"),
                           "j2H": self.load_result("j2H", len(frequency_list), "real"),
                           "p_hyst": self.load_result("p_hyst", len(frequency_list), "real")}
                return results

        else:
            if return_results:
                return {"FEM_results": "invalid"}

    def simulate(self):
        """
        Initializes a onelab client. Provides the GetDP based solver with the created mesh file.
        """
        ff.femmt_print(f"\n---\n"
              f"Initialize ONELAB API\n"
              f"Run Simulation\n")

        # -- Simulation --
        # create a new onelab client

        # Initial Clearing of gmsh data
        gmsh.clear()

        # get model file names with correct path
        solver = os.path.join(self.file_data.electro_magnetic_folder_path, "ind_axi_python_controlled.pro")

        os.chdir(self.file_data.working_directory)

        verbose = ""
        if ff.silent:
            verbose = "-verbose 1"
        else:
            verbose = "-verbose 5"

        # Run simulations as sub clients (non blocking??)
        mygetdp = os.path.join(self.file_data.onelab_folder_path, "getdp")
        self.onelab_client.runSubClient("myGetDP", mygetdp + " " + solver + " -msh " + self.file_data.e_m_mesh_file + " -solve Analysis -v2 " + verbose)

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   -  -  -  -  -  -  -  -  -  -  -
    # Miscellaneous
    def calculate_core_volume_with_air(self) -> float:
        """Calculates the volume of the core including air.

        :return: Volume of the core
        :rtype: float
        """
        core_height = self.core.window_h + self.core.core_inner_diameter / 2
        core_width = self.core.r_outer

        return np.pi * core_width**2 * core_height

    def calculate_core_volume(self) -> float:
        """Calculates the volume of the core excluding air.

        :return: Volume of the core.
        :rtype: float
        """
        core_height = self.core.window_h + self.core.core_inner_diameter / 2
        core_width = self.core.r_outer

        winding_height = self.core.window_h
        winding_width = self.core.window_w

        air_gap_volume = 0
        inner_leg_width = self.core.r_inner - winding_width
        for leg_position, position_value, height in self.air_gaps.midpoints:
            width = 0

            if leg_position == AirGapLegPosition.LeftLeg.value:
                # left leg
                # TODO this is wrong since the airgap is not centered on the y axis 
                width = core_width - self.core.r_inner
            elif leg_position == AirGapLegPosition.CenterLeg.value:
                # center leg
                width = inner_leg_width
            elif leg_position == AirGapLegPosition.RightLeg.value:
                # right leg
                # TODO this is wrong since the airgap is not centered on the y axis
                width = core_width - self.core.r_inner
            else:
                raise Exception(f"Unvalid leg position tag {leg_position} used for an air gap.")

            air_gap_volume += np.pi * width**2 * height

        return np.pi*(core_width**2 * core_height - (inner_leg_width+winding_width)**2 * winding_height + inner_leg_width**2 * winding_height) - air_gap_volume

    def get_wire_distances(self) -> List[List[float]]:
        """Helper function which returns the distance (radius) of each conductor to the y-axis 

        :return: Wire distances
        :rtype: List[List[float]]
        """
        wire_distance = []
        for winding in self.two_d_axi.p_conductor:
            # 5 points are for 1 wire
            num_points = len(winding)
            num_windings = num_points//5
            winding_list = []
            for i in range(num_windings):
                winding_list.append(winding[i*5][0])
            wire_distance.append(winding_list)

        return wire_distance

    def calculate_wire_lengths(self) -> List[float]:
        distances = self.get_wire_distances()
        lengths = []
        for winding in distances:
            lengths.append(sum([2 * np.pi * turn for turn in winding]))

        return lengths

    def calculate_wire_volumes(self) -> List[float]:
        wire_volumes = []
        wire_lenghts = self.calculate_wire_lengths()
        for index, winding in enumerate(self.windings):
            cross_section_area = 0
            if winding.conductor_type == ConductorType.RoundLitz or winding.conductor_type == ConductorType.RoundSolid:
                # For round wire its always the same
                cross_section_area = np.pi * winding.conductor_radius ** 2
            elif winding.conductor_type == ConductorType.RectangularSolid:
                # Since the foil sizes also depends on the winding scheme, conductor_arrangement and wrap_para_type
                # the volume calculation is different.
                for vww_index, vww in enumerate(self.virtual_winding_windows):
                    for vww_winding in vww.windings:
                        winding_type = vww_winding.winding_type
                        winding_scheme = vww_winding.winding_scheme
                        if vww_winding.winding_number == index:
                            if winding_type == WindingType.Single:
                                if winding_scheme == WindingScheme.Full:
                                    cross_section_area = self.core.window_h * self.core.window_w
                                elif winding_scheme == WindingScheme.FoilHorizontal:
                                    cross_section_area = self.core.window_w * winding.thickness
                                elif winding_scheme == WindingScheme.FoilVertical:
                                    wrap_para_type = winding.wrap_para
                                    if wrap_para_type == WrapParaType.FixedThickness:
                                        cross_section_area = self.core.window_h * winding.thickness
                                    elif wrap_para_type == WrapParaType.Interpolate:
                                        cross_section_area = self.core.window_h * self.core.window_w / vww.turns[vww_index]
                                    else:
                                        raise Exception(f"Unknown wrap para type {wrap_para_type}")
                                else:
                                    raise Exception(f"Unknown winding scheme {winding_scheme}")
                            elif winding_type == WindingType.Interleaved:
                                # Since interleaved winding type currently only supports round conductors this can be left empty.
                                pass
                            else:
                                raise Exception(f"Unknown winding type {winding_type}")
            else:
                raise Exception(f"Unknown conductor type {winding.conductor_type}")

            wire_volumes.append(cross_section_area * wire_lenghts[index])

        return wire_volumes

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    # GetDP Interaction / Simulation / Excitation
    def excitation(self, frequency: float, amplitude_list: List, phase_deg_list: List = None, ex_type: str = 'current', imposed_red_f=0):
        """
        - excitation of the electromagnetic problem
        - current, voltage or current density
        - frequency or reduced frequency

        :param frequency: Frequency
        :type frequency: float
        :param amplitude_list: Current amplitudes according to windings
        :type amplitude_list: List
        :param phase_deg_list: Current phases in degree according to the current amplitudes (according to windings)
        :type phase_deg_list: List
        :param ex_type: Excitation type. 'Current' implemented only. Future use may: 'voltage' and 'current_density'
        :type ex_type: str
        :param imposed_red_f:
        :type imposed_red_f:
        """
        ff.femmt_print(f"\n---\n"
              f"Excitation: \n"
              f"Frequency: {frequency}\n"
              f"Current(s): {amplitude_list}\n"
              f"Phase(s): {phase_deg_list}\n")


        # -- Excitation --
        self.flag_excitation_type = ex_type  # 'current', 'current_density', 'voltage'

        # Has the user provided a list of phase angles?
        phase_deg_list = phase_deg_list or []
        # phase_deg_list = np.asarray(phase_deg_list)

        for num in range(len(self.windings)):

            # Imposed current
            if self.flag_excitation_type == 'current':
                if len(phase_deg_list) == 0:
                    if self.component_type == ComponentType.Inductor:
                        # Define complex current phasor as real value
                        self.current[num] = complex(amplitude_list[num], 0)
                        phase_deg_list.append(0)  # set to zero

                    else:
                        raise ValueError
                else:
                    self.phase_deg = phase_deg_list
                    # Define complex current phasor as excitation
                    self.current[num] = complex(amplitude_list[num]*np.cos(np.deg2rad(phase_deg_list[num])),
                                                amplitude_list[num]*np.sin(np.deg2rad(phase_deg_list[num])))

        # Imposed current density
        if self.flag_excitation_type == 'current_density':
            raise NotImplementedError

        # Imposed voltage
        if self.flag_excitation_type == 'voltage':
            raise NotImplementedError


        # -- Frequency --

        self.frequency = frequency  # in Hz

        # Define reduced frequency (used for homogenization technique)
        self.red_freq = np.empty(2)

        if self.frequency != 0:
            self.delta = np.sqrt(2 / (2 * self.frequency * np.pi * self.windings[0].cond_sigma * self.mu0)) #TODO: distingish between material conductivities
            for num in range(len(self.windings)):
                if self.windings[num].conductor_type == ConductorType.RoundLitz:
                    self.red_freq[num] = self.windings[num].strand_radius / self.delta
                elif self.windings[num].conductor_type == ConductorType.RoundSolid:
                    self.red_freq[num] = self.windings[num].conductor_radius / self.delta
                else:
                    ff.femmt_print("Reduced Frequency does not have a physical value here")
                    ff.femmt_print(self.windings[num].conductor_type)
                    self.red_freq[num] = 1  # TODO: doesn't make sense like this -> rewrite fore conductor windings shape
        else:
            # DC case
            self.delta = 1e20  # random huge value
            self.red_freq[num] = 0

    def file_communication(self):
        """
        Interaction between python and Prolog files.
        """
        # --------------------------------- File Communication --------------------------------
        # All shared control variables and parameters are passed to a temporary Prolog file
        ff.femmt_print(f"\n---\n"
              f"File Communication\n")

        # Write initialization parameters for simulation in .pro file
        self.write_electro_magnetic_parameter_pro()

        # Write postprocessing parameters in .pro file
        self.write_electro_magnetic_post_pro()


    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    # Post-Processing
    def get_inductances(self, I0: float, op_frequency: float = 0, skin_mesh_factor: float = 1, visualize: bool = False):
        """

        :param visualize:
        :param skin_mesh_factor:
        :param I0:
        :param op_frequency:
        """

        # Remove "old" Inductance Logs
        try:
            os.remove(os.path.join(self.file_data.e_m_values_folder_path, "L_11.dat"))
            os.remove(os.path.join(self.file_data.e_m_values_folder_path, "L_22.dat"))
        except:
            # TODO: Find better way for exception
            pass

        # -- Inductance Estimation --
        self.mesh.mesh(frequency=op_frequency, skin_mesh_factor=skin_mesh_factor)
        # self.high_level_geo_gen(frequency=op_frequency, skin_mesh_factor=skin_mesh_factor)
        # self.mesh.generate_mesh()

        if self.valid:
            frequencies = [op_frequency] * 2
            currents = [[I0, 0], [0, I0]]
            phases = [[0, 180], [0, 180]]

            # self.excitation_sweep(frequencies=frequencies, currents=currents, phi=phases, show_last=visualize, meshing=False)
            self.excitation_sweep_old(frequencies=frequencies, currents=currents, phi=phases, show_last=visualize)

            ff.femmt_print(f"\n"
                  f"                             == Inductances ==                             \n")

            # Read the logged Flux_Linkages
            with open(os.path.join(self.file_data.e_m_values_folder_path, "Flux_Linkage_1.dat")) as fd:
                line = fd.readlines()[-2:]
                # Fluxes induced in Winding 1
                Phi_11 = float(line[0].split(sep=' ')[2])
                Phi_12 = float(line[1].split(sep=' ')[2])

            with open(os.path.join(self.file_data.e_m_values_folder_path, "Flux_Linkage_2.dat")) as fd:
                line = fd.readlines()[-2:]
                # Fluxes induced in Winding 2
                Phi_21 = float(line[0].split(sep=' ')[2])
                Phi_22 = float(line[1].split(sep=' ')[2])

            ff.femmt_print(f"\n"
                  f"Fluxes: \n"
                  f"Phi_11 = {Phi_11}     Induced by I_1 in Winding1 \n"
                  f"Phi_21 = {Phi_21}     Induced by I_1 in Winding2 \n"
                  f"Phi_12 = {Phi_12}     Induced by I_2 in Winding1 \n"
                  f"Phi_22 = {Phi_22}     Induced by I_2 in Winding2 \n")

            """
            # Old way
            # Calculation of inductance matrix
            L_s1 = 0.5*(L_11-sum(self.turns[0])**2/sum(self.turns)[1]**2*L_22+L_k1)
            L_m1 = L_11 - L_s1
            L_m = sum(self.turns[1])/sum(self.turns[0])*L_m1
            L_m2 = sum(self.turns[1])/sum(self.turns[0])*L_m
            L_s2 = L_22 - L_m2
            """

            # Turns Ratio n=N1/N2 with relative winding sense
            if phases[0][0] != phases[0][1]:
                n = -1 * sum(self.windings[0].turns) / sum(self.windings[1].turns)
            else:
                n = sum(self.windings[0].turns) / sum(self.windings[1].turns)

            ff.femmt_print(f"\n"
                  f"Turns Ratio:\n"
                  f"n = {n}\n"
                  )

            # Coupling Factors
            K_21 = Phi_21 / Phi_11
            K_12 = Phi_12 / Phi_22
            k = n / np.abs(n) * (K_21 * K_12) ** 0.5
            ff.femmt_print(f"Coupling Factors:\n"
                  f"K_12 = Phi_12 / Phi_22 = {K_12}\n"
                  f"K_21 = Phi_21 / Phi_11 = {K_21}\n"
                  f"k = Sqrt(K_12 * K_21) = M / Sqrt(L_11 * L_22) = {k}\n"
                  )

            # Read the logged inductance values
            with open(os.path.join(self.file_data.e_m_values_folder_path, "L_11.dat")) as fd:
                line = fd.readlines()[-1]
                words = line.split(sep=' ')
                self.L_11 = float(words[2])
            with open(os.path.join(self.file_data.e_m_values_folder_path, "L_22.dat")) as fd:
                line = fd.readlines()[-1]
                words = line.split(sep=' ')
                self.L_22 = float(words[2])
            ff.femmt_print(f"\n"
                  f"Self Inductances:\n"
                  f"L_11 = {self.L_11}\n"
                  f"L_22 = {self.L_22}\n"
                  )

            # Main/Counter Inductance
            self.M = k * (self.L_11 * self.L_22) ** 0.5
            M_ = self.L_11 * K_21  # Only to proof correctness - ideally: M = M_ = M__
            M__ = self.L_22 * K_12  # Only to proof correctness - ideally: M = M_ = M__
            ff.femmt_print(f"\n"
                  f"Main/Counter Inductance:\n"
                  f"M = k * Sqrt(L_11 * L_22) = {self.M}\n"
                  f"M_ = L_11 * K_21 = {M_}\n"
                  f"M__ = L_22 * K_12 = {M__}\n"
                  )

            # Stray Inductance with 'Turns Ratio' n as 'Transformation Ratio' n
            L_s1 = self.L_11 - self.M * n
            L_s2 = self.L_22 - self.M / n
            L_h = self.M * n
            ff.femmt_print(f"\n"
                  f"T-ECD (primary side transformed):\n"
                  f"[Underdetermined System: 'Transformation Ratio' := 'Turns Ratio']\n"
                  f"    - Transformation Ratio: n\n"
                  f"    - Primary Side Stray Inductance: L_s1\n"
                  f"    - Secondary Side Stray Inductance: L_s2\n"
                  f"    - Primary Side Main Inductance: L_h\n"
                  f"n := n = {n}\n"
                  f"L_s1 = L_11 - M * n = {L_s1}\n"
                  f"L_s2 = L_22 - M / n = {L_s2}\n"
                  f"L_h = M * n = {L_h}\n"
                  )

            # Stray Inductance concentrated on Primary Side
            self.n_conc = self.M / self.L_22
            self.L_s_conc = (1 - k ** 2) * self.L_11
            self.L_h_conc = self.M ** 2 / self.L_22

            ff.femmt_print(f"\n"
                  f"T-ECD (primary side concentrated):\n"
                  f"[Underdetermined System: n := M / L_22  -->  L_s2 = L_22 - M / n = 0]\n"
                  f"    - Transformation Ratio: n\n"
                  f"    - (Primary) Stray Inductance: L_s1\n"
                  f"    - Primary Side Main Inductance: L_h\n"
                  f"n := M / L_22 = k * Sqrt(L_11 / L_22) = {self.n_conc}\n"
                  f"L_s1 = (1 - k^2) * L_11 = {self.L_s_conc}\n"
                  f"L_h = M^2 / L_22 = k^2 * L_11 = {self.L_h_conc}\n"
                  )
            self.visualize()

        else:
            ff.femmt_print(f"Invalid Geometry Data!")

    def get_steinmetz_loss(self, Ipeak: float = None, ki: float = 1, alpha: float = 1.2, beta: float = 2.2, t_rise: float = 3e-6, t_fall: float = 3e-6,
                            f_switch: float = 100000, skin_mesh_factor: float = 0.5):

        """

        :param skin_mesh_factor:
        :param Ipeak:
        :param ki:
        :param alpha:
        :param beta:
        :param t_rise:
        :param t_fall:
        :param f_switch:
        """
        Ipeak = Ipeak or [10, 10]

        self.core.steinmetz_loss = 1
        self.Ipeak = Ipeak
        self.ki = ki
        self.alpha = alpha
        self.beta = beta
        self.t_rise = t_rise
        self.t_fall = t_fall
        self.f_switch = f_switch

        # TODO:
        #  - ki calculation
        #  - all piecewise linear functions (t_1, t_2, ... with for loop in .pro file)

        # Call Simulation
        # self.high_level_geo_gen(frequency=0, skin_mesh_factor=skin_mesh_factor)
        # self.mesh.generate_mesh()
        self.excitation(frequency=f_switch, amplitude_list=Ipeak, phase_deg_list=[0, 180])  # frequency and current
        self.file_com

    def write_electro_magnetic_parameter_pro(self):
        """
        Write materials and other parameters to the "Parameter.pro" file.
        This file is generated by python and is read by gmsh to hand over some parameters.

        Source for the parameters is the MagneticComponent object.
        """
        text_file = open(os.path.join(self.file_data.electro_magnetic_folder_path, "Parameter.pro"), "w")

        # Magnetic Component Type
        if self.component_type == ComponentType.Inductor:
            text_file.write(f"Flag_Transformer = 0;\n")
        if self.component_type == ComponentType.Transformer:
            text_file.write(f"Flag_Transformer = 1;\n")
        if self.component_type == ComponentType.IntegratedTransformer:
            text_file.write(f"Flag_Transformer = 1;\n")

        # Frequency
        text_file.write("Freq = %s;\n" % self.frequency)
        text_file.write(f"delta = {self.delta};\n")

        # Core Loss
        text_file.write(f"Flag_Steinmetz_loss = {self.core.steinmetz_loss};\n")
        text_file.write(f"Flag_Generalized_Steinmetz_loss = {self.core.generalized_steinmetz_loss};\n")

        if self.core.sigma != 0 and self.core.sigma is not None:
            text_file.write(f"Flag_Conducting_Core = 1;\n")
            if isinstance(self.core.sigma, str):
                # TODO: Make following definition general
                # self.core.sigma = 2 * np.pi * self.frequency * self.e0 * f_N95_er_imag(f=self.frequency) + 1 / 6
                self.core.sigma = 1 / 6
            text_file.write(f"sigma_core = {self.core.sigma};\n")
        else:
            text_file.write(f"Flag_Conducting_Core = 0;\n")

        if self.core.steinmetz_loss:
            text_file.write(f"ki = {self.core.ki};\n")
            text_file.write(f"alpha = {self.core.alpha};\n")
            text_file.write(f"beta = {self.core.beta};\n")
        if self.core.generalized_steinmetz_loss:
            text_file.write(f"t_rise = {self.t_rise};\n")
            text_file.write(f"t_fall = {self.t_fall};\n")
            text_file.write(f"f_switch = {self.f_switch};\n")

        # Conductor specific definitions
        for num in range(len(self.windings)):
            # -- Control Flags --
            if self.flag_excitation_type == 'current':
                text_file.write(f"Flag_ImposedVoltage = 0;\n")
            if self.flag_excitation_type == 'voltage':
                text_file.write(f"Flag_ImposedVoltage = 1;\n")
            if self.windings[num].conductor_type == ConductorType.RoundLitz:
                text_file.write(f"Flag_HomogenisedModel{num + 1} = 1;\n")
            else:
                text_file.write(f"Flag_HomogenisedModel{num + 1} = 0;\n")

            # -- Geometry --
            # Number of turns per conductor
            turns = 0
            for vww in self.virtual_winding_windows:
                for index, winding in enumerate(vww.windings):
                    if winding.winding_number == num:
                        turns += vww.turns[index]
            text_file.write(f"NbrCond{num + 1} = {turns};\n")

            # For stranded Conductors:
            # text_file.write(f"NbrstrandedCond = {self.turns};\n")  # redundant
            if self.windings[num].conductor_type == ConductorType.RoundLitz:
                text_file.write(f"NbrStrands{num + 1} = {self.windings[num].n_strands};\n")
                text_file.write(f"Fill{num + 1} = {self.windings[num].ff};\n")
                # ---
                # Litz Approximation Coefficients were created with 4 layers
                # That's why here a hard-coded 4 is implemented
                # text_file.write(f"NbrLayers{num+1} = {self.n_layers[num]};\n")
                text_file.write(f"NbrLayers{num + 1} = 4;\n")
            text_file.write(f"AreaCell{num + 1} = {self.windings[num].a_cell};\n")
            text_file.write(f"Rc{num + 1} = {self.windings[num].conductor_radius};\n")

            # -- Excitation --
            # Imposed current, current density or voltage
            if self.flag_excitation_type == 'current':
                text_file.write(f"Val_EE_{num + 1} = {abs(self.current[num])};\n")
                text_file.write(f"Phase_{num + 1} = {np.deg2rad(self.phase_deg[num])};\n")
                text_file.write(f"Parallel_{num + 1} = {self.windings[num].parallel};\n")

            if self.flag_excitation_type == 'current_density':
                text_file.write(f"Val_EE_{num + 1} = {self.current_density[num]};\n")
                raise NotImplementedError

            if self.flag_excitation_type == 'voltage':
                text_file.write(f"Val_EE_{num + 1} = {self.voltage[num]};\n")
                raise NotImplementedError

            ff.femmt_print(f"Cell surface area: {self.windings[num].a_cell} \n"
                  f"Reduced frequency: {self.red_freq[num]}")

            if self.red_freq[num] > 1.25 and self.windings[num].conductor_type == ConductorType.RoundLitz:
                # TODO: Allow higher reduced frequencies
                ff.femmt_print(f"Litz Coefficients only implemented for X<=1.25")
                raise Warning
            # Reduced Frequency
            text_file.write(f"Rr{num + 1} = {self.red_freq[num]};\n")

            # Material Properties
            # Conductor Material
            text_file.write(f"sigma_winding_{num + 1} = {self.windings[num].cond_sigma};\n")

        # -- Materials --

        # Nature Constants
        text_file.write(f"mu0 = 4.e-7 * Pi;\n")
        text_file.write(f"nu0 = 1 / mu0;\n")
        text_file.write(f"e0 = {self.e0};\n")

        # Material Properties

        # Core Material
        # if self.frequency == 0:
        if self.core.non_linear:
            text_file.write(f"Flag_NL = 1;\n")
            text_file.write(f"Core_Material = {self.core.material};\n")  # relative permeability is defined at simulation runtime            text_file.write(f"Flag_NL = 0;\n")
        else:
            text_file.write(f"Flag_NL = 0;\n")
            text_file.write(f"mur = {self.core.mu_rel};\n")  # mur is predefined to a fixed value
            if self.core.permeability_type == PermeabilityType.FromData:
                text_file.write(f"Flag_Permeability_From_Data = 1;\n")  # mur is predefined to a fixed value
            else:
                text_file.write(f"Flag_Permeability_From_Data = 0;\n")  # mur is predefined to a fixed value
            if self.core.permeability_type == PermeabilityType.FixedLossAngle:
                text_file.write(f"phi_mu_deg = {self.core.phi_mu_deg};\n")  # loss angle for complex representation of hysteresis loss
                text_file.write(f"mur_real = {self.core.mu_rel * np.cos(np.deg2rad(self.core.phi_mu_deg))};\n")  # Real part of complex permeability
                text_file.write(f"mur_imag = {self.core.mu_rel * np.sin(np.deg2rad(self.core.phi_mu_deg))};\n")  # Imaginary part of complex permeability
                text_file.write(f"Flag_Fixed_Loss_Angle = 1;\n")  # loss angle for complex representation of hysteresis loss
            else:
                text_file.write(f"Flag_Fixed_Loss_Angle = 0;\n")  # loss angle for complex representation of hysteresis loss

        # if self.frequency != 0:
        #    text_file.write(f"Flag_NL = 0;\n")
        #    text_file.write(f"mur = {self.core.mu_rel};\n")

        text_file.close()

    def write_electro_magnetic_post_pro(self):
        """
        """
        text_file = open(os.path.join(self.file_data.electro_magnetic_folder_path, "postquantities.pro"), "w") 

        # This is needed because the f string cant contain a \ in {}
        backslash = "\\"

        text_file.write(f"DirRes = \"{self.file_data.results_folder_path.replace(backslash, '/')}/\";\n")
        text_file.write(f"DirResFields = \"{self.file_data.e_m_fields_folder_path.replace(backslash, '/')}/\";\n")
        text_file.write(f"DirResVals = \"{self.file_data.e_m_values_folder_path.replace(backslash, '/')}/\";\n")
        text_file.write(f"DirResValsPrimary = \"{self.file_data.e_m_values_folder_path.replace(backslash, '/')}/Primary/\";\n")
        text_file.write(f"DirResValsSecondary = \"{self.file_data.e_m_values_folder_path.replace(backslash, '/')}/Secondary/\";\n")
        text_file.write(f"DirResCirc = \"{self.file_data.e_m_circuit_folder_path.replace(backslash, '/')}/\";\n")
        text_file.write(f"OptionPos = \"{self.file_data.results_folder_path.replace(backslash, '/')}/option.pos\";\n")

        # Visualisation
        if self.plot_fields == "standard":
            text_file.write(f"Flag_show_standard_fields = 1;\n")
        else:
            text_file.write(f"Flag_show_standard_fields = 0;\n")

        text_file.close()

    def calculate_and_write_log(self, sweep_number: int = 1, currents: List = None, frequencies: List = None):
        """
        Method reads back the results from the .dat result files created by the ONELAB simulation client and stores
        them in a dictionary. Additionally, the input settings which are used in order to create the simulation are also printed.
        From this data type a JSON log file is created.
        :param sweep_number: Number of sweep iterations that were done before. For a single simulation sweep_number = 1
        :param currents: Current values of the sweep iterations. Not needed for single simulation
        :param frequencies: frequencies values of the sweep iterations. Not needed for single simulation
        """
        # ---- Print simulation results ----

        fundamental_index = 0  # index of the fundamental frequency

        # create the dictionary used to log the result data
        #   - 'single_sweeps' contains a list of n=sweep_number sweep_dicts
        #   - 'total_losses' contains the sums of all sweeps
        #       - in case complex core parameters are used:
        #           - hysteresis losses are taken only from fundamental frequency
        #               - fundamental frequency is the smallest frequency != 0
        #   - 'simulation_settings' contains the simulation configuration (for more info see encode_settings())

        log_dict = {"single_sweeps": [], "total_losses": {}, "simulation_settings": {}}

        for sweep_run in range(0, sweep_number):
            # create dictionary sweep_dict with 'Winding' as a list of m=n_windings winding_dicts.
            # Single values, received during one of the n=sweep_number simulations, are added as 'core_eddy_losses' etc.
            sweep_dict = {}

            # Frequencies
            if sweep_number > 1:
                # sweep_simulation -> get frequencies from passed currents
                sweep_dict["f"] = frequencies[sweep_run]
                fundamental_index = np.argmin(np.ma.masked_where(np.array(frequencies) == 0, np.array(frequencies)))
            else:
                # single_simulation -> get frequency from instance variable
                sweep_dict["f"] = self.frequency

            # Winding names are needed to find the logging path
            winding_name = ["Primary", "Secondary", "Tertiary"]

            for winding in range(len(self.windings)):

                # Create empty winding dictionary
                # create dictionary winding_dict with 'turn_losses' as list of the j=number_turns turn losses.
                # Single values related to one winding are added as 'winding_losses' etc.
                winding_dict = {"turn_losses": [],
                                "flux": [],
                                "self_inductance": [],
                                "mag_field_energy": [],
                                "V": []}

                # Number turns
                turns = 0
                for vww in self.virtual_winding_windows:
                    for index, conductor in enumerate(vww.windings):
                        if conductor.winding_number == winding:
                            turns += vww.turns[index]
                winding_dict["number_turns"] = turns

                # Currents
                if sweep_number > 1:
                    # sweep_simulation -> get currents from passed currents
                    complex_current_phasor = currents[sweep_run][winding]
                else:
                    # single_simulation -> get current from instance variable
                    complex_current_phasor = self.current[winding]

                # Store complex value as list in json (because json isnt natively capable of complex values)
                winding_dict["I"] = [complex_current_phasor.real, complex_current_phasor.imag]


                # Case litz: Load homogenized results
                if self.windings[winding].conductor_type == ConductorType.RoundLitz:
                    winding_dict["winding_losses"] = self.load_result(res_name=f"j2H_{winding + 1}", last_n=sweep_number)[sweep_run]
                    for turn in range(0, winding_dict["number_turns"]):
                        winding_dict["turn_losses"].append(self.load_result(res_name=winding_name[winding] + f"/Losses_turn_{turn + 1}", last_n=sweep_number)[sweep_run])

                # Case litz: Load homogenized results
                else:
                    winding_dict["winding_losses"] = self.load_result(res_name=f"j2F_{winding + 1}", last_n=sweep_number)[sweep_run]
                    for turn in range(0, winding_dict["number_turns"]):
                        winding_dict["turn_losses"].append(self.load_result(res_name=winding_name[winding] + f"/Losses_turn_{turn + 1}", last_n=sweep_number)[sweep_run])

                # Flux
                winding_dict["flux"].append(self.load_result(res_name=f"Flux_Linkage_{winding + 1}", last_n=sweep_number)[sweep_run])
                winding_dict["flux"].append(self.load_result(res_name=f"Flux_Linkage_{winding + 1}", part="imaginary", last_n=sweep_number)[sweep_run])

                # Inductance
                winding_dict["self_inductance"].append(self.load_result(res_name=f"L_{winding + 1}{winding + 1}", part="real", last_n=sweep_number)[sweep_run])
                winding_dict["self_inductance"].append(self.load_result(res_name=f"L_{winding + 1}{winding + 1}", part="imaginary", last_n=sweep_number)[sweep_run])

                # Magnetic Field Energy
                winding_dict["mag_field_energy"].append(self.load_result(res_name=f"ME", last_n=sweep_number)[sweep_run])
                winding_dict["mag_field_energy"].append(self.load_result(res_name=f"ME", part="imaginary", last_n=sweep_number)[sweep_run])

                # Voltage
                winding_dict["V"].append(self.load_result(res_name=f"Voltage_{winding + 1}", part="real", last_n=sweep_number)[sweep_run])
                winding_dict["V"].append(self.load_result(res_name=f"Voltage_{winding + 1}", part="imaginary", last_n=sweep_number)[sweep_run])
                complex_voltage_phasor = complex(winding_dict["V"][0], winding_dict["V"][1])

                # Power
                # using 'winding_dict["V"][0]' to get first element (real part) of V. Use winding_dict["I"][0] to avoid typeerror
                winding_dict["P"] = (complex_voltage_phasor * complex_current_phasor.conjugate() / 2).real
                winding_dict["Q"] = (complex_voltage_phasor * complex_current_phasor.conjugate() / 2).imag
                winding_dict["S"] = np.sqrt(winding_dict["P"] ** 2 + winding_dict["Q"] ** 2)


                sweep_dict[f"winding{winding+1}"] = winding_dict


            # Core losses TODO: Choose between Steinmetz or complex core losses
            sweep_dict["core_eddy_losses"] = self.load_result(res_name="CoreEddyCurrentLosses", last_n=sweep_number)[sweep_run]
            sweep_dict["core_hyst_losses"] = self.load_result(res_name="p_hyst", last_n=sweep_number)[sweep_run]

            # Sum losses of all windings of one single run
            sweep_dict["all_winding_losses"] = sum(sweep_dict[f"winding{d+1}"]["winding_losses"] for d in range(len(self.windings)))

            log_dict["single_sweeps"].append(sweep_dict)


        # Total losses of excitation sweep
        # Sum losses of all sweep runs. For core losses just use hyst_losses of the fundamental frequency.
        # Also needed as excitation for steady state thermal simulations

        # Single Windings
        for winding in range(len(self.windings)):
            # Number of turns per conductor
            turns = 0
            for vww in self.virtual_winding_windows:
                for index, conductor in enumerate(vww.windings):
                    if conductor.winding_number == winding:
                        turns += vww.turns[index]

            log_dict["total_losses"][f"winding{winding + 1}"] = {
                "total": sum(sum(log_dict["single_sweeps"][d][f"winding{winding+1}"]["turn_losses"]) for d in range(len(log_dict["single_sweeps"]))),
                "turns": []
            }
            for turn in range(0, turns):
                log_dict["total_losses"][f"winding{winding + 1}"]["turns"].append(sum(log_dict["single_sweeps"][d][f"winding{winding+1}"]["turn_losses"][turn] for d in range(len(log_dict["single_sweeps"]))))

        # Winding (all windings)
        log_dict["total_losses"]["all_windings"] = sum(log_dict["single_sweeps"][d]["all_winding_losses"] for d in range(len(log_dict["single_sweeps"])))

        # Core
        log_dict["total_losses"]["eddy_core"] = sum(log_dict["single_sweeps"][d]["core_eddy_losses"] for d in range(len(log_dict["single_sweeps"])))
        # For core losses just use hyst_losses of the fundamental frequency. When using single_simulation, the fundamental frquency is at [0]
        # => just an approximation for excitation sweeps!
        log_dict["total_losses"]["hyst_core_fundamental_freq"] = log_dict["single_sweeps"][fundamental_index]["core_hyst_losses"]

        # Total losses of inductive component according to single or sweep simulation
        log_dict["total_losses"]["core"] = log_dict["total_losses"]["hyst_core_fundamental_freq"] + log_dict["total_losses"]["eddy_core"]

        # ---- Miscellaneous ----
        log_dict["misc"] = {
            "core_2daxi_volume": self.calculate_core_volume(),
            "core_2daxi_total_volume":self.calculate_core_volume_with_air(),
            "core_2daxi_weight": -1,
            "wire_lengths": self.calculate_wire_lengths(),
            "wire_volumes": self.calculate_wire_volumes(),
            "wire_weight": -1
        }

        # ---- Print current configuration ----
        log_dict["simulation_settings"] = MagneticComponent.encode_settings(self)

        # ====== save data as JSON ======
        with open(self.file_data.e_m_results_log_path, "w+", encoding='utf-8') as outfile:
            json.dump(log_dict, outfile, indent=2, ensure_ascii=False)

    def read_log(self):
        """
        """
        log = {}
        with open(self.file_data.e_m_results_log_path, "r") as fd:
            content = json.loads(fd.read())
            log = content

        return log

    def visualize(self):
        """
        - a post simulation viewer
        - allows to open ".pos"-files in gmsh
        - For example current density, ohmic losses or the magnetic field density can be visualized
        """
        # ---------------------------------------- Visualization in gmsh ---------------------------------------
        ff.femmt_print(f"\n---\n"
              f"Visualize fields in GMSH front end:\n")

        # gmsh.initialize()
        epsilon = 1e-9

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Colors
        # gmsh.option.setColor([(2, 1)], 50, 50, 50)
        # gmsh.option.setColor([(2, 2)], 0, 0, 230)
        # gmsh.option.setColor([(2, 3)], 150, 150, 0)
        # gmsh.option.setColor([(2, 4)], 150, 150, 0)
        # gmsh.option.setColor([(2, 5)], 150, 150, 0)
        # gmsh.option.setColor([(2, 6)], 150, 150, 0)
        # gmsh.option.setColor([(2, 7)], 150, 150, 0)
        # gmsh.option.setNumber("Mesh.SurfaceFaces", 1)
        # gmsh.option.setNumber("Mesh.ColorCarousel", 1)
        # gmsh.option.setNumber("Mesh.LabelType", 1)

        # Mesh
        gmsh.option.setNumber("Mesh.SurfaceEdges", 0)
        view = 0

        if any(self.windings[i].conductor_type != ConductorType.RoundLitz for i in range(len(self.windings))):
            # Ohmic losses (weighted effective value of current density)
            gmsh.open(os.path.join(self.file_data.e_m_fields_folder_path, "j2F.pos"))
            gmsh.option.setNumber(f"View[{view}].ScaleType", 2)
            gmsh.option.setNumber(f"View[{view}].RangeType", 2)
            gmsh.option.setNumber(f"View[{view}].SaturateValues", 1)
            gmsh.option.setNumber(f"View[{view}].CustomMin", gmsh.option.getNumber(f"View[{view}].Min") + epsilon)
            gmsh.option.setNumber(f"View[{view}].CustomMax", gmsh.option.getNumber(f"View[{view}].Max"))
            gmsh.option.setNumber(f"View[{view}].ColormapNumber", 1)
            gmsh.option.setNumber(f"View[{view}].IntervalsType", 2)
            gmsh.option.setNumber(f"View[{view}].NbIso", 40)
            gmsh.option.setNumber(f"View[{view}].ShowTime", 0)
            ff.femmt_print(gmsh.option.getNumber(f"View[{view}].Max"))
            view += 1

        if any(self.windings[i].conductor_type == ConductorType.RoundLitz for i in range(len(self.windings))):
            # Ohmic losses (weighted effective value of current density)
            gmsh.open(os.path.join(self.file_data.e_m_fields_folder_path, "jH.pos"))
            gmsh.option.setNumber(f"View[{view}].ScaleType", 2)
            gmsh.option.setNumber(f"View[{view}].RangeType", 2)
            gmsh.option.setNumber(f"View[{view}].SaturateValues", 1)
            gmsh.option.setNumber(f"View[{view}].CustomMin", gmsh.option.getNumber(f"View[{view}].Min") + epsilon)
            gmsh.option.setNumber(f"View[{view}].CustomMax", gmsh.option.getNumber(f"View[{view}].Max"))
            gmsh.option.setNumber(f"View[{view}].ColormapNumber", 1)
            gmsh.option.setNumber(f"View[{view}].IntervalsType", 2)
            gmsh.option.setNumber(f"View[{view}].ShowTime", 0)
            gmsh.option.setNumber(f"View[{view}].NbIso", 40)
            view += 1

        # Magnetic flux density
        gmsh.open(os.path.join(self.file_data.e_m_fields_folder_path, "Magb.pos"))
        gmsh.option.setNumber(f"View[{view}].ScaleType", 1)
        gmsh.option.setNumber(f"View[{view}].RangeType", 1)
        gmsh.option.setNumber(f"View[{view}].CustomMin", gmsh.option.getNumber(f"View[{view}].Min") + epsilon)
        gmsh.option.setNumber(f"View[{view}].CustomMax", gmsh.option.getNumber(f"View[{view}].Max"))
        gmsh.option.setNumber(f"View[{view}].ColormapNumber", 1)
        gmsh.option.setNumber(f"View[{view}].IntervalsType", 2)
        gmsh.option.setNumber(f"View[{view}].ShowTime", 0)
        gmsh.option.setNumber(f"View[{view}].NbIso", 40)
        view += 1

        """
        # Vector Potential
        gmsh.open(os.path.join(self.file_data.e_m_fields_folder_path, "raz.pos"))
        gmsh.option.setNumber(f"View[{view}].ScaleType", 1)
        gmsh.option.setNumber(f"View[{view}].RangeType", 1)
        gmsh.option.setNumber(f"View[{view}].CustomMin", gmsh.option.getNumber(f"View[{view}].Min") + epsilon)
        gmsh.option.setNumber(f"View[{view}].CustomMax", gmsh.option.getNumber(f"View[{view}].Max"))
        gmsh.option.setNumber(f"View[{view}].ColormapNumber", 1)
        gmsh.option.setNumber(f"View[{view}].IntervalsType", 2)
        gmsh.option.setNumber(f"View[{view}].NbIso", 40)
        view += 1
        """

        gmsh.fltk.run()
        # gmsh.finalize()

    def get_loss_data(self, last_n_values: int, loss_type: str = 'litz_loss'):
        """
        Returns the last n values from the chosen loss type logged in the result folder.

        :param last_n_values: _description_
        :type last_n_values: int
        :param loss_type: _description_, defaults to 'litz_loss'
        :type loss_type: str, optional
        """
        loss_file = None
        # Loss file location
        if loss_type == 'litz_loss':
            loss_file = 'j2H.dat'
        if loss_type == 'solid_loss':
            loss_file = 'j2F.dat'
        # Read the logged losses corresponding to the frequencies
        with open(os.path.join(self.file_data.e_m_values_folder_path, loss_file), newline='') as fd:
            reader = csv.reader(fd)
            data = list(reader)
        return data[-last_n_values:-1] + [data[-1]]

    def load_result(self, res_name: str, res_type: str = "value", last_n: int = 1, part:str = "real", position: int = 0):
        """
        Loads the "last_n" parameters from a result file of the scalar quantity "res_name".
        Either the real or imaginary part can be chosen.
        :param part: "real" or "imaginary" part can be chosen
        :param res_name: name of the quantity
        :param res_type: type of the quantity: "value" or "circuit"
        :param last_n: Number of parameters to be loaded
        :return: last_n entries of the chosen result file
        :rtype: list
        """
        if res_type=="value":
            res_path=self.file_data.e_m_values_folder_path
        if res_type=="circuit":
            res_path=self.file_data.e_m_circuit_folder_path


        with open(os.path.join(res_path, f"{res_name}.dat")) as fd:
            lines = fd.readlines()[-last_n:]

            if part == "real":
                result = [float(line.split(sep=' ')[1 + 2*position + 1]) for n, line in enumerate(lines)]
            if part == "imaginary":
                result = [float(line.split(sep=' ')[2 + 2*position + 1]) for n, line in enumerate(lines)]

            return result

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    # Litz Approximation [internal methods]
    def pre_simulate(self):
        """
        Used to determine the litz-approximation coefficients.
        """
        for num in range(len(self.windings)):
            if self.windings[num].conductor_type == ConductorType.RoundLitz:
                # ---
                # Litz Approximation Coefficients were created with 4 layers
                # That's why here a hard-coded 4 is implemented
                # if os.path.isfile(self.path +
                # f"/Strands_Coefficients/coeff/pB_RS_la{self.windings[num].ff}_{self.n_layers[num]}layer.dat"):
                if os.path.exists(os.path.join(self.file_data.e_m_strands_coefficients_folder_path, "coeff", f"pB_RS_la{self.windings[num].ff}_4layer.dat")):
                    ff.femmt_print("Coefficients for stands approximation are found.")

                else:
                    # Rounding X to fit it with corresponding parameters from the database
                    X = self.red_freq[num]
                    X = np.around(X, decimals=3)
                    ff.femmt_print(f"Rounded Reduced frequency X = {X}")
                    self.create_strand_coeff(num)

    def create_strand_coeff(self, num: int):
        """TODO Doc
        """
        ff.femmt_print(f"\n"
              f"Pre-Simulation\n"
              f"-----------------------------------------\n"
              f"Create coefficients for strands approximation\n")
        # Create a new onelab client
        # -- Pre-Simulation Settings --
        text_file = open(os.path.join(self.file_data.e_m_strands_coefficients_folder_path, "PreParameter.pro"), "w")
        # ---
        # Litz Approximation Coefficients are created with 4 layers
        # That's why here a hard-coded 4 is implemented
        # text_file.write(f"NbrLayers = {self.n_layers[num]};\n")
        text_file.write(f"NbrLayers = 4;\n")
        text_file.write(f"Fill = {self.windings[num].ff};\n")
        ff.femmt_print("Here")
        text_file.write(f"Rc = {self.windings[num].strand_radius};\n")  # double named!!! must be changed
        text_file.close()
        cell_geo = os.path.join(self.file_data.e_m_strands_coefficients_folder_path, "cell.geo")

        # Run gmsh as a sub client
        mygmsh = os.path.join(self.file_data.onelab_folder_path, "gmsh")
        self.onelab_client.runSubClient("myGmsh", mygmsh + " " + cell_geo + " -2 -v 2")

        modes = [1, 2]  # 1 = "skin", 2 = "proximity"
        reduced_frequencies = np.linspace(0, 1.25, 6)  # must be even
        for mode in modes:
            for rf in reduced_frequencies:
                # -- Pre-Simulation Settings --
                text_file = open(os.path.join(self.file_data.e_m_strands_coefficients_folder_path, "PreParameter.pro"), "w")
                text_file.write(f"Rr_cell = {rf};\n")
                text_file.write(f"Mode = {mode};\n")
                # Litz Approximation Coefficients are created with 4 layers
                # That's why here a hard-coded 4 is implemented
                # text_file.write(f"NbrLayers = {self.n_layers[num]};\n")
                text_file.write(f"NbrLayers = 4;\n")

                text_file.write(f"Fill = {self.windings[num].ff};\n")
                text_file.write(f"Rc = {self.windings[num].strand_radius};\n")  # double named!!! must be changed
                text_file.close()

                # get model file names with correct path
                input_file = os.path.join(self.file_data.e_m_strands_coefficients_folder_path, "cell_dat.pro")
                cell = os.path.join(self.file_data.e_m_strands_coefficients_folder_path, "cell.pro")

                # Run simulations as sub clients
                mygetdp = os.path.join(self.file_data.onelab_folder_path, "getdp")
                self.onelab_client.runSubClient("myGetDP", mygetdp + " " + cell + " -input " + input_file + " -solve MagDyn_a -v2")

        # Formatting stuff
        # Litz Approximation Coefficients are created with 4 layers
        # That's why here a hard-coded 4 is implemented
        # files = [self.path + f"/Strands_Coefficients/coeff/pB_RS_la{self.windings[num].ff}_{self.n_layers[num]}layer.dat",
        #         self.path + f"/Strands_Coefficients/coeff/pI_RS_la{self.windings[num].ff}_{self.n_layers[num]}layer.dat",
        #         self.path + f"/Strands_Coefficients/coeff/qB_RS_la{self.windings[num].ff}_{self.n_layers[num]}layer.dat",
        #         self.path + f"/Strands_Coefficients/coeff/qI_RS_la{self.windings[num].ff}_{self.n_layers[num]}layer.dat"]
        
        coeff_folder = os.path.join(self.file_data.e_m_strands_coefficients_folder_path, "coeff") 
        if not os.path.isdir(coeff_folder):
            os.mkdir(coeff_folder)

        files = [os.path.join(coeff_folder, f"pB_RS_la{self.windings[num].ff}_4layer.dat"),
                 os.path.join(coeff_folder, f"pI_RS_la{self.windings[num].ff}_4layer.dat"),
                 os.path.join(coeff_folder, f"qB_RS_la{self.windings[num].ff}_4layer.dat"),
                 os.path.join(coeff_folder, f"qI_RS_la{self.windings[num].ff}_4layer.dat")]
        for i in range(0, 4):
            with fileinput.FileInput(files[i], inplace=True) as file:
                for line in file:
                    ff.femmt_print(line.replace(' 0\n', '\n'), end='')

        # Corrects pB coefficient error at 0Hz
        # Must be changed in future in cell.pro
        for i in range(0, 4):
            with fileinput.FileInput(files[i], inplace=True) as file:
                for line in file:
                    ff.femmt_print(line.replace(' 0\n', ' 1\n'), end='')

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    # FEMM [alternative Solver]
    def femm_reference(self, freq: float, current: float, sign: bool = None, non_visualize: int = 0):
        """
        Allows reference simulations with the 2D open source electromagnetic FEM tool FEMM.
        Helpful to validate changes (especially in the Prolog Code).


        Blockprop <--> Group Convention:
                                            Ferrite := 0
                                            Air := 1
                                            Winding 1 := 2
                                            Winding 2 := 3
                                            ...
                                            Winding n := n+1

        :param sign:
        :param non_visualize:
        :param freq:
        :param current:
        """
        if os.name == 'nt':
            ff.install_pyfemm_if_missing()
            if not self.femm_is_imported:
                globals()["femm"] = __import__("femm")
                self.femm_is_imported = True
        else:
            raise Exception('You are using a computer that is not running windows. '
                            'This command is only executable on Windows computers.')


        self.create_folders(self.file_data.femm_folder_path)

        sign = sign or [1]

        if sign is None:
            sign = [1, 1]

        # == Pre Geometry ==
        self.high_level_geo_gen()

        # == Init ==
        femm.openfemm(non_visualize)
        femm.newdocument(0)
        femm.mi_probdef(freq, 'meters', 'axi', 1.e-8, 0, 30)

        # == Materials ==
        if self.core.sigma != 0:
            if isinstance(self.core.sigma, str):
                # TODO: Make following definition general
                # self.core.sigma = 2 * np.pi * self.frequency * self.e0 * f_N95_er_imag(f=self.frequency) + 1 / 6
                self.core.sigma = 1 / 6


        ff.femmt_print(f"{self.core.permeability_type=}, {self.core.sigma=}")
        if self.core.permeability_type == PermeabilityType.FixedLossAngle:
            femm.mi_addmaterial('Ferrite', self.core.mu_rel, self.core.mu_rel, 0, 0, self.core.sigma/1e6, 0, 0, 1, 0, self.core.phi_mu_deg, self.core.phi_mu_deg)
        elif self.core.permeability_type == PermeabilityType.RealValue:
            femm.mi_addmaterial('Ferrite', self.core.mu_rel, self.core.mu_rel, 0, 0, self.core.sigma/1e6, 0, 0, 1, 0, self.core.phi_mu_deg, self.core.phi_mu_deg)
        else:
            femm.mi_addmaterial('Ferrite', self.core.mu_rel, self.core.mu_rel, 0, 0, self.core.sigma/1e6, 0, 0, 1, 0, 0, 0)
        femm.mi_addmaterial('Air', 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0)
        if self.windings[0].conductor_type == ConductorType.RoundLitz:
            femm.mi_addmaterial('Copper', 1, 1, 0, 0, self.windings[0].cond_sigma/1e6, 0, 0, 1, 5, 0, 0, self.windings[0].n_strands,
                                2 * 1000 * self.windings[0].strand_radius)  # type := 5. last argument
            ff.femmt_print(f"Number of strands: {self.windings[0].n_strands}")
            ff.femmt_print(f"Diameter of strands in mm: {2 * 1000 * self.windings[0].strand_radius}")
        if self.windings[0].conductor_type == ConductorType.RoundSolid:
            femm.mi_addmaterial('Copper', 1, 1, 0, 0, self.windings[0].cond_sigma/1e6, 0, 0, 1, 0, 0, 0, 0, 0)

        # == Circuit ==
        # coil as seen from the terminals.
        femm.mi_addcircprop('Primary', current[0] * sign[0], 1)
        if self.component_type == (ComponentType.Transformer or ComponentType.IntegratedTransformer):
            femm.mi_addcircprop('Secondary', current[1] * sign[1], 1)

        # == Geometry ==
        # Add core
        if self.air_gaps.number == 0:
            femm.mi_drawline(self.two_d_axi.p_window[4, 0],
                            self.two_d_axi.p_window[4, 1],
                            self.two_d_axi.p_window[5, 0],
                            self.two_d_axi.p_window[5, 1])
            femm.mi_drawline(self.two_d_axi.p_window[5, 0],
                            self.two_d_axi.p_window[5, 1],
                            self.two_d_axi.p_window[7, 0],
                            self.two_d_axi.p_window[7, 1])
            femm.mi_drawline(self.two_d_axi.p_window[7, 0],
                            self.two_d_axi.p_window[7, 1],
                            self.two_d_axi.p_window[6, 0],
                            self.two_d_axi.p_window[6, 1])
            femm.mi_drawline(self.two_d_axi.p_window[6, 0],
                            self.two_d_axi.p_window[6, 1],
                            self.two_d_axi.p_window[4, 0],
                            self.two_d_axi.p_window[4, 1])
            femm.mi_drawline(0,
                            self.two_d_axi.p_outer[2, 1],
                            self.two_d_axi.p_outer[3, 0],
                            self.two_d_axi.p_outer[3, 1])
            femm.mi_drawline(self.two_d_axi.p_outer[3, 0],
                            self.two_d_axi.p_outer[3, 1],
                            self.two_d_axi.p_outer[1, 0],
                            self.two_d_axi.p_outer[1, 1])
            femm.mi_drawline(self.two_d_axi.p_outer[1, 0],
                            self.two_d_axi.p_outer[1, 1],
                            0,
                            self.two_d_axi.p_outer[0, 1])
            femm.mi_drawline(0,
                            self.two_d_axi.p_outer[0, 1],
                            0,
                            self.two_d_axi.p_outer[2, 1])
        elif self.air_gaps.number > 0:
            femm.mi_drawline(0,
                            self.two_d_axi.p_air_gaps[0, 1],
                            self.two_d_axi.p_air_gaps[1, 0],
                            self.two_d_axi.p_air_gaps[1, 1])
            femm.mi_drawline(self.two_d_axi.p_air_gaps[1, 0],
                            self.two_d_axi.p_air_gaps[1, 1],
                            self.two_d_axi.p_window[4, 0],
                            self.two_d_axi.p_window[4, 1])
            femm.mi_drawline(self.two_d_axi.p_window[4, 0],
                            self.two_d_axi.p_window[4, 1],
                            self.two_d_axi.p_window[5, 0],
                            self.two_d_axi.p_window[5, 1])
            femm.mi_drawline(self.two_d_axi.p_window[5, 0],
                            self.two_d_axi.p_window[5, 1],
                            self.two_d_axi.p_window[7, 0],
                            self.two_d_axi.p_window[7, 1])
            femm.mi_drawline(self.two_d_axi.p_window[7, 0],
                            self.two_d_axi.p_window[7, 1],
                            self.two_d_axi.p_window[6, 0],
                            self.two_d_axi.p_window[6, 1])
            femm.mi_drawline(self.two_d_axi.p_window[6, 0],
                            self.two_d_axi.p_window[6, 1],
                            self.two_d_axi.p_air_gaps[3, 0],
                            self.two_d_axi.p_air_gaps[3, 1])
            femm.mi_drawline(self.two_d_axi.p_air_gaps[3, 0],
                            self.two_d_axi.p_air_gaps[3, 1],
                            0,
                            self.two_d_axi.p_air_gaps[2, 1])
            femm.mi_drawline(0,
                            self.two_d_axi.p_air_gaps[2, 1],
                            0,
                            self.two_d_axi.p_outer[2, 1])
            femm.mi_drawline(0,
                            self.two_d_axi.p_outer[2, 1],
                            self.two_d_axi.p_outer[3, 0],
                            self.two_d_axi.p_outer[3, 1])
            femm.mi_drawline(self.two_d_axi.p_outer[3, 0],
                            self.two_d_axi.p_outer[3, 1],
                            self.two_d_axi.p_outer[1, 0],
                            self.two_d_axi.p_outer[1, 1])
            femm.mi_drawline(self.two_d_axi.p_outer[1, 0],
                            self.two_d_axi.p_outer[1, 1],
                            0,
                            self.two_d_axi.p_outer[0, 1])
            femm.mi_drawline(0,
                            self.two_d_axi.p_outer[0, 1],
                            0,
                            self.two_d_axi.p_air_gaps[0, 1])
        else:
            raise Exception("Negative air gap number is not allowed")
        # Add Coil

        # femm.mi_drawrectangle(self.two_d_axi.p_window[4, 0]+self.insulation.core_cond,
        # self.two_d_axi.p_window[4, 1]+self.component.insulation.core_cond,
        # self.two_d_axi.p_window[7, 0]-self.insulation.core_cond, self.two_d_axi.p_window[7, 1]-self.insulation.core_cond)
        #
        # femm.mi_addblocklabel(self.two_d_axi.p_window[7, 0]-2*self.insulation.core_cond,
        # self.two_d_axi.p_window[7, 1]-2*self.insulation.core_cond)
        #
        # femm.mi_selectlabel(self.two_d_axi.p_window[7, 0]-2*self.insulation.core_cond,
        # self.two_d_axi.p_window[7, 1]-2*self.insulation.core_cond)
        #
        # femm.mi_setblockprop('Copper', 0, 1, 'icoil', 0, 0, 1)
        #
        # femm.mi_clearselected()

        for num in range(len(self.windings)):
            if self.windings[num].conductor_type in [ConductorType.RoundLitz, ConductorType.RoundSolid]:
                for i in range(0, int(self.two_d_axi.p_conductor[num].shape[0] / 5)):
                    # 0: center | 1: left | 2: top | 3: right | 4.bottom
                    femm.mi_drawarc(self.two_d_axi.p_conductor[num][5 * i + 1][0],
                                    self.two_d_axi.p_conductor[num][5 * i + 1][1],
                                    self.two_d_axi.p_conductor[num][5 * i + 3][0],
                                    self.two_d_axi.p_conductor[num][5 * i + 3][1], 180, 2.5)
                    femm.mi_addarc(self.two_d_axi.p_conductor[num][5 * i + 3][0],
                                self.two_d_axi.p_conductor[num][5 * i + 3][1],
                                self.two_d_axi.p_conductor[num][5 * i + 1][0],
                                self.two_d_axi.p_conductor[num][5 * i + 1][1], 180, 2.5)
                    femm.mi_addblocklabel(self.two_d_axi.p_conductor[num][5 * i][0],
                                        self.two_d_axi.p_conductor[num][5 * i][1])
                    femm.mi_selectlabel(self.two_d_axi.p_conductor[num][5 * i][0], self.two_d_axi.p_conductor[num][5 * i][1])
                    if num == 0:
                        femm.mi_setblockprop('Copper', 1, 0, 'Primary', 0, 2, 1)
                    if num == 1:
                        # femm.mi_setblockprop('Copper', 0, 1e-4, 'Secondary', 0, 3, 1)
                        femm.mi_setblockprop('Copper', 1, 0, 'Secondary', 0, 3, 1)
                    femm.mi_clearselected()

        # Define an "open" boundary condition using the built-in function:
        femm.mi_makeABC()
        """
        # Alternative BC
        region_add = 1.1

        femm.mi_drawrectangle(0, region_add*self.two_d_axi.p_outer[0][1], region_add*self.two_d_axi.p_outer[3][0], 
        region_add*self.two_d_axi.p_outer[3][1])
        # mi_addboundprop('Asymptotic', 0, 0, 0, 0, 0, 0, 1 / (para.mu_0 * bound.width), 0, 2); % Mixed
        femm.mi_addboundprop('Asymptotic', 0, 0, 0, 0, 1, 50, 0, 0, 1)
        femm.mi_selectsegment(region_add*self.two_d_axi.p_outer[3][0], region_add*self.two_d_axi.p_outer[3][1])
        femm.mi_setsegmentprop('Asymptotic', 1, 1, 0, 0)
        """

        # == Labels/Designations ==

        # Label for core
        femm.mi_addblocklabel(self.two_d_axi.p_outer[3, 0] - 0.001, self.two_d_axi.p_outer[3, 1] - 0.001)
        femm.mi_selectlabel(self.two_d_axi.p_outer[3, 0] - 0.001, self.two_d_axi.p_outer[3, 1] - 0.001)
        femm.mi_setblockprop('Ferrite', 1, 0, '<None>', 0, 0, 0)
        femm.mi_clearselected()

        # Labels for air
        if self.air_gaps.number == 0:
            femm.mi_addblocklabel(self.two_d_axi.r_inner - 0.0001, 0)
            femm.mi_selectlabel(self.two_d_axi.r_inner - 0.001, 0)
            femm.mi_setblockprop('Air', 1, 0, '<None>', 0, 1, 0)
            femm.mi_clearselected()
        else:
            femm.mi_addblocklabel(0.001, 0)
            femm.mi_selectlabel(0.001, 0)
            femm.mi_setblockprop('Air', 1, 0, '<None>', 0, 1, 0)
            femm.mi_clearselected()
        femm.mi_addblocklabel(self.two_d_axi.p_outer[3, 0] + 0.001, self.two_d_axi.p_outer[3, 1] + 0.001)
        femm.mi_selectlabel(self.two_d_axi.p_outer[3, 0] + 0.001, self.two_d_axi.p_outer[3, 1] + 0.001)
        femm.mi_setblockprop('Air', 1, 0, '<None>', 0, 1, 0)
        femm.mi_clearselected()

        # Now, the finished input geometry can be displayed.
        femm.mi_zoomnatural()
        femm.mi_saveas(os.path.join(self.file_data.femm_folder_path, 'coil.fem'))
        femm.mi_analyze()
        femm.mi_loadsolution()

        # == Log results ==
        self.write_femm_log()

        """
        # If we were interested in the flux density at specific positions,
        # we could inquire at specific points directly:
        b0 = femm.mo_getb(0, 0)
        ff.femmt_print('Flux density at the center of the bar is %g T' % np.abs(b0[1]))
        b1 = femm.mo_getb(0.01, 0.05)
        ff.femmt_print(f"Flux density at r=1cm, z=5cm is {np.abs(b1[1])} T")

        # The program will report the terminal properties of the circuit:
        # current, voltage, and flux linkage
        vals = femm.mo_getcircuitproperties('icoil')


        # [i, v, \[Phi]] = MOGetCircuitProperties["icoil"]

        # If we were interested in inductance, it could be obtained by
        # dividing flux linkage by current
        L = 1000 * np.abs(vals[2]) / np.abs(vals[0])
        ff.femmt_print('The self-inductance of the coil is %g mH' % L)

        # Or we could, for example, plot the results along a line using
        zee = []
        bee = []
        for n in range(-100, 101):
            b = femm.mo_getb(0.01, n)
            zee.append(n)
            bee.append(b[1])

        plt.plot(zee, bee)
        plt.ylabel('Flux Density, Tesla')
        plt.xlabel('Distance along the z-axis, mm')
        plt.title('Plot of flux density along the axis')
        plt.show()
        """
        # When the analysis is completed, FEMM can be shut down.
        # femm.closefemm()

    def write_femm_log(self):
        """TODO Doc
        """
        if os.name == 'nt':
            ff.install_pyfemm_if_missing()
            if not self.femm_is_imported:
                globals()["femm"] = __import__("femm")
                self.femm_is_imported = True
        else:
            raise Exception('You are using a computer that is not running windows. '
                            'This command is only executable on Windows computers.')

        file = open(self.file_data.femm_results_log_path, 'w+', encoding='utf-8')

        log = {}

        # ff.femmt_print(hyst_loss)
        # tmp = femm.mo_getcircuitproperties('Primary')
        # ff.femmt_print(tmp)
        # self.tot_loss_femm = 0.5 * tmp[0] * tmp[1]
        # ff.femmt_print(self.tot_loss_femm)

        # Write Circuit Properties
        # log["Circuit Properties"] = femm.mo_getcircuitproperties('Primary')

        # Write Hysteresis Losses
        femm.mo_groupselectblock(0)
        log["Hysteresis Losses"] = femm.mo_blockintegral(3)
        femm.mo_clearblock()

        # Primary Winding Ciruit Properties
        circuit_properties_primary = femm.mo_getcircuitproperties('Primary')
        log["Primary Current"] = circuit_properties_primary[0]
        log["Primary Voltage"] = [circuit_properties_primary[1].real, circuit_properties_primary[1].imag]
        log["Primary Flux"] = [circuit_properties_primary[2].real, circuit_properties_primary[2].imag]
        log["Primary Self Inductance"] = [circuit_properties_primary[2].real / circuit_properties_primary[0],
                                        circuit_properties_primary[2].imag / circuit_properties_primary[0]]
        log["Primary Mean Power"] = [0.5*circuit_properties_primary[1].real*circuit_properties_primary[0],
                                    0.5*circuit_properties_primary[1].imag*circuit_properties_primary[0]]

        # Primary Winding Losses (with group n=2) by field intergation
        femm.mo_groupselectblock(2)
        log["Primary Winding Losses"] = femm.mo_blockintegral(6).real
        femm.mo_clearblock()

        if self.component_type == (ComponentType.Transformer or ComponentType.IntegratedTransformer):
            # secondary Winding Ciruit Properties
            circuit_properties_secondary = femm.mo_getcircuitproperties('Secondary')
            log["Secondary Current"] = circuit_properties_secondary[0]
            log["Secondary Voltage"] = [circuit_properties_secondary[1].real, circuit_properties_secondary[1].imag]
            log["Secondary Flux"] = [circuit_properties_secondary[2].real, circuit_properties_secondary[2].imag]
            log["Secondary Self Inductance"] = [circuit_properties_secondary[2].real / circuit_properties_secondary[0],
                                            circuit_properties_secondary[2].imag / circuit_properties_secondary[0]]
            log["Secondary Mean Power"] = [0.5 * circuit_properties_secondary[1].real * circuit_properties_secondary[0],
                                        0.5 * circuit_properties_secondary[1].imag * circuit_properties_secondary[0]]

            # secondary Winding Losses (with group n=2) by field intergation
            femm.mo_groupselectblock(3)
            log["Secondary Winding Losses"] = femm.mo_blockintegral(6).real
            femm.mo_clearblock()


        json.dump(log, file, indent=2, ensure_ascii=False)
        file.close()


    @staticmethod
    def calculate_point_average(x1: float, y1: float, x2: float, y2: float) -> List[float]:
        """Calculates the middle point between two given points.

        :param x1: Point1 x
        :type x1: float
        :param y1: Point1 y
        :type y1: float
        :param x2: Point2 x
        :type x2: float
        :param y2: Point2 y
        :type y2: float
        :return: Average Point x, y
        :rtype: float, float
        """
        # TODO Move to femmt_functions
        return (x1 + x2) / 2, (y1 + y2) / 2

    def femm_thermal_validation(self, thermal_conductivity_dict: Dict, boundary_temperature: Dict, case_gap_top: float, case_gap_right: float, case_gap_bot: float):
        """Creates a thermal model in femm and simulates it with the given thermal conductivities

        :param thermal_conductivity_dict: Dict containing conductivities for air, winding, case, core
        :type thermal_conductivity_dict: Dict
        :param boundary_temperature: Dict containing temperatures on boundary lines
        :type boundary_temperature: Dict
        :param case_gap_top: Length top case
        :type case_gap_top: float
        :param case_gap_right: Length right case
        :type case_gap_right: float
        :param case_gap_bot: Length bot case
        :type case_gap_bot: float
        """
        # Optional usage of FEMM tool by David Meeker
        # 2D Mesh and FEM interfaces (only for windows machines)
        if os.name == 'nt':
            ff.install_pyfemm_if_missing()
            if not self.femm_is_imported:
                globals()["femm"] = __import__("femm")
                self.femm_is_imported = True
        else:
            raise Exception('You are using a computer that is not running windows. '
                'This command is only executable on Windows computers.')

        # Get paths
        femm_model_file_path = os.path.join(self.file_data.femm_folder_path, "thermal-validation.FEH")

        self.file_data.create_folders(self.file_data.femm_folder_path)

        # Extract losses
        losses = read_log(self.file_data.e_m_results_log_path)

        # Extract wire_radii
        wire_radii = [winding.conductor_radius for winding in self.windings]

        # == Init ==
        femm.openfemm(0)
        femm.newdocument(2)
        femm.hi_probdef("meters", "axi", 1.e-8, 0)

        # == Materials ==
        # Core
        k_core = thermal_conductivity_dict["core"]
        q_vol_core = losses["core"] / self.calculate_core_volume()
        # c_core = 0.007
        c_core = 0

        # Air
        k_air = thermal_conductivity_dict["air"]
        q_vol_air = 0
        # c_air = 1004.7
        c_air = 0

        # Wire
        k_wire = thermal_conductivity_dict["winding"]
        # c_wire = 385
        c_wire = 0

        # Case
        k_case = thermal_conductivity_dict["case"]["top"] # Does not matter when the regions all have the same thermal coductivity.
        q_vol_case = 0
        # c_case = 0.01
        c_case = 0

        # Air gap
        k_air_gap = thermal_conductivity_dict["air_gaps"]
        q_vol_air_gap = 0
        c_air_gap = 0

        # Setup winding list
        winding_losses_list = []
        for i in range(1, 3):
            key = f"winding{i}"
            inner_winding_list = []
            if key in losses:
                for winding in losses[key]["turns"]:
                    inner_winding_list.append(winding)
            winding_losses_list.append(inner_winding_list)

        # Setup materials
        femm.hi_addmaterial('Core', k_core, k_core, q_vol_core, c_core)
        femm.hi_addmaterial('Air', k_air, k_air, q_vol_air, c_air)
        femm.hi_addmaterial('Air Gaps', k_air_gap, k_air_gap, q_vol_air_gap, c_air_gap)
        wire_distances = self.get_wire_distances()
        for winding_index, winding in enumerate(winding_losses_list):
            for i in range(len(winding)):
                femm.hi_addmaterial(f'Wire_{winding_index}_{i}', k_wire, k_wire, calculate_heat_flux_round_wire(winding[i], wire_radii[winding_index], wire_distances[winding_index][i]), c_wire)
        femm.hi_addmaterial('Case', k_case, k_case, q_vol_case, c_case)

        # Add boundary condition
        femm.hi_addboundprop("Boundary", 0, boundary_temperature, 0, 0, 0, 0)
        femm.hi_addboundprop("NeumannBoundary", 1, 0, 0, 0, 0, 0)

        # == Geometry ==
        self.high_level_geo_gen()
        # Add core
        if self.air_gaps.number == 0:
            femm.hi_drawline(self.two_d_axi.p_window[4, 0],
                            self.two_d_axi.p_window[4, 1],
                            self.two_d_axi.p_window[5, 0],
                            self.two_d_axi.p_window[5, 1])
            femm.hi_drawline(self.two_d_axi.p_window[5, 0],
                            self.two_d_axi.p_window[5, 1],
                            self.two_d_axi.p_window[7, 0],
                            self.two_d_axi.p_window[7, 1])
            femm.hi_drawline(self.two_d_axi.p_window[7, 0],
                            self.two_d_axi.p_window[7, 1],
                            self.two_d_axi.p_window[6, 0],
                            self.two_d_axi.p_window[6, 1])
            femm.hi_drawline(self.two_d_axi.p_window[6, 0],
                            self.two_d_axi.p_window[6, 1],
                            self.two_d_axi.p_window[4, 0],
                            self.two_d_axi.p_window[4, 1])
            femm.hi_drawline(0,
                            self.two_d_axi.p_outer[2, 1],
                            self.two_d_axi.p_outer[3, 0],
                            self.two_d_axi.p_outer[3, 1])
            femm.hi_drawline(self.two_d_axi.p_outer[3, 0],
                            self.two_d_axi.p_outer[3, 1],
                            self.two_d_axi.p_outer[1, 0],
                            self.two_d_axi.p_outer[1, 1])
            femm.hi_drawline(self.two_d_axi.p_outer[1, 0],
                            self.two_d_axi.p_outer[1, 1],
                            0,
                            self.two_d_axi.p_outer[0, 1])
            femm.hi_drawline(0,
                            self.two_d_axi.p_outer[0, 1],
                            0,
                            self.two_d_axi.p_outer[2, 1])
        elif self.air_gaps.number > 0:
            femm.hi_drawline(0,
                            self.two_d_axi.p_air_gaps[0, 1],
                            self.two_d_axi.p_air_gaps[1, 0],
                            self.two_d_axi.p_air_gaps[1, 1])
            femm.hi_drawline(self.two_d_axi.p_air_gaps[1, 0],
                            self.two_d_axi.p_air_gaps[1, 1],
                            self.two_d_axi.p_window[4, 0],
                            self.two_d_axi.p_window[4, 1])
            femm.hi_drawline(self.two_d_axi.p_window[4, 0],
                            self.two_d_axi.p_window[4, 1],
                            self.two_d_axi.p_window[5, 0],
                            self.two_d_axi.p_window[5, 1])
            femm.hi_drawline(self.two_d_axi.p_window[5, 0],
                            self.two_d_axi.p_window[5, 1],
                            self.two_d_axi.p_window[7, 0],
                            self.two_d_axi.p_window[7, 1])
            femm.hi_drawline(self.two_d_axi.p_window[7, 0],
                            self.two_d_axi.p_window[7, 1],
                            self.two_d_axi.p_window[6, 0],
                            self.two_d_axi.p_window[6, 1])
            femm.hi_drawline(self.two_d_axi.p_window[6, 0],
                            self.two_d_axi.p_window[6, 1],
                            self.two_d_axi.p_air_gaps[3, 0],
                            self.two_d_axi.p_air_gaps[3, 1])
            femm.hi_drawline(self.two_d_axi.p_air_gaps[3, 0],
                            self.two_d_axi.p_air_gaps[3, 1],
                            0,
                            self.two_d_axi.p_air_gaps[2, 1])
            femm.hi_drawline(0,
                            self.two_d_axi.p_air_gaps[2, 1],
                            0,
                            self.two_d_axi.p_outer[2, 1])
            femm.hi_drawline(0,
                            self.two_d_axi.p_outer[2, 1],
                            self.two_d_axi.p_outer[3, 0],
                            self.two_d_axi.p_outer[3, 1])
            femm.hi_drawline(self.two_d_axi.p_outer[3, 0],
                            self.two_d_axi.p_outer[3, 1],
                            self.two_d_axi.p_outer[1, 0],
                            self.two_d_axi.p_outer[1, 1])
            femm.hi_drawline(self.two_d_axi.p_outer[1, 0],
                            self.two_d_axi.p_outer[1, 1],
                            0,
                            self.two_d_axi.p_outer[0, 1])
            femm.hi_drawline(0,
                            self.two_d_axi.p_outer[0, 1],
                            0,
                            self.two_d_axi.p_air_gaps[0, 1])

            # In order for the simulation to work the air_gap must be closed:
            femm.hi_drawline(0, self.two_d_axi.p_air_gaps[0, 1], 0, self.two_d_axi.p_air_gaps[2, 1])

            # Close air gap to seperate from air
            femm.hi_drawline(self.two_d_axi.p_air_gaps[1, 0], self.two_d_axi.p_air_gaps[1, 1],
                                self.two_d_axi.p_air_gaps[3, 0], self.two_d_axi.p_air_gaps[3, 1])
        else:
            raise Exception("Negative air gap number is not allowed")

        # Add case
        femm.hi_drawline(0, self.two_d_axi.p_outer[2, 1], 0, self.two_d_axi.p_outer[2, 1] + case_gap_top)  # Top left line
        femm.hi_drawline(0, self.two_d_axi.p_outer[2, 1] + case_gap_top, self.two_d_axi.p_outer[3, 0] + case_gap_right, self.two_d_axi.p_outer[3, 1] + case_gap_top)  # Top line
        femm.hi_drawline(self.two_d_axi.p_outer[3, 0] + case_gap_right, self.two_d_axi.p_outer[3, 1] + case_gap_top, self.two_d_axi.p_outer[1, 0] + case_gap_right,
                        self.two_d_axi.p_outer[1, 1] - case_gap_bot)  # Right line
        femm.hi_drawline(self.two_d_axi.p_outer[1, 0] + case_gap_right, self.two_d_axi.p_outer[1, 1] - case_gap_bot, 0, self.two_d_axi.p_outer[0, 1] - case_gap_bot)  # Bottom line
        femm.hi_drawline(0, self.two_d_axi.p_outer[0, 1] - case_gap_bot, 0, self.two_d_axi.p_outer[0, 1])  # Bottom right line

        # Create boundary
        # femm.hi_selectsegment(*self.calculatePointAverage(0, self.two_d_axi.p_outer[2, 1], 0, self.two_d_axi.p_outer[2, 1] + caseGapTop))
        femm.hi_selectsegment(*MagneticComponent.calculate_point_average(0, self.two_d_axi.p_outer[2, 1] + case_gap_top, self.two_d_axi.p_outer[3, 0] + case_gap_right,
                                                                        self.two_d_axi.p_outer[3, 1] + case_gap_top))
        femm.hi_setsegmentprop("NeumannBoundary", 0, 1, 0, 2, "<None>")
        femm.hi_clearselected()

        femm.hi_selectsegment(*MagneticComponent.calculate_point_average(self.two_d_axi.p_outer[3, 0] + case_gap_right, self.two_d_axi.p_outer[3, 1] + case_gap_top,
                                                                        self.two_d_axi.p_outer[1, 0] + case_gap_right, self.two_d_axi.p_outer[1, 1] - case_gap_bot))
        femm.hi_selectsegment(*MagneticComponent.calculate_point_average(self.two_d_axi.p_outer[1, 0] + case_gap_right, self.two_d_axi.p_outer[1, 1] - case_gap_bot, 0,
                                                                        self.two_d_axi.p_outer[0, 1] - case_gap_bot))
        # femm.hi_selectsegment(*self.calculatePointAverage(0, self.two_d_axi.p_outer[0, 1] - caseGapBot, 0, self.two_d_axi.p_outer[0, 1]))
        femm.hi_setsegmentprop("Boundary", 0, 1, 0, 2, "<None>")
        femm.hi_clearselected()

        # Add case material
        material_x, material_y = self.calculate_point_average(0, self.two_d_axi.p_outer[2, 1], 0, self.two_d_axi.p_outer[2, 1] + case_gap_top)
        femm.hi_addblocklabel(material_x + 0.001, material_y)
        femm.hi_selectlabel(material_x + 0.001, material_y)
        femm.hi_setblockprop('Case', 1, 0, 0)
        femm.hi_clearselected()

        # Add Coil
        for num in range(len(self.windings)):
            for i in range(0, int(self.two_d_axi.p_conductor[num].shape[0] / 5)):
                # 0: center | 1: left | 2: top | 3: right | 4.bottom
                femm.hi_drawarc(self.two_d_axi.p_conductor[num][5 * i + 1][0],
                                self.two_d_axi.p_conductor[num][5 * i + 1][1],
                                self.two_d_axi.p_conductor[num][5 * i + 3][0],
                                self.two_d_axi.p_conductor[num][5 * i + 3][1], 180, 2.5)
                femm.hi_addarc(self.two_d_axi.p_conductor[num][5 * i + 3][0],
                            self.two_d_axi.p_conductor[num][5 * i + 3][1],
                            self.two_d_axi.p_conductor[num][5 * i + 1][0],
                            self.two_d_axi.p_conductor[num][5 * i + 1][1], 180, 2.5)
                femm.hi_addblocklabel(self.two_d_axi.p_conductor[num][5 * i][0],
                                    self.two_d_axi.p_conductor[num][5 * i][1])
                femm.hi_selectlabel(self.two_d_axi.p_conductor[num][5 * i][0], self.two_d_axi.p_conductor[num][5 * i][1])
                if num == 0:
                    femm.hi_setblockprop(f'Wire_{num}_{i}', 1, 0, 1)
                if num == 1:
                    femm.hi_setblockprop(f'Wire_{num}_{i}', 1, 0, 1)
                femm.hi_clearselected()

        # Define an "open" boundary condition using the built-in function:
        # femm.hi_makeABC()

        # == Labels/Designations ==
        # Label for air and air gap
        if self.air_gaps.number == 0:
            femm.hi_addblocklabel(self.two_d_axi.r_inner - 0.0001, 0)
            femm.hi_selectlabel(self.two_d_axi.r_inner - 0.001, 0)
            femm.hi_setblockprop('Air', 1, 0, 0)
            femm.hi_clearselected()
        else:
            femm.hi_addblocklabel(self.two_d_axi.r_inner - 0.0001, 0)
            femm.hi_selectlabel(self.two_d_axi.r_inner - 0.001, 0)
            femm.hi_setblockprop('Air', 1, 0, 0)
            femm.hi_clearselected()
            femm.hi_addblocklabel(0.001, 0)
            femm.hi_selectlabel(0.001, 0)
            femm.hi_setblockprop('Air Gaps', 1, 0, 0)
            femm.hi_clearselected()

        # Label for core
        femm.hi_addblocklabel(self.two_d_axi.p_outer[3, 0] - 0.001, self.two_d_axi.p_outer[3, 1] - 0.001)
        femm.hi_selectlabel(self.two_d_axi.p_outer[3, 0] - 0.001, self.two_d_axi.p_outer[3, 1] - 0.001)
        femm.hi_setblockprop('Core', 1, 0, 0)
        femm.hi_clearselected()

        # Not needed when the core is the boundary
        # femm.hi_addblocklabel(self.two_d_axi.p_outer[3, 0] + 0.001, self.two_d_axi.p_outer[3, 1] + 0.001)
        # femm.hi_selectlabel(self.two_d_axi.p_outer[3, 0] + 0.001, self.two_d_axi.p_outer[3, 1] + 0.001)
        # femm.hi_setblockprop('Air', 1, 0, 0)
        # femm.hi_clearselected()

        # Now, the finished input geometry can be displayed.
        femm.hi_zoomnatural()
        femm.hi_saveas(femm_model_file_path)
        femm.hi_analyze()
        femm.hi_loadsolution()
        input()  # So the window stays open
        # femm.closefemm()

    @staticmethod
    def encode_settings(o) -> Dict:
        """Encoes the magnetic component in a dictionary.

        :param o: Magnetic component containing the model.
        :type o: MagneticComponent
        :return: Model encodede as dictionary
        :rtype: Dict
        """
        content = {
            "date": datetime.today().strftime('%Y-%m-%d %H:%M:%S'),
            "component_type": o.component_type.name,
            "working_directory": o.file_data.working_directory,
            "core": o.core.to_dict(),
            "virtual_winding_windows": [vww.to_dict() for vww in o.virtual_winding_windows],
        }

        if o.air_gaps is not None:
            content["air_gaps"] = o.air_gaps.to_dict()
        
        if o.insulation is not None:
            content["insulation"] = o.insulation.to_dict()

        if o.stray_path is not None:
            content["stray_path"] = o.stray_path.__dict__
        
        return content

    @staticmethod    
    def decode_settings_from_log(log_file_path: str, working_directory: str = None):
        """Reads the given log and returns the magnetic component from th elog.

        :param log_file_path: Path to the log file
        :type log_file_path: str
        :param working_directory: If the working directory shall be a different than from the log file enter a new one here, defaults to None
        :type working_directory: str, optional
        :return: Magnetic component containing the model
        :rtype: MagneticComponent
        """
        if not os.path.isfile(log_file_path):
            raise Exception(f"File {log_file_path} does not exists or is not a file!")

        settings = None
        with open(log_file_path, "r") as fd:
            content = json.load(fd)
            settings = content["simulation_settings"]

        if settings is not None:
            cwd = working_directory if working_directory is not None else settings["working_directory"]
            geo = MagneticComponent(component_type=ComponentType[settings["component_type"]], working_directory=cwd)

            settings["core"]["loss_approach"] = LossApproach[settings["core"]["loss_approach"]]
            core = Core(**settings["core"])
            geo.set_core(core)

            if "air_gaps" in settings:
                air_gaps = AirGaps(AirGapMethod[settings["air_gaps"]["method"]], core)
                for air_gap in settings["air_gaps"]["air_gaps"]:
                    air_gaps.add_air_gap(AirGapLegPosition[air_gap["leg_position"]], air_gap["height"], air_gap["position_value"],)
                geo.set_air_gaps(air_gaps)

            if "insulation" in settings:
                insulation = Insulation()
                insulation.add_core_insulations(*settings["insulation"]["core_insulations"])
                insulation.add_winding_insulations(settings["insulation"]["inner_winding_insulations"], settings["insulation"]["vww_insulation"])
                geo.set_insulation(insulation)

            if "stray_path" in settings:
                stray_path = StrayPath(*settings["stray_path"])
                geo.set_stray_path(stray_path)

            virtual_winding_windows = settings["virtual_winding_windows"]
            new_virtual_winding_windows = []
            for vww in virtual_winding_windows:
                turns = vww["turns"]
                conductors = []
                for winding in vww["windings"]:
                    conductor = Conductor(winding["winding_number"], Conductivity[winding["conductivity"]])
                    conductor_type = ConductorType[winding["conductor_type"]]
                    if conductor_type == ConductorType.RectangularSolid:
                        conductor.set_rectangular_conductor(winding["thickness"])
                    elif conductor_type == ConductorType.RoundLitz:
                        conductor.set_litz_round_conductor(winding["conductor_radius"], winding["number_strands"], 
                        winding["strand_radius"], winding["fill_factor"], ConductorArrangement[winding["conductor_arrangement"]])
                    elif conductor_type == ConductorType.RoundSolid:
                        conductor.set_solid_round_conductor(winding["conductor_radius"], ConductorArrangement[winding["conductor_arrangement"]])
                    else:
                        raise Exception(f"Unknown conductor type {conductor_type.name}")
                    
                    conductors.append(conductor)

                new_vww = VirtualWindingWindow(vww["bot_bound"], vww["top_bound"], vww["left_bound"], vww["right_bound"])
                winding_type = WindingType[vww["winding_type"]]
                if winding_type == WindingType.Single:
                    winding_scheme = WindingScheme[vww["winding_scheme"]] if vww["winding_scheme"] is not None else None
                    wrap_para_type = WrapParaType[vww["wrap_para"]] if vww["wrap_para"] is not None else None
                    new_vww.set_winding(conductors[0], turns[0], winding_scheme, wrap_para_type)
                elif winding_type == WindingType.Interleaved:
                    new_vww.set_interleaved_winding(conductors[0], turns[0], conductors[1], turns[1], InterleavedWindingScheme[vww["winding_scheme"]], vww["winding_insulation"])
                else:
                    raise Exception(f"Winding type {winding_type} is not implemented")

                new_virtual_winding_windows.append(new_vww)

            winding_window = WindingWindow(core, insulation)
            winding_window.virtual_winding_windows = new_virtual_winding_windows
            geo.set_winding_window(winding_window)

            return geo

        raise Exception(f"Couldn't extract settings from file {log_file_path}")