# Python standard libraries
import csv
import fileinput
import os
import gmsh
import json
import warnings
import inspect
import time
import logging
from typing import List, Dict, Optional, Tuple
from datetime import datetime
import dataclasses
from matplotlib import pyplot as plt

# Third party libraries
from onelab import onelab
import materialdatabase as mdb
import numpy as np

# Local libraries
import femmt.functions as ff
from femmt.constants import *
from femmt.mesh import Mesh
from femmt.model import VirtualWindingWindow, WindingWindow, Core, Insulation, StrayPath, AirGaps, Conductor
from femmt.enumerations import *
from femmt.data import FileData, MeshData
from femmt.drawing import TwoDaxiSymmetric
from femmt.thermal import thermal_simulation, calculate_heat_flux_round_wire, read_results_log
from femmt.dtos import *
import femmt.functions_reluctance as fr


class MagneticComponent:
    """
    A MagneticComponent is the main object for all simulation purposes in femmt.

        - One or more "MagneticComponents" can be created
        - Each "MagneticComponent" owns its own instance variable values

    """
    # Initialization of all class variables
    # Common variables for all instances

    onelab_folder_path: str = None
    silent: bool = False

    def __init__(self, component_type: ComponentType = ComponentType.Inductor, working_directory: str = None,
                 verbosity: Verbosity = 2, is_gui: bool = False, simulation_name: Optional[str] = None):
        # TODO Add a enum? for the verbosity to combine silent and print_output_to_file variables
        """
        :param component_type: Available options:
                               - "inductor"
                               - "transformer"
                               - "integrated_transformer" (Transformer with included stray-path)
        :type component_type: ComponentType
        :param working_directory: Sets the working directory
        :type working_directory: string
        :param silent: True to reduce onelab simulation outputs in command line. False to see full command line output.
        :type silent: bool
        :param is_gui: Asks at first startup for onelab-path. Distinction between GUI and command line. Defaults to 'False' in command-line-mode.
        :type is_gui: bool
        :param simulation_name: name without any effect. Will just be displayed in the result-log file
        :type simulation_name: str
        """
        # Get caller filepath when no working_directory was set
        if working_directory is None:
            caller_filename = inspect.stack()[1].filename
            working_directory = os.path.join(os.path.dirname(caller_filename), "femmt")

        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # Create file paths class in order to handle all paths
        self.file_data = FileData(working_directory)
        # Clear result folder structure in case of missing
        self.file_data.clear_previous_simulation_results()

        # Variable to set silent mode
        self.verbosity = verbosity
        self.logger = logging.getLogger("FEMMTLogger")
        self.logger.setLevel(logging.INFO)
        if not gmsh.isInitialized():
            gmsh.initialize()

        if not verbosity == Verbosity.ToConsole:
            gmsh.option.setNumber("General.Terminal", 0)
            self.silent = True

        if verbosity == Verbosity.ToFile:
            fh = logging.FileHandler(self.file_data.femmt_log, mode="w")
            fh.setLevel(logging.INFO)
            self.logger.addHandler(fh)
            self.silent = True

        self.femmt_print(f"\n"
                         f"Initialized a new Magnetic Component of type {component_type.name}\n"
                         f"--- --- --- ---")

        # To make sure femm is only imported once
        self.femm_is_imported = False

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Component Geometry
        self.component_type = component_type  # "inductor", "transformer", "integrated_transformer" (or "three-phase-transformer")

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Components
        self.core = None  # Contains all information about the cores
        self.air_gaps = None  # Contains every air gap
        self.windings = None  # List of the different winding objects which the following structure: windings[0]: primary, windings[1]: secondary, windings[2]: tertiary ....
        self.insulation = None  # Contains information about the needed insulations
        self.winding_windows = None  # Contains a list of every winding_window which was created containing a list of virtual_winding_windows
        self.stray_path = None  # Contains information about the stray_path (only for integrated transformers)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Control Flags
        self.plot_fields = "standard"  # can be "standard" or False

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Excitation Parameters
        # Empty lists will be set when a winding window is added to the magnetic component
        self.imposed_reduced_frequency = None
        self.flag_excitation_type = None
        self.current = []  # Defined for every conductor
        self.current_density = []  # Defined for every conductor
        self.voltage = []  # Defined for every conductor
        self.frequency = None
        self.phase_deg = None  # Default is zero, Defined for every conductor
        self.red_freq = None  # [] * self.n_windings  # Defined for every conductor
        self.delta = None

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Steinmetz loss material coefficients and current waveform
        self.Ipeak = None
        self.ki = None
        self.alpha = None
        self.beta = None
        self.t_rise = None
        self.t_fall = None
        self.f_switch = None

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # MeshData to store the mesh size for different points
        # Object is added in set_core
        self.mesh_data = None
        self.mesh = None
        self.two_d_axi = None

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # -- Used for Litz Validation --
        self.sweep_frequencies = None

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        # 2 and 3 winding transformer inductance matrix
        self.L_1_1 = None
        self.L_2_2 = None
        self.L_3_3 = None
        self.M = None
        self.M_12 = None
        self.M_21 = None
        self.M_13 = None
        self.M_32 = None
        self.M_32 = None
        self.M_23 = None
        self.Pv = None
        # 2 and 3 winding transformer primary concentrated equivalent circuit
        self.n_conc = None
        self.L_s_conc = None
        self.L_h_conc = None
        self.L_s1 = None
        self.L_s2 = None
        self.L_s3 = None
        self.L_h = None
        self.L_s12 = None
        self.L_s13 = None
        self.L_s23 = None
        self.n_12 = None
        self.n_13 = None

        # -- FEMM variables --
        self.tot_loss_femm = None

        self.onelab_setup(is_gui)
        self.onelab_client = onelab.client(__file__)
        self.simulation_name = simulation_name

    def femmt_print(self, text: str):
        if self.verbosity != Verbosity.Silent:
            self.logger.info(text)

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   -  -  -  -  -  -  -  -  -  -  -
    # Thermal simulation
    def thermal_simulation(self, thermal_conductivity_dict: Dict, boundary_temperatures_dict: Dict,
                           boundary_flags_dict: Dict, case_gap_top: float,
                           case_gap_right: float, case_gap_bot: float, show_thermal_simulation_results: bool = True,
                           pre_visualize_geometry: bool = False, color_scheme: Dict = ff.colors_femmt_default,
                           colors_geometry: Dict = ff.colors_geometry_femmt_default, flag_insulation: bool = True):
        """
        Starts the thermal simulation using thermal_simulation.py

        :param thermal_conductivity_dict: Contains the thermal conductivities for every region
        :type thermal_conductivity_dict: Dict
        :param boundary_temperatures_dict: Contains the temperatures at each boundary line
        :type boundary_temperatures_dict: Dict
        :param boundary_flags_dict: Sets the boundary type (dirichlet or von neumann) for each boundary line
        :type boundary_flags_dict: Dict
        :param case_gap_top: Size of the top case
        :type case_gap_top: float
        :param case_gap_right: Size of the right case
        :type case_gap_right: float
        :param case_gap_bot: Size of the bot case
        :type case_gap_bot: float
        :param show_thermal_simulation_results: Shows thermal results in gmsh, defaults to True
        :type show_thermal_simulation_results: bool, optional
        :param pre_visualize_geometry: Shows the thermal model before simulation, defaults to False
        :type pre_visualize_geometry: bool, optional
        :param color_scheme: Color scheme for visualization, defaults to ff.colors_femmt_default
        :type color_scheme: Dict, optional
        :param colors_geometry: Color geometry for visualization, defaults to ff.colors_geometry_femmt_default
        :type colors_geometry: Dict, optional
        """
        # Create necessary folders
        self.file_data.create_folders(self.file_data.thermal_results_folder_path)

        self.mesh.generate_thermal_mesh(case_gap_top, case_gap_right, case_gap_bot, color_scheme, colors_geometry,
                                        pre_visualize_geometry)

        # insulation_tag = self.mesh.ps_insulation if flag_insulation and len(self.insulation.core_cond) == 4 else None

        if not os.path.exists(self.file_data.e_m_results_log_path):
            # Simulation results file not created
            raise Exception(
                "Cannot run thermal simulation -> Magnetic simulation needs to run first (no results_log.json found")

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
            "core_tags": self.mesh.ps_core,
            "background_tag": self.mesh.ps_air,
            "winding_tags": self.mesh.ps_cond,
            "air_gaps_tag": self.mesh.ps_air_gaps if self.air_gaps.number > 0 else None,
            "boundary_regions": self.mesh.thermal_boundary_region_tags,
            "insulations_tag": self.mesh.ps_insulation if flag_insulation and len(
                self.insulation.core_cond) == 4 else None
        }

        # Core area -> Is needed to estimate the heat flux
        # Power density for volumes W/m^3
        #core_area = self.calculate_core_volume()
        core_area = self.calculate_core_parts_volume()




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
            "show_thermal_fem_results": show_thermal_simulation_results,
            "print_sensor_values": False,
            "silent": self.verbosity == Verbosity.Silent,  # Add verbosity for thermal simulation
            "flag_insulation": flag_insulation
        }

        thermal_simulation.run_thermal(**thermal_parameters)

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   -  -  -  -  -  -  -  -  -  -  -
    # Setup
    def onelab_setup(self, is_gui: bool):
        """
        Either reads ONELAB parent folder path from config.json or asks the user to provide the ONELAB path it.
        Creates a config.json inside the site-packages folder at first run.

        :param is_gui: set to True to avoid terminal output question for onelab file path at first run. Used especially in GUI
        :type is_gui: bool
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
            onelab_path = os.path.normpath(input(
                "Enter the path of onelab's parent folder (path to folder which contains getdp, onelab executables): "))

            if os.path.exists(onelab_path):
                onelab_path_wrong = False
                break
            else:
                self.femmt_print('onelab not found! Tool searches for onelab.py in the folder. Please re-enter path!')
        self.file_data.onelab_folder_path = onelab_path

        # Write the path to the config.json
        onelab_path_dict = {"onelab": onelab_path}
        with open(os.path.join(self.file_data.config_path), 'w', encoding='utf-8') as fd:
            json.dump(onelab_path_dict, fd, indent=2, ensure_ascii=False)

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   -  -  -  -  -  -  -  -  -  -  -
    # Geometry Parts
    def high_level_geo_gen(self, frequency: float = None, skin_mesh_factor: float = None):
        """ Updates the mesh data and creates the model and mesh objects

        :param frequency: Frequency used in the mesh density, defaults to None
        :type frequency: float, optional
        :param skin_mesh_factor: Used in the mesh density, defaults to None
        :type skin_mesh_factor: float, optional
        """
        # Update mesh data
        self.mesh_data.update_data(frequency, skin_mesh_factor)

        # Create model
        self.two_d_axi = TwoDaxiSymmetric(self.core, self.mesh_data, self.air_gaps, self.winding_windows,
                                          self.stray_path,
                                          self.insulation, self.component_type, len(self.windings), self.verbosity,
                                          self.logger)
        self.two_d_axi.draw_model()

        # Create mesh
        self.mesh = Mesh(self.two_d_axi, self.windings, self.core.correct_outer_leg, self.file_data, self.verbosity,
                         self.logger, None)
        # self.mesh = Mesh(self.two_d_axi, self.windings, self.core.correct_outer_leg, self.file_data, None, ff.silent)

    def mesh(self, frequency: float = None, skin_mesh_factor: float = None):
        """Generates model and mesh.

        :param frequency: Frequency used in the mesh density, defaults to None
        :type frequency: float, optional
        :param skin_mesh_factor: Used in the mesh density, defaults to None
        :type skin_mesh_factor: float, optional
        """
        self.high_level_geo_gen(frequency=frequency, skin_mesh_factor=skin_mesh_factor)
        self.mesh.generate_hybrid_mesh()  # create the mesh itself with gmsh
        self.mesh.generate_electro_magnetic_mesh()  # assign the physical entities/domains to the mesh

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   -  -  -  -  -  -  -  -  -  -  -
    # Create Model
    def set_insulation(self, insulation: Insulation):
        """Adds the insulation to the model

        :param insulation: insulation object
        :type insulation: Insulation
        """
        if insulation.cond_cond is None or not insulation.cond_cond:
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

    def set_winding_windows(self, winding_windows: List[WindingWindow], mesh_accuracy: float = 0.5):
        """
        Adds the winding windows to the model. Creates the windings list, which contains the conductors
        from the virtual winding windows but sorted by the winding_number (ascending).
        Sets empty lists for excitation parameters

        :param winding_windows: List of WindingWindow objects
        :type winding_windows: List[WindingWindow]
        :param mesh_accuracy: a mesh_accuracy of 0.5 is recommended. Do not change this parameter, except performing thousands of simulations, e.g. a Pareto optimization. In this case, the value can be set e.g. to 0.8
        :type mesh_accuracy: float
        """
        self.winding_windows = winding_windows
        windings = []
        for ww in winding_windows:
            for vww in ww.virtual_winding_windows:
                if not vww.winding_is_set:
                    raise Exception("Each virtual winding window needs to have a winding")
                for winding in vww.windings:
                    if winding not in windings:
                        windings.append(winding)

        # print(f"{winding_window.virtual_winding_windows = }")
        # print(f"{windings = }")
        self.windings = sorted(windings, key=lambda x: x.winding_number)

        # Print statement was moved here so the silence functionality is not needed in Conductors class.
        # TODO Can this be even removed?
        for winding in self.windings:
            if winding.conductor_type == ConductorType.RoundLitz:
                self.femmt_print(f"Updated Litz Configuration: \n"
                                 f" ff: {winding.ff} \n"
                                 f" Number of layers/strands: {winding.n_layers}/{winding.n_strands} \n"
                                 f" Strand radius: {winding.strand_radius} \n"
                                 f" Conductor radius: {winding.conductor_radius}\n"
                                 f"---")

        # Set excitation parameter lists
        self.current = [None] * len(windings)
        self.current_density = [None] * len(windings)
        self.voltage = [None] * len(windings)
        self.phase_deg = np.zeros(len(windings))

        # Correct the turns lists in each vww, so that they have the same length
        for ww in winding_windows:
            for vww in ww.virtual_winding_windows:
                zeros_to_append = (len(self.windings) - len(vww.turns))
                if zeros_to_append < 0:
                    for i in range(0, -zeros_to_append):
                        vww.turns.pop()
                else:
                    for i in range(0, zeros_to_append):
                        vww.turns.append(0)

        # Default values for global_accuracy and padding
        self.mesh_data = MeshData(mesh_accuracy, 1.5, mu_0, self.core.core_inner_diameter, self.core.window_w, self.windings)

    def set_core(self, core: Core):
        """Adds the core to the model

        :param core: Core object
        :type core: Core
        """
        self.core = core

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   -  -  -  -  -  -  -  -  -  -  -
    # Pre-Processing

    def create_model(self, freq: float, skin_mesh_factor: float = 0.5, pre_visualize_geometry: bool = False,
                     save_png: bool = False, color_scheme: Dict = ff.colors_femmt_default,
                     colors_geometry: Dict = ff.colors_geometry_femmt_default, benchmark: bool = False):
        """
        Create a model from the abstract geometry description inside onelab including optional mesh generation

        :param freq: Frequency [Hz]
        :type freq: float
        :param skin_mesh_factor: [default to 0.5]
        :type skin_mesh_factor: float
        :param pre_visualize_geometry: True for a pre-visualisation (e.g. check your geometry) and after this a simulation runs, False for a direct simulation
        :type pre_visualize_geometry: bool
        :param save_png: True to save a png-figure, false for no figure
        :type save_png: bool
        :param color_scheme: color file (definition for red, green, blue, ...)
        :type color_scheme: Dict
        :param colors_geometry: definition for e.g. core is grey, winding is orange, ...
        :type colors_geometry: Dict
        """
        if self.core is None:
            raise Exception("A core class needs to be added to the magnetic component")
        if self.air_gaps is None:
            self.air_gaps = AirGaps(None, None)
            self.femmt_print("No air gaps are added")
        if self.insulation is None:
            raise Exception("An insulation class need to be added to the magnetic component")
        if self.winding_windows is None:
            raise Exception("Winding windows are not set properly. Please check the winding creation")

        if benchmark:
            start_time = time.time()
            self.high_level_geo_gen(frequency=freq, skin_mesh_factor=skin_mesh_factor)
            high_level_geo_gen_time = time.time() - start_time
            start_time = time.time()
            self.mesh.generate_hybrid_mesh(visualize_before=pre_visualize_geometry, save_png=save_png,
                                           color_scheme=color_scheme, colors_geometry=colors_geometry)
            generate_hybrid_mesh_time = time.time() - start_time

            return high_level_geo_gen_time, generate_hybrid_mesh_time
        else:
            self.high_level_geo_gen(frequency=freq, skin_mesh_factor=skin_mesh_factor)
            self.mesh.generate_hybrid_mesh(visualize_before=pre_visualize_geometry, save_png=save_png,
                                           color_scheme=color_scheme, colors_geometry=colors_geometry)

    def get_single_complex_permeability(self):
        """
        Function returns the complex permeability.
        In case of amplitude dependent material definition, the initial permeability is used.
        :return: complex
        """
        if self.core.permeability_type == PermeabilityType.FromData:
            # take datasheet value from database
            complex_permeability = mu_0 * mdb.MaterialDatabase(
                self.verbosity == Verbosity.Silent).get_material_attribute(material_name=self.core.material,
                                                                           attribute="initial_permeability")
            self.femmt_print(f"{complex_permeability = }")
        if self.core.permeability_type == PermeabilityType.FixedLossAngle:
            complex_permeability = mu_0 * self.core.mu_r_abs * complex(np.cos(np.deg2rad(self.core.phi_mu_deg)),
                                                                       np.sin(np.deg2rad(self.core.phi_mu_deg)))
        if self.core.permeability_type == PermeabilityType.RealValue:
            complex_permeability = mu_0 * self.core.mu_r_abs
        return complex_permeability

    def check_model_mqs_condition(self) -> None:
        """
        Check the model for magneto-quasi-static condition for frequencies != 0

        Is called before a simulation.
        Loads the permittivity from the material database (measurement or datasheet) and calculates the
        resonance ratio = diameter_to_wavelength_ratio / diameter_to_wavelength_ratio_of_first_resonance

        :return: None
        """
        if self.frequency != 0:
            if self.core.permittivity["datasource"] == "measurements" or self.core.permittivity[
                "datasource"] == "datasheet":
                epsilon_r, epsilon_phi_deg = mdb.MaterialDatabase(self.verbosity == Verbosity.Silent).get_permittivity(
                    temperature=self.core.temperature, frequency=self.frequency,
                    material_name=self.core.material,
                    datasource=self.core.permittivity["datasource"],
                    datatype=self.core.permittivity["datatype"],
                    measurement_setup=self.core.permittivity["measurement_setup"],
                    interpolation_type="linear")

                complex_permittivity = epsilon_0 * epsilon_r * complex(np.cos(np.deg2rad(epsilon_phi_deg)),
                                                                       np.sin(np.deg2rad(epsilon_phi_deg)))
                self.femmt_print(f"{complex_permittivity = }\n"
                                 f"{epsilon_r = }\n"
                                 f"{epsilon_phi_deg = }")

                ff.check_mqs_condition(radius=self.core.core_inner_diameter / 2, frequency=self.frequency,
                                       complex_permeability=self.get_single_complex_permeability(),
                                       complex_permittivity=complex_permittivity, conductivity=self.core.sigma,
                                       relative_margin_to_first_resonance=0.5, silent=self.silent)

            else:
                ff.check_mqs_condition(radius=self.core.core_inner_diameter / 2, frequency=self.frequency,
                                       complex_permeability=self.get_single_complex_permeability(),
                                       complex_permittivity=0, conductivity=self.core.sigma,
                                       relative_margin_to_first_resonance=0.5, silent=self.silent)

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   -  -  -  -  -  -  -  -  -  -  -
    # Miscellaneous
    def calculate_core_volume_with_air(self) -> float:
        """Calculates the volume of the core including air.

        :return: Volume of the core
        :rtype: float
        """
        if self.core.core_type == CoreType.Single:
            core_height = self.core.window_h + self.core.core_inner_diameter / 2
        elif self.core.core_type == CoreType.Stacked:
            core_height = self.core.window_h_bot + self.core.window_h_top + self.core.core_inner_diameter * 3 / 4  # TODO: could also be done arbitrarily

        core_width = self.core.r_outer

        return np.pi * core_width ** 2 * core_height

    def calculate_core_volume(self) -> float:
        """Calculates the volume of the core excluding air.

        :return: Volume of the core.
        :rtype: float
        """
        core_height = None
        winding_height = None
        if self.core.core_type == CoreType.Single:
            core_height = self.core.window_h + self.core.core_inner_diameter / 2
            winding_height = self.core.window_h
        elif self.core.core_type == CoreType.Stacked:
             core_height = self.core.window_h_bot + self.core.window_h_top + self.core.core_inner_diameter * 3 / 4  # TODO: could also be done arbitrarily
             winding_height = self.core.window_h_bot + self.core.window_h_top  # TODO: could also be done arbitrarily

        core_width = self.core.r_outer

        winding_width = self.core.window_w

        air_gap_volume = 0
        inner_leg_width = self.core.r_inner - winding_width


        for leg_position, position_value, height in self.air_gaps.midpoints:
            width = 0

            if leg_position == AirGapLegPosition.LeftLeg.value:
                    # left leg
                    # TODO this is wrong since the air gap is not centered on the y axis
                    width = core_width - self.core.r_inner
            elif leg_position == AirGapLegPosition.CenterLeg.value:
                # center leg
                width = inner_leg_width
            elif leg_position == AirGapLegPosition.RightLeg.value:
                # right leg
                # TODO this is wrong since the air gap is not centered on the y axis
                width = core_width - self.core.r_inner
            else:
                raise Exception(f"Invalid leg position tag {leg_position} used for an air gap.")

            air_gap_volume += np.pi * width ** 2 * height

        return np.pi * (core_width ** 2 * core_height - (
                inner_leg_width + winding_width) ** 2 * winding_height + inner_leg_width ** 2 * winding_height) - air_gap_volume

    def calculate_core_parts_volume(self) -> list:

        """Calculates the volume of the part core excluding air.

                        :return: Volume of the core part.
                        :rtype: list
                        """
        # Extract heights from the midpoints of air gaps
        heights = [point[2] for point in self.air_gaps.midpoints]
        core_parts_volumes = []

        def get_width(part_number):
            """
            If there is a stray path, calculate width based on its starting index and part number.
            part_number is the core_part_i+2; means that if the start_index is 0, the stray path is in core_part_2
            if the start_index is 1, the stray path is in core_part_3 and so on

            """

            if self.stray_path and part_number == self.stray_path.start_index + 2:
                return self.stray_path.length
            return self.core.core_inner_diameter / 2

        # # Sorting air gaps from lower to upper
        sorted_midpoints = sorted(self.air_gaps.midpoints, key=lambda x: x[1])
        # Finding position of first airgap
        bottommost_airgap_position = sorted_midpoints[0][1]
        bottommost_airgap_height = sorted_midpoints[0][2]
        # Finding position of last airgap
        topmost_airgap_position = sorted_midpoints[-1][1]
        topmost_airgap_height = sorted_midpoints[-1][2]

        # if single core
        if self.core.core_type == CoreType.Single:

            # For single core and more than one core_part, volume for every core part is calculated


            # core_part_1 is divided into subparts cores
            # # subpart1: bottom left subpart
            subpart1_1_height = bottommost_airgap_position + self.core.window_h / 2 - bottommost_airgap_height / 2
            subpart1_1_width = self.core.core_inner_diameter / 2
            subpart1_1_volume = np.pi * subpart1_1_width ** 2 * subpart1_1_height

            # # subpart2: bottom mid subpart
            subpart1_2_height = self.core.core_inner_diameter / 4
            subpart1_2_width = self.core.r_outer
            subpart1_2_volume = np.pi * subpart1_2_width ** 2 * subpart1_2_height


            # subpart3: right subpart
            subpart1_3_height = self.core.window_h
            subpart1_3_width = self.core.r_outer
            subpart1_3_volume = np.pi * subpart1_3_width ** 2 * subpart1_3_height - (
                        np.pi * (self.core.window_w + self.core.core_inner_diameter / 2) ** 2 * self.core.window_h)

            # subpart4: top mid subpart
            subpart1_4_height = self.core.core_inner_diameter / 4
            subpart1_4_width = self.core.r_outer
            subpart1_4_volume = np.pi * subpart1_4_width ** 2 * subpart1_4_height

            # subpart5: top left subpart
            subpart1_5_height = self.core.window_h / 2 - topmost_airgap_position - topmost_airgap_height / 2
            subpart1_5_width = self.core.core_inner_diameter / 2
            subpart1_5_volume = np.pi * subpart1_5_width ** 2 * subpart1_5_height

            # Calculate the volume of core part 1 by summing up subpart volumes
            core_part_1_volume = subpart1_1_volume + subpart1_2_volume + subpart1_3_volume + subpart1_4_volume + subpart1_5_volume
            core_parts_volumes.append(core_part_1_volume)

            # Calculate the volumes of the core parts between the air gaps
            for i in range(len(sorted_midpoints) - 1):
                air_gap_1_position = sorted_midpoints[i][1]
                air_gap_1_height = sorted_midpoints[i][2]
                air_gap_2_position = sorted_midpoints[i + 1][1]
                air_gap_2_height = sorted_midpoints[i + 1][2]
                # calculate the height based on airgap positions and heights, and the width
                core_part_height = air_gap_2_position - air_gap_2_height / 2 - (
                        air_gap_1_position + air_gap_1_height / 2)
                core_part_width = get_width(i + 2)
                # calculate the volume
                core_part_volume = np.pi * core_part_width ** 2 * core_part_height
                core_parts_volumes.append(core_part_volume)

            # Return the total core part volume
            #return core_parts_volumes

        elif self.core.core_type == CoreType.Stacked:

            # For stacked core types, the volume is divided into different core  * parts, each of which is further
            # divided into subparts to calculate the total volume of each core part.

            # core_part_2 : core part between the bottom airgap and subpart_1 of core_part_1
            core_part_1_height = self.core.window_h_bot / 2 - heights[0] / 2
            core_part_1_width = self.core.core_inner_diameter / 2
            core_part_1_volume = np.pi * core_part_1_width ** 2 * core_part_1_height
            core_parts_volumes.append(core_part_1_volume)

            # Core Part 1 Calculation
            # Core part 1 is calculated as the sum of three different subparts
            # subpart_1: bottom left subpart
            subpart2_1_height = self.core.window_h_bot / 2 - heights[0] / 2
            subpart2_1_width = self.core.core_inner_diameter / 2
            subpart2_1_volume = np.pi * subpart2_1_width ** 2 * subpart2_1_height

            # subpart_2 : bottom mid subpart
            subpart2_2_height = self.core.core_inner_diameter / 4
            subpart2_2_width = self.core.r_outer
            subpart2_2_volume = np.pi * subpart2_2_width ** 2 * subpart2_2_height

            # subpart_3: bottom right subpart
            subpart2_3_height = self.core.window_h_bot
            subpart2_3_width = self.core.r_outer
            subpart2_3_volume = np.pi * (subpart2_3_width ** 2 * subpart2_3_height) - np.pi * ((self.core.window_w + self.core.core_inner_diameter/2) ** 2 * self.core.window_h_bot)


            # Summing up the volumes of the subparts to get the total volume of core part 1
            core_part_2_volume = subpart2_1_volume + subpart2_2_volume + subpart2_3_volume
            core_parts_volumes.append(core_part_2_volume)

            # core_part_3 : left mid core part (stacked)
            core_part_3_height = self.core.core_inner_diameter / 4
            core_part_3_width = self.core.r_inner
            core_part_3_volume = np.pi * core_part_3_width ** 2 * core_part_3_height
            core_parts_volumes.append(core_part_3_volume)

            # core_part_4: right mid core part
            core_part_4_height = self.core.core_inner_diameter / 4
            core_part_4_width = self.core.r_outer
            core_part_4_volume = np.pi * core_part_4_width ** 2 * core_part_4_height - core_part_3_volume
            core_parts_volumes.append(core_part_4_volume)

            # core_part_5
            # core_part_5 is divided into 3 parts
            # subpart_1: left top subpart
            subpart5_1_height = self.core.window_h_top - heights[1] / 2
            subpart5_1_width = self.core.core_inner_diameter / 2
            subpart5_1_volume = np.pi * subpart5_1_width ** 2 * subpart5_1_height

            # subpart_2: mid top subpart
            subpart5_2_height = self.core.core_inner_diameter / 4
            subpart5_2_width = self.core.r_outer
            subpart5_2_volume = np.pi * subpart5_2_width ** 2 * subpart5_2_height

            # subpart 3: top right subpart
            subpart5_3_height = self.core.window_h_top
            subpart5_3_width = self.core.r_outer
            subpart5_3_volume = np.pi * (subpart5_3_width ** 2 * subpart5_3_height) - np.pi * ((self.core.window_w + self.core.core_inner_diameter / 2) ** 2 * self.core.window_h_top)

            # Summing up the volumes of the subparts to get the total volume of core_part_5
            core_part_5_volume = subpart5_1_volume + subpart5_2_volume + subpart5_3_volume
            core_parts_volumes.append(core_part_5_volume)

        # Core Volume Consistency Check
        # Sum all the core part volumes
        total_parts_volume = sum(core_parts_volumes)

        # Calculate the whole core volume
        whole_core_volume = self.calculate_core_volume()

        # Define a margin of error
        margin_of_error = 1e-5

        # Check if the volumes are equal within the margin of error
        if not (abs(whole_core_volume - total_parts_volume) <= margin_of_error):
            error_message = (f"Sum of core parts ({total_parts_volume}) does not equal "
                             f"the whole core volume ({whole_core_volume}) within the margin of error ({margin_of_error}).")
            raise ValueError(error_message)

        # Returning the final list of core part volumes
        return core_parts_volumes

    def calculate_core_weight(self) -> float:
        """
        Calculates the weight of the core in kg.
        This method is using the core volume from for an ideal rotation-symmetric core and the volumetric mass density from the material database.
        """
        if self.core.material == 'custom':
            volumetric_mass_density = 0
            warnings.warn(
                "Volumetric mass density not implemented for custom cores. Returns '0' in log-file: Core cost will also result to 0.")
        else:
            volumetric_mass_density = self.core.material_database.get_material_attribute(
                material_name=self.core.material, attribute="volumetric_mass_density")
        return self.calculate_core_volume() * volumetric_mass_density




    def get_wire_distances(self) -> List[List[float]]:
        """Helper function which returns the distance (radius) of each conductor to the y-axis

        :return: Wire distances
        :rtype: List[List[float]]
        """
        # wire_distance = []
        # for winding in self.two_d_axi.p_conductor:
        #     # 5 points are for 1 wire
        #     num_points = len(winding)
        #     num_windings = num_points // 5
        #     winding_list = []
        #     for i in range(num_windings):
        #         winding_list.append(winding[i * 5][0])
        #     wire_distance.append(winding_list)
        #
        # return wire_distance

        wire_distance = []
        for num, conductor in enumerate(self.two_d_axi.p_conductor):
            # If the conductor is of type RectangularSolid, it is represented by 4 points (the corners of the rectangle)
            if self.windings[num].conductor_type == ConductorType.RectangularSolid:
                num_points = len(conductor)
                num_turns = num_points // 4
                point_increment = 4
            else:
                num_points = len(conductor)
                num_turns = num_points // 5
                point_increment = 5

            winding_list = []
            for i in range(num_turns):
                winding_list.append(conductor[i * point_increment][0])
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
        wire_lengths = self.calculate_wire_lengths()
        for index, winding in enumerate(self.windings):
            cross_section_area = 0
            if winding.conductor_type == ConductorType.RoundLitz or winding.conductor_type == ConductorType.RoundSolid:
                # For round wire it is always the same
                cross_section_area = np.pi * winding.conductor_radius ** 2
            elif winding.conductor_type == ConductorType.RectangularSolid:
                # Since the foil sizes also depends on the winding scheme, conductor_arrangement and wrap_para_type
                # the volume calculation is different.
                for ww in self.winding_windows:
                    for vww_index, vww in enumerate(ww.virtual_winding_windows):
                        winding_type = vww.winding_type
                        winding_scheme = vww.winding_scheme
                        wrap_para_type = vww.wrap_para
                        for vww_winding in vww.windings:
                            if vww_winding.winding_number == index:
                                if winding_type == WindingType.Single:
                                    if winding_scheme == WindingScheme.Full:
                                        cross_section_area = self.core.window_h * self.core.window_w
                                    elif winding_scheme == WindingScheme.FoilHorizontal:
                                        cross_section_area = self.core.window_w * winding.thickness
                                    elif winding_scheme == WindingScheme.FoilVertical:
                                        if wrap_para_type == WrapParaType.FixedThickness:
                                            cross_section_area = self.core.window_h * winding.thickness
                                        elif wrap_para_type == WrapParaType.Interpolate:
                                            cross_section_area = self.core.window_h * self.core.window_w / vww.turns[
                                                vww_index]
                                        else:
                                            raise Exception(f"Unknown wrap para type {wrap_para_type}")
                                    else:
                                        raise Exception(f"Unknown winding scheme {winding_scheme}")
                                elif winding_type == WindingType.TwoInterleaved:
                                    # Since interleaved winding type currently only supports round conductors this can be left empty.
                                    pass
                                elif winding_type == WindingType.CenterTappedGroup:
                                    cross_section_area = self.core.window_w * winding.thickness
                                else:
                                    raise Exception(f"Unknown winding type {winding_type}")
            else:
                raise Exception(f"Unknown conductor type {winding.conductor_type}")

            wire_volumes.append(cross_section_area * wire_lengths[index])

        return wire_volumes

    def calculate_wire_weight(self) -> List[float]:
        wire_material = ff.wire_material_database()

        wire_weight = []

        # TODO: distinguish between wire material. Only copper at the moment
        for wire_volume in self.calculate_wire_volumes():
            wire_weight.append(wire_volume * wire_material["Copper"].volumetric_mass_density)

        return wire_weight

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    # GetDP Interaction / Simulation / Excitation
    def excitation(self, frequency: float, amplitude_list: List, phase_deg_list: List = None, ex_type: str = 'current',
                   plot_interpolation: bool = False, imposed_red_f=0):
        """
        - excitation of the electromagnetic problem
        - current, voltage or current density
        - frequency or reduced frequency

        :param plot_interpolation:
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

        # negative currents are not allowed and lead to wrong simulation results. Check for this.
        # this message appears after meshing but before simulation
        for amplitude in amplitude_list:
            if amplitude < 0:
                raise ValueError(
                    "Negative currents are not allowed. Use the phase + 180 degree to generate a negative current.")

        self.femmt_print(f"\n---\n"
                         f"Excitation: \n"
                         f"Frequency: {frequency}\n"
                         f"Current(s): {amplitude_list}\n"
                         f"Phase(s): {phase_deg_list}\n")

        # -- Excitation --
        self.flag_excitation_type = ex_type  # 'current', 'current_density', 'voltage'
        if self.core.permeability["datasource"] != MaterialDataSource.Custom:
            self.core.update_core_material_pro_file(frequency, self.file_data.electro_magnetic_folder_path,
                                                    plot_interpolation)  # frequency update to core class
        if self.core.permittivity["datasource"] != MaterialDataSource.Custom:
            self.core.update_sigma(frequency)
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
                        raise ValueError("Missing phases inside excitation, e.g. 'phase_deg_list = [0, 180]'. ")
                else:
                    self.phase_deg = phase_deg_list
                    # Define complex current phasor as excitation
                    self.current[num] = complex(amplitude_list[num] * np.cos(np.deg2rad(phase_deg_list[num])),
                                                amplitude_list[num] * np.sin(np.deg2rad(phase_deg_list[num])))

        # Imposed current density
        if self.flag_excitation_type == 'current_density':
            raise NotImplementedError

        # Imposed voltage
        if self.flag_excitation_type == 'voltage':
            raise NotImplementedError

        # -- Frequency --

        self.frequency = frequency  # in Hz

        # Define reduced frequency (used for homogenization technique)
        # self.red_freq = np.empty(2)
        self.red_freq = []
        for num in range(len(self.windings)):
            self.red_freq.append([])

        if self.frequency != 0:
            self.delta = np.sqrt(2 / (2 * self.frequency * np.pi * self.windings[
                0].cond_sigma * mu_0))  # TODO: distinguish between material conductivities
            for num in range(len(self.windings)):
                if self.windings[num].conductor_type == ConductorType.RoundLitz:
                    self.red_freq[num] = self.windings[num].strand_radius / self.delta
                elif self.windings[num].conductor_type == ConductorType.RoundSolid:
                    self.red_freq[num] = self.windings[num].conductor_radius / self.delta
                else:
                    self.femmt_print("Reduced Frequency does not have a physical value here")
                    self.femmt_print(self.windings[num].conductor_type)
                    self.red_freq[
                        num] = 1  # TODO: doesn't make sense like this -> rewrite fore conductor windings shape
        else:
            # DC case
            self.delta = 1e20  # random huge value
            for num in range(len(self.windings)):
                self.red_freq[num] = 0

    def simulate(self):
        """
        Initializes an onelab client. Provides the GetDP based solver with the created mesh file.
        """
        self.femmt_print(f"\n---\n"
                         f"Initialize ONELAB API\n"
                         f"Run Simulation\n")

        # -- Simulation --
        # create a new onelab client

        # Initial Clearing of gmsh data
        gmsh.clear()

        # get model file names with correct path
        solver = os.path.join(self.file_data.electro_magnetic_folder_path, "ind_axi_python_controlled.pro")

        os.chdir(self.file_data.working_directory)

        if self.verbosity == Verbosity.Silent:
            verbose = "-verbose 1"
        else:
            verbose = "-verbose 5"

        to_file_str = ""
        if self.verbosity == Verbosity.ToFile:
            to_file_str = " > " + self.file_data.getdp_log

        # Run simulations as sub clients (non-blocking??)
        getdp_filepath = os.path.join(self.file_data.onelab_folder_path, "getdp")
        self.onelab_client.runSubClient("myGetDP",
                                        getdp_filepath + " " + solver + " -msh " + self.file_data.e_m_mesh_file + " -solve Analysis -v2 " + verbose + to_file_str)

    def write_simulation_parameters_to_pro_files(self):
        """
        Interaction between python and Prolog files.

        Writes the simulation parameters to .pro-files

        Parameter.pro: includes material properties, currents, phases, ...
        postquantities.pro: includes directions to store the raw results from the FEM simulation

        """
        # All shared control variables and parameters are passed to a temporary Prolog file
        self.femmt_print(f"\n---\n"
                         f"Write simulation parameters to .pro files (file communication).\n")

        # Write initialization parameters for simulation in 'Parameter.pro' file
        self.write_electro_magnetic_parameter_pro()

        # Write postprocessing parameters in 'postquantities.pro' file
        self.write_electro_magnetic_post_pro()

    def overwrite_conductors_with_air(self, physical_surfaces_to_overwrite: list):
        """

        EXPERIMENTAL



        :return:
        """
        if True:
            with open(os.path.join(os.path.join(self.file_data.e_m_mesh_file)), "r") as mesh_file:
                mesh_data = mesh_file.read()

            for ps in physical_surfaces_to_overwrite:
                # mesh_data = mesh_data.replace(f'1 {ps} 4', f'1 {ps+1000000} 4')
                mesh_data = mesh_data.replace(f'1 {ps} 4', f'1 {ps + 1000000} 4')

            with open(os.path.join(os.path.join(self.file_data.e_m_mesh_file)), "w") as mesh_file:
                mesh_file.write(mesh_data)

    def overwrite_air_conductors_with_conductors(self, physical_surfaces_to_overwrite: list):
        """

        EXPERIMENTAL



        :return:
        """
        if True:
            with open(os.path.join(os.path.join(self.file_data.e_m_mesh_file)), "r") as mesh_file:
                mesh_data = mesh_file.read()

            for ps in physical_surfaces_to_overwrite:
                mesh_data = mesh_data.replace(f'1 {ps} 4', f'1 {ps - 1000000} 4')

            with open(os.path.join(os.path.join(self.file_data.e_m_mesh_file)), "w") as mesh_file:
                mesh_file.write(mesh_data)

    def single_simulation(self, freq: float, current: List[float], phi_deg: List[float] = None,
                          plot_interpolation: bool = False, show_fem_simulation_results: bool = True,
                          benchmark: bool = False):
        """
        Start a _single_ electromagnetic ONELAB simulation.

        :param plot_interpolation:
        :param freq: frequency to simulate
        :type freq: float
        :param current: current to simulate
        :param phi_deg: phase angle in degree
        :type phi_deg: List[float]
        :param show_fem_simulation_results: Set to True to show the simulation results after the simulation has finished
        :type show_fem_simulation_results: bool
        """
        # negative currents are not allowed and lead to wrong simulation results. Check for this.
        # this message appears before meshing and before simulation
        # there is another ValueError rising inside excitation()-method for safety (but after meshing).
        if type(current) is not list:
            raise Exception("The current must be given in a list.")
        for current_value in current:
            if current_value < 0:
                raise ValueError(
                    "Negative currents are not allowed. Use the phase + 180 degree to generate a negative current.")

        phi_deg = phi_deg or []
        if benchmark:
            start_time = time.time()
            self.mesh.generate_electro_magnetic_mesh()
            generate_electro_magnetic_mesh_time = time.time() - start_time

            start_time = time.time()
            self.excitation(frequency=freq, amplitude_list=current, phase_deg_list=phi_deg,
                            plot_interpolation=plot_interpolation)  # frequency and current
            self.check_model_mqs_condition()
            self.write_simulation_parameters_to_pro_files()
            self.generate_load_litz_approximation_parameters()

            prepare_simulation_time = time.time() - start_time

            start_time = time.time()
            self.simulate()
            real_simulation_time = time.time() - start_time

            start_time = time.time()
            self.calculate_and_write_log()  # TODO: reuse center tapped
            logging_time = time.time() - start_time
            if show_fem_simulation_results:
                self.visualize()

            return generate_electro_magnetic_mesh_time, prepare_simulation_time, real_simulation_time, logging_time
        else:
            self.mesh.generate_electro_magnetic_mesh()
            self.excitation(frequency=freq, amplitude_list=current, phase_deg_list=phi_deg,
                            plot_interpolation=plot_interpolation)  # frequency and current
            self.check_model_mqs_condition()
            self.write_simulation_parameters_to_pro_files()
            self.generate_load_litz_approximation_parameters()
            self.simulate()
            self.calculate_and_write_log()  # TODO: reuse center tapped
            if show_fem_simulation_results:
                self.visualize()

    def excitation_sweep(self, frequency_list: List, current_list_list: List, phi_deg_list_list: List,
                         show_last_fem_simulation: bool = False,
                         excitation_meshing_type: ExcitationMeshingType = None, skin_mesh_factor: float = 0.5,
                         visualize_before: bool = False, save_png: bool = False,
                         color_scheme: Dict = ff.colors_femmt_default,
                         colors_geometry: Dict = ff.colors_geometry_femmt_default,
                         inductance_dict: Dict = None, core_hyst_loss: float = None) -> None:
        """
        Performs a sweep simulation for frequency-current pairs. Both values can
        be passed in lists of the same length. The mesh is only created ones (fast sweep)!

        :Example Code for Inductor:

        >>> import femmt as fmt
        >>> fs_list = [0, 10000, 30000, 60000, 100000, 150000]
        >>> amplitude_list_list = [[10], [2], [1], [0.5], [0.2], [0.1]]
        >>> phase_list_list = [[0], [10], [20], [30], [40], [50]]
        >>> geo.excitation_sweep(frequency_list=fs_list, current_list_list=amplitude_list_list, phi_deg_list_list=phase_list_list)

        :Example Code for Transformer with 2 windings:

        >>> import femmt as fmt
        >>> fs_list = [0, 10000, 30000, 60000, 100000, 150000]
        >>> amplitude_list_list = [[10, 2], [2, 1], [1, 0.5], [0.5, 0.25], [0.2, 0.1], [0.1, 0.05]]
        >>> phase_list_list = [[0, 170], [10, 180], [20, 190], [30, 200], [40, 210], [50, 220]]
        >>> geo.excitation_sweep(frequency_list=fs_list, current_list_list=amplitude_list_list, phi_deg_list_list=phase_list_list)

        :param frequency_list: Frequency in a list
        :type frequency_list: List
        :param current_list_list: current amplitude, must be a list in a list, see example!
        :type current_list_list: List
        :param phi_deg_list_list: phase in degree, must be a list in a list, see example!
        :type phi_deg_list_list: List
        :param show_last_fem_simulation: shows last simulation in gmsh if set to True
        :type show_last_fem_simulation: bool
        :param visualize_before: show genarated mesh before the simulation is run
        :type visualize_before: bool
        :param color_scheme: colorfile (definition for red, green, blue, ...)
        :type color_scheme: Dict
        :param colors_geometry: definition for e.g. core is grey, winding is orange, ...
        :type colors_geometry: Dict
        :param save_png: True to save a .png
        :type save_png: bool

        """
        # negative currents are not allowed and lead to wrong simulation results. Check for this.
        # this message appears before meshing and before simulation
        # there is another ValueError rising inside excitation()-method for safety (but after meshing).
        for current_list in current_list_list:
            for current in current_list:
                if current < 0:
                    raise ValueError(
                        "Negative currents are not allowed. Use the phase + 180 degree to generate a negative current.")

        # frequencies = frequencies or []
        # currents = currents or []
        # phi = phi or []
        if show_last_fem_simulation:
            self.plot_fields = "standard"
        else:
            self.plot_fields = False

        # If one conductor is solid and no meshing type is given then change the meshing type to MeshEachFrequency
        # In case of litz wire, only the lowest frequency is meshed (frequency indecent due to litz-approximation)
        if excitation_meshing_type is None:
            for winding in self.windings:
                if winding.conductor_type == ConductorType.RoundSolid:
                    excitation_meshing_type = ExcitationMeshingType.MeshEachFrequency
                    break
                if winding.conductor_type == ConductorType.RoundLitz:
                    excitation_meshing_type = ExcitationMeshingType.MeshOnlyLowestFrequency

        if excitation_meshing_type == ExcitationMeshingType.MeshEachFrequency:
            for count_frequency, _ in enumerate(frequency_list):
                self.high_level_geo_gen(frequency=frequency_list[count_frequency], skin_mesh_factor=skin_mesh_factor)
                self.mesh.generate_hybrid_mesh(color_scheme, colors_geometry, visualize_before=visualize_before,
                                               save_png=save_png)
                self.mesh.generate_electro_magnetic_mesh()

                self.excitation(frequency=frequency_list[count_frequency],
                                amplitude_list=current_list_list[count_frequency],
                                phase_deg_list=phi_deg_list_list[count_frequency])  # frequency and current
                if count_frequency == 0: self.check_model_mqs_condition()
                self.write_simulation_parameters_to_pro_files()
                self.generate_load_litz_approximation_parameters()
                self.simulate()
        else:
            if excitation_meshing_type == ExcitationMeshingType.MeshOnlyHighestFrequency:
                self.high_level_geo_gen(frequency=max(frequency_list), skin_mesh_factor=skin_mesh_factor)
            elif excitation_meshing_type == ExcitationMeshingType.MeshOnlyLowestFrequency:
                self.high_level_geo_gen(frequency=min(frequency_list), skin_mesh_factor=skin_mesh_factor)
            else:
                raise Exception(f"Unknown excitation meshing type {excitation_meshing_type}")
            self.mesh.generate_hybrid_mesh(color_scheme, colors_geometry, visualize_before=visualize_before,
                                           save_png=save_png)
            self.mesh.generate_electro_magnetic_mesh()

            check_model_mqs_condition_already_performerd = False
            for count_frequency, value_frequency in enumerate(range(0, len(frequency_list))):
                self.excitation(frequency=frequency_list[count_frequency],
                                amplitude_list=current_list_list[count_frequency],
                                phase_deg_list=phi_deg_list_list[count_frequency])  # frequency and current
                if value_frequency != 0 and not check_model_mqs_condition_already_performerd:
                    self.check_model_mqs_condition()
                    check_model_mqs_condition_already_performerd = True
                self.write_simulation_parameters_to_pro_files()
                self.generate_load_litz_approximation_parameters()
                self.simulate()
                # self.visualize()

        self.calculate_and_write_log(sweep_number=len(frequency_list), currents=current_list_list,
                                     frequencies=frequency_list, inductance_dict=inductance_dict,
                                     core_hyst_losses=core_hyst_loss)

        if show_last_fem_simulation:
            self.write_simulation_parameters_to_pro_files()
            self.visualize()

    def component_study(self, time_current_vectors: List[List[List[float]]], fft_filter_value_factor: float = 0.01):
        """
        :param fft_filter_value_factor: Factor to filter frequencies from the fft. E.g. 0.01 [default] removes all amplitudes below 1 % of the maximum amplitude from the result-frequency list
        :type fft_filter_value_factor: float
        """

        # winding losses
        frequency_current_phase_deg_list = []

        # collect simulation input parameters from time_current_vectors
        hyst_loss_amplitudes = []
        hyst_loss_phases_deg = []
        hyst_frequency = 1 / (time_current_vectors[0][0][-1])
        for time_current_vector in time_current_vectors:
            # collect winding losses simulation input parameters
            [frequency_list, amplitude, phi_rad] = ff.fft(time_current_vector, mode='time', filter_value_factor=fft_filter_value_factor)
            phi_deg = np.rad2deg(phi_rad)
            frequency_current_phase_deg_list.append([frequency_list, amplitude, phi_deg])

            # collect hysteresis loss simulation input parameters
            hyst_loss_amplitudes.append(fr.max_value_from_value_vec(time_current_vector[1])[0])
            hyst_loss_phases_deg.append(
                fr.phases_deg_from_time_current(time_current_vector[0], time_current_vector[1])[0])

        # check if all frequency vectors include the same frequencies
        for count in range(len(frequency_current_phase_deg_list) - 1):
            if not np.array_equal(frequency_current_phase_deg_list[count][0],
                                  frequency_current_phase_deg_list[count + 1][0]):
                raise ValueError("Frequency vectors for different currents are not the same!")

        # transfer format from fft()-output to excitation_sweep()-input
        current_list_list = []
        phi_deg_list_list = []
        for count_frequency, frequency in enumerate(frequency_list):
            currents_single_frequency = []
            phi_deg_single_frequency = []
            for count_current, _ in enumerate(time_current_vectors):
                currents_single_frequency.append(frequency_current_phase_deg_list[count_current][1][count_frequency])
                phi_deg_single_frequency.append(frequency_current_phase_deg_list[count_current][2][count_frequency])
            current_list_list.append(currents_single_frequency)
            phi_deg_list_list.append(phi_deg_single_frequency)

        # get the inductance
        inductance_dict = self.get_inductances(I0=1, op_frequency=hyst_frequency, skin_mesh_factor=1)

        # calculate hysteresis losses
        # use a single simulation
        self.generate_load_litz_approximation_parameters()
        self.excitation(frequency=hyst_frequency, amplitude_list=hyst_loss_amplitudes, phase_deg_list=hyst_loss_phases_deg, plot_interpolation=False)  # frequency and current
        self.check_model_mqs_condition()
        self.write_simulation_parameters_to_pro_files()
        self.generate_load_litz_approximation_parameters()
        self.simulate()
        self.calculate_and_write_log()  # TODO: reuse center tapped
        [p_hyst] = self.load_result(res_name="p_hyst")

        # calculate the winding losses
        self.excitation_sweep(frequency_list, current_list_list, phi_deg_list_list, inductance_dict=inductance_dict,
                              core_hyst_loss=float(p_hyst))

    def center_tapped_pre_study(self, time_current_vectors: List[List[List[float]]], plot_waveforms: bool = False,
                                fft_filter_value_factor: float = 0.01) -> Dict:
        """
        As magnetizing currents are often non-sinusoidal, some corrections in the simulation current waveforms
        are needed. This function calculates the new current waveforms for the center tapped study to get
        inductance values and so on.

        :param time_current_vectors: time-current vectors for primary and secondary
        :type time_current_vectors: List[List[List[float]]]
        :param plot_waveforms: True to watch the pre-calculated waveforms
        :type plot_waveforms: bool
        :param fft_filter_value_factor: Factor to filter frequencies from the fft. E.g. 0.01 [default] removes all amplitudes below 1 % of the maximum amplitude from the result-frequency list
        :type fft_filter_value_factor: float
        :return: new current waveform vector
        :rtype: Dict

        return dict:
        center_tapped_study_excitation = {
            "hysteresis": {
                "frequency": None,
                "transformer": {
                    "current_amplitudes": None,
                    "current_phases_deg": None
                },
                "choke": {
                    "current_amplitudes": None,
                    "current_phases_deg": None
                }
            },
            "linear_losses": {
                "frequencies": None,
                "current_amplitudes": None,
                "current_phases_deg": None
            }
        }
        """

        def hysteresis_loss_excitation(input_time_current_vectors):
            # collect simulation input parameters from time_current_vectors
            hyst_loss_amplitudes = []
            hyst_loss_phases_deg = []
            hyst_frequency = 1 / (input_time_current_vectors[0][0][-1])
            for time_current_vector in input_time_current_vectors:
                # collect hysteresis loss simulation input parameters
                hyst_loss_amplitudes.append(fr.max_value_from_value_vec(time_current_vector[1])[0])
                hyst_loss_phases_deg.append(
                    fr.phases_deg_from_time_current(time_current_vector[0], time_current_vector[1])[0])
            return hyst_frequency, hyst_loss_amplitudes, hyst_loss_phases_deg

        def split_hysteresis_loss_excitation_center_tapped(hyst_frequency, hyst_loss_amplitudes, hyst_loss_phases_deg):
            # print(f"{hyst_frequency, hyst_loss_amplitudes, hyst_loss_phases_deg = }")
            hyst_loss_amplitudes[-1] = hyst_loss_amplitudes[-1] / 2
            hyst_loss_amplitudes.append(hyst_loss_amplitudes[-1])
            hyst_loss_phases_deg.append(hyst_loss_phases_deg[-1])
            # print(f"{hyst_frequency, hyst_loss_amplitudes, hyst_loss_phases_deg = }")
            return hyst_frequency, hyst_loss_amplitudes, hyst_loss_phases_deg

        def split_time_current_vectors_center_tapped(time_current_vectors: List[List[List[float]]]):
            # print(f"{time_current_vectors = }")
            positive_secondary_current = np.copy(time_current_vectors[1][1])
            positive_secondary_current[positive_secondary_current < 0] = 0
            # print(f"{positive_secondary_current = }")
            negative_secondary_current = np.copy(time_current_vectors[1][1])
            negative_secondary_current[negative_secondary_current > 0] = 0
            # print(f"{negative_secondary_current = }")

            center_tapped_time_current_vectors = [time_current_vectors[0],
                                                  [time_current_vectors[1][0], positive_secondary_current],
                                                  [time_current_vectors[1][0], negative_secondary_current]]

            if plot_waveforms:
                plt.plot(time_current_vectors[1][0], negative_secondary_current, label="negative_secondary_current")
                plt.plot(time_current_vectors[1][0], positive_secondary_current, label="positive_secondary_current")
                plt.plot(time_current_vectors[0][0], time_current_vectors[0][1], label="primary_current")
                plt.xlabel("time / s")
                plt.ylabel("current / A")
                plt.grid()
                plt.legend()
                plt.show()

            return center_tapped_time_current_vectors

        def linear_loss_excitation(time_current_vectors: List[List[List[float]]], fft_filter_value_factor: float = 0.01):
            """
            Perform FFT to get the primary and secondary currents to calculate the wire losses.
            These losses can be 'linear added' to get the total winding losses.

            :param time_current_vectors: primary and secondary current waveforms over time
            :type time_current_vectors: List[List[List[float]]]
            :param fft_filter_value_factor: Factor to filter frequencies from the fft. E.g. 0.01 [default] removes all amplitudes below 1 % of the maximum amplitude from the result-frequency list
            :type fft_filter_value_factor: float
            """
            # winding losses
            frequency_current_phase_deg_list = []
            # collect winding losses simulation input parameters
            for time_current_vector in time_current_vectors:
                [frequency_list, amplitude, phi_rad] = ff.fft(time_current_vector, mode='time', filter_value_factor=fft_filter_value_factor)
                phi_deg = np.rad2deg(phi_rad)
                frequency_current_phase_deg_list.append([frequency_list, amplitude, phi_deg])

            # check if all frequency vectors include the same frequencies
            # WORKAROUND: if any frequency is not included in one of the vectors it is added with amplitude  = 0 and phase = 0
            # TODO: recalculate the fft at the "missing frequencies and add their values...
            all_frequencies = set()
            for count in range(len(frequency_current_phase_deg_list) - 1):
                if not np.array_equal(frequency_current_phase_deg_list[count][0],
                                      frequency_current_phase_deg_list[count + 1][0]):
                    # raise ValueError("Frequency vectors for different currents are not the same!")
                    # print("Frequency vectors for different currents are not the same!")
                    all_frequencies = all_frequencies | set(frequency_current_phase_deg_list[count][0]) | set(
                        frequency_current_phase_deg_list[count + 1][0])
                    # print(f"Original: {frequency_current_phase_deg_list = }")

            for frequency in list(all_frequencies):
                for count in range(0, len(frequency_current_phase_deg_list)):
                    if frequency not in frequency_current_phase_deg_list[count][0]:
                        ii = np.searchsorted(frequency_current_phase_deg_list[count][0], frequency)
                        frequency_current_phase_deg_list[count][0] = np.insert(
                            frequency_current_phase_deg_list[count][0], ii, frequency)
                        frequency_current_phase_deg_list[count][1] = np.insert(
                            frequency_current_phase_deg_list[count][1], ii, 0)
                        frequency_current_phase_deg_list[count][2] = np.insert(
                            frequency_current_phase_deg_list[count][2], ii, 0)
            # print(f"Corrected: {frequency_current_phase_deg_list = }")
            return frequency_list, frequency_current_phase_deg_list

        center_tapped_study_excitation = {
            "hysteresis": {
                "frequency": None,
                "transformer": {
                    "current_amplitudes": None,
                    "current_phases_deg": None
                },
                "choke": {
                    "current_amplitudes": None,
                    "current_phases_deg": None
                }
            },
            "linear_losses": {
                "frequencies": None,
                "current_amplitudes": None,
                "current_phases_deg": None
            }
        }

        # Hysteresis Loss Excitation
        time_current_vectors[1][1] = time_current_vectors[1][1] * (-1)
        hyst_frequency, hyst_loss_amplitudes, hyst_loss_phases_deg = hysteresis_loss_excitation(time_current_vectors)
        hyst_frequency, hyst_loss_amplitudes, hyst_loss_phases_deg = split_hysteresis_loss_excitation_center_tapped(
            hyst_frequency, hyst_loss_amplitudes, hyst_loss_phases_deg)
        center_tapped_study_excitation["hysteresis"]["frequency"] = hyst_frequency

        if plot_waveforms:
            i_1 = hyst_loss_amplitudes[0] * np.cos(
                time_current_vectors[0][0] * 2 * np.pi * hyst_frequency - np.deg2rad(hyst_loss_phases_deg[0]))
            i_2 = hyst_loss_amplitudes[1] * np.cos(
                time_current_vectors[0][0] * 2 * np.pi * hyst_frequency - np.deg2rad(hyst_loss_phases_deg[1]))
            plt.plot(time_current_vectors[0][0], i_1, label="i_1")
            plt.plot(time_current_vectors[0][0], i_2, "-", label="i_2")
            plt.xlabel("time / s")
            plt.ylabel("current / A")
            plt.grid()
            plt.legend()
            plt.show()

        # calculate hysteresis losses in the xfmr
        xfmr_scale = 1.7
        center_tapped_study_excitation["hysteresis"]["transformer"]["current_amplitudes"] = list(
            np.array(hyst_loss_amplitudes) * xfmr_scale)
        center_tapped_study_excitation["hysteresis"]["transformer"]["current_phases_deg"] = hyst_loss_phases_deg

        # calculate hysteresis losses in the choke
        choke_hyst_loss_amplitudes = hyst_loss_amplitudes
        choke_hyst_loss_amplitudes[1] = choke_hyst_loss_amplitudes[0] * 7
        choke_hyst_loss_amplitudes[2] = choke_hyst_loss_amplitudes[0] * 7
        center_tapped_study_excitation["hysteresis"]["choke"]["current_amplitudes"] = choke_hyst_loss_amplitudes
        center_tapped_study_excitation["hysteresis"]["choke"]["current_phases_deg"] = [0, 180, 180]

        # Linear Loss Excitation
        time_current_vectors = split_time_current_vectors_center_tapped(time_current_vectors)
        frequency_list, frequency_current_phase_deg_list = linear_loss_excitation(time_current_vectors, fft_filter_value_factor)

        if plot_waveforms:
            i_1 = hyst_loss_amplitudes[0] * np.cos(
                time_current_vectors[0][0] * 2 * np.pi * hyst_frequency - np.deg2rad(hyst_loss_phases_deg[0]))
            i_2 = hyst_loss_amplitudes[1] * np.cos(
                time_current_vectors[0][0] * 2 * np.pi * hyst_frequency - np.deg2rad(hyst_loss_phases_deg[1]))
            i_3 = hyst_loss_amplitudes[2] * np.cos(
                time_current_vectors[0][0] * 2 * np.pi * hyst_frequency - np.deg2rad(hyst_loss_phases_deg[2]))
            plt.plot(time_current_vectors[0][0], i_1, label="i_1")
            plt.plot(time_current_vectors[0][0], i_2, "-", label="i_2")
            plt.plot(time_current_vectors[0][0], i_3, "--", label="i_3")
            plt.xlabel("time / s")
            plt.ylabel("current / A")
            plt.grid()
            plt.legend()
            plt.show()

        # transfer format from fft()-output to excitation_sweep()-input
        current_list_list = []
        phi_deg_list_list = []
        for count_frequency, frequency in enumerate(frequency_list):
            currents_single_frequency = []
            phi_deg_single_frequency = []
            for count_current, _ in enumerate(time_current_vectors):
                currents_single_frequency.append(frequency_current_phase_deg_list[count_current][1][count_frequency])
                phi_deg_single_frequency.append(frequency_current_phase_deg_list[count_current][2][count_frequency])
            current_list_list.append(currents_single_frequency)
            phi_deg_list_list.append(phi_deg_single_frequency)
        center_tapped_study_excitation["linear_losses"]["frequencies"] = list(frequency_list)
        center_tapped_study_excitation["linear_losses"]["current_amplitudes"] = current_list_list
        center_tapped_study_excitation["linear_losses"]["current_phases_deg"] = phi_deg_list_list

        return center_tapped_study_excitation

    def stacked_core_center_tapped_study(self, center_tapped_study_excitation, number_primary_coil_turns: int = None,
                                         non_sine_hysteresis_correction: bool = False):
        """
        Comprehensive component analysis for center tapped transformers with dedicated choke.

        :param non_sine_hysteresis_correction: True to enable the non-sinusoidal hysteresis correction factor
        :param center_tapped_study_excitation:
        :param number_primary_coil_turns: number of primary coil turns. Needed due to a special trick to get the transformer losses without effect of the choke
        :type number_primary_coil_turns: int
        :return:
        """

        def factor_triangular_hysteresis_loss_iGSE(D, alpha):
            nominator = 2 * (D ** (1 - alpha) + (1 - D) ** (1 - alpha))
            theta = np.linspace(0, 2 * np.pi, 100)
            integrant = np.abs(np.cos(theta)) ** alpha
            denominator = np.pi ** (alpha - 1) * np.trapz(integrant, x=theta)
            return nominator / denominator

        # get the inductance
        inductance_dict = self.get_inductances(I0=1, skin_mesh_factor=1,
                                               op_frequency=center_tapped_study_excitation["hysteresis"]["frequency"],
                                               silent=self.silent)

        # Initialize the hysteresis losses with zero
        p_hyst = 0
        # print(f"{p_hyst = }")

        ps_primary_coil_turns = [150000 + i for i in range(number_primary_coil_turns)]
        self.overwrite_conductors_with_air(ps_primary_coil_turns)
        self.excitation(frequency=center_tapped_study_excitation["hysteresis"]["frequency"],
                        amplitude_list=center_tapped_study_excitation["hysteresis"]["transformer"][
                            "current_amplitudes"],
                        phase_deg_list=center_tapped_study_excitation["hysteresis"]["transformer"][
                            "current_phases_deg"],
                        plot_interpolation=False)
        self.write_simulation_parameters_to_pro_files()
        self.generate_load_litz_approximation_parameters()
        self.simulate()
        self.calculate_and_write_log()  # TODO: reuse center tapped

        log = self.read_log()
        for i in [1, 2, 3, 4]:
            res = log['single_sweeps'][0]['core_parts'][f'core_part_{i}']['hyst_losses']
            # print(f"core_part_{i} = {res}")
            p_hyst += res
        # print(f"{p_hyst = }")

        if non_sine_hysteresis_correction:
            # Correct the hysteresis loss for the triangular shaped flux density waveform
            # alpha_from_db, beta_from_db, k_from_db = mdb.MaterialDatabase(ff.silent).get_steinmetz(temperature=self.core.temperature, material_name=self.core.material, datasource="measurements",
            #                                                               datatype=mdb.MeasurementDataType.Steinmetz, measurement_setup="LEA_LK",interpolation_type="linear")
            # p_hyst = factor_triangular_hysteresis_loss_iGSE(D=0.5, alpha=alpha_from_db) * p_hyst
            # print(f"{p_hyst = }")

            ps_primary_coil_turns = [150000 + i for i in range(number_primary_coil_turns)]

        self.overwrite_air_conductors_with_conductors(list(np.array(ps_primary_coil_turns) + 1000000))
        self.excitation(frequency=center_tapped_study_excitation["hysteresis"]["frequency"],
                        amplitude_list=center_tapped_study_excitation["hysteresis"]["choke"]["current_amplitudes"],
                        phase_deg_list=center_tapped_study_excitation["hysteresis"]["choke"]["current_phases_deg"],
                        plot_interpolation=False)
        self.write_simulation_parameters_to_pro_files()
        self.generate_load_litz_approximation_parameters()
        self.simulate()
        self.calculate_and_write_log()  # TODO: reuse center tapped

        log = self.read_log()
        for i in [3, 4, 5]:
            res = log['single_sweeps'][0]['core_parts'][f'core_part_{i}']['hyst_losses']
            # print(f"core_part_{i} = {res}")
            p_hyst += res

        # print(f"{p_hyst = }")

        # calculate the winding losses # TODO: avoid meshing twice
        self.excitation_sweep(center_tapped_study_excitation["linear_losses"]["frequencies"],
                              center_tapped_study_excitation["linear_losses"]["current_amplitudes"],
                              center_tapped_study_excitation["linear_losses"]["current_phases_deg"],
                              inductance_dict=inductance_dict, core_hyst_loss=float(p_hyst))

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    # Post-Processing
    def get_inductances(self, I0: float, op_frequency: float = 0, skin_mesh_factor: float = 1,
                        visualize_last_fem_simulation: bool = False, silent: bool = False):
        """
        Performs 'open' simulations with input current for each winding and calculates the inductance matrix and
        the primary concentrated equivalent circuit.

        * 2 winding transformer:
            - inductance matrix: L_1_1, L_2_2 and M
            - primary concentrated equivalent circuits: n_conc, L_s_conc and L_h_conc

        * 3 winding transformer
            - inductance matrix: L_1_1, L_2_2, L_3_3, M_12, M_13 and M_23
            - Stray Inductance with 'Turns Ratio' n as 'Transformation Ratio' n2 and n3: L_s1, L_s2, L_s3 and L_h

        :param visualize_last_fem_simulation: Show last FEM simulation result. Defaults to False.
        :type visualize_last_fem_simulation: bool
        :param skin_mesh_factor:
        :type skin_mesh_factor: bool
        :param I0: exciting peak current in A. For all transformers, this value is used in all windings.
        :type I0: float
        :param op_frequency: operating frequency in Hz
        :type op_frequency: float
        :param silent: True for not terminal output
        :type silent: bool
        """
        if len(self.windings) == 1:
            raise NotImplementedError(
                "For inductor, this function will not be implemented. See 'flux_over_current' in 'log_electro_magnetic.json' ")

        else:
            # Data-generation with FEM simulations
            self.high_level_geo_gen(frequency=op_frequency, skin_mesh_factor=skin_mesh_factor)
            frequencies, currents, phases = ff.create_open_circuit_excitation_sweep(I0, len(self.windings),
                                                                                    op_frequency)
            self.excitation_sweep(frequency_list=frequencies, current_list_list=currents, phi_deg_list_list=phases,
                                  show_last_fem_simulation=visualize_last_fem_simulation)

            # Post-Processing
            log = self.read_log()
            self_inductances = ff.get_self_inductances_from_log(log)
            flux_linkages = ff.get_flux_linkages_from_log(log)
            coupling_matrix = ff.get_coupling_matrix(flux_linkages)
            mean_coupling_factors = ff.get_mean_coupling_factors(coupling_matrix)
            inductance_matrix = ff.get_inductance_matrix(self_inductances, mean_coupling_factors, coupling_matrix)

            if not silent:
                ff.visualize_self_inductances(self_inductances, flux_linkages, silent=self.silent)
                ff.visualize_self_resistances(self_inductances, flux_linkages, op_frequency, silent=self.silent)
                ff.visualize_flux_linkages(flux_linkages, silent=self.silent)
                # visualize_coupling_factors(coupling_matrix)
                ff.visualize_mean_coupling_factors(mean_coupling_factors, silent=self.silent)
                # visualize_mean_mutual_inductances(inductance_matrix)
                # visualize_mutual_inductances(self_inductances, coupling_matrix)
                ff.visualize_inductance_matrix(inductance_matrix, silent=self.silent)
                # print(np.array(inductance_matrix).real)
                # print("")
                # print(np.array(inductance_matrix).imag)

        if len(self.windings) == 2:
            # Self inductances
            self.L_1_1 = self_inductances[0].real
            self.L_2_2 = self_inductances[1].real

            # Main/Counter Inductance
            self.M = inductance_matrix[0][1].real

            # Mean coupling factor
            k = mean_coupling_factors[0][1]

            # # 2 winding transformer
            # # Turns Ratio n=N1/N2
            # primary_turns = 0
            # secondary_turns = 0
            # for ww in self.winding_windows:
            #     for vww in ww.virtual_winding_windows:
            #         primary_turns += vww.turns[0]
            #         secondary_turns += vww.turns[1]
            # n = primary_turns / secondary_turns
            #
            # self.femmt_print(f"\n"
            #                f"Turns Ratio:\n"
            #                f"n = {n}\n"
            #                )
            #
            #
            # l_s1 = self.L_1_1 - self.M * n
            # l_s2 = self.L_2_2 - self.M / n
            # l_h = self.M * n
            # self.femmt_print(f"\n"
            #                f"T-ECD (primary side transformed):\n"
            #                f"[Underdetermined System: 'Transformation Ratio' := 'Turns Ratio']\n"
            #                f"    - Transformation Ratio: n\n"
            #                f"    - Primary Side Stray Inductance: L_s1\n"
            #                f"    - Secondary Side Stray Inductance: L_s2\n"
            #                f"    - Primary Side Main Inductance: L_h\n"
            #                f"n := n = {n}\n"
            #                f"L_s1 = L_1_1 - M * n = {l_s1}\n"
            #                f"L_s2 = L_2_2 - M / n = {l_s2}\n"
            #                f"L_h = M * n = {l_h}\n"
            #                )

            # Stray Inductance concentrated on Primary Side
            self.n_conc = self.M / self.L_2_2
            self.L_s_conc = (1 - k ** 2) * self.L_1_1
            self.L_h_conc = self.M ** 2 / self.L_2_2

            self.femmt_print(f"\n"
                             f"T-ECD (primary side concentrated):\n"
                             f"[Underdetermined System: n := M / L_2_2  -->  L_s2 = L_2_2 - M / n = 0]\n"
                             f"    - Transformation Ratio: n\n"
                             f"    - (Primary) Stray Inductance: L_s1\n"
                             f"    - Primary Side Main Inductance: L_h\n"
                             f"n := M / L_2_2 = k * Sqrt(L_1_1 / L_2_2) = {self.n_conc}\n"
                             f"L_s1 = (1 - k^2) * L_1_1 = {self.L_s_conc}\n"
                             f"L_h = M^2 / L_2_2 = k^2 * L_1_1 = {self.L_h_conc}\n"
                             )

            inductance = TransformerInductance(
                l_h_conc=self.L_h_conc,
                l_s_conc=self.L_s_conc,
                n_conc=self.n_conc,
                M=self.M,
                L_1_1=self.L_1_1,
                L_2_2=self.L_2_2,
            )

            return dataclasses.asdict(inductance)

        if len(self.windings) == 3:
            # Self inductances
            self.L_1_1 = self_inductances[0].real
            self.L_2_2 = self_inductances[1].real
            self.L_3_3 = self_inductances[2].real

            # Main/Counter Inductance
            self.M_12 = inductance_matrix[0][1].real
            self.M_13 = inductance_matrix[0][2].real
            self.M_23 = inductance_matrix[1][2].real

            # Stray Inductance with 'Turns Ratio' n as 'Transformation Ratio' n2 and n3
            self.L_s1 = self.L_1_1 - (self.M_12 * self.M_13) / self.M_23
            self.L_s2 = self.L_2_2 - (self.M_12 * self.M_23) / self.M_13
            self.L_s3 = self.L_3_3 - (self.M_13 * self.M_23) / self.M_12
            self.L_h = (self.M_12 * self.M_13) / self.M_23
            self.n_12 = np.sqrt(self.L_1_1 / self.L_2_2)  # self.M_13 / self.M_23
            self.n_13 = np.sqrt(self.L_1_1 / self.L_3_3)  # self.M_12 / self.M_23
            self.n_23 = np.sqrt(self.L_2_2 / self.L_3_3)

            # Shortcut Inductances
            self.L_s12 = self.L_s1 + self.n_12 ** 2 * self.L_s2
            self.L_s13 = self.L_s1 + self.n_13 ** 2 * self.L_s3
            self.L_s23 = self.L_s2 + (self.n_13 / self.n_12) ** 2 * self.L_s3

            self.femmt_print(f"\n"
                             f"T-ECD (Lh on primary side):\n"
                             f"    - Primary Side Stray Inductance: L_s1\n"
                             f"    - Secondary Side Stray Inductance: L_s2\n"
                             f"    - Tertiary Side Stray Inductance: L_s3\n"
                             f"    - Transformation Ratio with respect to the primary and the Secondary: n2\n"
                             f"    - Transformation Ratio with respect to the primary and the Tertiary: n3\n"
                             f"    - Primary Side Main Inductance: L_h\n"
                             f"L_s1 = L_1_1 - M_12 * M_13 / M_23 = {self.L_s1}\n"
                             f"L_s2 = L_2_2 - M_12 * M_23 / M_13 = {self.L_s2}\n"
                             f"L_s3 = L_3_3 - M_13 * M_23 / M_12 = {self.L_s3}\n"
                             f"n_12 = np.sqrt(self.L_1_1/self.L_2_2) = {self.n_12}\n"
                             f"n_13 = np.sqrt(self.L_1_1/self.L_3_3) = {self.n_13}\n"
                             f"n_23 = np.sqrt(self.L_2_2/self.L_3_3) = {self.n_23}\n"
                             f"L_h = M_12 * M_13 / M_23 = {self.L_h}\n\n"
                             f"Shortcut Inductances L_snm measured on winding n with short applied to winding m\n"
                             f"L_s12 = L_s1 + n_12**2 * L_s2 = {self.L_s12}\n"
                             f"L_s13 = L_s1 + n_13**2 * L_s3 = {self.L_s13}\n"
                             f"L_s23 = L_s2 + (n_13/n_12)**2 * L_s3 = {self.L_s23}\n"
                             )
            """
            # Stray Inductance concentrated on Primary Side
            self.n_conc = self.M / self.L_2_2
            self.L_s_conc = (1 - k ** 2) * self.L_1_1
            self.L_h_conc = self.M ** 2 / self.L_2_2

            self.femmt_print(f"\n"
                f"T-ECD (primary side concentrated):\n"
                f"[Underdetermined System: n := M / L_2_2  -->  L_s2 = L_2_2 - M / n = 0]\n"
                f"    - Transformation Ratio: n\n"
                f"    - (Primary) Stray Inductance: L_s1\n"
                f"    - Primary Side Main Inductance: L_h\n"
                f"n := M / L_2_2 = k * Sqrt(L_1_1 / L_2_2) = {self.n_conc}\n"
                f"L_s1 = (1 - k^2) * L_1_1 = {self.L_s_conc}\n"
                f"L_h = M^2 / L_2_2 = k^2 * L_1_1 = {self.L_h_conc}\n"
                )
            """

            inductances = ThreeWindingTransformerInductance(
                M_12=self.M_12,
                M_13=self.M_13,
                M_23=self.M_23,
                L_s1=self.L_s1,
                L_s2=self.L_s2,
                L_s3=self.L_s3,
                L_h=self.L_h,
                n_12=self.n_12,
                n_13=self.n_13,
                n_23=self.n_23,
                L_s12=self.L_s12,
                L_s13=self.L_s13,
                L_s23=self.L_s23
            )
            return dataclasses.asdict(inductances)

        # self.visualize()

    def get_steinmetz_loss(self, peak_current: float = None, ki: float = 1, alpha: float = 1.2, beta: float = 2.2,
                           t_rise: float = 3e-6, t_fall: float = 3e-6,
                           f_switch: float = 100000, skin_mesh_factor: float = 0.5):

        """

        :param skin_mesh_factor:
        :param peak_current:
        :param ki:
        :param alpha:
        :param beta:
        :param t_rise:
        :param t_fall:
        :param f_switch:
        """
        peak_current = peak_current or [10, 10]

        self.core.steinmetz_loss = 1
        self.Ipeak = peak_current
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
        self.excitation(frequency=f_switch, amplitude_list=peak_current,
                        phase_deg_list=[0, 180])  # frequency and current
        self.check_model_mqs_condition()
        self.write_simulation_parameters_to_pro_files()

    def write_electro_magnetic_parameter_pro(self):
        """
        Write materials and other parameters to the "Parameter.pro" file.
        This file is generated by python and is read by gmsh to hand over some parameters.

        Source for the parameters is the MagneticComponent object.
        """
        text_file = open(os.path.join(self.file_data.electro_magnetic_folder_path, "Parameter.pro"), "w")

        # Magnetic Component Type
        if self.component_type == ComponentType.Inductor:
            text_file.write(f"Number_of_Windings = {len(self.windings)};\n")

        if self.component_type == ComponentType.Transformer:
            text_file.write(f"Number_of_Windings = {len(self.windings)};\n")

        if self.component_type == ComponentType.IntegratedTransformer:
            text_file.write(f"Number_of_Windings = {len(self.windings)};\n")

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
                # self.core.sigma = 2 * np.pi * self.frequency * epsilon_0 * f_N95_er_imag(f=self.frequency) + 1 / 6
                self.core.sigma = 1 / 6
            text_file.write(f"sigma_core = {self.core.sigma.real};\n")
            text_file.write(f"sigma_core_imag = {self.core.sigma.imag};\n")
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
        for winding_number in range(len(self.windings)):
            # -- Control Flags --
            if self.flag_excitation_type == 'current':
                text_file.write(f"Flag_ImposedVoltage = 0;\n")
            if self.flag_excitation_type == 'voltage':
                text_file.write(f"Flag_ImposedVoltage = 1;\n")
            if self.windings[winding_number].conductor_type == ConductorType.RoundLitz:
                text_file.write(f"Flag_HomogenisedModel_{winding_number + 1} = 1;\n")
            else:
                text_file.write(f"Flag_HomogenisedModel_{winding_number + 1} = 0;\n")

            # -- Geometry --
            # Core Parts
            text_file.write(f"nCoreParts = {len(self.mesh.plane_surface_core)};\n")

            turns = ff.get_number_of_turns_of_winding(winding_windows=self.winding_windows, windings=self.windings,
                                                      winding_number=winding_number)

            if self.windings[winding_number].parallel:
                text_file.write(f"NbrCond_{winding_number + 1} = 1;\n")
                text_file.write(f"AreaCell_{winding_number + 1} = {self.windings[winding_number].a_cell * turns};\n")
            else:
                text_file.write(f"NbrCond_{winding_number + 1} = {turns};\n")
                text_file.write(f"AreaCell_{winding_number + 1} = {self.windings[winding_number].a_cell};\n")

            # For stranded Conductors:
            # text_file.write(f"NbrstrandedCond = {self.turns};\n")  # redundant
            if self.windings[winding_number].conductor_type == ConductorType.RoundLitz:
                text_file.write(f"NbrStrands_{winding_number + 1} = {self.windings[winding_number].n_strands};\n")
                text_file.write(f"Fill_{winding_number + 1} = {self.windings[winding_number].ff};\n")
                # ---
                # Litz Approximation Coefficients were created with 4 layers
                # That's why here a hard-coded 4 is implemented
                # text_file.write(f"NbrLayers{num+1} = {self.n_layers[num]};\n")
                text_file.write(f"NbrLayers_{winding_number + 1} = 4;\n")

            text_file.write(f"Rc_{winding_number + 1} = {self.windings[winding_number].conductor_radius};\n")

            # -- Excitation --
            # Imposed current, current density or voltage
            if self.flag_excitation_type == 'current':
                text_file.write(f"Val_EE_{winding_number + 1} = {abs(self.current[winding_number])};\n")
                text_file.write(f"Phase_{winding_number + 1} = {np.deg2rad(self.phase_deg[winding_number])};\n")
                text_file.write(
                    f"Parallel_{winding_number + 1} = {int(self.windings[winding_number].parallel == True)};\n")

            if self.flag_excitation_type == 'current_density':
                text_file.write(f"Val_EE_{winding_number + 1} = {self.current_density[winding_number]};\n")
                raise NotImplementedError

            if self.flag_excitation_type == 'voltage':
                text_file.write(f"Val_EE_{winding_number + 1} = {self.voltage[winding_number]};\n")
                raise NotImplementedError

            self.femmt_print(f"Cell surface area: {self.windings[winding_number].a_cell} \n"
                             f"Reduced frequency: {self.red_freq[winding_number]}")

            if self.red_freq[winding_number] > 1.25 and self.windings[
                winding_number].conductor_type == ConductorType.RoundLitz:
                # TODO: Allow higher reduced frequencies
                self.femmt_print(f"Litz Coefficients only implemented for X<=1.25")
                raise Warning
            # Reduced Frequency
            text_file.write(f"Rr_{winding_number + 1} = {self.red_freq[winding_number]};\n")

            # Material Properties
            # Conductor Material
            text_file.write(f"sigma_winding_{winding_number + 1} = {self.windings[winding_number].cond_sigma};\n")

        # -- Materials --

        # Nature Constants
        text_file.write(f"mu0 = 4.e-7 * Pi;\n")
        text_file.write(f"nu0 = 1 / mu0;\n")
        text_file.write(f"e0 = {epsilon_0};\n")

        # Material Properties

        # Core Material
        # if self.frequency == 0:
        if self.core.non_linear:
            text_file.write(f"Flag_NL = 1;\n")
            text_file.write(
                f"Core_Material = {self.core.material};\n")  # relative permeability is defined at simulation runtime            text_file.write(f"Flag_NL = 0;\n")
        else:
            text_file.write(f"Flag_NL = 0;\n")
            text_file.write(f"mur = {self.core.mu_r_abs};\n")  # mur is predefined to a fixed value
            if self.core.permeability_type == PermeabilityType.FromData:
                text_file.write(f"Flag_Permeability_From_Data = 1;\n")  # mur is predefined to a fixed value
            else:
                text_file.write(f"Flag_Permeability_From_Data = 0;\n")  # mur is predefined to a fixed value
            if self.core.permeability_type == PermeabilityType.FixedLossAngle:
                text_file.write(
                    f"phi_mu_deg = {self.core.phi_mu_deg};\n")  # loss angle for complex representation of hysteresis loss
                text_file.write(
                    f"mur_real = {self.core.mu_r_abs * np.cos(np.deg2rad(self.core.phi_mu_deg))};\n")  # Real part of complex permeability
                text_file.write(
                    f"mur_imag = {self.core.mu_r_abs * np.sin(np.deg2rad(self.core.phi_mu_deg))};\n")  # Imaginary part of complex permeability
                text_file.write(
                    f"Flag_Fixed_Loss_Angle = 1;\n")  # loss angle for complex representation of hysteresis loss
            else:
                text_file.write(
                    f"Flag_Fixed_Loss_Angle = 0;\n")  # loss angle for complex representation of hysteresis loss

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
        text_file.write(
            f"DirResValsCore = \"{self.file_data.e_m_values_folder_path.replace(backslash, '/')}/core_parts/\";\n")
        for i in range(1, len(self.windings) + 1):
            text_file.write(
                f"DirResValsWinding_{i} = \"{self.file_data.e_m_values_folder_path.replace(backslash, '/')}/Winding_{i}/\";\n")
            # text_file.write(f"DirResValsSecondary = \"{self.file_data.e_m_values_folder_path.replace(backslash, '/')}/Secondary/\";\n")
            # text_file.write(f"DirResValsTertiary = \"{self.file_data.e_m_values_folder_path.replace(backslash, '/')}/Tertiary/\";\n")
        text_file.write(f"DirResCirc = \"{self.file_data.e_m_circuit_folder_path.replace(backslash, '/')}/\";\n")
        text_file.write(f"OptionPos = \"{self.file_data.results_folder_path.replace(backslash, '/')}/option.pos\";\n")
        text_file.write(
            f"DirStrandCoeff = \"{self.file_data.e_m_strands_coefficients_folder_path.replace(backslash, '/')}/\";\n")

        # Visualisation
        if self.plot_fields == "standard":
            text_file.write(f"Flag_show_standard_fields = 1;\n")
        else:
            text_file.write(f"Flag_show_standard_fields = 0;\n")

        text_file.close()

    def calculate_and_write_log(self, sweep_number: int = 1, currents: List = None, frequencies: List = None,
                                inductance_dict: dict = None, core_hyst_losses: float = None):
        """
        Method reads back the results from the .dat result files created by the ONELAB simulation client and stores
        them in a result dictionary (JSON log file).

        This file includes:
         * results (losses, ...) of the simulation
         * geometry parameters of the given simulation
         * optional parameters can be added to the log, like the inductance values or the core hysteresis losses from
                    external simulations

        :param sweep_number: Number of sweep iterations that were done before. For a single simulation sweep_number = 1
        :type sweep_number: int
        :param currents: Current values of the sweep iterations. Not needed for single simulation
        :type currents: list
        :param frequencies: frequencies values of the sweep iterations. Not needed for single simulation
        :type frequencies: list
        :param inductance_dict: Optional inductance dict to include in the logfile.
        :type inductance_dict: dict
        :param core_hyst_losses: Optional core hysteresis losses value from another simulation. If a value is given, the external value is used in the result-log. Otherwise, the hysteresis losses of the fundamental frequency is used
        :type core_hyst_losses: float
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
            winding_name = ["Winding_" + str(i) for i in range(1, len(self.windings) + 1)]

            for winding_number in range(len(self.windings)):

                # Create empty winding dictionary
                # create dictionary winding_dict with 'turn_losses' as list of the j=number_turns turn losses.
                # Single values related to one winding are added as 'winding_losses' etc.
                winding_dict = {"turn_losses": [],
                                "flux": [],
                                "flux_over_current": [],
                                # "mag_field_energy": [],
                                "V": []}

                # Number turns
                turns = ff.get_number_of_turns_of_winding(winding_windows=self.winding_windows, windings=self.windings,
                                                          winding_number=winding_number)

                winding_dict["number_turns"] = turns

                # Currents
                if sweep_number > 1:
                    # sweep_simulation -> get currents from passed currents
                    complex_current_phasor = currents[sweep_run][winding_number]
                else:
                    # single_simulation -> get current from instance variable
                    complex_current_phasor = self.current[winding_number]

                # Store complex value as list in json (because json is not natively capable of complex values)
                winding_dict["I"] = [complex_current_phasor.real, complex_current_phasor.imag]

                # Case litz: Load homogenized results
                if self.windings[winding_number].conductor_type == ConductorType.RoundLitz:
                    winding_dict["winding_losses"] = \
                    self.load_result(res_name=f"j2H_{winding_number + 1}", last_n=sweep_number)[sweep_run]
                    for turn in range(0, winding_dict["number_turns"]):
                        winding_dict["turn_losses"].append(
                            self.load_result(res_name=winding_name[winding_number] + f"/Losses_turn_{turn + 1}",
                                             last_n=sweep_number)[sweep_run])

                # Case Solid: Load results, (pitfall for parallel windings results are only stored in one turn!)
                else:
                    winding_dict["winding_losses"] = \
                    self.load_result(res_name=f"j2F_{winding_number + 1}", last_n=sweep_number)[sweep_run]
                    if self.windings[winding_number].parallel:
                        winding_dict["turn_losses"].append(
                            self.load_result(res_name=winding_name[winding_number] + f"/Losses_turn_{1}",
                                             last_n=sweep_number)[sweep_run])
                    else:
                        for turn in range(0, winding_dict["number_turns"]):
                            winding_dict["turn_losses"].append(
                                self.load_result(res_name=winding_name[winding_number] + f"/Losses_turn_{turn + 1}",
                                                 last_n=sweep_number)[sweep_run])

                # Magnetic Field Energy
                # winding_dict["mag_field_energy"].append(self.load_result(res_name=f"ME", last_n=sweep_number)[sweep_run])
                # winding_dict["mag_field_energy"].append(self.load_result(res_name=f"ME", part="imaginary", last_n=sweep_number)[sweep_run])

                # Voltage
                winding_dict["V"].append(
                    self.load_result(res_name=f"Voltage_{winding_number + 1}", part="real", last_n=sweep_number)[
                        sweep_run])
                winding_dict["V"].append(
                    self.load_result(res_name=f"Voltage_{winding_number + 1}", part="imaginary", last_n=sweep_number)[
                        sweep_run])
                complex_voltage_phasor = complex(winding_dict["V"][0], winding_dict["V"][1])

                # Inductance
                # winding_dict["self_inductance"].append(self.load_result(res_name=f"L_{winding + 1}{winding + 1}", part="real", last_n=sweep_number)[sweep_run])
                # winding_dict["self_inductance"].append(self.load_result(res_name=f"L_{winding + 1}{winding + 1}", part="imaginary", last_n=sweep_number)[sweep_run])
                # Inductance from voltage

                if complex_current_phasor == 0 or sweep_dict["f"] == 0:  # if-statement to avoid div by zero error
                    winding_dict["flux_over_current"] = [0, 0]
                else:
                    winding_dict["flux_over_current"].append((complex_voltage_phasor / (
                                complex(0, 1) * 2 * np.pi * complex_current_phasor * sweep_dict["f"])).real)
                    winding_dict["flux_over_current"].append((complex_voltage_phasor / (
                                complex(0, 1) * 2 * np.pi * complex_current_phasor * sweep_dict["f"])).imag)

                # Flux
                winding_dict["flux"].append(
                    self.load_result(res_name=f"Flux_Linkage_{winding_number + 1}", last_n=sweep_number)[sweep_run])
                winding_dict["flux"].append(
                    self.load_result(res_name=f"Flux_Linkage_{winding_number + 1}", part="imaginary",
                                     last_n=sweep_number)[sweep_run])
                # Flux from voltage
                # winding_dict["flux"].append((complex(winding_dict["self_inductance"][-2], winding_dict["self_inductance"][-1])*self.current[winding]).real)  # (L*I).real
                # winding_dict["flux"].append((complex(winding_dict["self_inductance"][-2], winding_dict["self_inductance"][-1])*self.current[winding]).imag)  # (L*I).imag

                # Power
                # using 'winding_dict["V"][0]' to get first element (real part) of V. Use winding_dict["I"][0] to avoid typeerror
                winding_dict["P"] = (complex_voltage_phasor * complex_current_phasor.conjugate() / 2).real
                winding_dict["Q"] = (complex_voltage_phasor * complex_current_phasor.conjugate() / 2).imag
                winding_dict["S"] = np.sqrt(winding_dict["P"] ** 2 + winding_dict["Q"] ** 2)

                sweep_dict[f"winding{winding_number + 1}"] = winding_dict

            # Core losses TODO: Choose between Steinmetz or complex core losses
            sweep_dict["core_eddy_losses"] = self.load_result(res_name="CoreEddyCurrentLosses", last_n=sweep_number)[
                sweep_run]
            sweep_dict["core_hyst_losses"] = self.load_result(res_name="p_hyst", last_n=sweep_number)[sweep_run]

            # Core Part losses
            # If there are multiple core parts, calculate losses for each part individually
            # the core losses are caluclated by summing the eddy and hyst losses
            if len(self.mesh.plane_surface_core) > 1:
                sweep_dict["core_parts"] = {}
                for i in range(0, len(self.mesh.plane_surface_core)):
                    sweep_dict["core_parts"][f"core_part_{i + 1}"] = {}
                    sweep_dict["core_parts"][f"core_part_{i + 1}"]["eddy_losses"] = \
                    self.load_result(res_name=f"core_parts/CoreEddyCurrentLosses_{i + 1}", last_n=sweep_number)[sweep_run]
                    sweep_dict["core_parts"][f"core_part_{i + 1}"]["hyst_losses"] = \
                    self.load_result(res_name=f"core_parts/p_hyst_{i + 1}", last_n=sweep_number)[sweep_run]
                    # finding the total losses for every core_part
                    eddy = sweep_dict["core_parts"][f"core_part_{i + 1}"]["eddy_losses"]
                    hyst = sweep_dict["core_parts"][f"core_part_{i + 1}"]["hyst_losses"]
                    sweep_dict["core_parts"][f"core_part_{i + 1}"][f"total_core_part_{i + 1}"] = eddy + hyst

            # For a single core part, total losses are simply the sum of eddy and hysteresis losses
            else:
                # if I have only one core_part
                sweep_dict["core_parts"] = {}  # parent dictionary is initialized first
                sweep_dict["core_parts"]["core_part_1"] = {}
                sweep_dict["core_parts"]["core_part_1"]["total_core_part_1"] = sweep_dict["core_eddy_losses"] + sweep_dict["core_hyst_losses"]

            # Sum losses of all windings of one single run
            sweep_dict["all_winding_losses"] = sum(
                sweep_dict[f"winding{d + 1}"]["winding_losses"] for d in range(len(self.windings)))

            log_dict["single_sweeps"].append(sweep_dict)

        # Total losses of excitation sweep
        # Sum losses of all sweep runs. For core losses just use hyst_losses of the fundamental frequency.
        # Also needed as excitation for steady state thermal simulations

        # Single Windings
        for winding_number in range(len(self.windings)):
            # Number of turns per conductor
            turns = 0
            for ww in self.winding_windows:
                for vww in ww.virtual_winding_windows:
                    for conductor in vww.windings:
                        if conductor.winding_number == winding_number:
                            turns += vww.turns[conductor.winding_number]

            log_dict["total_losses"][f"winding{winding_number + 1}"] = {
                "total": sum(sum(log_dict["single_sweeps"][d][f"winding{winding_number + 1}"]["turn_losses"]) for d in
                             range(len(log_dict["single_sweeps"]))),
                "turns": []
            }
            if self.windings[winding_number].parallel:
                log_dict["total_losses"][f"winding{winding_number + 1}"]["turns"].append(
                    sum(log_dict["single_sweeps"][d][f"winding{winding_number + 1}"]["turn_losses"][0] for d in
                        range(len(log_dict["single_sweeps"]))))
            else:
                for turn in range(0, turns):
                    log_dict["total_losses"][f"winding{winding_number + 1}"]["turns"].append(
                        sum(log_dict["single_sweeps"][d][f"winding{winding_number + 1}"]["turn_losses"][turn] for d in
                            range(len(log_dict["single_sweeps"]))))

        # Winding (all windings)
        log_dict["total_losses"]["all_windings"] = sum(
            log_dict["single_sweeps"][d]["all_winding_losses"] for d in range(len(log_dict["single_sweeps"])))

        # Core
        log_dict["total_losses"]["eddy_core"] = sum(
            log_dict["single_sweeps"][d]["core_eddy_losses"] for d in range(len(log_dict["single_sweeps"])))
        # In the total losses calculation, individual core part losses are also accounted for
        if len(self.mesh.plane_surface_core) > 1:
            for i in range(len(self.mesh.plane_surface_core)):
                log_dict["total_losses"][f"total_core_part_{i + 1}"] = sweep_dict["core_parts"][f"core_part_{i + 1}"][
                    f"total_core_part_{i + 1}"]
        # If there is only one core part, set its total losses directly
        if len(self.mesh.plane_surface_core) == 1:

            log_dict["total_losses"]["total_core_part_1"] = sweep_dict["core_parts"]["core_part_1"]["total_core_part_1"]


        if isinstance(core_hyst_losses, float) or isinstance(core_hyst_losses, int):
            log_dict["total_losses"]["hyst_core_fundamental_freq"] = core_hyst_losses
        elif core_hyst_losses is not None:
            raise ValueError("External core hysteresis losses are given in non-float format")
        else:
            log_dict["total_losses"]["hyst_core_fundamental_freq"] = log_dict["single_sweeps"][fundamental_index][
                "core_hyst_losses"]

        # Total losses of inductive component according to single or sweep simulation
        log_dict["total_losses"]["core"] = log_dict["total_losses"]["hyst_core_fundamental_freq"] + \
                                           log_dict["total_losses"]["eddy_core"]

        log_dict["total_losses"]["total_losses"] = log_dict["total_losses"]["hyst_core_fundamental_freq"] + \
                                                   log_dict["total_losses"]["eddy_core"] + log_dict["total_losses"][
                                                       "all_windings"]

        # ---- Introduce calculations for writing the misc-dict into the result-log ----
        wire_type_list = []
        for winding_number in self.windings:
            wire_type_list.append(winding_number.conductor_type.name)

        single_strand_cross_section_list = []
        for winding_number in self.windings:
            if winding_number.strand_radius:
                single_strand_cross_section = winding_number.strand_radius ** 2 * np.pi
                single_strand_cross_section_list.append(single_strand_cross_section)
            else:
                single_strand_cross_section_list.append(None)

        wire_weight_list = self.calculate_wire_weight()
        core_weight = self.calculate_core_weight()

        # ---- Miscellaneous ----
        log_dict["misc"] = {
            "core_2daxi_volume": self.calculate_core_volume(), #self.calculate_core_volume(),
            "core_2daxi_total_volume": self.calculate_core_volume_with_air(),
            "core_2daxi_weight": core_weight,
            "wire_lengths": self.calculate_wire_lengths(),
            "wire_volumes": self.calculate_wire_volumes(),
            "wire_weight": wire_weight_list,
            "core_cost": ff.cost_function_core(core_weight, core_type="ferrite"),
            "winding_cost": ff.cost_function_winding(wire_weight_list=wire_weight_list, wire_type_list=wire_type_list,
                                                     single_strand_cross_section_list=single_strand_cross_section_list),
            "total_cost_incl_margin": ff.cost_function_total(core_weight, core_type="ferrite",
                                                             wire_weight_list=wire_weight_list,
                                                             wire_type_list=wire_type_list,
                                                             single_strand_cross_section_list=single_strand_cross_section_list)
        }

        # ---- Print current configuration ----
        log_dict["simulation_settings"] = MagneticComponent.encode_settings(self)

        if isinstance(inductance_dict, Dict):
            log_dict["inductances"] = inductance_dict

        # ====== save data as JSON ======
        with open(self.file_data.e_m_results_log_path, "w+", encoding='utf-8') as outfile:
            json.dump(log_dict, outfile, indent=2, ensure_ascii=False)

    def read_log(self) -> Dict:
        """
        Read results from electromagnetic simulation
        :return: Logfile as a dictionary
        :rtype: Dict
        """
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
        self.femmt_print(f"\n---\n"
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
            gmsh.open(os.path.join(self.file_data.e_m_fields_folder_path, "j2F_density.pos"))
            gmsh.option.setNumber(f"View[{view}].ScaleType", 2)
            gmsh.option.setNumber(f"View[{view}].RangeType", 2)
            gmsh.option.setNumber(f"View[{view}].SaturateValues", 1)
            gmsh.option.setNumber(f"View[{view}].CustomMin", gmsh.option.getNumber(f"View[{view}].Min") + epsilon)
            gmsh.option.setNumber(f"View[{view}].CustomMax", gmsh.option.getNumber(f"View[{view}].Max"))
            gmsh.option.setNumber(f"View[{view}].ColormapNumber", 1)
            gmsh.option.setNumber(f"View[{view}].IntervalsType", 2)
            gmsh.option.setNumber(f"View[{view}].NbIso", 40)
            gmsh.option.setNumber(f"View[{view}].ShowTime", 0)
            self.femmt_print(gmsh.option.getNumber(f"View[{view}].Max"))
            view += 1

        if any(self.windings[i].conductor_type == ConductorType.RoundLitz for i in range(len(self.windings))):
            # Ohmic losses (weighted effective value of current density)
            gmsh.open(os.path.join(self.file_data.e_m_fields_folder_path, "jH_density.pos"))
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

    def load_result(self, res_name: str, res_type: str = "value", last_n: int = 1, part: str = "real",
                    position: int = 0):
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
        if res_type == "value":
            res_path = self.file_data.e_m_values_folder_path
        if res_type == "circuit":
            res_path = self.file_data.e_m_circuit_folder_path


        with open(os.path.join(res_path, f"{res_name}.dat")) as fd:
            lines = fd.readlines()[-last_n:]

            if part == "real":
                result = [float(line.split(sep=' ')[1 + 2 * position + 1]) for n, line in enumerate(lines)]
            if part == "imaginary":
                result = [float(line.split(sep=' ')[2 + 2 * position + 1]) for n, line in enumerate(lines)]

            return result

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    # Litz Approximation [internal methods]
    def generate_load_litz_approximation_parameters(self):
        """
        Used to determine the litz-approximation coefficients.

        Checks if litz-approximation parameters exists. In case of non-existing litz-parameters for the certain
        litz, the litz parameters are generated directly
        """
        for num in range(len(self.windings)):
            if self.windings[num].conductor_type == ConductorType.RoundLitz:
                # ---
                # Litz Approximation Coefficients were created with 4 layers
                # That's why here a hard-coded 4 is implemented
                # if os.path.isfile(self.path +
                # f"/Strands_Coefficients/coeff/pB_RS_la{self.windings[num].ff}_{self.n_layers[num]}layer.dat"):
                if os.path.exists(os.path.join(self.file_data.e_m_strands_coefficients_folder_path, "coeff",
                                               f"pB_RS_la{self.windings[num].ff}_4layer.dat")):
                    self.femmt_print("Coefficients for stands approximation are found.")

                else:
                    # Rounding X to fit it with corresponding parameters from the database
                    X = self.red_freq[num]
                    X = np.around(X, decimals=3)
                    self.femmt_print(f"Rounded Reduced frequency X = {X}")
                    self.create_strand_coeff(num)

    def create_strand_coeff(self, winding_number: int) -> None:
        """
        Creates the initial strand coefficients for a certain litz wire.
        This function comes into account in case of the litz-coefficients are not known so far.
        After generating the litz-coefficients, the results are stored to avoid a second time-consuming
        strand coefficients generation.

        This function sends commands via text-file to gmsh to generate the strand coefficients. The onelab (gmsh) client
        is generated in a new instance.

        :param winding_number: Winding number
        :type winding_number: int
        """
        self.femmt_print(f"\n"
                         f"Pre-Simulation\n"
                         f"-----------------------------------------\n"
                         f"Create coefficients for strands approximation\n")

        text_file = open(os.path.join(self.file_data.e_m_strands_coefficients_folder_path, "PreParameter.pro"), "w")

        # Litz Approximation Coefficients are created with 4 layers
        # That's why here a hard-coded 4 is implemented
        # text_file.write(f"NbrLayers = {self.n_layers[num]};\n")
        text_file.write(f"NbrLayers = 4;\n")
        text_file.write(f"Fill = {self.windings[winding_number].ff};\n")
        text_file.write(f"Rc = {self.windings[winding_number].strand_radius};\n")  # double named!!! must be changed
        text_file.close()
        cell_geo = os.path.join(self.file_data.e_m_strands_coefficients_folder_path, "cell.geo")

        verbose = "-verbose 1" if self.silent else "-verbose 5"

        # Run gmsh as a sub client
        gmsh_client = os.path.join(self.file_data.onelab_folder_path, "gmsh")
        self.onelab_client.runSubClient("myGmsh", gmsh_client + " " + cell_geo + " -2 " + verbose)

        modes = [1, 2]  # 1 = "skin", 2 = "proximity"
        reduced_frequencies = np.linspace(0, 1.25, 6)  # must be even
        for mode in modes:
            for rf in reduced_frequencies:
                # -- Pre-Simulation Settings --
                text_file = open(os.path.join(self.file_data.e_m_strands_coefficients_folder_path, "PreParameter.pro"),
                                 "w")
                text_file.write(f"Rr_cell = {rf};\n")
                text_file.write(f"Mode = {mode};\n")
                # Litz Approximation Coefficients are created with 4 layers
                # That's why here a hard-coded 4 is implemented
                # text_file.write(f"NbrLayers = {self.n_layers[num]};\n")
                text_file.write(f"NbrLayers = 4;\n")
                text_file.write(f"Fill = {self.windings[winding_number].ff};\n")
                text_file.write(
                    f"Rc = {self.windings[winding_number].strand_radius};\n")  # double named!!! must be changed
                text_file.close()

                # get model file names with correct path
                input_file = os.path.join(self.file_data.e_m_strands_coefficients_folder_path, "cell_dat.pro")
                cell = os.path.join(self.file_data.e_m_strands_coefficients_folder_path, "cell.pro")

                # Run simulations as sub clients
                mygetdp = os.path.join(self.file_data.onelab_folder_path, "getdp")
                self.onelab_client.runSubClient("myGetDP",
                                                mygetdp + " " + cell + " -input " + input_file + " -solve MagDyn_a " + verbose)

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

        files = [os.path.join(coeff_folder, f"pB_RS_la{self.windings[winding_number].ff}_4layer.dat"),
                 os.path.join(coeff_folder, f"pI_RS_la{self.windings[winding_number].ff}_4layer.dat"),
                 os.path.join(coeff_folder, f"qB_RS_la{self.windings[winding_number].ff}_4layer.dat"),
                 os.path.join(coeff_folder, f"qI_RS_la{self.windings[winding_number].ff}_4layer.dat")]
        for i in range(0, 4):
            with fileinput.FileInput(files[i], inplace=True) as file:
                for line in file:
                    print(line.replace(' 0\n', '\n'), end='')

        # Corrects pB coefficient error at 0Hz
        # Must be changed in future in cell.pro
        for i in range(0, 4):
            with fileinput.FileInput(files[i], inplace=True) as file:
                for line in file:
                    print(line.replace(' 0\n', ' 1\n'), end='')

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

        self.file_data.create_folders(self.file_data.femm_folder_path)

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
                # self.core.sigma = 2 * np.pi * self.frequency * epsilon_0 * f_N95_er_imag(f=self.frequency) + 1 / 6
                self.core.sigma = 1 / 6

        self.femmt_print(f"{self.core.permeability_type=}, {self.core.sigma=}")
        if self.core.permeability_type == PermeabilityType.FixedLossAngle:
            femm.mi_addmaterial('Ferrite', self.core.mu_r_abs, self.core.mu_r_abs, 0, 0, self.core.sigma / 1e6, 0, 0, 1,
                                0, self.core.phi_mu_deg, self.core.phi_mu_deg)
        elif self.core.permeability_type == PermeabilityType.RealValue:
            femm.mi_addmaterial('Ferrite', self.core.mu_r_abs, self.core.mu_r_abs, 0, 0, self.core.sigma / 1e6, 0, 0, 1,
                                0, self.core.phi_mu_deg, self.core.phi_mu_deg)
        else:
            femm.mi_addmaterial('Ferrite', self.core.mu_r_abs, self.core.mu_r_abs, 0, 0, self.core.sigma / 1e6, 0, 0, 1,
                                0, 0, 0)
        femm.mi_addmaterial('Air', 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0)

        for i in range(0, len(self.windings)):
            if self.windings[i].conductor_type == ConductorType.RoundLitz:
                femm.mi_addmaterial('Litz', 1, 1, 0, 0, self.windings[i].cond_sigma / 1e6, 0, 0, 1, 5, 0, 0,
                                    self.windings[i].n_strands,
                                    2 * 1000 * self.windings[i].strand_radius)  # type := 5. last argument
                self.femmt_print(f"Number of strands: {self.windings[i].n_strands}")
                self.femmt_print(f"Diameter of strands in mm: {2 * 1000 * self.windings[i].strand_radius}")
            if self.windings[i].conductor_type == ConductorType.RoundSolid:
                femm.mi_addmaterial('Copper', 1, 1, 0, 0, self.windings[i].cond_sigma / 1e6, 0, 0, 1, 0, 0, 0, 0, 0)

        # == Circuit ==
        # coil as seen from the terminals.
        for i in range(len(self.windings)):
            femm.mi_addcircprop('Winding' + str(i + 1), current[i] * sign[i], int(not self.windings[i].parallel))

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
                    femm.mi_selectlabel(self.two_d_axi.p_conductor[num][5 * i][0],
                                        self.two_d_axi.p_conductor[num][5 * i][1])

                    winding_name = 'Winding' + str(num + 1)
                    if self.windings[num].conductor_type == ConductorType.RoundLitz:
                        femm.mi_setblockprop('Litz', 1, 0, winding_name, 0, num + 2, 1)
                    else:
                        femm.mi_setblockprop('Copper', 1, 0, winding_name, 0, num + 2, 1)

                    femm.mi_clearselected()

        # Define an "open" boundary condition using the built-in function:
        femm.mi_makeABC()
        """
        # Alternative BC
        region_add = 1.1

        femm.mi_drawrectangle(0, region_add*self.two_d_axi.p_outer[0][1], region_add*self.two_d_axi.p_outer[3][0], 
        region_add*self.two_d_axi.p_outer[3][1])
        # mi_addboundprop('Asymptotic', 0, 0, 0, 0, 0, 0, 1 / (mu_0 * bound.width), 0, 2); % Mixed
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
        self.femmt_print('Flux density at the center of the bar is %g T' % np.abs(b0[1]))
        b1 = femm.mo_getb(0.01, 0.05)
        self.femmt_print(f"Flux density at r=1cm, z=5cm is {np.abs(b1[1])} T")

        # The program will report the terminal properties of the circuit:
        # current, voltage, and flux linkage
        vals = femm.mo_getcircuitproperties('icoil')


        # [i, v, \[Phi]] = MOGetCircuitProperties["icoil"]

        # If we were interested in inductance, it could be obtained by
        # dividing flux linkage by current
        L = 1000 * np.abs(vals[2]) / np.abs(vals[0])
        self.femmt_print('The self-inductance of the coil is %g mH' % L)

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

        # self.femmt_print(hyst_loss)
        # tmp = femm.mo_getcircuitproperties('Primary')
        # self.femmt_print(tmp)
        # self.tot_loss_femm = 0.5 * tmp[0] * tmp[1]
        # self.femmt_print(self.tot_loss_femm)

        # Write Circuit Properties
        # log["Circuit Properties"] = femm.mo_getcircuitproperties('Primary')

        # Write Hysteresis Losses
        femm.mo_groupselectblock(0)
        log["Hysteresis Losses"] = femm.mo_blockintegral(3)
        femm.mo_clearblock()

        # Primary Winding circuit Properties
        circuit_properties_primary = femm.mo_getcircuitproperties('Primary')
        log["Primary Current"] = circuit_properties_primary[0]
        log["Primary Voltage"] = [circuit_properties_primary[1].real, circuit_properties_primary[1].imag]
        log["Primary Flux"] = [circuit_properties_primary[2].real, circuit_properties_primary[2].imag]
        log["Primary Self Inductance"] = [circuit_properties_primary[2].real / circuit_properties_primary[0],
                                          circuit_properties_primary[2].imag / circuit_properties_primary[0]]
        log["Primary Mean Power"] = [0.5 * circuit_properties_primary[1].real * circuit_properties_primary[0],
                                     0.5 * circuit_properties_primary[1].imag * circuit_properties_primary[0]]
        for i in range(len(self.windings)):
            circuit_properties = femm.mo_getcircuitproperties('Winding' + str(i + 1))
            log["Winding" + str(i + 1) + " Current"] = circuit_properties[0]
            log["Winding" + str(i + 1) + " Voltage"] = [circuit_properties[1].real, circuit_properties[1].imag]
            log["Winding" + str(i + 1) + " Flux"] = [circuit_properties[2].real, circuit_properties[2].imag]
            if circuit_properties[0] != 0:
                log["Winding" + str(i + 1) + " Self Inductance"] = [circuit_properties[2].real / circuit_properties[0],
                                                                    circuit_properties[2].imag / circuit_properties[0]]
            else:
                log["Winding" + str(i + 1) + " Self Inductance"] = [0, 0]
            log["Winding" + str(i + 1) + " Mean Power"] = [0.5 * circuit_properties[1].real * circuit_properties[0],
                                                           0.5 * circuit_properties[1].imag * circuit_properties[0]]

            # Obtain winding losses for each winding
            femm.mo_groupselectblock(i + 2)
            log["Winding" + str(i + 1) + " Losses"] = femm.mo_blockintegral(6).real
            femm.mo_clearblock()

        json.dump(log, file, indent=2, ensure_ascii=False)
        file.close()

    @staticmethod
    def calculate_point_average(x1: float, y1: float, x2: float, y2: float) -> Tuple[float, float]:
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

    def femm_thermal_validation(self, thermal_conductivity_dict: Dict, boundary_temperature: Dict, case_gap_top: float,
                                case_gap_right: float, case_gap_bot: float):
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
        # 2D Mesh and FEM interfaces (only for Windows machines)
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
        losses = read_results_log(self.file_data.e_m_results_log_path)

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
        k_case = thermal_conductivity_dict["case"][
            "top"]  # Does not matter when the regions all have the same thermal conductivity.
        q_vol_case = 0
        # c_case = 0.01
        c_case = 0

        # Air gap
        k_air_gap = thermal_conductivity_dict["air_gaps"]
        q_vol_air_gap = 0
        c_air_gap = 0

        # Setup winding list
        winding_losses_list = []
        for i in range(1, 10):
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

        # Check the type of conductor for this winding
        for winding_index, winding in enumerate(winding_losses_list):
            for i in range(len(winding)):
                # Add a new material to FEMM with Litz conductor properties
                femm.hi_addmaterial(f'Wire_{winding_index}_{i}', k_wire, k_wire,
                                    calculate_heat_flux_round_wire(winding[i], wire_radii[winding_index],
                                                                   wire_distances[winding_index][i]), c_wire)
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

            # Close air gap to separate from air
            femm.hi_drawline(self.two_d_axi.p_air_gaps[1, 0], self.two_d_axi.p_air_gaps[1, 1],
                             self.two_d_axi.p_air_gaps[3, 0], self.two_d_axi.p_air_gaps[3, 1])
        else:
            raise Exception("Negative air gap number is not allowed")

        # Add case
        femm.hi_drawline(0, self.two_d_axi.p_outer[2, 1], 0,
                         self.two_d_axi.p_outer[2, 1] + case_gap_top)  # Top left line
        femm.hi_drawline(0, self.two_d_axi.p_outer[2, 1] + case_gap_top, self.two_d_axi.p_outer[3, 0] + case_gap_right,
                         self.two_d_axi.p_outer[3, 1] + case_gap_top)  # Top line
        femm.hi_drawline(self.two_d_axi.p_outer[3, 0] + case_gap_right, self.two_d_axi.p_outer[3, 1] + case_gap_top,
                         self.two_d_axi.p_outer[1, 0] + case_gap_right,
                         self.two_d_axi.p_outer[1, 1] - case_gap_bot)  # Right line
        femm.hi_drawline(self.two_d_axi.p_outer[1, 0] + case_gap_right, self.two_d_axi.p_outer[1, 1] - case_gap_bot, 0,
                         self.two_d_axi.p_outer[0, 1] - case_gap_bot)  # Bottom line
        femm.hi_drawline(0, self.two_d_axi.p_outer[0, 1] - case_gap_bot, 0,
                         self.two_d_axi.p_outer[0, 1])  # Bottom right line

        # Create boundary
        # femm.hi_selectsegment(*self.calculatePointAverage(0, self.two_d_axi.p_outer[2, 1], 0, self.two_d_axi.p_outer[2, 1] + caseGapTop))
        femm.hi_selectsegment(*MagneticComponent.calculate_point_average(0, self.two_d_axi.p_outer[2, 1] + case_gap_top,
                                                                         self.two_d_axi.p_outer[3, 0] + case_gap_right,
                                                                         self.two_d_axi.p_outer[3, 1] + case_gap_top))
        femm.hi_setsegmentprop("NeumannBoundary", 0, 1, 0, 2, "<None>")
        femm.hi_clearselected()

        femm.hi_selectsegment(*MagneticComponent.calculate_point_average(self.two_d_axi.p_outer[3, 0] + case_gap_right,
                                                                         self.two_d_axi.p_outer[3, 1] + case_gap_top,
                                                                         self.two_d_axi.p_outer[1, 0] + case_gap_right,
                                                                         self.two_d_axi.p_outer[1, 1] - case_gap_bot))
        femm.hi_selectsegment(*MagneticComponent.calculate_point_average(self.two_d_axi.p_outer[1, 0] + case_gap_right,
                                                                         self.two_d_axi.p_outer[1, 1] - case_gap_bot, 0,
                                                                         self.two_d_axi.p_outer[0, 1] - case_gap_bot))
        # femm.hi_selectsegment(*self.calculatePointAverage(0, self.two_d_axi.p_outer[0, 1] - caseGapBot, 0, self.two_d_axi.p_outer[0, 1]))
        femm.hi_setsegmentprop("Boundary", 0, 1, 0, 2, "<None>")
        femm.hi_clearselected()

        # Add case material
        material_x, material_y = self.calculate_point_average(0, self.two_d_axi.p_outer[2, 1], 0,
                                                              self.two_d_axi.p_outer[2, 1] + case_gap_top)
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
                femm.hi_selectlabel(self.two_d_axi.p_conductor[num][5 * i][0],
                                    self.two_d_axi.p_conductor[num][5 * i][1])
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
        """Encodes the magnetic component in a dictionary.

        :param o: Magnetic component containing the model.
        :type o: MagneticComponent
        :return: Model encodes as dictionary
        :rtype: Dict
        """
        content = {
            "simulation_name": o.simulation_name,
            "date": datetime.today().strftime('%Y-%m-%d %H:%M:%S'),
            "component_type": o.component_type.name,
            "working_directory": o.file_data.working_directory,
            "core": o.core.to_dict(),
            "winding_windows": [ww.to_dict() for ww in o.winding_windows],
        }

        if o.air_gaps is not None:
            content["air_gaps"] = o.air_gaps.to_dict()

        if o.insulation is not None:
            content["insulation"] = o.insulation.to_dict()

        if o.stray_path is not None:
            content["stray_path"] = o.stray_path.__dict__

        return content

    @staticmethod
    def decode_settings_from_log(log_file_path: str, working_directory: str = None, verbosity: bool = False):
        """
        Reads the given log and returns the magnetic component from th elog.

        :param log_file_path: Path to the log file
        :type log_file_path: str
        :param working_directory: If the working directory shall be a different from the log file enter a new one here, defaults to None
        :type working_directory: str, optional
        :param silent: True to avoid terminal outputs. Helps to speed up simulations.
        :type silent: bool
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
            geo = MagneticComponent(component_type=ComponentType[settings["component_type"]], working_directory=cwd,
                                    verbosity = 2)

            settings["core"]["loss_approach"] = LossApproach[settings["core"]["loss_approach"]]
            core_type = settings["core"]["core_type"]
            #print(core_type)
            if core_type == CoreType.Single:
                core_dimensions = SingleCoreDimensions(core_inner_diameter=settings["core"]["core_inner_diameter"],
                                                                window_w=settings["core"]["window_w"],
                                                                window_h=settings["core"]["window_h"],
                                                                core_h=settings["core"]["core_h"])

            elif core_type == CoreType.Stacked:
                core_dimensions = StackedCoreDimensions(core_inner_diameter=settings["core"]["core_inner_diameter"],
                                                                window_w=settings["core"]["window_w"],
                                                                window_h_bot=settings["core"]["window_h_bot"],
                                                                window_h_top=settings["core"]["window_h_top"])
                                                                #ToDo: core_h not implemented yet.
                                                                #core_h=settings["core"]["core_h"])
            else:
                raise ValueError("unknown core_type for decoding from result_log.")


            if isinstance(settings["core"]["sigma"], List):
                # in case of sigma is a complex number, it is given as a list and needs to translated to complex.
                settings["core"]["sigma"] = complex(settings["core"]["sigma"][0], settings["core"]["sigma"][1])

            if settings["core"]["material"] != 'custom':
                # a custom core does not need a material, measurement_setup and _datatype
                settings["core"]["material"] = Material(settings["core"]["material"])
                settings["core"]["permeability_measurement_setup"] = MeasurementSetup(settings["core"]["permeability_measurement_setup"])
                settings["core"]["permeability_datatype"] = MeasurementDataType(settings["core"]["permeability_datatype"])
                settings["core"]["permittivity_measurement_setup"] = MeasurementSetup(settings["core"]["permittivity_measurement_setup"])
                settings["core"]["permittivity_datatype"] = MeasurementDataType(settings["core"]["permittivity_datatype"])

            settings["core"]["permeability_datasource"] = MaterialDataSource(settings["core"]["permeability_datasource"])
            settings["core"]["permittivity_datasource"] = MaterialDataSource(settings["core"]["permittivity_datasource"])

            core = Core(core_dimensions=core_dimensions, **settings["core"])
            geo.set_core(core)

            if "air_gaps" in settings:
                air_gaps = AirGaps(AirGapMethod[settings["air_gaps"]["method"]], core)
                for air_gap in settings["air_gaps"]["air_gaps"]:
                    air_gaps.add_air_gap(AirGapLegPosition[air_gap["leg_position"]], air_gap["height"],
                                         air_gap["position_value"], )
                geo.set_air_gaps(air_gaps)

            if "insulation" in settings:
                insulation = Insulation()
                insulation.add_core_insulations(*settings["insulation"]["core_insulations"])
                insulation.add_winding_insulations(settings["insulation"]["inner_winding_insulations"])
                geo.set_insulation(insulation)

            if "stray_path" in settings:
                stray_path = StrayPath(**settings["stray_path"])
                geo.set_stray_path(stray_path)

            new_virtual_winding_windows = []
            winding_windows = settings["winding_windows"]
            for winding_window in winding_windows:
                virtual_winding_windows = winding_window["virtual_winding_windows"]
                for vww in virtual_winding_windows:
                    turns = vww["turns"]
                    conductors = []
                    for winding in vww["windings"]:
                        winding_number = winding["winding_number"]
                        conductor = Conductor(winding["winding_number"], Conductivity[winding["conductivity"]])
                        conductor_type = ConductorType[winding["conductor_type"]]
                        if conductor_type == ConductorType.RectangularSolid:
                            conductor.set_rectangular_conductor(winding["thickness"])
                        elif conductor_type == ConductorType.RoundLitz:
                            # 3 of 4 wire preferences are allowed, so fill-factor is set to None, even the value is known from the log.
                            conductor.set_litz_round_conductor(winding["conductor_radius"], winding["number_strands"],
                                                               winding["strand_radius"], None,
                                                               ConductorArrangement[winding["conductor_arrangement"]])
                        elif conductor_type == ConductorType.RoundSolid:
                            conductor.set_solid_round_conductor(winding["conductor_radius"],
                                                                ConductorArrangement[winding["conductor_arrangement"]])
                        else:
                            raise Exception(f"Unknown conductor type {conductor_type.name}")

                        conductors.append(conductor)

                    new_vww = VirtualWindingWindow(vww["bot_bound"], vww["top_bound"], vww["left_bound"],
                                                   vww["right_bound"])
                    winding_type = WindingType[vww["winding_type"]]
                    if winding_type == WindingType.Single:
                        print(f"Winding Type Single")
                        winding_scheme = WindingScheme[vww["winding_scheme"]] if vww["winding_scheme"] is not None else None
                        wrap_para_type = WrapParaType[vww["wrap_para"]] if vww["wrap_para"] is not None else None
                        new_vww.set_winding(conductors[0], turns[winding_number], winding_scheme, wrap_para_type)
                        print(turns[0])
                    elif winding_type == WindingType.TwoInterleaved:
                        new_vww.set_interleaved_winding(conductors[0], turns[0], conductors[1], turns[1],
                                                        InterleavedWindingScheme[vww["winding_scheme"]],
                                                        vww["winding_insulation"])
                    else:
                        raise Exception(f"Winding type {winding_type} is not implemented")
                    new_virtual_winding_windows.append(new_vww)

            winding_window = WindingWindow(core, insulation)
            winding_window.virtual_winding_windows = new_virtual_winding_windows
            geo.set_winding_windows([winding_window])

            return geo

        raise Exception(f"Couldn't extract settings from file {log_file_path}")




