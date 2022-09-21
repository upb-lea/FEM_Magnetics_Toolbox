# Usual Python libraries
import csv
import fileinput
import numpy as np
import os
import sys
import gmsh
import json
import warnings
import inspect

from onelab import onelab
from matplotlib import pyplot as plt
from scipy.optimize import brentq
from scipy.integrate import quad
from scipy.interpolate import interp1d
from datetime import datetime
from typing import List, Union, Optional, Dict
from .thermal.thermal_simulation import *
from .thermal.thermal_functions import *
from .femmt_functions import *
from .electro_magnetic.Analytical_Core_Data import *
from .femmt_model_classes import *
from .femmt_enumerations import *

#  ===== Main Class  =====
class MagneticComponent:
    """
    A MagneticComponent is the main object for all simulation purposes in femmt.

        - One or more "MagneticComponents" can be created
        - Each "MagneticComponent" owns its own instance variable values

    """
    # Initialization of all class variables
    # Common variables for all instances

    onelab_folder_path = None

    def __init__(self, component_type: ComponentType = ComponentType.Inductor, working_directory = None):
        """
        :param component_type: Available options:
                               - "inductor"
                               - "transformer"
                               - "integrated_transformer" (Transformer with included stray-path)
        :type component_type: string
        :param working_directory: Sets the working directory
        :type working_directory: string
        """
        print(f"\n"
              f"Initialized a new Magnetic Component of type {component_type}\n"
              f"--- --- --- ---")

        # Get caller filepath when no working_directpry was set
        if working_directory is None:
            caller_filename = inspect.stack()[1].filename 
            working_directory = os.path.join(os.path.dirname(caller_filename), "femmt")

        if not os.path.exists(working_directory):
            os.mkdir(working_directory)
 
        self.update_paths(working_directory)

        self.correct_outer_leg = False
        
        # Initialization of all instance variables

        # Breaking variable
        self.valid = True

        # To make sure femm is only imported once
        self.femm_is_imported = False

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Component Geometry
        self.component_type = component_type  # "inductor", "transformer", "integrated_transformer" (or "three-phase-transformer")
        self.dimensionality = "2D"  # "2D" or "3D" # TODO As Enum? Is this even needed?
        self.symmetry = "radial"  # "radial", "linear", None # TODO As Enum? Is this even needed?

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Components
        self.core = None
        self.air_gaps = None
        self.windings = None
        self.isolation = None
        self.virtual_winding_windows = None
        self.stray_path = None

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Control Flags
        self.region = None  # Apply an outer Region or directly apply a constraint on the Core Boundary
        self.plot_fields = "standard"  # can be "standard" or False

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Windings
        self.n_windings = None  # Number of conductors/windings
        if component_type == ComponentType.Inductor:
            self.n_windings = 1
        elif component_type == ComponentType.Transformer or component_type == ComponentType.IntegratedTransformer:
            self.n_windings = 2
        else:
            raise Exception("Unknown component type")
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Virtual Winding Windows
        self.virtual_winding_windows = None
        self.vw_type = None  # "center" and "full_window" are the only cases implemented yet; #TODO: replace

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Geometric Parameters/Coordinates
        self.n_windows = None

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Excitation Parameters
        self.imposed_reduced_frequency = None
        self.flag_excitation_type = None
        self.current = [None] * self.n_windings  # Defined for every conductor
        self.current_density = [None] * self.n_windings  # Defined for every conductor
        self.voltage = [None] * self.n_windings  # Defined for every conductor
        self.frequency = None
        self.phase_deg = np.zeros(self.n_windings)  # Default is zero, Defined for every conductor
        self.red_freq = None  # [] * self.n_windings  # Defined for every conductor
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
        # Meshing
        self.mesh = self.Mesh(self)

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

        self.onelab_setup()
        self.onelab_client = onelab.client(__file__)

    def create_folders(self, *args):
        for folder in list(args):
            if not os.path.exists(folder):
                os.mkdir(folder)

    def update_paths(self, working_directory):
        # Setup folder paths 
        self.working_directory = working_directory
        self.femmt_folder_path = os.path.dirname(__file__)
        self.mesh_folder_path = os.path.join(self.working_directory, "mesh")
        self.electro_magnetic_folder_path = os.path.join(self.femmt_folder_path, "electro_magnetic")
        self.results_folder_path = os.path.join(self.working_directory, "results")
        self.e_m_values_folder_path = os.path.join(self.results_folder_path, "values")
        self.e_m_fields_folder_path = os.path.join(self.results_folder_path, "fields")
        self.e_m_circuit_folder_path = os.path.join(self.results_folder_path, "circuit")
        self.e_m_strands_coefficients_folder_path = os.path.join(self.electro_magnetic_folder_path, "Strands_Coefficients")
        self.femm_folder_path = os.path.join(self.working_directory, "femm")
        self.reluctance_model_folder_path = os.path.join(self.working_directory, "reluctance_model")
        self.thermal_results_folder_path = os.path.join(self.results_folder_path, "thermal")

        # Setup file paths
        self.e_m_results_log_path = os.path.join(self.results_folder_path, "log_electro_magnetic.json")
        self.femm_results_log_path = os.path.join(self.femm_folder_path, "result_log_femm.json")
        self.config_path = os.path.join(self.femmt_folder_path, "config.json")
        self.e_m_mesh_file = os.path.join(self.mesh_folder_path, "electro_magnetic.msh")
        self.model_geo_file = os.path.join(self.mesh_folder_path, "model.geo_unrolled")
        self.hybrid_mesh_file = os.path.join(self.mesh_folder_path, "hybrid.msh")
        self.hybrid_color_mesh_file = os.path.join(self.mesh_folder_path, "hybrid_color.msh")
        self.hybrid_color_visualize_file = os.path.join(self.mesh_folder_path, "hybrid_color.png")
        self.thermal_mesh_file = os.path.join(self.mesh_folder_path, "thermal.msh")

        # Create necessary folders
        self.create_folders(self.femmt_folder_path, self.mesh_folder_path, self.electro_magnetic_folder_path, 
            self.results_folder_path, self.e_m_values_folder_path, self.e_m_fields_folder_path, 
            self.e_m_circuit_folder_path, self.e_m_strands_coefficients_folder_path)

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   -  -  -  -  -  -  -  -  -  -  -
    # Thermal simulation
    def calculate_core_volume(self) -> float:
        # TODO core_h and core_w should always be set
        if self.core.core_h is not None and self.core.core_w is not None:
            core_height = self.core.core_h
            core_width = self.core.core_w
        else:
            core_height = self.core.window_h + self.core.core_w / 2
            core_width = self.two_d_axi.r_outer
        winding_height = self.core.window_h
        winding_width = self.core.window_w

        air_gap_volume = 0
        inner_leg_width = self.two_d_axi.r_inner - winding_width
        for leg_position, position_value, height in self.air_gaps.midpoints:
            width = 0

            if leg_position == AirGapLegPosition.LeftLeg.value:
                # left leg
                # TODO this is wrong since the airgap is not centered on the y axis 
                width = core_width - self.two_d_axi.r_inner
            elif leg_position == AirGapLegPosition.CenterLeg.value:
                # center leg
                width = inner_leg_width
            elif leg_position == AirGapLegPosition.RightLeg.value:
                # right leg
                # TODO this is wrong since the airgap is not centered on the y axis
                width = core_width - self.two_d_axi.r_inner
            else:
                raise Exception(f"Unvalid leg position tag {leg_position} used for an air gap.")

            air_gap_volume += np.pi * width**2 * height

        return np.pi*(core_width**2 * core_height - (inner_leg_width+winding_width)**2 * winding_height + inner_leg_width**2 * winding_height) - air_gap_volume

    def get_wire_distances(self):
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

    # Start thermal simulation
    def thermal_simulation(self, thermal_conductivity, boundary_temperatures, boundary_flags, case_gap_top,
                           case_gap_right, case_gap_bot, show_results=True, visualize_before=False, color_scheme: Dict = colors_femmt_default,
                           colors_geometry: Dict = colors_geometry_femmt_default) -> None:
        """
        
        Starts the thermal simulation using thermal_simulation.py

        :param color_scheme: colorfile (definition for red, green, blue, ...)
        :type color_scheme: Dict
        :param colors_geometry: definition for e.g. core is grey, winding is orange, ...
        :type colors_geometry: Dict

        :return: -
        """
        # Create necessary folders
        self.create_folders(self.thermal_results_folder_path)

        self.mesh.generate_thermal_mesh(case_gap_top, case_gap_right, case_gap_bot, color_scheme, colors_geometry, visualize_before)

        if not os.path.exists(self.e_m_results_log_path):
            # Simulation results file not created
            raise Exception("Cannot run thermal simulation -> Magnetic simulation needs to run first (no results_log.json found")

        # Check if the results log path simulation settings fit the current simulation settings
        current_settings = encode_settings(self)
        del current_settings["working_directory"]
        del current_settings["date"]

        log_settings = None
        with open(self.e_m_results_log_path, "r") as fd:
            content = json.load(fd)
            log_settings = content["simulation_settings"]
        del log_settings["working_directory"]
        del log_settings["date"]
        
        if current_settings != log_settings:
            raise Exception(f"The settings from the log file {self.e_m_results_log_path} do not match the current simulation settings. \
                                Please re-run the magnetic simulation.")

        tags = {
            "core_tag": self.mesh.ps_core,
            "background_tag": self.mesh.ps_air,
            "winding_tags": self.mesh.ps_cond,
            "air_gaps_tag": self.mesh.ps_air_gaps if self.air_gaps.number > 0 else None,
            "boundary_regions": self.mesh.thermal_boundary_region_tags,
            "isolations_tag": self.mesh.ps_isolation
        }

        # Core area -> Is needed to estimate the heat flux
        # Power density for volumes W/m^3
        core_area = self.calculate_core_volume()

        # Set wire radii
        wire_radii = [winding.conductor_radius for winding in self.windings]

        thermal_parameters = {
            "onelab_folder_path": self.onelab_folder_path,
            "model_mesh_file_path": self.thermal_mesh_file,
            "results_log_file_path": self.e_m_results_log_path,
            "results_folder_path": self.thermal_results_folder_path,
            "tags_dict": tags,
            "thermal_conductivity_dict": thermal_conductivity,
            "boundary_temperatures": boundary_temperatures,
            "boundary_flags": boundary_flags,
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
            "show_results": show_results,
            "print_sensor_values": False
        }

        run_thermal(**thermal_parameters)

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   -  -  -  -  -  -  -  -  -  -  -
    # Setup
    def onelab_setup(self) -> None:
        """

        Either reads ONELAB parent folder path from config.json or asks the user to provide the ONELAB path it.
        Creates a config.json inside the site-packages folder at first run.

        :return: -
        """
        # check if config.json is available and not empty
        if os.path.isfile(self.config_path) and os.stat(self.config_path).st_size != 0:
            onelab_path = ""
            with open(self.config_path, "r") as fd:
                loaded_dict = json.loads(fd.read())
                onelab_path = loaded_dict['onelab']

            if os.path.exists(onelab_path) and os.path.isfile(os.path.join(onelab_path, "onelab.py")):
                # Path found
                self.onelab_folder_path = onelab_path
                return

        # Let the user enter the onelab_path:
        # Find out the onelab_path of installed module, or in case of running directly from git, find the onelab_path of git repository
        # loop until path is correct
        onelab_path_wrong = True
        path_wrong = True
        while onelab_path_wrong:
            onelab_path = os.path.normpath(input("Enter the path of onelabs parent folder (path to folder which contains getdp, onelab executables): "))

            if os.path.exists(onelab_path):
                onelab_path_wrong = False
                break
            else:
                print('onelab not found! Tool searches for onelab.py in the folder. Please re-enter path!')
        self.onelab_folder_path = onelab_path

        # Write the path to the config.json
        onelab_path_dict = {"onelab": onelab_path}
        with open(os.path.join(self.config_path), 'w', encoding='utf-8') as fd:
            json.dump(onelab_path_dict, fd, indent=2, ensure_ascii=False)

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   -  -  -  -  -  -  -  -  -  -  -
    # Geometry Parts
    def high_level_geo_gen(self, dimensionality="2D", frequency=None, skin_mesh_factor=None, isolation_deltas = None):
        """
        - high level geometry generation
        - based on chosen core and conductor types and simulation mode
        - calls "low level" methods, that create all points needed for mesh generation
     
        :return:
     
        """
        # Always reset the to valid
        self.valid = True

        # Mesh-Parameters must be updated depending on geometry size
        self.mesh.c_core = self.core.core_w / 10. * self.mesh.global_accuracy
        self.mesh.c_window = self.core.window_w / 30 * self.mesh.global_accuracy

        # Update Skin Depth (needed for meshing)
        self.mesh.skin_mesh_factor = skin_mesh_factor
        if frequency is not None:
            if frequency == 0:
                self.delta = 1e9
            else:
                self.delta = np.sqrt(2 / (2 * frequency * np.pi * self.windings[0].cond_sigma * self.mu0))
            for i in range(0, self.n_windings):
                if self.windings[i].conductor_type == ConductorType.Solid:
                    self.mesh.c_conductor[i] = min([self.delta * self.mesh.skin_mesh_factor, self.windings[i].conductor_radius / 4 * self.mesh.global_accuracy]) #* self.mesh.skin_mesh_factor])
                    self.mesh.c_center_conductor[i] = self.windings[i].conductor_radius / 4 * self.mesh.global_accuracy  # * self.mesh.skin_mesh_factor
                elif self.windings[i].conductor_type == ConductorType.Litz:
                    self.mesh.c_conductor[i] = self.windings[i].conductor_radius / 4 * self.mesh.global_accuracy
                    self.mesh.c_center_conductor[i] = self.windings[i].conductor_radius / 4 * self.mesh.global_accuracy
                else:
                    self.mesh.c_conductor[i] = 0.0001  # TODO: dynamic implementation

        # -- Core-type --
        self.n_windows = 2

        # -- Symmetry -- [Choose between asymmetric, symmetric and axi symmetric]
        self.dimensionality = dimensionality
        if self.dimensionality == "2D":
            self.two_d_axi = self.TwoDaxiSymmetric(self)
            self.two_d_axi.update(isolation_deltas)

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   -  -  -  -  -  -  -  -  -  -  -
    # Create Model
    def set_isolation(self, isolation: Isolation):
        if isolation.cond_cond is None or not isolation.cond_cond:
            raise Exception("Isolations between the conductors must be set")

        if isolation.core_cond is None or not isolation.core_cond:
            raise Exception("Isolations between the core and the conductors must be set")

        self.isolation = isolation

    def set_stray_path(self, stray_path: StrayPath):
        self.stray_path = stray_path

    def set_air_gaps(self, air_gaps: AirGaps):
        # Sorting air gaps from lower to upper
        air_gaps.midpoints.sort(key=lambda x: x[1])

        self.air_gaps = air_gaps

    def set_windings(self, windings: List[Winding]):
        self.windings = windings

        self.virtual_winding_windows = []

        if self.component_type == ComponentType.Inductor:
            if len(windings) != 1:
                raise Exception("Inductor was set but the number of windings is not 1")

            self.vw_type = VirtualWindingType.FullWindow

        if self.component_type == ComponentType.Transformer:
            if len(self.windings) != 2:
                raise Exception("Only transformers with 2 windings are allowed")

            if self.windings[0].winding_type == WindingType.Interleaved and self.windings[1].winding_type == WindingType.Interleaved:
                self.vw_type = VirtualWindingType.FullWindow
            elif self.windings[0].winding_type == WindingType.Interleaved and self.windings[1].winding_type != WindingType.Interleaved or \
                    self.windings[0].winding_type != WindingType.Interleaved and self.windings[1].winding_type == WindingType.Interleaved:
                raise Exception("When interleaved scheme is needed. Please set both winding types to interleaved")
            else:
                self.vw_type = VirtualWindingType.Split2

        if self.component_type == ComponentType.IntegratedTransformer:
            if len(self.windings) != 2:
                raise Exception("Only integrated transformers with 2 windings are allowed")

            if self.windings[0].winding_type == WindingType.Interleaved and self.windings[1].winding_type == WindingType.Interleaved:
                self.windings[0].turns = [self.windings[0].turns_primary, self.windings[0].turns_secondary]
                self.windings[1].turns = [self.windings[1].turns_primary, self.windings[1].turns_secondary]

            self.vw_type = VirtualWindingType.Split2

        if self.vw_type == VirtualWindingType.FullWindow:
            self.virtual_winding_windows = [VirtualWindingWindow(self.windings[0].winding_type, self.windings[0].winding_scheme)]
        elif self.vw_type == VirtualWindingType.Split2:
            vww1 = VirtualWindingWindow(self.windings[0].winding_type, self.windings[0].winding_scheme)
            vww2 = VirtualWindingWindow(self.windings[1].winding_type, self.windings[1].winding_scheme)
            self.virtual_winding_windows = [vww1, vww2]
        else:
            # Currently there is no third or default VirtualWindingType
            pass

    def set_core(self, core: Core):
        self.core = core

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   -  -  -  -  -  -  -  -  -  -  -
    # Full Geometry
    class TwoDaxiSymmetric:
        """
        - creates all points needed for the radial axi-symmetric core topology

        :return:

        """

        def __init__(self, component):
            self.component = component
            # -- Arrays for geometry data -- [all points with (x, y, z, mesh_accuracy)]
            self.p_outer = None  # np.zeros((4, 4))
            self.p_region_bound = None  # np.zeros((4, 4))
            self.p_window = None  # np.zeros((4 * self.component.n_windows, 4))
            self.p_air_gaps = None  # np.zeros((4 * self.component.air_gaps.number, 4))
            self.p_conductor = []
            self.p_iso_core = []
            self.p_iso_pri_sec = []
            for i in range(0, self.component.n_windings):
                self.p_conductor.insert(i, [])

            self.r_inner = None
            self.r_outer = None

        def draw_outer(self):
            """
            Draws the outer points

            :return:
            """
            # Outer Core
            # (A_zyl=2pi*r*h => h=0.5r=0.25core_w <=> ensure A_zyl=A_core on the tiniest point)
            # TODO Case core_h is not None
            if self.component.core.core_h is None:
                self.p_outer[0][:] = [-self.r_outer,
                                    -(self.component.core.window_h / 2 + self.component.core.core_w / 4),
                                    0,
                                    self.component.mesh.c_core]

                self.p_outer[1][:] = [self.r_outer,
                                    -(self.component.core.window_h / 2 + self.component.core.core_w / 4),
                                    0,
                                    self.component.mesh.c_core]

                self.p_outer[2][:] = [-self.r_outer,
                                    (self.component.core.window_h / 2 + self.component.core.core_w / 4),
                                    0,
                                    self.component.mesh.c_core]

                self.p_outer[3][:] = [self.r_outer,
                                    (self.component.core.window_h / 2 + self.component.core.core_w / 4),
                                    0,
                                    self.component.mesh.c_core]
            else:
                self.p_outer[0][:] = [-self.r_outer,
                                    -self.component.core.core_h/2,
                                    0,
                                    self.component.mesh.c_core]

                self.p_outer[1][:] = [self.r_outer,
                                    -self.component.core.core_h/2,
                                    0,
                                    self.component.mesh.c_core]

                self.p_outer[2][:] = [-self.r_outer,
                                    self.component.core.core_h/2,
                                    0,
                                    self.component.mesh.c_core]

                self.p_outer[3][:] = [self.r_outer,
                                    self.component.core.core_h/2,
                                    0,
                                    self.component.mesh.c_core]

        def draw_window(self):
            # Window
            # At this point both windows (in a cut) are modeled
            # print(f"win: c_window: {self.component.mesh.c_window}")
            self.p_window[0] = [-self.r_inner,
                                -self.component.core.window_h / 2,
                                0,
                                self.component.mesh.c_window]

            self.p_window[1] = [-self.component.core.core_w / 2,
                                -self.component.core.window_h / 2,
                                0,
                                self.component.mesh.c_window]

            self.p_window[2] = [-self.r_inner,
                                self.component.core.window_h / 2,
                                0,
                                self.component.mesh.c_window]

            self.p_window[3] = [-self.component.core.core_w / 2,
                                self.component.core.window_h / 2,
                                0,
                                self.component.mesh.c_window]

            self.p_window[4] = [self.component.core.core_w / 2,
                                -self.component.core.window_h / 2,
                                0,
                                self.component.mesh.c_window]

            self.p_window[5] = [self.r_inner,
                                -self.component.core.window_h / 2,
                                0,
                                self.component.mesh.c_window]

            self.p_window[6] = [self.component.core.core_w / 2,
                                self.component.core.window_h / 2,
                                0,
                                self.component.mesh.c_window]

            self.p_window[7] = [self.r_inner,
                                self.component.core.window_h / 2,
                                0,
                                self.component.mesh.c_window]

        def draw_air_gaps(self):
            # Air gaps
            # "air_gaps" is a list with [position_tag, air_gap_position, air_gap_h, c_air_gap]
            #   - position_tag: specifies the gapped "leg"
            #   - air_gap_position: specifies the coordinate of the air gap's center point along the specified leg
            #   - air_gap_h: height/length of the air gap
            #   - c_air_gap: mesh accuracy factor
            # at this point the 4 corner points of each air gap are generated out of "air_gaps"

            mesh_accuracy = self.component.core.window_w / 20 * self.component.mesh.global_accuracy

            for i in range(0, self.component.air_gaps.number):

                # # Left leg (-1)
                # if self.component.air_gaps.midpoints[i][0] == -1:
                #     self.p_air_gaps[i * 4] = [-(self.component.core.core_w + self.component.core.window_w),
                #                               self.component.air_gaps.midpoints[i][1] -
                #     self.component.air_gaps.midpoints[i][2] / 2, 0, self.component.air_gaps.midpoints[i][3]]
                #     self.p_air_gaps[i * 4 + 1] = [-(self.component.core.core_w / 2 + self.component.core.window_w),
                #                                   self.component.air_gaps.midpoints[i][1] -
                #     self.component.air_gaps.midpoints[i][2] / 2, 0, self.component.air_gaps.midpoints[i][3]]
                #     self.p_air_gaps[i * 4 + 2] = [-(self.component.core.core_w + self.component.core.window_w),
                #                                   self.component.air_gaps.midpoints[i][1] +
                #     self.component.air_gaps.midpoints[i][2] / 2, 0, self.component.air_gaps.midpoints[i][3]]
                #     self.p_air_gaps[i * 4 + 3] = [-(self.component.core.core_w / 2 + self.component.core.window_w),
                #                                   self.component.air_gaps.midpoints[i][1] +
                #     self.component.air_gaps.midpoints[i][2] / 2, 0, self.component.air_gaps.midpoints[i][3]]
                #
                # # Right leg (+1)
                # if self.component.air_gaps.midpoints[i][0] == 1:
                #     self.p_air_gaps[i * 4] = [self.component.core.core_w / 2 + self.component.core.window_w,
                #                               self.component.air_gaps.midpoints[i][1] -
                #     self.component.air_gaps.midpoints[i][2] / 2, 0, self.component.air_gaps.midpoints[i][3]]
                #     self.p_air_gaps[i * 4 + 1] = [self.component.core.core_w + self.component.core.window_w,
                #                                   self.component.air_gaps.midpoints[i][1] -
                #     self.component.air_gaps.midpoints[i][2] / 2, 0, self.component.air_gaps.midpoints[i][3]]
                #     self.p_air_gaps[i * 4 + 2] = [self.component.core.core_w / 2 + self.component.core.window_w,
                #                                   self.component.air_gaps.midpoints[i][1] +
                #     self.component.air_gaps.midpoints[i][2] / 2, 0, self.component.air_gaps.midpoints[i][3]]
                #     self.p_air_gaps[i * 4 + 3] = [self.component.core.core_w + self.component.core.window_w,
                #                                   self.component.air_gaps.midpoints[i][1] +
                #     self.component.air_gaps.midpoints[i][2] / 2, 0, self.component.air_gaps.midpoints[i][3]]

                # Center leg (0)
                if self.component.air_gaps.midpoints[i][0] == 0:
                    # The center points are transformed each into 4 corner points

                    air_gap_y_position = self.component.air_gaps.midpoints[i][1]
                    air_gap_height = self.component.air_gaps.midpoints[i][2]
                    air_gap_length_top = self.component.core.core_w / 2
                    air_gap_length_bot = self.component.core.core_w / 2

                    # Check for stray_paths in integrated transformers
                    if self.component.component_type == ComponentType.IntegratedTransformer:
                        if self.component.stray_path.start_index == i:
                            # Stray path is above current air_gap
                            air_gap_length_top = self.component.stray_path.length
                        elif self.component.stray_path.start_index + 1 == i:
                            # Stray path is below current air_gap
                            air_gap_length_bot = self.component.stray_path.length

                    # Bottom left
                    self.p_air_gaps[i * 4 + 0] = [0,
                                                  air_gap_y_position -
                                                  air_gap_height / 2,
                                                  0,
                                                  self.component.mesh.c_core]

                    # Bottom right
                    self.p_air_gaps[i * 4 + 1] = [air_gap_length_bot,
                                                  air_gap_y_position -
                                                  air_gap_height / 2,
                                                  0,
                                                  mesh_accuracy]

                    # Top left
                    self.p_air_gaps[i * 4 + 2] = [0,
                                                  air_gap_y_position +
                                                  air_gap_height / 2,
                                                  0,
                                                  self.component.mesh.c_core]

                    # Top right
                    self.p_air_gaps[i * 4 + 3] = [air_gap_length_top,
                                                  air_gap_y_position +
                                                  air_gap_height / 2,
                                                  0,
                                                  mesh_accuracy]

            # In order to close the air gap when a stray_path is added, additional points need to be added
            if self.component.component_type == ComponentType.IntegratedTransformer:
                top_point = [self.component.core.core_w / 2,
                                self.component.air_gaps.midpoints[self.component.stray_path.start_index+1][1] -
                                air_gap_height / 2,
                                0,
                                mesh_accuracy]
                bot_point = [self.component.core.core_w / 2,
                                self.component.air_gaps.midpoints[self.component.stray_path.start_index][1] +
                                air_gap_height / 2,
                                0,
                                mesh_accuracy]
                self.p_close_air_gaps = [top_point, bot_point]

        def draw_virtual_winding_windows(self):
            # Virtual Windows
            # TODO: make this part of the class VWW...
            #  self.component.vw_type must be an input or so to that class
            separation_hor = 0  # self.component.core.window_h * 0.5
            separation_vert = self.component.core.window_w * 0.5

            if not self.component.component_type == ComponentType.IntegratedTransformer:
                # Some examples for virtual windows
                # Concentrated windings

                if self.component.vw_type == VirtualWindingType.FullWindow:
                    """
                    The winding window is completely used by one VWW.
                    In case of a transformer, an interleaved winding scheme must be used.
                    """
                    # top window
                    minimum = -self.component.core.window_h / 2 + self.component.isolation.core_cond[0]  # bottom
                    maximum = self.component.core.window_h / 2 - self.component.isolation.core_cond[1]  # top
                    left = self.component.core.core_w / 2 + self.component.isolation.core_cond[2]
                    right = self.r_inner - self.component.isolation.core_cond[3]

                    # Sum the windows up in a list
                    # self.virtual_windows = [[min, max, left, right]]
                    self.component.virtual_winding_windows[0].bot_bound = minimum
                    self.component.virtual_winding_windows[0].top_bound = maximum
                    self.component.virtual_winding_windows[0].left_bound = left
                    self.component.virtual_winding_windows[0].right_bound = right

                if self.component.vw_type == VirtualWindingType.Split2:
                    """
                    The winding window is split into two VWWs.
                    The primary winding is placed in the upper half,
                    the secondary winding is placed in the lower half of the winding window
                    """
                    separation_hor = 0

                    # top window
                    min21 = -separation_hor + self.component.isolation.cond_cond[-1] / 2  # separation_hor
                    max21 = self.component.core.window_h / 2 - self.component.isolation.core_cond[0]  # top
                    left21 = self.component.core.core_w / 2 + self.component.isolation.core_cond[2]
                    right21 = self.r_inner - self.component.isolation.core_cond[3]

                    # bottom window
                    min11 = -self.component.core.window_h / 2 + self.component.isolation.core_cond[1]  # bottom
                    max11 = -separation_hor - self.component.isolation.cond_cond[-1] / 2  # separation_hor
                    left11 = self.component.core.core_w / 2 + self.component.isolation.core_cond[2]
                    right11 = self.r_inner - self.component.isolation.core_cond[3]

                    # Sum the windows up in a list
                    virtual_windows = [[min11, max11, left11, right11],
                                       [min21, max21, left21, right21]]
                    for vww in range(0, len(virtual_windows)):
                        self.component.virtual_winding_windows[vww].bot_bound = virtual_windows[vww][0]
                        self.component.virtual_winding_windows[vww].top_bound = virtual_windows[vww][1]
                        self.component.virtual_winding_windows[vww].left_bound = virtual_windows[vww][2]
                        self.component.virtual_winding_windows[vww].right_bound = virtual_windows[vww][3]

                # Because of the switch to Enumerations this is currently not used
                if self.component.vw_type == "something_else":
                    # bottom left window
                    min11 = -self.component.core.window_h / 2 + self.component.isolation.core_cond[0]  # bottom
                    max11 = -separation_hor - self.component.isolation.cond_cond[-1] / 2  # separation_hor
                    left11 = self.component.core.core_w / 2 + self.component.isolation.core_cond[0]
                    right11 = self.r_inner - self.component.isolation.cond_cond[0] - separation_vert

                    # bottom right window
                    min12 = -self.component.core.window_h / 2 + self.component.isolation.core_cond[0]  # bottom
                    max12 = -separation_hor - self.component.isolation.cond_cond[-1] / 2  # separation_hor
                    left12 = self.r_inner + self.component.isolation.cond_cond[0] - separation_vert
                    right12 = self.r_inner - self.component.isolation.core_cond[0]

                    # top window
                    min21 = -separation_hor + self.component.isolation.cond_cond[-1] / 2  # separation_hor
                    max21 = self.component.core.window_h / 2 - self.component.isolation.core_cond[0]  # top
                    left21 = self.component.core.core_w / 2 + self.component.isolation.core_cond[0]
                    right21 = self.r_inner - self.component.isolation.core_cond[0]

                    # Sum the windows up in a list
                    virtual_windows = [[min11, max11, left11, right11],
                                       [min12, max12, left12, right12],
                                       [min21, max21, left21, right21]]
                    # TODO: More flexible virtual winging windows

            # With dedicated stray path:
            if self.component.component_type == ComponentType.IntegratedTransformer:
                """
                If dedicated stray path is the chosen typology the are two winding windows
                These can either be split up into more virtual windows or (in case of bifilar windings) not
        
                """
                # TODO: Separation in more Virtual Winding Windows

                # bot window
                min11 = -self.component.core.window_h / 2 + self.component.isolation.core_cond[1]  # bottom
                max11 = self.p_air_gaps[self.component.stray_path.start_index*4+2][1] - self.component.isolation.core_cond[0]  # sep_hor
                left11 = self.component.core.core_w / 2 + self.component.isolation.core_cond[2]
                right11 = self.r_inner - self.component.isolation.core_cond[3]

                # top window
                min21 = self.p_air_gaps[(self.component.stray_path.start_index+1)*4][1] + self.component.isolation.core_cond[1]
                max21 = self.component.core.window_h / 2 - self.component.isolation.core_cond[0]  # top
                left21 = self.component.core.core_w / 2 + self.component.isolation.core_cond[2]
                right21 = self.r_inner - self.component.isolation.core_cond[3]

                # Store the window boarders in the VWW objects
                virtual_windows = [[min21, max21, left21, right21], [min11, max11, left11, right11]]
                for vww in range(0, len(virtual_windows)):
                    self.component.virtual_winding_windows[vww].bot_bound = virtual_windows[vww][0]
                    self.component.virtual_winding_windows[vww].top_bound = virtual_windows[vww][1]
                    self.component.virtual_winding_windows[vww].left_bound = virtual_windows[vww][2]
                    self.component.virtual_winding_windows[vww].right_bound = virtual_windows[vww][3]

        def draw_conductors(self):
            # Conductors
            for n_win in range(0, len(self.component.virtual_winding_windows)):
                """
                - Work through the virtual winding windows
                - Usual cases are one window for classical transformers or two windows for transformers 
                   with a dedicated stray path
                - There can be as many virtual winding windows as the user wants to define...
                - To automatically fill a virtual winding window with windings, #TODO: self.interleaving[n_win] 
                   can be chosen to "bifilar", "vertical", "horizontal", ["hexa", "square"] or completely with 
                   one of the windings by "primary" or "secondary"
        
                """
                # Boarders of the VWW:
                # start with the lower one
                bot_bound = self.component.virtual_winding_windows[n_win].bot_bound
                top_bound = self.component.virtual_winding_windows[n_win].top_bound
                left_bound = self.component.virtual_winding_windows[n_win].left_bound
                right_bound = self.component.virtual_winding_windows[n_win].right_bound

                if self.component.virtual_winding_windows[n_win].winding == WindingType.Interleaved:

                    if self.component.virtual_winding_windows[n_win].scheme == WindingScheme.Bifilar:
                        """
                        - Bifilar interleaving means a uniform winding scheme of two conductors (prim. and sec.)
                        - Can only be used for conductors of identical radius (in terms of litz radius for 
                          stranded wires)
                        - Excess windings are placed below the bifilar ones
        
                        """
                        if self.component.windings[0].conductor_radius != self.component.windings[1].conductor_radius:
                            print("For bifilar winding scheme both conductors must be of the same radius!")
                        else:
                            print("Bifilar winding scheme is applied")
                            # n1 = self.n_turns[0]/self.n_turns[1]
                            # n2 = self.n_stray_turns[0]/self.n_stray_turns[1]

                            # for
                            #     if self.component.virtual_winding_windows[num].scheme == "hexa":
                            #         y = bot_bound + self.component.windings[num].conductor_radius
                            #         x = left_bound + self.component.windings[num].conductor_radius
                            #         i = 0
                            #         base_line = True
                            #         # Case n_conductors higher that "allowed" is missing
                            #         while x < right_bound - self.component.windings[num].conductor_radius and
                            #         i < self.component.windings[num].turns:
                            #             while y < top_bound - self.component.windings[num].conductor_radius and
                            #             i < self.component.windings[num].turns:
                            #                 self.p_conductor[num].append([x, y, 0,
                            #                                           self.component.mesh.c_center_conductor[num]])
                            #                 self.p_conductor[num].append(
                            #                     [x - self.component.windings[num].conductor_radius, y, 0,
                            #                 self.component.mesh.c_conductor[num]])
                            #                 self.p_conductor[num].append(
                            #                     [x, y + self.component.windings[num].conductor_radius, 0,
                            #                 self.component.mesh.c_conductor[num]])
                            #                 self.p_conductor[num].append(
                            #                     [x + self.component.windings[num].conductor_radius, y, 0,
                            #                 self.component.mesh.c_conductor[num]])
                            #                 self.p_conductor[num].append(
                            #                     [x, y - self.component.windings[num].conductor_radius, 0,
                            #                 self.component.mesh.c_conductor[num]])
                            #                 i += 1
                            #                 y += self.component.windings[num].conductor_radius * 2 + \
                            #                      self.component.isolation.cond_cond[num]
                            #                 # from bottom to top
                            #             x += 2 * np.cos(np.pi / 6) * (self.component.windings[num].conductor_radius +
                            #             self.component.isolation.cond_cond[
                            #                 num] / 2)  # * np.sqrt(2 / 3 * np.pi / np.sqrt(3))  # one step from
                            #     # left to
                            #                 # right
                            #             # depending on what line, hexa scheme starts shifted
                            #             # reset y to "new" bottom
                            #             base_line = (not base_line)
                            #             if base_line:
                            #                 y = bot_bound + self.component.windings[num].conductor_radius
                            #             else:
                            #                 y = bot_bound + 2 * self.component.windings[num].conductor_radius +
                            #                 self.component.isolation.cond_cond[num] / 2

                    if self.component.virtual_winding_windows[n_win].scheme == WindingScheme.Vertical:
                        """
                        - Vertical interleaving means a winding scheme where the two conductors are alternating 
                           in vertical (y-)direction
                        - This is practically uncommon
                        - If the turns ratio is != 1, the scheme always begins with the "higher-turns-number's" 
                           conductor
                        """

                    if self.component.virtual_winding_windows[n_win].scheme == WindingScheme.Horizontal:
                        """
                        - Horizontal interleaving means a winding scheme where the two conductors are alternating in 
                        horizontal  (x-)direction (Tonnenwicklung)
                        - This is practically most common
                        - If the turns ratio is != 1, the scheme always begins with the "higher-turns-number's" 
                          conductor
                          
                        """

                        # assume 2 winding transformer and dedicated stray path:
                        if self.component.component_type == ComponentType.IntegratedTransformer or self.component.n_windings == 2:

                            # top window
                            if n_win == 0:
                                # Initialize the list, that counts the already placed conductors
                                N_completed = [0, 0]

                                # Initialize the starting conductor
                                if self.component.windings[0].turns[n_win] >= self.component.windings[1].turns[n_win]:
                                    # Primary starts first
                                    self.col_cond_start = 0
                                else:
                                    # Secondary starts fist
                                    self.col_cond_start = 1

                                col_cond = self.col_cond_start

                                # Initialize the x and y coordinate
                                x = left_bound + self.component.windings[col_cond].conductor_radius
                                y = top_bound - self.component.windings[col_cond].conductor_radius
                                self.top_window_iso_counter = 0
                                # Continue placing as long as not all conductors have been placed
                                while (self.component.windings[0].turns[n_win] - N_completed[0] != 0) or \
                                        (self.component.windings[1].turns[n_win] - N_completed[1] != 0):
                                    if self.component.windings[col_cond].turns[n_win] - N_completed[col_cond] != 0:
                                        # is this winding not already finished?
                                        if x < right_bound - self.component.windings[col_cond].conductor_radius:
                                            while y > bot_bound + self.component.windings[col_cond].conductor_radius and \
                                                    N_completed[col_cond] < self.component.windings[col_cond].turns[n_win]:
                                                self.p_conductor[col_cond].append([x, y, 0, self.component.mesh.c_center_conductor[col_cond]])

                                                self.p_conductor[col_cond].append(
                                                    [x - self.component.windings[col_cond].conductor_radius,
                                                     y,
                                                     0,
                                                     self.component.mesh.c_conductor[col_cond]])

                                                self.p_conductor[col_cond].append([x,
                                                                                   y + self.component.windings[col_cond].conductor_radius,
                                                                                   0,
                                                                                   self.component.mesh.c_conductor[col_cond]])

                                                self.p_conductor[col_cond].append(
                                                    [x + self.component.windings[col_cond].conductor_radius,
                                                     y,
                                                     0,
                                                     self.component.mesh.c_conductor[col_cond]])

                                                self.p_conductor[col_cond].append([x,
                                                                                   y - self.component.windings[col_cond].conductor_radius,
                                                                                   0,
                                                                                   self.component.mesh.c_conductor[col_cond]])

                                                N_completed[col_cond] += 1

                                                y -= self.component.windings[col_cond].conductor_radius * 2 + self.component.isolation.cond_cond[col_cond]  # one from bot to top

                                            x += self.component.windings[col_cond].conductor_radius + \
                                                 self.component.windings[(col_cond + 1) % 2].conductor_radius + \
                                                 self.component.isolation.cond_cond[2]  # from left to right

                                            # Reset y
                                            col_cond = (col_cond + 1) % 2
                                            y = top_bound - self.component.windings[col_cond].conductor_radius
                                            self.top_window_iso_counter += 1
                                        else:
                                            break

                                    else:
                                        # is this winding already finished? - continue with the other one
                                        col_cond = (col_cond + 1) % 2

                                        # Correct the reset of y and correct x displacement
                                        x += self.component.windings[col_cond].conductor_radius - \
                                             self.component.windings[(col_cond + 1) % 2].conductor_radius \
                                             - self.component.isolation.cond_cond[2] + self.component.isolation.cond_cond[
                                                 col_cond]

                                        y = top_bound - self.component.windings[col_cond].conductor_radius
                                        self.top_window_iso_counter -= 1
                            #  bottom window
                            if n_win == 1:
                                # Initialize the list, that counts the already placed conductors
                                N_completed = [0, 0]

                                # Initialize the starting conductor
                                if self.component.windings[0].turns[n_win] >= self.component.windings[1].turns[n_win]:
                                    col_cond = 0
                                else:
                                    col_cond = 1

                                # Initialize the x and y coordinate
                                x = left_bound + self.component.windings[col_cond].conductor_radius
                                y = bot_bound + self.component.windings[col_cond].conductor_radius

                                # Continue placing as long as not all conductors have been placed
                                while (self.component.windings[0].turns[n_win] - N_completed[0] != 0) or \
                                        (self.component.windings[1].turns[n_win] - N_completed[1] != 0):
                                    if self.component.windings[col_cond].turns[n_win] - N_completed[col_cond] != 0:
                                        # is this winding not already finished?
                                        if x < right_bound - self.component.windings[col_cond].conductor_radius:
                                            while y < top_bound - self.component.windings[col_cond].conductor_radius and \
                                                    N_completed[col_cond] < self.component.windings[col_cond].turns[
                                                n_win]:
                                                self.p_conductor[col_cond].append(
                                                    [x, y, 0, self.component.mesh.c_center_conductor[col_cond]])

                                                self.p_conductor[col_cond].append(
                                                    [x - self.component.windings[col_cond].conductor_radius,
                                                     y,
                                                     0,
                                                     self.component.mesh.c_conductor[col_cond]])

                                                self.p_conductor[col_cond].append([x,
                                                                                   y + self.component.windings[
                                                                                       col_cond].conductor_radius,
                                                                                   0,
                                                                                   self.component.mesh.c_conductor[
                                                                                       col_cond]])

                                                self.p_conductor[col_cond].append(
                                                    [x + self.component.windings[col_cond].conductor_radius,
                                                     y,
                                                     0,
                                                     self.component.mesh.c_conductor[col_cond]])

                                                self.p_conductor[col_cond].append([x,
                                                                                   y - self.component.windings[
                                                                                       col_cond].conductor_radius,
                                                                                   0,
                                                                                   self.component.mesh.c_conductor[
                                                                                       col_cond]])

                                                N_completed[col_cond] += 1

                                                y += self.component.windings[col_cond].conductor_radius * 2 + \
                                                     self.component.isolation.cond_cond[col_cond]  # one from bot to top

                                            x += self.component.windings[col_cond].conductor_radius + \
                                                 self.component.windings[(col_cond + 1) % 2].conductor_radius + \
                                                 self.component.isolation.cond_cond[2]  # from left to right

                                            # Reset y
                                            col_cond = (col_cond + 1) % 2
                                            y = bot_bound + self.component.windings[col_cond].conductor_radius

                                        else:
                                            break

                                    else:
                                        # is this winding already finished? - continue with the other one
                                        col_cond = (col_cond + 1) % 2

                                        # Correct the reset of y and correct x displacement
                                        x += self.component.windings[col_cond].conductor_radius - \
                                             self.component.windings[(col_cond + 1) % 2].conductor_radius \
                                             - self.component.isolation.cond_cond[2] + \
                                             self.component.isolation.cond_cond[
                                                 col_cond]

                                        y = bot_bound + self.component.windings[col_cond].conductor_radius

                    """Blockwise concentrated"""
                    if isinstance(self.component.virtual_winding_windows[n_win].scheme, list):
                        """
                        - interleaving with a list means a concentrated winding scheme of ("hexagonal", "square" or 
                          mixed) in virtual winding window
                        - only valid for two winding case (=transformer) 
                        - vertical stacking
                        - block winding

                        how many turns fit in arrow?
                        from top to bot
                        while (not placed all cond.):
                            1. start with the primary winding from bot / left
                            2. continue with the secondary from top / right
                            3.CHECK solution conditions
                       
                        """
                        # CHECK for two winding transformer
                        if len(self.component.virtual_winding_windows[n_win].scheme) != 2:
                            print(f"Interleaving with a list is only valid for the two winding case.\n"
                                  f"Therefore the scheme must be a list of length 2 but is of length "
                                  f"{len(self.component.virtual_winding_windows[n_win].scheme)}")
                            raise Warning

                        for num in range(0, len(self.component.virtual_winding_windows[n_win].scheme)):

                            y = None

                            # Cases
                            if num == 0:
                                y = bot_bound + self.component.windings[num].conductor_radius
                            if num == 1:
                                y = top_bound - self.component.windings[num].conductor_radius

                            # Initialization
                            x = left_bound + self.component.windings[num].conductor_radius
                            i = 0

                            if self.component.virtual_winding_windows[n_win].scheme[num] == WindingScheme.Square:

                                # Primary winding from bottom to top
                                if num == 0:
                                    while y < top_bound - self.component.windings[num].conductor_radius and \
                                            i < self.component.windings[num].turns[n_win]:
                                        while x < right_bound - self.component.windings[num].conductor_radius and \
                                                i < self.component.windings[num].turns[n_win]:
                                            self.p_conductor[num].append(
                                                [x, y, 0, self.component.mesh.c_center_conductor[num]])
                                            self.p_conductor[num].append(
                                                [x - self.component.windings[num].conductor_radius, y, 0,
                                                 self.component.mesh.c_conductor[num]])
                                            self.p_conductor[num].append(
                                                [x, y + self.component.windings[num].conductor_radius, 0,
                                                 self.component.mesh.c_conductor[num]])
                                            self.p_conductor[num].append(
                                                [x + self.component.windings[num].conductor_radius, y, 0,
                                                 self.component.mesh.c_conductor[num]])
                                            self.p_conductor[num].append(
                                                [x, y - self.component.windings[num].conductor_radius, 0,
                                                 self.component.mesh.c_conductor[num]])
                                            i += 1
                                            x += self.component.windings[num].conductor_radius * 2 + \
                                                 self.component.isolation.cond_cond[
                                                     num]  # from left to right
                                        y += self.component.windings[num].conductor_radius * 2 + \
                                             self.component.isolation.cond_cond[
                                                 num]  # one step from bot to top
                                        x = left_bound + self.component.windings[
                                            num].conductor_radius  # always the same

                                # Secondary winding from top to bottom
                                if num == 1:
                                    while y > bot_bound + self.component.windings[num].conductor_radius and i < \
                                            self.component.windings[num].turns[n_win]:
                                        while x < right_bound - self.component.windings[num].conductor_radius and i < \
                                                self.component.windings[num].turns[n_win]:
                                            self.p_conductor[num].append(
                                                [x, y, 0, self.component.mesh.c_center_conductor[num]])
                                            self.p_conductor[num].append(
                                                [x - self.component.windings[num].conductor_radius, y, 0,
                                                 self.component.mesh.c_conductor[num]])
                                            self.p_conductor[num].append(
                                                [x, y + self.component.windings[num].conductor_radius, 0,
                                                 self.component.mesh.c_conductor[num]])
                                            self.p_conductor[num].append(
                                                [x + self.component.windings[num].conductor_radius, y, 0,
                                                 self.component.mesh.c_conductor[num]])
                                            self.p_conductor[num].append(
                                                [x, y - self.component.windings[num].conductor_radius, 0,
                                                 self.component.mesh.c_conductor[num]])
                                            i += 1

                                            x += self.component.windings[num].conductor_radius * 2 + \
                                                 self.component.isolation.cond_cond[
                                                     num]  # from left to right
                                        y += -(self.component.windings[num].conductor_radius * 2) - \
                                             self.component.isolation.cond_cond[
                                                 num]  # one step from bot to top
                                        x = left_bound + self.component.windings[
                                            num].conductor_radius  # always the same

                            if self.component.virtual_winding_windows[n_win].scheme[num] == WindingScheme.Hexagonal:

                                # Primary winding from bottom to top
                                if num == 0:

                                    base_line = True

                                    while y < top_bound - self.component.windings[num].conductor_radius and \
                                            i < self.component.windings[num].turns[n_win]:
                                        while x < right_bound - self.component.windings[num].conductor_radius and \
                                                i < self.component.windings[num].turns[n_win]:

                                            self.p_conductor[num].append(
                                                [x, y, 0, self.component.mesh.c_center_conductor[num]])
                                            self.p_conductor[num].append(
                                                [x - self.component.windings[num].conductor_radius, y, 0,
                                                 self.component.mesh.c_conductor[num]])
                                            self.p_conductor[num].append(
                                                [x, y + self.component.windings[num].conductor_radius, 0,
                                                 self.component.mesh.c_conductor[num]])
                                            self.p_conductor[num].append(
                                                [x + self.component.windings[num].conductor_radius, y, 0,
                                                 self.component.mesh.c_conductor[num]])
                                            self.p_conductor[num].append(
                                                [x, y - self.component.windings[num].conductor_radius, 0,
                                                 self.component.mesh.c_conductor[num]])
                                            i += 1

                                            x += 2 * np.cos(np.pi / 6) * (
                                                    self.component.windings[num].conductor_radius +
                                                    self.component.isolation.cond_cond[num] / 2)

                                            # depending on what line, hexa scheme starts shifted
                                            # reset y to "new" bottom
                                            base_line = (not base_line)
                                            if base_line:
                                                y -= (self.component.windings[num].conductor_radius +
                                                      self.component.isolation.cond_cond[num])
                                            else:
                                                y += (self.component.windings[num].conductor_radius +
                                                      self.component.isolation.cond_cond[num])

                                        # Undo last base_line reset
                                        if base_line:
                                            y += (self.component.windings[num].conductor_radius +
                                                  self.component.isolation.cond_cond[num])
                                        else:
                                            y -= (self.component.windings[num].conductor_radius +
                                                  self.component.isolation.cond_cond[num])

                                        base_line = True
                                        x = left_bound + self.component.windings[num].conductor_radius
                                        y += self.component.windings[num].conductor_radius + \
                                             self.component.isolation.cond_cond[num]

                                # Secondary winding from top to bottom
                                if num == 1:

                                    base_line = True

                                    while y > bot_bound + self.component.windings[num].conductor_radius and \
                                            i < self.component.windings[num].turns[n_win]:
                                        while x < right_bound - self.component.windings[num].conductor_radius and \
                                                i < self.component.windings[num].turns[n_win]:
                                            print(f"i: {i} "
                                                  f"x: {x} "
                                                  f"y: {y} ")

                                            self.p_conductor[num].append(
                                                [x, y, 0, self.component.mesh.c_center_conductor[num]])
                                            self.p_conductor[num].append(
                                                [x - self.component.windings[num].conductor_radius, y, 0,
                                                 self.component.mesh.c_conductor[num]])
                                            self.p_conductor[num].append(
                                                [x, y + self.component.windings[num].conductor_radius, 0,
                                                 self.component.mesh.c_conductor[num]])
                                            self.p_conductor[num].append(
                                                [x + self.component.windings[num].conductor_radius, y, 0,
                                                 self.component.mesh.c_conductor[num]])
                                            self.p_conductor[num].append(
                                                [x, y - self.component.windings[num].conductor_radius, 0,
                                                 self.component.mesh.c_conductor[num]])

                                            i += 1
                                            x += 2 * np.cos(np.pi / 6) * (
                                                    self.component.windings[num].conductor_radius +
                                                    self.component.isolation.cond_cond[num] / 2)

                                            # depending on what line, hexa scheme starts shifted
                                            # reset y to "new" bottom
                                            base_line = (not base_line)
                                            if base_line:
                                                y += (self.component.windings[num].conductor_radius +
                                                      self.component.isolation.cond_cond[num])
                                            else:
                                                y -= (self.component.windings[num].conductor_radius +
                                                      self.component.isolation.cond_cond[num])

                                        # Undo last base_line reset
                                        if base_line:
                                            y -= (self.component.windings[num].conductor_radius +
                                                  self.component.isolation.cond_cond[num])
                                        else:
                                            y += (self.component.windings[num].conductor_radius +
                                                  self.component.isolation.cond_cond[num])

                                        base_line = True
                                        x = left_bound + self.component.windings[num].conductor_radius
                                        # from top to bottom
                                        y -= (self.component.windings[num].conductor_radius +
                                              self.component.isolation.cond_cond[num])

                else:
                    # other case is non-interleaved
                    if self.component.virtual_winding_windows[n_win].winding == WindingType.Primary:
                        num = 0
                    elif self.component.virtual_winding_windows[n_win].winding == WindingType.Secondary:
                        num = 1
                    else:
                        raise Exception(f"Unusable winding type with this settings: {self.component.virtual_winding_windows[n_win].winding}" )

                    if self.component.windings[num].conductor_type == ConductorType.Full:
                        if sum(self.component.windings[num].turns) != 1:
                            print(f"For a \"full\" conductor you must choose 1 turn for each conductor!")
                        # full window conductor
                        self.p_conductor[num].append([left_bound, bot_bound, 0, self.component.mesh.c_conductor[num]])
                        self.p_conductor[num].append([right_bound, bot_bound, 0, self.component.mesh.c_conductor[num]])
                        self.p_conductor[num].append([left_bound, top_bound, 0, self.component.mesh.c_conductor[num]])
                        self.p_conductor[num].append([right_bound, top_bound, 0, self.component.mesh.c_conductor[num]])

                    if self.component.windings[num].conductor_type == ConductorType.Stacked:
                        # Stack defined number of turns and chosen thickness
                        for i in range(0, self.component.windings[num].turns[n_win]):
                            # CHECK if top bound is reached
                            if (bot_bound + (i + 1) * self.component.windings[num].thickness +
                                i * self.component.isolation.cond_cond[num]) <= top_bound:
                                # stacking from the ground
                                self.p_conductor[num].append(
                                    [left_bound, bot_bound + i * self.component.windings[num].thickness + i *
                                     self.component.isolation.cond_cond[num], 0, self.component.mesh.c_conductor[num]])
                                self.p_conductor[num].append([right_bound,
                                                              bot_bound + i * self.component.windings[
                                                                  num].thickness + i *
                                                              self.component.isolation.cond_cond[num], 0,
                                                              self.component.mesh.c_conductor[num]])
                                self.p_conductor[num].append([left_bound,
                                                              bot_bound + (i + 1) * self.component.windings[
                                                                  num].thickness + i *
                                                              self.component.isolation.cond_cond[num], 0,
                                                              self.component.mesh.c_conductor[num]])
                                self.p_conductor[num].append([right_bound,
                                                              bot_bound + (i + 1) * self.component.windings[
                                                                  num].thickness + i *
                                                              self.component.isolation.cond_cond[num], 0,
                                                              self.component.mesh.c_conductor[num]])

                    if self.component.windings[num].conductor_type == ConductorType.Foil:
                        # Wrap defined number of turns and chosen thickness
                        if self.component.windings[num].wrap_para == WrapParaType.Fixed_Thickness:
                            for i in range(0, self.component.windings[num].turns[n_win]):
                                # CHECK if right bound is reached
                                if (left_bound + (i + 1) * self.component.windings[num].thickness +
                                    i * self.component.isolation.cond_cond[num]) <= right_bound:
                                    # Foils
                                    self.p_conductor[num].append(
                                        [left_bound + i * self.component.windings[num].thickness + i *
                                         self.component.isolation.cond_cond[num],
                                         bot_bound, 0, self.component.mesh.c_conductor[num]])
                                    self.p_conductor[num].append(
                                        [left_bound + (i + 1) * self.component.windings[num].thickness + i *
                                         self.component.isolation.cond_cond[num], bot_bound, 0,
                                         self.component.mesh.c_conductor[num]])
                                    self.p_conductor[num].append(
                                        [left_bound + i * self.component.windings[num].thickness + i *
                                         self.component.isolation.cond_cond[num],
                                         top_bound, 0, self.component.mesh.c_conductor[num]])
                                    self.p_conductor[num].append(
                                        [left_bound + (i + 1) * self.component.windings[num].thickness + i *
                                         self.component.isolation.cond_cond[num], top_bound, 0,
                                         self.component.mesh.c_conductor[num]])

                        # Fill the allowed space in the Winding Window with a chosen number of turns
                        if self.component.windings[num].wrap_para == WrapParaType.Interpolate:
                            x_interpol = np.linspace(left_bound, right_bound + self.component.isolation.cond_cond[num],
                                                     self.component.windings[num].turns + 1)
                            for i in range(0, self.component.windings[num].turns):
                                # Foils
                                self.p_conductor[num].append(
                                    [x_interpol[i], bot_bound, 0, self.component.mesh.c_conductor[num]])
                                self.p_conductor[num].append(
                                    [x_interpol[i + 1] - self.component.isolation.cond_cond[num], bot_bound, 0,
                                     self.component.mesh.c_conductor[num]])
                                self.p_conductor[num].append(
                                    [x_interpol[i], top_bound, 0, self.component.mesh.c_conductor[num]])
                                self.p_conductor[num].append(
                                    [x_interpol[i + 1] - self.component.isolation.cond_cond[num], top_bound, 0,
                                     self.component.mesh.c_conductor[num]])

                    # Round Conductors:
                    if self.component.windings[num].conductor_type in [ConductorType.Litz, ConductorType.Solid]:

                        if self.component.virtual_winding_windows[num].scheme == WindingScheme.Square:
                            y = bot_bound + self.component.windings[num].conductor_radius
                            x = left_bound + self.component.windings[num].conductor_radius
                            i = 0
                            # Case n_conductors higher that "allowed" is missing
                            while x < right_bound - self.component.windings[num].conductor_radius and i < self.component.windings[num].turns[n_win]:
                                while y < top_bound - self.component.windings[num].conductor_radius and i < self.component.windings[num].turns[n_win]:
                                    self.p_conductor[num].append([x, y, 0, self.component.mesh.c_center_conductor[num]])
                                    self.p_conductor[num].append(
                                        [x - self.component.windings[num].conductor_radius, y, 0,
                                         self.component.mesh.c_conductor[num]])
                                    self.p_conductor[num].append(
                                        [x, y + self.component.windings[num].conductor_radius, 0,
                                         self.component.mesh.c_conductor[num]])
                                    self.p_conductor[num].append(
                                        [x + self.component.windings[num].conductor_radius, y, 0,
                                         self.component.mesh.c_conductor[num]])
                                    self.p_conductor[num].append(
                                        [x, y - self.component.windings[num].conductor_radius, 0,
                                         self.component.mesh.c_conductor[num]])
                                    i += 1
                                    y += self.component.windings[num].conductor_radius * 2 + self.component.isolation.cond_cond[num]  # one step from left to right
                                x += self.component.windings[num].conductor_radius * 2 + self.component.isolation.cond_cond[num]  # from left to top
                                y = bot_bound + self.component.windings[num].conductor_radius

                        if self.component.virtual_winding_windows[num].scheme == WindingScheme.Square_Full_Width:
                            y = bot_bound + self.component.windings[num].conductor_radius
                            x = left_bound + self.component.windings[num].conductor_radius
                            i = 0
                            # Case n_conductors higher that "allowed" is missing
                            while y < top_bound - self.component.windings[num].conductor_radius \
                                    and i < self.component.windings[num].turns[n_win]:
                                while x < right_bound - self.component.windings[num].conductor_radius \
                                        and i < self.component.windings[num].turns[n_win]:
                                    self.p_conductor[num].append([x, y, 0, self.component.mesh.c_center_conductor[num]])
                                    self.p_conductor[num].append(
                                        [x - self.component.windings[num].conductor_radius, y, 0,
                                         self.component.mesh.c_conductor[num]])
                                    self.p_conductor[num].append(
                                        [x, y + self.component.windings[num].conductor_radius, 0,
                                         self.component.mesh.c_conductor[num]])
                                    self.p_conductor[num].append(
                                        [x + self.component.windings[num].conductor_radius, y, 0,
                                         self.component.mesh.c_conductor[num]])
                                    self.p_conductor[num].append(
                                        [x, y - self.component.windings[num].conductor_radius, 0,
                                         self.component.mesh.c_conductor[num]])
                                    i += 1
                                    x += self.component.windings[num].conductor_radius * 2 + \
                                         self.component.isolation.cond_cond[
                                             num]  # from left to top
                                y += self.component.windings[num].conductor_radius * 2 + \
                                     self.component.isolation.cond_cond[
                                         num]  # one step from left to right
                                x = left_bound + self.component.windings[num].conductor_radius  # always the same

                        if self.component.virtual_winding_windows[num].scheme == WindingScheme.Hexagonal:
                            y = bot_bound + self.component.windings[num].conductor_radius
                            x = left_bound + self.component.windings[num].conductor_radius
                            i = 0
                            base_line = True
                            # Case n_conductors higher that "allowed" is missing
                            while x < right_bound - self.component.windings[num].conductor_radius \
                                    and i < self.component.windings[num].turns[n_win]:
                                while y < top_bound - self.component.windings[num].conductor_radius and \
                                        i < self.component.windings[num].turns[n_win]:
                                    self.p_conductor[num].append([x, y, 0, self.component.mesh.c_center_conductor[num]])
                                    self.p_conductor[num].append(
                                        [x - self.component.windings[num].conductor_radius, y, 0,
                                         self.component.mesh.c_conductor[num]])
                                    self.p_conductor[num].append(
                                        [x, y + self.component.windings[num].conductor_radius, 0,
                                         self.component.mesh.c_conductor[num]])
                                    self.p_conductor[num].append(
                                        [x + self.component.windings[num].conductor_radius, y, 0,
                                         self.component.mesh.c_conductor[num]])
                                    self.p_conductor[num].append(
                                        [x, y - self.component.windings[num].conductor_radius, 0,
                                         self.component.mesh.c_conductor[num]])
                                    i += 1
                                    y += self.component.windings[num].conductor_radius * 2 + \
                                         self.component.isolation.cond_cond[
                                             num]  # from bottom to top
                                x += 2 * np.cos(np.pi / 6) * (
                                        self.component.windings[num].conductor_radius +
                                        self.component.isolation.cond_cond[num] / 2)
                                # * np.sqrt(2 / 3 * np.pi / np.sqrt(3))  # one step from left to right
                                # depending on what line, hexa scheme starts shifted
                                # reset y to "new" bottom
                                base_line = (not base_line)
                                if base_line:
                                    y = bot_bound + self.component.windings[num].conductor_radius
                                else:
                                    y = bot_bound + 2 * self.component.windings[num].conductor_radius + \
                                        self.component.isolation.cond_cond[
                                            num] / 2

            # Checking the Conductors
            for num in range(0, self.component.n_windings):
                # Convert to numpy
                # Check if all Conductors could be resolved
                self.p_conductor[num] = np.asarray(self.p_conductor[num])

                # TODO:CHECKS for rect. conductors
                """ CHECK: rectangle conductors with 4 points
                if self.component.windings[num].conductor_type == "full" or 
                        self.component.windings[num].conductor_type == "stacked" or \
                        self.component.windings[num].conductor_type == "foil":
                    if int(self.p_conductor[num].shape[0]/4) < self.component.windings[num].turns:
                        warnings.warn("Too many turns that do not fit in the winding window.")
                        # self.component.windings[num].turns = int(self.p_conductor[num].shape[0]/4)
                        self.component.valid = None
                """

                # CHECK: round conductors with 5 points
                if self.component.windings[num].conductor_type in [ConductorType.Solid, ConductorType.Litz]:
                    if int(self.p_conductor[num].shape[0] / 5) < sum(self.component.windings[num].turns):
                        # Warning: warnings.warn("Too many turns that do not fit in the winding window.")
                        # Correct: self.component.windings[num].turns = int(self.p_conductor[num].shape[0]/5)
                        # TODO: break, but remove warning. valid bit should be set to False
                        #  Code must go to the next parameter-iteration step for geometric sweep
                        self.component.valid = False
                        warnings.warn("Too many turns that do not fit in the winding window.")

            # Region for Boundary Condition
            self.p_region_bound[0][:] = [-self.r_outer * self.component.mesh.padding,
                                         -(self.component.core.window_h / 2 + self.component.core.core_w / 4)
                                         * self.component.mesh.padding,
                                         0,
                                         self.component.mesh.c_core * self.component.mesh.padding]
            self.p_region_bound[1][:] = [self.r_outer * self.component.mesh.padding,
                                         -(self.component.core.window_h / 2 + self.component.core.core_w / 4)
                                         * self.component.mesh.padding,
                                         0,
                                         self.component.mesh.c_core * self.component.mesh.padding]
            self.p_region_bound[2][:] = [-self.r_outer * self.component.mesh.padding,
                                         (self.component.core.window_h / 2 + self.component.core.core_w / 4)
                                         * self.component.mesh.padding,
                                         0,
                                         self.component.mesh.c_core * self.component.mesh.padding]
            self.p_region_bound[3][:] = [self.r_outer * self.component.mesh.padding,
                                         (self.component.core.window_h / 2 + self.component.core.core_w / 4)
                                         * self.component.mesh.padding,
                                         0,
                                         self.component.mesh.c_core * self.component.mesh.padding]

        def draw_isolations(self, isolation_deltas):
            """
            DISCLAIMER
            Because the distance from the core to the winding is set by
            iso.core_cond, a delta, which is used in order for no overlapping lines will cause
            the "real" isolation to be slightly smaller than set by the user.
            """

            window_h = self.component.core.window_h
            iso = self.component.isolation
            mesh = self.component.mesh

            # Since there are many cases in which alternating conductors would lead to slightly different
            # mesh densities a simplification is made: Just use the lowest mesh density to be safe all the time.
            mesh_density_to_winding = min(mesh.c_conductor)

            mesh_density_to_core = mesh.c_window
             

            # Using the delta the lines and points from the isolation and the core/windings are not overlapping
            # which makes creating the mesh more simpler
            # Isolation between winding and core
            iso_core_delta_left = isolation_deltas["core_left"]
            iso_core_delta_top = isolation_deltas["core_top"]
            iso_core_delta_right = isolation_deltas["core_right"]
            iso_core_delta_bot = isolation_deltas["core_bot"]
            iso_iso_delta = isolation_deltas["iso_iso"]
            iso_winding_delta_left = isolation_deltas["winding_left"]
            iso_winding_delta_top = isolation_deltas["winding_top"]
            iso_winding_delta_right = isolation_deltas["winding_right"]
            iso_winding_delta_bot = isolation_deltas["winding_bot"]

            self.p_iso_core = [] # Order: Left, Top, Right, Bot
            self.p_iso_pri_sec = []

            if self.component.component_type == ComponentType.IntegratedTransformer:
                # TODO implement for integrated_transformers
                warnings.warn("Isolations are not set because they are not implemented for integrated transformers.")
            else:
                # Useful points for isolation creation
                left_x = self.component.core.core_w / 2 + iso_core_delta_left
                top_y = window_h / 2 - iso_core_delta_top
                right_x = self.r_inner - iso_core_delta_right
                bot_y = -window_h / 2 + iso_core_delta_bot

                # Useful lengths
                left_iso_width = iso.core_cond[2] - iso_core_delta_left - iso_winding_delta_left
                top_iso_height = iso.core_cond[0] - iso_core_delta_top - iso_winding_delta_top
                right_iso_width = iso.core_cond[3] - iso_core_delta_right - iso_winding_delta_right
                bot_iso_height = iso.core_cond[1] - iso_core_delta_bot - iso_winding_delta_bot

                # Core to Pri isolation
                iso_core_left = [
                    [
                        left_x,
                        top_y - top_iso_height - iso_iso_delta,
                        0,
                        mesh_density_to_core
                    ],
                    [
                        left_x + left_iso_width,
                        top_y - top_iso_height - iso_iso_delta,
                        0,
                        mesh_density_to_winding
                    ],
                    [
                        left_x + left_iso_width,
                        bot_y + bot_iso_height + iso_iso_delta,
                        0,
                        mesh_density_to_winding
                    ],
                    [
                        left_x,
                        bot_y + bot_iso_height + iso_iso_delta,
                        0,
                        mesh_density_to_core
                    ]
                ]
                iso_core_top = [
                    [
                        left_x,
                        top_y,
                        0,
                        mesh_density_to_core
                    ],
                    [
                        right_x,
                        top_y,
                        0,
                        mesh_density_to_core
                    ],
                    [
                        right_x,
                        top_y - top_iso_height,
                        0,
                        mesh_density_to_winding
                    ],
                    [
                        left_x,
                        top_y - top_iso_height,
                        0,
                        mesh_density_to_winding
                    ]
                ]
                iso_core_right = [
                    [
                        right_x - right_iso_width,
                        top_y - top_iso_height - iso_iso_delta,
                        0,
                        mesh_density_to_winding
                    ],
                    [
                        right_x,
                        top_y - top_iso_height - iso_iso_delta,
                        0,
                        mesh_density_to_core
                    ],
                    [
                        right_x,
                        bot_y + bot_iso_height + iso_iso_delta,
                        0,
                        mesh_density_to_core
                    ],
                    [
                        right_x - right_iso_width,
                        bot_y + bot_iso_height + iso_iso_delta,
                        0,
                        mesh_density_to_winding
                    ]
                ]
                iso_core_bot = [
                    [
                        left_x,
                        bot_y + bot_iso_height,
                        0,
                        mesh_density_to_winding
                    ],
                    [
                        right_x,
                        bot_y + bot_iso_height,
                        0,
                        mesh_density_to_winding
                    ],
                    [
                        right_x,
                        bot_y,
                        0,
                        mesh_density_to_core
                    ],
                    [
                        left_x,
                        bot_y,
                        0,
                        mesh_density_to_core
                    ]
                ]

                self.p_iso_core = [iso_core_left, iso_core_top, iso_core_right, iso_core_bot]

                # Isolation between virtual winding windows
                if self.component.vw_type == VirtualWindingType.FullWindow:
                    # Only one vww -> The winding can be interleaved 
                    # -> still an isolation between pri and sec necessary
                    if len(self.component.virtual_winding_windows) == 1 and self.component.virtual_winding_windows[0].winding == WindingType.Interleaved:
                        # vertical isolations needed between the layers
                        # bifilar and vertical do not exist yet
                        vww = self.component.virtual_winding_windows[0]
                        winding_0 = self.component.windings[0]
                        winding_1 = self.component.windings[1]
                        current_x = vww.left_bound + 2 * winding_0.conductor_radius

                        mesh_density_left_side = 0
                        mesh_density_right_side = 0
                        if winding_0.turns[0] > winding_1.turns[0]:
                            mesh_density_left_side = mesh.c_conductor[0]
                            mesh_density_right_side = mesh.c_conductor[1]
                        else:
                            mesh_density_left_side = mesh.c_conductor[1]
                            mesh_density_right_side = mesh.c_conductor[0]

                        if vww.scheme == WindingScheme.Horizontal:
                            self.p_iso_pri_sec = []
                            for index in range(self.top_window_iso_counter-1):
                                self.p_iso_pri_sec.append([
                                    [
                                        current_x + iso_winding_delta_left,
                                        top_y - top_iso_height - iso_iso_delta,
                                        0,
                                        mesh_density_left_side
                                    ],
                                    [
                                        current_x + iso.cond_cond[2] - iso_winding_delta_right,
                                        top_y - top_iso_height - iso_iso_delta,
                                        0,
                                        mesh_density_right_side
                                    ],
                                    [
                                        current_x + iso.cond_cond[2] - iso_winding_delta_right,
                                        bot_y + bot_iso_height + iso_iso_delta,
                                        0,
                                        mesh_density_right_side
                                    ],
                                    [
                                        current_x + iso_winding_delta_left,
                                        bot_y + bot_iso_height + iso_iso_delta,
                                        0,
                                        mesh_density_left_side
                                    ]
                                ])
                                # The sec winding can start first when it has a higher number of turns
                                # col_cond
                                if index % 2 == self.col_cond_start:
                                    current_x += iso.cond_cond[2] + 2 * winding_1.conductor_radius
                                else:
                                    current_x += iso.cond_cond[2] + 2 * winding_0.conductor_radius
                        elif vww.scheme == WindingScheme.Vertical:
                            raise Exception("Vertical scheme not implemented yet!")
                        elif vww.scheme == WindingScheme.Bifilar:
                            raise Exception("Bifilar scheme not implemented yet!")
                        else:
                            raise Exception(f"The winding scheme {vww.scheme} is unknown.")
                elif self.component.vw_type == VirtualWindingType.Split2:
                    # Two vwws -> a horizontal isolation is needed
                    vww_bot = self.component.virtual_winding_windows[0]
                    vww_top = self.component.virtual_winding_windows[1]
                    self.p_iso_pri_sec.append([
                        [
                            vww_top.left_bound,
                            vww_top.bot_bound - iso_winding_delta_top,
                            0,
                            mesh_density_to_winding
                        ],
                        [
                            vww_top.right_bound,
                            vww_top.bot_bound - iso_winding_delta_top,
                            0,
                            mesh_density_to_winding
                        ],
                        [
                            vww_top.right_bound,
                            vww_bot.top_bound + iso_winding_delta_bot,
                            0,
                            mesh_density_to_winding
                        ],
                        [
                            vww_top.left_bound,
                            vww_bot.top_bound + iso_winding_delta_bot,
                            0,
                            mesh_density_to_winding
                        ]
                    ])
                else:
                    warnings.warn(f"Isolations are not implemented for components with type {self.component.vw_type}")

        def update(self, isolation_deltas = None):

            # Preallocate the arrays, in which the geometries' point coordinates will be stored
            self.p_outer = np.zeros((4, 4))
            self.p_region_bound = np.zeros((4, 4))
            self.p_window = np.zeros((4 * self.component.n_windows, 4))
            self.p_air_gaps = np.zeros((4 * self.component.air_gaps.number, 4)) 

            # Fitting the outer radius to ensure surface area
            self.r_inner = self.component.core.window_w + self.component.core.core_w / 2
            if self.component.correct_outer_leg:
                A_out = 200 * 10 ** -6
                self.r_outer = np.sqrt(A_out / np.pi + self.r_inner ** 2)  # Hardcode for PQ 40/40
            else:
                self.r_outer = np.sqrt((self.component.core.core_w / 2) ** 2 + self.r_inner ** 2)

            # np.sqrt(window_w**2 + window_w * core_w + core_w**2/2)
            #
            self.draw_outer()

            self.draw_window()

            self.draw_air_gaps()

            self.draw_virtual_winding_windows()

            self.draw_conductors()

            if isolation_deltas is None:
                isolation_deltas = {
                    "core_left": 0.00001,
                    "core_top": 0.00001,
                    "core_bot": 0.00001,
                    "core_right": 0.00001,
                    "iso_iso" : 0.00001,
                    "winding_left": 0.00001,
                    "winding_top": 0.00001,
                    "winding_right": 0.00001,
                    "winding_bot": 0.00001
                }

            self.draw_isolations(isolation_deltas)

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

            self.component.create_folders(self.component.reluctance_model_folder_path)

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

                    # print(f"{b_peak=}")
                    # print(f"{A_part=}")

                    if b_peak > self.b_max:
                        air_gap_lengths[air_gap_name] = "saturated"
                        # print("saturated")
                        # break the outer for loop
                        break

                    # Go through the distributed air gaps (of each reluctance)
                    for n_distributed in range(0, self.n_ag_per_rel[n_reluctance]):

                        if self.air_gap_types[n_reluctance][n_distributed] == "round-inf":
                            # print(self.component.core.window_h,
                            #       self.component.stray_path.midpoint,
                            #       self.component.stray_path.width)

                            def r_sct(length):
                                return r_round_inf(l=length,
                                                   r=self.component.core.core_w / 2,
                                                   sigma=sigma(l=length,
                                                               w=self.component.core.core_w / 2,
                                                               R_equivalent=r_basis(l=length,
                                                                                    w=self.component.core.core_w,
                                                                                    h=h_basis))
                                                   ) - R_0

                            if self.visualize_airgaps:
                                def r_sct_ideal(length):
                                    return r_round_inf(l=length,
                                                       r=self.component.core.core_w / 2,
                                                       sigma=1
                                                       ) - R_0

                        if self.air_gap_types[n_reluctance][n_distributed] == "round-round":
                            # def r_sct(length):
                            #    return r_round_round(length)
                            pass

                        if self.air_gap_types[n_reluctance][n_distributed] == "cyl-cyl":
                            def r_sct(length):
                                return r_cyl_cyl(l=length,
                                                 sigma=sigma(l=length,
                                                             w=self.component.stray_path.width,
                                                             R_equivalent=r_basis(l=length,
                                                                                  w=self.component.stray_path.width,
                                                                                  h=self.component.core.window_w
                                                                                    - length) / 2),
                                                 w=self.component.stray_path.width,
                                                 r_o=self.component.core.window_w + self.component.core.core_w / 2
                                                 ) - R_0

                            def r_sct_real(length):
                                return r_cyl_cyl_real(l=length,
                                                      sigma=sigma(l=length,
                                                                  w=self.component.stray_path.width,
                                                                  R_equivalent=r_basis(l=length,
                                                                                       w=self.component.stray_path.width,
                                                                                       h=self.component.core.window_w
                                                                                         - length) / 2),
                                                      w=self.component.stray_path.width,
                                                      r_o=self.component.core.window_w + self.component.core.core_w / 2,
                                                      h_real_core=self.real_core_width
                                                      ) - R_0

                        # print(f"\n  {air_gap_name}")
                        # print(f"n_reluctance {n_reluctance}")
                        # print(f"self.component.stray_path.width {self.component.stray_path.width}")
                        # print(f"max_length[n_reluctance] {self.max_length[n_reluctance]}")
                        # print(f"R_0 {R_0}")
                        # print(f"r_sct(a) {r_sct(1e-6)}")
                        # print(f"r_sct(b) {r_sct(self.max_length[n_reluctance])}")

                        # Check for different signs (zero crossing)
                        if r_sct(1e-6) * r_sct(self.max_length[n_reluctance]) > 0:
                            air_gap_lengths[air_gap_name] = "out of bounds"
                            # print("out of bounds")

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

            # print(f"{air_gap_lengths=}")
            return air_gap_lengths, b_peaks

        def get_core_loss(self):
            """
            Calculates the hysteresis loss corresponding to complex core parameters.
            :return:
            """
            # print(f"Calculate Analytical Core Losses:")

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
            # print(p_top_1st, p_bot_1st, p_stray_1st)
            self.p_hyst_nom_1st = sum([p_top_1st, p_bot_1st, p_stray_1st])

            p_top, p_bot, p_stray = self.hysteresis_loss(Phi_top_nom_peak, Phi_bot_nom_peak, Phi_stray_nom_peak)
            # print(p_top, p_bot, p_stray)
            self.p_hyst_nom = sum([p_top, p_bot, p_stray])

            self.visualize_phi_core_loss(Phi_top_nom, Phi_bot_nom, Phi_stray_nom,
                                         phase_top, phase_bot, phase_stray)

            # print(f"Analytical Core Losses = \n\n")

        def visualize_phi_core_loss(self, Phi_top_init, Phi_bot_init, Phi_stray_init,
                                    phase_top, phase_bot, phase_stray):
            """
            Visualization of the fluxes used for the core loss calculation.
            :return:
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
            """

            :return:
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
            :return:
            """
            length_corner = self.A_core / self.component.core.core_w / np.pi
            ri = self.component.core.core_w / 2
            ro = ri + self.component.core.window_w

            # Top Part
            b_top = Phi_top / self.A_core
            Vol_top = self.A_core * (100 - self.component.stray_path.midpoint) / 100 * self.component.core.window_h + \
                      self.A_core * ((100 - self.component.stray_path.midpoint) / 100 * self.component.core.window_h - \
                                     self.air_gap_lengths["R_top"])

            p_top = 0.5 * self.component.mu0 * f_N95_mu_imag(self.f_1st, b_top) * Vol_top * 2 * np.pi * self.f_1st * \
                    (b_top / self.component.core.mu_rel / self.component.mu0) ** 2 + \
                    self.p_loss_cyl(Phi_top, length_corner, ri, ro)[0]

            # Bot Part
            b_bot = Phi_bot / self.A_core
            Vol_bot = self.A_core * self.component.stray_path.midpoint / 100 * self.component.core.window_h + \
                      self.A_core * (self.component.stray_path.midpoint / 100 * self.component.core.window_h + \
                                     self.air_gap_lengths["R_bot"])
            p_bot = 0.5 * self.component.mu0 * f_N95_mu_imag(self.f_1st, b_bot) * Vol_bot * 2 * np.pi * self.f_1st * \
                    (b_bot / self.component.core.mu_rel / self.component.mu0) ** 2 + \
                    self.p_loss_cyl(Phi_bot, length_corner, ri, ro)[0]

            # Stray Path
            p_stray = self.p_loss_cyl(Phi_stray, self.component.stray_path.width, ri,
                                      ro - self.air_gap_lengths["R_stray"])[0]

            return p_top, p_bot, p_stray

        def p_loss_cyl(self, Phi, w, ri, ro):
            """

            :return:
            """

            def b(r, Phi, w):
                return Phi / (2 * np.pi * r * w)

            def p_loss_density(r, Phi, w):
                return 2 * np.pi * r * w * \
                       np.pi * self.f_1st * \
                       self.component.mu0 * f_N95_mu_imag(self.f_1st, b(r, Phi, w)) * \
                       (b(r, Phi, w) / self.component.core.mu_rel / self.component.mu0) ** 2

            return quad(p_loss_density, ri, ro, args=(Phi, w), epsabs=1e-4)

        @staticmethod
        def max_phi_from_time_phi(Phi_top, Phi_bot, Phi_stray):
            """

            :param Phi_stray:
            :param Phi_bot:
            :param Phi_top:
            :return:
            """
            Phi_top_peak = max([abs(ele) for ele in Phi_top])
            Phi_bot_peak = max([abs(ele) for ele in Phi_bot])
            Phi_stray_peak = max([abs(ele) for ele in Phi_stray])

            # print(f"{Phi_top_peak=}\n"
            #       f"{Phi_bot_peak=}\n"
            #       f"{Phi_stray_peak=}\n")

            return Phi_top_peak, Phi_bot_peak, Phi_stray_peak

        def phases_from_time_phi(self, Phi_top, Phi_bot, Phi_stray):
            """
            Returns the phases of the peaks
            :param Phi_stray:
            :param Phi_bot:
            :param Phi_top:
            :return:
            """
            # print(np.array(Phi_top))
            # print(np.array(Phi_top).argmax(axis=0))
            # print(self.nom_current[0][0])

            phase_top = self.nom_current[0][0][np.array(Phi_top).argmax(axis=0)]
            phase_bot = self.nom_current[0][0][np.array(Phi_bot).argmax(axis=0)]
            phase_stray = self.nom_current[0][0][np.array(Phi_stray).argmax(axis=0)]

            # print(np.array(Phi_top).argmax(axis=0))
            # print(np.array(Phi_bot).argmax(axis=0))
            # print(np.array(Phi_stray).argmax(axis=0))

            return phase_top, phase_bot, phase_stray

        def phi_from_time_currents(self, time_current_1, time_current_2, visualize=False):
            """

            :param visualize:
            :param time_current_2:
            :param time_current_1:
            :return:
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
            """
            Method defines instance variables.
            :return:
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
                # print(self.max_phi)

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
                # print(N)
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
                    R_matrix = calculate_reluctances(N=self.N, L=self.L_goal)

                    # R_goal = [R_top, R_bot, R_stray]
                    R_goal = [R_matrix[0, 0] + R_matrix[0, 1], R_matrix[1, 1] + R_matrix[0, 1], -R_matrix[0, 1]]

                    # Stray path specific parameters
                    # print(R_matrix)
                    # print(R_goal)
                    # print(self.L_goal)
                    self.stray_path_parametrization_two_d_axi()

                    # Check for negative Reluctances
                    if all(R >= 0 for R in R_goal) and self.singularity == False:

                        # Calculate the air gap lengths with the help of SCT
                        # air_gap_lengths is a dictionary with air gap names and the associated length
                        self.air_gap_lengths, self.b_peaks = self.calculate_air_gap_lengths_with_sct(reluctances=R_goal)

                        # print(air_gap_lengths.values())

                        # Check for invalid data
                        if self.air_gap_lengths.values() in ['saturated', 'out of bounds']:

                            results = None
                            # print("Invalid Data\n\n")

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

                # print(f"{results=}")
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
            print(f"\n"
                  f"--- ---\n"
                  f"Perform reluctance calculations\n\n"
                  # f"Goal inductance values                   : {L_goal}\n\n"
                  f"Number of initial reluctance parameters  : {len(parameters_init)}\n")

            # Update the core to use its internal core parameter calculation functionality
            # Set attributes of core with given keywords
            for index, parameters in enumerate(parameters_init):
                for key, value in parameters.items():
                    # print(key, value)
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

            print(f"Number of singularities: {self.n_singularities}\n")

            # Save Results including invalid parameters
            # print(f"{parameter_results=}")

            np.save(self.component.path + '/Reluctance_Model/parameter_results.npy', parameter_results)

            # Filter all entries, that are None
            # Elements of list reluctance_parameters are either a dictionary or None
            valid_parameter_results = [x for x in parameter_results if x is not None]

            # Save resulting valid parameters
            np.save(self.component.path + '/Reluctance_Model/valid_parameter_results.npy', valid_parameter_results)

            print(f"Number of valid parameters: {len(valid_parameter_results)}\n\n"
                  f"Ready with reluctance calculations\n"
                  f"--- ---\n")

            return valid_parameter_results

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    # Meshing
    class Mesh:
        """

        """

        def __init__(self, component):
            self.component = component

            # Initialize gmsh once
            if not gmsh.isInitialized():
                gmsh.initialize()

            # Characteristic lengths [for mesh sizes]
            self.global_accuracy = 0.5  # Parameter for mesh-accuracy
            self.padding = 1.5  # ... > 1

            self.skin_mesh_factor = None
            self.c_core = None
            self.c_window = None
            self.c_conductor = [None] * self.component.n_windings
            self.c_center_conductor = [None] * self.component.n_windings  # used for the mesh accuracy in the conductors
            # self.c_air_gap = None  # TODO: make usable again ?

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Pre-Definitions

        def set_empty_lists(self):
            # Points
            self.p_core = []
            self.p_island = []
            self.p_cond = [[], []]
            self.p_region = []
            self.p_iso_core = []
            self.p_iso_pri_sec = []
            # Curves
            self.l_bound_core = []
            self.l_bound_air = []
            self.l_core_air = []
            self.l_cond = [[], []]
            self.l_region = []
            self.l_air_gaps_air = []
            self.l_iso_core = []
            self.l_iso_pri_sec = []
            # Curve Loops
            self.curve_loop_cond = [[], []]
            self.curve_loop_island = []
            self.curve_loop_air = []
            self.curve_loop_air_gaps = []
            self.curve_loop_iso_core = []
            self.curve_loop_iso_pri_sec = []
            # curve_loop_outer_air = []
            # curve_loop_bound = []
            # Plane Surfaces
            self.plane_surface_core = []
            self.plane_surface_cond = [[], []]
            self.plane_surface_air = []
            self.plane_surface_outer_air = []
            self.plane_surface_air_gaps = []
            self.plane_surface_iso_core = []
            self.plane_surface_iso_pri_sec = []

        def generate_hybrid_mesh(self, color_scheme: Dict = colors_femmt_default, colors_geometry: Dict = colors_geometry_femmt_default,
                                 visualize_before: bool = False,
                                 save_png: bool = True, refine=0, alternative_error=0):
            """
            - interaction with gmsh
            - mesh generation
                - Skin depth based forward meshing
                - adaptive refinement [future TODO]
                    - with the help of mesh-size-fields/background meshes
                    - with an appropriate local error metric

            :return:

            """
            # Initialization
            self.set_empty_lists()

            print("Hybrid Mesh Generation in Gmsh")

            """
            if refine == 1:
                # TODO: Adaptive Meshing
                # Choose applied Error Function
                if alternative_error == 1:
                    local_error = self.alternative_local_error()  # something like current density
                else:
                    local_error = self.local_error()  # Here a "real" numeric error should be applied

                self.create_background_mesh(local_error)
            """

            # gmsh.option.setNumber("General.Terminal", 1)
            gmsh.model.add(os.path.join(self.component.e_m_mesh_file, "geometry"))


            # ------------------------------------------ Geometry -------------------------------------------
            # Core generation
            # --------------------------------------- Points --------------------------------------------
            if self.component.dimensionality == "2D":
                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                # Geometry definitions: points -> lines -> curve loops -> surfaces
                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                # Main Core

                # Points
                # (index refers to sketch)

                # First point (left point of lowest air gap)
                if self.component.air_gaps.number > 0:
                    self.p_core.append(gmsh.model.geo.addPoint(0,
                                                               self.component.two_d_axi.p_air_gaps[0][1],
                                                               0,
                                                               self.c_core))
                if self.component.air_gaps.number == 0:
                    self.p_core.append(None)  # dummy filled for no air gap special case

                # Go down and counter-clockwise
                # Four points around the core
                if self.component.correct_outer_leg:
                    correction_of_outer_points = 0.00
                else:
                    correction_of_outer_points = 0

                self.p_core.append(gmsh.model.geo.addPoint(0,
                                                           self.component.two_d_axi.p_outer[1][1],
                                                           self.component.two_d_axi.p_outer[1][2],
                                                           self.component.two_d_axi.p_outer[1][3]))

                self.p_core.append(gmsh.model.geo.addPoint(self.component.two_d_axi.p_outer[1][0],
                                                           self.component.two_d_axi.p_outer[1][1] + correction_of_outer_points,
                                                           self.component.two_d_axi.p_outer[1][2],
                                                           self.component.two_d_axi.p_outer[1][3]))

                self.p_core.append(gmsh.model.geo.addPoint(self.component.two_d_axi.p_outer[3][0],
                                                           self.component.two_d_axi.p_outer[3][1] - correction_of_outer_points,
                                                           self.component.two_d_axi.p_outer[3][2],
                                                           self.component.two_d_axi.p_outer[3][3]))

                self.p_core.append(gmsh.model.geo.addPoint(0,
                                                           self.component.two_d_axi.p_outer[3][1],
                                                           self.component.two_d_axi.p_outer[3][2],
                                                           self.component.two_d_axi.p_outer[3][3]))

                # Two points of highest air gap
                if self.component.air_gaps.number > 0:
                    self.p_core.append(gmsh.model.geo.addPoint(0,
                                                               self.component.two_d_axi.p_air_gaps[-2][1],
                                                               self.component.two_d_axi.p_air_gaps[-2][2],
                                                               self.c_core))

                    self.p_core.append(gmsh.model.geo.addPoint(self.component.two_d_axi.p_air_gaps[-1][0],
                                                               self.component.two_d_axi.p_air_gaps[-1][1],
                                                               self.component.two_d_axi.p_air_gaps[-1][2],
                                                               self.c_window))

                if self.component.air_gaps.number == 0:
                    self.p_core.append(None)  # dummy filled for no air gap special case
                    self.p_core.append(None)  # dummy filled for no air gap special case

                # Clockwise
                # Four points of the window
                self.p_core.append(gmsh.model.geo.addPoint(self.component.two_d_axi.p_window[6][0],
                                                           self.component.two_d_axi.p_window[6][1],
                                                           self.component.two_d_axi.p_window[6][2],
                                                           self.component.two_d_axi.p_window[6][3]))

                self.p_core.append(gmsh.model.geo.addPoint(self.component.two_d_axi.p_window[7][0],
                                                           self.component.two_d_axi.p_window[7][1],
                                                           self.component.two_d_axi.p_window[7][2],
                                                           self.component.two_d_axi.p_window[7][3]))

                self.p_core.append(gmsh.model.geo.addPoint(self.component.two_d_axi.p_window[5][0],
                                                           self.component.two_d_axi.p_window[5][1],
                                                           self.component.two_d_axi.p_window[5][2],
                                                           self.component.two_d_axi.p_window[5][3]))

                self.p_core.append(gmsh.model.geo.addPoint(self.component.two_d_axi.p_window[4][0],
                                                           self.component.two_d_axi.p_window[4][1],
                                                           self.component.two_d_axi.p_window[4][2],
                                                           self.component.two_d_axi.p_window[4][3]))

                # Last point of lowest air gap
                if self.component.air_gaps.number > 0:
                    self.p_core.append(gmsh.model.geo.addPoint(self.component.two_d_axi.p_air_gaps[1][0],
                                                               self.component.two_d_axi.p_air_gaps[1][1],
                                                               self.component.two_d_axi.p_air_gaps[1][2],
                                                               self.c_window))

                if self.component.air_gaps.number == 0:
                    self.p_core.append(None)  # dummy filled for no air gap special case

                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                # Curves
                # (index refers to sketch)
                # To be added: Case Air Gaps directly on outer leg

                # Curves: Boundary - Core
                if self.component.air_gaps.number > 0:
                    self.l_bound_core.append(gmsh.model.geo.addLine(self.p_core[0],
                                                                    self.p_core[1]))
                self.l_bound_core.append(gmsh.model.geo.addLine(self.p_core[1],
                                                                self.p_core[2]))
                self.l_bound_core.append(gmsh.model.geo.addLine(self.p_core[2],
                                                                self.p_core[3]))
                self.l_bound_core.append(gmsh.model.geo.addLine(self.p_core[3],
                                                                self.p_core[4]))
                if self.component.air_gaps.number > 0:
                    self.l_bound_core.append(gmsh.model.geo.addLine(self.p_core[4],
                                                                    self.p_core[5]))
                if self.component.air_gaps.number == 0:
                    self.l_bound_core.append(gmsh.model.geo.addLine(self.p_core[4],
                                                                    self.p_core[1]))
                # Curves: Core - Air
                if self.component.air_gaps.number > 0:
                    self.l_core_air.append(gmsh.model.geo.addLine(self.p_core[5],
                                                                  self.p_core[6]))
                    self.l_core_air.append(gmsh.model.geo.addLine(self.p_core[6],
                                                                  self.p_core[7]))
                self.l_core_air.append(gmsh.model.geo.addLine(self.p_core[7],
                                                              self.p_core[8]))
                self.l_core_air.append(gmsh.model.geo.addLine(self.p_core[8],
                                                              self.p_core[9]))
                self.l_core_air.append(gmsh.model.geo.addLine(self.p_core[9],
                                                              self.p_core[10]))

                if self.component.air_gaps.number > 0:
                    self.l_core_air.append(gmsh.model.geo.addLine(self.p_core[10],
                                                                  self.p_core[11]))
                    self.l_core_air.append(gmsh.model.geo.addLine(self.p_core[11],
                                                                  self.p_core[0]))
                if self.component.air_gaps.number == 0:
                    self.l_core_air.append(gmsh.model.geo.addLine(self.p_core[10],
                                                                  self.p_core[7]))

                # Plane: Main Core --> plane_surface_core[0]
                if self.component.air_gaps.number > 0:
                    curve_loop_core = gmsh.model.geo.addCurveLoop(self.l_bound_core + self.l_core_air)
                    self.plane_surface_core.append(gmsh.model.geo.addPlaneSurface([-curve_loop_core]))
                if self.component.air_gaps.number == 0:
                    curve_loop_bound_core = gmsh.model.geo.addCurveLoop(self.l_bound_core)
                    curve_loop_core_air = gmsh.model.geo.addCurveLoop(self.l_core_air)
                    self.plane_surface_core.append(gmsh.model.geo.addPlaneSurface([-curve_loop_bound_core, curve_loop_core_air]))

                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                # - Core parts between Air Gaps
                # Points of Core Islands (index refers to sketch)


                stray_path_air_gap_top_point = None
                stray_path_air_gap_bot_point = None
                l_core_air_air_gap = []

                if self.component.air_gaps.number > 1:
                    for point in self.component.two_d_axi.p_air_gaps[2:-2]:
                        self.p_island.append(gmsh.model.geo.addPoint(*point))

                    # Add two more points for closing of the air gap for a stray_path
                    if self.component.component_type == ComponentType.IntegratedTransformer:
                        stray_path_air_gap_top_point = gmsh.model.geo.addPoint(*self.component.two_d_axi.p_close_air_gaps[0])
                        stray_path_air_gap_bot_point = gmsh.model.geo.addPoint(*self.component.two_d_axi.p_close_air_gaps[1])

                    # Curves of Core Islands (index refers to sketch)
                    for i in range(0, int(len(self.p_island) / 4)):
                        if self.component.component_type == ComponentType.IntegratedTransformer and self.component.stray_path.start_index == i:
                            l_core_air_air_gap.append(gmsh.model.geo.addLine(self.p_island[4 * i + 0], stray_path_air_gap_bot_point))
                            self.l_core_air.append(gmsh.model.geo.addLine(stray_path_air_gap_bot_point, self.p_island[4 * i + 1]))
                            self.l_core_air.append(gmsh.model.geo.addLine(self.p_island[4 * i + 1], self.p_island[4 * i + 3]))
                            self.l_core_air.append(gmsh.model.geo.addLine(self.p_island[4 * i + 3], stray_path_air_gap_top_point))
                            l_core_air_air_gap.append(gmsh.model.geo.addLine(stray_path_air_gap_top_point, self.p_island[4 * i + 2]))
                            self.l_bound_core.append(gmsh.model.geo.addLine(self.p_island[4 * i + 2], self.p_island[4 * i + 0]))

                            self.curve_loop_island.append(gmsh.model.geo.addCurveLoop(
                                [l_core_air_air_gap[-2], self.l_core_air[-3], self.l_core_air[-2], l_core_air_air_gap[-1], self.l_core_air[-1], self.l_bound_core[-1]]))
                        else:
                            # Default
                            self.l_core_air.append(gmsh.model.geo.addLine(self.p_island[4 * i + 0], self.p_island[4 * i + 1]))
                            self.l_core_air.append(gmsh.model.geo.addLine(self.p_island[4 * i + 1], self.p_island[4 * i + 3]))
                            self.l_core_air.append(gmsh.model.geo.addLine(self.p_island[4 * i + 3], self.p_island[4 * i + 2]))
                            self.l_bound_core.append(gmsh.model.geo.addLine(self.p_island[4 * i + 2], self.p_island[4 * i + 0]))
                            
                            self.curve_loop_island.append(gmsh.model.geo.addCurveLoop(
                                [self.l_core_air[-3], self.l_core_air[-2], self.l_core_air[-1], self.l_bound_core[-1]]))
                        
                        # Iterative plane creation
                        self.plane_surface_core.append(gmsh.model.geo.addPlaneSurface([-self.curve_loop_island[-1]]))

                # Curves: Boundary - Air
                if self.component.air_gaps.number == 1:
                    self.l_bound_air.append(gmsh.model.geo.addLine(self.p_core[0], self.p_core[5]))
                else:
                    for i in range(0, int(len(self.p_island) / 4)):
                        if i == 0:  # First Line
                            self.l_bound_air.append(gmsh.model.geo.addLine(self.p_core[0], self.p_island[0]))
                        else:  # Middle Lines
                            self.l_bound_air.append(
                                gmsh.model.geo.addLine(self.p_island[4 * (i - 1) + 2], self.p_island[4 * i + 0]))
                        if i == int(len(self.p_island) / 4) - 1:  # Last Line
                            self.l_bound_air.append(gmsh.model.geo.addLine(self.p_island[-2], self.p_core[5]))
                
                # Curves: Close air gaps
                if self.component.air_gaps.number > 0:
                    for i in range(0, self.component.air_gaps.number):
                        bottom_point = self.p_core[11] if i == 0 else self.p_island[(i-1)*4+3]
                        top_point = self.p_core[6] if i == self.component.air_gaps.number-1 else self.p_island[i*4+1]

                        if self.component.component_type == ComponentType.IntegratedTransformer:
                            if self.component.stray_path.start_index == i:
                                # Stray path is above current air_gap
                                top_point = stray_path_air_gap_bot_point
                            elif self.component.stray_path.start_index + 1 == i: 
                                # Stray path is below current air_gap
                                bottom_point = stray_path_air_gap_top_point

                        self.l_air_gaps_air.append(gmsh.model.geo.addLine(bottom_point, top_point))

                for i in range(0, self.component.air_gaps.number):
                    left = self.l_bound_air[i]
                    right  = self.l_air_gaps_air[i]

                    if i == self.component.air_gaps.number-1:
                        top = self.l_core_air[0]
                    elif self.component.component_type == ComponentType.IntegratedTransformer and self.component.stray_path.start_index == i:
                        top = l_core_air_air_gap[0]
                    else:
                        top = self.l_core_air[7+3*i]
                        
                    if i == 0:
                        bottom = self.l_core_air[6] 
                    elif self.component.component_type == ComponentType.IntegratedTransformer and self.component.stray_path.start_index + 1 == i:
                        bottom = l_core_air_air_gap[1]
                    else:
                        bottom = self.l_core_air[6+3*i]

                    curve_loop = gmsh.model.geo.addCurveLoop([left, top, bottom, right], -1, True)
                    self.curve_loop_air_gaps.append(curve_loop)
                    self.plane_surface_air_gaps.append(gmsh.model.geo.addPlaneSurface([curve_loop]))
                
                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                # Conductors
                # Points of Conductors
                for num in range(0, self.component.n_windings):
                    for i in range(0, self.component.two_d_axi.p_conductor[num].shape[0]):
                        self.p_cond[num].append(
                            gmsh.model.geo.addPoint(
                                self.component.two_d_axi.p_conductor[num][i][0],  # x
                                self.component.two_d_axi.p_conductor[num][i][1],  # y
                                0,  # z
                                self.component.two_d_axi.p_conductor[num][i][3]))

                    # Curves of Conductors
                    if self.component.windings[num].conductor_type in [ConductorType.Litz, ConductorType.Solid]:
                        for i in range(0, int(len(self.p_cond[num]) / 5)):
                            self.l_cond[num].append(gmsh.model.geo.addCircleArc(
                                self.p_cond[num][5 * i + 1],
                                self.p_cond[num][5 * i + 0],
                                self.p_cond[num][5 * i + 2]))
                            self.l_cond[num].append(gmsh.model.geo.addCircleArc(
                                self.p_cond[num][5 * i + 2],
                                self.p_cond[num][5 * i + 0],
                                self.p_cond[num][5 * i + 3]))
                            self.l_cond[num].append(gmsh.model.geo.addCircleArc(
                                self.p_cond[num][5 * i + 3],
                                self.p_cond[num][5 * i + 0],
                                self.p_cond[num][5 * i + 4]))
                            self.l_cond[num].append(gmsh.model.geo.addCircleArc(
                                self.p_cond[num][5 * i + 4],
                                self.p_cond[num][5 * i + 0],
                                self.p_cond[num][5 * i + 1]))
                            # Iterative plane creation
                            self.curve_loop_cond[num].append(gmsh.model.geo.addCurveLoop(
                                [self.l_cond[num][i * 4 + 0],
                                 self.l_cond[num][i * 4 + 1],
                                 self.l_cond[num][i * 4 + 2],
                                 self.l_cond[num][i * 4 + 3]]))
                            self.plane_surface_cond[num].append(
                                gmsh.model.geo.addPlaneSurface([self.curve_loop_cond[num][i]]))
                    else:
                        # Rectangle conductor cut
                        for i in range(0, int(len(self.p_cond[num]) / 4)):
                            self.l_cond[num].append(gmsh.model.geo.addLine(self.p_cond[num][4 * i + 0],
                                                                           self.p_cond[num][4 * i + 2]))
                            self.l_cond[num].append(gmsh.model.geo.addLine(self.p_cond[num][4 * i + 2],
                                                                           self.p_cond[num][4 * i + 3]))
                            self.l_cond[num].append(gmsh.model.geo.addLine(self.p_cond[num][4 * i + 3],
                                                                           self.p_cond[num][4 * i + 1]))
                            self.l_cond[num].append(gmsh.model.geo.addLine(self.p_cond[num][4 * i + 1],
                                                                           self.p_cond[num][4 * i + 0]))
                            # Iterative plane creation
                            self.curve_loop_cond[num].append(gmsh.model.geo.addCurveLoop([self.l_cond[num][i * 4 + 0],
                                                                                          self.l_cond[num][i * 4 + 1],
                                                                                          self.l_cond[num][i * 4 + 2],
                                                                                          self.l_cond[num][i * 4 + 3]]))
                            self.plane_surface_cond[num].append(
                                gmsh.model.geo.addPlaneSurface([self.curve_loop_cond[num][i]]))


                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                # Isolations
                # Core to Pri
                if self.component.two_d_axi.p_iso_core: # Check if list is not empty
                    # Points
                    for iso in self.component.two_d_axi.p_iso_core:
                        p_iso = []
                        for i in iso:
                            p_iso.append(gmsh.model.geo.addPoint(i[0], i[1], i[2], i[3]))
                        self.p_iso_core.append(p_iso)
                    # Lines
                    self.l_iso_core = [[gmsh.model.geo.addLine(iso[i], iso[(i+1)%4]) for i in range(4)] for iso in self.p_iso_core]

                    # Curve loop and surface
                    self.curve_loop_iso_core = []
                    self.plane_surface_iso_core = []
                    for iso in self.l_iso_core:
                        cl = gmsh.model.geo.addCurveLoop(iso)
                        self.curve_loop_iso_core.append(cl)
                        self.plane_surface_iso_core.append(gmsh.model.geo.addPlaneSurface([cl]))

                # Pri to Sec
                if self.component.two_d_axi.p_iso_pri_sec: # Check if list is not empty
                    # Points
                    for iso in self.component.two_d_axi.p_iso_pri_sec:
                        p_iso = []
                        for i in iso:
                            p_iso.append(gmsh.model.geo.addPoint(i[0], i[1], i[2], i[3]))
                        self.p_iso_pri_sec.append(p_iso)

                    # Lines
                    # iso_pri_sec contains a list of multiple isolations
                    self.l_iso_pri_sec = [[gmsh.model.geo.addLine(iso[i], iso[(i+1)%4]) for i in range(4)] for iso in self.p_iso_pri_sec]
                        
                    # Curve loops and surfaces
                    self.curve_loop_iso_pri_sec = []
                    self.plane_surface_iso_pri_sec = []
                    for iso in self.l_iso_pri_sec:
                        cl = gmsh.model.geo.addCurveLoop(iso)
                        self.curve_loop_iso_pri_sec.append(cl)
                        self.plane_surface_iso_pri_sec.append(gmsh.model.geo.addPlaneSurface([cl]))

                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                # Air
                # Points are partwise double designated

                """ This is without the split of the air gaps
                l_air_tmp = self.l_core_air[:7]
                for i in range(0, len(self.l_bound_air)):
                    l_air_tmp.append(self.l_bound_air[i])
                    if i < len(self.l_bound_air) - 1:
                        l_air_tmp.append(self.l_core_air[7 + 3 * i])
                        l_air_tmp.append(self.l_core_air[7 + 3 * i + 1])
                        l_air_tmp.append(self.l_core_air[7 + 3 * i + 2])

                self.curve_loop_air.append(gmsh.model.geo.addCurveLoop(l_air_tmp))
                """

                # With closed air gaps
                l_air_tmp = []
                if self.component.air_gaps.number == 0:
                    l_air_tmp = self.l_core_air
                elif self.component.air_gaps.number > 0:
                    l_air_tmp = self.l_core_air[1:6] + self.l_air_gaps_air

                    for i in range(self.component.air_gaps.number - 1):
                        if self.component.component_type == ComponentType.IntegratedTransformer and i == self.component.stray_path.start_index:
                            l_air_tmp.append(self.l_core_air[7+3*i])
                            l_air_tmp.append(self.l_core_air[8+3*i])
                            l_air_tmp.append(self.l_core_air[9+3*i])
                        else:
                            l_air_tmp.append(self.l_core_air[8+3*i])
                            
                #for i in range(0, self.component.air_gaps.number):
                #    l_air_tmp.append(self.l_air_gaps_air[i])
                #    l_air_tmp.append(self.l_air_gaps_air[i+1])

                #self.curve_loop_air.append(gmsh.model.geo.addCurveLoop(l_air_tmp, -1, True))
                #for i in range(0, self.component.air_gaps.number):
                #    l_air_tmp.append(self.l_air_gaps_air[i])
                #    l_air_tmp.append(self.l_air_gaps_air[i+1])

                self.curve_loop_air.append(gmsh.model.geo.addCurveLoop(l_air_tmp, -1, True))

                # Need flatten list of all! conductors
                flatten_curve_loop_cond = [j for sub in self.curve_loop_cond for j in sub]

                # The first curve loop represents the outer bounds: self.curve_loop_air (should only contain one element)
                # The other curve loops represent holes in the surface -> For each conductor as well as each isolation
                self.plane_surface_air.append(
                    gmsh.model.geo.addPlaneSurface(self.curve_loop_air + flatten_curve_loop_cond + self.curve_loop_iso_core + self.curve_loop_iso_pri_sec))

                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                # Boundary
                if self.component.region is None:
                    self.l_bound_tmp = self.l_bound_core[:5]
                    for i in range(0, len(self.l_bound_air)):
                        self.l_bound_tmp.append(self.l_bound_air[-i - 1])
                        if i != len(self.l_bound_air) - 1:  # last run
                            self.l_bound_tmp.append(self.l_bound_core[-i - 1])

                else:
                    # Generate Lines of Region
                    # start top left and go clockwise
                    self.p_region.append(gmsh.model.geo.addPoint(0,
                                                                 self.component.two_d_axi.p_region_bound[2][1],
                                                                 self.component.two_d_axi.p_region_bound[2][2],
                                                                 self.component.two_d_axi.p_region_bound[2][3]))

                    self.p_region.append(gmsh.model.geo.addPoint(self.component.two_d_axi.p_region_bound[3][0],
                                                                 self.component.two_d_axi.p_region_bound[3][1],
                                                                 self.component.two_d_axi.p_region_bound[3][2],
                                                                 self.component.two_d_axi.p_region_bound[3][3]))

                    self.p_region.append(gmsh.model.geo.addPoint(self.component.two_d_axi.p_region_bound[1][0],
                                                                 self.component.two_d_axi.p_region_bound[1][1],
                                                                 self.component.two_d_axi.p_region_bound[1][2],
                                                                 self.component.two_d_axi.p_region_bound[1][3]))

                    self.p_region.append(gmsh.model.geo.addPoint(0,
                                                                 self.component.two_d_axi.p_region_bound[0][1],
                                                                 self.component.two_d_axi.p_region_bound[0][2],
                                                                 self.component.two_d_axi.p_region_bound[0][3]))

                    # Outer Region Lines
                    self.l_region.append(gmsh.model.geo.addLine(self.p_core[4],
                                                                self.p_region[0]))
                    self.l_region.append(gmsh.model.geo.addLine(self.p_region[0],
                                                                self.p_region[1]))
                    self.l_region.append(gmsh.model.geo.addLine(self.p_region[1],
                                                                self.p_region[2]))
                    self.l_region.append(gmsh.model.geo.addLine(self.p_region[2],
                                                                self.p_region[3]))
                    self.l_region.append(gmsh.model.geo.addLine(self.p_region[3],
                                                                self.p_core[1]))

                    # Boundary Line
                    self.l_bound_tmp = [self.l_bound_core[4]]

                    for i in range(0, len(self.l_region)):
                        self.l_bound_tmp.append(self.l_region[i])

                    self.l_bound_tmp.append(self.l_bound_core[0])

                    for i in range(0, len(self.l_bound_air)):
                        self.l_bound_tmp.append(self.l_bound_air[-i - 1])
                        if i != len(self.l_bound_air) - 1:  # last run
                            self.l_bound_tmp.append(self.l_bound_core[-i - 1])

                    # Outer Air Surface
                    curve_loop_outer_air = gmsh.model.geo.addCurveLoop(self.l_region + self.l_bound_core[1:4])
                    self.plane_surface_outer_air.append(gmsh.model.geo.addPlaneSurface([curve_loop_outer_air]))

            gmsh.model.geo.synchronize()
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if visualize_before or save_png:
                color_scheme = colors_femmt_default
                colors_geometry = colors_geometry_femmt_default

                # core color
                for i in range(0, len(self.plane_surface_core)):
                    gmsh.model.setColor([(2, self.plane_surface_core[i])], color_scheme[colors_geometry["core"]][0], 
                    color_scheme[colors_geometry["core"]][1], color_scheme[colors_geometry["core"]][2], recursive=True)

                # air gap color
                if self.plane_surface_air_gaps:
                    # only colorize air-gap in case of air gaps
                    for air_gap in self.plane_surface_air_gaps:
                        gmsh.model.setColor([(2, air_gap)], color_scheme[colors_geometry["air_gap"]][0], 
                            color_scheme[colors_geometry["air_gap"]][1], color_scheme[colors_geometry["air_gap"]][2], recursive=True)

                # air/potting-material inside core window
                gmsh.model.setColor([(2, self.plane_surface_air[0])], color_scheme[colors_geometry["potting_inner"]][0], 
                color_scheme[colors_geometry["potting_inner"]][1], color_scheme[colors_geometry["potting_inner"]][2], recursive=True)

                # winding colors
                for winding_number in range(0, self.component.n_windings):
                    for turn_number in range(0, len(self.plane_surface_cond[winding_number])):
                        gmsh.model.setColor([(2, self.plane_surface_cond[winding_number][turn_number])], 
                        color_scheme[colors_geometry["winding"][winding_number]][0], color_scheme[colors_geometry["winding"][winding_number]][1], 
                        color_scheme[colors_geometry["winding"][winding_number]][2], recursive=True)

                # isolation color (inner isolation / bobbin)
                gmsh.model.setColor([(2, iso) for iso in self.plane_surface_iso_core + self.plane_surface_iso_pri_sec], 
                    color_scheme[colors_geometry["isolation"]][0], color_scheme[colors_geometry["isolation"]][1], color_scheme[colors_geometry["isolation"]][2], recursive=True)

                if visualize_before:
                    gmsh.fltk.run()

                # Output .msh file
                # gmsh.option.setNumber("Mesh.SaveAll", 1)
                if save_png:
                    if '-nopopup' not in sys.argv:
                        gmsh.fltk.initialize()

                    gmsh.write(self.component.hybrid_color_visualize_file)  # save png

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # TODO The following algorithms try to modify the mesh in order to reduce the runtime. But maybe the synchronize() calls
            # have a high runtime. Check if thats true and when it does try to reduce the number of synchronize() calls by adding all points first and
            # embed them later together:
            # This is added here therefore the additional points are not seen in the pictures and views
            self.forward_meshing()

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # No mesh is generated here, because generating a mesh, saving it as *.msh, loading it, appending more geometry data
            # and then mesh again can cause bugs in the mesh
            # Therefore only the model geometry is saved and the mesh will be generated later
            # -> Save file as geo: File extension must be *.geo_unrolled
            gmsh.write(self.component.model_geo_file)

        def generate_electro_magnetic_mesh(self, refine = 0):
            print("Electro Magnetic Mesh Generation in Gmsh (write physical entities)")

            gmsh.open(self.component.model_geo_file)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Define physical Surfaces and Curves
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Core
            self.ps_core = gmsh.model.geo.addPhysicalGroup(2, self.plane_surface_core, tag=2000)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Conductors
            self.ps_cond = [[], []]
            for num in range(0, self.component.n_windings):

                if self.component.windings[num].conductor_type in [ConductorType.Foil, ConductorType.Solid, ConductorType.Full, ConductorType.Stacked]:
                    for i in range(0, sum(self.component.windings[num].turns)):
                        self.ps_cond[num].append(
                            gmsh.model.geo.addPhysicalGroup(2, [self.plane_surface_cond[num][i]], tag=4000 + 1000 * num + i))
                if self.component.windings[num].conductor_type == ConductorType.Litz:
                    for i in range(0, sum(self.component.windings[num].turns)):
                        self.ps_cond[num].append(
                            gmsh.model.geo.addPhysicalGroup(2, [self.plane_surface_cond[num][i]], tag=6000 + 1000 * num + i))

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Air, air_gaps and iso (since isolation is handled as air, as well as the air gaps)
            air_and_air_gaps = self.plane_surface_air + self.plane_surface_air_gaps + self.plane_surface_iso_core + self.plane_surface_iso_pri_sec
            self.ps_air = gmsh.model.geo.addPhysicalGroup(2, air_and_air_gaps, tag=1000)
            # ps_air_ext = gmsh.model.geo.addPhysicalGroup(2, plane_surface_outer_air, tag=1001)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Boundary
            self.pc_bound = gmsh.model.geo.addPhysicalGroup(1, self.l_bound_tmp, tag=1111)
            # print(f"Physical Conductor Surfaces: {ps_cond}")

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Set names [optional]
            gmsh.model.setPhysicalName(2, self.ps_core, "CORE")
            for num in range(0, self.component.n_windings):
                for i in range(0, len(self.ps_cond[num])):
                    gmsh.model.setPhysicalName(2, self.ps_cond[num][i], f"COND{num + 1}")
            gmsh.model.setPhysicalName(2, self.ps_air, "AIR")
            gmsh.model.setPhysicalName(1, self.pc_bound, "BOUND")

            # mshopt # Explicit stray path air gap optimization
            # mshopt if not self.component.component_type == "integrated_transformer":
            # mshopt     stray_path_mesh_optimizer = []
            # mshopt     stray_path_mesh_optimizer.append(gmsh.model.geo.addPoint(stray_path_gap[0][0],
            # stray_path_gap[0][1]+0.0001, stray_path_gap[0][2], 0.5*stray_path_gap[0][3]))
            # mshopt     stray_path_mesh_optimizer.append(gmsh.model.geo.addPoint(stray_path_gap[1][0],
            # stray_path_gap[1][1]-0.0001, stray_path_gap[1][2], 0.5*stray_path_gap[1][3]))
            # mshopt     print(f"plane_surface_core: {plane_surface_core}")
            # mshopt     print(f"stray_path_mesh_optimizer: {stray_path_mesh_optimizer}")
            # mshopt     print(f"stray_path_mesh_optimizer coordinates: {stray_path_gap[0][0], stray_path_gap[0][1],
            # stray_path_gap[0][2], stray_path_gap[0][3]}\n"
            # mshopt           f"{stray_path_gap[1][0], stray_path_gap[1][1], stray_path_gap[1][2],
            # stray_path_gap[1][3]}")

            # Synchronize
            gmsh.model.geo.synchronize()

            # Output .msh file
            # TODO: What are these flags about???
            gmsh.option.setNumber("Mesh.SaveAll", 1)
            gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
            gmsh.option.setNumber("Mesh.SurfaceFaces", 0)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # TODO: Adaptive Meshing
            if refine == 1:
                print("\n ------- \nRefined Mesh Creation ")
                # mesh the new gmsh.model using the size field
                bg_field = gmsh.model.mesh.field.add("PostView")
                # TODO: gmsh.model.mesh.field.setNumber(bg_field, "ViewTag", sf_view)
                gmsh.model.mesh.field.setAsBackgroundMesh(bg_field)
                print("\nMeshing...\n")
                gmsh.model.mesh.generate(2)
            else:
                print("\nMeshing...\n")
                gmsh.model.mesh.generate(2)

            if not os.path.exists(self.component.mesh_folder_path):
                os.mkdir(self.component.mesh_folder_path)

            gmsh.write(self.component.e_m_mesh_file)

        def generate_thermal_mesh(self, case_gap_top, case_gap_right, case_gap_bot, color_scheme, colors_geometry, visualize_before):
            print("Thermal Mesh Generation in Gmsh (write physical entities)")

            gmsh.open(self.component.model_geo_file)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Create case around the core

            # Core line tags
            core_line_tags = []
            if self.component.air_gaps.number == 0:
                core_line_tags = [3, 2, 1]
            else:
                core_line_tags = [4, 3, 2]

            # Core point tags
            core_point_tags = []
            if self.component.air_gaps.number == 0:
                core_point_tags = [4, 3, 2, 1]
            else:
                core_point_tags = [5, 4, 3, 2]

            tl_point = core_point_tags[0] # Top left - default 5
            tr_point = core_point_tags[1] # Top right - default 4
            br_point = core_point_tags[2] # Bottom right - default 3
            bl_point = core_point_tags[3] # Bottom left - default 2

            top_line = core_line_tags[0] # default 4
            right_line = core_line_tags[1] # default 3 
            bottom_line = core_line_tags[2] # default 2

            # Get positions from points
            tl_point_pos  = gmsh.model.getValue(0, tl_point, [])
            tr_point_pos  = gmsh.model.getValue(0, tr_point, [])
            br_point_pos  = gmsh.model.getValue(0, br_point, [])
            bl_point_pos  = gmsh.model.getValue(0, bl_point, [])

            mesh = self.c_core*4 # It typically does not need to be the same size as c_core, but it shouldn't be too big either

            # Create 5 new areas: top, top right, right, bottom right, bottom
            # top
            top_case_left_point = gmsh.model.geo.addPoint(tl_point_pos[0], tl_point_pos[1] + case_gap_top, tl_point_pos[2], mesh)
            top_case_right_point = gmsh.model.geo.addPoint(tr_point_pos[0], tr_point_pos[1] + case_gap_top, tr_point_pos[2], mesh)
            top_case_left_line = gmsh.model.geo.addLine(tl_point, top_case_left_point)
            top_case_top_line = gmsh.model.geo.addLine(top_case_left_point, top_case_right_point)
            top_case_right_line = gmsh.model.geo.addLine(top_case_right_point, tr_point)
            top_case_curve_loop = gmsh.model.geo.addCurveLoop([top_case_left_line, top_case_top_line, top_case_right_line, top_line])
            top_case_surface = gmsh.model.geo.addPlaneSurface([top_case_curve_loop])

            # top right
            top_right_case_top_right_point = gmsh.model.geo.addPoint(tr_point_pos[0] + case_gap_right, tr_point_pos[1] + case_gap_top, tr_point_pos[2], mesh)
            top_right_case_right_point = gmsh.model.geo.addPoint(tr_point_pos[0] + case_gap_right, tr_point_pos[1], tr_point_pos[2], mesh)
            top_right_case_bottom_line = gmsh.model.geo.addLine(tr_point, top_right_case_right_point)
            top_right_case_right_line = gmsh.model.geo.addLine(top_right_case_right_point, top_right_case_top_right_point)
            top_right_case_top_line = gmsh.model.geo.addLine(top_right_case_top_right_point, top_case_right_point)
            top_right_case_curve_loop = gmsh.model.geo.addCurveLoop([top_case_right_line, top_right_case_bottom_line, top_right_case_right_line, top_right_case_top_line])
            top_right_case_surface = gmsh.model.geo.addPlaneSurface([top_right_case_curve_loop])

            # right
            right_case_bottom_point = gmsh.model.geo.addPoint(br_point_pos[0] + case_gap_right, br_point_pos[1], br_point_pos[2], mesh)
            right_case_right_line = gmsh.model.geo.addLine(top_right_case_right_point, right_case_bottom_point)
            right_case_bottom_line = gmsh.model.geo.addLine(right_case_bottom_point, br_point)
            right_case_curve_loop = gmsh.model.geo.addCurveLoop([top_right_case_bottom_line, right_case_right_line, right_case_bottom_line, right_line])
            right_case_surface = gmsh.model.geo.addPlaneSurface([right_case_curve_loop])

            # bottom right
            bottom_right_case_bottom_right_point = gmsh.model.geo.addPoint(br_point_pos[0] + case_gap_right, br_point_pos[1] - case_gap_bot, br_point_pos[2], mesh)
            bottom_right_case_bottom_point = gmsh.model.geo.addPoint(br_point_pos[0], br_point_pos[1] - case_gap_bot, br_point_pos[2], mesh)
            bottom_right_case_left_line = gmsh.model.geo.addLine(br_point, bottom_right_case_bottom_point)
            bottom_right_case_bottom_line = gmsh.model.geo.addLine(bottom_right_case_bottom_point, bottom_right_case_bottom_right_point)
            bottom_right_case_right_line = gmsh.model.geo.addLine(bottom_right_case_bottom_right_point, right_case_bottom_point)
            bottom_right_case_curve_loop = gmsh.model.geo.addCurveLoop([right_case_bottom_line, bottom_right_case_left_line, bottom_right_case_bottom_line, bottom_right_case_right_line])
            bottom_right_case_surface = gmsh.model.geo.addPlaneSurface([bottom_right_case_curve_loop])

            # bottom
            bottom_case_bottom_left_point = gmsh.model.geo.addPoint(bl_point_pos[0], bl_point_pos[1] - case_gap_bot, bl_point_pos[2], mesh)
            bottom_case_bottom_line = gmsh.model.geo.addLine(bottom_right_case_bottom_point, bottom_case_bottom_left_point)
            bottom_case_left_line = gmsh.model.geo.addLine(bottom_case_bottom_left_point, bl_point)
            bottom_case_curve_loop = gmsh.model.geo.addCurveLoop([bottom_case_bottom_line, bottom_case_left_line, bottom_line, bottom_right_case_left_line])
            bottom_case_surface = gmsh.model.geo.addPlaneSurface([bottom_case_curve_loop])

            gmsh.model.geo.synchronize()

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Define physical Surfaces and Curves
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Core
            self.ps_core = gmsh.model.geo.addPhysicalGroup(2, self.plane_surface_core, tag=2000)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Conductors
            self.ps_cond = [[], []]
            for num in range(0, self.component.n_windings):
                if self.component.windings[num].conductor_type in [ConductorType.Foil, ConductorType.Solid, ConductorType.Full, ConductorType.Stacked]:
                    for i in range(0, sum(self.component.windings[num].turns)):
                        self.ps_cond[num].append(
                            gmsh.model.geo.addPhysicalGroup(2, [self.plane_surface_cond[num][i]], tag=4000 + 1000 * num + i))
                if self.component.windings[num].conductor_type == ConductorType.Litz:
                    for i in range(0, sum(self.component.windings[num].turns)):
                        self.ps_cond[num].append(
                            gmsh.model.geo.addPhysicalGroup(2, [self.plane_surface_cond[num][i]], tag=6000 + 1000 * num + i))

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Air
            self.ps_air = gmsh.model.geo.addPhysicalGroup(2, self.plane_surface_air, tag=1000)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Air gaps
            self.ps_air_gaps = gmsh.model.geo.addPhysicalGroup(2, self.plane_surface_air_gaps, tag=1001)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Isolations
            # TODO Currently isolations can only have the same material
            self.ps_isolation = gmsh.model.geo.addPhysicalGroup(2, self.plane_surface_iso_pri_sec + self.plane_surface_iso_core)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Boundary
            self.thermal_boundary_region_tags = {
                "BOUNDARY_TOP"            : top_case_top_line,
                "BOUNDARY_TOP_RIGHT"      : top_right_case_top_line,
                "BOUNDARY_RIGHT_TOP"      : top_right_case_right_line,
                "BOUNDARY_RIGHT"          : right_case_right_line,
                "BOUNDARY_RIGHT_BOTTOM"   : bottom_right_case_right_line,
                "BOUNDARY_BOTTOM_RIGHT"   : bottom_right_case_bottom_line,
                "BOUNDARY_BOTTOM"         : bottom_case_bottom_line
            }

            for key in self.thermal_boundary_region_tags:
                self.thermal_boundary_region_tags[key] = create_physical_group(1, [self.thermal_boundary_region_tags[key]], key)

            # Add surface physical groups
            # INFO: The physical groups are not created in the createRectWithPhysicalGroup because it causes a bug with the index counter when
            # 1D physical groups (lines) are added after 2D physical groups (surfaces)
            top_surface_physical_group = create_physical_group(2, [top_case_surface], "TopCase")
            top_right_surface_physical_group = create_physical_group(2, [top_right_case_surface], "TopRightCase")
            right_surface_physical_group = create_physical_group(2, [right_case_surface], "RightCase")
            bottom_right_surface_physical_group = create_physical_group(2, [bottom_right_case_surface], "BottomRightCase")
            bottom_surface_physical_group = create_physical_group(2, [bottom_case_surface], "BottomCase")

            self.thermal_boundary_ps_groups = [top_surface_physical_group, top_right_surface_physical_group, 
                right_surface_physical_group, bottom_right_surface_physical_group, bottom_surface_physical_group]

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Set names [optional]
            gmsh.model.setPhysicalName(2, self.ps_core, "CORE")
            for num in range(0, self.component.n_windings):
                for i in range(0, len(self.ps_cond[num])):
                    gmsh.model.setPhysicalName(2, self.ps_cond[num][i], f"COND{num + 1}")
            gmsh.model.setPhysicalName(2, self.ps_air, "AIR")
            gmsh.model.setPhysicalName(2, self.ps_air_gaps, "AIR_GAPS")
            gmsh.model.setPhysicalName(2, self.ps_isolation, "ISOLATIONS")
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # - Forward Meshing -
            p_inter = None
            # Inter Conductors
            for n_win in range(0, len(self.component.virtual_winding_windows)):
                if self.component.virtual_winding_windows[n_win].winding != WindingType.Interleaved:
                    for num in range(0, self.component.n_windings):
                        p_inter = []
                        x_inter = []
                        y_inter = []
                        j = 0

                        if self.component.windings[num].conductor_type == ConductorType.Solid and \
                                self.component.windings[num].turns[n_win] > 1:
                            while self.component.two_d_axi.p_conductor[num][5 * j][1] == \
                                    self.component.two_d_axi.p_conductor[num][5 * j + 5][1]:
                                x_inter.append(
                                    0.5 * (self.component.two_d_axi.p_conductor[num][5 * j][0] +
                                           self.component.two_d_axi.p_conductor[num][5 * j + 5][0]))
                                j += 1
                                if j == self.component.windings[num].turns[n_win] - 1:
                                    break
                            j += 1
                            if int(self.component.windings[num].turns[n_win] / j) > 1:
                                for i in range(0, int(self.component.windings[num].turns[n_win] / j)):
                                    if 5 * j * i + 5 * j >= len(self.component.two_d_axi.p_conductor[num][:]):
                                        break
                                    y_inter.append(0.5 * (self.component.two_d_axi.p_conductor[num][5 * j * i][1] +
                                                          self.component.two_d_axi.p_conductor[num][5 * j * i + 5 * j][1]))
                                for x in x_inter:
                                    for y in y_inter:
                                        p_inter.append(gmsh.model.geo.addPoint(x,
                                                                               y,
                                                                               0,
                                                                               self.c_center_conductor[num]))

            # Synchronize
            gmsh.model.geo.synchronize()

            # Conductor Center
            for num in range(0, self.component.n_windings):
                for i in range(0, int(len(self.p_cond[num]) / 5)):
                    gmsh.model.mesh.embed(0, [self.p_cond[num][5 * i + 0]], 2, self.plane_surface_cond[num][i])

            # Embed points for mesh refinement
            # Inter Conductors
            for n_win in range(0, len(self.component.virtual_winding_windows)):
                if self.component.virtual_winding_windows[n_win].winding != WindingType.Interleaved:
                    gmsh.model.mesh.embed(0, p_inter, 2, self.plane_surface_air[0])

            # Synchronize again
            gmsh.model.geo.synchronize()

            # Set colors:
            color_case = color_scheme[colors_geometry["potting_outer"]]
            color_background = color_scheme[colors_geometry["potting_inner"]]
            color_core = color_scheme[colors_geometry["core"]]
            color_air_gap = color_scheme[colors_geometry["air_gap"]]
            color_isolation = color_scheme[colors_geometry["isolation"]]

            # Set case color to core color
            for sf in [top_case_surface, top_right_case_surface, right_case_surface, bottom_right_case_surface, bottom_case_surface]:
                gmsh.model.setColor([(2, sf)], color_case[0], color_case[1], color_case[2], recursive=True)

            # core color
            for i in range(0, len(self.plane_surface_core)):
                gmsh.model.setColor([(2, self.plane_surface_core[i])], color_core[0], color_core[1], color_core[2], recursive=True)

            # air gap color
            if self.plane_surface_air_gaps:
                # only colorize air-gap in case of air gaps
                gmsh.model.setColor([(2, self.plane_surface_air[0]), (2, self.plane_surface_air_gaps[0])], color_air_gap[0], color_air_gap[1],
                                    color_air_gap[2], recursive=True)

            # air/potting-material inside core window
            gmsh.model.setColor([(2, self.plane_surface_air[0])], color_background[0], color_background[1], color_background[2], recursive=True)

            # winding colors
            for winding_number in range(0, self.component.n_windings):
                for turn_number in range(0, len(self.plane_surface_cond[winding_number])):
                    gmsh.model.setColor([(2, self.plane_surface_cond[winding_number][turn_number])], color_scheme[colors_geometry["winding"][winding_number]][0], 
                        color_scheme[colors_geometry["winding"][winding_number]][1], color_scheme[colors_geometry["winding"][winding_number]][2], recursive=True)

            # isolation color (inner isolation / bobbin)
            gmsh.model.setColor([(2, iso) for iso in self.plane_surface_iso_core + self.plane_surface_iso_pri_sec], color_isolation[0], color_isolation[1], 
                color_isolation[2], recursive=True)

            if visualize_before:
                gmsh.fltk.run()
            
            # Output .msh file
            gmsh.model.mesh.generate(2)

            if not os.path.exists(self.component.mesh_folder_path):
                os.mkdir(self.component.mesh_folder_path)

            gmsh.write(self.component.thermal_mesh_file)

        def mesh(self, frequency=None, skin_mesh_factor=None):
            self.component.high_level_geo_gen(frequency=frequency, skin_mesh_factor=skin_mesh_factor)
            if self.component.valid:
                self.generate_hybrid_mesh()
                self.generate_electro_magnetic_mesh()

        def forward_meshing(self):
            """
            In this function multiple techniques in order to raise the mesh density at certain points are applied.

            :return:
            """

            p_inter = None
            # Inter Conductors
            for n_win in range(0, len(self.component.virtual_winding_windows)):
                if self.component.virtual_winding_windows[n_win].winding != WindingType.Interleaved:
                    for num in range(0, self.component.n_windings):
                        p_inter = []
                        x_inter = []
                        y_inter = []
                        j = 0

                        if self.component.windings[num].conductor_type == ConductorType.Solid and \
                                self.component.windings[num].turns[n_win] > 1:
                            while self.component.two_d_axi.p_conductor[num][5 * j][1] == \
                                    self.component.two_d_axi.p_conductor[num][5 * j + 5][1]:
                                x_inter.append(
                                    0.5 * (self.component.two_d_axi.p_conductor[num][5 * j][0] +
                                           self.component.two_d_axi.p_conductor[num][5 * j + 5][0]))
                                j += 1
                                if j == self.component.windings[num].turns[n_win] - 1:
                                    break
                            j += 1
                            if int(self.component.windings[num].turns[n_win] / j) > 1:
                                for i in range(0, int(self.component.windings[num].turns[n_win] / j)):
                                    if 5 * j * i + 5 * j >= len(self.component.two_d_axi.p_conductor[num][:]):
                                        break
                                    y_inter.append(0.5 * (self.component.two_d_axi.p_conductor[num][5 * j * i][1] +
                                                          self.component.two_d_axi.p_conductor[num][5 * j * i + 5 * j][1]))
                                for x in x_inter:
                                    for y in y_inter:
                                        p_inter.append(gmsh.model.geo.addPoint(x,
                                                                               y,
                                                                               0,
                                                                               self.c_center_conductor[num]))

            # TODO: Inter conductor meshing!
            if all(winding.conductor_type == ConductorType.Solid for winding in self.component.windings):
                print(f"Making use of skin based meshing\n")
                for num in range(0, self.component.n_windings):
                    for i in range(0, int(len(self.p_cond[num]) / 5)):
                        gmsh.model.mesh.embed(0, [self.p_cond[num][5 * i + 0]], 2, self.plane_surface_cond[num][i])

                # Embed points for mesh refinement
                # Inter Conductors
                for n_win in range(0, len(self.component.virtual_winding_windows)):
                    if self.component.virtual_winding_windows[n_win].winding != WindingType.Interleaved:
                        gmsh.model.mesh.embed(0, p_inter, 2, self.plane_surface_air[0])
                # Stray path
                # mshopt gmsh.model.mesh.embed(0, stray_path_mesh_optimizer, 2, plane_surface_core[2])

            # Winding window rasterization:
            # In order adjust the mesh density in empty parts of the winding window a grid of possible points
            # is put on the winding window. Every point that is too close to the conductors is removed.
            # Every remaining point is added to the mesh with a higher mesh density

            min_distance = max([winding.conductor_radius for winding in self.component.windings]) + max(self.component.isolation.cond_cond)
            left_bound = self.component.core.core_w / 2
            right_bound = self.component.two_d_axi.r_inner
            top_bound = self.component.core.window_h / 2
            bot_bound = -self.component.core.window_h / 2

            width = right_bound - left_bound
            height = top_bound - bot_bound

            number_cols = 16 # Can be changed. More points equal higher raster density 
            number_rows = int(number_cols*height/width) # Assumption: number_cols/number_rows = width/height

            cell_width = width/(number_cols+1)
            cell_height =  height/(number_rows+1)

            # Get all possible points
            possible_points = []
            x = left_bound + cell_width/2
            y = bot_bound + cell_height/2
            for i in range(number_cols+1):
                for j in range(number_rows+1):
                    possible_points.append([x + i * cell_width, y + j * cell_height])

            fixed_points = []
            conductors = self.component.two_d_axi.p_conductor
            for winding in range(self.component.n_windings):
                for i in range(len(conductors[winding])//5):
                    point = conductors[winding][i*5]
                    fixed_points.append([point[0], point[1]])

            # Because the points need to be embed into the right surface. The points now will be split between different isolations and the air in the winding window.
            # TODO Currently primary secondary isolation is not considered
            left_iso = []
            right_iso = []
            top_iso = []
            bot_iso = []
            primary_secondary_iso = []
            air = []

            for i in range(len(self.component.two_d_axi.p_iso_pri_sec)):
                primary_secondary_iso.append([])

            # Isolations are currently not implemented for integrated transformers
            if self.component.component_type is not ComponentType.IntegratedTransformer:
                iso_core_left = self.component.two_d_axi.p_iso_core[0]
                iso_core_top = self.component.two_d_axi.p_iso_core[1]
                iso_core_right = self.component.two_d_axi.p_iso_core[2]
                iso_core_bot = self.component.two_d_axi.p_iso_core[3]

            # Extract all free_points
            for i in range(len(possible_points)):
                x = possible_points[i][0]
                y = possible_points[i][1]

                # Check collision with fixed points
                valid = True
                for fixed_point in fixed_points:
                    dist = np.sqrt((fixed_point[0]-x)**2 + (fixed_point[1]-y)**2)
                    if dist < min_distance:
                        valid = False
                        break

                if not valid:
                    continue

                # Check if point is in stray_path
                if self.component.component_type == ComponentType.IntegratedTransformer:
                    start_index = self.component.stray_path.start_index
                    stray_path_top_bound = self.component.air_gaps.midpoints[start_index+1][1] - self.component.air_gaps.midpoints[start_index+1][2] / 2
                    stray_path_bot_bound = self.component.air_gaps.midpoints[start_index][1] + self.component.air_gaps.midpoints[start_index][2] / 2
                    stray_path_right_bound = self.component.stray_path.length
                    stray_path_left_bound = left_bound

                    if x > stray_path_left_bound and x < stray_path_right_bound and y > stray_path_bot_bound and y < stray_path_top_bound:
                        continue

                # Point seems to be valid. Now find out in which surface the point belongs
                point = gmsh.model.geo.addPoint(x, y, 0, 2*self.c_window)

                if self.component.component_type is not ComponentType.IntegratedTransformer:
                    if point_is_in_rect(x, y, iso_core_left):
                        # Left iso
                        left_iso.append(point)
                    elif point_is_in_rect(x, y, iso_core_top):
                        # Top iso
                        top_iso.append(point)
                    elif point_is_in_rect(x, y, iso_core_right):
                        # Right iso
                        right_iso.append(point)
                    elif point_is_in_rect(x, y, iso_core_bot):
                        # Bot iso
                        bot_iso.append(point)
                    else:
                        # Check if point is in pri_to_sec isolation
                        added = False
                        for i in range(len(primary_secondary_iso)):
                            if point_is_in_rect(x, y, self.component.two_d_axi.p_iso_pri_sec[i]):
                                primary_secondary_iso[i].append(point)
                                added = True
                                break

                        if not added:
                            # Air
                            air.append(point)
                else:
                    air.append(point)
            # Call synchronize so the points will be added to the model            
            gmsh.model.geo.synchronize()

            # Embed points into surfaces
            if self.component.component_type is not ComponentType.IntegratedTransformer:
                gmsh.model.mesh.embed(0, left_iso, 2, self.plane_surface_iso_core[0])
                gmsh.model.mesh.embed(0, top_iso, 2, self.plane_surface_iso_core[1])
                gmsh.model.mesh.embed(0, right_iso, 2, self.plane_surface_iso_core[2])
                gmsh.model.mesh.embed(0, bot_iso, 2, self.plane_surface_iso_core[3])

                for i in range(len(primary_secondary_iso)):
                    gmsh.model.mesh.embed(0, primary_secondary_iso[i], 2, self.plane_surface_iso_pri_sec[i])
            
            gmsh.model.mesh.embed(0, air, 2, self.plane_surface_air[0])

            # Synchronize again
            gmsh.model.geo.synchronize()

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    # GetDP Interaction / Simulation / Excitation
    def excitation(self, frequency: float, amplitude_list: List, phase_deg_list: List = None, ex_type: str = 'current', imposed_red_f=0) -> None:
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

        :return: None
        :rtype: None

        """
        print(f"\n---\n"
              f"Excitation: \n"
              f"Frequency: {frequency}\n"
              f"Current(s): {amplitude_list}\n"
              f"Phase(s): {phase_deg_list}\n")


        # -- Excitation --
        self.flag_excitation_type = ex_type  # 'current', 'current_density', 'voltage'

        # Has the user provided a list of phase angles?
        phase_deg_list = phase_deg_list or []
        # phase_deg_list = np.asarray(phase_deg_list)

        for num in range(0, self.n_windings):

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
            for num in range(0, self.n_windings):
                if self.windings[num].conductor_type == ConductorType.Litz:
                    self.red_freq[num] = self.windings[num].strand_radius / self.delta
                elif self.windings[num].conductor_type == ConductorType.Solid:
                    self.red_freq[num] = self.windings[num].conductor_radius / self.delta
                else:
                    print("Reduced Frequency does not have a physical value here")
                    print(self.windings[num].conductor_type)
                    self.red_freq[num] = 1  # TODO: doesn't make sense like this -> rewrite fore conductor windings shape
        else:
            # DC case
            self.delta = 1e20  # random huge value
            self.red_freq[num] = 0

    def file_communication(self):
        """
        Interaction between python and Prolog files.

        :return:

        """
        # --------------------------------- File Communication --------------------------------
        # All shared control variables and parameters are passed to a temporary Prolog file
        print(f"\n---\n"
              f"File Communication\n")

        # Write initialization parameters for simulation in .pro file
        self.write_electro_magnetic_parameter_pro()

        # Write postprocessing parameters in .pro file
        self.write_electro_magnetic_post_pro()

    def write_electro_magnetic_parameter_pro(self):
        """
        Write materials and other parameters to the "Parameter.pro" file.
        This file is generated by python and is read by gmsh to hand over some parameters.

        Source for the parameters is the MagneticComponent object.
        :return: None
        :rtype: None
        """
        text_file = open(os.path.join(self.electro_magnetic_folder_path, "Parameter.pro"), "w")

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
        for num in range(0, self.n_windings):
            # -- Control Flags --
            if self.flag_excitation_type == 'current':
                text_file.write(f"Flag_ImposedVoltage = 0;\n")
            if self.flag_excitation_type == 'voltage':
                text_file.write(f"Flag_ImposedVoltage = 1;\n")
            if self.windings[num].conductor_type == ConductorType.Litz:
                text_file.write(f"Flag_HomogenisedModel{num + 1} = 1;\n")
            else:
                text_file.write(f"Flag_HomogenisedModel{num + 1} = 0;\n")

            # -- Geometry --
            # Number of conductors
            text_file.write(f"NbrCond{num + 1} = {sum(self.windings[num].turns)};\n")

            # For stranded Conductors:
            # text_file.write(f"NbrstrandedCond = {self.turns};\n")  # redundant
            if self.windings[num].conductor_type == ConductorType.Litz:
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

            print(f"Cell surface area: {self.windings[num].a_cell} \n"
                  f"Reduced frequency: {self.red_freq[num]}")

            if self.red_freq[num] > 1.25 and self.windings[num].conductor_type == ConductorType.Litz:
                # TODO: Allow higher reduced frequencies
                print(f"Litz Coefficients only implemented for X<=1.25")
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

        :return:
        """
        text_file = open(os.path.join(self.electro_magnetic_folder_path, "postquantities.pro"), "w") 

        # This is needed because the f string cant contain a \ in {}
        backslash = "\\"

        text_file.write(f"DirRes = \"{self.results_folder_path.replace(backslash, '/')}/\";\n")
        text_file.write(f"DirResFields = \"{self.e_m_fields_folder_path.replace(backslash, '/')}/\";\n")
        text_file.write(f"DirResVals = \"{self.e_m_values_folder_path.replace(backslash, '/')}/\";\n")
        text_file.write(f"DirResValsPrimary = \"{self.e_m_values_folder_path.replace(backslash, '/')}/Primary/\";\n")
        text_file.write(f"DirResValsSecondary = \"{self.e_m_values_folder_path.replace(backslash, '/')}/Secondary/\";\n")
        text_file.write(f"DirResCirc = \"{self.e_m_circuit_folder_path.replace(backslash, '/')}/\";\n")
        text_file.write(f"OptionPos = \"{self.results_folder_path.replace(backslash, '/')}/option.pos\";\n")

        # Visualisation
        if self.plot_fields == "standard":
            text_file.write(f"Flag_show_standard_fields = 1;\n")
        else:
            text_file.write(f"Flag_show_standard_fields = 0;\n")

        text_file.close()

    def simulate(self):
        """
        Initializes a onelab client. Provides the GetDP based solver with the created mesh file.

        :return:

        """
        print(f"\n---\n"
              f"Initialize ONELAB API\n"
              f"Run Simulation\n")

        # -- Simulation --
        # create a new onelab client

        # Initial Clearing of gmsh data
        gmsh.clear()

        # get model file names with correct path
        solver = os.path.join(self.electro_magnetic_folder_path, "ind_axi_python_controlled.pro")

        os.chdir(self.working_directory)

        # Run simulations as sub clients (non blocking??)
        mygetdp = os.path.join(self.onelab_folder_path, "getdp")
        self.onelab_client.runSubClient("myGetDP", mygetdp + " " + solver + " -msh " + self.e_m_mesh_file + " -solve Analysis -v2")

    def calculate_and_write_log(self, sweep_number: int = 1, currents: List = None, frequencies: List = None):
        """
        Method reads back the results from the .dat result files created by the ONELAB simulation client and stores
        them in a dictionary. Additionally, the input settings which are used in order to create the simulation are also printed.
        From this data type a JSON log file is created.
        :param sweep_number: Number of sweep iterations that were done before. For a single simulation sweep_number = 1
        :param currents: Current values of the sweep iterations. Not needed for single simulation
        :param frequencies: frequencies values of the sweep iterations. Not needed for single simulation
        :return:
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

            for winding in range(0, self.n_windings):

                # Create empty winding dictionary
                # create dictionary winding_dict with 'turn_losses' as list of the j=number_turns turn losses.
                # Single values related to one winding are added as 'winding_losses' etc.
                winding_dict = {"turn_losses": [],
                                "flux": [],
                                "self_inductance": [],
                                "mag_field_energy": [],
                                "V": []}

                # Number turns
                # TODO Since both windings can have a primary and a secondary turn, currently
                # only the turns corresponding to the winding number are printed:
                #   Winding1 -> primary turns
                #   Winding2 -> secondary turns
                # But technically Winding1 can have secondary turns and Winding2 can have primary turns
                if self.windings[winding].winding_type == WindingType.Interleaved:
                    winding_dict["number_turns"] = self.windings[winding].turns[0]
                else:
                    winding_dict["number_turns"] = self.windings[winding].turns[winding]

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
                if self.windings[winding].conductor_type == ConductorType.Litz:
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
            sweep_dict["all_winding_losses"] = sum(sweep_dict[f"winding{d+1}"]["winding_losses"] for d in range(0, self.n_windings))

            log_dict["single_sweeps"].append(sweep_dict)


        # Total losses of excitation sweep
        # Sum losses of all sweep runs. For core losses just use hyst_losses of the fundamental frequency.
        # Also needed as excitation for steady state thermal simulations

        # Single Windings
        for winding in range(0, self.n_windings):
            # Same as above
            turns = 0
            if self.windings[winding].winding_type == WindingType.Interleaved:
                turns = self.windings[winding].turns[0]
            else:
                turns = self.windings[winding].turns[winding]
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

        # ---- Print current configuration ----
        log_dict["simulation_settings"] = encode_settings(self)

        # ====== save data as JSON ======
        with open(self.e_m_results_log_path, "w+", encoding='utf-8') as outfile:
            json.dump(log_dict, outfile, indent=2, ensure_ascii=False)

    def read_log(self):
        log = {}
        with open(self.e_m_results_log_path, "r") as fd:
            content = json.loads(fd.read())
            log = content

        return log

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    # Post-Processing
    def visualize(self):
        """
        - a post simulation viewer
        - allows to open ".pos"-files in gmsh
        - For example current density, ohmic losses or the magnetic field density can be visualized

        :return:

        """
        # ---------------------------------------- Visualization in gmsh ---------------------------------------
        print(f"\n---\n"
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

        if any(self.windings[i].conductor_type != ConductorType.Litz for i in range(0, self.n_windings)):
            # Ohmic losses (weighted effective value of current density)
            gmsh.open(os.path.join(self.e_m_fields_folder_path, "j2F.pos"))
            gmsh.option.setNumber(f"View[{view}].ScaleType", 2)
            gmsh.option.setNumber(f"View[{view}].RangeType", 2)
            gmsh.option.setNumber(f"View[{view}].SaturateValues", 1)
            gmsh.option.setNumber(f"View[{view}].CustomMin", gmsh.option.getNumber(f"View[{view}].Min") + epsilon)
            gmsh.option.setNumber(f"View[{view}].CustomMax", gmsh.option.getNumber(f"View[{view}].Max"))
            gmsh.option.setNumber(f"View[{view}].ColormapNumber", 1)
            gmsh.option.setNumber(f"View[{view}].IntervalsType", 2)
            gmsh.option.setNumber(f"View[{view}].NbIso", 40)
            gmsh.option.setNumber(f"View[{view}].ShowTime", 0)
            print(gmsh.option.getNumber(f"View[{view}].Max"))
            view += 1

        if any(self.windings[i].conductor_type == ConductorType.Litz for i in range(0, self.n_windings)):
            # Ohmic losses (weighted effective value of current density)
            gmsh.open(os.path.join(self.e_m_fields_folder_path, "jH.pos"))
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
        gmsh.open(os.path.join(self.e_m_fields_folder_path, "Magb.pos"))
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
        gmsh.open(os.path.join(self.e_m_fields_folder_path, "raz.pos"))
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

    def get_loss_data(self, last_n_values, loss_type='litz_loss'):
        """
        Returns the last n values from the chosen loss type logged in the result folder.

        :param last_n_values:
        :param loss_type:

        :return:

        """
        loss_file = None
        # Loss file location
        if loss_type == 'litz_loss':
            loss_file = 'j2H.dat'
        if loss_type == 'solid_loss':
            loss_file = 'j2F.dat'
        # Read the logged losses corresponding to the frequencies
        with open(os.path.join(self.e_m_values_folder_path, loss_file), newline='') as fd:
            reader = csv.reader(fd)
            data = list(reader)
        return data[-last_n_values:-1] + [data[-1]]

    def load_result(self, res_name, res_type="value", last_n: int = 1, part="real", position: int = 0):
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
            res_path=self.e_m_values_folder_path
        if res_type=="circuit":
            res_path=self.e_m_circuit_folder_path


        with open(os.path.join(res_path, f"{res_name}.dat")) as fd:
            lines = fd.readlines()[-last_n:]

            if part == "real":
                result = [float(line.split(sep=' ')[1 + 2*position + 1]) for n, line in enumerate(lines)]
            if part == "imaginary":
                result = [float(line.split(sep=' ')[2 + 2*position + 1]) for n, line in enumerate(lines)]

            return result

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    # FEMM [alternative Solver]
    def femm_reference(self, freq, current, sign=None, non_visualize=0):
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

        :return:

        """
        if os.name == 'nt':
            install_pyfemm_if_missing()
            if not self.femm_is_imported:
                globals()["femm"] = __import__("femm")
                self.femm_is_imported = True
        else:
            raise Exception('You are using a computer that is not running windows. '
                            'This command is only executable on Windows computers.')


        self.create_folders(self.femm_folder_path)

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


        print(f"{self.core.permeability_type=}, {self.core.sigma=}")
        if self.core.permeability_type == PermeabilityType.FixedLossAngle:
            femm.mi_addmaterial('Ferrite', self.core.mu_rel, self.core.mu_rel, 0, 0, self.core.sigma/1e6, 0, 0, 1, 0, self.core.phi_mu_deg, self.core.phi_mu_deg)
        elif self.core.permeability_type == PermeabilityType.RealValue:
            femm.mi_addmaterial('Ferrite', self.core.mu_rel, self.core.mu_rel, 0, 0, self.core.sigma/1e6, 0, 0, 1, 0, self.core.phi_mu_deg, self.core.phi_mu_deg)
        else:
            femm.mi_addmaterial('Ferrite', self.core.mu_rel, self.core.mu_rel, 0, 0, self.core.sigma/1e6, 0, 0, 1, 0, 0, 0)
        femm.mi_addmaterial('Air', 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0)
        if self.windings[0].conductor_type == ConductorType.Litz:
            femm.mi_addmaterial('Copper', 1, 1, 0, 0, self.windings[0].cond_sigma/1e6, 0, 0, 1, 5, 0, 0, self.windings[0].n_strands,
                                2 * 1000 * self.windings[0].strand_radius)  # type := 5. last argument
            print(f"Number of strands: {self.windings[0].n_strands}")
            print(f"Diameter of strands in mm: {2 * 1000 * self.windings[0].strand_radius}")
        if self.windings[0].conductor_type == ConductorType.Solid:
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

        # femm.mi_drawrectangle(self.two_d_axi.p_window[4, 0]+self.isolation.core_cond,
        # self.two_d_axi.p_window[4, 1]+self.component.isolation.core_cond,
        # self.two_d_axi.p_window[7, 0]-self.isolation.core_cond, self.two_d_axi.p_window[7, 1]-self.isolation.core_cond)
        #
        # femm.mi_addblocklabel(self.two_d_axi.p_window[7, 0]-2*self.isolation.core_cond,
        # self.two_d_axi.p_window[7, 1]-2*self.isolation.core_cond)
        #
        # femm.mi_selectlabel(self.two_d_axi.p_window[7, 0]-2*self.isolation.core_cond,
        # self.two_d_axi.p_window[7, 1]-2*self.isolation.core_cond)
        #
        # femm.mi_setblockprop('Copper', 0, 1, 'icoil', 0, 0, 1)
        #
        # femm.mi_clearselected()

        for num in range(0, self.n_windings):
            if self.windings[num].conductor_type in [ConductorType.Litz, ConductorType.Solid]:
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
        femm.mi_saveas(os.path.join(self.femm_folder_path, 'coil.fem'))
        femm.mi_analyze()
        femm.mi_loadsolution()

        # == Log results ==
        self.write_femm_log()

        """
        # If we were interested in the flux density at specific positions,
        # we could inquire at specific points directly:
        b0 = femm.mo_getb(0, 0)
        print('Flux density at the center of the bar is %g T' % np.abs(b0[1]))
        b1 = femm.mo_getb(0.01, 0.05)
        print(f"Flux density at r=1cm, z=5cm is {np.abs(b1[1])} T")

        # The program will report the terminal properties of the circuit:
        # current, voltage, and flux linkage
        vals = femm.mo_getcircuitproperties('icoil')


        # [i, v, \[Phi]] = MOGetCircuitProperties["icoil"]

        # If we were interested in inductance, it could be obtained by
        # dividing flux linkage by current
        L = 1000 * np.abs(vals[2]) / np.abs(vals[0])
        print('The self-inductance of the coil is %g mH' % L)

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
        """

        :return:
        """
        if os.name == 'nt':
            install_pyfemm_if_missing()
            if not self.femm_is_imported:
                globals()["femm"] = __import__("femm")
                self.femm_is_imported = True
        else:
            raise Exception('You are using a computer that is not running windows. '
                            'This command is only executable on Windows computers.')

        file = open(self.femm_results_log_path, 'w+', encoding='utf-8')

        log = {}

        # print(hyst_loss)
        # tmp = femm.mo_getcircuitproperties('Primary')
        # print(tmp)
        # self.tot_loss_femm = 0.5 * tmp[0] * tmp[1]
        # print(self.tot_loss_femm)

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
    def calculate_point_average(x1, y1, x2, y2):
        """
        Used for femm_thermaL_validation
        """
        # TODO Move to femmt_functions
        return (x1 + x2) / 2, (y1 + y2) / 2

    def femm_thermal_validation(self, thermal_conductivity_dict, boundary_temperature, case_gap_top, case_gap_right, case_gap_bot):
        """
        Creates a thermal model in femm and simulates it with the given thermal conductivities

        :thermal_conductivity_dict: Dict containing conductivities for air, winding, case, core

        :return:
        """
        # Optional usage of FEMM tool by David Meeker
        # 2D Mesh and FEM interfaces (only for windows machines)
        if os.name == 'nt':
            install_pyfemm_if_missing()
            if not self.femm_is_imported:
                globals()["femm"] = __import__("femm")
                self.femm_is_imported = True
        else:
            raise Exception('You are using a computer that is not running windows. '
                'This command is only executable on Windows computers.')

        # Get paths
        femm_model_file_path = os.path.join(self.femm_folder_path, "thermal-validation.FEH")

        self.create_folders(self.femm_folder_path)

        # Extract losses
        losses = read_results_log(self.e_m_results_log_path)

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
        for num in range(0, self.n_windings):
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

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    # Litz Approximation [internal methods]
    def pre_simulate(self):
        """
        Used to determine the litz-approximation coefficients.

        :return:

        """
        for num in range(0, self.n_windings):
            if self.windings[num].conductor_type == ConductorType.Litz:
                # ---
                # Litz Approximation Coefficients were created with 4 layers
                # That's why here a hard-coded 4 is implemented
                # if os.path.isfile(self.path +
                # f"/Strands_Coefficients/coeff/pB_RS_la{self.windings[num].ff}_{self.n_layers[num]}layer.dat"):
                if os.path.exists(os.path.join(self.e_m_strands_coefficients_folder_path, "coeff", f"pB_RS_la{self.windings[num].ff}_4layer.dat")):
                    print("Coefficients for stands approximation are found.")

                else:
                    # Rounding X to fit it with corresponding parameters from the database
                    X = self.red_freq[num]
                    X = np.around(X, decimals=3)
                    print(f"Rounded Reduced frequency X = {X}")
                    self.create_strand_coeff(num)

    def create_strand_coeff(self, num):
        """

        :return:

        """
        print(f"\n"
              f"Pre-Simulation\n"
              f"-----------------------------------------\n"
              f"Create coefficients for strands approximation\n")
        # Create a new onelab client
        # -- Pre-Simulation Settings --
        text_file = open(os.path.join(self.e_m_strands_coefficients_folder_path, "PreParameter.pro"), "w")
        # ---
        # Litz Approximation Coefficients are created with 4 layers
        # That's why here a hard-coded 4 is implemented
        # text_file.write(f"NbrLayers = {self.n_layers[num]};\n")
        text_file.write(f"NbrLayers = 4;\n")
        text_file.write(f"Fill = {self.windings[num].ff};\n")
        print("Here")
        text_file.write(f"Rc = {self.windings[num].strand_radius};\n")  # double named!!! must be changed
        text_file.close()
        cell_geo = os.path.join(self.e_m_strands_coefficients_folder_path, "cell.geo")

        # Run gmsh as a sub client
        mygmsh = os.path.join(self.onelab_folder_path, "gmsh")
        self.onelab_client.runSubClient("myGmsh", mygmsh + " " + cell_geo + " -2 -v 2")

        modes = [1, 2]  # 1 = "skin", 2 = "proximity"
        reduced_frequencies = np.linspace(0, 1.25, 6)  # must be even
        for mode in modes:
            for rf in reduced_frequencies:
                # -- Pre-Simulation Settings --
                text_file = open(os.path.join(self.e_m_strands_coefficients_folder_path, "PreParameter.pro"), "w")
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
                input_file = os.path.join(self.e_m_strands_coefficients_folder_path, "cell_dat.pro")
                cell = os.path.join(self.e_m_strands_coefficients_folder_path, "cell.pro")

                # Run simulations as sub clients
                mygetdp = os.path.join(self.onelab_folder_path, "getdp")
                self.onelab_client.runSubClient("myGetDP", mygetdp + " " + cell + " -input " + input_file + " -solve MagDyn_a -v2")

        # Formatting stuff
        # Litz Approximation Coefficients are created with 4 layers
        # That's why here a hard-coded 4 is implemented
        # files = [self.path + f"/Strands_Coefficients/coeff/pB_RS_la{self.windings[num].ff}_{self.n_layers[num]}layer.dat",
        #         self.path + f"/Strands_Coefficients/coeff/pI_RS_la{self.windings[num].ff}_{self.n_layers[num]}layer.dat",
        #         self.path + f"/Strands_Coefficients/coeff/qB_RS_la{self.windings[num].ff}_{self.n_layers[num]}layer.dat",
        #         self.path + f"/Strands_Coefficients/coeff/qI_RS_la{self.windings[num].ff}_{self.n_layers[num]}layer.dat"]
        
        coeff_folder = os.path.join(self.e_m_strands_coefficients_folder_path, "coeff") 
        if not os.path.isdir(coeff_folder):
            os.mkdir(coeff_folder)

        files = [os.path.join(coeff_folder, f"pB_RS_la{self.windings[num].ff}_4layer.dat"),
                 os.path.join(coeff_folder, f"pI_RS_la{self.windings[num].ff}_4layer.dat"),
                 os.path.join(coeff_folder, f"qB_RS_la{self.windings[num].ff}_4layer.dat"),
                 os.path.join(coeff_folder, f"qI_RS_la{self.windings[num].ff}_4layer.dat")]
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

    def pre_simulation(self):
        """
        - Complete "pre-simulation" call

        :return:

        """
        self.high_level_geo_gen()
        self.excitation(frequency=100000, amplitude_list=1)  # arbitrary values: frequency and current
        self.file_communication()
        self.pre_simulate()

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    # Standard Simulations
    def create_model(self, freq: float, skin_mesh_factor: float = 0.5, visualize_before: bool = False,
                     save_png: bool = False, color_scheme: Dict = colors_femmt_default,
                     colors_geometry: Dict = colors_geometry_femmt_default):
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
            print("No air gaps are added")
        if self.windings is None or not self.windings:
            raise Exception("Winding classes need to be added to the magnetic component")
        if self.isolation is None:
            raise Exception("An isolation class need to be added to the magnetic component")
        if self.virtual_winding_windows is None:
            raise Exception("Virtual winding windows are not set properly. Please check the winding creation")

        self.high_level_geo_gen(frequency=freq, skin_mesh_factor=skin_mesh_factor)
        if self.valid:
            self.mesh.generate_hybrid_mesh(visualize_before=visualize_before, save_png=save_png)
        else:
            raise Exception("The model is not valid. The simulation won't start.")

    def single_simulation(self, freq: float, current: List[float], phi_deg: List[float] = None, show_results = True) -> None:
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

        :return: None

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

    def get_inductances(self, I0, op_frequency=0, skin_mesh_factor=1, visualize=False):
        """

        :param visualize:
        :param skin_mesh_factor:
        :param I0:
        :param op_frequency:

        :return:

        """

        # Remove "old" Inductance Logs
        try:
            os.remove(os.path.join(self.e_m_values_folder_path, "L_11.dat"))
            os.remove(os.path.join(self.e_m_values_folder_path, "L_22.dat"))
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

            print(f"\n"
                  f"                             == Inductances ==                             \n")

            # Read the logged Flux_Linkages
            with open(os.path.join(self.e_m_values_folder_path, "Flux_Linkage_1.dat")) as fd:
                line = fd.readlines()[-2:]
                # Fluxes induced in Winding 1
                Phi_11 = float(line[0].split(sep=' ')[2])
                Phi_12 = float(line[1].split(sep=' ')[2])

            with open(os.path.join(self.e_m_values_folder_path, "Flux_Linkage_2.dat")) as fd:
                line = fd.readlines()[-2:]
                # Fluxes induced in Winding 2
                Phi_21 = float(line[0].split(sep=' ')[2])
                Phi_22 = float(line[1].split(sep=' ')[2])

            print(f"\n"
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

            print(f"\n"
                  f"Turns Ratio:\n"
                  f"n = {n}\n"
                  )

            # Coupling Factors
            K_21 = Phi_21 / Phi_11
            K_12 = Phi_12 / Phi_22
            k = n / np.abs(n) * (K_21 * K_12) ** 0.5
            print(f"Coupling Factors:\n"
                  f"K_12 = Phi_12 / Phi_22 = {K_12}\n"
                  f"K_21 = Phi_21 / Phi_11 = {K_21}\n"
                  f"k = Sqrt(K_12 * K_21) = M / Sqrt(L_11 * L_22) = {k}\n"
                  )

            # Read the logged inductance values
            with open(os.path.join(self.e_m_values_folder_path, "L_11.dat")) as fd:
                line = fd.readlines()[-1]
                words = line.split(sep=' ')
                self.L_11 = float(words[2])
            with open(os.path.join(self.e_m_values_folder_path, "L_22.dat")) as fd:
                line = fd.readlines()[-1]
                words = line.split(sep=' ')
                self.L_22 = float(words[2])
            print(f"\n"
                  f"Self Inductances:\n"
                  f"L_11 = {self.L_11}\n"
                  f"L_22 = {self.L_22}\n"
                  )

            # Main/Counter Inductance
            self.M = k * (self.L_11 * self.L_22) ** 0.5
            M_ = self.L_11 * K_21  # Only to proof correctness - ideally: M = M_ = M__
            M__ = self.L_22 * K_12  # Only to proof correctness - ideally: M = M_ = M__
            print(f"\n"
                  f"Main/Counter Inductance:\n"
                  f"M = k * Sqrt(L_11 * L_22) = {self.M}\n"
                  f"M_ = L_11 * K_21 = {M_}\n"
                  f"M__ = L_22 * K_12 = {M__}\n"
                  )

            # Stray Inductance with 'Turns Ratio' n as 'Transformation Ratio' n
            L_s1 = self.L_11 - self.M * n
            L_s2 = self.L_22 - self.M / n
            L_h = self.M * n
            print(f"\n"
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

            print(f"\n"
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
            print(f"Invalid Geometry Data!")

    def excitation_sweep_old(self, frequencies, currents, phi, show_last=False):
        """
        Performs a sweep simulation for frequency-current pairs. Both values can
        be passed in lists of the same length. The mesh is only created ones (fast sweep)!

        :Example Code:
            #. geo = MagneticComponent()
            #. geo.mesh()
            #. fs = np.linspace(0, 250000, 6)
            #. cs = [10, 2, 1, 0.5, 0.2, 0.1]
            #. geo.excitation_sweep(frequencies=fs, currents=cs)

        :param phi:
        :param currents:
        :param frequencies:
        :param show_last:

        :return:

        """
        # frequencies = frequencies or []
        # currents = currents or []
        # phi = phi or []
        print(frequencies, currents, phi)
        for i in range(0, len(frequencies)):
            self.excitation(frequency=frequencies[i], amplitude_list=currents[i], phase_deg_list=phi[i])  # frequency and current
            self.file_communication()
            self.pre_simulate()
            self.simulate()
        if show_last:
            self.visualize()

    def excitation_sweep(self, frequency_list: List, current_list_list: List, phi_deg_list_list: List,
                         show_last: bool = False, return_results: bool = False, 
                         excitation_meshing_type: ExcitationMeshingType = None, skin_mesh_factor: float = 0.5, visualize_before: bool = False, save_png: bool = False,
                         color_scheme: Dict = colors_femmt_default, colors_geometry: Dict = colors_geometry_femmt_default) -> Dict:
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
                if winding.conductor_type == ConductorType.Solid:
                    excitation_meshing_type = ExcitationMeshingType.MeshEachFrequency
                    break
                if winding.conductor_type == ConductorType.Litz:
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

    def get_steinmetz_loss(self, Ipeak=None, ki=1, alpha=1.2, beta=2.2, t_rise=3e-6, t_fall=3e-6, f_switch=100000,
                           skin_mesh_factor=0.5):
        """

        :param skin_mesh_factor:
        :param Ipeak:
        :param ki:
        :param alpha:
        :param beta:
        :param t_rise:
        :param t_fall:
        :param f_switch:

        :return:

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


#  ===== Additional functions =====

def encode_settings(o: MagneticComponent):
    content = {
        "date": datetime.today().strftime('%Y-%m-%d %H:%M:%S'),
        "component_type": o.component_type.name,
        "working_directory": o.working_directory,
        "core": o.core.to_dict(),
        "air_gaps": o.air_gaps.to_dict(),
        "windings": [winding.to_dict() for winding in o.windings],
        "isolation": o.isolation.to_dict(),
        "virtual_winding_windows": [vww.to_dict() for vww in o.virtual_winding_windows],
    }

    if o.stray_path is not None:
        content["stray_path"] = o.stray_path.__dict__

    return content

def decode_settings_from_log(log_file_path: str, working_directory: str = None):
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

        air_gaps = AirGaps(AirGapMethod[settings["air_gaps"]["method"]], core)
        for air_gap in settings["air_gaps"]["air_gaps"]:
            air_gaps.add_air_gap(AirGapLegPosition[air_gap["leg_position"]], air_gap["position_value"], air_gap["height"])
        geo.set_air_gaps(air_gaps)

        windings = []
        for settings_winding in settings["windings"]:
            settings_cond = settings_winding["conductor_settings"]
            winding = Winding(settings_winding["turns_primary"], settings_winding["turns_secondary"], 
                                Conductivity[settings_winding["conductivity"]], WindingType[settings_winding["winding_type"]], 
                                WindingScheme[settings_winding["winding_scheme"]])
            conductor_type = ConductorType[settings_cond["conductor_type"]]
            if conductor_type == ConductorType.Foil:
                winding.set_foil_conductor(settings_cond["thickness"], settings_cond["wrap_para"])
            elif conductor_type == ConductorType.Full: 
                winding.set_full_conductor(settings_cond["thickness"], settings_cond["wrap_para"])
            elif conductor_type == ConductorType.Stacked:
                winding.set_stacked_conductor(settings_cond["thickness"], settings_cond["wrap_para"])
            elif conductor_type == ConductorType.Litz:
                winding.set_litz_conductor(settings_cond["conductor_radius"], settings_cond["n_strands"], settings_cond["strand_radius"],
                                            settings_cond["ff"])
            elif conductor_type == ConductorType.Solid:
                winding.set_solid_conductor(settings_cond["conductor_radius"])
            else:
                raise Exception(f"Unknown conductor type {conductor_type}")

            windings.append(winding)

        geo.set_windings(windings)

        isolation = Isolation()
        isolation.add_core_isolations(*settings["isolation"]["core_isolations"])
        isolation.add_winding_isolations(*settings["isolation"]["winding_isolations"])
        geo.set_isolation(isolation)

        return geo

    raise Exception(f"Couldn't extract settings from file {log_file_path}")