# Usual Python libraries
import csv
import fileinput
import numpy as np
import os
import pathlib
import sys
from matplotlib import pyplot as plt
from scipy.optimize import brentq
import gmsh
from onelab import onelab
import json
from scipy.integrate import quad
from scipy.interpolate import interp1d
import warnings
from thermal import thermal as th
from thermal import thermal_functions as th_functions
import shutil
# import pandas as pd
# import re
# import time
from typing import List, Union, Optional
from .femmt_functions import inner_points, min_max_inner_points, call_for_path, id_generator, NbrStrands, NbrLayers, \
    fft, compare_fft_list, r_basis, sigma, r_round_inf, r_round_round, r_cyl_cyl, r_cheap_cyl_cyl, \
    install_femm_if_missing, calculate_reluctances, r_cyl_cyl_real
from .Analytical_Core_Data import f_N95_mu_imag, f_N95_er_imag

# Optional usage of FEMM tool by David Meeker
# 2D Mesh and FEM interfaces (only for windows machines)
if os.name == 'nt':
    install_femm_if_missing()
    import femm


#  ===== Main Class  =====
class MagneticComponent:
    """
    A MagneticComponent is the main object for all simulation purposes in femmt.

        - One or more "MagneticComponents" can be created
        - Each "MagneticComponent" owns its own instance variable values

    """
    # Initialization of all class variables
    # Common variables for all instances

    # -- Parent folder path --
    path = str(pathlib.Path(__file__).parent.absolute())  # Path of FEMMT files
    onelab = None  # Path to Onelab installation folder
    path_mesh = "mesh/"
    path_res = "results/"
    path_res_fields = "results/fields/"
    path_res_values = "results/values/"
    path_res_circuit = "results/circuit/"
    path_res_FEMM = "FEMM/"

    def __init__(self, component_type="inductor"):
        """
        :param component_type: Available options:
                               - "inductor"
                               - "transformer"
                               - "integrated_transformer" (Transformer with included stray-path)
        :type component_type: string
        """
        print(f"\n"
              f"Initialized a new Magnetic Component of type {component_type}\n"
              f"--- --- --- ---")

        # Initialization of all instance variables

        # Breaking variable
        self.valid = True
        self.mu0 = np.pi * 4e-7

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Geometry

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Control Flags
        self.visualize_before = False
        self.region = None  # Apply an outer Region or directly apply a constraint on the Core Boundary
        self.padding = 1.5  # ... > 1
        self.y_symmetric = 1  # Mirror-symmetry across y-axis
        self.dimensionality = "2D axi"  # Axial-symmetric model (idealized full-cylindrical)
        self.s = 0.5  # Parameter for mesh-accuracy
        self.component_type = component_type  # "inductor", "transformer", "integrated_transformer" or
        # "three-phase-transformer"
        self.plot_fields = "standard"

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Core
        self.core = self.Core(self)
        self.core.update(type="EI", core_w=0.02, window_w=0.01, window_h=0.03)  # some initial values

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Air Gaps
        self.air_gaps = self.AirGaps(self)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Windings
        self.n_windings = None  # Number of conductors/windings
        if component_type == "inductor":
            self.n_windings = 1
            self.windings = [self.Winding()]

        if component_type == "transformer":
            self.n_windings = 2
            self.windings = [self.Winding(), self.Winding()]

        if component_type == "integrated_transformer":
            self.n_windings = 2
            self.windings = [self.Winding(), self.Winding()]
            self.stray_path = self.StrayPath()

        if component_type == "three_phase_transformer":
            self.n_windings = 3
            self.windings = [self.Winding(), self.Winding(), self.Winding()]
            raise NotImplemented

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Virtual Winding Windows
        self.virtual_winding_windows = None
        self.vw_type = None  # "center" and "full_window" are the only cases implemented yet; #TODO: replace

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Isolation
        self.isolation = self.Isolation()

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Geometric Parameters/Coordinates
        self.n_windows = None

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Excitation Parameters
        self.flag_imposed_reduced_frequency = None
        self.flag_excitation_type = None
        self.current = [None] * self.n_windings  # Defined for every conductor
        self.current_density = [None] * self.n_windings  # Defined for every conductor
        self.voltage = [None] * self.n_windings  # Defined for every conductor
        self.frequency = None
        self.phase_tmp = np.zeros(self.n_windings)  # Default is zero, Defined for every conductor
        self.red_freq = None  # [] * self.n_windings  # Defined for every conductor
        self.delta = None

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Materials
        self.core.steinmetz_loss = 0
        self.core.generalized_steinmetz_loss = 0
        self.Ipeak = None
        self.ki = None
        self.alpha = None
        self.beta = None
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

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   -  -  -  -  -  -  -  -  -  -  -
    # Setup
    def onelab_setup(self) -> None:
        """

        Either reads ONELAB parent folder path from config.json or asks the user to provide the ONELAB path it.
        Creates a config.json inside the site-packages folder at first run.

        :return: -
        """
        # find out path of femmt (installed module or directly opened in git)?
        module_file_path = os.path.dirname(__file__)
        config_file_path = os.path.join(module_file_path, 'config.json')

        # check if config.json is available and not empty
        if os.path.isfile(config_file_path) and os.stat(config_file_path).st_size != 0:
            path = ""
            with open(config_file_path, "r") as fd:
                loaded_dict = json.loads(fd.read())
                path = loaded_dict['onelab']

            if os.path.exists(path):
                self.onelab = path
            else:
                self.onelab = call_for_path("onelab")
        else:
            self.onelab = call_for_path("onelab")

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   -  -  -  -  -  -  -  -  -  -  -
    # Geometry Parts
    def high_level_geo_gen(self, core_type="EI", dimensionality="2D axi", frequency=None, skin_mesh_factor=1.0):
        """
        - high level geometry generation
        - based on chosen core and conductor types and simulation mode
        - calls "low level" methods, that creates all points needed for mesh generation
     
        :return:
     
        """
        # Always reset the to valid
        self.valid = True

        # ==============================
        # High-Level Geometry Generation
        # ==============================

        # Mesh-Parameters must be updated depending on geometry size
        self.mesh.c_core = self.core.core_w / 10. * self.s
        self.mesh.c_window = self.core.window_w / 30 * self.s

        # Update Skin Depth (needed for meshing)
        self.mesh.skin_mesh_factor = skin_mesh_factor
        if frequency is not None:
            if frequency == 0:
                self.delta = 1e9
            else:
                self.delta = np.sqrt(2 / (2 * frequency * np.pi * self.windings[0].cond_sigma * self.mu0))
            for i in range(0, self.n_windings):
                if self.windings[i].conductor_type == "solid":
                    self.mesh.c_conductor[i] = min([self.delta * self.mesh.skin_mesh_factor,
                                                    self.windings[i].conductor_radius / 4 * self.mesh.skin_mesh_factor])
                    self.mesh.c_center_conductor[i] = self.windings[i].conductor_radius / 4 * self.mesh.skin_mesh_factor
                elif self.windings[i].conductor_type == "litz":
                    self.mesh.c_conductor[i] = self.windings[i].conductor_radius / 4 * self.mesh.skin_mesh_factor
                    self.mesh.c_center_conductor[i] = self.windings[i].conductor_radius / 4 * self.mesh.skin_mesh_factor
                else:
                    self.mesh.c_conductor[i] = 0.0001  # revisit

        # -- Core-type --
        if self.core.type == core_type:
            self.n_windows = 2

        # -- Symmetry -- [Choose between asymmetric, symmetric and axi symmetric]
        self.dimensionality = dimensionality
        if self.dimensionality == "2D axi":
            self.y_symmetric = 1

        if self.core.type == "EI" and self.dimensionality == "2D axi":
            self.ei_axi = self.EIaxi(self)
            self.ei_axi.update()

    class VirtualWindingWindow:
        """
        A virtual winding window is the area, where either some kind of interleaved conductors or a one winding
        (primary, secondary,...) is placed in a certain way.
        """

        def __init__(self):
            # Rectangular frame:
            self.bot_bound = None
            self.top_bound = None
            self.left_bound = None
            self.right_bound = None

            # Arrangement of the Conductors in the virtual winding window
            # Obviously depends on the chosen conductor type
            self.winding = None  # "interleaved" | "primary"/"secondary"
            self.scheme = None  # "bifilar", "vertical", "horizontal", ["hexa", "square"] | "hexa", "square"
            self.turns = [None]  # Must be an ordered list with the number of turns of the conductors in the window
            # [N_primary, N_secondary,...]

    class Winding:
        """
        A winding defines a conductor which is wound around a magnetic component such as transformer or inductance.
        The winding is defined by its conductor and the way it is placed in the magnetic component. To allow different
        arrangements of the conductors in several winding windows (hexagonal or square packing, interleaved, ...) in
        this class only the conductor parameters are specified. Then, by calling class:Winding in
        class:VirtualWindingWindow the arrangement of the conductors is specified.
        """

        def __init__(self):
            self.cond_sigma = 5.8e7  # perfect copper
            self.turns = None
            self.conductor_type = None  # List of possible conductor types
            self.ff = None
            self.n_layers = None
            self.n_strands = None
            self.strand_radius = None
            self.conductor_radius = None
            self.a_cell = None
            self.thickness = None
            self.wrap_para = None

    class Core:
        """
        frequency = 0: mu_rel only used if non_linear == False
        frequency > 0: mu_rel is used

        """

        def __init__(self, component, re_mu_rel=3000):
            """
            :param re_mu_rel:
            
            """
            self.component = component  # Convention: parent

            # Complex Loss
            # TDK N95 as standard material:
            self.re_mu_rel = re_mu_rel  # Real part of relative Core Permeability  [B-Field and frequency-dependent]
            # self.im_mu_rel = 1750 * np.sin(10 * np.pi / 180)  # Imaginary part of relative Core Permeability
            # B-Field and frequency-dependent:
            self.im_mu_rel = 2500 * np.sin(20 * np.pi / 180)  # Imaginary part of relative Core Permeability
            self.im_epsilon_rel = 6e+4 * np.sin(
                20 * np.pi / 180)  # Imaginary part of complex equivalent permeability  [only frequency-dependent]
            self.material = 95_100  # 95 := TDK-N95 | Currently only works with Numbers corresponding to BH.pro

            # Steinmetz Loss
            self.steinmetz_loss = 0
            self.Ipeak = None
            self.ki = None
            self.alpha = None
            self.beta = None
            self.t_rise = None
            self.t_fall = None
            self.f_switch = None

            # Dimensions
            self.type = "EI"  # Basic shape of magnetic conductor
            self.core_w = None  # Axi symmetric case | core_w := core radius
            self.window_w = None  # Winding window width
            self.window_h = None  # Winding window height

        def update(self, type: str,
                   re_mu_rel: float = 3000,
                   im_mu_rel: float = 2500 * np.sin(20 * np.pi / 180),
                   im_epsilon_rel: float = 6e+4 * np.sin(20 * np.pi / 180),
                   material=95,
                   non_linear: bool = False,
                   **kwargs) -> None:
            """
            # TODO: Steinmetz - parameters
            Updates the core structure.

            - One positional parameter: type
            - All other core parameters can be passed or adjusted by keyword calling

                - Allows single parameter changing, depending on core-type
                - Strict keyword usage!
                - All dimensions are in meters

            :param im_epsilon_rel:
            :type im_epsilon_rel: float
            :param non_linear: True/False
            :type non_linear: bool
            :param material: specified in BH.pro. Currently available: '95_100' #ToDo: Add Values!
            :type material: str
            :param re_mu_rel: material's real mu_rel part, found in the datasheet (used if non_linear == False)
            :type re_mu_rel: float
            :param im_mu_rel: material's imaginary mu_rel part, found in the datasheet (used if non_linear == False)
            :type im_mu_rel: float
            :param type: There is actually only one type: "EI"
            :type type: str
            :param kwargs:
                - Case "2D, axisym., EI": 'core_w', 'window_w', 'window_h'
                - Case "3D, EI": ...tba...

            :return: None

            """

            print(f"Update the magnetic Core to {self.type}-type with following parameters: {kwargs}\n"
                  f"---")

            # Material Properties
            self.non_linear = non_linear
            self.re_mu_rel = re_mu_rel
            self.im_mu_rel = im_mu_rel
            self.im_epsilon_rel = im_epsilon_rel
            self.material = material
            self.type = type

            # Set attributes of core with given keywords
            for key, value in kwargs.items():
                setattr(self, key, value)

    class StrayPath:
        """

        """

        def __init__(self):
            # Dedicated Stray Path
            self.start_index = None  # lower air gap that characterizes the stray path
            self.radius = None
            self.width = None
            self.midpoint = None
            # TODO: Thickness of the stray path must be fitted for the real Tablet (effective area of the
            #  "stray air gap" is different in axi-symmetric approximation

        def update(self, **kwargs):
            # Set attributes of core with given keywords
            for key, value in kwargs.items():
                setattr(self, key, value)

    class AirGaps:
        """

        """

        def __init__(self, component):
            """

            """
            self.component = component
            # self.number = 1
            self.number = None  # Number of air gaps [==1: air gap in center | >1: random air gaps]

            # self.midpoints = np.empty((self.number, 4))
            self.midpoints = None  # list: [position_tag, air_gap_position, air_gap_h, c_air_gap]

        def update(self, method: str = None, n_air_gaps: Union[int, None] = None,
                   position_tag: List[float] = None,
                   air_gap_position: List[float] = None,
                   air_gap_h: List[float] = None, **kwargs) -> None:
            """
            Updates the air gap structure.

            - Strict keyword usage!
            - All dimensions are in meters
            - all parameters are in lists!
            - first chose the method, second, transfer parameters:

                - "center": ONE air gap exactly in core's middle
                - "random": random count of air gaps in the inner/outer leg
                - "percent": Easy way to split air gaps over the inner/outer leg
                - "manually": Place air gaps manually

            - "self.midpoints" is a list with [position_tag, air_gap_position, air_gap_h, c_air_gap]
            -  c_air_gap: mesh accuracy factor
            - "EI 2D axi": position_tag = 0  # '-1': left leg | '0': center leg | '1': right leg

            :param n_air_gaps: number of air gaps
            :type n_air_gaps: int
            :param position_tag: specifies the gapped "leg"
            :type position_tag: float
            :param air_gap_h: List of air gap high, list length depending on total air gap count
            :type air_gap_h: float
            :param air_gap_position: specifies the coordinate of the air gap's center point along the specified leg
            :type air_gap_position: float
            :param method: "random", "center", "percent", "manually"
            :type method: str

            :return: None

            :Example:

            >>> geo.air_gaps.update(method="center", n_air_gaps=1, air_gap_h=[0.002]) # 'center': single air gap in the middle
            >>> geo.air_gaps.update(method="percent", n_air_gaps=2, position_tag=[0, 0], air_gap_h=[0.003, 0.001], air_gap_position=[10, 80]) # 'percent': Place air gaps manually using percentages
            >>> geo.air_gaps.update(method="percent", n_air_gaps=2, position_tag=[0, 0], air_gap_h=[0.003, 0.001], air_gap_position=[10, 80]) # random
            >>> geo.air_gaps.update(method="manually", n_air_gaps=2, position_tag=[0, 0], air_gap_h=[0.003, 0.001], air_gap_position=[0.000, 0.003]) # manually

            """
            self.number = n_air_gaps
            position_tag = position_tag
            air_gap_position = air_gap_position
            air_gap_h = air_gap_h

            print(f"Update the air gaps.\n"
                  f"---")

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - -

            # Update the mesh accuracy of the window
            self.component.mesh.c_window = self.component.core.window_w / 20 * self.component.s

            # Rewrite variables
            self.midpoints = np.empty((self.number, 4))
            self.component.mesh.c_air_gap = [None] * self.number

            # Update air gaps with chosen method

            # Center
            if method == "center" and self.component.dimensionality == "2D axi":
                if self.number > 1:
                    print(f"{self.number} are too many air gaps for the 'center' option!")
                    raise Warning
                else:
                    self.component.mesh.c_air_gap[0] = air_gap_h[0] * self.component.s
                    self.midpoints[0, :] = np.array([0, 0, air_gap_h[0], self.component.mesh.c_air_gap[0]])

            # Deterministic
            if (method == "manually" or method == "percent") and self.component.dimensionality == "2D axi":
                for i in range(0, self.number):
                    if method == "percent":
                        # air_gap_position[i] = air_gap_position[i] / 100 * \
                        #                       (self.component.core.window_h - air_gap_h[i]) \
                        #                       - (self.component.core.window_h / 2 - air_gap_h[i] / 2)
                        air_gap_position[i] = air_gap_position[i] / 100 * self.component.core.window_h \
                                              - self.component.core.window_h / 2

                    # Overlapping Control
                    for j in range(0, self.midpoints.shape[0]):
                        if self.midpoints[j, 1] + self.midpoints[j, 2] / 2 > air_gap_position[i] > self.midpoints[
                            j, 1] - \
                                self.midpoints[j, 2] / 2:
                            if position_tag[i] == self.midpoints[j, 0]:
                                print(f"Overlapping Air Gap")
                                # raise Warning
                        else:
                            # self.component.mesh.c_air_gap[i] = air_gap_h[i] * self.component.s
                            self.component.mesh.c_air_gap[i] = self.component.mesh.c_window
                            # print(f"c_window: {self.component.mesh.c_window}")
                            self.midpoints[i, :] = np.array([position_tag[i],
                                                             air_gap_position[i],
                                                             air_gap_h[i],
                                                             self.component.mesh.c_air_gap[i]])

            #  TODO: Proof whether air gaps are valid

            # Random
            if method == "random" and self.component.dimensionality == "2D axi":
                position_tag = [0] * self.number

                i = 0
                while i in range(0, self.number):
                    height = np.random.rand(1) * 0.001 + 0.001
                    position = np.random.rand(1) * (self.component.core.window_h - height) - (
                            self.component.core.window_h / 2 - height / 2)
                    self.component.mesh.c_air_gap[i] = height * self.component.s
                    # Overlapping Control
                    for j in range(0, self.midpoints.shape[0]):
                        if self.midpoints[j, 1] + self.midpoints[j, 2] / 2 > position > self.midpoints[j, 1] - \
                                self.midpoints[j, 2] / 2:
                            if position_tag[i] == self.midpoints[j, 0]:
                                print(f"Overlapping air Gaps have been corrected")
                    else:
                        self.midpoints[i, :] = np.array([position_tag[i],
                                                         position,
                                                         height,
                                                         self.component.mesh.c_air_gap[i]])
                        i += 1

    class Isolation:
        """

        """

        def __init__(self, cond_cond=None, core_cond=None):
            """
            - Isolation
                - Between two turns of common conductors: first n_conductor arguments of cond_cond
                - Between two neighboured conductors: last n_conductor-1 arguments

            """
            self.cond_cond = cond_cond or []
            self.core_cond = core_cond or []

    # Update Methods
    def update_conductors(self, n_turns=None, conductor_type=None, winding=None, scheme=None, conductor_radii=None,
                          litz_para_type=None, ff=None, strands_numbers=None, strand_radii=None, thickness=None,
                          wrap_para=None, cond_cond_isolation=None, core_cond_isolation=None) -> None:
        """
        This Method allows the user to easily update/initialize the Windings in terms of their conductors and
        arrangement.

        Note: set all parameters in a list!

        :param wrap_para:
        :param ff: fill-factor, values between [0....1]
        :param winding:
                - "interleaved"
                - ""
        :param scheme:
        :param litz_para_type: there is one degree of freedom.
                - "implicit_litz_radius":
                - "implicit_ff":
                - "implicit_strands_number":

        :param thickness:
        :param cond_cond_isolation:
        :param core_cond_isolation:
        :param n_turns:
        :param strand_radii:
        :param strands_numbers:
        :param conductor_radii:
        :param conductor_type:  - "stacked"  # Vertical packing of conductors
                                - "full"     # One massive Conductor in each window
                                - "foil"     # Horizontal packing of conductors
                                - "solid"    # Massive wires
                                - "litz"     # Litz wires
        :return:

        """
        n_turns = n_turns or []
        conductor_type = conductor_type or []
        winding = winding or []
        scheme = scheme or []
        conductor_radii = conductor_radii or []
        litz_para_type = litz_para_type or []
        ff = ff or []
        strands_numbers = strands_numbers or []
        strand_radii = strand_radii or []
        thickness = thickness or []
        wrap_para = wrap_para or []
        cond_cond_isolation = cond_cond_isolation or []
        core_cond_isolation = core_cond_isolation or []

        print(f"Update the conductors...\n"
              f"---")

        # - - - - - - - - - - - - - - - - - - - - - CHECK Input Parameters - - - - - - - - - - - - - - - - - - - - - - -
        if self.component_type == "inductor":
            if len(n_turns) != 1 or len(conductor_type) != 1:
                print(f"Wrong number of conductor parameters passed to inductor model!")
                raise Warning

        if self.component_type == ("transformer" or "integrated_transformer"):
            if len(n_turns) != 2 or len(conductor_type) != 2:
                raise Warning(f"Wrong number of conductor parameters passed to transformer model!")

        # - - - - - - - - - - - - - - - - - Definition of Virtual Winding Windows - - - - - - - - - - - - - - - - - - -
        # Inductor: One virtual winding window
        if self.component_type == "inductor":
            self.vw_type = "full_window"  # One Virtual Winding Window

        # Two winding transformer
        # One virtual winding window
        if not self.component_type == "integrated_transformer":
            if len(winding) == 1 and winding[0] == "interleaved":
                self.vw_type = "full_window"  # One Virtual Winding Window
            if len(winding) == 2:
                self.vw_type = "center"  # Two Virtual Winding Windows #TODO: Adjust center (no longer just zero line)
        else:
            pass
            # Dedicated stray path: Two virtual winding windows
            # TODO: Subdivision into further Virtual Winding Windows

        self.virtual_winding_windows = []
        for windows in range(0, len(winding)):
            # If there is a dedicated stray path, always two Virtual Winding Windows are used!
            self.virtual_winding_windows.append(self.VirtualWindingWindow())

        # Fill the Virtual Winding Windows with the given data
        for vww in range(0, len(self.virtual_winding_windows)):
            self.virtual_winding_windows[vww].winding = winding[vww]
            self.virtual_winding_windows[vww].scheme = scheme[vww]
            self.virtual_winding_windows[vww].turns = list(map(list, zip(*n_turns)))[vww]
            # Need number of turns per VWW but given is a form
            # of [list_of_primary_turns, list_of_secondary_turns]

        # - - - - - - - - - - - - - - - - - Definition of the Isolation Parameters - - - - - - - - - - - - - - - - - - -
        # isolation parameters as lists:
        self.isolation.core_cond = core_cond_isolation
        self.isolation.cond_cond = cond_cond_isolation

        # - - - - - - - - - - - - - - - - - Definition of the Conductor Parameters - - - - - - - - - - - - - - - - - - -
        # self.n_windings is implied by component type

        for i in range(0, self.n_windings):
            # general conductor parameters
            self.windings[i].turns = n_turns[i]  # turns in a list which corresponds to the Virtual Winding Windows
            self.windings[i].conductor_type = conductor_type[i]

            if self.windings[i].conductor_type == ('stacked' or 'full' or 'foil'):
                # foil/stacked parameters
                self.windings[i].thickness = thickness[i]
                self.windings[i].wrap_para = wrap_para[i]

            # round conductors
            if self.windings[i].conductor_type == 'solid':
                self.windings[i].conductor_radius = conductor_radii[i]
                self.windings[i].a_cell = np.pi * self.windings[
                    i].conductor_radius ** 2  # Cross section of the solid conductor

            if self.windings[i].conductor_type == 'litz':
                if litz_para_type[i] == 'implicit_ff':
                    self.update_litz_configuration(num=i,
                                                   litz_parametrization_type=litz_para_type[i],
                                                   conductor_radius=conductor_radii[i],
                                                   n_strands=strands_numbers[i],
                                                   strand_radius=strand_radii[i])
                if litz_para_type[i] == 'implicit_litz_radius':
                    self.update_litz_configuration(num=i,
                                                   litz_parametrization_type=litz_para_type[i],
                                                   ff=ff[i],
                                                   n_strands=strands_numbers[i],
                                                   strand_radius=strand_radii[i])
                if litz_para_type[i] == 'implicit_strands_number':
                    self.update_litz_configuration(num=i,
                                                   litz_parametrization_type=litz_para_type[i],
                                                   ff=ff[i],
                                                   conductor_radius=conductor_radii[i],
                                                   strand_radius=strand_radii[i])

            if self.windings[i].conductor_type == ('stacked' or 'full' or 'foil'):
                self.windings[i].a_cell = 1  # TODO: Surface size needed?
                self.windings[i].conductor_radius = 1  # revisit
                # Surface of the litz approximated hexagonal cell
                # self.a_cell = np.pi * self.conductor_radius**2  # * self.ff

    def update_litz_configuration(self, num=0, litz_parametrization_type='implicit_ff', strand_radius=None, ff=None,
                                  conductor_radius=None, n_strands=None):
        """
        - updates the conductor #num (simple transformer: num=0 -> primary winding,
                                                                     num=1 -> secondary winding)
        - used to change litz configurations
        - also used at first initialisation of the geometry
        - needed to always make sure that the relation between litz parameters (strand radius, fill factor, number of
          layers/strands and conductor/litz radius) is valid and consistent
        - 4 parameters, 1 degree of freedom (dof)
        - all parameters are list parameters!

        :param num: internal counter for primary/secondary winding. Do not change!
        :param litz_parametrization_type:
        :param strand_radius: radius of a single strand in [m]
        :type strand_radius: float
        :param ff: in 0....1
        :type ff: float
        :param conductor_radius: radius of conductor in [m], in a list
        :type conductor_radius: list[float]
        :param n_strands: number of strands for one conductor in a list
        :type: n_strands: list[float]

        :return:

        """
        # Choose one of three different parametrization types
        if litz_parametrization_type == 'implicit_ff':
            ff_exact = n_strands * strand_radius ** 2 / conductor_radius ** 2
            print(f"Exact fill factor: {ff_exact}")
            ff = np.around(ff_exact, decimals=2)  # TODO: Interpolation instead of rounding

        if litz_parametrization_type == 'implicit_litz_radius':
            conductor_radius = np.sqrt(n_strands * strand_radius ** 2 / ff)

        if litz_parametrization_type == 'implicit_strands_number':
            n_strands = conductor_radius ** 2 / strand_radius ** 2 * ff

        # Save parameters in Winding objects
        self.windings[num].conductor_radius = conductor_radius
        self.windings[num].n_layers = NbrLayers(n_strands)
        self.windings[num].n_strands = n_strands
        self.windings[num].strand_radius = strand_radius
        self.windings[num].ff = ff
        self.windings[num].a_cell = self.windings[num].n_strands * self.windings[num].strand_radius ** 2 * np.pi \
                                    / self.windings[num].ff

        # Print updated Litz Data
        print(f"Updated Litz Configuration: \n"
              f" ff: {self.windings[num].ff} \n"
              f" Number of layers/strands: {self.windings[num].n_layers}/{self.windings[num].n_strands} \n"
              f" Strand radius: {self.windings[num].strand_radius} \n"
              f" Conductor radius: {self.windings[num].conductor_radius}\n"
              f"---")

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   -  -  -  -  -  -  -  -  -  -  -
    # Full Geometry
    class EIaxi:
        """
        - creates all points needed for the radial axi-symmetric EI core typology

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
            for i in range(0, self.component.n_windings):
                self.p_conductor.insert(i, [])

            self.r_inner = None
            self.r_outer = None

        def draw_outer(self):
            """
            Draws the

            :return:

            """
            # Outer Core
            # (A_zyl=2pi*r*h => h=0.5r=0.25core_w <=> ensure A_zyl=A_core on the tiniest point)
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
                    # The center points are passed by air_gaps.update() and at this point transformed each into 4
                    # corner points
                    self.p_air_gaps[i * 4 + 0] = [-self.component.core.core_w / 2,
                                                  self.component.air_gaps.midpoints[i][1] -
                                                  self.component.air_gaps.midpoints[i][2] / 2,
                                                  0,
                                                  self.component.air_gaps.midpoints[i][3]]

                    self.p_air_gaps[i * 4 + 1] = [self.component.core.core_w / 2,
                                                  self.component.air_gaps.midpoints[i][1] -
                                                  self.component.air_gaps.midpoints[i][2] / 2,
                                                  0,
                                                  self.component.air_gaps.midpoints[i][3]]

                    self.p_air_gaps[i * 4 + 2] = [-self.component.core.core_w / 2,
                                                  self.component.air_gaps.midpoints[i][1] +
                                                  self.component.air_gaps.midpoints[i][2] / 2,
                                                  0,
                                                  self.component.air_gaps.midpoints[i][3]]

                    self.p_air_gaps[i * 4 + 3] = [self.component.core.core_w / 2,
                                                  self.component.air_gaps.midpoints[i][1] +
                                                  self.component.air_gaps.midpoints[i][2] / 2,
                                                  0,
                                                  self.component.air_gaps.midpoints[i][3]]

        def draw_virtual_winding_windows(self):
            # Virtual Windows
            # TODO: make this part of the class VWW...
            #  self.component.vw_type must be an input or so to that class
            separation_hor = 0  # self.component.core.window_h * 0.5
            separation_vert = self.component.core.window_w * 0.5

            if not self.component.component_type == "integrated_transformer":
                # Some examples for virtual windows
                # Concentrated windings

                if self.component.vw_type == "full_window":
                    """
                    The winding window is completely used by one VWW.
                    In case of a transformer, an interleaved winding scheme must be used.
                    """
                    # top window
                    minimum = -self.component.core.window_h / 2 + self.component.isolation.core_cond[0] / 2  # bottom
                    maximum = self.component.core.window_h / 2 - self.component.isolation.core_cond[0]  # top
                    left = self.component.core.core_w / 2 + self.component.isolation.core_cond[0]
                    right = self.r_inner - self.component.isolation.core_cond[0]

                    # Sum the windows up in a list
                    # self.virtual_windows = [[min, max, left, right]]
                    self.component.virtual_winding_windows[0].bot_bound = minimum
                    self.component.virtual_winding_windows[0].top_bound = maximum
                    self.component.virtual_winding_windows[0].left_bound = left
                    self.component.virtual_winding_windows[0].right_bound = right

                if self.component.vw_type == "center":
                    """
                    The winding window is split into two VWWs.
                    The primary winding is placed in the upper half,
                    the secondary winding is placed in the lower half of the winding window
                    """
                    separation_hor = 0

                    # top window
                    min21 = -separation_hor + self.component.isolation.cond_cond[-1] / 2  # separation_hor
                    max21 = self.component.core.window_h / 2 - self.component.isolation.core_cond[0]  # top
                    left21 = self.component.core.core_w / 2 + self.component.isolation.core_cond[0]
                    right21 = self.r_inner - self.component.isolation.core_cond[0]

                    # bottom window
                    min11 = -self.component.core.window_h / 2 + self.component.isolation.core_cond[0] / 2  # bottom
                    max11 = -separation_hor - self.component.isolation.cond_cond[-1] / 2  # separation_hor
                    left11 = self.component.core.core_w / 2 + self.component.isolation.core_cond[0]
                    right11 = self.r_inner - self.component.isolation.core_cond[0]

                    # Sum the windows up in a list
                    virtual_windows = [[min11, max11, left11, right11],
                                       [min21, max21, left21, right21]]
                    for vww in range(0, len(virtual_windows)):
                        self.component.virtual_winding_windows[vww].bot_bound = virtual_windows[vww][0]
                        self.component.virtual_winding_windows[vww].top_bound = virtual_windows[vww][1]
                        self.component.virtual_winding_windows[vww].left_bound = virtual_windows[vww][2]
                        self.component.virtual_winding_windows[vww].right_bound = virtual_windows[vww][3]

                if self.component.vw_type == "something_else":
                    # bottom left window
                    min11 = -self.component.core.window_h / 2 + self.component.isolation.core_cond[0] / 2  # bottom
                    max11 = -separation_hor - self.component.isolation.cond_cond[-1] / 2  # separation_hor
                    left11 = self.component.core.core_w / 2 + self.component.isolation.core_cond[0]
                    right11 = self.r_inner - self.component.isolation.cond_cond[0] - separation_vert

                    # bottom right window
                    min12 = -self.component.core.window_h / 2 + self.component.isolation.core_cond[0] / 2  # bottom
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
            if self.component.component_type == "integrated_transformer":
                """
                If dedicated stray path is the chosen typology the are two winding windows
                These can either be split up into more virtual windows or (in case of bifilar windings) not
        
                """
                # TODO: Separation in more Virtual Winding Windows

                # bot window
                island_right_tmp = inner_points(self.p_window[4], self.p_window[6], self.p_air_gaps)
                min11 = -self.component.core.window_h / 2 + self.component.isolation.core_cond[1]  # bottom
                max11 = island_right_tmp[(self.component.stray_path.start_index - 1) * 2][1] - \
                        self.component.isolation.core_cond[0]  # sep_hor
                left11 = self.component.core.core_w / 2 + self.component.isolation.core_cond[0]
                right11 = self.r_inner - self.component.isolation.core_cond[1]

                # top window
                min21 = island_right_tmp[(self.component.stray_path.start_index - 1) * 2 + 1][1] + \
                        self.component.isolation.core_cond[0] * 1.5 # TODO:remove 1.5 # sep_hor
                max21 = self.component.core.window_h / 2 - self.component.isolation.core_cond[1]  # top
                left21 = self.component.core.core_w / 2 + self.component.isolation.core_cond[0]
                right21 = self.r_inner - self.component.isolation.core_cond[1]

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

                if self.component.virtual_winding_windows[n_win].winding == "interleaved":

                    if self.component.virtual_winding_windows[n_win].scheme == "bifilar":
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

                    if self.component.virtual_winding_windows[n_win].scheme == "vertical":
                        """
                        - Vertical interleaving means a winding scheme where the two conductors are alternating 
                           in vertical (y-)direction
                        - This is practically uncommon
                        - If the turns ratio is != 1, the scheme always begins with the "higher-turns-number's" 
                           conductor
                        """

                    if self.component.virtual_winding_windows[n_win].scheme == "horizontal":
                        """
                        - Horizontal interleaving means a winding scheme where the two conductors are alternating in 
                        horizontal  (x-)direction (Tonnenwicklung)
                        - This is practically most common
                        - If the turns ratio is != 1, the scheme always begins with the "higher-turns-number's" 
                          conductor
                          
                        """

                        # assume 2 winding transformer and dedicated stray path:
                        if self.component.component_type == "integrated_transformer" or self.component.n_windings == 2:

                            # top window
                            if n_win == 0:
                                # Initialize the list, that counts the already placed conductors
                                N_completed = [0, 0]

                                # Initialize the starting conductor
                                if self.component.windings[0].turns[n_win] >= self.component.windings[1].turns[n_win]:
                                    col_cond = 0
                                else:
                                    col_cond = 1

                                # Initialize the x and y coordinate
                                x = left_bound + self.component.windings[col_cond].conductor_radius
                                y = top_bound - self.component.windings[col_cond].conductor_radius

                                # Continue placing as long as not all conductors have been placed
                                while (self.component.windings[0].turns[n_win] - N_completed[0] != 0) or \
                                        (self.component.windings[1].turns[n_win] - N_completed[1] != 0):
                                    if self.component.windings[col_cond].turns[n_win] - N_completed[col_cond] != 0:
                                        # is this winding not already finished?
                                        if x < right_bound - self.component.windings[col_cond].conductor_radius:
                                            while y > bot_bound + self.component.windings[col_cond].conductor_radius and \
                                                    N_completed[col_cond] < self.component.windings[col_cond].turns[n_win]:
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

                                                y -= self.component.windings[col_cond].conductor_radius * 2 + \
                                                     self.component.isolation.cond_cond[col_cond]  # one from bot to top

                                            x += self.component.windings[col_cond].conductor_radius + \
                                                 self.component.windings[(col_cond + 1) % 2].conductor_radius + \
                                                 self.component.isolation.cond_cond[2]  # from left to right

                                            # Reset y
                                            col_cond = (col_cond + 1) % 2
                                            y = top_bound - self.component.windings[col_cond].conductor_radius

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

                            if self.component.virtual_winding_windows[n_win].scheme[num] == "square":

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

                            if self.component.virtual_winding_windows[n_win].scheme[num] == "hexa":

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
                    if self.component.virtual_winding_windows[n_win].winding == "primary":
                        num = 0
                    if self.component.virtual_winding_windows[n_win].winding == "secondary":
                        num = 1

                    if self.component.windings[num].conductor_type == "full":
                        if sum(self.component.windings[num].turns) != 1:
                            print(f"For a \"full\" conductor you must choose 1 turn for each conductor!")
                        # full window conductor
                        self.p_conductor[num].append([left_bound, bot_bound, 0, self.component.mesh.c_conductor[num]])
                        self.p_conductor[num].append([right_bound, bot_bound, 0, self.component.mesh.c_conductor[num]])
                        self.p_conductor[num].append([left_bound, top_bound, 0, self.component.mesh.c_conductor[num]])
                        self.p_conductor[num].append([right_bound, top_bound, 0, self.component.mesh.c_conductor[num]])

                    if self.component.windings[num].conductor_type == "stacked":
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

                    if self.component.windings[num].conductor_type == "foil":
                        # Wrap defined number of turns and chosen thickness
                        if self.component.wrap_para[num] == "fixed_thickness":
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
                        if self.component.wrap_para[num] == "interpolate":
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
                    if self.component.windings[num].conductor_type == "litz" or \
                            self.component.windings[num].conductor_type == "solid":

                        if self.component.virtual_winding_windows[num].scheme == "square":
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

                        if self.component.virtual_winding_windows[num].scheme == "hexa":
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
                if self.component.windings[num].conductor_type == "solid" or \
                        self.component.windings[num].conductor_type == "litz":
                    if int(self.p_conductor[num].shape[0] / 5) < sum(self.component.windings[num].turns):
                        # Warning: warnings.warn("Too many turns that do not fit in the winding window.")
                        # Correct: self.component.windings[num].turns = int(self.p_conductor[num].shape[0]/5)
                        # TODO: break, but remove warning. valid bit should be set to False
                        #  Code must go to the next parameter-iteration step for geometric sweep
                        self.component.valid = False

            # Region for Boundary Condition
            self.p_region_bound[0][:] = [-self.r_outer * self.component.padding,
                                         -(self.component.core.window_h / 2 + self.component.core.core_w / 4)
                                         * self.component.padding,
                                         0,
                                         self.component.mesh.c_core * self.component.padding]
            self.p_region_bound[1][:] = [self.r_outer * self.component.padding,
                                         -(self.component.core.window_h / 2 + self.component.core.core_w / 4)
                                         * self.component.padding,
                                         0,
                                         self.component.mesh.c_core * self.component.padding]
            self.p_region_bound[2][:] = [-self.r_outer * self.component.padding,
                                         (self.component.core.window_h / 2 + self.component.core.core_w / 4)
                                         * self.component.padding,
                                         0,
                                         self.component.mesh.c_core * self.component.padding]
            self.p_region_bound[3][:] = [self.r_outer * self.component.padding,
                                         (self.component.core.window_h / 2 + self.component.core.core_w / 4)
                                         * self.component.padding,
                                         0,
                                         self.component.mesh.c_core * self.component.padding]

        def update(self):

            # Preallocate the arrays, in which the geometries' point coordinates will be stored
            self.p_outer = np.zeros((4, 4))
            self.p_region_bound = np.zeros((4, 4))
            self.p_window = np.zeros((4 * self.component.n_windows, 4))
            self.p_air_gaps = np.zeros((4 * self.component.air_gaps.number, 4))

            # Fitting the outer radius to ensure surface area
            self.r_inner = self.component.core.window_w + self.component.core.core_w / 2
            self.r_outer = np.sqrt((self.component.core.core_w / 2) ** 2 + self.r_inner ** 2)
            # np.sqrt(window_w**2 + window_w * core_w + core_w**2/2)

            #
            self.draw_outer()

            self.draw_window()

            self.draw_air_gaps()

            self.draw_virtual_winding_windows()

            self.draw_conductors()

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

        def calculate_air_gap_lengths_idealized(self, reluctances, types):
            """

            :param reluctances:
            :param types:

            :return:

            """
            air_gap_lengths = [None] * len(reluctances)

            # Go through reluctances
            for n_reluctance, R_0 in enumerate(reluctances):
                if self.component.core.type == "EI" and self.component.dimensionality == "2D axi":

                    if types[n_reluctance] == ("round-round" or "round-inf"):
                        A_core = (self.component.core.core_w / 2) ** 2 * np.pi
                        length = R_0 * self.component.mu0 * A_core

                    if types[n_reluctance] == "cyl-cyl":
                        # return R_0 * self.component.mu0 * w * np.pi * (r_o + r_i)
                        length = None  # TODO: Find explicit formula of idealized air gap length
                air_gap_lengths[n_reluctance] = length

            return air_gap_lengths

        def calculate_air_gap_lengths_with_sct(self, reluctances):
            """
            Method calculates air gap lengths according to the given reluctances.
            Method uses several instance variables of the Reluctance Model.

            .. todo:: List with lists of air gap lengths important to always keep the order of elements [or use a dictionary] future use case integrated_transformer:  [[l_top_1, l_top_2, ...],
                                                      [l_bot_1, l_bot_2, ...],
                                                      [l_stray_1, l_stray_2, ...]]
                                                      use case now to be implemented:  [[l_top_1], [l_bot_1], [l_stray_1]]

            :param reluctances:

            :return: Dictionary with air gap names and the associated lengths

            """
            air_gap_lengths = {}
            b_peaks = {}

            # Go through reluctances
            for n_reluctance, R_0 in enumerate(reluctances):

                # Define the Reluctance function to be solved for the air gap length
                if self.component.core.type == "EI" and self.component.dimensionality == "2D axi":

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
                                length_vec = np.linspace(1e-6, self.max_length[n_reluctance]/4)
                                R_res = np.array([r_sct(length)+R_0 for length in length_vec])
                                R_res_ideal = np.array([r_sct_ideal(length)+R_0 for length in length_vec])
                                R_0_vec = np.array([R_0 for length in length_vec])
                                length_vec = length_vec * 1000

                                plt.figure(figsize=(5, 3))
                                plt.plot(length_vec, R_res*1e-6, label=r"$R_{\mathrm{fringing}}$")
                                plt.plot(length_vec, R_res_ideal*1e-6, label=r"$R_{\mathrm{ideal}}$")
                                plt.plot(length_vec, R_0_vec*1e-6, label=r"$R_{\mathrm{goal}}$")
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
            length_corner = self.A_core/self.component.core.core_w/np.pi
            ri = self.component.core.core_w/2
            ro = ri + self.component.core.window_w

            # Top Part
            b_top = Phi_top / self.A_core
            Vol_top = self.A_core * (100-self.component.stray_path.midpoint) / 100 * self.component.core.window_h + \
                      self.A_core * ((100-self.component.stray_path.midpoint) / 100 * self.component.core.window_h - \
                      self.air_gap_lengths["R_top"])

            p_top = 0.5 * self.component.mu0 * f_N95_mu_imag(self.f_1st, b_top) * Vol_top * 2 * np.pi * self.f_1st * \
                    (b_top / self.component.core.re_mu_rel / self.component.mu0)**2 + \
                    self.p_loss_cyl(Phi_top, length_corner, ri, ro)[0]

            # Bot Part
            b_bot = Phi_bot / self.A_core
            Vol_bot = self.A_core * self.component.stray_path.midpoint / 100 * self.component.core.window_h + \
                      self.A_core * (self.component.stray_path.midpoint / 100 * self.component.core.window_h + \
                      self.air_gap_lengths["R_bot"])
            p_bot = 0.5 * self.component.mu0 * f_N95_mu_imag(self.f_1st, b_bot) * Vol_bot * 2 * np.pi * self.f_1st * \
                    (b_bot / self.component.core.re_mu_rel / self.component.mu0)**2 + \
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
                       self.component.mu0 * f_N95_mu_imag(self.f_1st, b(r, Phi, w)) *\
                       (b(r, Phi, w) / self.component.core.re_mu_rel / self.component.mu0)**2

            return quad(p_loss_density, ri, ro, args=(Phi, w), epsabs=1e-4)

        def max_phi_from_time_phi(self, Phi_top, Phi_bot, Phi_stray):
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

        def visualize_current_and_flux(self, t, Phi_top, Phi_bot, Phi_stray, I1, I2):
            figure, axis = plt.subplots(2, figsize=(4, 4))

            t = np.array(t) * 10**6 / 200000 / 2 / np.pi

            axis[0].plot(t, 1000*np.array(Phi_top), label=r"$\mathcal{\phi}_{\mathrm{top}}$")
            axis[0].plot(t, 1000*np.array(Phi_bot), label=r"$\mathcal{\phi}_{\mathrm{bot}}$")
            axis[0].plot(t, 1000*np.array(Phi_stray), label=r"$\mathcal{\phi}_{\mathrm{stray}}$")

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

            plt.savefig(f"C:/Users/tillp/sciebo/Exchange Till/04_Documentation/Reluctance_Model_Current_Shapes/{I1[0]}.pdf",
                        bbox_inches="tight")

            # plt.show()

        def stray_path_parametrization_ei_axi(self):
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
                    b_stray_mean = (np.abs(Phi_top_max_peak) + np.abs(Phi_bot_max_peak)) / 2 / np.pi / (self.component.core.core_w / 2)**2
                    self.component.stray_path.width = Phi_stray_max_peak / b_stray_mean / (np.pi * self.component.core.core_w)
                if self.stray_path_parametrization == "max_flux":
                    phi_abs = np.array([np.abs(Phi_top_max_peak), np.abs(Phi_bot_max_peak)])
                    b_stray_max = phi_abs.max(0) / np.pi / (self.component.core.core_w / 2)**2
                    self.component.stray_path.width = Phi_stray_max_peak / b_stray_max / self.b_stray_rel_overshoot\
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
            if self.component.core.type == "EI" and self.component.dimensionality == "2D axi":

                # Check winding matrix

                # Inductor
                if self.N.shape == (1,) and self.component.component_type == "inductor":
                    raise NotImplemented

                # Transformer
                if self.N.shape == (2,) and self.component.component_type == "transformer":
                    raise NotImplemented

                # Dedicated stray path
                if self.N.shape == (2, 2) and self.component.component_type == "integrated_transformer":

                    # Calculate goal reluctance matrix
                    R_matrix = calculate_reluctances(N=self.N, L=self.L_goal)

                    # R_goal = [R_top, R_bot, R_stray]
                    R_goal = [R_matrix[0, 0] + R_matrix[0, 1], R_matrix[1, 1] + R_matrix[0, 1], -R_matrix[0, 1]]

                    # Stray path specific parameters
                    # print(R_matrix)
                    # print(R_goal)
                    # print(self.L_goal)
                    self.stray_path_parametrization_ei_axi()

                    # Check for negative Reluctances
                    if all(R >= 0 for R in R_goal) and self.singularity == False:

                        # Calculate the air gap lengths with the help of SCT
                        # air_gap_lengths is a dictionary with air gap names and the associated length
                        self.air_gap_lengths, self.b_peaks = self.calculate_air_gap_lengths_with_sct(reluctances=R_goal)

                        # print(air_gap_lengths.values())

                        # Check for invalid data
                        if "saturated" in self.air_gap_lengths.values() or \
                                "out of bounds" in self.air_gap_lengths.values():
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

            if not os.path.isdir(self.component.path + "/" + "Reluctance_Model"):
                os.mkdir(self.component.path + "/" + "Reluctance_Model")

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

            # Characteristic lengths [for mesh sizes]
            self.skin_mesh_factor = None
            self.c_core = None
            self.c_window = None
            self.c_conductor = [None] * self.component.n_windings
            self.c_center_conductor = [None] * self.component.n_windings  # used for the mesh accuracy in the conductors
            # self.c_air_gap = None  # TODO: make usable again ?

        def generate_mesh(self, refine=0, alternative_error=0):
            """
            - interaction with gmsh
            - mesh generation
                - Skin depth based forward meshing
                - adaptive refinement
                    - with the help of mesh-size-fields/background meshes
                    - with an appropriate local error metric

            :return:

            """
            # Initialization
            gmsh.initialize()

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

            gmsh.option.setNumber("General.Terminal", 1)
            gmsh.model.add(self.component.path_mesh + "geometry")
            # ------------------------------------------ Geometry -------------------------------------------
            # Core generation
            if self.component.core.type == "EI":
                # --------------------------------------- Points --------------------------------------------
                if self.component.y_symmetric == 1:
                    if self.component.dimensionality == "2D axi":

                        # Find points of air gaps (used later)
                        if self.component.air_gaps.number > 0:
                            # Top and bottom point
                            center_right = min_max_inner_points(self.component.ei_axi.p_window[4],
                                                                self.component.ei_axi.p_window[6],
                                                                self.component.ei_axi.p_air_gaps)
                            island_right = inner_points(self.component.ei_axi.p_window[4],
                                                        self.component.ei_axi.p_window[6],
                                                        self.component.ei_axi.p_air_gaps)

                            # Dedicated stray path:
                            if self.component.component_type == "integrated_transformer":
                                # mshopt stray_path_gap = [[], []]
                                # mshopt stray_path_gap[0][:] = island_right[(self.stray_path.start_index-1)*2][:]
                                # mshopt stray_path_gap[1][:] = island_right[(self.stray_path.start_index-1)*2+1][:]
                                island_right[(self.component.stray_path.start_index - 1) * 2][0] = \
                                    self.component.stray_path.radius
                                island_right[(self.component.stray_path.start_index - 1) * 2 + 1][0] = \
                                    self.component.stray_path.radius

                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        # Pre-Definitions
                        # Points
                        p_core = []
                        p_island = []
                        p_cond = [[], []]
                        p_region = []
                        # Curves
                        l_bound_core = []
                        l_bound_air = []
                        l_core_air = []
                        l_cond = [[], []]
                        l_region = []
                        curve_loop_cond = [[], []]
                        # Curve Loops
                        curve_loop_island = []
                        curve_loop_air = []
                        # curve_loop_outer_air = []
                        # curve_loop_bound = []
                        # Plane Surfaces
                        plane_surface_core = []
                        plane_surface_cond = [[], []]
                        plane_surface_air = []
                        plane_surface_outer_air = []

                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        # Geometry definitions: points -> lines -> curve loops -> surfaces
                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        # Main Core

                        # Points
                        # (index refers to sketch)

                        # First point (left point of lowest air gap)
                        if self.component.air_gaps.number > 0:
                            p_core.append(gmsh.model.geo.addPoint(0,
                                                                  center_right[0][1],
                                                                  center_right[0][2],
                                                                  self.c_core))
                        if self.component.air_gaps.number == 0:
                            p_core.append(None)  # dummy filled for no air gap special case

                        # Go down and counter-clockwise
                        # Four points around the core
                        p_core.append(gmsh.model.geo.addPoint(0,
                                                              self.component.ei_axi.p_outer[1][1],
                                                              self.component.ei_axi.p_outer[1][2],
                                                              self.component.ei_axi.p_outer[1][3]))

                        p_core.append(gmsh.model.geo.addPoint(self.component.ei_axi.p_outer[1][0],
                                                              self.component.ei_axi.p_outer[1][1],
                                                              self.component.ei_axi.p_outer[1][2],
                                                              self.component.ei_axi.p_outer[1][3]))

                        p_core.append(gmsh.model.geo.addPoint(self.component.ei_axi.p_outer[3][0],
                                                              self.component.ei_axi.p_outer[3][1],
                                                              self.component.ei_axi.p_outer[3][2],
                                                              self.component.ei_axi.p_outer[3][3]))

                        p_core.append(gmsh.model.geo.addPoint(0,
                                                              self.component.ei_axi.p_outer[3][1],
                                                              self.component.ei_axi.p_outer[3][2],
                                                              self.component.ei_axi.p_outer[3][3]))

                        # Two points of highest air gap
                        if self.component.air_gaps.number > 0:
                            p_core.append(gmsh.model.geo.addPoint(0,
                                                                  center_right[1][1],
                                                                  center_right[1][2],
                                                                  self.c_core))

                            p_core.append(gmsh.model.geo.addPoint(center_right[1][0],
                                                                  center_right[1][1],
                                                                  center_right[1][2],
                                                                  self.c_window))

                        if self.component.air_gaps.number == 0:
                            p_core.append(None)  # dummy filled for no air gap special case
                            p_core.append(None)  # dummy filled for no air gap special case

                        # Clockwise
                        # Four points of the window
                        p_core.append(gmsh.model.geo.addPoint(self.component.ei_axi.p_window[6][0],
                                                              self.component.ei_axi.p_window[6][1],
                                                              self.component.ei_axi.p_window[6][2],
                                                              self.component.ei_axi.p_window[6][3]))

                        p_core.append(gmsh.model.geo.addPoint(self.component.ei_axi.p_window[7][0],
                                                              self.component.ei_axi.p_window[7][1],
                                                              self.component.ei_axi.p_window[7][2],
                                                              self.component.ei_axi.p_window[7][3]))

                        p_core.append(gmsh.model.geo.addPoint(self.component.ei_axi.p_window[5][0],
                                                              self.component.ei_axi.p_window[5][1],
                                                              self.component.ei_axi.p_window[5][2],
                                                              self.component.ei_axi.p_window[5][3]))

                        p_core.append(gmsh.model.geo.addPoint(self.component.ei_axi.p_window[4][0],
                                                              self.component.ei_axi.p_window[4][1],
                                                              self.component.ei_axi.p_window[4][2],
                                                              self.component.ei_axi.p_window[4][3]))

                        # Last point of lowest air gap
                        if self.component.air_gaps.number > 0:
                            p_core.append(gmsh.model.geo.addPoint(center_right[0][0],
                                                                  center_right[0][1],
                                                                  center_right[0][2],
                                                                  self.c_window))

                        if self.component.air_gaps.number == 0:
                            p_core.append(None)  # dummy filled for no air gap special case

                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        # Curves
                        # (index refers to sketch)
                        # To be added: Case Air Gaps directly on outer leg

                        # Curves: Boundary - Core
                        if self.component.air_gaps.number > 0:
                            l_bound_core.append(gmsh.model.geo.addLine(p_core[0],
                                                                       p_core[1]))
                        l_bound_core.append(gmsh.model.geo.addLine(p_core[1],
                                                                   p_core[2]))
                        l_bound_core.append(gmsh.model.geo.addLine(p_core[2],
                                                                   p_core[3]))
                        l_bound_core.append(gmsh.model.geo.addLine(p_core[3],
                                                                   p_core[4]))
                        if self.component.air_gaps.number > 0:
                            l_bound_core.append(gmsh.model.geo.addLine(p_core[4],
                                                                       p_core[5]))
                        if self.component.air_gaps.number == 0:
                            l_bound_core.append(gmsh.model.geo.addLine(p_core[4],
                                                                       p_core[1]))
                        # Curves: Core - Air
                        if self.component.air_gaps.number > 0:
                            l_core_air.append(gmsh.model.geo.addLine(p_core[5],
                                                                     p_core[6]))
                            l_core_air.append(gmsh.model.geo.addLine(p_core[6],
                                                                     p_core[7]))
                        l_core_air.append(gmsh.model.geo.addLine(p_core[7],
                                                                 p_core[8]))
                        l_core_air.append(gmsh.model.geo.addLine(p_core[8],
                                                                 p_core[9]))
                        l_core_air.append(gmsh.model.geo.addLine(p_core[9],
                                                                 p_core[10]))

                        if self.component.air_gaps.number > 0:
                            l_core_air.append(gmsh.model.geo.addLine(p_core[10],
                                                                     p_core[11]))
                            l_core_air.append(gmsh.model.geo.addLine(p_core[11],
                                                                     p_core[0]))
                        if self.component.air_gaps.number == 0:
                            l_core_air.append(gmsh.model.geo.addLine(p_core[10],
                                                                     p_core[7]))
                        # Plane: Main Core --> plane_surface_core[0]
                        if self.component.air_gaps.number > 0:
                            curve_loop_core = gmsh.model.geo.addCurveLoop(l_bound_core + l_core_air)
                            plane_surface_core.append(gmsh.model.geo.addPlaneSurface([curve_loop_core]))
                        if self.component.air_gaps.number == 0:
                            curve_loop_bound_core = gmsh.model.geo.addCurveLoop(l_bound_core)
                            curve_loop_core_air = gmsh.model.geo.addCurveLoop(l_core_air)
                            plane_surface_core.append(
                                gmsh.model.geo.addPlaneSurface([curve_loop_bound_core, curve_loop_core_air]))

                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        # - Core parts between Air Gaps
                        # Points of Core Islands (index refers to sketch)
                        if self.component.air_gaps.number != 0:
                            while island_right.shape[0] > 0:
                                # take two points with lowest y-coordinates
                                min_island_right = np.argmin(island_right[:, 1])
                                p_island.append(gmsh.model.geo.addPoint(0,
                                                                        island_right[min_island_right, 1],
                                                                        island_right[min_island_right, 2],
                                                                        self.c_core))
                                p_island.append(gmsh.model.geo.addPoint(island_right[min_island_right, 0],
                                                                        island_right[min_island_right, 1],
                                                                        island_right[min_island_right, 2],
                                                                        self.c_window))

                                island_right = np.delete(island_right, min_island_right, 0)

                            # Curves of Core Islands (index refers to sketch)
                            for i in range(0, int(len(p_island) / 4)):
                                l_core_air.append(gmsh.model.geo.addLine(p_island[4 * i + 0], p_island[4 * i + 1]))
                                l_core_air.append(gmsh.model.geo.addLine(p_island[4 * i + 1], p_island[4 * i + 3]))
                                l_core_air.append(gmsh.model.geo.addLine(p_island[4 * i + 3], p_island[4 * i + 2]))
                                l_bound_core.append(gmsh.model.geo.addLine(p_island[4 * i + 2], p_island[4 * i + 0]))
                                # Iterative plane creation
                                curve_loop_island.append(gmsh.model.geo.addCurveLoop(
                                    [l_core_air[-3], l_core_air[-2], l_core_air[-1], l_bound_core[-1]]))
                                plane_surface_core.append(gmsh.model.geo.addPlaneSurface([curve_loop_island[-1]]))

                        # Curves: Boundary - Air
                        if self.component.air_gaps.number == 1:
                            l_bound_air.append(gmsh.model.geo.addLine(p_core[0], p_core[5]))
                        else:
                            for i in range(0, int(len(p_island) / 4)):
                                if i == 0:  # First Line
                                    l_bound_air.append(gmsh.model.geo.addLine(p_core[0], p_island[0]))
                                else:  # Middle Lines
                                    l_bound_air.append(
                                        gmsh.model.geo.addLine(p_island[4 * (i - 1) + 2], p_island[4 * i + 0]))
                                if i == int(len(p_island) / 4) - 1:  # Last Line
                                    l_bound_air.append(gmsh.model.geo.addLine(p_island[-2], p_core[5]))

                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        # Conductors
                        # Points of Conductors
                        for num in range(0, self.component.n_windings):
                            for i in range(0, self.component.ei_axi.p_conductor[num].shape[0]):
                                p_cond[num].append(
                                    gmsh.model.geo.addPoint(
                                        self.component.ei_axi.p_conductor[num][i][0],  # x
                                        self.component.ei_axi.p_conductor[num][i][1],  # y
                                        0,  # z
                                        self.component.ei_axi.p_conductor[num][i][3]))

                            # Curves of Conductors
                            if self.component.windings[num].conductor_type == "litz" or \
                                    self.component.windings[num].conductor_type == "solid":
                                for i in range(0, int(len(p_cond[num]) / 5)):
                                    l_cond[num].append(gmsh.model.geo.addCircleArc(
                                        p_cond[num][5 * i + 1],
                                        p_cond[num][5 * i + 0],
                                        p_cond[num][5 * i + 2]))
                                    l_cond[num].append(gmsh.model.geo.addCircleArc(
                                        p_cond[num][5 * i + 2],
                                        p_cond[num][5 * i + 0],
                                        p_cond[num][5 * i + 3]))
                                    l_cond[num].append(gmsh.model.geo.addCircleArc(
                                        p_cond[num][5 * i + 3],
                                        p_cond[num][5 * i + 0],
                                        p_cond[num][5 * i + 4]))
                                    l_cond[num].append(gmsh.model.geo.addCircleArc(
                                        p_cond[num][5 * i + 4],
                                        p_cond[num][5 * i + 0],
                                        p_cond[num][5 * i + 1]))
                                    # Iterative plane creation
                                    curve_loop_cond[num].append(gmsh.model.geo.addCurveLoop(
                                        [l_cond[num][i * 4 + 0],
                                         l_cond[num][i * 4 + 1],
                                         l_cond[num][i * 4 + 2],
                                         l_cond[num][i * 4 + 3]]))
                                    plane_surface_cond[num].append(
                                        gmsh.model.geo.addPlaneSurface([curve_loop_cond[num][i]]))
                            else:
                                # Rectangle conductor cut
                                for i in range(0, int(len(p_cond[num]) / 4)):
                                    l_cond[num].append(gmsh.model.geo.addLine(p_cond[num][4 * i + 0],
                                                                              p_cond[num][4 * i + 2]))
                                    l_cond[num].append(gmsh.model.geo.addLine(p_cond[num][4 * i + 2],
                                                                              p_cond[num][4 * i + 3]))
                                    l_cond[num].append(gmsh.model.geo.addLine(p_cond[num][4 * i + 3],
                                                                              p_cond[num][4 * i + 1]))
                                    l_cond[num].append(gmsh.model.geo.addLine(p_cond[num][4 * i + 1],
                                                                              p_cond[num][4 * i + 0]))
                                    # Iterative plane creation
                                    curve_loop_cond[num].append(gmsh.model.geo.addCurveLoop([l_cond[num][i * 4 + 0],
                                                                                             l_cond[num][i * 4 + 1],
                                                                                             l_cond[num][i * 4 + 2],
                                                                                             l_cond[num][i * 4 + 3]]))
                                    plane_surface_cond[num].append(
                                        gmsh.model.geo.addPlaneSurface([curve_loop_cond[num][i]]))

                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        # Air
                        # Points are partwise double designated

                        l_air_tmp = l_core_air[:7]
                        for i in range(0, len(l_bound_air)):
                            l_air_tmp.append(l_bound_air[i])
                            if i < len(l_bound_air) - 1:
                                l_air_tmp.append(l_core_air[7 + 3 * i])
                                l_air_tmp.append(l_core_air[7 + 3 * i + 1])
                                l_air_tmp.append(l_core_air[7 + 3 * i + 2])

                        curve_loop_air.append(gmsh.model.geo.addCurveLoop(l_air_tmp))

                        # Need flatten list of all! conductors
                        flatten_curve_loop_cond = [j for sub in curve_loop_cond for j in sub]
                        plane_surface_air.append(
                            gmsh.model.geo.addPlaneSurface(curve_loop_air + flatten_curve_loop_cond))

                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        # Boundary
                        if self.component.region is None:
                            l_bound_tmp = l_bound_core[:5]
                            for i in range(0, len(l_bound_air)):
                                l_bound_tmp.append(l_bound_air[-i - 1])
                                if i != len(l_bound_air) - 1:  # last run
                                    l_bound_tmp.append(l_bound_core[-i - 1])

                        else:
                            # Generate Lines of Region
                            # start top left and go clockwise
                            p_region.append(gmsh.model.geo.addPoint(0,
                                                                    self.component.ei_axi.p_region_bound[2][1],
                                                                    self.component.ei_axi.p_region_bound[2][2],
                                                                    self.component.ei_axi.p_region_bound[2][3]))

                            p_region.append(gmsh.model.geo.addPoint(self.component.ei_axi.p_region_bound[3][0],
                                                                    self.component.ei_axi.p_region_bound[3][1],
                                                                    self.component.ei_axi.p_region_bound[3][2],
                                                                    self.component.ei_axi.p_region_bound[3][3]))

                            p_region.append(gmsh.model.geo.addPoint(self.component.ei_axi.p_region_bound[1][0],
                                                                    self.component.ei_axi.p_region_bound[1][1],
                                                                    self.component.ei_axi.p_region_bound[1][2],
                                                                    self.component.ei_axi.p_region_bound[1][3]))

                            p_region.append(gmsh.model.geo.addPoint(0,
                                                                    self.component.ei_axi.p_region_bound[0][1],
                                                                    self.component.ei_axi.p_region_bound[0][2],
                                                                    self.component.ei_axi.p_region_bound[0][3]))

                            # Outer Region Lines
                            l_region.append(gmsh.model.geo.addLine(p_core[4],
                                                                   p_region[0]))
                            l_region.append(gmsh.model.geo.addLine(p_region[0],
                                                                   p_region[1]))
                            l_region.append(gmsh.model.geo.addLine(p_region[1],
                                                                   p_region[2]))
                            l_region.append(gmsh.model.geo.addLine(p_region[2],
                                                                   p_region[3]))
                            l_region.append(gmsh.model.geo.addLine(p_region[3],
                                                                   p_core[1]))

                            # Boundary Line
                            l_bound_tmp = [l_bound_core[4]]

                            for i in range(0, len(l_region)):
                                l_bound_tmp.append(l_region[i])

                            l_bound_tmp.append(l_bound_core[0])

                            for i in range(0, len(l_bound_air)):
                                l_bound_tmp.append(l_bound_air[-i - 1])
                                if i != len(l_bound_air) - 1:  # last run
                                    l_bound_tmp.append(l_bound_core[-i - 1])

                            # Outer Air Surface
                            curve_loop_outer_air = gmsh.model.geo.addCurveLoop(l_region + l_bound_core[1:4])
                            plane_surface_outer_air.append(gmsh.model.geo.addPlaneSurface([curve_loop_outer_air]))

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if self.component.visualize_before:
                # synchronize must be placed before!!! the color properties are set

                gmsh.model.geo.synchronize()
                # Colors
                for i in range(0, len(plane_surface_core)):
                    gmsh.model.setColor([(2, plane_surface_core[i])], 50, 50, 50)
                gmsh.model.setColor([(2, plane_surface_air[0])], 255, 255, 255)
                for num in range(0, self.component.n_windings):
                    for i in range(0, len(plane_surface_cond[num])):
                        gmsh.model.setColor([(2, plane_surface_cond[num][i])], 200*num, 200*(1-num), 0)

                # Output .msh file
                gmsh.option.setNumber("Mesh.SaveAll", 1)
                # gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)

                # gmsh.model.geo.synchronize()
                # gmsh.write(self.component.path_mesh + "color.msh")

                gmsh.write(
                    self.component.path + "/" + self.component.path_mesh + "color.msh")  # Win10 can handle slash

                gmsh.model.mesh.generate(2)
                gmsh.fltk.run()

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Define physical Surfaces and Curves
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Core
            ps_core = gmsh.model.geo.addPhysicalGroup(2, plane_surface_core, tag=2000)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Conductors
            ps_cond = [[], []]
            for num in range(0, self.component.n_windings):

                if self.component.windings[num].conductor_type == "foil" or \
                        self.component.windings[num].conductor_type == "solid" or \
                        self.component.windings[num].conductor_type == "full" or \
                        self.component.windings[num].conductor_type == "stacked":
                    for i in range(0, sum(self.component.windings[num].turns)):
                        ps_cond[num].append(
                            gmsh.model.geo.addPhysicalGroup(2, [plane_surface_cond[num][i]], tag=4000 + 1000 * num + i))
                if self.component.windings[num].conductor_type == "litz":
                    for i in range(0, sum(self.component.windings[num].turns)):
                        ps_cond[num].append(
                            gmsh.model.geo.addPhysicalGroup(2, [plane_surface_cond[num][i]], tag=6000 + 1000 * num + i))

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Air
            ps_air = gmsh.model.geo.addPhysicalGroup(2, plane_surface_air, tag=1000)
            # ps_air_ext = gmsh.model.geo.addPhysicalGroup(2, plane_surface_outer_air, tag=1001)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Boundary
            pc_bound = gmsh.model.geo.addPhysicalGroup(1, l_bound_tmp, tag=1111)
            # print(f"Physical Conductor Surfaces: {ps_cond}")

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Set names [optional]
            gmsh.model.setPhysicalName(2, ps_core, "CORE")
            for num in range(0, self.component.n_windings):
                for i in range(0, len(ps_cond[num])):
                    gmsh.model.setPhysicalName(2, ps_cond[num][i], f"COND{num + 1}")
            gmsh.model.setPhysicalName(2, ps_air, "AIR")
            gmsh.model.setPhysicalName(1, pc_bound, "BOUND")

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # - Forward Meshing -
            p_inter = None
            # Inter Conductors
            for n_win in range(0, len(self.component.virtual_winding_windows)):
                if self.component.virtual_winding_windows[n_win].winding != "interleaved":
                    for num in range(0, self.component.n_windings):
                        p_inter = []
                        x_inter = []
                        y_inter = []
                        j = 0

                        if self.component.windings[num].conductor_type == "solid" and \
                                self.component.windings[num].turns[n_win] > 1:
                            while self.component.ei_axi.p_conductor[num][5 * j][1] == \
                                    self.component.ei_axi.p_conductor[num][5 * j + 5][1]:
                                x_inter.append(
                                    0.5 * (self.component.ei_axi.p_conductor[num][5 * j][0] +
                                           self.component.ei_axi.p_conductor[num][5 * j + 5][0]))
                                j += 1
                                if j == self.component.windings[num].turns[n_win] - 1:
                                    break
                            j += 1
                            if int(self.component.windings[num].turns[n_win] / j) > 1:
                                for i in range(0, int(self.component.windings[num].turns[n_win] / j)):
                                    if 5 * j * i + 5 * j >= len(self.component.ei_axi.p_conductor[num][:]):
                                        break
                                    y_inter.append(0.5 * (self.component.ei_axi.p_conductor[num][5 * j * i][1] +
                                                          self.component.ei_axi.p_conductor[num][5 * j * i + 5 * j][1]))
                                for x in x_inter:
                                    for y in y_inter:
                                        p_inter.append(gmsh.model.geo.addPoint(x,
                                                                               y,
                                                                               0,
                                                                               self.c_center_conductor[num]))

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

            # Conductor Center
            for num in range(0, self.component.n_windings):
                for i in range(0, int(len(p_cond[num]) / 5)):
                    gmsh.model.mesh.embed(0, [p_cond[num][5 * i + 0]], 2, plane_surface_cond[num][i])

            # Embed points for mesh refinement
            # Inter Conductors
            for n_win in range(0, len(self.component.virtual_winding_windows)):
                if self.component.virtual_winding_windows[n_win].winding != "interleaved":
                    gmsh.model.mesh.embed(0, p_inter, 2, plane_surface_air[0])
            # Stray path
            # mshopt gmsh.model.mesh.embed(0, stray_path_mesh_optimizer, 2, plane_surface_core[2])

            # Synchronize again
            gmsh.model.geo.synchronize()

            # Output .msh file
            # TODO: What are these flags about???
            gmsh.option.setNumber("Mesh.SaveAll", 1)
            gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
            gmsh.option.setNumber("Mesh.SurfaceFaces", 1)

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
                gmsh.write(self.component.path_mesh + "geometry.msh")
            else:
                print("\nMeshing...\n")
                gmsh.model.mesh.generate(2)

            # Mesh direction
            if not os.path.isdir(self.component.path + "/" + self.component.path_mesh):
                os.mkdir(self.component.path + "/" + self.component.path_mesh)

            # Check operating system
            if sys.platform == "linux" or sys.platform == "linux2":
                gmsh.write(self.component.path + "/" + self.component.path_mesh + "geometry.msh")
            elif sys.platform == "darwin":  # OS X
                gmsh.write(self.component.path + "/" + self.component.path_mesh + "geometry.msh")
            elif sys.platform == "win32":
                gmsh.write(
                    self.component.path + "/" + self.component.path_mesh + "geometry.msh")  # Win10 can handle slash

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Open gmsh GUI for visualization
            # gmsh.fltk.run()
            # Terminate gmsh
            gmsh.finalize()

        def mesh(self, frequency=None, skin_mesh_factor=1):
            self.component.high_level_geo_gen(frequency=frequency, skin_mesh_factor=skin_mesh_factor)
            if self.component.valid:
                self.generate_mesh()

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    # GetDP Interaction / Simulation / Excitation
    def excitation(self, f, i, phases=None, ex_type='current', imposed_red_f=0):
        """
        - excitation of the electromagnetic problem
        - current, voltage or current density
        - frequency or reduced frequency

        :param phases:
        :param f:
        :param i:
        :param ex_type:
        :param imposed_red_f:

        :return:

        """
        phases = phases or []

        print(f"\n---\n"
              f"Excitation: \n"
              f"Frequency: {f}\n"
              f"Current(s): {i}\n"
              f"Phase(s): {phases}\n")

        # -- Excitation --
        self.flag_imposed_reduced_frequency = imposed_red_f  # if == 0 --> impose frequency f
        self.flag_excitation_type = ex_type  # 'current', 'current_density', 'voltage'

        phases = np.asarray(phases)
        self.red_freq = np.empty(2)

        for num in range(0, self.n_windings):
            # Imposed current, current density or voltage
            if self.flag_excitation_type == 'current':
                self.current[num] = i[num]
                if len(phases) != 0:
                    self.phase_tmp = phases / 180
            if self.flag_excitation_type == 'current_density':
                raise NotImplementedError
            if self.flag_excitation_type == 'voltage':
                raise NotImplementedError
                # self.voltage = 2

            # -- Frequency --
            self.frequency = f  # in Hz
            if self.flag_imposed_reduced_frequency == 1:
                self.red_freq[num] = 4
            else:
                if self.frequency != 0:
                    self.delta = np.sqrt(2 / (2 * self.frequency * np.pi * self.windings[0].cond_sigma * self.mu0))

                    if self.windings[num].conductor_type == "litz":
                        self.red_freq[num] = self.windings[num].strand_radius / self.delta
                    elif self.windings[num].conductor_type == "solid":
                        self.red_freq[num] = self.windings[num].conductor_radius / self.delta
                    else:
                        print("Wrong???")
                        print(self.windings[num].conductor_type)
                        self.red_freq[num] = 1  # TODO: doesn't make sense like this
                else:
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

        text_file = open(self.path + "/Parameter.pro", "w")

        text_file.write(f"DirRes = \"{self.path_res}\";\n")
        text_file.write(f"DirResFields = \"{self.path_res_fields}\";\n")
        text_file.write(f"DirResVals = \"{self.path_res_values}\";\n")
        text_file.write(f"DirResValsPrimary = \"{self.path_res_values}Primary/\";\n")
        text_file.write(f"DirResValsSecondary = \"{self.path_res_values}Secondary/\";\n")
        text_file.write(f"DirResCirc = \"{self.path_res_circuit}\";\n")

        # Visualisation
        if self.plot_fields == "standard":
            text_file.write(f"Flag_show_standard_fields = 1;\n")
        else:
            text_file.write(f"Flag_show_standard_fields = 0;\n")

        # Magnetic Component Type
        if self.component_type == 'inductor':
            text_file.write(f"Flag_Transformer = 0;\n")
        if self.component_type == 'transformer':
            text_file.write(f"Flag_Transformer = 1;\n")
        if self.component_type == 'integrated_transformer':
            text_file.write(f"Flag_Transformer = 1;\n")

        # Frequency
        text_file.write("Freq = %s;\n" % self.frequency)
        text_file.write(f"delta = {self.delta};\n")

        # Core Loss
        text_file.write(f"Flag_Steinmetz_loss = {self.core.steinmetz_loss};\n")
        text_file.write(f"Flag_Generalized_Steinmetz_loss = {self.core.generalized_steinmetz_loss};\n")

        # TODO: Make following definitions general
        self.core.im_epsilon_rel = f_N95_er_imag(f=self.frequency)
        text_file.write(f"e_r_imag = {self.core.im_epsilon_rel};\n")
        text_file.write(f"mu_r_imag = {self.core.im_mu_rel};\n")

        if self.core.steinmetz_loss:
            text_file.write(f"ki = {self.ki};\n")
            text_file.write(f"alpha = {self.alpha};\n")
            text_file.write(f"beta = {self.beta};\n")
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
            if self.windings[num].conductor_type == 'litz':
                text_file.write(f"Flag_HomogenisedModel{num + 1} = 1;\n")
            else:
                text_file.write(f"Flag_HomogenisedModel{num + 1} = 0;\n")
            text_file.write("Flag_imposedRr = %s;\n" % self.flag_imposed_reduced_frequency)

            # -- Geometry --
            # Number of conductors
            text_file.write(f"NbrCond{num + 1} = {sum(self.windings[num].turns)};\n")

            # For stranded Conductors:
            # text_file.write(f"NbrstrandedCond = {self.turns};\n")  # redundant
            if self.windings[num].conductor_type == "litz":
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
                text_file.write(f"Val_EE_{num + 1} = {self.current[num]};\n")
                text_file.write(f"Phase_{num + 1} = Pi*{self.phase_tmp[num]};\n")
            if self.flag_excitation_type == 'current_density':
                text_file.write(f"Val_EE_{num + 1} = {self.current_density[num]};\n")
            if self.flag_excitation_type == 'voltage':
                text_file.write(f"Val_EE_{num + 1} = {self.voltage[num]};\n")
            print(f"Cell surface area: {self.windings[num].a_cell} \n"
                  f"Reduced frequency: {self.red_freq[num]}")
            if self.red_freq[num] > 1.25 and self.windings[num].conductor_type == "litz":
                # TODO: Allow higher reduced frequencies
                print(f"Litz Coefficients only implemented for X<=1.25")
                raise Warning
            # Reduced Frequency
            text_file.write(f"Rr{num + 1} = {self.red_freq[num]};\n")

        """
        # Coordinates of the rectangular winding window
        if self.dimensionality == "2D axi":
            text_file.write("Xw1 = %s;\n" % self.ei_axi.p_window[4, 0])
            text_file.write("Xw2 = %s;\n" % self.ei_axi.p_window[5, 0])
        else:
            raise NotImplementedError("Only axi-symmetric case implemented :(")
        """
        # -- Materials --

        # Nature Constants
        text_file.write(f"mu0 = 4.e-7 * Pi;\n")
        text_file.write(f"nu0 = 1 / mu0;\n")

        # Material Properties
        # Conductor Material
        text_file.write(f"SigmaCu = {self.windings[0].cond_sigma};\n")

        # Core Material
        # if self.frequency == 0:
        if self.core.non_linear:
            text_file.write(f"Flag_NL = 1;\n")
            text_file.write(f"Core_Material = {self.core.material};\n")
        else:
            text_file.write(f"Flag_NL = 0;\n")
            text_file.write(f"mur = {self.core.re_mu_rel};\n")
        # if self.frequency != 0:
        #    text_file.write(f"Flag_NL = 0;\n")
        #    text_file.write(f"mur = {self.core.re_mu_rel};\n")

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

        # get model file names with correct path
        msh_file = self.onelab_client.getPath(self.path_mesh + 'geometry.msh')
        solver = self.onelab_client.getPath('ind_axi_python_controlled' + '.pro')

        # Run simulations as sub clients (non blocking??)
        mygetdp = self.onelab + 'getdp'
        self.onelab_client.runSubClient('myGetDP', mygetdp + ' ' + solver + ' -msh ' + msh_file + ' -solve Analysis -v2')

    def write_log(self):
        """
        Method reads back the results from the .dat result files created by the ONELAB simulation client and stores
        them in a JSON log file.
        :return:
        :rtype: None
        """

        log_path = self.path + "/" + self.path_res + 'result_log.json'
        file = open(log_path, 'w', encoding='utf-8')

        winding_result_path = ["Primary", "Secondary", "Tertiary"]

        log = {}
        log["Losses"] = {}

        # Write Core Losses
        log["Losses"]["Core_Eddy_Current"] = self.load_result(res_name="CoreEddyCurrentLosses")[0]

        # Write Winding Losses
        total_winding_losses = 0
        for winding in range(0, self.n_windings):
            log["Losses"][f"Winding_{winding + 1}"] = {}
            log["Losses"][f"Winding_{winding + 1}"][f"Turns"] = []
            if self.windings[winding].conductor_type == "litz":
                losses_winding = self.load_result(res_name=f"j2H_{winding+1}")
                for turn in range(0, self.windings[winding].turns[0]):
                    log["Losses"][f"Winding_{winding + 1}"][f"Turns"].append(self.load_result(res_name=winding_result_path[winding]+f"/Losses_turn_{turn + 1}")[0])
            else:
                losses_winding = self.load_result(res_name=f"j2F_{winding+1}")
                for turn in range(0, self.windings[winding].turns[0]):
                    log["Losses"][f"Winding_{winding + 1}"][f"Turns"].append(self.load_result(res_name=winding_result_path[winding]+f"/Losses_turn_{turn + 1}")[0])

            log["Losses"][f"Winding_{winding+1}"]["total_winding"] = losses_winding[0]

            total_winding_losses += losses_winding[0]
        log["Losses"]["all_windings"] = total_winding_losses




        # loss_log["Winding_Total"] = self.load_result(res_name="WindingTotal")



        json.dump(log, file, indent=2, ensure_ascii=False)
        file.close()

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   -  -  -  -  -  -  -  -  -  -  -
    # Thermal simulation

    # Start thermal simulation
    def thermal_simulation(self, thermal_conductivity) -> None:
        """

        Starts the thermal simulation using thermal.py

        :return: -
        """

        # Set necessary paths
        onelab_folder_path = self.onelab
        model_mesh_file_path_old = os.path.join(self.mesh.component.path, self.mesh.component.path_mesh, "geometry.msh")
        results_log_file_path = os.path.join(self.path, self.path_res, "result_log.json")
        model_mesh_file_path = os.path.join(os.path.dirname(th.__file__), "thermal_mesh.msh")

        # Create copy of the current mesh, because it will be changed for the thermal simulation
        shutil.copy(model_mesh_file_path_old, model_mesh_file_path)

        # Set tags
        # TODO The tags need to be set dynamically based by the model
        winding_tags = [list(range(4000, 4036)), list(range(7000, 7011)), None]
        tags = {
            "core_tag": 2000,
            "background_tag": 1000,
            "winding_tags": winding_tags,
            "core_point_tags": [5, 4, 3, 2],  # Order: top left, top right, bottom right, bottom left
        }

        # Mesh size -> Used when creating the case
        # TODO Currently fixed.. Can be changed dynamically?
        mesh_size = 0.001

        # Core area -> Is needed to estimate the heat flux
        # TODO Needs to be calculated dynamically
        core_area = 0.00077
        # Use th_functions.calculate_heat_flux_core(losses["Core_Eddy_Current"], )
        # Set wire radii
        wire_radii = [winding.conductor_radius for winding in self.windings]

        # When a gmsh window should open showing the simulation results
        show_results = True

        th.thermal(onelab_folder_path, model_mesh_file_path, results_log_file_path, tags, thermal_conductivity,
                   mesh_size, core_area, wire_radii, show_results)

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

        gmsh.initialize()
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

        if any(self.windings[i].conductor_type != 'litz' for i in range(0, self.n_windings)):
            # Ohmic losses (weighted effective value of current density)
            gmsh.open(self.path + "/" + self.path_res_fields + "j2F.pos")
            gmsh.option.setNumber(f"View[{view}].ScaleType", 2)
            gmsh.option.setNumber(f"View[{view}].RangeType", 2)
            gmsh.option.setNumber(f"View[{view}].SaturateValues", 1)
            gmsh.option.setNumber(f"View[{view}].CustomMin", gmsh.option.getNumber(f"View[{view}].Min") + epsilon)
            gmsh.option.setNumber(f"View[{view}].CustomMax", gmsh.option.getNumber(f"View[{view}].Max"))
            gmsh.option.setNumber(f"View[{view}].ColormapNumber", 1)
            gmsh.option.setNumber(f"View[{view}].IntervalsType", 2)
            gmsh.option.setNumber(f"View[{view}].NbIso", 40)
            print(gmsh.option.getNumber(f"View[{view}].Max"))
            view += 1

        if any(self.windings[i].conductor_type == 'litz' for i in range(0, self.n_windings)):
            # Ohmic losses (weighted effective value of current density)
            gmsh.open(self.path + "/" + self.path_res_fields + "jH.pos")
            gmsh.option.setNumber(f"View[{view}].ScaleType", 2)
            gmsh.option.setNumber(f"View[{view}].RangeType", 2)
            gmsh.option.setNumber(f"View[{view}].SaturateValues", 1)
            gmsh.option.setNumber(f"View[{view}].CustomMin", gmsh.option.getNumber(f"View[{view}].Min") + epsilon)
            gmsh.option.setNumber(f"View[{view}].CustomMax", gmsh.option.getNumber(f"View[{view}].Max"))
            gmsh.option.setNumber(f"View[{view}].ColormapNumber", 1)
            gmsh.option.setNumber(f"View[{view}].IntervalsType", 2)
            gmsh.option.setNumber(f"View[{view}].NbIso", 40)
            view += 1

        # Magnetic flux density
        gmsh.open(self.path + "/" + self.path_res_fields + "Magb.pos")
        gmsh.option.setNumber(f"View[{view}].ScaleType", 1)
        gmsh.option.setNumber(f"View[{view}].RangeType", 1)
        gmsh.option.setNumber(f"View[{view}].CustomMin", gmsh.option.getNumber(f"View[{view}].Min") + epsilon)
        gmsh.option.setNumber(f"View[{view}].CustomMax", gmsh.option.getNumber(f"View[{view}].Max"))
        gmsh.option.setNumber(f"View[{view}].ColormapNumber", 1)
        gmsh.option.setNumber(f"View[{view}].IntervalsType", 2)
        gmsh.option.setNumber(f"View[{view}].NbIso", 40)
        view += 1

        """
        # Vector Potential
        gmsh.open(self.path + "/" + self.path_res_fields + "raz.pos")
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
        gmsh.finalize()

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
        with open(self.path + "/" + self.path_res_values + loss_file, newline='') as f:
            reader = csv.reader(f)
            data = list(reader)
        return data[-last_n_values:-1] + [data[-1]]

    def load_result(self, res_name, last_n=1, part="real"):
        """
        Loads the "last_n" parameters from a result file of the scalar quantity "res_name".
        Either the real or imaginary part can be chosen.
        :param part: "real" or "imaginary" part can be chosen
        :param res_name: name of the quantity
        :param last_n: Number of parameters to be loaded
        :return: last_n entries of the chosen result file
        :rtype: list
        """
        with open(self.path + "/" + self.path_res_values + res_name + ".dat") as f:
            lines = f.readlines()[-last_n:]

            if part == "real":
                result = [float(line.split(sep=' ')[2]) for n, line in enumerate(lines)]
            if part == "imaginary":
                result = [float(line.split(sep=' ')[3]) for n, line in enumerate(lines)]

            return result

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    # FEMM [alternative Solver]
    def femm_reference(self, freq, current, sigma_cu, sign=None, non_visualize=0):
        """
        Allows reference simulations with the 2D open source electromagnetic FEM tool FEMM.
        Helpful to validate changes (especially in the Prolog Code).

        :param sign:
        :param non_visualize:
        :param freq:
        :param current:
        :param sigma_cu:

        :return:

        """
        if os.name != 'nt':
            raise Exception('You are using a computer that is not running windows. '
                          'This command is only executable on Windows computers.')

        sign = sign or [1]

        if sign is None:
            sign = [1, 1]
        if not os.path.isdir(self.path + "/" + self.path_res_FEMM):
            os.mkdir(self.path + "/" + self.path_res_FEMM)

        # == Pre Geometry ==
        self.high_level_geo_gen()

        if self.air_gaps.number != 1:
            raise NotImplementedError

        # == Init ==
        femm.openfemm(non_visualize)
        femm.newdocument(0)
        femm.mi_probdef(freq, 'meters', 'axi', 1.e-8, 0, 30)

        # == Materials ==
        femm.mi_addmaterial('Ferrite', self.core.re_mu_rel, self.core.re_mu_rel, 0, 0, 0, 0, 0, 1, 0, 0, 0)
        femm.mi_addmaterial('Air', 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0)
        if self.windings[0].conductor_type == "litz":
            femm.mi_addmaterial('Copper', 1, 1, 0, 0, sigma_cu, 0, 0, 1, 5, 0, 0, self.windings[0].n_strands,
                                2 * 1000 * self.windings[0].strand_radius)  # type := 5. last argument
            print(f"Number of strands: {self.windings[0].n_strands}")
            print(f"Diameter of strands in mm: {2 * 1000 * self.windings[0].strand_radius}")
        if self.windings[0].conductor_type == "solid":
            femm.mi_addmaterial('Copper', 1, 1, 0, 0, sigma_cu, 0, 0, 1, 0, 0, 0, 0, 0)

        # == Circuit ==
        # coil as seen from the terminals.
        femm.mi_addcircprop('Primary', current[0] * sign[0], 1)
        if self.component_type == ('transformer' or 'integrated_transformer'):
            femm.mi_addcircprop('Secondary', current[1] * sign[1], 1)

        # == Geometry ==
        # Add core
        femm.mi_drawline(0,
                         self.ei_axi.p_air_gaps[0, 1],
                         self.ei_axi.p_air_gaps[1, 0],
                         self.ei_axi.p_air_gaps[1, 1])
        femm.mi_drawline(self.ei_axi.p_air_gaps[1, 0],
                         self.ei_axi.p_air_gaps[1, 1],
                         self.ei_axi.p_window[4, 0],
                         self.ei_axi.p_window[4, 1])
        femm.mi_drawline(self.ei_axi.p_window[4, 0],
                         self.ei_axi.p_window[4, 1],
                         self.ei_axi.p_window[5, 0],
                         self.ei_axi.p_window[5, 1])
        femm.mi_drawline(self.ei_axi.p_window[5, 0],
                         self.ei_axi.p_window[5, 1],
                         self.ei_axi.p_window[7, 0],
                         self.ei_axi.p_window[7, 1])
        femm.mi_drawline(self.ei_axi.p_window[7, 0],
                         self.ei_axi.p_window[7, 1],
                         self.ei_axi.p_window[6, 0],
                         self.ei_axi.p_window[6, 1])
        femm.mi_drawline(self.ei_axi.p_window[6, 0],
                         self.ei_axi.p_window[6, 1],
                         self.ei_axi.p_air_gaps[3, 0],
                         self.ei_axi.p_air_gaps[3, 1])
        femm.mi_drawline(self.ei_axi.p_air_gaps[3, 0],
                         self.ei_axi.p_air_gaps[3, 1],
                         0,
                         self.ei_axi.p_air_gaps[2, 1])
        femm.mi_drawline(0,
                         self.ei_axi.p_air_gaps[2, 1],
                         0,
                         self.ei_axi.p_outer[2, 1])
        femm.mi_drawline(0,
                         self.ei_axi.p_outer[2, 1],
                         self.ei_axi.p_outer[3, 0],
                         self.ei_axi.p_outer[3, 1])
        femm.mi_drawline(self.ei_axi.p_outer[3, 0],
                         self.ei_axi.p_outer[3, 1],
                         self.ei_axi.p_outer[1, 0],
                         self.ei_axi.p_outer[1, 1])
        femm.mi_drawline(self.ei_axi.p_outer[1, 0],
                         self.ei_axi.p_outer[1, 1],
                         0,
                         self.ei_axi.p_outer[0, 1])
        femm.mi_drawline(0,
                         self.ei_axi.p_outer[0, 1],
                         0,
                         self.ei_axi.p_air_gaps[0, 1])

        # Add Coil

        # femm.mi_drawrectangle(self.ei_axi.p_window[4, 0]+self.isolation.core_cond,
        # self.ei_axi.p_window[4, 1]+self.component.isolation.core_cond,
        # self.ei_axi.p_window[7, 0]-self.isolation.core_cond, self.ei_axi.p_window[7, 1]-self.isolation.core_cond)
        #
        # femm.mi_addblocklabel(self.ei_axi.p_window[7, 0]-2*self.isolation.core_cond,
        # self.ei_axi.p_window[7, 1]-2*self.isolation.core_cond)
        #
        # femm.mi_selectlabel(self.ei_axi.p_window[7, 0]-2*self.isolation.core_cond,
        # self.ei_axi.p_window[7, 1]-2*self.isolation.core_cond)
        #
        # femm.mi_setblockprop('Copper', 0, 1, 'icoil', 0, 0, 1)
        #
        # femm.mi_clearselected()

        for num in range(0, self.n_windings):
            if self.windings[num].conductor_type == "litz" or self.windings[num].conductor_type == "solid":
                for i in range(0, int(self.ei_axi.p_conductor[num].shape[0] / 5)):
                    # 0: center | 1: left | 2: top | 3: right | 4.bottom
                    femm.mi_drawarc(self.ei_axi.p_conductor[num][5 * i + 1][0],
                                    self.ei_axi.p_conductor[num][5 * i + 1][1],
                                    self.ei_axi.p_conductor[num][5 * i + 3][0],
                                    self.ei_axi.p_conductor[num][5 * i + 3][1], 180, 2.5)
                    femm.mi_addarc(self.ei_axi.p_conductor[num][5 * i + 3][0],
                                   self.ei_axi.p_conductor[num][5 * i + 3][1],
                                   self.ei_axi.p_conductor[num][5 * i + 1][0],
                                   self.ei_axi.p_conductor[num][5 * i + 1][1], 180, 2.5)
                    femm.mi_addblocklabel(self.ei_axi.p_conductor[num][5 * i][0],
                                          self.ei_axi.p_conductor[num][5 * i][1])
                    femm.mi_selectlabel(self.ei_axi.p_conductor[num][5 * i][0], self.ei_axi.p_conductor[num][5 * i][1])
                    if num == 0:
                        femm.mi_setblockprop('Copper', 1, 0, 'Primary', 0, 0, 1)
                    if num == 1:
                        # femm.mi_setblockprop('Copper', 0, 1e-4, 'Secondary', 0, 0, 1)
                        femm.mi_setblockprop('Copper', 1, 0, 'Secondary', 0, 0, 1)
                    femm.mi_clearselected()

        # Define an "open" boundary condition using the built-in function:
        femm.mi_makeABC()
        """
        # Alternative BC
        region_add = 1.1

        femm.mi_drawrectangle(0, region_add*self.ei_axi.p_outer[0][1], region_add*self.ei_axi.p_outer[3][0], 
        region_add*self.ei_axi.p_outer[3][1])
        # mi_addboundprop('Asymptotic', 0, 0, 0, 0, 0, 0, 1 / (para.mu_0 * bound.width), 0, 2); % Mixed
        femm.mi_addboundprop('Asymptotic', 0, 0, 0, 0, 1, 50, 0, 0, 1)
        femm.mi_selectsegment(region_add*self.ei_axi.p_outer[3][0], region_add*self.ei_axi.p_outer[3][1])
        femm.mi_setsegmentprop('Asymptotic', 1, 1, 0, 0)
        """

        # == Labels/Designations ==

        # Label for core
        femm.mi_addblocklabel(self.ei_axi.p_outer[3, 0] - 0.001, self.ei_axi.p_outer[3, 1] - 0.001)
        femm.mi_selectlabel(self.ei_axi.p_outer[3, 0] - 0.001, self.ei_axi.p_outer[3, 1] - 0.001)
        femm.mi_setblockprop('Ferrite', 1, 0, '<None>', 0, 0, 0)
        femm.mi_clearselected()

        # Labels for air
        femm.mi_addblocklabel(0.001, 0)
        femm.mi_selectlabel(0.001, 0)
        femm.mi_setblockprop('Air', 1, 0, '<None>', 0, 0, 0)
        femm.mi_clearselected()
        femm.mi_addblocklabel(self.ei_axi.p_outer[3, 0] + 0.001, self.ei_axi.p_outer[3, 1] + 0.001)
        femm.mi_selectlabel(self.ei_axi.p_outer[3, 0] + 0.001, self.ei_axi.p_outer[3, 1] + 0.001)
        femm.mi_setblockprop('Air', 1, 0, '<None>', 0, 0, 0)
        femm.mi_clearselected()

        # Now, the finished input geometry can be displayed.
        femm.mi_zoomnatural()
        femm.mi_saveas(self.path + "/" + self.path_res_FEMM + 'coil.fem')
        femm.mi_analyze()
        femm.mi_loadsolution()

        # == Losses ==
        tmp = femm.mo_getcircuitproperties('Primary')
        self.tot_loss_femm = 0.5 * tmp[0] * tmp[1]
        print(self.tot_loss_femm)

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

    @staticmethod
    def calculate_point_average(x1, y1, x2, y2):
        """
        Used for femm_thermaL_validation
        """
        return (x1+x2)/2, (y1+y2)/2

    def femm_thermal_validation(self, thermal_conductivity_dict):
        # Get paths
        femm_model_file_path = os.path.join(os.path.dirname(th.__file__), "femm", "validation.FEH")
        results_log_file_path = os.path.join(self.path, self.path_res, "result_log.json")

        # Extract losses
        losses = th_functions.read_results_log(results_log_file_path)

        # Extract wire_radii
        wire_radii = [winding.conductor_radius for winding in self.windings]

        # TODO Needs to be calculated by the model
        core_area = 0.00077

        # == Init ==
        femm.openfemm(0)
        femm.newdocument(2)
        femm.hi_probdef("meters", "axi", 1.e-8, 0)

        # == Materials ==
        # Core
        k_core = thermal_conductivity_dict["core"]
        #q_vol_core = th_functions.calculate_heat_flux_core(losses["Core_Eddy_Current"], )
        q_vol_core = losses["Core_Eddy_Current"]/core_area
        #c_core = 0.007
        c_core = 0

        # Air
        k_air = thermal_conductivity_dict["air"]
        q_vol_air = 0
        #c_air = 1004.7
        c_air = 0

        # Wire
        k_wire = thermal_conductivity_dict["winding"]
        #c_wire = 385
        c_wire = 0

        # Case
        k_case = thermal_conductivity_dict["case"]
        q_vol_case = 0 
        #c_case = 0.01
        c_case = 0

        # Setup winding list
        winding_losses_list = []
        for i in range(1, 3):
            key = f"Winding_{i}"
            inner_winding_list = []
            if key in losses:
                for loss in losses[key]["Turns"]:
                    inner_winding_list.append(loss)
            winding_losses_list.append(inner_winding_list)

        # Setup materials
        femm.hi_addmaterial('Core', k_core, k_core, q_vol_core, c_core)
        femm.hi_addmaterial('Air', k_air, k_air, q_vol_air, c_air)
        for winding_index, winding in enumerate(winding_losses_list):
            for i in range(len(winding)):
                femm.hi_addmaterial(f'Wire_{winding_index}_{i}', k_wire, k_wire, th_functions.calculate_heat_flux_round_wire(winding[i], wire_radii[winding_index]), c_wire)
        femm.hi_addmaterial('Case', k_case, k_case, q_vol_case, c_case) 

        # Add boundary condition
        femm.hi_addboundprop("Boundary", 0, 273, 0, 0, 0, 0)
        femm.hi_setsegmentprop("Boundary", 0, 1, 0, 2, "<None>")

        # == Geometry ==
        self.high_level_geo_gen()
        # Add core
        femm.hi_drawline(0,
                        self.ei_axi.p_air_gaps[0, 1],
                        self.ei_axi.p_air_gaps[1, 0],
                        self.ei_axi.p_air_gaps[1, 1])
        femm.hi_drawline(self.ei_axi.p_air_gaps[1, 0],
                        self.ei_axi.p_air_gaps[1, 1],
                        self.ei_axi.p_window[4, 0],
                        self.ei_axi.p_window[4, 1])
        femm.hi_drawline(self.ei_axi.p_window[4, 0],
                        self.ei_axi.p_window[4, 1],
                        self.ei_axi.p_window[5, 0],
                        self.ei_axi.p_window[5, 1])
        femm.hi_drawline(self.ei_axi.p_window[5, 0],
                        self.ei_axi.p_window[5, 1],
                        self.ei_axi.p_window[7, 0],
                        self.ei_axi.p_window[7, 1])
        femm.hi_drawline(self.ei_axi.p_window[7, 0],
                        self.ei_axi.p_window[7, 1],
                        self.ei_axi.p_window[6, 0],
                        self.ei_axi.p_window[6, 1])
        femm.hi_drawline(self.ei_axi.p_window[6, 0],
                        self.ei_axi.p_window[6, 1],
                        self.ei_axi.p_air_gaps[3, 0],
                        self.ei_axi.p_air_gaps[3, 1])
        femm.hi_drawline(self.ei_axi.p_air_gaps[3, 0],
                        self.ei_axi.p_air_gaps[3, 1],
                        0,
                        self.ei_axi.p_air_gaps[2, 1])
        femm.hi_drawline(0,
                        self.ei_axi.p_air_gaps[2, 1],
                        0,
                        self.ei_axi.p_outer[2, 1])
        femm.hi_drawline(0,
                        self.ei_axi.p_outer[2, 1],
                        self.ei_axi.p_outer[3, 0],
                        self.ei_axi.p_outer[3, 1])
        femm.hi_drawline(self.ei_axi.p_outer[3, 0],
                        self.ei_axi.p_outer[3, 1],
                        self.ei_axi.p_outer[1, 0],
                        self.ei_axi.p_outer[1, 1])
        femm.hi_drawline(self.ei_axi.p_outer[1, 0],
                        self.ei_axi.p_outer[1, 1],
                        0,
                        self.ei_axi.p_outer[0, 1])      
        femm.hi_drawline(0,
                        self.ei_axi.p_outer[0, 1],
                        0,
                        self.ei_axi.p_air_gaps[0, 1])

        # In order for the simulation to work the air_gap must be closed:
        femm.hi_drawline(0, self.ei_axi.p_air_gaps[0, 1], 0, self.ei_axi.p_air_gaps[2, 1])

        # Add case
        # Case size
        case_gap_top = 0.0015
        case_gap_right = 0.0025
        case_gap_bot = 0.002

        femm.hi_drawline(0, self.ei_axi.p_outer[2, 1], 0, self.ei_axi.p_outer[2, 1] + case_gap_top) # Top left line
        femm.hi_drawline(0, self.ei_axi.p_outer[2, 1] + case_gap_top, self.ei_axi.p_outer[3, 0] + case_gap_right, self.ei_axi.p_outer[3, 1] + case_gap_top) # Top line
        femm.hi_drawline(self.ei_axi.p_outer[3, 0] + case_gap_right, self.ei_axi.p_outer[3, 1] + case_gap_top, self.ei_axi.p_outer[1, 0] + case_gap_right, self.ei_axi.p_outer[1, 1] - case_gap_bot) # Right line
        femm.hi_drawline(self.ei_axi.p_outer[1, 0] + case_gap_right, self.ei_axi.p_outer[1, 1] - case_gap_bot, 0, self.ei_axi.p_outer[0, 1] - case_gap_bot) # Bottom line
        femm.hi_drawline(0, self.ei_axi.p_outer[0, 1] - case_gap_bot, 0, self.ei_axi.p_outer[0, 1]) # Bottom right line

        # Create boundary
        #femm.hi_selectsegment(*self.calculatePointAverage(0, self.ei_axi.p_outer[2, 1], 0, self.ei_axi.p_outer[2, 1] + caseGapTop))
        femm.hi_selectsegment(*MagneticComponent.calculate_point_average(0, self.ei_axi.p_outer[2, 1] + case_gap_top, self.ei_axi.p_outer[3, 0] + case_gap_right, self.ei_axi.p_outer[3, 1] + case_gap_top))
        femm.hi_selectsegment(*MagneticComponent.calculate_point_average(self.ei_axi.p_outer[3, 0] + case_gap_right, self.ei_axi.p_outer[3, 1] + case_gap_top, self.ei_axi.p_outer[1, 0] + case_gap_right, self.ei_axi.p_outer[1, 1] - case_gap_bot))
        femm.hi_selectsegment(*MagneticComponent.calculate_point_average(self.ei_axi.p_outer[1, 0] + case_gap_right, self.ei_axi.p_outer[1, 1] - case_gap_bot, 0, self.ei_axi.p_outer[0, 1] - case_gap_bot))
        #femm.hi_selectsegment(*self.calculatePointAverage(0, self.ei_axi.p_outer[0, 1] - caseGapBot, 0, self.ei_axi.p_outer[0, 1]))
        femm.hi_setsegmentprop("Boundary", 0, 1, 0, 2, "<None>")

        # Add case material
        material_x, material_y = self.calculate_point_average(0, self.ei_axi.p_outer[2, 1], 0, self.ei_axi.p_outer[2, 1] + case_gap_top)
        femm.hi_addblocklabel(material_x + 0.001, material_y)
        femm.hi_selectlabel(material_x + 0.001, material_y)
        femm.hi_setblockprop('Case', 1, 0, 0)
        femm.hi_clearselected()

        # Add Coil
        for num in range(0, self.n_windings):
            for i in range(0, int(self.ei_axi.p_conductor[num].shape[0] / 5)):
                # 0: center | 1: left | 2: top | 3: right | 4.bottom
                femm.hi_drawarc(self.ei_axi.p_conductor[num][5 * i + 1][0],
                                self.ei_axi.p_conductor[num][5 * i + 1][1],
                                self.ei_axi.p_conductor[num][5 * i + 3][0],
                                self.ei_axi.p_conductor[num][5 * i + 3][1], 180, 2.5)
                femm.hi_addarc(self.ei_axi.p_conductor[num][5 * i + 3][0],
                            self.ei_axi.p_conductor[num][5 * i + 3][1],
                            self.ei_axi.p_conductor[num][5 * i + 1][0],
                            self.ei_axi.p_conductor[num][5 * i + 1][1], 180, 2.5)
                femm.hi_addblocklabel(self.ei_axi.p_conductor[num][5 * i][0],
                                    self.ei_axi.p_conductor[num][5 * i][1])
                femm.hi_selectlabel(self.ei_axi.p_conductor[num][5 * i][0], self.ei_axi.p_conductor[num][5 * i][1])
                if num == 0:
                    femm.hi_setblockprop(f'Wire_{num}_{i}', 1, 0, 1)
                if num == 1:
                    femm.hi_setblockprop(f'Wire_{num}_{i}', 1, 0, 1)
                femm.hi_clearselected()

        # Define an "open" boundary condition using the built-in function:
        #femm.hi_makeABC()


        # == Labels/Designations ==
        # Label for core
        femm.hi_addblocklabel(self.ei_axi.p_outer[3, 0] - 0.001, self.ei_axi.p_outer[3, 1] - 0.001)
        femm.hi_selectlabel(self.ei_axi.p_outer[3, 0] - 0.001, self.ei_axi.p_outer[3, 1] - 0.001)
        femm.hi_setblockprop('Core', 1, 0, 0)
        femm.hi_clearselected()

        # Labels for air
        femm.hi_addblocklabel(0.001, 0)
        femm.hi_selectlabel(0.001, 0)
        femm.hi_setblockprop('Air', 1, 0, 0)
        femm.hi_clearselected()

        # Not needed when the core is the boundary
        #femm.hi_addblocklabel(self.ei_axi.p_outer[3, 0] + 0.001, self.ei_axi.p_outer[3, 1] + 0.001)
        #femm.hi_selectlabel(self.ei_axi.p_outer[3, 0] + 0.001, self.ei_axi.p_outer[3, 1] + 0.001)
        #femm.hi_setblockprop('Air', 1, 0, 0)
        #femm.hi_clearselected()

        # Now, the finished input geometry can be displayed.
        femm.hi_zoomnatural()
        femm.hi_saveas(femm_model_file_path)
        femm.hi_analyze()
        femm.hi_loadsolution()
        input() # So the window stays open
        # femm.closefemm()

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    # Litz Approximation [internal methods]
    def pre_simulate(self):
        """
        Used to determine the litz-approximation coefficients.

        :return:

        """
        for num in range(0, self.n_windings):
            if self.windings[num].conductor_type == 'litz':
                # ---
                # Litz Approximation Coefficients were created with 4 layers
                # That's why here a hard-coded 4 is implemented
                # if os.path.isfile(self.path +
                # f"/Strands_Coefficents/coeff/pB_RS_la{self.windings[num].ff}_{self.n_layers[num]}layer.dat"):
                if os.path.isfile(self.path + f"/Strands_Coefficents/coeff/pB_RS_la{self.windings[num].ff}_4layer.dat"):
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
        text_file = open(self.path + "/Strands_Coefficents/PreParameter.pro", "w")
        # ---
        # Litz Approximation Coefficients are created with 4 layers
        # That's why here a hard-coded 4 is implemented
        # text_file.write(f"NbrLayers = {self.n_layers[num]};\n")
        text_file.write(f"NbrLayers = 4;\n")
        text_file.write(f"Fill = {self.windings[num].ff};\n")
        print("Here")
        text_file.write(f"Rc = {self.windings[num].strand_radius};\n")  # double named!!! must be changed
        text_file.close()
        cell_geo = self.onelab_client.getPath('Strands_Coefficents/cell.geo')

        # Run gmsh as a sub client
        mygmsh = self.onelab + 'gmsh'
        self.onelab_client.runSubClient('myGmsh', mygmsh + ' ' + cell_geo + ' -2 -v 2')

        modes = [1, 2]  # 1 = "skin", 2 = "proximity"
        reduced_frequencies = np.linspace(0, 1.25, 6)  # must be even
        for mode in modes:
            for rf in reduced_frequencies:
                # -- Pre-Simulation Settings --
                text_file = open(self.path + "/Strands_Coefficents/PreParameter.pro", "w")
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
                input_file = self.onelab_client.getPath('Strands_Coefficents/cell_dat.pro')
                cell = self.onelab_client.getPath('Strands_Coefficents/cell.pro')

                # Run simulations as sub clients
                mygetdp = self.onelab + 'getdp'
                self.onelab_client.runSubClient('myGetDP', mygetdp + ' ' + cell + ' -input ' + input_file + ' -solve MagDyn_a -v2')

        # Formatting stuff
        # Litz Approximation Coefficients are created with 4 layers
        # That's why here a hard-coded 4 is implemented
        # files = [self.path + f"/Strands_Coefficents/coeff/pB_RS_la{self.windings[num].ff}_{self.n_layers[num]}layer.dat",
        #         self.path + f"/Strands_Coefficents/coeff/pI_RS_la{self.windings[num].ff}_{self.n_layers[num]}layer.dat",
        #         self.path + f"/Strands_Coefficents/coeff/qB_RS_la{self.windings[num].ff}_{self.n_layers[num]}layer.dat",
        #         self.path + f"/Strands_Coefficents/coeff/qI_RS_la{self.windings[num].ff}_{self.n_layers[num]}layer.dat"]
        files = [self.path + f"/Strands_Coefficents/coeff/pB_RS_la{self.windings[num].ff}_4layer.dat",
                 self.path + f"/Strands_Coefficents/coeff/pI_RS_la{self.windings[num].ff}_4layer.dat",
                 self.path + f"/Strands_Coefficents/coeff/qB_RS_la{self.windings[num].ff}_4layer.dat",
                 self.path + f"/Strands_Coefficents/coeff/qI_RS_la{self.windings[num].ff}_4layer.dat"]
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
        self.excitation(f=100000, i=1)  # arbitrary values: frequency and current
        self.file_communication()
        self.pre_simulate()

    #  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    # Standard Simulations
    def single_simulation(self, freq: float, current: List[float], phi_deg: List[float] = None,
                          skin_mesh_factor: float = 1, NL_core=0) -> None:
        """

        Start a _single_ ONELAB simulation.
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

        self.high_level_geo_gen(frequency=freq, skin_mesh_factor=skin_mesh_factor)
        if self.valid:
            self.mesh.generate_mesh()
            self.excitation(f=freq, i=current, phases=phi_deg)  # frequency and current
            self.file_communication()
            self.pre_simulate()
            self.simulate()
            self.write_log()
            self.visualize()

        #results =

        #return results

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
            os.remove(self.path + "/" + self.path_res_values + "L_11.dat")
            os.remove(self.path + "/" + self.path_res_values + "L_22.dat")
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
            with open(self.path + "/" + self.path_res_values + "Flux_Linkage_1.dat") as f:
                line = f.readlines()[-2:]
                # Fluxes induced in Winding 1
                Phi_11 = float(line[0].split(sep=' ')[2])
                Phi_12 = float(line[1].split(sep=' ')[2])

            with open(self.path + "/" + self.path_res_values + "Flux_Linkage_2.dat") as f:
                line = f.readlines()[-2:]
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
            with open(self.path + "/" + self.path_res_values + "L_11.dat") as f:
                line = f.readlines()[-1]
                words = line.split(sep=' ')
                self.L_11 = float(words[2])
            with open(self.path + "/" + self.path_res_values + "L_22.dat") as f:
                line = f.readlines()[-1]
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
            self.excitation(f=frequencies[i], i=currents[i], phases=phi[i])  # frequency and current
            self.file_communication()
            self.pre_simulate()
            self.simulate()
        if show_last:
            self.visualize()

    def excitation_sweep(self, frequencies, currents, phi, show_last=False, return_results=False, meshing=True):
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
        if show_last:
            self.plot_fields = "standard"
        else:
            self.plot_fields = False

        if meshing:
            self.high_level_geo_gen(frequency=frequencies[0])  # TODO: Must be changed for solid sim.
            if self.valid:
                self.mesh.generate_mesh()

        if self.valid:

            for i in range(0, len(frequencies)):
                self.excitation(f=frequencies[i], i=currents[i], phases=phi[i])  # frequency and current
                self.file_communication()
                self.pre_simulate()
                self.simulate()

            if show_last:
                self.visualize()

            if return_results:
                # Return a Dictionary with the results
                results = {"j2F": self.load_result("j2F", len(frequencies), "real"),
                           "j2H": self.load_result("j2H", len(frequencies), "real"),
                           "p_hyst": self.load_result("p_hyst", len(frequencies), "real")}
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
        self.excitation(f=f_switch, i=Ipeak, phases=[0, 180])  # frequency and current
        self.file_communication()
        self.pre_simulate()
        self.simulate()
        self.visualize()
