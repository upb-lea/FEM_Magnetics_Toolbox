# Usual Python libraries
import csv
import fileinput
import numpy.typing as npt
import numpy as np
import os
import pandas as pd
import pathlib
import re
import sys
import time
import warnings
from matplotlib import pyplot as plt
import gmsh
# import adapt_mesh
from onelab import onelab
import json
import random
import string
import subprocess
import pkg_resources
from typing import Union
# Self written functions. It is necessary to write a . before the function, due to handling
# this package also as a pip-package
# from .femmt_functions import id_generator, inner_points, min_max_inner_points, call_for_path, NbrStrands


def install_femm_if_missing():
    required = {'pyfemm'}
    installed = {pkg.key for pkg in pkg_resources.working_set}
    missing = required - installed

    if missing:
        print("Missing 'pyfemm' installation.")
        print("Installing 'pyfemm' ...")
        python = sys.executable
        subprocess.check_call([python, '-m', 'pip', 'install', *missing], stdout=subprocess.DEVNULL)
        print("'pyfemm' is now installed!")


# Optional usage of FEMM tool by David Meeker
# 2D Mesh and FEM interfaces (only for windows machines)
if os.name == 'nt':
    install_femm_if_missing()
    import femm


#  ===== Main Class  =====
class MagneticComponent:
    """
    - Initialization of all class variables
    - Common variables for all instances
    """
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
        - Initialization of all instance variables
        - One or more "MagneticComponents" can be created
        - Each "MagneticComponent" owns its own instance variable values
        :param component_type: Available options are "inductor" and "transformer"
        """
        print(f"\n"
              f"Initialized a new Magnetic Component of type {component_type}\n"
              f"--- --- --- ---")

        # Breaking variable
        self.valid = True

        # ==============================
        # Geometry
        # ==============================
        # -- Control Flags --
        self.region = None # Apply an outer Region or directly apply a constraint on the Core Boundary
        self.padding = 1.5 # ... > 1
        self.y_symmetric = 1  # Mirror-symmetry across y-axis
        self.dimensionality = "2D axi"  # Axial-symmetric model (idealized full-cylindrical)
        self.s = 0.5  # Parameter for mesh-accuracy
        self.component_type = component_type  # "inductor" or "transformer"

        # -- Core --
        self.core_type = "EI"  # Basic shape of magnetic conductor
        self.core_w = None  # Axi symmetric case | core_w := core radius
        self.window_w = None  # Winding window width
        self.window_h = None  # Winding window height
        self.core_update_count = 0
        self.update_core(core_type="EI", core_w=0.02, window_w=0.01, window_h=0.03)  # some initial values

        # -- Air gaps --
        self.n_air_gaps = 1  # Number of air gaps [==1: air gap in center | >1: random air gaps]
        self.air_gaps = np.empty((self.n_air_gaps, 4))  # list: [position_tag, air_gap_position, air_gap_h, c_air_gap]

        # -- Dedicated Stray Path --
        self.dedicated_stray_path = False
        self.start_i = None
        #TODO: Thickness of the stray path must be fitted for the real Tablet (effective area of the "stray air gap" is
        # different in axialsymmetric approximation
        self.end_i = None
        self.added_bar = None

        # -- Windings --
        self.n_windings = None  # Number of conductors/windings

        if component_type == "inductor":
            self.n_windings = 1
            self.windings = [self.Winding()]

        if component_type == "transformer":
            self.n_windings = 2
            self.windings = [self.Winding(), self.Winding()]

        if component_type == "three_phase_transformer":
            self.n_windings = 3
            self.windings = [self.Winding(), self.Winding(), self.Winding()]
            raise NotImplemented


        # -- Virtual Winding Windows --
        self.virtual_winding_windows = None
        self.vw_type = None  # "center" and "full_window" are the only cases implemented yet; #TODO: ersetzen

        # -- Isolation ---
        self.core_cond_isolation = [None] * self.n_windings  # gap between Core and each Conductor
        self.cond_cond_isolation = [None] * (self.n_windings * 2 - 1)  # \n
        # first n_conductor arguments: isolation gap between two turns of common conductors
        # last n_conductor-1 arguments: gap between two neighboured conductors
        # 12, 13, 23
        #TODO: (n-1)! = (n-1)*(n-2)...*1 + 1

        # -- Geometric Parameters/Coordinates --
        self.n_windows = None
        self.p_outer = None
        self.p_window = None
        self.p_conductor = []
        for i in range(0, self.n_windings):
            self.p_conductor.insert(i, [])
        self.p_air_gaps = None

        # ==============================
        # Materials
        # ==============================
        # frequency = 0: mu_rel only used if flag_non_linear_core == 0
        # frequency > 0: mu_rel is used
        self.cond_sigma = 5.8e7  # perfect copper
        self.mu0 = np.pi*4e-7
        # TDk N95 as standard material:
        self.core_re_mu_rel = 3000  # Real part of relative Core Permeability  [B-Field and frequency-dependend]
        self.core_im_mu_rel = 1750 * np.sin(10 *np.pi/180)   # Imaginary part of relative Core Permeability  [B-Field and frequency-dependend]
        self.core_im_epsilon_rel = 6e+4 * np.sin(20 *np.pi/180)  # Imaginary part of complex equivalent permeability  [only frequency-dependend]
        self.core_material = 95_100  # 95  # 95 := TDK-N95 | Currently only works with Numbers corresponding to BH.pro

        # ==============================
        # Problem Definition
        # ==============================
        # -- Excitation Parameters --
        self.flag_imposed_reduced_frequency = None
        self.flag_excitation_type = None
        self.flag_non_linear_core = None
        self.current = [None] * self.n_windings  # Defined for every conductor
        self.current_density = [None] * self.n_windings  # Defined for every conductor
        self.voltage = [None] * self.n_windings  # Defined for every conductor
        self.frequency = None
        self.phase_tmp = np.zeros(self.n_windings)  # Default is zero, Defined for every conductor
        self.red_freq = [None] * self.n_windings  # Defined for every conductor
        self.delta = None

        # ==============================
        # Meshing
        # ==============================
        # -- Characteristic lengths -- [for mesh sizes]
        self.skin_mesh_factor = None
        self.c_core = None
        # self.c_core =  self.core_w/10. * self.s
        self.c_window = None
        self.c_conductor = [None] * self.n_windings  # self.delta  # self.s /20 #self.window_w/30 * self.s
        self.c_center_conductor = [None] * self.n_windings  # used for the mesh accuracy in the conductors
        self.c_air_gap = []
        self.island_mesh_acc = 1  # will be multiplied with self.c_window later
        self.zero_mesh_acc = self.island_mesh_acc  # points that are pulled to x=0

        # -- Used for Litz Validation --
        self.sweep_frequencies = None

        # -- Core Loss --
        self.core_loss_simlation = 0
        self.Ipeak = None
        self.ki = None
        self.alpha = None
        self.beta = None
        self.t_rise = None
        self.t_fall = None
        self.f_switch = None

        # -- Results --
        self.L_11 = None
        self.L_22 = None
        self.M = None
        self.Pv = None
        # - Primary Concentrated
        self.ü_conc = None
        self.L_s_conc = None
        self.L_h_conc = None

        # -- FEMM variables --
        self.tot_loss_femm = None

    # ==== Back-End Methods ====
    # Methods that are usually not directly called from outside

    # === Setup ===
    def onelab_setup(self):
        """
        Either reads onelab parent folder path from config.p or asks the user to provide it.
        Creates a config.p at first run.
        :return: -
        """
        # find out path of femmt (installed module or directly opened in git)?
        module_file_path = pathlib.Path(__file__).parent.absolute()
        config_file_path = module_file_path / 'config.json'

        # check if config.json is available and not empty
        if pathlib.Path.is_file(config_file_path) and pathlib.Path.stat(config_file_path).st_size != 0:
            json_file = config_file_path.open('rb')
            loaded_dict = json.load(json_file)
            json_file.close()
            path = loaded_dict['onelab']

            if os.path.exists(path):
                self.onelab = path
            else:
                self.onelab = call_for_path("onelab")
        else:
            self.onelab = call_for_path("onelab")

    # === Geometry ===
    def high_level_geo_gen(self, core_type="EI", dimensionality="2D axi", frequency=None, skin_mesh_factor=1):
        """
        - high level geometry generation
        - based on chosen core and conductor types and simulation mode
        - calls "low level" methods, that creates all points needed for mesh generation
        :return:
        """
        # ==============================
        # High-Level Geometry Generation
        # ==============================

        # Mesh-Parameters must be updated depending on geometry size
        self.c_core =  self.core_w/10. * self.s
        self.c_window = self.window_w/30 * self.s

        # Update Skin Depth (needed for meshing)
        self.skin_mesh_factor = skin_mesh_factor
        if frequency != None:
            if frequency == 0:
                self.delta = 1e9
            else:
                self.delta = np.sqrt(2 / (2 * frequency * np.pi * self.cond_sigma * self.mu0))
            for i in range(0, self.n_windings):
                if self.windings[i].conductor_type == "solid":
                    self.c_conductor[i] = min([self.delta * self.skin_mesh_factor,
                                               self.windings[i].conductor_radius / 4 * self.skin_mesh_factor])
                    self.c_center_conductor[i] = self.windings[i].conductor_radius / 4 * self.skin_mesh_factor
                    #print(f"Werte Leiter: {self.c_conductor, self.c_center_conductor}")
                elif self.windings[i].conductor_type == "litz":
                    self.c_conductor[i] = self.windings[i].conductor_radius / 4 * self.skin_mesh_factor
                    self.c_center_conductor[i] = self.windings[i].conductor_radius / 4 * self.skin_mesh_factor
                    #print(f"Werte Leiter: {self.c_conductor, self.c_center_conductor}")
                else:
                    self.c_conductor[i] = 0.0001  # revisit

        # -- Core-type --
        if self.core_type == core_type:
            self.n_windows = 2

        # -- Symmetry -- [Choose between asymmetric, symmetric and axi symmetric]
        self.dimensionality = dimensionality
        if self.dimensionality == "2D axi":
            self.y_symmetric = 1

        if self.core_type == "EI" and self.dimensionality == "2D axi":
            self.ei_axi()

    class VirtualWindingWindow:
        """
        A virtual winding window is the area, where either some kind of interleaved conductors or a one winding
        (primary, secondary,...) is placed in a certain way.
        """
        def __init__(self):
            ## Rectangular frame:
            self.bot_bound = None
            self.top_bound = None
            self.left_bound = None
            self.right_bound = None

            # Arrangement of the Conductors in the virtual winding window
            # Obviously depends on the chosen conductor type
            self.winding = None  # "interleaved" | "primary"/"secondary"
            self.scheme = None   # "bifilar", "vertical", "horizontal", ["hexa", "square"] | "hexa", "square"
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

    def ei_axi(self):
        """
        - creates all points needed for the radial axi-symmetric EI core typology
        :return:
        """
        # -- Air Gap Data -- [random air gap generation]
        #self.update_air_gaps()

        # -- Arrays for geometry data -- [all points with (x, y, z, mesh_accuracy)]
        self.p_outer = np.zeros((4, 4))
        self.p_region_bound = np.zeros((4, 4))
        self.p_window = np.zeros((4*self.n_windows, 4))
        self.p_air_gaps = np.zeros((4*self.n_air_gaps, 4))

        # -- Geometry data --

        # Fitting the outer radius to ensure surface area
        r_inner = self.window_w + self.core_w/2
        r_outer = np.sqrt((self.core_w/2)**2 + r_inner**2)  # np.sqrt(window_w**2 + window_w * core_w + core_w**2/2)

        # Outer Core
        # (A_zyl=2pi*r*h => h=0.5r=0.25core_w <=> ensure A_zyl=A_core on the tiniest point)
        self.p_outer[0][:] = [-r_outer, -(self.window_h / 2 + self.core_w/4), 0, self.c_core]
        self.p_outer[1][:] = [r_outer, -(self.window_h / 2 + self.core_w/4), 0, self.c_core]
        self.p_outer[2][:] = [-r_outer, (self.window_h / 2 + self.core_w/4), 0, self.c_core]
        self.p_outer[3][:] = [r_outer, (self.window_h / 2 + self.core_w/4), 0, self.c_core]

        # Window
        # At this point both windows (in a cut) are modeled
        # print(f"win: c_window: {self.c_window}")
        self.p_window[0] = [-r_inner, -self.window_h/2, 0, self.c_window]
        self.p_window[1] = [-self.core_w/2, -self.window_h/2, 0, self.c_window]
        self.p_window[2] = [-r_inner, self.window_h/2, 0, self.c_window]
        self.p_window[3] = [-self.core_w/2, self.window_h/2, 0, self.c_window]
        self.p_window[4] = [self.core_w/2, -self.window_h/2, 0, self.c_window]
        self.p_window[5] = [r_inner, -self.window_h/2, 0, self.c_window]
        self.p_window[6] = [self.core_w/2, self.window_h/2, 0, self.c_window]
        self.p_window[7] = [r_inner, self.window_h/2, 0, self.c_window]


        # - Air gaps -
        # "air_gaps" is a list with [position_tag, air_gap_position, air_gap_h, c_air_gap]
        #   - position_tag: specifies the gapped "leg"
        #   - air_gap_position: specifies the coordinate of the air gap's center point along the specified leg
        #   - air_gap_h: height/length of the air gap
        #   - c_air_gap: mesh accuracy factor
        # at this point the 4 corner points of each air gap are generated out of "air_gaps"
        for i in range(0, self.n_air_gaps):
            """
            # Left leg (-1)
            if self.air_gaps[i][0] == -1:
                self.p_air_gaps[i * 4] = [-(self.core_w + self.window_w), self.air_gaps[i][1] - self.air_gaps[i][2] / 2, 0, self.air_gaps[i][3]]
                self.p_air_gaps[i * 4 + 1] = [-(self.core_w / 2 + self.window_w), self.air_gaps[i][1] - self.air_gaps[i][2] / 2, 0, self.air_gaps[i][3]]
                self.p_air_gaps[i * 4 + 2] = [-(self.core_w + self.window_w), self.air_gaps[i][1] + self.air_gaps[i][2] / 2, 0, self.air_gaps[i][3]]
                self.p_air_gaps[i * 4 + 3] = [-(self.core_w / 2 + self.window_w), self.air_gaps[i][1] + self.air_gaps[i][2] / 2, 0, self.air_gaps[i][3]]

            # Right leg (+1)
            if self.air_gaps[i][0] == 1:
                self.p_air_gaps[i * 4] = [self.core_w / 2 + self.window_w, self.air_gaps[i][1] - self.air_gaps[i][2] / 2, 0, self.air_gaps[i][3]]
                self.p_air_gaps[i * 4 + 1] = [self.core_w + self.window_w, self.air_gaps[i][1] - self.air_gaps[i][2] / 2, 0, self.air_gaps[i][3]]
                self.p_air_gaps[i * 4 + 2] = [self.core_w / 2 + self.window_w, self.air_gaps[i][1] + self.air_gaps[i][2] / 2, 0, self.air_gaps[i][3]]
                self.p_air_gaps[i * 4 + 3] = [self.core_w + self.window_w, self.air_gaps[i][1] + self.air_gaps[i][2] / 2, 0, self.air_gaps[i][3]]
            """
            # Center leg (0)
            if self.air_gaps[i][0] == 0:
                print(self.air_gaps[i][3])
                #TODO: sadly the center points are passed by update_air_gaps() and at this point transformed into 4 corner points
                self.p_air_gaps[i * 4 + 0] = [-self.core_w/2, self.air_gaps[i][1] - self.air_gaps[i][2] / 2, 0, self.air_gaps[i][3]]
                self.p_air_gaps[i * 4 + 1] = [ self.core_w/2, self.air_gaps[i][1] - self.air_gaps[i][2] / 2, 0, self.air_gaps[i][3]]
                self.p_air_gaps[i * 4 + 2] = [-self.core_w/2, self.air_gaps[i][1] + self.air_gaps[i][2] / 2, 0, self.air_gaps[i][3]]
                self.p_air_gaps[i * 4 + 3] = [ self.core_w/2, self.air_gaps[i][1] + self.air_gaps[i][2] / 2, 0, self.air_gaps[i][3]]


        # - Virtual Windows -
        #TODO: make this part of the class VWW...
        # self.vw_type must be an input or so to that class
        separation_hor = 0 #self.window_h * 0.5
        separation_vert = self.window_w * 0.5

        if self.dedicated_stray_path == False:
            # Some examples for virtual windows
            # Concentrated windings

            if self.vw_type == "full_window":
                """
                The winding window is completely used by one VWW.
                In case of a transformer, an interleaved winding scheme must be used.
                """
                # top window
                min =-self.window_h / 2 + self.core_cond_isolation[0] / 2  # bottom
                max = self.window_h / 2 - self.core_cond_isolation[0]  # top
                left = self.core_w / 2 + self.core_cond_isolation[0]
                right = r_inner - self.core_cond_isolation[0]

                # Sum the windows up in a list
                #self.virtual_windows = [[min, max, left, right]]
                self.virtual_winding_windows[0].bot_bound = min
                self.virtual_winding_windows[0].top_bound = max
                self.virtual_winding_windows[0].left_bound = left
                self.virtual_winding_windows[0].right_bound = right

            if self.vw_type == "center":
                """
                The winding window is split into two VWWs.
                The primary winding is placed in the upper half,
                the secondary winding is placed in the lower half of the winding window
                """
                separation_hor = 0

                # top window
                min21 = -separation_hor + self.cond_cond_isolation[-1] / 2  # separation_hor
                max21 = self.window_h / 2 - self.core_cond_isolation[0]  # top
                left21 = self.core_w / 2 + self.core_cond_isolation[0]
                right21 = r_inner - self.core_cond_isolation[0]

                # bottom window
                min11 =-self.window_h / 2 + self.core_cond_isolation[0] / 2  # bottom
                max11 = -separation_hor - self.cond_cond_isolation[-1] / 2  # separation_hor
                left11 = self.core_w / 2 + self.core_cond_isolation[0]
                right11 = r_inner - self.core_cond_isolation[0]

                # Sum the windows up in a list
                virtual_windows = [[min11, max11, left11, right11], [min21, max21, left21, right21]]
                for vww in range (0, len(virtual_windows)):
                    self.virtual_winding_windows[vww].bot_bound = virtual_windows[vww][0]
                    self.virtual_winding_windows[vww].top_bound = virtual_windows[vww][1]
                    self.virtual_winding_windows[vww].left_bound = virtual_windows[vww][2]
                    self.virtual_winding_windows[vww].right_bound = virtual_windows[vww][3]

            if self.vw_type == "something_else":
                # bottom left window
                min11 = -self.window_h / 2 + self.core_cond_isolation[0] / 2  # bottom
                max11 = -separation_hor - self.cond_cond_isolation[-1] / 2  # separation_hor
                left11 = self.core_w / 2 + self.core_cond_isolation[0]
                right11 = r_inner - self.cond_cond_isolation[0] - separation_vert

                # bottom right window
                min12 =-self.window_h / 2 + self.core_cond_isolation[0] / 2  # bottom
                max12 = -separation_hor - self.cond_cond_isolation[-1] / 2  # separation_hor
                left12 = r_inner + self.cond_cond_isolation[0] - separation_vert
                right12 = r_inner - self.core_cond_isolation[0]

                # top window
                min21 = -separation_hor + self.cond_cond_isolation[-1] / 2  # separation_hor
                max21 = self.window_h / 2 - self.core_cond_isolation[0]  # top
                left21 = self.core_w / 2 + self.core_cond_isolation[0]
                right21 = r_inner - self.core_cond_isolation[0]

                # Sum the windows up in a list
                virtual_windows = [[min11, max11, left11, right11], [min12, max12, left12, right12],
                                        [min21, max21, left21, right21]]
                #TODO:More flexible virtual winging windows


        # With dedicated stray path:
        if self.dedicated_stray_path == True:
            """
            If dedicated stray path is the chosen typology the are two winding windows
            These can either be split up into more virtual windows or (in case of bifilar windings) not
            """
            #TODO: Unterteilung in weitere virtuelle Fenster

            # top window
            island_right_tmp = inner_points(self.p_window[4], self.p_window[6], self.p_air_gaps)
            min11 = -self.window_h / 2 + self.core_cond_isolation[0] / 2  # bottom
            max11 = island_right_tmp[(self.start_i - 1) * 2][1] - self.core_cond_isolation[0] / 2  # separation_hor
            left11 = self.core_w / 2 + self.core_cond_isolation[0]
            right11 = r_inner - self.core_cond_isolation[0]

            # bot window
            min21 = island_right_tmp[(self.start_i - 1) * 2 + 1][1] + self.core_cond_isolation[0] / 2  # separation_hor
            max21 = self.window_h / 2 - self.core_cond_isolation[0]  # top
            left21 = self.core_w / 2 + self.core_cond_isolation[0]
            right21 = r_inner - self.core_cond_isolation[0]

            # Store the window boarders in the VWW objects
            virtual_windows = [[min21, max21, left21, right21], [min11, max11, left11, right11]]
            for vww in range(0, len(virtual_windows)):
                self.virtual_winding_windows[vww].bot_bound = virtual_windows[vww][0]
                self.virtual_winding_windows[vww].top_bound = virtual_windows[vww][1]
                self.virtual_winding_windows[vww].left_bound = virtual_windows[vww][2]
                self.virtual_winding_windows[vww].right_bound = virtual_windows[vww][3]


        # - Conductors -
        for n_win in range(0, len(self.virtual_winding_windows)):
            """
            - Work through the virtual winding windows
            - Usual cases are one window for classical transformers or two windows for transformers with a dedicated 
              stray path
            - There can be as many virtual winding windows as the user wants to define...
            - To automatically fill a virtual winding window with windings, #TODO: self.interleaving[n_win] can be chosen
              to "bifilar", "vertical", "horizontal", ["hexa", "square"] or completely with one of the windings by "primary"
              or "secondary"
            """
            # Boarders of the VWW:
            # start with the lower one
            bot_bound = self.virtual_winding_windows[n_win].bot_bound
            top_bound = self.virtual_winding_windows[n_win].top_bound
            left_bound = self.virtual_winding_windows[n_win].left_bound
            right_bound = self.virtual_winding_windows[n_win].right_bound

            if self.virtual_winding_windows[n_win].winding == "interleaved":

                if self.virtual_winding_windows[n_win].scheme == "bifilar":
                    """
                    - Bifilar interleaving means a uniform winding scheme of two conductors (prim. and sec.)
                    - Can only be used for conductors of identical radius (in terms of litz radius for stranded wires)
                    - Excess windings are placed below the bifilar ones
                    """
                    if self.windings[0].conductor_radius != self.windings[1].conductor_radius:
                        print("For bifilar winding scheme both conductors must be of the same radius!")
                    else:
                        print("Bifilar winding scheme is applied")
                        #ü1 = self.n_turns[0]/self.n_turns[1]
                        #ü2 = self.n_stray_turns[0]/self.n_stray_turns[1]
                        """
                        for
                            if self.virtual_winding_windows[num].scheme == "hexa":
                                y = bot_bound + self.windings[num].conductor_radius
                                x = left_bound + self.windings[num].conductor_radius
                                i = 0
                                base_line = True
                                # Case n_conductors higher that "allowed" is missing
                                while x < right_bound - self.windings[num].conductor_radius and i < self.windings[num].turns:
                                    while y < top_bound - self.windings[num].conductor_radius and i < self.windings[num].turns:
                                        self.p_conductor[num].append([x, y, 0, self.c_center_conductor[num]])
                                        self.p_conductor[num].append([x - self.windings[num].conductor_radius, y, 0, self.c_conductor[num]])
                                        self.p_conductor[num].append([x, y + self.windings[num].conductor_radius, 0, self.c_conductor[num]])
                                        self.p_conductor[num].append([x + self.windings[num].conductor_radius, y, 0, self.c_conductor[num]])
                                        self.p_conductor[num].append([x, y - self.windings[num].conductor_radius, 0, self.c_conductor[num]])
                                        i += 1
                                        y += self.windings[num].conductor_radius * 2 + self.cond_cond_isolation[num]  # from bottom to top
                                    x += 2 * np.cos(np.pi / 6) * (self.windings[num].conductor_radius + self.cond_cond_isolation[
                                        num] / 2)  # * np.sqrt(2 / 3 * np.pi / np.sqrt(3))  # one step from left to right
                                    # depending on what line, hexa scheme starts shifted
                                    # reset y to "new" bottom
                                    base_line = (not base_line)
                                    if base_line:
                                        y = bot_bound + self.windings[num].conductor_radius
                                    else:
                                        y = bot_bound + 2 * self.windings[num].conductor_radius + self.cond_cond_isolation[num] / 2
                        """

                if self.virtual_winding_windows[n_win].scheme == "vertical":
                    """
                    - Vertical interleaving means a winding scheme where the two conductors are alternating in vertical 
                      (y-)direction
                    - This is practically uncommon
                    - If the turns ratio is != 1, the scheme always begins with the "higher-turns-number's" conductor
                    """

                if self.virtual_winding_windows[n_win].scheme == "horizontal":
                    #- Horizontal interleaving means a winding scheme where the two conductors are alternating in horizontal
                    #  (x-)direction (Tonnenwicklung)
                    #- This is practically most common
                    #- If the turns ratio is != 1, the scheme always begins with the "higher-turns-number's" conductor

                    # assume 2 winding transformer and dedicated stray path:
                    if (self.dedicated_stray_path == True) or (self.n_windings == 2):
                        # Initialize the list, that counts the already placed conductors
                        N_completed=[0, 0]

                        # Initialize the starting conductor
                        if self.windings[0].turns[n_win] >= self.windings[1].turns[n_win]:
                            col_cond = 0
                        else:
                            col_cond = 1

                        # Initialize the x and y coordinate
                        x = left_bound + self.windings[col_cond].conductor_radius
                        y = bot_bound + self.windings[col_cond].conductor_radius


                        # Continue placing as long as not all conductors have been placed
                        while (self.windings[0].turns[n_win]-N_completed[0] != 0) or (self.windings[1].turns[n_win]-N_completed[1] != 0):
                            #print(N_completed[col_cond])
                            if self.windings[col_cond].turns[n_win]-N_completed[col_cond] != 0:
                                # is this winding not already finished?
                                if x < right_bound - self.windings[col_cond].conductor_radius:
                                    while y < top_bound - self.windings[col_cond].conductor_radius and \
                                            N_completed[col_cond] < self.windings[col_cond].turns[n_win]:

                                        self.p_conductor[col_cond].append([x, y, 0, self.c_center_conductor[col_cond]])
                                        self.p_conductor[col_cond].append([x - self.windings[col_cond].conductor_radius, y, 0, self.c_conductor[col_cond]])
                                        self.p_conductor[col_cond].append([x, y + self.windings[col_cond].conductor_radius, 0, self.c_conductor[col_cond]])
                                        self.p_conductor[col_cond].append([x + self.windings[col_cond].conductor_radius, y, 0, self.c_conductor[col_cond]])
                                        self.p_conductor[col_cond].append([x, y - self.windings[col_cond].conductor_radius, 0, self.c_conductor[col_cond]])
                                        N_completed[col_cond] += 1
                                        #print(self.p_conductor)
                                        #print(N_completed)
                                        #print(col_cond)
                                        print(x, right_bound-self.windings[col_cond].conductor_radius)

                                        y += self.windings[col_cond].conductor_radius * 2 + self.cond_cond_isolation[col_cond]  # one from bot to top
                                        #print(y)
                                    x += self.windings[col_cond].conductor_radius + self.windings[(col_cond + 1) % 2].conductor_radius + self.cond_cond_isolation[2]  # from left to right
                                    # Reset y
                                    col_cond = (col_cond + 1) % 2
                                    y = bot_bound + self.windings[col_cond].conductor_radius
                                    #print(col_cond)
                                else:
                                    break
                            else:
                                # is this winding already finished? - continue with the other one
                                col_cond = (col_cond + 1) % 2
                                # Correct the reset of y and correct x displacement
                                x += self.windings[col_cond].conductor_radius - self.windings[(col_cond + 1) % 2].conductor_radius \
                                     - self.cond_cond_isolation[2] + self.cond_cond_isolation[col_cond]
                                y = bot_bound + self.windings[col_cond].conductor_radius
                                #print(col_cond)

                """Blockwise concentrated"""
                if isinstance(self.virtual_winding_windows[n_win].scheme, list):
                    """
                    - interleaving with a list means a concentrated winding scheme of ("hexagonal", "square" or mixed) 
                      in virtual winding window
                    - only valid for two winding case 
                    - vertical stacking
                    - block winding

                    how many turns fit in arow?
                    von top nach bot
                    while (not placed all cond.):
                        1. start with the primary winding from bot / left
                        2. continue with the secondary from top / right
                        3.CHECK solation conditions
                    """
                    # CHECK for two winding transformer
                    if len(self.virtual_winding_windows[n_win].scheme) != 2:
                        print(f"Interleaving with a list is only valid for the two winding case.\n"
                              f"Therefore the scheme must be a list of lentgh 2 but is of length "
                              f"{len(self.virtual_winding_windows[n_win].scheme)}")
                        raise Warning


                    for num in range(0, len(self.virtual_winding_windows[n_win].scheme)):

                        # Cases
                        if num == 0:
                            y = bot_bound + self.windings[num].conductor_radius
                        if num == 1:
                            y = top_bound - self.windings[num].conductor_radius

                        # Initialization
                        x = left_bound + self.windings[num].conductor_radius
                        i = 0

                        if self.virtual_winding_windows[n_win].scheme[num] == "square":
                                                       
                            # Primary winding from bottom to top
                            if num == 0:
                                while y < top_bound - self.windings[num].conductor_radius and \
                                        i < self.windings[num].turns[n_win]:
                                    while x < right_bound - self.windings[num].conductor_radius and \
                                            i < self.windings[num].turns[n_win]:
                                        self.p_conductor[num].append([x, y, 0, self.c_center_conductor[num]])
                                        self.p_conductor[num].append(
                                            [x - self.windings[num].conductor_radius, y, 0, self.c_conductor[num]])
                                        self.p_conductor[num].append(
                                            [x, y + self.windings[num].conductor_radius, 0, self.c_conductor[num]])
                                        self.p_conductor[num].append(
                                            [x + self.windings[num].conductor_radius, y, 0, self.c_conductor[num]])
                                        self.p_conductor[num].append(
                                            [x, y - self.windings[num].conductor_radius, 0, self.c_conductor[num]])
                                        i += 1
                                        x += self.windings[num].conductor_radius * 2 + self.cond_cond_isolation[
                                            num]  # from left to right
                                    y += self.windings[num].conductor_radius * 2 + self.cond_cond_isolation[
                                        num]  # one step from bot to top
                                    x = left_bound + self.windings[num].conductor_radius  # always the same

                            # Secondary winding from top to bottom
                            if num == 1:
                                while y > bot_bound + self.windings[num].conductor_radius and i < self.windings[num].turns[n_win]:
                                    while x < right_bound - self.windings[num].conductor_radius and i < self.windings[num].turns[n_win]:
                                        self.p_conductor[num].append([x, y, 0, self.c_center_conductor[num]])
                                        self.p_conductor[num].append(
                                            [x - self.windings[num].conductor_radius, y, 0, self.c_conductor[num]])
                                        self.p_conductor[num].append(
                                            [x, y + self.windings[num].conductor_radius, 0, self.c_conductor[num]])
                                        self.p_conductor[num].append(
                                            [x + self.windings[num].conductor_radius, y, 0, self.c_conductor[num]])
                                        self.p_conductor[num].append(
                                            [x, y - self.windings[num].conductor_radius, 0, self.c_conductor[num]])
                                        i += 1

                                        x += self.windings[num].conductor_radius * 2 + self.cond_cond_isolation[
                                            num]  # from left to right
                                    y += -(self.windings[num].conductor_radius * 2) - self.cond_cond_isolation[
                                        num]  # one step from bot to top
                                    x = left_bound + self.windings[num].conductor_radius  # always the same

                        if self.virtual_winding_windows[n_win].scheme[num] == "hexa":

                            # Primary winding from bottom to top
                            if num == 0:

                                base_line = True

                                while y < top_bound - self.windings[num].conductor_radius and i < self.windings[num].turns[n_win]:
                                    while x < right_bound - self.windings[num].conductor_radius and i < self.windings[num].turns[n_win]:
                                        print(f"i: {i}")
    
                                        self.p_conductor[num].append([x, y, 0, self.c_center_conductor[num]])
                                        self.p_conductor[num].append(
                                            [x - self.windings[num].conductor_radius, y, 0, self.c_conductor[num]])
                                        self.p_conductor[num].append(
                                            [x, y + self.windings[num].conductor_radius, 0, self.c_conductor[num]])
                                        self.p_conductor[num].append(
                                            [x + self.windings[num].conductor_radius, y, 0, self.c_conductor[num]])
                                        self.p_conductor[num].append(
                                            [x, y - self.windings[num].conductor_radius, 0, self.c_conductor[num]])
                                        i += 1

                                        x += 2 * np.cos(np.pi / 6) * (
                                                    self.windings[num].conductor_radius + self.cond_cond_isolation[
                                                num] / 2)
                                        # depending on what line, hexa scheme starts shifted
                                        # reset y to "new" bottom
                                        base_line = (not base_line)
                                        if base_line:
                                            y -=  (self.windings[num].conductor_radius + self.cond_cond_isolation[num])
                                        else:
                                            y +=  (self.windings[num].conductor_radius + self.cond_cond_isolation[num])

                                    # Undo last base_line reset
                                    if base_line:
                                        y +=  (self.windings[num].conductor_radius + self.cond_cond_isolation[num])
                                    else:
                                        y -=  (self.windings[num].conductor_radius + self.cond_cond_isolation[num])

                                    base_line = True
                                    x = left_bound + self.windings[num].conductor_radius
                                    y += self.windings[num].conductor_radius + self.cond_cond_isolation[num]

                            # Secondary winding from top to bottom
                            if num == 1:

                                base_line = True

                                while y > bot_bound + self.windings[num].conductor_radius and i < self.windings[num].turns[n_win]:
                                    while x < right_bound - self.windings[num].conductor_radius and i < self.windings[num].turns[
                                        n_win]:
                                        print(f"i: {i} "
                                              f"x: {x} "
                                              f"y: {y} ")

                                        self.p_conductor[num].append([x, y, 0, self.c_center_conductor[num]])
                                        self.p_conductor[num].append(
                                            [x - self.windings[num].conductor_radius, y, 0, self.c_conductor[num]])
                                        self.p_conductor[num].append(
                                            [x, y + self.windings[num].conductor_radius, 0, self.c_conductor[num]])
                                        self.p_conductor[num].append(
                                            [x + self.windings[num].conductor_radius, y, 0, self.c_conductor[num]])
                                        self.p_conductor[num].append(
                                            [x, y - self.windings[num].conductor_radius, 0, self.c_conductor[num]])
                                        i += 1


                                        x += 2 * np.cos(np.pi / 6) * (
                                                self.windings[num].conductor_radius + self.cond_cond_isolation[num] / 2)

                                        # depending on what line, hexa scheme starts shifted
                                        # reset y to "new" bottom
                                        base_line = (not base_line)
                                        if base_line:
                                            y += (self.windings[num].conductor_radius + self.cond_cond_isolation[num])
                                        else:
                                            y -= (self.windings[num].conductor_radius + self.cond_cond_isolation[num])

                                    # Undo last base_line reset
                                    if base_line:
                                        y -= (self.windings[num].conductor_radius + self.cond_cond_isolation[num])
                                    else:
                                        y += (self.windings[num].conductor_radius + self.cond_cond_isolation[num])

                                    base_line = True
                                    x = left_bound + self.windings[num].conductor_radius
                                    # from top to bottom
                                    y -= (self.windings[num].conductor_radius + self.cond_cond_isolation[num])


            else:
                # other case is non-interleaved
                if self.virtual_winding_windows[n_win].winding == "primary":
                    num = 0
                if self.virtual_winding_windows[n_win].winding == "secondary":
                    num = 1

                if self.windings[num].conductor_type == "full":
                    if sum(self.windings[num].turns) != 1:
                        print(f"For a \"full\" conductor you must choose 1 turn for each conductor!")
                    # full window conductor
                    self.p_conductor[num].append([left_bound, bot_bound, 0, self.c_conductor[num]])
                    self.p_conductor[num].append([right_bound, bot_bound, 0, self.c_conductor[num]])
                    self.p_conductor[num].append([left_bound, top_bound, 0, self.c_conductor[num]])
                    self.p_conductor[num].append([right_bound, top_bound, 0, self.c_conductor[num]])


                if self.windings[num].conductor_type == "stacked":
                    # Stack defined number of turns and chosen thickness
                    for i in range(0, self.windings[num].turns[n_win]):
                        # CHECK if top bound is reached
                        if (bot_bound + (i+1)*self.windings[num].thickness + i*self.cond_cond_isolation[num]) <= top_bound:
                            # stacking from the ground
                            self.p_conductor[num].append([left_bound, bot_bound+ i*self.windings[num].thickness + i*self.cond_cond_isolation[num], 0, self.c_conductor[num]])
                            self.p_conductor[num].append([right_bound, bot_bound + i*self.windings[num].thickness + i*self.cond_cond_isolation[num], 0, self.c_conductor[num]])
                            self.p_conductor[num].append([left_bound, bot_bound+ (i+1)*self.windings[num].thickness + i*self.cond_cond_isolation[num], 0, self.c_conductor[num]])
                            self.p_conductor[num].append([right_bound, bot_bound+ (i+1)*self.windings[num].thickness + i*self.cond_cond_isolation[num], 0, self.c_conductor[num]])


                if self.windings[num].conductor_type == "foil":
                    # Wrap defined number of turns and chosen thickness
                    if self.wrap_para[num] == "fixed_thickness":
                        for i in range(0, self.windings[num].turns[n_win]):
                            # CHECK if right bound is reached
                            if (left_bound + (i + 1) * self.windings[num].thickness + i * self.cond_cond_isolation[num]) <= right_bound:
                                # Foils
                                self.p_conductor[num].append([left_bound + i    *self.windings[num].thickness + i*self.cond_cond_isolation[num], bot_bound, 0, self.c_conductor[num]])
                                self.p_conductor[num].append([left_bound + (i+1)*self.windings[num].thickness + i*self.cond_cond_isolation[num], bot_bound, 0, self.c_conductor[num]])
                                self.p_conductor[num].append([left_bound + i    *self.windings[num].thickness + i*self.cond_cond_isolation[num], top_bound, 0, self.c_conductor[num]])
                                self.p_conductor[num].append([left_bound + (i+1)*self.windings[num].thickness + i*self.cond_cond_isolation[num], top_bound, 0, self.c_conductor[num]])


                    # Fill the allowed space in the Winding Window with a chosen number of turns
                    if self.wrap_para[num] == "interpolate":
                        x_interpol = np.linspace(left_bound, right_bound+self.cond_cond_isolation[num], self.windings[num].turns+1)
                        print(x_interpol)
                        for i in range(0, self.windings[num].turns):
                            # Foils
                            self.p_conductor[num].append([x_interpol[i], bot_bound, 0, self.c_conductor[num]])
                            self.p_conductor[num].append([x_interpol[i+1] - self.cond_cond_isolation[num], bot_bound, 0, self.c_conductor[num]])
                            self.p_conductor[num].append([x_interpol[i], top_bound, 0, self.c_conductor[num]])
                            self.p_conductor[num].append([x_interpol[i+1] - self.cond_cond_isolation[num], top_bound, 0, self.c_conductor[num]])


                # Round Conductors:
                if self.windings[num].conductor_type == "litz" or self.windings[num].conductor_type == "solid":

                    if self.virtual_winding_windows[num].scheme == "square":
                        y = bot_bound + self.windings[num].conductor_radius
                        x = left_bound + self.windings[num].conductor_radius
                        i = 0
                        # Case n_conductors higher that "allowed" is missing
                        while y < top_bound-self.windings[num].conductor_radius and i < self.windings[num].turns[n_win]:
                            while x < right_bound-self.windings[num].conductor_radius and i < self.windings[num].turns[n_win]:
                                self.p_conductor[num].append([x, y, 0, self.c_center_conductor[num]])
                                self.p_conductor[num].append([x-self.windings[num].conductor_radius, y, 0, self.c_conductor[num]])
                                self.p_conductor[num].append([x, y+self.windings[num].conductor_radius, 0, self.c_conductor[num]])
                                self.p_conductor[num].append([x+self.windings[num].conductor_radius, y, 0, self.c_conductor[num]])
                                self.p_conductor[num].append([x, y-self.windings[num].conductor_radius, 0, self.c_conductor[num]])
                                i += 1
                                x += self.windings[num].conductor_radius * 2 + self.cond_cond_isolation[num]  # from left to top
                            y += self.windings[num].conductor_radius * 2 + self.cond_cond_isolation[num]  # one step from left to right
                            x = left_bound + self.windings[num].conductor_radius  # always the same

                    if self.virtual_winding_windows[num].scheme == "hexa":
                        y = bot_bound + self.windings[num].conductor_radius
                        x = left_bound + self.windings[num].conductor_radius
                        i = 0
                        base_line = True
                        # Case n_conductors higher that "allowed" is missing
                        while x < right_bound - self.windings[num].conductor_radius and i < self.windings[num].turns[n_win]:
                            while y < top_bound - self.windings[num].conductor_radius and i < self.windings[num].turns[n_win]:
                                self.p_conductor[num].append([x, y, 0, self.c_center_conductor[num]])
                                self.p_conductor[num].append([x - self.windings[num].conductor_radius, y, 0, self.c_conductor[num]])
                                self.p_conductor[num].append([x, y + self.windings[num].conductor_radius, 0, self.c_conductor[num]])
                                self.p_conductor[num].append([x + self.windings[num].conductor_radius, y, 0, self.c_conductor[num]])
                                self.p_conductor[num].append([x, y - self.windings[num].conductor_radius, 0, self.c_conductor[num]])
                                i += 1
                                y += self.windings[num].conductor_radius * 2 + self.cond_cond_isolation[num]  # from bottom to top
                            x += 2 * np.cos(np.pi/6) * (self.windings[num].conductor_radius + self.cond_cond_isolation[num]/2) #* np.sqrt(2 / 3 * np.pi / np.sqrt(3))  # one step from left to right
                            # depending on what line, hexa scheme starts shifted
                            # reset y to "new" bottom
                            base_line = (not base_line)
                            if base_line:
                                y = bot_bound + self.windings[num].conductor_radius
                            else:
                                y = bot_bound + 2*self.windings[num].conductor_radius + self.cond_cond_isolation[num]/2


        # Checking the Conductors
        for num in range(0, self.n_windings):
            # Convert to numpy
            # Check if all Conductors could be resolved
            self.p_conductor[num] = np.asarray(self.p_conductor[num])


            #TODO:CHECKS for rect. conductors
            """ CHECK: rectangle conductors with 4 points
            if self.windings[num].conductor_type == "full" or self.windings[num].conductor_type == "stacked" or \
                    self.windings[num].conductor_type == "foil":
                if int(self.p_conductor[num].shape[0]/4) < self.windings[num].turns:
                    print("Could not resolve all conductors.")
                    # self.windings[num].turns = int(self.p_conductor[num].shape[0]/4)
                    self.valid = None
            """

            # CHECK: round conductors with 5 points
            if self.windings[num].conductor_type == "solid" or self.windings[num].conductor_type == "litz":
                if int(self.p_conductor[num].shape[0]/5) < sum(self.windings[num].turns):
                    print("Could not resolve all conductors.")
                    # self.windings[num].turns = int(self.p_conductor[num].shape[0]/5)
                    # TODO: break and warning
                    self.valid = None

        # Region for Boundary Condition
        self.p_region_bound[0][:] = [-r_outer*self.padding, -(self.window_h / 2 + self.core_w/4)*self.padding, 0, self.c_core*self.padding]
        self.p_region_bound[1][:] = [r_outer*self.padding, -(self.window_h / 2 + self.core_w/4)*self.padding, 0, self.c_core*self.padding]
        self.p_region_bound[2][:] = [-r_outer*self.padding, (self.window_h / 2 + self.core_w/4)*self.padding, 0, self.c_core*self.padding]
        self.p_region_bound[3][:] = [r_outer*self.padding, (self.window_h / 2 + self.core_w/4)*self.padding, 0, self.c_core*self.padding]

    # === Meshing ===
    def generate_mesh(self, refine=0, alternative_error=0):
        """
        - interaction with gmsh
        - mesh generation
            - Skin depth based forward meshing
            [- adaptive refinement
                - with the help of mesh-size-fields/background meshes
                - with an appropriate local error metric ]
        :return:
        """
        # Initialization
        gmsh.initialize()

        if refine == 1:
            # Choose applied Error Function
            if alternative_error == 1:
                local_error = self.alternative_local_error()  # something like current density
            else:
                local_error = self.local_error()  # Here a "real" numeric error should be applied

            self.create_background_mesh(local_error)

        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.model.add(self.path_mesh + "geometry")
        # ------------------------------------------ Geometry -------------------------------------------
        # Core generation
        if self.core_type == "EI":
            # --------------------------------------- Points --------------------------------------------
            if self.y_symmetric == 1:
                if self.dimensionality == "2D axi":

                    # Find points of air gaps (used later)
                    if self.n_air_gaps > 0:
                        # Top and bottom point
                        center_right = min_max_inner_points(self.p_window[4], self.p_window[6], self.p_air_gaps)
                        island_right = inner_points(self.p_window[4], self.p_window[6], self.p_air_gaps)
                        # print(f"Air gap island points: {center_right}\n{island_right}")
                        # Explicite stray path:
                        if self.dedicated_stray_path == True:
                            # mshopt stray_path_gap = [[], []]
                            # mshopt stray_path_gap[0][:] = island_right[(self.start_i-1)*2][:]
                            # mshopt stray_path_gap[1][:] = island_right[(self.start_i-1)*2+1][:]
                            island_right[(self.start_i-1)*2][0] = self.added_bar
                            island_right[(self.start_i-1)*2+1][0] = self.added_bar


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
                    curve_loop_outer_air = []
                    curve_loop_bound = []
                    # Plane Surfaces
                    plane_surface_core = []
                    plane_surface_cond = [[], []]
                    plane_surface_air = []
                    plane_surface_outer_air = []

                    # =====================
                    # Main Core

                    # """ Points """
                    # (index refers to sketch)

                    # First point (left point of lowest air gap)
                    if self.n_air_gaps > 0:
                        #p_core.append(gmsh.model.geo.addPoint(0, center_right[0][1], center_right[0][2],
                        #                                        center_right[0][3]*self.zero_mesh_acc))
                        p_core.append(gmsh.model.geo.addPoint(0, center_right[0][1], center_right[0][2],
                                                                self.c_core))
                    if self.n_air_gaps == 0:
                        p_core.append(None)  # dummy filled for no air gap special case

                    # Go down and counter-clockwise
                    # Four points around the core
                    p_core.append(gmsh.model.geo.addPoint(0, self.p_outer[1][1], self.p_outer[1][2],
                                                          self.p_outer[1][3]))
                    p_core.append(gmsh.model.geo.addPoint(self.p_outer[1][0], self.p_outer[1][1], self.p_outer[1][2],
                                                          self.p_outer[1][3]))
                    p_core.append(gmsh.model.geo.addPoint(self.p_outer[3][0], self.p_outer[3][1], self.p_outer[3][2],
                                                          self.p_outer[3][3]))
                    p_core.append(gmsh.model.geo.addPoint(0, self.p_outer[3][1], self.p_outer[3][2],
                                                          self.p_outer[3][3]))

                    # Two points of highest air gap
                    if self.n_air_gaps > 0:
                        #p_core.append(gmsh.model.geo.addPoint(0, center_right[1][1], center_right[1][2],
                        #                                        center_right[1][3]*self.zero_mesh_acc))
                        p_core.append(gmsh.model.geo.addPoint(0, center_right[1][1], center_right[1][2],
                                                                self.c_core))
                        #p_core.append(gmsh.model.geo.addPoint(center_right[1][0], center_right[1][1],
                        #                                      center_right[1][2], center_right[1][3]))
                        p_core.append(gmsh.model.geo.addPoint(center_right[1][0], center_right[1][1],
                                                              center_right[1][2], self.c_window))
                    if self.n_air_gaps == 0:
                        p_core.append(None)  # dummy filled for no air gap special case
                        p_core.append(None)  # dummy filled for no air gap special case

                    # Clockwise
                    # Four points of window
                    p_core.append(gmsh.model.geo.addPoint(self.p_window[6][0], self.p_window[6][1], self.p_window[6][2],
                                                          self.p_window[6][3]))
                    p_core.append(gmsh.model.geo.addPoint(self.p_window[7][0], self.p_window[7][1], self.p_window[7][2],
                                                          self.p_window[7][3]))
                    p_core.append(gmsh.model.geo.addPoint(self.p_window[5][0], self.p_window[5][1], self.p_window[5][2],
                                                          self.p_window[5][3]))
                    p_core.append(gmsh.model.geo.addPoint(self.p_window[4][0], self.p_window[4][1], self.p_window[4][2],
                                                          self.p_window[4][3]))

                    # Last point of lowest air gap
                    if self.n_air_gaps > 0:
                        p_core.append(gmsh.model.geo.addPoint(center_right[0][0], center_right[0][1],
                                                              center_right[0][2], self.c_window))
                    if self.n_air_gaps == 0:
                        p_core.append(None)  # dummy filled for no air gap special case


                    # """ Curves """
                    # (index refers to sketch)
                    # To be added: Case Air Gaps directly on outer leg
                    # Curves: Boundary - Core
                    if self.n_air_gaps > 0:
                        l_bound_core.append(gmsh.model.geo.addLine(p_core[0], p_core[1]))
                    l_bound_core.append(gmsh.model.geo.addLine(p_core[1], p_core[2]))
                    l_bound_core.append(gmsh.model.geo.addLine(p_core[2], p_core[3]))
                    l_bound_core.append(gmsh.model.geo.addLine(p_core[3], p_core[4]))
                    if self.n_air_gaps > 0:
                        l_bound_core.append(gmsh.model.geo.addLine(p_core[4], p_core[5]))
                    if self.n_air_gaps == 0:
                        l_bound_core.append(gmsh.model.geo.addLine(p_core[4], p_core[1]))
                    # Curves: Core - Air
                    if self.n_air_gaps > 0:
                        l_core_air.append(gmsh.model.geo.addLine(p_core[5], p_core[6]))
                        l_core_air.append(gmsh.model.geo.addLine(p_core[6], p_core[7]))
                    l_core_air.append(gmsh.model.geo.addLine(p_core[7], p_core[8]))
                    l_core_air.append(gmsh.model.geo.addLine(p_core[8], p_core[9]))
                    l_core_air.append(gmsh.model.geo.addLine(p_core[9], p_core[10]))
                    if self.n_air_gaps > 0:
                        l_core_air.append(gmsh.model.geo.addLine(p_core[10], p_core[11]))
                        l_core_air.append(gmsh.model.geo.addLine(p_core[11], p_core[0]))
                    if self.n_air_gaps == 0:
                        l_core_air.append(gmsh.model.geo.addLine(p_core[10], p_core[7]))
                    # Plane: Main Core --> plane_surface_core[0]
                    if self.n_air_gaps > 0:
                        curve_loop_core = gmsh.model.geo.addCurveLoop(l_bound_core + l_core_air)
                        plane_surface_core.append(gmsh.model.geo.addPlaneSurface([curve_loop_core]))
                    if self.n_air_gaps == 0:
                        curve_loop_bound_core = gmsh.model.geo.addCurveLoop(l_bound_core)
                        curve_loop_core_air = gmsh.model.geo.addCurveLoop(l_core_air)
                        plane_surface_core.append(
                            gmsh.model.geo.addPlaneSurface([curve_loop_bound_core, curve_loop_core_air]))
                    # =====================
                    # Core Islands (Core parts between Air Gaps)
                    # Points of Core Islands (index refers to sketch)
                    if self.n_air_gaps != 0:
                        while island_right.shape[0] > 0:
                            # take two points with lowest y-coordinates
                            min_island_right = np.argmin(island_right[:, 1])
                            #p_island.append(gmsh.model.geo.addPoint(0, island_right[min_island_right, 1],
                            #                                        island_right[min_island_right, 2],
                            #                                        island_right[min_island_right, 3]*self.zero_mesh_acc))
                            #p_island.append(gmsh.model.geo.addPoint(island_right[min_island_right, 0],
                            #                                        island_right[min_island_right, 1],
                            #                                        island_right[min_island_right, 2],
                            #                                        island_right[min_island_right, 3]))
                            p_island.append(gmsh.model.geo.addPoint(0, island_right[min_island_right, 1],
                                                                    island_right[min_island_right, 2],
                                                                    self.c_core))
                            p_island.append(gmsh.model.geo.addPoint(island_right[min_island_right, 0],
                                                                    island_right[min_island_right, 1],
                                                                    island_right[min_island_right, 2],
                                                                    self.c_window))
                            print(f"isl: {island_right[min_island_right, 3]}")
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
                    if self.n_air_gaps == 1:
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
                    # =====================
                    # Conductors
                    # Points of Conductors
                    for num in range(0, self.n_windings):
                        for i in range(0, self.p_conductor[num].shape[0]):
                            """
                            print(f"x = {np.round(self.p_conductor[num][i][0], decimals=3)}  "
                                  f"y = {np.round(self.p_conductor[num][i][1], decimals=3)}  "
                                  f"z = 0  "
                                  f"mesh_accuracy = {np.round(self.p_conductor[num][i][3], decimals=4)}  \n"
                                  )
                            """
                            p_cond[num].append(gmsh.model.geo.addPoint(self.p_conductor[num][i][0],  # x
                                                                       self.p_conductor[num][i][1],  # y
                                                                       0,                            # z
                                                                       self.p_conductor[num][i][3])) # mesh_accuracy
                        # Curves of Conductors
                        if self.windings[num].conductor_type == "litz" or self.windings[num].conductor_type == "solid":
                            for i in range(0, int(len(p_cond[num]) / 5)):
                                l_cond[num].append(gmsh.model.geo.addCircleArc(
                                    p_cond[num][5 * i + 1], p_cond[num][5 * i + 0], p_cond[num][5 * i + 2]))
                                l_cond[num].append(gmsh.model.geo.addCircleArc(
                                    p_cond[num][5 * i + 2], p_cond[num][5 * i + 0], p_cond[num][5 * i + 3]))
                                l_cond[num].append(gmsh.model.geo.addCircleArc(
                                    p_cond[num][5 * i + 3], p_cond[num][5 * i + 0], p_cond[num][5 * i + 4]))
                                l_cond[num].append(gmsh.model.geo.addCircleArc(
                                    p_cond[num][5 * i + 4], p_cond[num][5 * i + 0], p_cond[num][5 * i + 1]))
                                # Iterative plane creation
                                curve_loop_cond[num].append(gmsh.model.geo.addCurveLoop(
                                    [l_cond[num][i * 4 + 0], l_cond[num][i * 4 + 1], l_cond[num][i * 4 + 2], l_cond[num][i * 4 + 3]]))
                                plane_surface_cond[num].append(gmsh.model.geo.addPlaneSurface([curve_loop_cond[num][i]]))
                        else:
                            # Rectangle conductor cut
                            for i in range(0, int(len(p_cond[num]) / 4)):
                                l_cond[num].append(gmsh.model.geo.addLine(p_cond[num][4 * i + 0], p_cond[num][4 * i + 2]))
                                l_cond[num].append(gmsh.model.geo.addLine(p_cond[num][4 * i + 2], p_cond[num][4 * i + 3]))
                                l_cond[num].append(gmsh.model.geo.addLine(p_cond[num][4 * i + 3], p_cond[num][4 * i + 1]))
                                l_cond[num].append(gmsh.model.geo.addLine(p_cond[num][4 * i + 1], p_cond[num][4 * i + 0]))
                                # Iterative plane creation
                                curve_loop_cond[num].append(gmsh.model.geo.addCurveLoop(
                                    [l_cond[num][i * 4 + 0], l_cond[num][i * 4 + 1], l_cond[num][i * 4 + 2], l_cond[num][i * 4 + 3]]))
                                plane_surface_cond[num].append(gmsh.model.geo.addPlaneSurface([curve_loop_cond[num][i]]))
                    print(f"plane_surface_cond {plane_surface_cond}")
                    # =====================

                    # Air (Points are partwise double designated)
                    l_air_tmp = l_core_air[:7]
                    for i in range(0, len(l_bound_air)):
                        l_air_tmp.append(l_bound_air[i])
                        if i < len(l_bound_air) - 1:
                            l_air_tmp.append(l_core_air[7 + 3 * i])
                            l_air_tmp.append(l_core_air[7 + 3 * i + 1])
                            l_air_tmp.append(l_core_air[7 + 3 * i + 2])
                    curve_loop_air.append(gmsh.model.geo.addCurveLoop(l_air_tmp))
                    flatten_curve_loop_cond = [j for sub in curve_loop_cond for j in sub]
                    plane_surface_air.append(gmsh.model.geo.addPlaneSurface(curve_loop_air + flatten_curve_loop_cond))  # xfmr

                    # =====================
                    # Bound
                    if self.region == None:
                        l_bound_tmp = l_bound_core[:5]
                        for i in range(0, len(l_bound_air)):
                            l_bound_tmp.append(l_bound_air[-i - 1])
                            if i != len(l_bound_air) - 1:  # last run
                                l_bound_tmp.append(l_bound_core[-i - 1])
                    else:
                        # Generate Lines of Region
                        # start top left and go clockwise
                        p_region.append(gmsh.model.geo.addPoint(0, self.p_region_bound[2][1],
                                                                self.p_region_bound[2][2], self.p_region_bound[2][3]))
                        p_region.append(gmsh.model.geo.addPoint(self.p_region_bound[3][0], self.p_region_bound[3][1],
                                                                self.p_region_bound[3][2], self.p_region_bound[3][3]))
                        p_region.append(gmsh.model.geo.addPoint(self.p_region_bound[1][0], self.p_region_bound[1][1],
                                                                self.p_region_bound[1][2], self.p_region_bound[1][3]))
                        p_region.append(gmsh.model.geo.addPoint(0, self.p_region_bound[0][1],
                                                                self.p_region_bound[0][2], self.p_region_bound[0][3]))

                        # Outer Region Lines
                        l_region.append(gmsh.model.geo.addLine(p_core[4], p_region[0]))
                        l_region.append(gmsh.model.geo.addLine(p_region[0], p_region[1]))
                        l_region.append(gmsh.model.geo.addLine(p_region[1], p_region[2]))
                        l_region.append(gmsh.model.geo.addLine(p_region[2], p_region[3]))
                        l_region.append(gmsh.model.geo.addLine(p_region[3], p_core[1]))

                        # Boundary Line
                        l_bound_tmp = []
                        l_bound_tmp.append(l_bound_core[4])

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

                    #print(l_bound_tmp)
                    #curve_loop_bound.append(gmsh.model.geo.addCurveLoop(l_bound_tmp, reorient=True))

        # Define physical Surfaces and Curves
        # Core
        ps_core = gmsh.model.geo.addPhysicalGroup(2, plane_surface_core, tag=2000)
        # Conductors
        ps_cond = [[], []]  # xfmr
        for num in range(0, self.n_windings):

            if self.windings[num].conductor_type == "foil" or self.windings[num].conductor_type == "solid" or \
                    self.windings[num].conductor_type == "full" or self.windings[num].conductor_type == "stacked":
                for i in range(0, sum(self.windings[num].turns)):
                    print(plane_surface_cond)
                    ps_cond[num].append(
                        gmsh.model.geo.addPhysicalGroup(2, [plane_surface_cond[num][i]], tag=4000 + 1000 * num + i))
            if self.windings[num].conductor_type == "litz":
                for i in range(0, sum(self.windings[num].turns)):
                    ps_cond[num].append(
                        gmsh.model.geo.addPhysicalGroup(2, [plane_surface_cond[num][i]], tag=6000 + 1000 * num + i))


        # Air
        ps_air = gmsh.model.geo.addPhysicalGroup(2, plane_surface_air, tag=1000)
        ps_air_ext = gmsh.model.geo.addPhysicalGroup(2, plane_surface_outer_air, tag=1001)
        # Boundary
        pc_bound = gmsh.model.geo.addPhysicalGroup(1, l_bound_tmp, tag=1111)
        #print(f"Physical Conductor Surfaces: {ps_cond}")  # xfmr
        gmsh.model.setPhysicalName(2, ps_core, "CORE")
        for num in range(0, self.n_windings):
            for i in range(0, len(ps_cond[num])):
                gmsh.model.setPhysicalName(2, ps_cond[num][i], f"COND{num+1}")
        gmsh.model.setPhysicalName(2, ps_air, "AIR")
        gmsh.model.setPhysicalName(1, pc_bound, "BOUND")

        # Remove Points from Model
        #for i in range(9,39):
        #    gmsh.model.geo.remove(dimTags=[(0, i)])
        #print(f"P_cond: {p_cond}")


        # - Forward Meshing -
        # Inter Conductors
        for n_win in range(0, len(self.virtual_winding_windows)):
            if self.virtual_winding_windows[n_win].winding != "interleaved":
                for num in range(0, self.n_windings):
                    p_inter = []
                    x_inter = []
                    y_inter = []
                    j = 0

                    if self.windings[num].conductor_type == "solid" and self.windings[num].turns[n_win] > 1:
                        while self.p_conductor[num][5*j][1] == self.p_conductor[num][5*j+5][1]:
                            x_inter.append(0.5*(self.p_conductor[num][5*j][0]+self.p_conductor[num][5*j+5][0]))
                            j += 1
                            if j == self.windings[num].turns[n_win]-1:
                                break
                        j += 1
                        #print(f"j = {j}")
                        if int(self.windings[num].turns[n_win]/j) > 1:
                            for i in range(0, int(self.windings[num].turns[n_win]/j)):
                                if 5*j*i+5*j >= len(self.p_conductor[num][:]):
                                    break
                                y_inter.append(0.5*(self.p_conductor[num][5*j*i][1]+self.p_conductor[num][5*j*i+5*j][1]))
                            for x in x_inter:
                                for y in y_inter:
                                    p_inter.append(gmsh.model.geo.addPoint(x, y, 0, self.c_center_conductor[num]))
                    #print(f"x_inter = {x_inter}")
                    #print(f"y_inter = {y_inter}")
                    #print(f"p_inter = {p_inter}")



        # mshopt # Explicite stray path air gap optimization
        # mshopt if self.dedicated_stray_path != None:
        # mshopt     stray_path_mesh_optimizer = []
            # mshopt     stray_path_mesh_optimizer.append(gmsh.model.geo.addPoint(stray_path_gap[0][0], stray_path_gap[0][1]+0.0001, stray_path_gap[0][2], 0.5*stray_path_gap[0][3]))
            # mshopt     stray_path_mesh_optimizer.append(gmsh.model.geo.addPoint(stray_path_gap[1][0], stray_path_gap[1][1]-0.0001, stray_path_gap[1][2], 0.5*stray_path_gap[1][3]))
            # mshopt     print(f"plane_surface_core: {plane_surface_core}")
            # mshopt     print(f"stray_path_mesh_optimizer: {stray_path_mesh_optimizer}")
            # mshopt     print(f"stray_path_mesh_optimizer coordinates: {stray_path_gap[0][0], stray_path_gap[0][1], stray_path_gap[0][2], stray_path_gap[0][3]}\n"
            # mshopt           f"{stray_path_gap[1][0], stray_path_gap[1][1], stray_path_gap[1][2], stray_path_gap[1][3]}")


        # Synchronize
        gmsh.model.geo.synchronize()

        # Conductor Center
        for num in range(0, self.n_windings):
            for i in range(0, int(len(p_cond[num]) / 5)):
                gmsh.model.mesh.embed(0, [p_cond[num][5 * i + 0]], 2, plane_surface_cond[num][i])

        # Embed ponts for mesh refinement
        # Inter Conductors
        for n_win in range(0, len(self.virtual_winding_windows)):
            if self.virtual_winding_windows[n_win].winding != "interleaved":
                gmsh.model.mesh.embed(0, p_inter, 2, plane_surface_air[0])
        # Stray path
        # mshopt gmsh.model.mesh.embed(0, stray_path_mesh_optimizer, 2, plane_surface_core[2])


        # Synchronize again
        gmsh.model.geo.synchronize()

        # Output .msh file
        gmsh.option.setNumber("Mesh.SaveAll", 0)
        gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        gmsh.option.setNumber("Mesh.SurfaceFaces", 1)

        # Colors
        for i in range(0, len(plane_surface_core)):
            gmsh.model.setColor([(2, plane_surface_core[i])], 50, 50, 50)
        gmsh.model.setColor([(2, plane_surface_air[0])], 0, 0, 230)
        for num in range(0, self.n_windings):
            for i in range(0, len(plane_surface_cond[num])):
                gmsh.model.setColor([(2, plane_surface_cond[num][i])], 150, 150, 0)
        # -----------------------------------------
        if refine == 1:
            print("\n ------- \nRefined Mesh Creation ")
            # mesh the new gmsh.model using the size field
            bg_field = gmsh.model.mesh.field.add("PostView")
            # TODO: gmsh.model.mesh.field.setNumber(bg_field, "ViewTag", sf_view)
            gmsh.model.mesh.field.setAsBackgroundMesh(bg_field)
            print("\nMeshing...\n")
            gmsh.model.mesh.generate(2)
            gmsh.write(self.path_mesh + "geometry.msh")
        else:
            print("\nMeshing...\n")
            gmsh.model.mesh.generate(2)

        # Mesh direction
        if not os.path.isdir(self.path + "/" + self.path_mesh):
            os.mkdir(self.path + "/" + self.path_mesh)

        # Check operating system
        if sys.platform == "linux" or sys.platform == "linux2":
            gmsh.write(self.path + "/" + self.path_mesh + "geometry.msh")
        elif sys.platform == "darwin":  # OS X
            gmsh.write(self.path + "/" + self.path_mesh + "geometry.msh")
        elif sys.platform == "win32":
            gmsh.write(self.path + "/" + self.path_mesh + "geometry.msh")  # Win10 can handle slash

        # Open gmsh GUI for visualization
        # gmsh.fltk.run()

        # Terminate gmsh
        gmsh.finalize()

    def mesh(self, frequency=None, skin_mesh_factor=1):
        self.high_level_geo_gen(frequency=frequency, skin_mesh_factor=skin_mesh_factor)
        if self.valid:
            self.generate_mesh()

    # == Adaptive Meshing ==
    # !!! Experimental Code !!!
    def triangle_max_edge(self, x):
        a = np.sum((x[:, 0, :] - x[:, 1, :]) ** 2, 1) ** 0.5
        b = np.sum((x[:, 0, :] - x[:, 2, :]) ** 2, 1) ** 0.5
        c = np.sum((x[:, 1, :] - x[:, 2, :]) ** 2, 1) ** 0.5
        return np.maximum(a, np.maximum(b, c))

    def triangle_mean_edge(self, x):
        a = np.sum((x[:, 0, :] - x[:, 1, :]) ** 2, 1) ** 0.5
        b = np.sum((x[:, 0, :] - x[:, 2, :]) ** 2, 1) ** 0.5
        c = np.sum((x[:, 1, :] - x[:, 2, :]) ** 2, 1) ** 0.5
        return (a + b + c)/3

    def compute_size_field_old(self, nodes, triangles, err, N):
        x = nodes[triangles]
        a = 2.
        d = 2.
        fact = (a ** ((2. + a) / (1. + a)) + a ** (1. / (1. + a))) * np.sum(err ** (2. / (1. + a)))
        ri = err ** (2. / (2. * (1 + a))) * a ** (1. / (d * (1. + a))) * ((1. + a) * N / fact) ** (1. / d)
        return self.triangle_max_edge(x) / ri

    def compute_size_field(self, nodes, triangles, err, N):
        x = nodes[triangles]
        threshold = 0.01
        print(f"Error{err}")
        err[err > 0.1] = 0.5
        err[err < 0.1] = 0.9
        print(f"Error{err}")
        return self.triangle_max_edge(x) * err
        # return self.triangle_max_edge(x) - err**(0.1) * self.triangle_max_edge(x)

    class Mesh:
        """
        Currently unused except for experimental adaptive Meshing.
        #TODO: Make the mesh an object for increased reusability
        """
        def __init__(self):
            self.vtags, vxyz, _ = gmsh.model.mesh.getNodes()
            self.vxyz = vxyz.reshape((-1, 3))
            vmap = dict({j: i for i, j in enumerate(self.vtags)})
            self.triangles_tags, evtags = gmsh.model.mesh.getElementsByType(2)
            evid = np.array([vmap[j] for j in evtags])
            self.triangles = evid.reshape((self.triangles_tags.shape[-1], -1))

    def refine_mesh(self, local=0):
        """

        :return:
        """

        # --------------------------------------
        if local == 1:
            self.generate_mesh(refine=1)

        # --------------------------------------
        if local == 0:
            # Refine current mesh
            gmsh.model.mesh.refine()
            print("\nMeshing...\n")
            gmsh.model.mesh.generate(2)
            # --------------------------------------
            # Mesh generation
            #gmsh.model.mesh.generate(2)
            # Check operating system
            if sys.platform == "linux" or sys.platform == "linux2":
                gmsh.write(self.path + "/" + self.path_mesh + "geometry.msh")
            elif sys.platform == "darwin":
                # OS X
                gmsh.write(self.path + "/" + self.path_mesh + "geometry.msh")
            elif sys.platform == "win32":
                gmsh.write(self.path + "/" + self.path_mesh + "geometry.msh")  # Win10 can handle slash

            # Terminate gmsh
            gmsh.finalize()

    def find_neighbours(self, file="results/J_rms.pos"):
        # Open loss/error results
        dest_file = open(self.path + "error.dat", "w")
        # Read the logged losses corresponding to the frequencies
        with open(self.path + "/" + self.path_res_fields + 'J_rms.pos') as f:
            read = 0
            for line in f:
                if line == "$ElementNodeData\n":
                    read = (read + 1) % 2
                    print(line, read)
                if read == 1 and ' ' in line:
                    # words = line.split(sep=' ')
                    words = line.replace(' ', ', ')
                    for word in words:
                        dest_file.write(word)
        dest_file.close()

    def alternative_local_error(self, loss_file='/results/fields/J_rms.pos'):
        """

        :return:
        """
        # Open loss/error results
        error_file = open(self.path + "/mesh_error.dat", "w")
        # Read the logged losses corresponding to the frequencies
        with open(self.path + loss_file) as f:
            read = 0
            for line in f:
                if line == "$ElementNodeData\n":
                    read = (read + 1) % 2
                    print(line, read)
                if read == 1 and ' ' in line:
                    # words = line.split(sep=' ')
                    words = line.replace(' ', ', ')
                    for word in words:
                        error_file.write(word)
        error_file.close()

        # Load local Error values
        data = pd.read_csv(self.path + "/mesh_error.dat")
        local_error = data.iloc[:, 2].to_numpy()
        local_error = np.insert(local_error, 0, 0., axis=0)

        local_error = local_error / np.max(local_error) + 0.001
        print(f"Längen: {len(local_error), local_error} ")  # first of the 3 loss values

        # ---- Open post-processing results ----

        # Elements ----
        elements = []
        values = []
        # error_file = open(self.path + "/mesh_error.dat", "w")
        # Read the logged losses corresponding to the frequencies
        with open(self.path + "/" + self.path_res_fields + 'J_rms.pos') as file:
            read = 0
            for line in file:
                if line == "$Elements\n" or line == "$EndElements\n":
                    read = (read + 1) % 2
                    print(line, read)
                if read == 1 and ' ' in line:
                    words = line.split(sep=' ')
                    elements.append(words)

        # Convert Elements to Dataframe
        element_frame = pd.DataFrame(elements)
        # Dropout not needed columns
        element_frame.drop(element_frame.columns[[1, 2, 3, 4, 8]], axis=1, inplace=True)
        element_frame.columns = ['NumElement', 'Node1', 'Node2', 'Node3']
        print(f"Number of Elements: {len(elements)}\n"
              # f"Elements: {elements}\n"
              f"Element Dataframe: {element_frame}")
        element_frame.to_csv(path_or_buf="elements.txt")

        # Values ----
        with open(self.path + "/" + self.path_res_fields + 'J_rms.pos') as file:
            read = 0
            for line in file:
                if line == "$ElementNodeData\n":
                    read = (read + 1) % 2
                    print(line, read)
                if read == 1 and ' ' in line:
                    words = re.split(' |\n', line)
                    values.append(words)

        # Convert Values to Dataframe
        value_frame = pd.DataFrame(values)
        # Dropout not needed columns
        value_frame.drop(value_frame.columns[[1, 5]], axis=1, inplace=True)
        value_frame.columns = ['NumElement', 'Node1', 'Node2', 'Node3']
        # value_frame['NumElement', 'Node1', 'Node2', 'Node3'] = pd.to_numeric(value_frame['NumElement', 'Node1', 'Node2', 'Node3'], downcast="float")
        value_frame['NumElement'] = pd.to_numeric(value_frame['NumElement'], downcast="float")
        value_frame['Node1'] = pd.to_numeric(value_frame['Node1'], downcast="float")
        value_frame['Node2'] = pd.to_numeric(value_frame['Node2'], downcast="float")
        value_frame['Node3'] = pd.to_numeric(value_frame['Node3'], downcast="float")
        print(f"Number of Values: {len(values)}\n"
              # f"Values: {values}\n"
              f"Values Dataframe: {value_frame}")

        # ---- Neighbour algorithm || Error calculation ----
        local_error = np.zeros(len(element_frame.index))

        nodes = ['Node1', 'Node2', 'Node3']
        for i in value_frame.index:
            mean_cell = 0
            # Mean loss per cell
            for node in nodes:
                mean_cell += value_frame[node][i] / len(nodes)
            # Local Variance
            for node in nodes:
                local_error[i] += (value_frame[node][i] - mean_cell) ** 2

        # Distribute Local Error on neighbour cells
        nodes_neighbours_found = []
        distribution_factor = 2

        """
        local_error_copy = local_error
        print(local_error)

        # every element
        for i in element_frame.index:
            # every element's node
            for node in nodes:
                # Node already considered?
                if not element_frame[node][i] in nodes_neighbours_found:
                    nodes_neighbours_found += element_frame[node][i]
                    # Value[Node] == 0 ? : skip
                    if not value_frame[node][i] == 0:
                        # search in every element for Node
                        for j in element_frame.index:
                            for node in nodes:
                                if value_frame[node][j] == 0:
                                    # add some error in neighboured Nodes
                                    if element_frame[node][j] == element_frame[node][i]:
                                        local_error_copy[j] += local_error[i] * distribution_factor


        local_error = local_error_copy
        """
        print(local_error)
        # Error Normalization
        return local_error / np.max(local_error) + 0.0001

        # print(f"Local Error: {local_error[3387]}\n"
        #      f"Length of Local Error: {len(local_error)}")
        """
        # Load local Error values
        data = pd.read_csv(self.path + "/mesh_error.dat")
        local_error = data.iloc[:, 2].to_numpy()
        local_error = np.insert(local_error, 0, 0., axis=0)

        local_error = local_error/np.max(local_error) + 0.001
        print(len(local_error), local_error)  # first of the 3 loss values
        """

    def local_error(self, loss_file='/results/fields/error.pos'):
        """
        - Method shall return the normalized numeric local error of the last adaptive simulation step
        - Local error can be used to optimize the mesh in the next iteration step
        :return:
        """
        # Open loss/error results
        error_file = open(self.path + "/mesh_error.dat", "w")
        # Read the logged losses corresponding to the frequencies
        with open(self.path + loss_file) as f:
            read = 0
            for line in f:
                if line == "$ElementNodeData\n":
                    read = (read + 1) % 2
                    print(line, read)
                if read == 1 and ' ' in line:
                    # words = line.split(sep=' ')
                    words = line.replace(' ', ', ')
                    for word in words:
                        error_file.write(word)
        error_file.close()

        # Load local Error values
        data = pd.read_csv(self.path + "/mesh_error.dat")
        print(f"Data: {data}")
        local_error = data.iloc[:, 2].to_numpy()
        local_error = np.insert(local_error, 0, 1e-5, axis=0)

        return local_error / np.max(local_error)
        print(f"Längen: {len(local_error), local_error} ")

    def create_background_mesh(self, local_error):
        gmsh.open(self.path + "/" + self.path_mesh + "geometry.msh")  # Open current mesh
        N = 50000  # Number of elements after remeshing
        mesh = self.Mesh()  # Create virtual mesh
        """
        print(f"Mesh nodes: {mesh.vxyz} \n "
              f"Mesh nodes.shape: {mesh.vxyz.shape} \n "
              f"Mesh node tags: {mesh.vtags} \n"
              f"Mesh triangles: {mesh.triangles} \n"
              f"Mesh triangles.shape: {mesh.triangles.shape} \n"
              f"Mesh triangle tags: {mesh.triangles_tags} \n"
              )
        """  # some printing options

        err_view = gmsh.view.add("element-wise error")
        gmsh.view.addModelData(err_view, 0, self.path_mesh + "/" + self.path_mesh + "geometry", "ElementData", mesh.triangles_tags, local_error[:, None])
        gmsh.view.write(err_view, "err.pos")

        # Refinement
        sf_ele = self.compute_size_field(mesh.vxyz, mesh.triangles, local_error, N)
        np.savetxt("sf_ele.txt", sf_ele)
        sf_view = gmsh.view.add("mesh size field")
        gmsh.view.addModelData(sf_view, 0, self.path_mesh + "/" + self.path_mesh + "geometry", "ElementData", mesh.triangles_tags, sf_ele[:, None])
        gmsh.view.write(sf_view, "sf.pos")

    # === GetDP Interaction / Simulation / Exciation ===
    def excitation(self, f, i, phases=[], ex_type='current', imposed_red_f=0):
        """
        - excitation of the electromagnetic problem
        - current, voltage or current density
        - frequency or reduced frequency
        :param phases:
        :param f:
        :param i:
        :param nonlinear:
        :param ex_type:
        :param imposed_red_f:
        :return:
        """
        print(f"\n---\n"
              f"Excitation: \n"
              f"Frequency: {f}\n"
              f"Current(s): {i}\n"
              f"Phase(s): {phases}\n")

        # -- Excitation --
        self.flag_imposed_reduced_frequency = imposed_red_f  # if == 0 --> impose frequency f
        self.flag_excitation_type = ex_type  # 'current', 'current_density', 'voltage'

        phases = np.asarray(phases)
        for num in range(0, self.n_windings):
            # Imposed current, current density or voltage
            if self.flag_excitation_type == 'current':
                self.current[num] = i[num]
                if len(phases) != 0:
                    self.phase_tmp = phases/180
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
                    self.delta = np.sqrt(2 / (2 * self.frequency * np.pi * self.cond_sigma * self.mu0))

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
        text_file.write(f"DirResCirc = \"{self.path_res_circuit}\";\n")

        # Magnetic Component Type
        if self.component_type == 'inductor':
            text_file.write(f"Flag_Transformer = 0;\n")
        if self.component_type == 'transformer':
            text_file.write(f"Flag_Transformer = 1;\n")

        # Frequency
        text_file.write("Freq = %s;\n" % self.frequency)
        text_file.write(f"delta = {self.delta};\n")

        # Core Loss
        text_file.write(f"Flag_Core_Loss = {self.core_loss_simlation};\n")
        text_file.write(f"e_r_imag = {self.core_im_epsilon_rel};\n")
        text_file.write(f"mu_r_imag = {self.core_im_mu_rel};\n")

        if self.core_loss_simlation:
            text_file.write(f"ki = {self.ki};\n")
            text_file.write(f"alpha = {self.alpha};\n")
            text_file.write(f"beta = {self.beta};\n")
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
            if self.windings[num].conductor_type == 'litz':  # xfmr
                text_file.write(f"Flag_HomogenisedModel{num+1} = 1;\n")
            else:
                text_file.write(f"Flag_HomogenisedModel{num+1} = 0;\n")
            text_file.write("Flag_imposedRr = %s;\n" % self.flag_imposed_reduced_frequency)

            # -- Geometry --
            # Number of conductors
            text_file.write(f"NbrCond{num+1} = {sum(self.windings[num].turns)};\n")

            # For stranded Conductors:
            # text_file.write(f"NbrstrandedCond = {self.turns};\n")  # redundant
            if self.windings[num].conductor_type == "litz":
                text_file.write(f"NbrStrands{num+1} = {self.windings[num].n_strands};\n")
                text_file.write(f"Fill{num+1} = {self.windings[num].ff};\n")
                # ---
                # Litz Approximation Coefficients were created with 4 layers
                # Thats why here a hard-coded 4 is implemented
                # text_file.write(f"NbrLayers{num+1} = {self.n_layers[num]};\n")
                text_file.write(f"NbrLayers{num+1} = 4;\n")
            text_file.write(f"AreaCell{num+1} = {self.windings[num].a_cell};\n")
            text_file.write(f"Rc{num+1} = {self.windings[num].conductor_radius};\n")

            # -- Excitation --
            # Imposed current, current density or voltage
            if self.flag_excitation_type == 'current':
                text_file.write(f"Val_EE_{num+1} = {self.current[num]};\n")
                text_file.write(f"Phase_{num+1} = Pi*{self.phase_tmp[num]};\n")
            if self.flag_excitation_type == 'current_density':
                text_file.write(f"Val_EE_{num+1} = {self.current_density[num]};\n")
            if self.flag_excitation_type == 'voltage':
                text_file.write(f"Val_EE_{num+1} = {self.voltage[num]};\n")
            print(f"Cell surface area: {self.windings[num].a_cell} \n"
                  f"Reduced frequency: {self.red_freq[num]}")
            if self.red_freq[num] > 1.25 and self.windings[num].conductor_type == "litz":
                #TODO: Allow higher reduced frequencies
                print(f"Litz Coefficients only implemented for X<=1.25")
                raise Warning
            # Reduced Frequency
            text_file.write(f"Rr{num+1} = {self.red_freq[num]};\n")

        """
        # Coordinates of the rectangular winding window
        if self.dimensionality == "2D axi":
            text_file.write("Xw1 = %s;\n" % self.p_window[4, 0])
            text_file.write("Xw2 = %s;\n" % self.p_window[5, 0])
        else:
            raise NotImplementedError("Only axi-symmetric case implemented :(")
        """
        # -- Materials --

        # Nature Constants
        text_file.write(f"mu0 = 4.e-7 * Pi;\n")
        text_file.write(f"nu0 = 1 / mu0;\n")

        # Material Properties
        # Conductor Material
        text_file.write(f"SigmaCu = {self.cond_sigma};\n")

        # Core Material
        #if self.frequency == 0:
        if self.flag_non_linear_core == 1:
            text_file.write(f"Flag_NL = 1;\n")
            text_file.write(f"Core_Material = {self.core_material};\n")
        else:
            text_file.write(f"Flag_NL = 0;\n")
            text_file.write(f"mur = {self.core_re_mu_rel};\n")
        #if self.frequency != 0:
        #    text_file.write(f"Flag_NL = 0;\n")
        #    text_file.write(f"mur = {self.core_re_mu_rel};\n")

        text_file.close()

    def simulate(self):
        """
        Initializes a onelab client. Provides the GetDP based solver with the created mesh file.
        :return:
        """
        print(f"\n---\n"
              f"Inititalize ONELAB API\n"
              f"Run Simulation\n")

        self.onelab_setup()
        # -- Simulation --
        # create a new onelab client
        c = onelab.client(__file__)

        # get model file names with correct path
        msh_file = c.getPath(self.path_mesh + 'geometry.msh')
        solver = c.getPath('ind_axi_python_controlled' + '.pro')

        # Run simulations as sub clients (non blocking??)
        mygetdp = self.onelab + 'getdp'
        c.runSubClient('myGetDP', mygetdp + ' ' + solver + ' -msh ' + msh_file + ' -solve Analysis -v2')

    # === Post-Processing ===
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
        # Mesh
        gmsh.option.setNumber("Mesh.SurfaceEdges", 0)
        view = 0

        #if self.conductor_type[0] != 'litz':
        if any(self.windings[i].conductor_type != 'litz' for i in range(0, self.n_windings)):
            # Ohmic losses (weightend effective value of current density)
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
            # Ohmic losses (weightend effective value of current density)
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

    def data_logging(self, sim_choice):
        """
        !!! not finally implemented !!!

        This method shall do the saving and loading of results! with date and time
        :return:
        """
        frequencies = None
        # Data Logging with date and time
        datetime = time.strftime("%Y%m%d_%H%M%S")

        target_femm = 'Pv_FEMM_' + datetime + '.json'
        target_femmt = 'Pv_FEMMT_' + datetime + '.json'

        # Either read old data or create
        if sim_choice != 'show':
          # pseudo 2D dataframe: ['strand number, strand radius'], 'frequency'  | const. litz radius --> ff
          if not os.path.isfile(target_femm):
              df_pv_femm = pd.DataFrame([], index=frequencies, columns=[])
              df_pv_femmt = pd.DataFrame([], index=frequencies, columns=[])
          else:
              # Read Loss Data
              df_pv_femm = pd.read_json(target_femm)
              df_pv_femmt = pd.read_json(target_femmt)
          # print(df_pv_femm, df_pv_femmt)

    def get_loss_data(self, last_n_values, loss_type='litz_loss'):
        """
        Returns the last n values from the chosen loss type logged in the result folder.
        :param last_n_values:
        :param loss_type:
        :return:
        """
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

    # === Alternative FEMM Solver ===
    def femm_reference(self, freq, current, sigma, sign=[1], non_visualize=0):
        """
        Allows reference simulations with the 2D open source electromagnetic FEM tool FEMM.
        Helpful to validate changes (especially in the Prolog Code).
        :param sign:
        :param non_visualize:
        :param freq:
        :param current:
        :param sigma:
        :return:
        """
        if sign is None:
            sign = [1, 1]
        if not os.path.isdir(self.path + "/" + self.path_res_FEMM):
            os.mkdir(self.path + "/" + self.path_res_FEMM)

        # == Pre Geometry ==
        self.high_level_geo_gen()
        #self.ei_axi()

        if self.n_air_gaps != 1:
            raise NotImplementedError

        # == Init ==
        femm.openfemm(non_visualize)
        femm.newdocument(0)
        femm.mi_probdef(freq, 'meters', 'axi', 1.e-8, 0, 30)

        # == Materials ==
        femm.mi_addmaterial('Ferrite', self.core_re_mu_rel, self.core_re_mu_rel, 0, 0, 0, 0, 0, 1, 0, 0, 0)
        femm.mi_addmaterial('Air', 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0)
        if self.conductor_type[0] == "litz":
            femm.mi_addmaterial('Copper', 1, 1, 0, 0, sigma, 0, 0, 1, 5, 0, 0, self.windings[0].n_strands, 2*1000*self.windings[0].strand_radius)  # type := 5. last argument
            print(f"Strandsnumber: {self.windings[0].n_strands}")
            print(f"Strandsdiameter in mm: {2 * 1000 * self.strand_radius[0]}")
        if self.conductor_type[0] == "solid":
            femm.mi_addmaterial('Copper', 1, 1, 0, 0, sigma, 0, 0, 1, 0, 0, 0, 0, 0)

        # == Circuit ==
        # coil as seen from the terminals.
        femm.mi_addcircprop('Primary', current[0]*sign[0], 1)
        if self.component_type == 'transformer':
            femm.mi_addcircprop('Secondary', current[1]*sign[1], 1)

        # == Geometry ==
        # Add core
        femm.mi_drawline(0, self.p_air_gaps[0, 1], self.p_air_gaps[1, 0], self.p_air_gaps[1, 1])
        femm.mi_drawline(self.p_air_gaps[1, 0], self.p_air_gaps[1, 1],  self.p_window[4, 0], self.p_window[4, 1])
        femm.mi_drawline(self.p_window[4, 0], self.p_window[4, 1], self.p_window[5, 0], self.p_window[5, 1])
        femm.mi_drawline(self.p_window[5, 0], self.p_window[5, 1], self.p_window[7, 0], self.p_window[7, 1])
        femm.mi_drawline(self.p_window[7, 0], self.p_window[7, 1], self.p_window[6, 0], self.p_window[6, 1])
        femm.mi_drawline(self.p_window[6, 0], self.p_window[6, 1], self.p_air_gaps[3, 0], self.p_air_gaps[3, 1])
        femm.mi_drawline(self.p_air_gaps[3, 0], self.p_air_gaps[3, 1], 0, self.p_air_gaps[2, 1])
        femm.mi_drawline(0, self.p_air_gaps[2, 1], 0, self.p_outer[2, 1])
        femm.mi_drawline(0, self.p_outer[2, 1], self.p_outer[3, 0], self.p_outer[3, 1])
        femm.mi_drawline(self.p_outer[3, 0], self.p_outer[3, 1], self.p_outer[1, 0], self.p_outer[1, 1])
        femm.mi_drawline(self.p_outer[1, 0], self.p_outer[1, 1], 0, self.p_outer[0, 1])
        femm.mi_drawline(0, self.p_outer[0, 1], 0, self.p_air_gaps[0, 1])
        # Add Coil
        """
        femm.mi_drawrectangle(self.p_window[4, 0]+self.core_cond_isolation, self.p_window[4, 1]+self.core_cond_isolation, self.p_window[7, 0]-self.core_cond_isolation, self.p_window[7, 1]-self.core_cond_isolation)
        femm.mi_addblocklabel(self.p_window[7, 0]-2*self.core_cond_isolation, self.p_window[7, 1]-2*self.core_cond_isolation)
        femm.mi_selectlabel(self.p_window[7, 0]-2*self.core_cond_isolation, self.p_window[7, 1]-2*self.core_cond_isolation)
        femm.mi_setblockprop('Copper', 0, 1, 'icoil', 0, 0, 1)
        femm.mi_clearselected()
        """
        for num in range(0, self.n_windings):
            if self.conductor_type[0] == "litz" or self.conductor_type[0] == "solid":
                for i in range(0, int(self.p_conductor[num].shape[0] / 5)):
                    # 0: center | 1: left | 2: top | 3: right | 4.bottom
                    femm.mi_drawarc(self.p_conductor[num][5*i+1][0], self.p_conductor[num][5*i+1][1], self.p_conductor[num][5*i+3][0], self.p_conductor[num][5*i+3][1], 180, 2.5)
                    femm.mi_addarc(self.p_conductor[num][5*i+3][0], self.p_conductor[num][5*i+3][1], self.p_conductor[num][5*i+1][0], self.p_conductor[num][5*i+1][1],  180, 2.5)
                    femm.mi_addblocklabel(self.p_conductor[num][5*i][0], self.p_conductor[num][5*i][1])
                    femm.mi_selectlabel(self.p_conductor[num][5*i][0], self.p_conductor[num][5*i][1])
                    if num == 0:
                        femm.mi_setblockprop('Copper', 1, 0, 'Primary', 0, 0, 1)
                    if num == 1:
                        #femm.mi_setblockprop('Copper', 0, 1e-4, 'Secondary', 0, 0, 1)
                        femm.mi_setblockprop('Copper', 1, 0, 'Secondary', 0, 0, 1)
                    femm.mi_clearselected

        # Define an "open" boundary condition using the built-in function:
        femm.mi_makeABC()
        """
        # Alternative BC
        region_add = 1.1

        femm.mi_drawrectangle(0, region_add*self.p_outer[0][1], region_add*self.p_outer[3][0], region_add*self.p_outer[3][1])
        # mi_addboundprop('Asymptotic', 0, 0, 0, 0, 0, 0, 1 / (para.mu_0 * bound.width), 0, 2); % Mixed
        femm.mi_addboundprop('Asymptotic', 0, 0, 0, 0, 1, 50, 0, 0, 1)
        femm.mi_selectsegment(region_add*self.p_outer[3][0], region_add*self.p_outer[3][1])
        femm.mi_setsegmentprop('Asymptotic', 1, 1, 0, 0)
        """

        # == Labels/Designations ==

        # Label for core
        femm.mi_addblocklabel(self.p_outer[3, 0]-0.001, self.p_outer[3, 1]-0.001)
        femm.mi_selectlabel(self.p_outer[3, 0]-0.001, self.p_outer[3, 1]-0.001)
        femm.mi_setblockprop('Ferrite', 1, 0, '<None>', 0, 0, 0)
        femm.mi_clearselected()

        # Labels for air
        femm.mi_addblocklabel(0.001, 0)
        femm.mi_selectlabel(0.001, 0)
        femm.mi_setblockprop('Air', 1, 0, '<None>', 0, 0, 0)
        femm.mi_clearselected()
        femm.mi_addblocklabel(self.p_outer[3, 0]+0.001, self.p_outer[3, 1]+0.001)
        femm.mi_selectlabel(self.p_outer[3, 0]+0.001, self.p_outer[3, 1]+0.001)
        femm.mi_setblockprop('Air', 1, 0, '<None>', 0, 0, 0)
        femm.mi_clearselected()

        # Now, the finished input geometry can be displayed.
        femm.mi_zoomnatural()
        femm.mi_saveas(self.path_res_FEMM + 'coil.fem')
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

    # === Litz Approximation ===

    def pre_simulate(self):
        """
        Used to determine the litz-approximation coefficients.
        :return:
        """
        for num in range(0, self.n_windings):
            if self.windings[num].conductor_type == 'litz':
                # ---
                # Litz Approximation Coefficients were created with 4 layers
                # Thats why here a hard-coded 4 is implemented
                # if os.path.isfile(self.path + f"/pre/coeff/pB_RS_la{self.windings[num].ff}_{self.n_layers[num]}layer.dat"):
                if os.path.isfile(self.path + f"/pre/coeff/pB_RS_la{self.windings[num].ff}_4layer.dat"):
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
        text_file = open(self.path + "/pre/PreParameter.pro", "w")
        # ---
        # Litz Approximation Coefficients are created with 4 layers
        # Thats why here a hard-coded 4 is implemented
        #text_file.write(f"NbrLayers = {self.n_layers[num]};\n")
        text_file.write(f"NbrLayers = 4;\n")
        text_file.write(f"Fill = {self.windings[num].ff};\n")
        print("Here")
        text_file.write(f"Rc = {self.strand_radius[num]};\n")  # double named!!! must be changed
        text_file.close()
        self.onelab_setup()
        c = onelab.client(__file__)
        cell_geo = c.getPath('pre/cell.geo')

        # Run gmsh as a sub client
        mygmsh = self.onelab + 'gmsh'
        c.runSubClient('myGmsh', mygmsh + ' ' + cell_geo + ' -2 -v 2')

        modes = [1, 2]  # 1 = "skin", 2 = "proximity"
        reduced_frequencies = np.linspace(0, 1.25, 6)  # must be even
        for mode in modes:
            for rf in reduced_frequencies:
                # -- Pre-Simulation Settings --
                text_file = open(self.path + "/pre/PreParameter.pro", "w")
                text_file.write(f"Rr_cell = {rf};\n")
                text_file.write(f"Mode = {mode};\n")
                # Litz Approximation Coefficients are created with 4 layers
                # Thats why here a hard-coded 4 is implemented
                # text_file.write(f"NbrLayers = {self.n_layers[num]};\n")
                text_file.write(f"NbrLayers = 4;\n")

                text_file.write(f"Fill = {self.windings[num].ff};\n")
                text_file.write(f"Rc = {self.strand_radius[num]};\n")  # double named!!! must be changed
                text_file.close()

                # get model file names with correct path
                input_file = c.getPath('pre/cell_dat.pro')
                cell = c.getPath('pre/cell.pro')

                # Run simulations as sub clients
                mygetdp = self.onelab + 'getdp'
                c.runSubClient('myGetDP', mygetdp + ' ' + cell + ' -input ' + input_file + ' -solve MagDyn_a -v2')

        # Formatting stuff
        # Litz Approximation Coefficients are created with 4 layers
        # Thats why here a hard-coded 4 is implemented
        #files = [self.path + f"/pre/coeff/pB_RS_la{self.windings[num].ff}_{self.n_layers[num]}layer.dat",
        #         self.path + f"/pre/coeff/pI_RS_la{self.windings[num].ff}_{self.n_layers[num]}layer.dat",
        #         self.path + f"/pre/coeff/qB_RS_la{self.windings[num].ff}_{self.n_layers[num]}layer.dat",
        #         self.path + f"/pre/coeff/qI_RS_la{self.windings[num].ff}_{self.n_layers[num]}layer.dat"]
        files = [self.path + f"/pre/coeff/pB_RS_la{self.windings[num].ff}_4layer.dat",
                 self.path + f"/pre/coeff/pI_RS_la{self.windings[num].ff}_4layer.dat",
                 self.path + f"/pre/coeff/qB_RS_la{self.windings[num].ff}_4layer.dat",
                 self.path + f"/pre/coeff/qI_RS_la{self.windings[num].ff}_4layer.dat"]
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

    # ==== Front-End Methods =====
    # === Geometry Definitions ===
    def update_core(self, core_type, core_re_mu_rel = 3000,
                                     core_im_mu_rel = 1750 * np.sin(10 *np.pi/180),
                                     core_im_epsilon_rel = 6e+4 * np.sin(20 *np.pi/180),
                                     core_material = 95,
                                     non_linear = 0,
                            **kwargs):
        """
        - One positional parameter core_type
        - All core parameters can be passed or adjusted by keyword calling
            - Allows single parameter changing
            - Depending on core_type
            - Strict keyword usage!
        :param core_type:
        :param kwargs:
            - Case "2D, axisym., EI": 'core_w', 'window_w', 'window_h'
            - Case "3D, EI": ...tba...
        :return:
        """

        if self.core_update_count == 0:
            print(f"Initialize the magnetic Core as {self.core_type}-type with some standard values.\n"
                  f"---")
        else:
            print(f"Update the magnetic Core to {self.core_type}-type with following parameters: {kwargs}\n"
                  f"---")
        self.core_update_count+=1

        # Material Properties
        self.flag_non_linear_core   = non_linear
        self.core_re_mu_rel         = core_re_mu_rel
        self.core_im_mu_rel         = core_im_mu_rel
        self.core_im_epsilon_rel    = core_im_epsilon_rel
        self.core_material          = core_material


        self.core_type = core_type
        if self.core_type == "EI":
            if self.dimensionality == "2D axi":
                for key, value in kwargs.items():
                    if key == 'core_w':
                        self.core_w = value
                    if key == 'window_w':
                        self.window_w = value
                    if key == 'window_h':
                        self.window_h = value
        if self.core_type == "EI":
            if self.dimensionality == "3D":
                # tba 3D Group
                None

    def update_air_gaps(self, method="center", n_air_gaps=[], position_tag=[0], air_gap_position=[0], air_gap_h=[0.001],
                        **kwargs):
        """
        - "self.air_gaps" is a list with [position_tag, air_gap_position, air_gap_h, c_air_gap]
           - position_tag: specifies the gapped "leg"
           - air_gap_position: specifies the coordinate of the air gap's center point along the specified leg
           - air_gap_h: height/length of the air gap
           - c_air_gap: mesh accuracy factor
        - "EI 2D axi": position_tag = 0  # '-1': left leg | '0': center leg | '1': right leg
        :param n_air_gaps:
        :param position_tag:
        :param air_gap_h:
        :param air_gap_position:
        :param method: "random", "center", "percent", "manually"
        :return:
        """
        print(f"Update the air gaps.\n"
              f"---")

        # Update the mesh accuracy of the window
        self.c_window = self.window_w/20 * self.s

        # Rewrite variables
        self.n_air_gaps = n_air_gaps
        self.air_gaps = np.empty((self.n_air_gaps, 4))
        self.c_air_gap = [None] * self.n_air_gaps


        # Update air gaps with chosen method

        # Center
        if method == "center" and self.dimensionality == "2D axi":
            if self.n_air_gaps > 1:
                print(f"{self.n_air_gaps} are too many air gaps for the 'center' option!")
                raise Warning
            else:
                self.c_air_gap[0] = air_gap_h[0] * self.s
                self.air_gaps[0, :] = np.array([0, 0, air_gap_h[0], self.c_air_gap[0]])

        # Deterministic
        if (method == "manually" or method == "percent") and self.dimensionality == "2D axi":
            for i in range(0, self.n_air_gaps):
                if method == "percent":
                    air_gap_position[i] = air_gap_position[i] / 100 * (self.window_h - air_gap_h[i]) - (
                                self.window_h / 2 - air_gap_h[i] / 2)
                # Overlapping Control
                for j in range(0, self.air_gaps.shape[0]):
                    if self.air_gaps[j, 1]+self.air_gaps[j, 2]/2 > air_gap_position[i] > self.air_gaps[j, 1]-self.air_gaps[j, 2]/2:
                        if position_tag[i] == self.air_gaps[j, 0]:
                            print(f"Overlapping Air Gap")
                            #raise Warning
                    else:
                        # self.c_air_gap[i] = air_gap_h[i] * self.s
                        self.c_air_gap[i] = self.c_window
                        #print(f"c_window: {self.c_window}")
                        self.air_gaps[i, :] = np.array([position_tag[i], air_gap_position[i], air_gap_h[i], self.c_air_gap[i]])

        # Random
        if method == "random" and self.dimensionality == "2D axi":
            position_tag = [0] * self.n_air_gaps

            i = 0
            while i in range(0, self.n_air_gaps):
                height = np.random.rand(1) * 0.001 + 0.001
                position = np.random.rand(1) * (self.window_h - height) - (self.window_h / 2 - height / 2)
                self.c_air_gap[i] = height * self.s
                # Overlapping Control
                for j in range(0, self.air_gaps.shape[0]):
                    if self.air_gaps[j, 1]+self.air_gaps[j, 2]/2 > position > self.air_gaps[j, 1]-self.air_gaps[j, 2]/2:
                        if position_tag[i] == self.air_gaps[j, 0]:
                            print(f"Overlapping air Gaps have been corrected")
                else:
                    self.air_gaps[i, :] = np.array([position_tag[i], position, height, self.c_air_gap[i]])
                    i += 1

        # Optional modelling of a stray path
        for key, value in kwargs.items():
            if key == 'stray_path':
                self.dedicated_stray_path = True
                self.start_i, self.end_i = value[0]
                self.added_bar = value[1]
                print(f"Stray Path arguments{self.start_i, self.end_i, self.added_bar}")

    def update_conductors(self, n_turns=[], conductor_type=[], winding=[], scheme=[], conductor_radii=[],
                          litz_para_type=[], ff=[], strands_numbers=[], strand_radii=[], thickness=[],
                          wrap_para=[], cond_cond_isolation=[], core_cond_isolation=[]):
        """
        This Method allows the user to easily update/initialize the Windings in terms of their conductors and
        arrangement.
        :param cond_cond_isolation:
        :param core_cond_isolation:
        :param n_turns:
        :param strand_radii:
        :param strands_numbers:
        :param conductor_radii:
        :param conductor_radius:
        :param n_conductors:
        :param conductor_type:  - "stacked"  # Vertical packing of conductors
                                - "full"     # One massive Conductor in each window
                                - "foil"     # Horizontal packing of conductors
                                - "solid"    # Massive wires
                                - "litz"     # Litz wires
        :return:
        """
        print(f"Update the conductors...\n"
                  f"---")

        # - - - - - - - - - - - - - - - - - - - - - CHECK Input Parameters - - - - - - - - - - - - - - - - - - - - - - -
        if self.component_type == "inductor":
            if len(n_turns) != 1 or len(conductor_type) != 1:
                print(f"Wrong number of conductor parameters passed to inductor model!")
                raise Warning

        if self.component_type == "transformer":
            if len(n_turns) != 2 or len(conductor_type) != 2:
                print(f"Wrong number of conductor parameters passed to transformer model!")
                raise Warning

        # - - - - - - - - - - - - - - - - - Definition of Virtual Winding Windows - - - - - - - - - - - - - - - - - - -
        # Inductor: One virtual winding window
        if self.component_type == "inductor":
            self.vw_type = "full_window"  # One Virtual Winding Window

        # Two winding transformer
        # One virtual winding window
        if not self.dedicated_stray_path:
            if len(winding) == 1 and winding[0] == "interleaved":
                self.vw_type = "full_window"  # One Virtual Winding Window
            if len(winding) == 2:
                self.vw_type = "center"  # Two Virtual Winding Windows #TODO: Center anpassen (nicht mehr nur Nulllinie)
        else:
            None
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
        self.core_cond_isolation = core_cond_isolation
        self.cond_cond_isolation = cond_cond_isolation

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
                self.windings[i].a_cell = np.pi * self.windings[i].conductor_radius ** 2  # Cross section of the solid conductor

            if self.windings[i].conductor_type == 'litz':
                if litz_para_type[i] == 'implicite_ff':
                    self.update_litz_configuration(num=i,
                                                   litz_parametrization_type=litz_para_type[i],
                                                   conductor_radius=conductor_radii[i],
                                                   n_strands=strands_numbers[i],
                                                   strand_radius=strand_radii[i])
                if litz_para_type[i] == 'implicite_litz_radius':
                    self.update_litz_configuration(num=i,
                                                   litz_parametrization_type=litz_para_type[i],
                                                   ff=ff[i],
                                                   n_strands=strands_numbers[i],
                                                   strand_radius=strand_radii[i])
                if litz_para_type[i] == 'implicite_strands_number':
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

    def update_litz_configuration(self, num=0, litz_parametrization_type='implicite_ff', strand_radius=None, ff=None,
                                  conductor_radius=None, n_strands=None):
        """
        - updates the conductor #num (simple transformer: num=0 -> primary winding,
                                                                     num=1 -> secondary winding)
        - used to change litz configurations
        - also used at first initialisation of the geometry
        - needed to always make sure that the relation between litz parameters (strand radius, fill factor, number of
          layers/strands and conductor/litz radius) is valid and consistent
        - 4 parameters, 1 degree of freedom (dof)
        :param num:
        :param litz_parametrization_type:
        :param strand_radius:
        :param ff:
        :param conductor_radius:
        :param n_strands:
        :return:
        """
        # Choose one of three different parametrization types
        if litz_parametrization_type == 'implicite_ff':
            ff_exact = n_strands * strand_radius ** 2 / conductor_radius ** 2
            print(f"Exact fill factor: {ff_exact}")
            ff = np.around(ff_exact, decimals=2)  # TODO: Interpolation instead of rounding

        if litz_parametrization_type == 'implicite_litz_radius':
            conductor_radius = np.sqrt(n_strands*strand_radius ** 2 / ff)

        if litz_parametrization_type == 'implicite_strands_number':
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

    # === Standard Simulations ===
    def single_simulation(self, freq, current, phi=[], skin_mesh_factor=1, NL_core=0):
        """
        - can be used for a single simulation
        - no sweeping at all
        :param freq:
        :param current:
        :return:
        """
        self.high_level_geo_gen(frequency=freq, skin_mesh_factor=skin_mesh_factor)
        self.generate_mesh()
        self.excitation(f=freq, i=current, phases=phi)  # frequency and current
        self.file_communication()
        self.pre_simulate()
        self.simulate()
        self.visualize()

    def get_inductances(self, I0, op_frequency=0, skin_mesh_factor=1, visualize=0):
        """

        :param mesh_accuracy:
        :param I0:
        :param op_frequency:
        :return:
        """


        # Remove "old" Inductance Logs
        try:
            os.remove(self.path + "/" + self.path_res_values + "L_11.dat")
            os.remove(self.path + "/" + self.path_res_values + "L_22.dat")
        except:
            print("Could not find Inductance logs")

        # -- Inductance Estimation --
        self.mesh(frequency=op_frequency, skin_mesh_factor=skin_mesh_factor)

        if self.valid:
            frequencies = [op_frequency] * 2
            currents = [[I0, 0], [0, I0]]
            phases = [[0, 180], [0, 180]]

            self.excitation_sweep(frequencies=frequencies, currents=currents, phi=phases[0], show_last=visualize)

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
            K_21 = Phi_21 / Phi_22
            K_12 = Phi_12 / Phi_11
            k = n/np.abs(n) * (K_21*K_12)**0.5
            print(f"Coupling Factors:\n"
                  f"K_12 = Phi_21 / Phi_22 = {K_12}\n"
                  f"K_21 = Phi_12 / Phi_11 = {K_21}\n"
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
            self.M = k * (self.L_11*self.L_22)**0.5
            M_ = self.L_11 * K_12  # Only to proof correctness - ideally: M = M_ = M__
            M__ = self.L_22 * K_21  # Only to proof correctness - ideally: M = M_ = M__
            print(f"\n"
                  f"Main/Counter Inductance:\n"  
                  f"M = k * Sqrt(L_11 * L_22) = {self.M}\n"
                  f"M_ = L_11 * K_12 = {M_}\n"
                  f"M__ = L_22 * K_21 = {M__}\n"
                  )

            # Stray Inductance with 'Turns Ratio' n as 'Transformation Ratio' ü
            L_s1 = self.L_11 - self.M * n
            L_s2 = self.L_22 - self.M / n
            L_h = self.M * n
            print(f"\n"
                  f"T-ECD (primary side transformed):\n"
                  f"[Underdetermined System: 'Transformation Ratio' := 'Turns Ratio']\n"
                  f"    - Transformation Ratio: ü\n"
                  f"    - Primary Side Stray Inductance: L_s1\n"
                  f"    - Secondary Side Stray Inductance: L_s2\n"
                  f"    - Primary Side Main Inductance: L_h\n"
                  f"ü := n = {n}\n"
                  f"L_s1 = L_11 - M * n = {L_s1}\n"
                  f"L_s2 = L_22 - M / n = {L_s2}\n"
                  f"L_h = M * n = {L_h}\n"
                  )

            # Stray Inductance concentrated on Primary Side
            self.ü_conc = self.M / self.L_22
            self.L_s_conc = (1 - k**2) * self.L_11
            self.L_h_conc = self.M**2 / self.L_22
            print(f"\n"
                  f"T-ECD (primary side concentrated):\n"
                  f"[Underdetermined System: ü := M / L_22  -->  L_s2 = L_22 - M / n = 0]\n"
                  f"    - Transformation Ratio: ü\n"
                  f"    - (Primary) Stray Inductance: L_s1\n"
                  f"    - Primary Side Main Inductance: L_h\n"
                  f"ü := M / L_22 = k * Sqrt(L_11 / L_22) = {self.ü_conc}\n"
                  f"L_s1 = (1 - k^2) * L_11 = {self.L_s_conc}\n"
                  f"L_h = M^2 / L_22 = k^2 * L_11 = {self.L_h_conc}\n"
                  )
            self.visualize()

        else:
            print(f"Invalid Geommetry Data!")

    def excitation_sweep(self, frequencies=[], currents=[], phi=[0, 180], show_last=False):
        """
        Performs a sweep simulation for frequency-current pairs. Both values can
        be passed in lists of the same length. The mesh is only created ones (fast sweep)!

        Example Code:
            1 geo = MagneticComponent()
            2 geo.mesh()
            3 fs = np.linspace(0, 250000, 6)
            4 cs = [10, 2, 1, 0.5, 0.2, 0.1]
            5 geo.excitation_sweep(frequencies=fs, currents=cs)
        :param currents:
        :param frequencies:
        :param show_last:
        :return:
        """
        #if len(frequencies) != len(currents):
        #    print('len(frequencies) != len(currents)')
        #    raise Exception
        for i in range(0, len(frequencies)):
            self.excitation(f=frequencies[i], i=currents[i], phases=phi)  # frequency and current
            self.file_communication()
            self.pre_simulate()
            self.simulate()
        if show_last:
            self.visualize()

    def get_Core_Loss(self, Ipeak=[10, 10], ki=1, alpha=1.2, beta=2.2, t_rise=3e-6, t_fall=3e-6, f_switch=100000,
                      skin_mesh_factor=0.5):
        """

        :param ki:
        :param alpha:
        :param beta:
        :param t_rise:
        :param t_fall:
        :param f_switch:
        :return:
        """
        self.core_loss_simlation = 1
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
        #self.high_level_geo_gen(frequency=0, skin_mesh_factor=skin_mesh_factor)
        #self.generate_mesh()
        self.excitation(f=f_switch, i=Ipeak, phases=[0, 180])  # frequency and current
        self.file_communication()
        self.pre_simulate()
        self.simulate()
        self.visualize()

    # === Special/Experimental Simulations ===
    def adaptive_single_simulation(self, freq, current, phi=[], max_iter=1, local=0):
        """
        - can be used for a single simulation
        - no sweeping at all
        :param freq:
        :param current:
        :return:
        """
        self.high_level_geo_gen(frequency=freq)
        self.generate_mesh()
        self.excitation(f=freq, i=current, phases=phi)  # frequency and current
        self.file_communication()
        self.pre_simulate()
        self.simulate()
        self.visualize()

        for i in range(0, max_iter):
            self.refine_mesh(local=local)
            self.excitation(f=freq, i=current, phases=phi)  # frequency and current
            self.file_communication()
            self.pre_simulate()
            self.simulate()
            self.visualize()

    def litz_loss_comparison(self, ff, n_strands, strand_radius, sim_choice, sweep_parameter='fill_factor', nametag=''):
        """

        :param nametag:
        :param sweep_parameter:
        :param ff:
        :param n_strands:
        :param strand_radius:
        :param sim_choice:
        :return:
        """
        # TODO: Bring it to newest Update (for example: strands_numbers instead of layer_number)

        # Excitation: fixed sweep parameters
        frequencies = np.linspace(0, 250000, 6)
        currents = [2, 2, 2, 2, 2, 2]

        # Find Json logfile
        target_femm = self.path + nametag + 'Pv_FEMM_' + sweep_parameter + '.json'
        target_femmt = self.path + nametag + 'Pv_FEMMT_' + sweep_parameter + '.json'

        # Update model to chosen litz parameters
        self.update_litz_configuration(litz_parametrization_type='implicite_litz_radius', n_strands=n_strands,
                                       strand_radius=strand_radius, ff=ff)
        for fill_factors in self.ff:
            if fill_factors < 0.4 or fill_factors > 0.9:
                print(f"Skip simulation with non realistic fill factor {fill_factors}")

            else:
                # Either read old data or create new dataframe
                if sim_choice != 'show':
                    if not os.path.isfile(target_femm):
                        df_pv_femm = pd.DataFrame([], index=None, columns=[])
                        df_pv_femm.insert(loc=0, column=f"Frequencies", value=frequencies)
                    else:
                        df_pv_femm = pd.read_json(target_femm)
                    if not os.path.isfile(target_femmt):
                        df_pv_femmt = pd.DataFrame([], index=None, columns=[])
                        df_pv_femmt.insert(loc=0, column=f"Frequencies", value=frequencies)
                    else:
                        df_pv_femmt = pd.read_json(target_femmt)

                # Column tags
                # TODO: Which winding is meant here? Only Inductor case?
                femmt_col_tag = f"onelab, {self.n_strands}, {self.strand_radius}, {self.ff}"
                femm_col_tag = f"femm, {self.n_strands}, {self.strand_radius}, {self.ff}"

                # Prevent from rewriting already used parametrization
                # rename duplicated tag by adding randomness | string of length 6
                femmt_col_tags = df_pv_femmt.columns.values.tolist()
                femm_col_tags = df_pv_femm.columns.values.tolist()
                random_id = id_generator()
                if any(femmt_col_tag in s for s in femmt_col_tags):
                    femmt_col_tag = femmt_col_tag + random_id
                if any(femm_col_tag in s for s in femm_col_tags):
                    femm_col_tag = femm_col_tag + random_id

                # -- Reference simulation with FEMM --
                if sim_choice == 'both' or sim_choice == 'femm':
                    pv_femm = []
                    # Iterate on frequency
                    for f in frequencies:
                        self.femm_reference(freq=f, current=2, sigma=58, non_visualize=1)
                        pv_femm.append(self.tot_loss_femm.real)
                    # Add new or rewrite old data // Maybe ask for replacement and wait for 30 s then go on...
                    df_pv_femm.insert(loc=1, column=femm_col_tag, value=np.transpose(pv_femm), allow_duplicates=True)
                    # Save dataframes in .json files
                    df_pv_femm.to_json(target_femm, date_format=None)

                # -- Onelab simulation with FEMMT --
                if sim_choice == 'both' or sim_choice == 'femmt':
                    # High level geometry generation + Generate Mesh
                    self.mesh()
                    # Iterate on frequency
                    self.excitation_sweep(currents=currents, frequencies=frequencies)

                    # Get losses from Onelab result file
                    data = self.get_loss_data(last_n_values=len(frequencies), loss_type='litz_loss')

                    pv_femmt = []
                    for lines in data:
                        print(re.findall(r"[-+]?\d*\.\d+|\d+", lines[0]))
                        fls = re.findall(r"[-+]?\d*\.\d+|\d+", lines[0])
                        if len(fls) == 3:
                            pv_femmt.append(float(fls[1]))
                        else:
                            pv_femmt.append(0.0)
                            warnings.warn("There is something wrong with the Loss data!")

                    # Add new or rewrite old data
                    df_pv_femmt.insert(loc=1, column=femmt_col_tag, value=np.transpose(pv_femmt))

                    # Save dataframes in .json files
                    df_pv_femmt.to_json(target_femmt, date_format=None)

    def load_litz_loss_logs(self, tag=''):
        """
        Used with litz_loss_comparison()
        :param tag:
        :return:
        """
        # TODO: Bring it to newest Update (for example: strands_numbers instead of layer_number)

        # Read Json to pandas dataframe
        df_femm = pd.read_json(self.path + tag + 'Pv_FEMM_fill_factor.json', convert_axes=False)
        df_femmt = pd.read_json(self.path + tag + 'Pv_FEMMT_fill_factor.json', convert_axes=False)

        # Make frequencies the index
        df_Pv_femm = df_femm.set_index('Frequencies')
        df_Pv_femmt = df_femmt.set_index('Frequencies')

        # Correct 0Hz error in femm
        print(df_Pv_femm)
        print(df_Pv_femm.iloc[:, 0])
        df_Pv_femm.iloc[0] = df_Pv_femm.iloc[0] * 2

        """
        # Error plotting
        error = pd.DataFrame([], columns=[])
        for col in range(0, len(df_Pv_femmt.columns)):
            error.insert(loc=0, column=f"ff{fillfactors[-col - 1]}",
                         value=(df_Pv_femmt.iloc[:, col] - df_Pv_femm.iloc[:, col]).div(df_Pv_femmt.iloc[:, col]))
        """

        # Error plotting
        error = pd.DataFrame([], columns=[])
        for col in range(0, len(df_Pv_femmt.columns)):
            error.insert(loc=0, column=df_Pv_femmt.columns[col].replace("onelab, ", ""),
                         value=(df_Pv_femmt.iloc[:, col] - df_Pv_femm.iloc[:, col]).div(df_Pv_femmt.iloc[:, col]))

        # print("Error: \n", error)
        # print("FEMM: \n", df_Pv_femm)
        # print("FEMMT: \n", df_Pv_femmt)

        # Single Plots
        # df_Pv_femm.plot(title="FEMM")
        # df_Pv_femmt.plot(title="FEMMT")
        error.plot(title=r"$\frac{P_{v,onelab}-P_{v,femm}}{P_{v,onelab}}$")

        # Two Plots in one
        ax = df_Pv_femm.plot(marker="p")
        df_Pv_femmt.plot(ax=ax, marker="o", title=r"Losses in stranded conductors", xlabel="Frequency in Hz", ylabel="Losses in W")

        # ax.text(25000, 1, 'FEMM', color='r', ha='right', rotation=0, wrap=True)
        # ax.text(25000, 0.5, 'FEMMT', color='g', ha='right', rotation=0, wrap=True)
        plt.show()

# ----------------------------------------------------------------------------------------------------------------------
#  ===== Static Functions  =====
#  Used somewhere in the Code


def inner_points(a, b, input_points):
    """
    Returns the input points that have a common coordinate as the two
    interval borders a and b
    :param a:
    :param b:
    :param input_points:
    :return:
    """
    [min, max] = [None, None]
    output = input_points
    dim = None
    # Find equal dimension
    for i in range(0, 3):
        if a[i] == b[i] and a[i] != 0:
            dim = i
    if dim == None:
        raise Exception("Given points do not have a common dimension")
    n = 0
    while n < output.shape[0]:
        if a[dim] != output[n, dim]:
            output = np.delete(output, n, 0)
        else:
            n += 1
    if output.shape[0] == 0:
        raise Exception("Not implemented Error: No air gaps between interval borders")
    if output.shape[0] % 2 == 1:
        raise Exception("Odd number of input points")
    if dim == 2:
        raise Exception("Not implemented Error: Only 2D is implemeted")
    dim2 = (dim+1) % 2
    if output.shape[0] >= 2:
        argmax = np.argmax(output[:, dim2])
        output = np.delete(output, argmax, 0)
        argmin = np.argmin(output[:, dim2])
        output = np.delete(output, argmin, 0)
        #if output.shape[0] == 0:
            #print("Only one air gap in this leg. No island needed.")
    return output


def min_max_inner_points(a, b, input_points):
    """
    Returns the input points that have a common coordinate and
    the minimum distance from the interval borders.
    :param a:
    :param b:
    :param input_points:
    :return:
    """

    [min, max] = [None, None]
    buffer = input_points
    dim = None
    # Find equal dimension
    for i in range(0, 3):
        if a[i] == b[i] and a[i] != 0:
            dim = i
    if dim == None:
        raise Exception("Given points do not have a common dimension")
    n = 0
    while n < buffer.shape[0]:
        if a[dim] != buffer[n, dim]:
            buffer = np.delete(buffer, n, 0)
        else:
            n += 1
    if buffer.shape[0] == 0:
        print("No air gaps between interval borders")
    if buffer.shape[0] % 2 == 1:
        raise Exception("Odd number of input points")
    if dim == 2:
        raise Exception("Not implemented Error: Only 2D is implemeted")
    dim2 = (dim+1) % 2
    if buffer.shape[0] >= 2:
        argmax = np.argmax(buffer[:, 1])
        max = buffer[argmax]
        argmin = np.argmin(buffer[:, 1])
        min = buffer[argmin]
    return [min, max]


def call_for_path(destination, config_file="config.json"):
    """
    asks the user to enter the filepath of a destinated file WITHOUT the suffix
    stores a the filepath as a python string declaration in a config file
    returns the filepath
    :param destination:
    :param config_file:
    :return:
    """
    # pickle_file = open(config_file, "w")
    # path = input(f"Please enter the parent folder path of {destination} in ways of 'C:.../onelab-Windows64/': ")
    # pickle.dumps(path, pickle_file) # f"{destination} = '{path}'\n")
    # pickle_file.close()

    # Find out the path of installed module, or in case of running directly from git, find the path of git repository
    module_file_path = pathlib.Path(__file__).parent.absolute()

    path = input(f"Please enter the parent folder path of {destination} in ways of 'C:.../onelab-Windows64/': ")
    dict = {"onelab": path}
    file = open(module_file_path / config_file, 'w', encoding='utf-8')
    json.dump(dict, file, ensure_ascii=False)
    file.close()

    return path


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


def NbrStrands(NbrLayers):
    """
    Returns the number of strands in a hexagonal litz winding with a
    specified number of layers (NbrLayers). CAUTION: Zero number of
    layers corresponds to a single strand.
    :param NbrLayers:
    :return:
    """
    return 3*(NbrLayers+1)**2 - 3*(NbrLayers+1) + 1


def NbrLayers(NbrStrands):
    """
    Returns the number of layers in a hexagonal litz winding with a
    specified number of strands (NbrStrands). CAUTION: Zero number of
    layers corresponds to a single strand.
    :param NbrStrands:
    :return:
    """
    return np.sqrt(0.25+(NbrStrands-1)/3)-0.5


def fft(period_vector_t_i: npt.ArrayLike, sample_factor: float = 1000, plot: str = 'no', rad: str ='no',
        f0: Union[float, None]=None, title: str='ffT') -> npt.NDArray[list]:
    """
    A fft for a input signal. Input signal is in vector format and should include one period.

    Output vector includes only frequencies with amplitudes > 1% of input signal

    Minimal example:
    example_waveform = np.array([[0, 1.34, 3.14, 4.48, 6.28],[-175.69, 103.47, 175.69, -103.47,-175.69]])
    out = fft(example_waveform, plot=True, rad='yes', f0=25000, title='ffT input current')

    :param period_vector_t_i: numpy-array [[time-vector[,[current-vector]]. One period only
    :param sample_factor: f_sampling/f_period, defaults to 1000
    :param plot: insert anything else than "no" or 'False' to show a plot to visualize input and output
    :param rad: 'no' for time domain input vector, anything else than 'no' for 2pi-time domain
    :param f0: set when rad != 'no' and rad != False
    :param title: plot window title, defaults to 'ffT'
    :return: numpy-array [[frequency-vector],[amplitude-vector],[phase-vector]]
    """

    t = period_vector_t_i[0]
    i = period_vector_t_i[1]

    if rad != 'no' and rad!=False:
        if f0 is None:
            raise ValueError("if rad!='no', a fundamental frequency f0 must be set")
        else:
            period_vector_t_i[0] = period_vector_t_i[0] / (2 * np.pi * f0)

    # time domain
    t_interp = np.linspace(0, t[-1], sample_factor)
    i_interp = np.interp(t_interp, t, i)

    f0 = round(1 / t[-1])
    Fs = f0 * sample_factor

    # frequency domain
    f = np.linspace(0, (sample_factor - 1) * f0, sample_factor)
    x = np.fft.fft(i_interp)
    x_mag = np.abs(x) / sample_factor
    phi_rad = np.angle(x)

    f_corrected = f[0:int(sample_factor / 2 + 1)]
    x_mag_corrected = 2 * x_mag[0:int(sample_factor / 2 + 1)]
    x_mag_corrected[0] = x_mag_corrected[0] / 2
    phi_rad_corrected = phi_rad[0:int(sample_factor / 2 + 1)]

    f_out = []
    x_out = []
    phi_rad_out = []
    for count, value in enumerate(x_mag_corrected):
        if x_mag_corrected[count] > 0.01 * max(i):
            f_out.append(f_corrected[count])
            x_out.append(x_mag_corrected[count])
            phi_rad_out.append(phi_rad_corrected[count])

    if plot != 'no' and plot != False:
        print(f"{title = }")
        print(f"{t[-1] = }")
        print(f"{f0 = }")
        print(f"{Fs = }")
        print(f"{sample_factor = }")
        print(f"f_out = {np.around(f_out, 0)}")
        print(f"x_out = {np.around(x_out, 1)}")
        print(f"phi_rad_out = {np.around(phi_rad_out, 1)}")

        reconstructed_signal = 0
        for i_range in range(len(f_out)):
            reconstructed_signal += x_out[i_range] * np.cos(
                2 * np.pi * f_out[i_range] * t_interp + phi_rad_out[i_range])

        fig, [ax1, ax2] = plt.subplots(num=title, nrows=2, ncols=1)
        ax1.plot(t, i, label='original signal')
        ax1.plot(t_interp, reconstructed_signal, label='reconstructed signal')
        ax1.grid()
        ax1.set_title('Signal')
        ax1.set_xlabel('time in s')
        ax1.set_ylabel('Amplitude')
        ax1.legend()
        ax2.stem(f_out, x_out)
        ax2.grid()
        ax2.set_title('ffT')
        ax2.set_xlabel('Frequency in Hz')
        ax2.set_ylabel('Amplitude')
        plt.tight_layout()
        plt.show()

    return np.array([f_out, x_out, phi_rad_out])


def compare_fft_list(list: list, rad: float = 'no', f0: Union[float,None] = None) -> None:
    """
    generate fft curves from input curves and compare them to each other

    minimal example:
    example_waveform = np.array([[0, 1.34, 3.14, 4.48, 6.28],[-175.69, 103.47, 175.69, -103.47,-175.69]])
    example_waveform_2 = np.array([[0, 0.55, 3.14, 3.69, 6.28],[-138.37, 257.58, 138.37, -257.58, -138.37]])
    compare_fft_list([example_waveform, example_waveform_2], rad='yes', f0=25000)

    :param list: list of fft-compatible numpy-arrays [element, element, ... ], each element format like [[time-vector[,[current-vector]]. One period only
    :param rad: 'no' for time domain input vector, anything else than 'no' for 2pi-time domain
    :param f0: set when rad != 'no'
    :return: plot
    """

    out = []
    for count, value in enumerate(list):
        out.append([fft(list[count], sample_factor=1000, plot='no', rad=rad, f0=f0)])

    fig, axs = plt.subplots(2, len(list), sharey=True)
    for count, value in enumerate(list):
        axs[0, count].plot(list[count][0], list[count][1], label='original signal')
        axs[0, count].grid()
        axs[0, count].set_xlabel('time in s')
        axs[0, count].set_ylabel('Amplitude')
        axs[1, count].stem(out[count][0][0], out[count][0][1])
        axs[1, count].grid()
        axs[1, count].set_xlabel('frequency in Hz')
        axs[1, count].set_ylabel('Amplitude')

        # ax1.plot(t_interp, reconstructed_signal, label='reconstructed signal')

    plt.tight_layout()
    plt.show()


# ----------------------------------------------------------------------------------------------------------------------
# Reluctance Model [with calculation]
mu0 = 4e-7*np.pi


def r_basis(l, w, h):
    """
    1-dim reluctance per-unit-of-length
    [according to "A Novel Approach for 3D Air Gap Reluctance Calculations" - J. Mühlethaler, J.W. Kolar, A. Ecklebe]
    :param w:
    :param l:
    :param h:
    :return:
    """
    if l <= 0:
        l = 0.0000001
    return 1 / (mu0 * (w/2/l + 2/np.pi * (1+np.log(np.pi*h/4/l))))


def sigma(l, w, R_equivalent):
    """
    1-dim fringing factor
    [according to "A Novel Approach for 3D Air Gap Reluctance Calculations" - J. Mühlethaler, J.W. Kolar, A. Ecklebe]
    :param w:
    :param l:
    :param R_equivalent:
    :return:
    """
    return R_equivalent / (l/mu0/w)


def r_round_inf(l, sigma, r):
    """
    3-dim reluctance for 2-dim axial-symmetric air gaps
    :param sigma:
    :param r:
    :param l:
    :return:
    """
    return sigma**2 * l/mu0/r**2/np.pi


def r_round_round(l, sigma, r):
    """
    3-dim reluctance for 2-dim axial-symmetric air gaps
    :param sigma:
    :param r:
    :param l:
    :return:
    """
    return sigma**2 * l/mu0/r**2/np.pi


def r_cyl_cyl(l, sigma, w, r_o):
    """

    :param l:
    :param sigma:
    :param w:
    :param r_o:
    :return:
    """
    return sigma * np.log(r_o/(r_o-l)) / 2/mu0/np.pi/w


def r_cheap_cyl_cyl(r_o, l, w):
    """

    :param r_o:
    :param l:
    :param w:
    :return:
    """
    r_i = r_o - l
    return (r_o-r_i) / mu0/w/np.pi/(r_o+r_i)

