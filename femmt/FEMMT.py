# Usual Python libraries
import csv
import fileinput
import numpy as np
import os
import pandas as pd
import pathlib
import re
import sys
import time
import warnings
from matplotlib import pyplot as plt
# FEM and Mesh interfaces
import femm
import gmsh
from onelab import onelab
# Self written functions
from functions import id_generator, inner_points, min_max_inner_points, call_for_path, NbrStrands


# Master Class
class MagneticComponent:
    def __init__(self, conductor_type="litz"):

        # ==============================
        # Settings
        # ==============================

        # -- Geometry Control Flags --
        self.y_symmetric = 1  # Mirror-symmetry across y-axis
        self.axi_symmetric = 1  # Axial-symmetric model (idealized full-cylindrical)
        self.s = 0.5  # Parameter for mesh-accuracy

        # -- Geometry --

        # - Core -
        self.core_type = "EI"  # Basic shape of magnetic conductor
        self.n_air_gaps = 1  # Number of air gaps [==1: air gap in center | >1: random air gaps]

        # - Dimensions -
        self.core_cond_isolation = 0.001  # gap between Core and Conductor
        self.cond_cond_isolation = 0.0002  # gap between Core and Conductor
        self.core_w = 0.02  # Axi symmetric case | core_w := core radius
        self.window_w = 0.01  # Winding window width
        self.window_h = 0.03  # Winding window height

        # - Conductor -
        self.n_conductors = 33  # Number of (homogenised) conductors in one window
        self.conductor_type = conductor_type  # Stranded wires
        """
        conductor_type = "stacked"  # Vertical packing of conductors
        conductor_type = "full"  # One massive Conductor in each window
        conductor_type = "foil"  # Horizontal packing of conductors
        conductor_type = "solid" # Massive wires
        conductor_type = "litz" # Litz wires
        """

        if self.conductor_type == 'solid':
            self.conductor_radius = 0.0011879
            self.A_cell = np.pi * self.conductor_radius**2  # * self.FF  # Surface of the litz approximated hexagonal cell

        if self.conductor_type == 'litz':
            self.update_litz_configuration(conductor_radius=0.0012, n_layers=10, strand_radius=0.05e-3)
            #self.A_cell = np.pi * self.conductor_radius**2  # * self.FF  # Surface of the litz approximated hexagonal cell


        # -- Materials --
        # frequency = 0: mu_rel only used if flag_non_linear_core == 0
        # frequency > 0: mu_rel is used
        self.mu0 = 4e-7*np.pi
        self.mu_rel = 3000   # relative Core Permeability
        self.core_material = 95  # 95 := TDK-N95 | Currently only works with Numbers corresponding to BH.pro
        self.sigma = 5.8e7

        # -- Characteristic lengths -- [for mesh sizes]
        self.c_core = self.core_w/10. * self.s
        self.c_window = self.window_w/10 * self.s
        self.c_conductor = self.window_w/30 * self.s

        # -- Parent folder path --
        self.path = str(pathlib.Path(__file__).parent.absolute())

        # -- Further Attributes --
        self.c_airgap = None
        self.n_windows = None
        self.p_outer = None
        self.p_window = None
        self.p_conductor = None
        self.p_air_gaps = None

        self.flag_imposed_reduced_frequency = None
        self.flag_excitation_type = None
        self.flag_non_linear_core = None
        self.current = None
        self.voltage = None
        self.frequency = None
        self.red_freq = None
        self.delta = None

        self.onelab = None
        self.sweep_frequencies = None

        # FEMM variables
        self.tot_loss_femm = None

    # ==== Back-End Methods =====
    def update_litz_configuration(self, litz_parametrization_type='implicite_FF', strand_radius=None, ff=None, conductor_radius=None, n_layers=None):
        """
        - used to change litz configurations
        - also used at first initialisation of the geometry
        - needed to always make sure that the relation between litz parameters (strand radius, fill factor, number of
          layers/strands and conductor/litz radius) is valid and consistent
        - 4 parameters, 1 degree of freedom (dof)
        :param litz_parametrization_type:
        :param strand_radius:
        :param ff:
        :param conductor_radius:
        :param n_layers:
        :return:
        """
        # Litz Approximation
        if litz_parametrization_type == 'implicite_FF':
            self.conductor_radius = conductor_radius
            self.n_layers = n_layers
            self.n_strands = NbrStrands(self.n_layers)
            self.strand_radius = strand_radius
            ff_exact = self.n_strands*self.strand_radius**2/self.conductor_radius**2 # hexagonal packing: ~90.7% are theoretical maximum
            print(f"Exact fill factor: {ff_exact}")
            self.FF = np.around(ff_exact, decimals=2)
            print(f"Updated Litz Configuration: \n"
                  f"FF: {self.FF} \n"
                  f"Number of layers/strands: {self.n_layers}/{self.n_strands} \n"
                  f"Strand radius: {self.strand_radius}"
                  f"Conductor radius: {self.conductor_radius}")
            #print(f"Rounded fill factor: {self.FF}")

        if litz_parametrization_type == 'implicite_litz_radius':
            self.FF = ff
            self.n_layers = n_layers
            self.n_strands = NbrStrands(self.n_layers)
            self.strand_radius = strand_radius
            self.conductor_radius = np.sqrt(self.n_strands*self.strand_radius**2/self.FF)
            print(f"Updated Litz Configuration: \n "
                  f"FF: {self.FF} \n "
                  f"Number of layers/strands: {self.n_layers}/{self.n_strands} \n"
                  f"Strand radius: {self.strand_radius}"
                  f"Conductor radius: {self.conductor_radius}")

        self.A_cell = self.n_strands * self.strand_radius**2 * np.pi / self.FF

    def onelab_setup(self):
        """
        Either reads onelab parent folder path from config.py or asks the user to provide it.
        Creates a config.py at first use.
        :return:
        """
        if os.path.isfile(self.path + "/config.py"):
            import config
            with open('config.py') as f:
                if 'onelab' in f.read():
                    if os.path.exists(config.onelab):
                        self.onelab = config.onelab
                else:
                    self.onelab = call_for_path("onelab")
        else:
            self.onelab = call_for_path("onelab")

    def high_level_geo_gen(self, core_type="EI", axi_symmetric=1):
        """
        - high level geometry generation
        - based on chosen core and conductor types and simulation mode
        - calls "low level" methods, that creates all points needed for mesh generation
        :return:
        """
        # ==============================
        # High-Level Geometry Generation
        # ==============================

        # -- Core-type --
        if self.core_type == core_type:
            self.n_windows = 2

        # -- Symmetry -- [Choose between asymmetric, symmetric and axi symmetric]
        if self.y_symmetric == 0:
            self.axi_symmetric = 0
        if self.axi_symmetric == 1:
            self.y_symmetric = axi_symmetric

        if self.core_type == "EI" and self.axi_symmetric == 1:
            self.ei_axi()

    def ei_axi(self):
        """
        - creates all points needed for the radial axi-symetric EI core typology
        :return:
        """
        # -- Air Gap Data -- [random air gap generation]
        air_gaps = np.empty((self.n_air_gaps, 4))
        i = 0
        while i in range(0, self.n_air_gaps):
            position_tag = 0  # '-1': left leg | '0': center leg | '1': right leg
            if self.n_air_gaps == 1:
                airgap_h = 0.001
                airgap_position = 0.5*(self.window_h-airgap_h)-(self.window_h/2-airgap_h/2)
            else:
                airgap_h = np.random.rand(1)*0.005 + 0.001
                airgap_position = np.random.rand(1)*(self.window_h-airgap_h)-(self.window_h/2-airgap_h/2)
            self.c_airgap = airgap_h / 3 * self.s
            # Overlapping Control
            for j in range(0, air_gaps.shape[0]):
                if position_tag == air_gaps[j, 0] and air_gaps[j, 1] + air_gaps[j, 2] / 2 > airgap_position > air_gaps[j, 1] - air_gaps[j, 2] / 2:
                    print("Overlapping air Gaps have been corrected")
            else:
                air_gaps[i, :] = np.array([position_tag, airgap_position, airgap_h, self.c_airgap])
                i += 1

        # -- Arrays for geometry data -- [all points with (x, y, z, mesh_accuracy)]
        self.p_outer = np.zeros((4, 4))
        self.p_window = np.zeros((4*self.n_windows, 4))
        self.p_air_gaps = np.zeros((4*self.n_air_gaps, 4))

        # -- Geometry data --

        """
        if self.y_symmetric == 0:
            # Outer
            self.p_outer[0][:] = [-(self.core_w + self.window_w), -(self.window_h / 2 + self.core_w), 0, self.c_core]
            self.p_outer[1][:] = [self.core_w + self.window_w, -(self.window_h / 2 + self.core_w), 0, self.c_core]
            self.p_outer[2][:] = [-(self.core_w + self.window_w), (self.window_h / 2 + self.core_w), 0, self.c_core]
            self.p_outer[3][:] = [self.core_w + self.window_w, (self.window_h / 2 + self.core_w), 0, self.c_core]
            # Window
            self.p_window[0] = [-(self.core_w/2+self.window_w), -self.window_h/2, 0, self.c_window]
            self.p_window[1] = [-self.core_w/2, -self.window_h/2, 0, self.c_window]
            self.p_window[2] = [-(self.core_w/2+self.window_w)/2, self.window_h, 0, self.c_window]
            self.p_window[3] = [-self.core_w/2, self.window_h/2, 0, self.c_window]
            self.p_window[4] = [self.core_w/2, -self.window_h/2, 0, self.c_window]
            self.p_window[5] = [(self.core_w/2+self.window_w), -self.window_h/2, 0, self.c_window]
            self.p_window[6] = [self.core_w/2, self.window_h/2, 0, self.c_window]
            self.p_window[7] = [(self.core_w/2+self.window_w), self.window_h/2, 0, self.c_window]
        """


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
        self.p_window[0] = [-r_inner, -self.window_h/2, 0, self.c_window]
        self.p_window[1] = [-self.core_w/2, -self.window_h/2, 0, self.c_window]
        self.p_window[2] = [-r_inner, self.window_h/2, 0, self.c_window]
        self.p_window[3] = [-self.core_w/2, self.window_h/2, 0, self.c_window]
        self.p_window[4] = [self.core_w/2, -self.window_h/2, 0, self.c_window]
        self.p_window[5] = [r_inner, -self.window_h/2, 0, self.c_window]
        self.p_window[6] = [self.core_w/2, self.window_h/2, 0, self.c_window]
        self.p_window[7] = [r_inner, self.window_h/2, 0, self.c_window]

        # - Conductors -

        # Case: no conductors [only theoretical]
        self.p_conductor = np.empty(0)

        if self.conductor_type == "full":
            # full window conductor
            self.p_conductor[0][:] = [self.core_cond_isolation + self.core_w/2, -self.window_h/2 + self.core_cond_isolation, 0, self.c_conductor]
            self.p_conductor[1][:] = [r_inner - self.core_cond_isolation, -self.window_h/2 + self.core_cond_isolation, 0, self.c_conductor]
            self.p_conductor[2][:] = [self.core_cond_isolation + self.core_w/2, self.window_h/2 - self.core_cond_isolation, 0, self.c_conductor]
            self.p_conductor[3][:] = [r_inner - self.core_cond_isolation, self.window_h/2 - self.core_cond_isolation, 0, self.c_conductor]

        if self.conductor_type == "stacked":
            # stacking from the ground
            self.p_conductor = np.empty((4*self.n_conductors, 4))
            for i in range(0, self.n_conductors):
                # two conductors above
                self.p_conductor[4*i+0][:] = [self.core_cond_isolation + self.core_w/2, (1-i)*self.core_cond_isolation + i*(-self.window_h/2 + self.core_cond_isolation), 0, self.c_conductor]
                self.p_conductor[4*i+1][:] = [r_inner - self.core_cond_isolation, (1-i)*self.core_cond_isolation + i*(-self.window_h/2 + self.core_cond_isolation), 0, self.c_conductor]
                self.p_conductor[4*i+2][:] = [self.core_cond_isolation + self.core_w/2, -i*self.core_cond_isolation + (1-i)*(self.window_h/2 - self.core_cond_isolation), 0, self.c_conductor]
                self.p_conductor[4*i+3][:] = [r_inner - self.core_cond_isolation, -i*self.core_cond_isolation + (1-i)*(self.window_h/2 - self.core_cond_isolation), 0, self.c_conductor]

        if self.conductor_type == "foil":
            self.p_conductor = np.empty((4*self.n_conductors, 4))
            left_bound = self.core_cond_isolation + self.core_w/2
            right_bound = r_inner - self.core_cond_isolation
            x_interpol = np.linspace(left_bound, right_bound, self.n_conductors+1)  # instead should FF window from inside to outside with fixed copper thickness
            for i in range(0, self.n_conductors):
                # Foils
                self.p_conductor[4 * i + 0][:] = [x_interpol[i] + self.cond_cond_isolation, -self.window_h / 2 + self.core_cond_isolation, 0, self.c_conductor]
                self.p_conductor[4 * i + 1][:] = [x_interpol[i+1] - self.cond_cond_isolation, -self.window_h / 2 + self.core_cond_isolation, 0, self.c_conductor]
                self.p_conductor[4 * i + 2][:] = [x_interpol[i] + self.cond_cond_isolation, self.window_h / 2 - self.core_cond_isolation, 0, self.c_conductor]
                self.p_conductor[4 * i + 3][:] = [x_interpol[i+1] - self.cond_cond_isolation, self.window_h / 2 - self.core_cond_isolation, 0, self.c_conductor]

        if self.conductor_type == "litz" or self.conductor_type == "solid":
            self.p_conductor = []  # center points are stored
            left_bound = self.core_w/2
            right_bound = r_inner - self.core_cond_isolation
            top_bound = self.window_h/2
            bot_bound = -self.window_h/2

            y = bot_bound + self.core_cond_isolation + self.conductor_radius
            x = left_bound + self.core_cond_isolation + self.conductor_radius
            i = 0
            # Case n_conductors higher that "allowed" is missing
            while y < top_bound-self.core_cond_isolation-self.conductor_radius and i < self.n_conductors:
                while x < right_bound-self.core_cond_isolation-self.conductor_radius and i < self.n_conductors:
                    self.p_conductor.append([x, y, 0, self.c_conductor])
                    self.p_conductor.append([x-self.conductor_radius, y, 0, self.c_conductor])
                    self.p_conductor.append([x, y+self.conductor_radius, 0, self.c_conductor])
                    self.p_conductor.append([x+self.conductor_radius, y, 0, self.c_conductor])
                    self.p_conductor.append([x, y-self.conductor_radius, 0, self.c_conductor])
                    i += 1
                    x += self.conductor_radius * 2 + self.cond_cond_isolation
                y += self.conductor_radius * 2 + self.cond_cond_isolation
                x = left_bound + self.core_cond_isolation + self.conductor_radius
            self.p_conductor = np.asarray(self.p_conductor)
            if int(self.p_conductor.shape[0]/5) < self.n_conductors:
                print("Could not resolve all conductors.")
                self.n_conductors = int(self.p_conductor.shape[0]/5)

        for i in range(0, self.n_air_gaps):
            # Left leg (-1)
            if air_gaps[i][0] == -1:
                self.p_air_gaps[i * 4] = [-(self.core_w + self.window_w), air_gaps[i][1] - air_gaps[i][2] / 2, 0, air_gaps[i][3]]
                self.p_air_gaps[i * 4 + 1] = [-(self.core_w / 2 + self.window_w), air_gaps[i][1] - air_gaps[i][2] / 2, 0, air_gaps[i][3]]
                self.p_air_gaps[i * 4 + 2] = [-(self.core_w + self.window_w), air_gaps[i][1] + air_gaps[i][2] / 2, 0, air_gaps[i][3]]
                self.p_air_gaps[i * 4 + 3] = [-(self.core_w / 2 + self.window_w), air_gaps[i][1] + air_gaps[i][2] / 2, 0, air_gaps[i][3]]

            # Center leg (0)
            if air_gaps[i][0] == 0:
                self.p_air_gaps[i * 4] = [-self.core_w/2, air_gaps[i][1] - air_gaps[i][2] / 2, 0, air_gaps[i][3]]
                self.p_air_gaps[i * 4 + 1] = [self.core_w/2, air_gaps[i][1] - air_gaps[i][2] / 2, 0, air_gaps[i][3]]
                self.p_air_gaps[i * 4 + 2] = [-self.core_w/2, air_gaps[i][1] + air_gaps[i][2] / 2, 0, air_gaps[i][3]]
                self.p_air_gaps[i * 4 + 3] = [self.core_w/2, air_gaps[i][1] + air_gaps[i][2] / 2, 0, air_gaps[i][3]]

            # Right leg (+1)
            if air_gaps[i][0] == 1:
                self.p_air_gaps[i * 4] = [self.core_w / 2 + self.window_w, air_gaps[i][1] - air_gaps[i][2] / 2, 0, air_gaps[i][3]]
                self.p_air_gaps[i * 4 + 1] = [self.core_w + self.window_w, air_gaps[i][1] - air_gaps[i][2] / 2, 0, air_gaps[i][3]]
                self.p_air_gaps[i * 4 + 2] = [self.core_w / 2 + self.window_w, air_gaps[i][1] + air_gaps[i][2] / 2, 0, air_gaps[i][3]]
                self.p_air_gaps[i * 4 + 3] = [self.core_w + self.window_w, air_gaps[i][1] + air_gaps[i][2] / 2, 0, air_gaps[i][3]]

    def generate_mesh(self):
        """
        - interaction with gmsh
        - mesh generation
        :return:
        """
        # Initialization
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.model.add("geometry")
        # ------------------------------------------ Geometry -------------------------------------------
        # Core generation
        if self.core_type == "EI":
            # --------------------------------------- Points --------------------------------------------
            if self.y_symmetric == 1:
                if self.axi_symmetric == 1:

                    # Find points of air gaps (used later)
                    if self.n_air_gaps > 0:
                        center_right = min_max_inner_points(self.p_window[4], self.p_window[6], self.p_air_gaps)
                        island_right = inner_points(self.p_window[4], self.p_window[6], self.p_air_gaps)

                    # Pre-Definitions
                    # Points
                    p_core = []
                    p_island = []
                    p_cond = []
                    # Curves
                    l_bound_core = []
                    l_bound_air = []
                    l_core_air = []
                    l_cond = []
                    curve_loop_cond = []
                    # Curve Loops
                    curve_loop_island = []
                    curve_loop_air = []
                    curve_loop_bound = []
                    # Plane Surfaces
                    plane_surface_core = []
                    plane_surface_cond = []
                    plane_surface_air = []

                    # =====================
                    # Main Core
                    # Points of Main Core (index refers to sketch)
                    if self.n_air_gaps > 0:
                        p_core.append(
                            gmsh.model.geo.addPoint(0, center_right[0][1], center_right[0][2], center_right[0][3]))
                    if self.n_air_gaps == 0:
                        p_core.append(None)  # dummy filled for no air gap special case
                    p_core.append(gmsh.model.geo.addPoint(0, self.p_outer[1][1], self.p_outer[1][2], self.p_outer[1][3]))
                    p_core.append(gmsh.model.geo.addPoint(self.p_outer[1][0], self.p_outer[1][1], self.p_outer[1][2],
                                                          self.p_outer[1][3]))
                    p_core.append(gmsh.model.geo.addPoint(self.p_outer[3][0], self.p_outer[3][1], self.p_outer[3][2],
                                                          self.p_outer[3][3]))
                    p_core.append(gmsh.model.geo.addPoint(0, self.p_outer[3][1], self.p_outer[3][2], self.p_outer[3][3]))
                    if self.n_air_gaps > 0:
                        p_core.append(
                            gmsh.model.geo.addPoint(0, center_right[1][1], center_right[1][2], center_right[1][3]))
                        p_core.append(
                            gmsh.model.geo.addPoint(center_right[1][0], center_right[1][1], center_right[1][2],
                                                    center_right[1][3]))
                    if self.n_air_gaps == 0:
                        p_core.append(None)  # dummy filled for no air gap special case
                        p_core.append(None)  # dummy filled for no air gap special case
                    p_core.append(gmsh.model.geo.addPoint(self.p_window[6][0], self.p_window[6][1], self.p_window[6][2],
                                                          self.p_window[6][3]))
                    p_core.append(gmsh.model.geo.addPoint(self.p_window[7][0], self.p_window[7][1], self.p_window[7][2],
                                                          self.p_window[7][3]))
                    p_core.append(gmsh.model.geo.addPoint(self.p_window[5][0], self.p_window[5][1], self.p_window[5][2],
                                                          self.p_window[5][3]))
                    p_core.append(gmsh.model.geo.addPoint(self.p_window[4][0], self.p_window[4][1], self.p_window[4][2],
                                                          self.p_window[4][3]))
                    if self.n_air_gaps > 0:
                        p_core.append(
                            gmsh.model.geo.addPoint(center_right[0][0], center_right[0][1], center_right[0][2],
                                                    center_right[0][3]))
                    if self.n_air_gaps == 0:
                        p_core.append(None)  # dummy filled for no air gap special case
                    # Curves of Main Core (index refers to sketch)
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
                            min_island_right = np.argmin(island_right[:, 1])
                            p_island.append(gmsh.model.geo.addPoint(0, island_right[min_island_right, 1],
                                                                    island_right[min_island_right, 2],
                                                                    island_right[min_island_right, 3]))
                            p_island.append(gmsh.model.geo.addPoint(island_right[min_island_right, 0],
                                                                    island_right[min_island_right, 1],
                                                                    island_right[min_island_right, 2],
                                                                    island_right[min_island_right, 3]))
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
                    for i in range(0, self.p_conductor.shape[0]):
                        p_cond.append(
                            gmsh.model.geo.addPoint(self.p_conductor[i][0], self.p_conductor[i][1], 0, self.c_conductor))
                    # Curves of Conductors
                    if self.conductor_type == "litz" or self.conductor_type == "solid":
                        for i in range(0, int(len(p_cond) / 5)):
                            l_cond.append(
                                gmsh.model.geo.addCircleArc(p_cond[5 * i + 1], p_cond[5 * i + 0], p_cond[5 * i + 2]))
                            l_cond.append(
                                gmsh.model.geo.addCircleArc(p_cond[5 * i + 2], p_cond[5 * i + 0], p_cond[5 * i + 3]))
                            l_cond.append(
                                gmsh.model.geo.addCircleArc(p_cond[5 * i + 3], p_cond[5 * i + 0], p_cond[5 * i + 4]))
                            l_cond.append(
                                gmsh.model.geo.addCircleArc(p_cond[5 * i + 4], p_cond[5 * i + 0], p_cond[5 * i + 1]))
                            # Iterative plane creation
                            curve_loop_cond.append(gmsh.model.geo.addCurveLoop(
                                [l_cond[i * 4 + 0], l_cond[i * 4 + 1], l_cond[i * 4 + 2], l_cond[i * 4 + 3]]))
                            plane_surface_cond.append(gmsh.model.geo.addPlaneSurface([curve_loop_cond[i]]))
                    else:
                        for i in range(0, int(len(p_cond) / 4)):
                            l_cond.append(gmsh.model.geo.addLine(p_cond[4 * i + 0], p_cond[4 * i + 2]))
                            l_cond.append(gmsh.model.geo.addLine(p_cond[4 * i + 2], p_cond[4 * i + 3]))
                            l_cond.append(gmsh.model.geo.addLine(p_cond[4 * i + 3], p_cond[4 * i + 1]))
                            l_cond.append(gmsh.model.geo.addLine(p_cond[4 * i + 1], p_cond[4 * i + 0]))
                            # Iterative plane creation
                            curve_loop_cond.append(gmsh.model.geo.addCurveLoop(
                                [l_cond[i * 4 + 0], l_cond[i * 4 + 1], l_cond[i * 4 + 2], l_cond[i * 4 + 3]]))
                            plane_surface_cond.append(gmsh.model.geo.addPlaneSurface([curve_loop_cond[i]]))
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
                    plane_surface_air.append(gmsh.model.geo.addPlaneSurface(curve_loop_air + curve_loop_cond))

                    # =====================
                    # Bound
                    l_bound_tmp = l_bound_core[:5]
                    for i in range(0, len(l_bound_air)):
                        l_bound_tmp.append(l_bound_air[-i - 1])
                        if i != len(l_bound_air) - 1:  # last run
                            l_bound_tmp.append(l_bound_core[-i - 1])
                    # curve_loop_bound.append(gmsh.model.geo.addCurveLoop(l_bound_tmp, reorient=True))

        # Define physical Surfaces and Curves
        # Core
        ps_core = gmsh.model.geo.addPhysicalGroup(2, plane_surface_core, tag=2000)
        # Conductors
        ps_cond = []
        if self.n_conductors == 2 and self.conductor_type == "stacked":
            ps_cond[0] = gmsh.model.geo.addPhysicalGroup(2, [plane_surface_cond[0]], tag=4000)
            ps_cond[1] = gmsh.model.geo.addPhysicalGroup(2, [plane_surface_cond[1]], tag=4001)
        if self.conductor_type == "foil" or self.conductor_type == "solid":
            for i in range(0, self.n_conductors):
                ps_cond.append(gmsh.model.geo.addPhysicalGroup(2, [plane_surface_cond[i]], tag=4000 + i))
        if self.conductor_type == "litz":
            for i in range(0, self.n_conductors):
                ps_cond.append(gmsh.model.geo.addPhysicalGroup(2, [plane_surface_cond[i]], tag=6000 + i))
        # Air
        ps_air = gmsh.model.geo.addPhysicalGroup(2, plane_surface_air, tag=1000)
        # Boundary
        pc_bound = gmsh.model.geo.addPhysicalGroup(1, l_bound_tmp, tag=1111)

        gmsh.model.setPhysicalName(2, ps_core, "CORE")
        gmsh.model.setPhysicalName(2, ps_cond[0], "COND1")
        gmsh.model.setPhysicalName(2, ps_air, "AIR")
        gmsh.model.setPhysicalName(1, pc_bound, "BOUND")

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
        for i in range(0, len(plane_surface_cond)):
            gmsh.model.setColor([(2, plane_surface_cond[i])], 150, 150, 0)

        # Mesh generation
        gmsh.model.mesh.generate(2)

        # Check operating system
        if sys.platform == "linux" or sys.platform == "linux2":
            gmsh.write(self.path + "/geometry.msh")
        elif sys.platform == "darwin":
            # OS X
            gmsh.write(self.path + "/geometry.msh")
        elif sys.platform == "win32":
            gmsh.write(self.path + "/geometry.msh")  # Win10 can handle slash

        # Open gmsh GUI for visualization
        # gmsh.fltk.run()

        # Terminate gmsh
        gmsh.finalize()

    def excitation(self, f, i, nonlinear=0, ex_type='current', imposed_red_f=0):
        """
        - excitation of the electromagnetic problem
        - current, voltage or current density
        - frequency or reduced frequency
        :param f:
        :param i:
        :param nonlinear:
        :param ex_type:
        :param imposed_red_f:
        :return:
        """
        # -- Excitation --
        self.flag_imposed_reduced_frequency = imposed_red_f  # if == 0 --> impose frequency f
        self.flag_excitation_type = ex_type  # 'current', 'current_density', 'voltage'
        self.flag_non_linear_core = nonlinear

        # Imposed current, current density or voltage
        if self.flag_excitation_type == 'current':
            self.current = i
        if self.flag_excitation_type == 'current_density':
            raise NotImplementedError
        if self.flag_excitation_type == 'voltage':
            raise NotImplementedError
            # self.voltage = 2

        # -- Frequency --
        self.frequency = f  # in Hz
        if self.flag_imposed_reduced_frequency == 1:
            self.red_freq = 4
        else:
            if self.frequency != 0:
                self.delta = np.sqrt(2 / (2 * self.frequency * np.pi * self.sigma * self.mu0))
                self.red_freq = self.strand_radius / self.delta
            else:
                self.delta = 1e20  # random huge value
                self.red_freq = 0

    def file_communication(self):
        """
        Interaction between python and Prolog files.
        :return:
        """
        # ------------------------------------- File Communication ----------------------------------
        # All shared control variables and parameters are passed to a temporary Prolog file
        text_file = open("Parameter.pro", "w")

        # -- Control Flags --
        if self.flag_excitation_type == 'current':
            text_file.write(f"Flag_ImposedVoltage = 0;\n")
        if self.flag_excitation_type == 'voltage':
            text_file.write(f"Flag_ImposedVoltage = 1;\n")
        if self.conductor_type == 'litz':
            text_file.write(f"Flag_HomogenisedModel = 1;\n")
        else:
            text_file.write(f"Flag_HomogenisedModel = 0;\n")
        text_file.write("Flag_imposedRr = %s;\n" % self.flag_imposed_reduced_frequency)
        # -- Geometry --
        # Number of conductors
        text_file.write(f"NbrCond = {self.n_conductors};\n")
        # For stranded Conductors:
        text_file.write(f"NbrstrandedCond = {self.n_conductors};\n")  # redundant
        text_file.write(f"NbrStrands = {self.n_strands};\n")
        text_file.write(f"Rc = {self.conductor_radius};\n")
        text_file.write(f"Fill = {self.FF};\n")
        text_file.write(f"NbrLayers = {self.n_layers};\n")
        #text_file.write(f"NbrLayers = 0;\n")
        #text_file.write(f"AreaCell = {self.A_cell};\n")
        text_file.write(f"AreaCell = {self.A_cell};\n")
        """
        # Coordinates of the rectangular winding window
        if self.axi_symmetric == 1:
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
        text_file.write(f"SigmaCu = {self.sigma};\n")

        # Core Material
        if self.frequency == 0:
            if self.flag_non_linear_core == 1:
                text_file.write(f"Flag_NL = 1;\n")
                text_file.write(f"Core_Material = {self.core_material};\n")
            else:
                text_file.write(f"Flag_NL = 0;\n")
                text_file.write(f"mur = {self.mu_rel};\n")
        if self.frequency != 0:
            text_file.write(f"Flag_NL = 0;\n")
            text_file.write(f"mur = {self.mu_rel};\n")

        # -- Excitation --
        # Imposed current, current density or voltage
        if self.flag_excitation_type == 'current':
            text_file.write(f"Val_EE = {self.current};\n")
        if self.flag_excitation_type == 'current_density':
            text_file.write(f"Val_EE = {self.current_density};\n")
        if self.flag_excitation_type == 'voltage':
            text_file.write(f"Val_EE = {self.voltage};\n")

        # Frequency and reduced Frequency
        text_file.write("Freq = %s;\n" % self.frequency)
        text_file.write(f"delta = {self.delta};\n")
        text_file.write(f"Rr = {self.red_freq};\n")

        text_file.close()

    def simulate(self):
        """
        Initializes a onelab client. Provides the GetDP based solver with the created mesh file.
        :return:
        """
        self.onelab_setup()
        # -- Simulation --
        # create a new onelab client
        c = onelab.client(__file__)

        # get model file names with correct path
        msh_file = c.getPath('geometry.msh')
        solver = c.getPath('ind_axi_python_controlled' + '.pro')

        # Run simulations as sub clients (non blocking??)
        mygetdp = self.onelab + 'getdp'
        c.runSubClient('myGetDP', mygetdp + ' ' + solver + ' -msh ' + msh_file + ' -solve Analysis -v2')

    def visualize(self):
        """
        - a post simulation viewer
        - allows to open ".pos"-files in gmsh
        - For example current density, ohmic losses or the magnetic field density can be visualized
        :return:
        """
        # ---------------------------------------- Visualization in gmsh ---------------------------------------
        gmsh.initialize()
        epsilon = 1e-9
        # Mesh
        gmsh.option.setNumber("Mesh.SurfaceEdges", 0)

        if self.conductor_type != 'litz':
            # Ohmic losses (weightend effective value of current density)
            gmsh.open("res/j2F.pos")
            gmsh.option.setNumber("View[0].ScaleType", 2)
            gmsh.option.setNumber("View[0].RangeType", 2)
            gmsh.option.setNumber("View[0].SaturateValues", 1)
            gmsh.option.setNumber("View[0].CustomMin", gmsh.option.getNumber("View[0].Min") + epsilon)
            gmsh.option.setNumber("View[0].CustomMax", gmsh.option.getNumber("View[0].Max"))
            gmsh.option.setNumber("View[0].ColormapNumber", 1)
            gmsh.option.setNumber("View[0].IntervalsType", 2)
            gmsh.option.setNumber("View[0].NbIso", 40)

            # Magnetic flux density
            gmsh.open("res/Magb.pos")
            gmsh.option.setNumber("View[1].ScaleType", 1)
            gmsh.option.setNumber("View[1].RangeType", 1)
            gmsh.option.setNumber("View[1].CustomMin", gmsh.option.getNumber("View[1].Min") + epsilon)
            gmsh.option.setNumber("View[1].CustomMax", gmsh.option.getNumber("View[1].Max"))
            gmsh.option.setNumber("View[1].ColormapNumber", 1)
            gmsh.option.setNumber("View[1].IntervalsType", 2)
            gmsh.option.setNumber("View[1].NbIso", 40)

            print(gmsh.option.getNumber("View[0].Max"))

        if self.conductor_type == 'litz':
            # Ohmic losses (weightend effective value of current density)
            gmsh.open("res/jH.pos")
            gmsh.option.setNumber("View[0].ScaleType", 1)
            # gmsh.option.setNumber("View[0].RangeType", 2)
            # gmsh.option.setNumber("View[0].SaturateValues", 1)
            # gmsh.option.setNumber("View[0].CustomMin", gmsh.option.getNumber("View[0].Min") + epsilon)
            # gmsh.option.setNumber("View[0].CustomMax", gmsh.option.getNumber("View[0].Max"))
            gmsh.option.setNumber("View[0].ColormapNumber", 1)
            gmsh.option.setNumber("View[0].IntervalsType", 2)
            gmsh.option.setNumber("View[0].NbIso", 40)

            # Magnetic flux density
            gmsh.open("res/MagbH.pos")
            gmsh.option.setNumber("View[1].ScaleType", 1)
            gmsh.option.setNumber("View[1].RangeType", 1)
            gmsh.option.setNumber("View[1].CustomMin", gmsh.option.getNumber("View[1].Min") + epsilon)
            gmsh.option.setNumber("View[1].CustomMax", gmsh.option.getNumber("View[1].Max"))
            gmsh.option.setNumber("View[1].ColormapNumber", 1)
            gmsh.option.setNumber("View[1].IntervalsType", 2)
            gmsh.option.setNumber("View[1].NbIso", 40)

        gmsh.fltk.run()
        gmsh.finalize()

    def femm_reference(self, freq, current, sigma, non_visualize=0):
        """
        Allows reference simulations with the 2D open source electromagnetic FEM tool FEMM.
        Helpful to validate changes (especially in the Prolog Code).
        :param non_visualize:
        :param freq:
        :param current:
        :param sigma:
        :return:
        """

        # == Pre Geometry ==
        self.high_level_geo_gen()
        self.ei_axi()

        if self.n_air_gaps != 1:
            raise NotImplementedError

        # == Init ==
        femm.openfemm(non_visualize)
        femm.newdocument(0)
        femm.mi_probdef(freq, 'meters', 'axi', 1.e-8, 0, 30)

        # == Materials ==
        femm.mi_addmaterial('Ferrite', 3000, 3000, 0, 0, 0, 0, 0, 1, 0, 0, 0)
        femm.mi_addmaterial('Air', 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0)
        if self.conductor_type == "litz":
            femm.mi_addmaterial('Copper', 1, 1, 0, 0, sigma, 0, 0, 1, 5, 0, 0, self.n_strands, 2*1000*self.strand_radius)  # type := 5. last argument
        if self.conductor_type == "solid":
            femm.mi_addmaterial('Copper', 1, 1, 0, 0, sigma, 0, 0, 1, 0, 0, 0, 0, 0)

        print(f"Strandsnumber: {self.n_strands}")
        print(f"Strandsdiameter in mm: {2*1000*self.strand_radius}")

        # == Circuit ==
        # coil as seen from the terminals.
        femm.mi_addcircprop('icoil', current, 1)

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
        if self.conductor_type == "litz" or self.conductor_type == "solid":
            for i in range(0, int(self.p_conductor.shape[0] / 5)):
                # 0: center | 1: left | 2: top | 3: right | 4.bottom
                femm.mi_drawarc(self.p_conductor[5*i+1][0], self.p_conductor[5*i+1][1], self.p_conductor[5*i+3][0], self.p_conductor[5*i+3][1], 180, 2.5)
                femm.mi_addarc(self.p_conductor[5*i+3][0], self.p_conductor[5*i+3][1], self.p_conductor[5*i+1][0], self.p_conductor[5*i+1][1],  180, 2.5)
                femm.mi_addblocklabel(self.p_conductor[5*i][0], self.p_conductor[5*i][1])
                femm.mi_selectlabel(self.p_conductor[5*i][0], self.p_conductor[5*i][1])
                femm.mi_setblockprop('Copper', 1, 1, 'icoil', 0, 0, 1)
                femm.mi_clearselected

        # Define an "open" boundary condition using the built-in function:
        femm.mi_makeABC()

        # == Labels/Designations ==

        # Label for core
        femm.mi_addblocklabel(self.p_outer[3, 0]-0.001, self.p_outer[3, 1]-0.001)
        femm.mi_selectlabel(self.p_outer[3, 0]-0.001, self.p_outer[3, 1]-0.001)
        femm.mi_setblockprop('Ferrite', 1, 1, '<None>', 0, 0, 0)
        femm.mi_clearselected()

        # Labels for air
        femm.mi_addblocklabel(0.001, 0)
        femm.mi_selectlabel(0.001, 0)
        femm.mi_setblockprop('Air', 1, 1, '<None>', 0, 0, 0)
        femm.mi_clearselected()
        femm.mi_addblocklabel(self.p_outer[3, 0]+0.001, self.p_outer[3, 1]+0.001)
        femm.mi_selectlabel(self.p_outer[3, 0]+0.001, self.p_outer[3, 1]+0.001)
        femm.mi_setblockprop('Air', 1, 1, '<None>', 0, 0, 0)
        femm.mi_clearselected()

        # Now, the finished input geometry can be displayed.
        femm.mi_zoomnatural()
        femm.mi_saveas('coil.fem')
        femm.mi_analyze()
        femm.mi_loadsolution()

        # == Losses ==
        tmp = femm.mo_getcircuitproperties('icoil')
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

    def data_logging(self, sim_choice):
        """
        !!! not ready implemented !!!

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
          # pseudo 2D dataframe: ['strand number, strand radius'], 'frequency'  | const. litz radius --> FF
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
        with open(self.path + '/res/' + loss_file, newline='') as f:
            reader = csv.reader(f)
            data = list(reader)
        return data[-last_n_values:-1] + [data[-1]]

    def pre_simulate(self):
        """
        Used to determine the litz-approximation coefficients.
        :return:
        """
        if self.conductor_type == 'litz':

            if os.path.isfile(self.path + f"/pre/coeff/pB_RS_la{self.FF}_{self.n_layers}layer.dat"):
            #if os.path.isfile(self.path + f"/pre/coeff/pB_RS_la{self.FF}_0layer.dat"):
                print("Coefficients for stands approximation are found.")

            else:
                print("Create coefficients for strands approximation")
                # need reduced frequency X for pre-simulation
                # Rounding X to fit it with corresponding parameters from the database
                X = self.red_freq
                X = np.around(X, decimals=3)
                print(f"Rounded Reduced frequency X = {X}")

                # Create new file with [0 1] for the first element entry

                # create a new onelab client
                # -- Pre-Simulation Settings --
                text_file = open("pre/PreParameter.pro", "w")
                text_file.write(f"NbrLayers = {self.n_layers};\n")
                #text_file.write(f"NbrLayers = 0;\n")
                text_file.write(f"Fill = {self.FF};\n")
                print("Here")
                text_file.write(f"Rc = {self.strand_radius};\n")  # double named!!! must be changed
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
                        text_file = open("pre/PreParameter.pro", "w")
                        text_file.write(f"Rr_cell = {rf};\n")
                        text_file.write(f"Mode = {mode};\n")
                        text_file.write(f"NbrLayers = {self.n_layers};\n")
                        #text_file.write(f"NbrLayers = 0;\n")
                        text_file.write(f"Fill = {self.FF};\n")
                        text_file.write(f"Rc = {self.strand_radius};\n")  # double named!!! must be changed
                        text_file.close()

                        # get model file names with correct path
                        input_file = c.getPath('pre/cell_dat.pro')
                        cell = c.getPath('pre/cell.pro')

                        # Run simulations as sub clients
                        mygetdp = self.onelab + 'getdp'
                        c.runSubClient('myGetDP', mygetdp + ' ' + cell + ' -input ' + input_file + ' -solve MagDyn_a -v2')

                # Formatting stuff
                files = [self.path + f"/pre/coeff/pB_RS_la{self.FF}_{self.n_layers}layer.dat",
                         self.path + f"/pre/coeff/pI_RS_la{self.FF}_{self.n_layers}layer.dat",
                         self.path + f"/pre/coeff/qB_RS_la{self.FF}_{self.n_layers}layer.dat",
                         self.path + f"/pre/coeff/qI_RS_la{self.FF}_{self.n_layers}layer.dat"]
                #files = [self.path + f"/pre/coeff/pB_RS_la{self.FF}_0layer.dat",
                #         self.path + f"/pre/coeff/pI_RS_la{self.FF}_0layer.dat",
                #         self.path + f"/pre/coeff/qB_RS_la{self.FF}_0layer.dat",
                #         self.path + f"/pre/coeff/qI_RS_la{self.FF}_0layer.dat"]
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
        self.excitation(f=100000, i=1)  # frequency and current
        self.file_communication()
        self.pre_simulate()

    # ==== Front-End Methods =====
    def single_simulation(self, freq, current):
        """
        - can be used for a single simulation
        - no sweeping at all
        :param freq:
        :param current:
        :return:
        """
        self.high_level_geo_gen()
        self.generate_mesh()
        self.excitation(f=freq, i=current)  # frequency and current
        self.file_communication()
        self.pre_simulate()
        self.simulate()
        self.visualize()

    def mesh(self):
        self.high_level_geo_gen()
        self.generate_mesh()

    def excitation_sweep(self, frequencies=[], currents=[], show_last=False):
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
        if len(frequencies) != len(currents):
            print('len(frequencies) != len(currents)')
            raise Exception
        for i in range(0, len(frequencies)):
            self.excitation(f=frequencies[i], i=currents[i])  # frequency and current
            self.file_communication()
            self.pre_simulate()
            self.simulate()
        if show_last:
            self.visualize()

    def litz_loss_comparison(self, FF, n_layers, strand_radius, sim_choice, sweep_parameter='fill_factor', nametag=''):
        """

        :param sweep_parameter:
        :param FF:
        :param n_layers:
        :param strand_radius:
        :param sim_choice:
        :return:
        """
        # Excitation: fixed sweep parameters
        frequencies = np.linspace(0, 250000, 6)
        currents = [2, 2, 2, 2, 2, 2]

        # Find Json logfile
        target_femm = self.path + nametag + 'Pv_FEMM_' + sweep_parameter + '.json'
        target_femmt = self.path + nametag + 'Pv_FEMMT_' + sweep_parameter + '.json'

        # Update model to chosen litz parameters
        self.update_litz_configuration(litz_parametrization_type='implicite_litz_radius', n_layers=n_layers,
                                       strand_radius=strand_radius, ff=FF)
        if self.FF < 0.4 or self.FF > 0.9:
            print(f"Skip simulation with non realistic fill factor {self.FF}")

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
            femmt_col_tag = f"onelab, {self.n_layers}, {self.strand_radius}, {self.FF}"
            femm_col_tag = f"femm, {self.n_layers}, {self.strand_radius}, {self.FF}"

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
            error.insert(loc=0, column=f"FF{fillfactors[-col - 1]}",
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
        df_Pv_femmt.plot(ax=ax, marker="o")
        # ax.text(25000, 1, 'FEMM', color='r', ha='right', rotation=0, wrap=True)
        # ax.text(25000, 0.5, 'FEMMT', color='g', ha='right', rotation=0, wrap=True)
        plt.show()