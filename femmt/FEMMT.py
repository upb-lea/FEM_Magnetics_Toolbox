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
import gmsh
# import adapt_mesh
from onelab import onelab
import json
import random
import string
# Self written functions. It is necessary to write a . before the function, du to handling this package also as a pip-package
#from .femmt_functions import id_generator, inner_points, min_max_inner_points, call_for_path, NbrStrands
# FEM and Mesh interfaces, only for windows machines
if os.name == 'nt':
    import femm


# Master Class
class MagneticComponent:
    def __init__(self, component_type="inductor"):
        """

        :param component_type:
        """
        # ==============================
        # Settings
        # ==============================
        # -- Parent folder path --
        self.path = str(pathlib.Path(__file__).parent.absolute())  # Path of FEMMT files
        self.onelab = None  # Path to Onelab installation folder

        # ==============================
        # Geometry
        # ==============================
        # -- Control Flags --
        self.y_symmetric = 1  # Mirror-symmetry across y-axis
        self.dimensionality = "2D axi"  # Axial-symmetric model (idealized full-cylindrical)
        self.s = 1  # Parameter for mesh-accuracy
        self.component_type = component_type  # "inductor" or "transformer"

        # -- Core --
        self.core_type = "EI"  # Basic shape of magnetic conductor
        self.core_w = None  # Axi symmetric case | core_w := core radius
        self.window_w = None  # Winding window width
        self.window_h = None  # Winding window height
        self.update_core(core_type="EI", core_w=0.02, window_w=0.01, window_h=0.03)  # some initial values

        # -- Air gaps --
        self.n_air_gaps = 1  # Number of air gaps [==1: air gap in center | >1: random air gaps]
        self.air_gaps = np.empty((self.n_air_gaps, 4))  # list with [position_tag, air_gap_position, air_gap_h, c_air_gap]

        # -- Isolation ---
        self.core_cond_isolation = 0.001  # gap between Core and Conductors
        self.cond_cond_isolation = 0.0002  # gap between two Conductors

        # -- Conductor --
        self.n_conductors = None  # Number of (homogenised) conductors in one window
        if component_type == "inductor":
            self.n_conductors = 1
        if component_type == "transformer":
            self.n_conductors = 2
        self.conductor_type = [None] * self.n_conductors  # List of possible conductor types
        self.turns = [None] * self.n_conductors
        self.FF = [None] * self.n_conductors
        self.n_layers = [None] * self.n_conductors
        self.n_strands = [None] * self.n_conductors
        self.strand_radius = [None] * self.n_conductors
        self.conductor_radius = [None] * self.n_conductors
        self.A_cell = [None] * self.n_conductors
        self.update_conductors(n_turns=[None] * self.n_conductors,
                               conductor_type=[None] * self.n_conductors,
                               conductor_radix=[None] * self.n_conductors,
                               layer_numbers=[None] * self.n_conductors,
                               strand_radix=[None] * self.n_conductors)

        # -- Geometric Parameters/Coordinates --
        self.n_windows = None
        self.p_outer = None
        self.p_window = None
        self.p_conductor = [[], []]
        self.p_air_gaps = None

        # ==============================
        # Materials
        # ==============================
        # frequency = 0: mu_rel only used if flag_non_linear_core == 0
        # frequency > 0: mu_rel is used
        self.mu0 = 4e-7*np.pi
        self.mu_rel = 3000   # relative Core Permeability
        self.core_material = 95  # 95 := TDK-N95 | Currently only works with Numbers corresponding to BH.pro
        self.sigma = 5.8e7

        # ==============================
        # Problem Definition
        # ==============================
        # -- Excitation Parameters --
        self.flag_imposed_reduced_frequency = None
        self.flag_excitation_type = None
        self.flag_non_linear_core = None
        self.current = [None] * self.n_conductors  # Defined for every conductor
        self.current_density = [None] * self.n_conductors  # Defined for every conductor
        self.voltage = [None] * self.n_conductors  # Defined for every conductor
        self.frequency = None
        self.phase_tmp = np.zeros(self.n_conductors)  # Default is zero, Defined for every conductor
        self.red_freq = [None] * self.n_conductors  # Defined for every conductor
        self.delta = None

        # ==============================
        # Meshing
        # ==============================
        # -- Characteristic lengths -- [for mesh sizes]
        self.skin_mesh_factor = None
        self.c_core = self.core_w/10. * self.s
        self.c_window = self.window_w/10 * self.s
        self.c_conductor = [None] * self.n_conductors  # self.delta  # self.s /20 #self.window_w/30 * self.s
        self.c_center_conductor = [None] * self.n_conductors  # used for the mesh accuracy in the conductors
        self.c_air_gap = []

        # -- Used for Litz Validation --
        self.sweep_frequencies = None

        # -- FEMM variables --
        self.tot_loss_femm = None

    # ==== Back-End Methods =====
    def update_core(self, core_type,  **kwargs):
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
        self.n_air_gaps = n_air_gaps
        self.air_gaps = np.empty((self.n_air_gaps, 4))
        self.c_air_gap = [None] * self.n_air_gaps

        """
        # Optional updating the number of air gaps
        for key, value in kwargs.items():
            if key == 'n_air_gaps':
                self.n_air_gaps = value
        """

        # Update air gaps with chosen method
        if method == "center" and self.dimensionality == "2D axi":
            if self.n_air_gaps > 1:
                print(f"{self.n_air_gaps} are too many air gaps for the 'center' option!")
                raise Warning
            else:
                self.c_air_gap[0] = air_gap_h[0] / 3 * self.s
                self.air_gaps[0, :] = np.array([0, 0, air_gap_h[0], self.c_air_gap[0]])

        if method == "random" and self.dimensionality == "2D axi":
            position_tag = [0] * self.n_air_gaps

            i = 0
            while i in range(0, self.n_air_gaps):
                height = np.random.rand(1) * 0.001 + 0.001
                position = np.random.rand(1) * (self.window_h - height) - (self.window_h / 2 - height / 2)
                self.c_air_gap[i] = height / 3 * self.s
                # Overlapping Control
                for j in range(0, self.air_gaps.shape[0]):
                    if self.air_gaps[j, 1]+self.air_gaps[j, 2]/2 > position > self.air_gaps[j, 1]-self.air_gaps[j, 2]/2:
                        if position_tag[i] == self.air_gaps[j, 0]:
                            print(f"Overlapping air Gaps have been corrected")
                else:
                    self.air_gaps[i, :] = np.array([position_tag[i], position, height, self.c_air_gap[i]])
                    i += 1

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
                            raise Warning
                else:
                    self.c_air_gap[i] = air_gap_h[i] / 3 * self.s
                    self.air_gaps[i, :] = np.array([position_tag[i], air_gap_position[i], air_gap_h[i], self.c_air_gap[i]])

    def update_conductors(self, n_turns=[], conductor_type=[], conductor_radix=[], layer_numbers=[], strand_radix=[]):
        """
        conductor_type = "stacked"  # Vertical packing of conductors
        conductor_type = "full"     # One massive Conductor in each window
        conductor_type = "foil"     # Horizontal packing of conductors
        conductor_type = "solid"    # Massive wires
        conductor_type = "litz"     # Litz wires
        :param n_turns:
        :param strand_radix:
        :param layer_numbers:
        :param conductor_radix:
        :param conductor_radius:
        :param n_conductors:
        :param conductor_type:
        :return:
        """
        if self.component_type == "inductor":
            if len(n_turns) != 1 or len(conductor_type) != 1:
                print(f"Wrong number of conductor parameters passed to inductor model!")
                raise Warning
        if self.component_type == "transformer":
            if len(n_turns) != 2 or len(conductor_type) != 2:
                print(f"Wrong number of conductor parameters passed to transformer model!")
                raise Warning

        self.turns = n_turns
        self.conductor_type = conductor_type

        for i in range(0, self.n_conductors):
            print(i)
            if self.conductor_type[i] == 'solid':
                self.conductor_radius[i] = conductor_radix[i]
                self.A_cell[i] = np.pi * self.conductor_radius[i]**2  # Area of the litz approximated hexagonal cell
            if self.conductor_type[i] == 'litz':
                self.update_litz_configuration(num=i,
                                               litz_parametrization_type='implicite_FF',
                                               ff=None,
                                               conductor_radius=conductor_radix[i],
                                               n_layers=layer_numbers[i],
                                               strand_radius=strand_radix[i])
                # Surface of the litz approximated hexagonal cell
                # self.A_cell = np.pi * self.conductor_radius**2  # * self.FF

    def update_litz_configuration(self, num=0, litz_parametrization_type='implicite_FF', strand_radius=None, ff=None,
                                  conductor_radius=None, n_layers=None):
        """
        - updates the conductor with number num (simple transformer: num=0 -> primary winding,
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
        :param n_layers:
        :return:
        """
        # Litz Approximation
        if litz_parametrization_type == 'implicite_FF':
            self.conductor_radius[num] = conductor_radius
            self.n_layers[num] = n_layers
            self.n_strands[num] = NbrStrands(self.n_layers[num])
            self.strand_radius[num] = strand_radius
            # hexagonal packing: ~90.7% are theoretical maximum FF
            ff_exact = self.n_strands[num]*self.strand_radius[num]**2/self.conductor_radius[num]**2
            print(f"Exact fill factor: {ff_exact}")
            self.FF[num] = np.around(ff_exact, decimals=2)
            print(f"Updated Litz Configuration: \n"
                  f"FF: {self.FF[num]} \n"
                  f"Number of layers/strands: {self.n_layers[num]}/{self.n_strands[num]} \n"
                  f"Strand radius: {self.strand_radius[num]} \n"
                  f"Conductor radius: {self.conductor_radius[num]}")
            # print(f"Rounded fill factor: {self.FF}")

        if litz_parametrization_type == 'implicite_litz_radius':
            self.FF[num] = ff
            self.n_layers[num] = n_layers
            self.n_strands[num] = NbrStrands(self.n_layers)
            self.strand_radius[num] = strand_radius
            self.conductor_radius[num] = np.sqrt(self.n_strands[num]*self.strand_radius[num]**2/self.FF[num])
            print(f"Updated Litz Configuration: \n "
                  f"FF: {self.FF[num]} \n "
                  f"Number of layers/strands: {self.n_layers[num]}/{self.n_strands[num]} \n"
                  f"Strand radius: {self.strand_radius[num]} \n"
                  f"Conductor radius: {self.conductor_radius[num]}")

        self.A_cell[num] = self.n_strands[num] * self.strand_radius[num]**2 * np.pi / self.FF[num]

    def onelab_setup(self):
        """
        Either reads onelab parent folder path from config.p or asks the user to provide it.
        Creates a config.p at first run.
        :return: -
        """
        if os.path.isfile(self.path + "/config.json") and os.stat(self.path + "/config.json") != 0:
            json_file = open(self.path + '/config.json','rb') #with open(self.path + '/config.p') as f:
            loaded_dict = json.load(json_file)
            json_file.close()
            path = loaded_dict['onelab']

            if os.path.exists(path):
                self.onelab = path
            else:
                self.onelab = call_for_path("onelab")
        else:
            self.onelab = call_for_path("onelab")

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

        # Update Skin Depth (needed for meshing)
        self.skin_mesh_factor = skin_mesh_factor
        if frequency != None:
            if frequency == 0:
                self.delta = 1e9
            else:
                self.delta = np.sqrt(2 / (2 * frequency * np.pi * self.sigma * self.mu0))
            for i in range(0, len(self.conductor_radius)):
                self.c_conductor[i] = min([self.delta * self.skin_mesh_factor,
                                    self.conductor_radius[i] / 4 * self.skin_mesh_factor])
        for i in range(0, len(self.conductor_radius)):
            self.c_center_conductor[i] = self.conductor_radius[i] / 4 * self.skin_mesh_factor

        print(f"Werte Leiter: {self.c_conductor, self.c_center_conductor}")

        # -- Core-type --
        if self.core_type == core_type:
            self.n_windows = 2

        # -- Symmetry -- [Choose between asymmetric, symmetric and axi symmetric]
        self.dimensionality = dimensionality
        if self.dimensionality == "2D axi":
            self.y_symmetric = 1

        if self.core_type == "EI" and self.dimensionality == "2D axi":
            self.ei_axi()

    def ei_axi(self):
        """
        - creates all points needed for the radial axi-symmetric EI core typology
        :return:
        """
        # -- Air Gap Data -- [random air gap generation]
        #self.update_air_gaps()

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
        for num in range(0, self.n_conductors):

            # Case: no conductors [only theoretical]
            # self.p_conductor = np.empty((num, 0))
            """
            if self.conductor_type == "full":
                # full window conductor
                self.p_conductor[0][:] = [self.core_cond_isolation + self.core_w/2, -self.window_h/2 + self.core_cond_isolation, 0, self.c_conductor]
                self.p_conductor[1][:] = [r_inner - self.core_cond_isolation, -self.window_h/2 + self.core_cond_isolation, 0, self.c_conductor]
                self.p_conductor[2][:] = [self.core_cond_isolation + self.core_w/2, self.window_h/2 - self.core_cond_isolation, 0, self.c_conductor]
                self.p_conductor[3][:] = [r_inner - self.core_cond_isolation, self.window_h/2 - self.core_cond_isolation, 0, self.c_conductor]
    
            if self.conductor_type == "stacked":
                # stacking from the ground
                self.p_conductor = np.empty((4*self.turns, 4))
                for i in range(0, self.turns):
                    # two conductors above
                    self.p_conductor[4*i+0][:] = [self.core_cond_isolation + self.core_w/2, (1-i)*self.core_cond_isolation + i*(-self.window_h/2 + self.core_cond_isolation), 0, self.c_conductor]
                    self.p_conductor[4*i+1][:] = [r_inner - self.core_cond_isolation, (1-i)*self.core_cond_isolation + i*(-self.window_h/2 + self.core_cond_isolation), 0, self.c_conductor]
                    self.p_conductor[4*i+2][:] = [self.core_cond_isolation + self.core_w/2, -i*self.core_cond_isolation + (1-i)*(self.window_h/2 - self.core_cond_isolation), 0, self.c_conductor]
                    self.p_conductor[4*i+3][:] = [r_inner - self.core_cond_isolation, -i*self.core_cond_isolation + (1-i)*(self.window_h/2 - self.core_cond_isolation), 0, self.c_conductor]
    
            if self.conductor_type == "foil":
                self.p_conductor = np.empty((4*self.turns, 4))
                left_bound = self.core_cond_isolation + self.core_w/2
                right_bound = r_inner - self.core_cond_isolation
                x_interpol = np.linspace(left_bound, right_bound, self.turns+1)  # instead should FF window from inside to outside with fixed copper thickness
                for i in range(0, self.turns):
                    # Foils
                    self.p_conductor[4 * i + 0][:] = [x_interpol[i] + self.cond_cond_isolation, -self.window_h / 2 + self.core_cond_isolation, 0, self.c_conductor]
                    self.p_conductor[4 * i + 1][:] = [x_interpol[i+1] - self.cond_cond_isolation, -self.window_h / 2 + self.core_cond_isolation, 0, self.c_conductor]
                    self.p_conductor[4 * i + 2][:] = [x_interpol[i] + self.cond_cond_isolation, self.window_h / 2 - self.core_cond_isolation, 0, self.c_conductor]
                    self.p_conductor[4 * i + 3][:] = [x_interpol[i+1] - self.cond_cond_isolation, self.window_h / 2 - self.core_cond_isolation, 0, self.c_conductor]
            """

            if self.conductor_type[num] == "litz" or self.conductor_type[num] == "solid":
                left_bound = self.core_w/2
                right_bound = r_inner - self.core_cond_isolation

                if self.component_type == "transformer":
                    # xfmr: oben prim - unten sek
                    if num == 0:
                        top_bound = self.window_h/2
                        bot_bound = 0
                    if num == 1:
                        top_bound = 0
                        bot_bound = -self.window_h/2

                if self.component_type == "inductor":
                    top_bound = self.window_h/2
                    bot_bound = -self.window_h/2


                y = bot_bound + self.core_cond_isolation + self.conductor_radius[num]
                x = left_bound + self.core_cond_isolation + self.conductor_radius[num]
                i = 0
                # Case n_conductors higher that "allowed" is missing
                while y < top_bound-self.core_cond_isolation-self.conductor_radius[num] and i < self.turns[num]:
                    while x < right_bound-self.core_cond_isolation-self.conductor_radius[num] and i < self.turns[num]:
                        self.p_conductor[num].append([x, y, 0, self.c_center_conductor[num]])
                        self.p_conductor[num].append([x-self.conductor_radius[num], y, 0, self.c_conductor[num]])
                        self.p_conductor[num].append([x, y+self.conductor_radius[num], 0, self.c_conductor[num]])
                        self.p_conductor[num].append([x+self.conductor_radius[num], y, 0, self.c_conductor[num]])
                        self.p_conductor[num].append([x, y-self.conductor_radius[num], 0, self.c_conductor[num]])
                        i += 1
                        x += self.conductor_radius[num] * 2 + self.cond_cond_isolation
                    y += self.conductor_radius[num] * 2 + self.cond_cond_isolation
                    x = left_bound + self.core_cond_isolation + self.conductor_radius[num]
                self.p_conductor[num] = np.asarray(self.p_conductor[num])
                if int(self.p_conductor[num].shape[0]/5) < self.turns[num]:
                    print("Could not resolve all conductors.")
                    self.turns[num] = int(self.p_conductor[num].shape[0]/5)
                print(f"Conductors: {self.p_conductor}")
        # - Air gaps -
        # "air_gaps" is a list with [position_tag, air_gap_position, air_gap_h, c_air_gap]
        #   - position_tag: specifies the gapped "leg"
        #   - air_gap_position: specifies the coordinate of the air gap's center point along the specified leg
        #   - air_gap_h: height/length of the air gap
        #   - c_air_gap: mesh accuracy factor
        # at this point the 4 corner points of each air gap are generated out of "air_gaps"
        for i in range(0, self.n_air_gaps):
            # Left leg (-1)
            if self.air_gaps[i][0] == -1:
                self.p_air_gaps[i * 4] = [-(self.core_w + self.window_w), self.air_gaps[i][1] - self.air_gaps[i][2] / 2, 0, self.air_gaps[i][3]]
                self.p_air_gaps[i * 4 + 1] = [-(self.core_w / 2 + self.window_w), self.air_gaps[i][1] - self.air_gaps[i][2] / 2, 0, self.air_gaps[i][3]]
                self.p_air_gaps[i * 4 + 2] = [-(self.core_w + self.window_w), self.air_gaps[i][1] + self.air_gaps[i][2] / 2, 0, self.air_gaps[i][3]]
                self.p_air_gaps[i * 4 + 3] = [-(self.core_w / 2 + self.window_w), self.air_gaps[i][1] + self.air_gaps[i][2] / 2, 0, self.air_gaps[i][3]]

            # Center leg (0)
            if self.air_gaps[i][0] == 0:
                self.p_air_gaps[i * 4] = [-self.core_w/2, self.air_gaps[i][1] - self.air_gaps[i][2] / 2, 0, self.air_gaps[i][3]]
                self.p_air_gaps[i * 4 + 1] = [self.core_w/2, self.air_gaps[i][1] - self.air_gaps[i][2] / 2, 0, self.air_gaps[i][3]]
                self.p_air_gaps[i * 4 + 2] = [-self.core_w/2, self.air_gaps[i][1] + self.air_gaps[i][2] / 2, 0, self.air_gaps[i][3]]
                self.p_air_gaps[i * 4 + 3] = [self.core_w/2, self.air_gaps[i][1] + self.air_gaps[i][2] / 2, 0, self.air_gaps[i][3]]

            # Right leg (+1)
            if self.air_gaps[i][0] == 1:
                self.p_air_gaps[i * 4] = [self.core_w / 2 + self.window_w, self.air_gaps[i][1] - self.air_gaps[i][2] / 2, 0, self.air_gaps[i][3]]
                self.p_air_gaps[i * 4 + 1] = [self.core_w + self.window_w, self.air_gaps[i][1] - self.air_gaps[i][2] / 2, 0, self.air_gaps[i][3]]
                self.p_air_gaps[i * 4 + 2] = [self.core_w / 2 + self.window_w, self.air_gaps[i][1] + self.air_gaps[i][2] / 2, 0, self.air_gaps[i][3]]
                self.p_air_gaps[i * 4 + 3] = [self.core_w + self.window_w, self.air_gaps[i][1] + self.air_gaps[i][2] / 2, 0, self.air_gaps[i][3]]

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
            gmsh.model.mesh.generate(2)
            # --------------------------------------
            # Mesh generation
            #gmsh.model.mesh.generate(2)
            # Check operating system
            if sys.platform == "linux" or sys.platform == "linux2":
                gmsh.write(self.path + "/geometry.msh")
            elif sys.platform == "darwin":
                # OS X
                gmsh.write(self.path + "/geometry.msh")
            elif sys.platform == "win32":
                gmsh.write(self.path + "/geometry.msh")  # Win10 can handle slash

            # Terminate gmsh
            gmsh.finalize()

    def find_neighbours(self, file="res/J_rms.pos"):
        # Open loss/error results
        dest_file = open(self.path + "error.dat", "w")
        # Read the logged losses corresponding to the frequencies
        with open(self.path + '/res/J_rms.pos') as f:
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

    def alternative_local_error(self, loss_file='/res/J_rms.pos'):
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
        with open(self.path + '/res/J_rms.pos') as file:
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
        with open(self.path + '/res/J_rms.pos') as file:
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

    def local_error(self, loss_file='/res/error.pos'):
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

    def create_background_mesh(self):
        gmsh.open(self.path + "/geometry.msh")  # Open current mesh
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
        gmsh.view.addModelData(err_view, 0, "geometry", "ElementData", mesh.triangles_tags, local_error[:, None])
        gmsh.view.write(err_view, "err.pos")

        # Refinement
        sf_ele = self.compute_size_field(mesh.vxyz, mesh.triangles, local_error, N)
        np.savetxt("sf_ele.txt", sf_ele)
        sf_view = gmsh.view.add("mesh size field")
        gmsh.view.addModelData(sf_view, 0, "geometry", "ElementData", mesh.triangles_tags, sf_ele[:, None])
        gmsh.view.write(sf_view, "sf.pos")

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

            self.create_background_mesh()

        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.model.add("geometry")
        # ------------------------------------------ Geometry -------------------------------------------
        # Core generation
        if self.core_type == "EI":
            # --------------------------------------- Points --------------------------------------------
            if self.y_symmetric == 1:
                if self.dimensionality == "2D axi":

                    # Find points of air gaps (used later)
                    if self.n_air_gaps > 0:
                        center_right = min_max_inner_points(self.p_window[4], self.p_window[6], self.p_air_gaps)
                        island_right = inner_points(self.p_window[4], self.p_window[6], self.p_air_gaps)

                    # Pre-Definitions
                    # Points
                    p_core = []
                    p_island = []
                    p_cond = [[], []]
                    # Curves
                    l_bound_core = []
                    l_bound_air = []
                    l_core_air = []
                    l_cond = [[], []]
                    curve_loop_cond = [[], []]
                    # Curve Loops
                    curve_loop_island = []
                    curve_loop_air = []
                    curve_loop_bound = []
                    # Plane Surfaces
                    plane_surface_core = []
                    plane_surface_cond = [[], []]
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
                    for num in range(0, self.n_conductors):
                        for i in range(0, self.p_conductor[num].shape[0]):
                            p_cond[num].append(gmsh.model.geo.addPoint(self.p_conductor[num][i][0],
                                                                  self.p_conductor[num][i][1], 0, self.p_conductor[num][i][3]))
                        # Curves of Conductors
                        if self.conductor_type[num] == "litz" or self.conductor_type[num] == "solid":
                            for i in range(0, int(len(p_cond[num]) / 5)):
                                l_cond[num].append(
                                    gmsh.model.geo.addCircleArc(p_cond[num][5 * i + 1], p_cond[num][5 * i + 0], p_cond[num][5 * i + 2]))
                                l_cond[num].append(
                                    gmsh.model.geo.addCircleArc(p_cond[num][5 * i + 2], p_cond[num][5 * i + 0], p_cond[num][5 * i + 3]))
                                l_cond[num].append(
                                    gmsh.model.geo.addCircleArc(p_cond[num][5 * i + 3], p_cond[num][5 * i + 0], p_cond[num][5 * i + 4]))
                                l_cond[num].append(
                                    gmsh.model.geo.addCircleArc(p_cond[num][5 * i + 4], p_cond[num][5 * i + 0], p_cond[num][5 * i + 1]))
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
        ps_cond = [[], []]  # xfmr
        for num in range(0, self.n_conductors):
            if self.turns == 2 and self.conductor_type == "stacked":
                ps_cond[0] = gmsh.model.geo.addPhysicalGroup(2, [plane_surface_cond[num][0]], tag=4000)  # ??? xfmr
                ps_cond[1] = gmsh.model.geo.addPhysicalGroup(2, [plane_surface_cond[num][1]], tag=4001)
            if self.conductor_type[num] == "foil" or self.conductor_type[num] == "solid":
                for i in range(0, self.turns[num]):
                    ps_cond[num].append(gmsh.model.geo.addPhysicalGroup(2, [plane_surface_cond[num][i]], tag=4000 + 1000*num + i))
            if self.conductor_type[num] == "litz":
                for i in range(0, self.turns[num]):
                    ps_cond[num].append(gmsh.model.geo.addPhysicalGroup(2, [plane_surface_cond[num][i]], tag=6000 + 1000*num + i))

        # Air
        ps_air = gmsh.model.geo.addPhysicalGroup(2, plane_surface_air, tag=1000)
        # Boundary
        pc_bound = gmsh.model.geo.addPhysicalGroup(1, l_bound_tmp, tag=1111)
        print(f"Physical Conductor Surfaces: {ps_cond}")  # xfmr
        gmsh.model.setPhysicalName(2, ps_core, "CORE")
        for num in range(0, self.n_conductors):
            for i in range(0, len(ps_cond[num])):
                gmsh.model.setPhysicalName(2, ps_cond[num][i], f"COND{num+1}")
        gmsh.model.setPhysicalName(2, ps_air, "AIR")
        gmsh.model.setPhysicalName(1, pc_bound, "BOUND")

        # Remove Points from Model
        #for i in range(9,39):
        #    gmsh.model.geo.remove(dimTags=[(0, i)])
        print(f"P_cond: {p_cond}")

        # Forward Meshing
        # Inter Conductors
        for num in range(0, self.n_conductors):
            p_inter = []
            x_inter = []
            y_inter = []
            j = 0

            if self.turns[num] > 1:
                while self.p_conductor[num][5*j][1] == self.p_conductor[num][5*j+5][1]:
                    x_inter.append(0.5*(self.p_conductor[num][5*j][0]+self.p_conductor[num][5*j+5][0]))
                    j += 1
                    if j == self.turns[num]-1:
                        break
                j += 1
                print(f"j = {j}")
                if int(self.turns[num]/j) > 1:
                    for i in range(0, int(self.turns[num]/j)):
                        if 5*j*i+5*j >= len(self.p_conductor[num][:]):
                            break
                        y_inter.append(0.5*(self.p_conductor[num][5*j*i][1]+self.p_conductor[num][5*j*i+5*j][1]))
                    for x in x_inter:
                        for y in y_inter:
                            p_inter.append(gmsh.model.geo.addPoint(x, y, 0, self.c_center_conductor[num]))
            print(f"x_inter = {x_inter}")
            print(f"y_inter = {y_inter}")
            print(f"p_inter = {p_inter}")

        # Synchronize
        gmsh.model.geo.synchronize()
        # Conductor Center
        for num in range(0, self.n_conductors):
            for i in range(0, int(len(p_cond[num]) / 5)):
                gmsh.model.mesh.embed(0, [p_cond[num][5 * i + 0]], 2, plane_surface_cond[num][i])

        # Inter Conductors
        gmsh.model.mesh.embed(0, p_inter, 2, plane_surface_air[0])

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
        for num in range(0, self.n_conductors):
            for i in range(0, len(plane_surface_cond[num])):
                gmsh.model.setColor([(2, plane_surface_cond[num][i])], 150, 150, 0)
        # -----------------------------------------
        if refine == 1:
            print("\n ------- \nRefined Mesh Creation ")
            # mesh the new gmsh.model using the size field
            bg_field = gmsh.model.mesh.field.add("PostView")
            gmsh.model.mesh.field.setNumber(bg_field, "ViewTag", sf_view)
            gmsh.model.mesh.field.setAsBackgroundMesh(bg_field)
            gmsh.model.mesh.generate(2)
            gmsh.write("geometry.msh")
        else:
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

    def excitation(self, f, i, phases=[], nonlinear=0, ex_type='current', imposed_red_f=0, ):
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
        # -- Excitation --
        self.flag_imposed_reduced_frequency = imposed_red_f  # if == 0 --> impose frequency f
        self.flag_excitation_type = ex_type  # 'current', 'current_density', 'voltage'
        self.flag_non_linear_core = nonlinear

        phases = np.asarray(phases)
        for num in range(0, self.n_conductors):
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
                    self.delta = np.sqrt(2 / (2 * self.frequency * np.pi * self.sigma * self.mu0))

                    if self.conductor_type[num] == "litz":
                        self.red_freq[num] = self.strand_radius[num] / self.delta
                    else:
                        self.red_freq[num] = self.conductor_radius[num] / self.delta
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
        text_file = open(self.path + "/Parameter.pro", "w")

        # Magnetic Component Type
        if self.component_type == 'inductor':
            text_file.write(f"Flag_Transformer = 0;\n")
        if self.component_type == 'transformer':
            text_file.write(f"Flag_Transformer = 1;\n")

        # Frequency
        text_file.write("Freq = %s;\n" % self.frequency)
        text_file.write(f"delta = {self.delta};\n")

        # Conductor specific definitions
        for num in range(0, self.n_conductors):
            # -- Control Flags --
            if self.flag_excitation_type == 'current':
                text_file.write(f"Flag_ImposedVoltage = 0;\n")
            if self.flag_excitation_type == 'voltage':
                text_file.write(f"Flag_ImposedVoltage = 1;\n")
            if self.conductor_type[num] == 'litz':  # xfmr
                text_file.write(f"Flag_HomogenisedModel{num+1} = 1;\n")
            else:
                text_file.write(f"Flag_HomogenisedModel{num+1} = 0;\n")
            text_file.write("Flag_imposedRr = %s;\n" % self.flag_imposed_reduced_frequency)

            # -- Geometry --
            # Number of conductors
            text_file.write(f"NbrCond{num+1} = {self.turns[num]};\n")
            # For stranded Conductors:
            # text_file.write(f"NbrstrandedCond = {self.turns};\n")  # redundant
            if self.conductor_type[num] == "litz":
                text_file.write(f"NbrStrands{num+1} = {self.n_strands[num]};\n")
                text_file.write(f"Fill{num+1} = {self.FF[num]};\n")
                text_file.write(f"NbrLayers{num+1} = {self.n_layers[num]};\n")
            text_file.write(f"AreaCell{num+1} = {self.A_cell[num]};\n")
            text_file.write(f"Rc{num+1} = {self.conductor_radius[num]};\n")

            # -- Excitation --
            # Imposed current, current density or voltage
            if self.flag_excitation_type == 'current':
                text_file.write(f"Val_EE_{num+1} = {self.current[num]};\n")
                text_file.write(f"Phase_{num+1} = Pi*{self.phase_tmp[num]};\n")
            if self.flag_excitation_type == 'current_density':
                text_file.write(f"Val_EE_{num+1} = {self.current_density[num]};\n")
            if self.flag_excitation_type == 'voltage':
                text_file.write(f"Val_EE_{num+1} = {self.voltage[num]};\n")
            print(f"Cell surface area: {self.A_cell[num]} \n"
                  f"Reduced frequency: {self.red_freq[num]}")
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

        view = 0

        #if self.conductor_type[0] != 'litz':
        if any(type != 'litz' for type in self.conductor_type):
            # Ohmic losses (weightend effective value of current density)
            gmsh.open(self.path + "/res/j2F.pos")
            gmsh.option.setNumber(f"View[{view}].ScaleType", 2)
            gmsh.option.setNumber(f"View[view].RangeType", 2)
            gmsh.option.setNumber(f"View[view].SaturateValues", 1)
            gmsh.option.setNumber(f"View[view].CustomMin", gmsh.option.getNumber(f"View[view].Min") + epsilon)
            gmsh.option.setNumber(f"View[view].CustomMax", gmsh.option.getNumber(f"View[view].Max"))
            gmsh.option.setNumber(f"View[view].ColormapNumber", 1)
            gmsh.option.setNumber(f"View[view].IntervalsType", 2)
            gmsh.option.setNumber(f"View[view].NbIso", 40)
            print(gmsh.option.getNumber("View[view].Max"))
            view += 1

        if any(type == 'litz' for type in self.conductor_type):
            # Ohmic losses (weightend effective value of current density)
            gmsh.open(self.path + "/res/jH.pos")
            gmsh.option.setNumber(f"View[view].ScaleType", 2)
            gmsh.option.setNumber("View[view].RangeType", 2)
            gmsh.option.setNumber("View[view].SaturateValues", 1)
            gmsh.option.setNumber("View[view].CustomMin", gmsh.option.getNumber("View[view].Min") + epsilon)
            gmsh.option.setNumber("View[view].CustomMax", gmsh.option.getNumber("View[view].Max"))
            gmsh.option.setNumber(f"View[view].ColormapNumber", 1)
            gmsh.option.setNumber(f"View[view].IntervalsType", 2)
            gmsh.option.setNumber(f"View[view].NbIso", 40)
            view += 1
        # Magnetic flux density
        gmsh.open(self.path + "/res/Magb.pos")
        gmsh.option.setNumber(f"View[view].ScaleType", 1)
        gmsh.option.setNumber(f"View[view].RangeType", 1)
        gmsh.option.setNumber(f"View[view].CustomMin", gmsh.option.getNumber(f"View[view].Min") + epsilon)
        gmsh.option.setNumber(f"View[view].CustomMax", gmsh.option.getNumber(f"View[view].Max"))
        gmsh.option.setNumber(f"View[view].ColormapNumber", 1)
        gmsh.option.setNumber(f"View[view].IntervalsType", 2)
        gmsh.option.setNumber(f"View[view].NbIso", 40)
        view += 1

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
        #self.ei_axi()

        if self.n_air_gaps != 1:
            raise NotImplementedError

        # == Init ==
        femm.openfemm(non_visualize)
        femm.newdocument(0)
        femm.mi_probdef(freq, 'meters', 'axi', 1.e-8, 0, 30)

        # == Materials ==
        femm.mi_addmaterial('Ferrite', 3000, 3000, 0, 0, 0, 0, 0, 1, 0, 0, 0)
        femm.mi_addmaterial('Air', 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0)
        if self.conductor_type[0] == "litz":
            femm.mi_addmaterial('Copper', 1, 1, 0, 0, sigma, 0, 0, 1, 5, 0, 0, self.n_strands[0], 2*1000*self.strand_radius[0])  # type := 5. last argument
            print(f"Strandsnumber: {self.n_strands[0]}")
            print(f"Strandsdiameter in mm: {2 * 1000 * self.strand_radius[0]}")
        if self.conductor_type[0] == "solid":
            femm.mi_addmaterial('Copper', 1, 1, 0, 0, sigma, 0, 0, 1, 0, 0, 0, 0, 0)



        # == Circuit ==
        # coil as seen from the terminals.
        femm.mi_addcircprop('icoil', current[0], 1)

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
        for num in range(0, self.n_conductors):
            if self.conductor_type[0] == "litz" or self.conductor_type[0] == "solid":
                for i in range(0, int(self.p_conductor[num].shape[0] / 5)):
                    # 0: center | 1: left | 2: top | 3: right | 4.bottom
                    femm.mi_drawarc(self.p_conductor[num][5*i+1][0], self.p_conductor[num][5*i+1][1], self.p_conductor[num][5*i+3][0], self.p_conductor[num][5*i+3][1], 180, 2.5)
                    femm.mi_addarc(self.p_conductor[num][5*i+3][0], self.p_conductor[num][5*i+3][1], self.p_conductor[num][5*i+1][0], self.p_conductor[num][5*i+1][1],  180, 2.5)
                    femm.mi_addblocklabel(self.p_conductor[num][5*i][0], self.p_conductor[num][5*i][1])
                    femm.mi_selectlabel(self.p_conductor[num][5*i][0], self.p_conductor[num][5*i][1])
                    femm.mi_setblockprop('Copper', 0, 1e-4, 'icoil', 0, 0, 1)
                    # femm.mi_setblockprop('Copper', 1, 0, 'icoil', 0, 0, 1)
                    femm.mi_clearselected

        # Define an "open" boundary condition using the built-in function:
        femm.mi_makeABC()

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
        for num in range(0, self.n_conductors):
            if self.conductor_type[num] == 'litz':

                if os.path.isfile(self.path + f"/pre/coeff/pB_RS_la{self.FF[num]}_{self.n_layers[num]}layer.dat"):
                #if os.path.isfile(self.path + f"/pre/coeff/pB_RS_la{self.FF}_0layer.dat"):
                    print("Coefficients for stands approximation are found.")

                else:
                    print("Create coefficients for strands approximation")
                    # need reduced frequency X for pre-simulation
                    # Rounding X to fit it with corresponding parameters from the database
                    X = self.red_freq[num]
                    X = np.around(X, decimals=3)
                    print(f"Rounded Reduced frequency X = {X}")

                    # Create new file with [0 1] for the first element entry

                    # create a new onelab client
                    # -- Pre-Simulation Settings --
                    text_file = open(self.path + "/pre/PreParameter.pro", "w")
                    text_file.write(f"NbrLayers = {self.n_layers[num]};\n")
                    #text_file.write(f"NbrLayers = 0;\n")
                    text_file.write(f"Fill = {self.FF[num]};\n")
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
                            text_file.write(f"NbrLayers = {self.n_layers[num]};\n")
                            #text_file.write(f"NbrLayers = 0;\n")
                            text_file.write(f"Fill = {self.FF[num]};\n")
                            text_file.write(f"Rc = {self.strand_radius[num]};\n")  # double named!!! must be changed
                            text_file.close()

                            # get model file names with correct path
                            input_file = c.getPath('pre/cell_dat.pro')
                            cell = c.getPath('pre/cell.pro')

                            # Run simulations as sub clients
                            mygetdp = self.onelab + 'getdp'
                            c.runSubClient('myGetDP', mygetdp + ' ' + cell + ' -input ' + input_file + ' -solve MagDyn_a -v2')

                    # Formatting stuff
                    files = [self.path + f"/pre/coeff/pB_RS_la{self.FF[num]}_{self.n_layers[num]}layer.dat",
                             self.path + f"/pre/coeff/pI_RS_la{self.FF[num]}_{self.n_layers[num]}layer.dat",
                             self.path + f"/pre/coeff/qB_RS_la{self.FF[num]}_{self.n_layers[num]}layer.dat",
                             self.path + f"/pre/coeff/qI_RS_la{self.FF[num]}_{self.n_layers[num]}layer.dat"]
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
        self.excitation(f=100000, i=1)  # arbitrary values: frequency and current
        self.file_communication()
        self.pre_simulate()

    # ==== Front-End Methods =====
    def single_simulation(self, freq, current, phi=[], skin_mesh_factor=1):
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

    def mesh(self, frequency=None, skin_mesh_factor=1):
        self.high_level_geo_gen(frequency=frequency, skin_mesh_factor=skin_mesh_factor)
        self.generate_mesh()

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

    def litz_loss_comparison(self, FF, n_layers, strand_radius, sim_choice, sweep_parameter='fill_factor', nametag=''):
        """

        :param nametag:
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
        for fill_factors in self.FF:
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

    def get_inductances(self, I0, op_frequency=0, mesh_accuracy=0.5):
        """

        :param I0:
        :param op_frequency:
        :return:
        """

        # Remove "old" Inductance Logs
        try:
            os.remove(self.path + "/res/L_11.dat")
            os.remove(self.path + "/res/L_22.dat")
        except:
            print("Could not find Inductance logs")

        # -- Inductance Estimation --
        self.mesh(frequency=op_frequency, skin_mesh_factor=mesh_accuracy)
        # 2nd-open
        frequencies = [op_frequency] * 4
        currents = [[1, 0], [0, 1], [I0, self.turns[0] / self.turns[1] * I0], [self.turns[1] / self.turns[0] * I0, I0]]
        self.excitation_sweep(frequencies=frequencies, currents=currents, show_last=1)

        # Open loss/error results
        inductance_matrix = open(self.path + "/res/Inductance_Matrix.dat", "w")
        # Read the logged inductance values
        with open(self.path + "/res/L_11.dat") as f:
            next_line = False
            count_lines = 0
            for line in f:
                if next_line:
                    next_line = False
                    count_lines += 1
                    words = line.split(sep=' ')
                    if count_lines == 1:
                        L_11 = float(words[2])
                    if count_lines == 2:
                        L_k1 = float(words[2])
                if "L_11" in line:
                    next_line = True
                if count_lines == 2:
                    break
        # Read the logged inductance values
        with open(self.path + "/res/L_22.dat") as f:
            next_line = False
            count_lines = 0
            for line in f:
                if next_line:
                    next_line = False
                    words = line.split(sep=' ')
                    L_22 = float(words[2])
                    count_lines += 1
                if "L_22" in line:
                    next_line = True
                if count_lines == 1:
                    break

        # Calculation of inductance matrix
        L_s1 = 0.5*(L_11-self.turns[0]**2/self.turns[1]**2*L_22+L_k1)
        L_m1 = L_11 - L_s1
        L_m = self.turns[1]/self.turns[0]*L_m1
        L_m2 = self.turns[1]/self.turns[0]*L_m
        L_s2 = L_22 - L_m2

        print(f"L_11 = {L_11}\n"
              f"L_22 = {L_22}\n"
              f"L_m = {L_m}\n"
              f"L_s1 = {L_s1}\n"
              f"L_s2 = {L_s2}\n")
        """
        if read == 1 and ' ' in line:
            # words = line.split(sep=' ')
            words = line.replace(' ', ', ')
            for word in words:
                inductance_matrix.write(word)
        """
        inductance_matrix.close()



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
        if output.shape[0] == 0:
            print("Only one air gap in this leg. No island needed.")
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

    path = input(f"Please enter the parent folder path of {destination} in ways of 'C:.../onelab-Windows64/': ")
    dict = {"onelab": path}
    file = open(config_file, 'w', encoding='utf-8')
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