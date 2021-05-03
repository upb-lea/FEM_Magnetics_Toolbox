import gmsh
import numpy as np
import os
import pathlib
import sys
from onelab import onelab
from functions import inner_points, min_max_inner_points, call_for_path, NbrStrands


class MagneticComponent:
    def __init__(self):

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
        self.conductor_type = "litz"  # Stranded wires
        """
        conductor_type = "stacked"  # Vertical packing of conductors
        conductor_type = "full"  # One massive Conductor in each window
        conductor_type = "foil"  # Horizontal packing of conductors
        conductor_type = "solid"  # Massive wires
        """

        if self.conductor_type == 'solid':
            self.conductor_radius = 0.0011879
        # Litz Approximation
        self.FF = 0.9  # hexagonal packing: ~90.7% are theoretical maximum
        self.n_layers = 6
        self.n_strands = NbrStrands(self.n_layers)
        self.strand_radius = 0.1e-3
        if self.conductor_type == 'litz':
            self.conductor_radius = np.sqrt(self.n_strands/self.FF)*self.strand_radius  # Must be calculated from strand
        self.A_cell = np.pi * self.conductor_radius**2  #* FF  # Surface of the litz approximated hexagonal cell

        # -- Materials --
        # frequency = 0: mu_rel only used if flag_non_linear_core == 0
        # frequency > 0: mu_rel is used
        self.mu0 = 4e-7*np.pi
        self.mu_rel = 3000   # relative Core Permeability
        self.core_material = 95  # 95 := TDK-N95 | Currently only works with Numbers corresponding to BH.pro
        self.sigma = 6e7

        # -- Characteristic lengths -- [for mesh sizes]
        self.c_core = self.core_w/10. * self.s
        self.c_window = self.window_w/10 * self.s
        self.c_conductor = self.window_w/10 * self.s

        # -- Further Attributes --
        self.path = None
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
        self.flag_imposed_reduced_frequency = None
        self.red_freq = None
        self.delta = None
        self.red_freq = None
        self.delta = None
        self.red_freq = None

    def high_level_geo_gen(self):

        # ==============================
        # High-Level Geometry Generation
        # ==============================

        # -- Core-type --
        if self.core_type == "EI":
            self.n_windows = 2

        # -- Symmetry -- [Choose between asymmetric, symmetric and axi symmetric]
        if self.y_symmetric == 0:
            self.axi_symmetric = 0
        if self.axi_symmetric == 1:
            self.y_symmetric = 1


        # -- Air Gap Data -- [random air gap generation]
        if self.y_symmetric == 0:
            raise(NotImplementedError, "Not up to date! Air Gap generation must be adopted from symmetric case")
            """
            air_gaps = np.empty((n_air_gaps, 4))
            for i in range(0, n_air_gaps):
                if i < n_air_gaps/2:
                    position_tag = -1
                else:
                    position_tag = 1
                airgap_position = np.random.rand(1)*window_h-window_h/2
                airgap_h = np.random.rand(1)*0.005
                c_airgap = airgap_h / 3 * s
                air_gaps[i, :] = np.array([position_tag, airgap_position, airgap_h, c_airgap])
                """
        if self.y_symmetric == 1:
            air_gaps = np.empty((self.n_air_gaps, 4))
            i = 0
            while i in range(0, self.n_air_gaps):
                position_tag = 0  # '-1': left leg | '0': center leg | '1': right leg
                if self.n_air_gaps == 1:
                    airgap_h = 0.0005
                    airgap_position = 0.5*(self.window_h-airgap_h)-(self.window_h/2-airgap_h/2)
                else:
                    np.random.rand(1)*0.005 + 0.001
                    np.random.rand(1)*(self.window_h-airgap_h)-(self.window_h/2-airgap_h/2)
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
        if self.core_type == "EI":

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

            if self.y_symmetric == 1:
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
                    # left window
                    # p_conductor[4*i+0][:] = [-r_inner + core_cond_isolation, -window_h/2 + core_cond_isolation, 0, c_conductor]
                    # p_conductor[4*i+1][:] = [-core_cond_isolation - core_w/2, -window_h/2 + core_cond_isolation, 0, c_conductor]
                    # p_conductor[4*i+2][:] = [-r_inner + core_cond_isolation, window_h/2 - core_cond_isolation, 0, c_conductor]
                    # p_conductor[4*i+3][:] = [-core_cond_isolation - core_w/2, window_h/2 - core_cond_isolation, 0, c_conductor]
                    # right window
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

                if self.conductor_type == "litz" or conductor_type == "solid":
                    self.p_conductor = []  # center points are stored
                    left_bound = self.core_cond_isolation + self.core_w/2
                    right_bound = r_inner - self.core_cond_isolation
                    top_bound = self.window_h/2
                    bot_bound = -self.window_h/2

                    y = bot_bound + self.core_cond_isolation + self.conductor_radius
                    x = left_bound + self.core_cond_isolation + self.conductor_radius
                    i = 0
                    # Case n_conductors higher that "allowed" is missing
                    while y < top_bound:
                        while x < right_bound and i < self.n_conductors:
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
                        print("Could not resolve all conductors")

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

        # Parent folder path
        self.path = str(pathlib.Path(__file__).parent.absolute())

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

    def excitation(self, f, I):

        # -- Excitation --
        self.flag_imposed_reduced_frequency = 0 # if == 0 --> impose frequency f
        self.flag_excitation_type = 'current'  # 'current', 'current_density', 'voltage'
        self.flag_non_linear_core = 0

        # Imposed current, current density or voltage
        if self.flag_excitation_type == 'current':
            self.current = I
        if self.flag_excitation_type == 'current_density':
            raise NotImplementedError
        if self.flag_excitation_type == 'voltage':
            self.voltage = 2

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
        text_file.write(f"AreaCell = {self.A_cell};\n")
        # Coordinates of the rectangular winding window
        if self.axi_symmetric == 1:
            text_file.write("Xw1 = %s;\n" % self.p_window[4, 0])
            text_file.write("Xw2 = %s;\n" % self.p_window[5, 0])
        else:
            raise NotImplementedError("Only axi-symmetric case implemented :(")

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
        # -------------------------------------- Onelab Setup -------------------------------------
        if os.path.isfile(self.path + "/config.py"):
            print("if")
            import config
            with open('config.py') as f:
                if 'mygetdp' in f.read():
                    if os.path.isfile(config.mygetdp + ".exe") or os.path.isfile(config.mygetdp):
                        mygetdp = config.mygetdp
                else:
                    mygetdp = call_for_path("mygetdp")
        else:
            print("else")
            mygetdp = call_for_path("mygetdp")

        # -------------------------------------- Pre-Simulation --------------------------------------
        if self.conductor_type == 'litz':
            if os.path.isfile(self.path + f"/pre/coeff/pB_RS_la{self.FF}_1layer.dat"):
                print("Coefficients for stands approximation are found.")

            else:
                print("Create coefficients for strands approximation")
                # need reduced frequency X for pre-simulation
                X = self.red_freq
                print(f"Exact Reduced frequency X = {X}")

                # Rounding X to fit it with corresbonding parameters from the database
                X = np.around(X, decimals=3)
                print(f"Rounded Reduced frequency X = {X}")

                """
                # create a new onelab client
                c = onelab.client(__file__)


                # get model file names with correct path
                input_file = c.getPath('cell_dat.pro')
                solver = c.getPath('ind_axi_python_controlled' + '.pro')

                # Run simulations as sub clients (non blocking??)
                c.runSubClient('myGetDP', mygetdp + ' ' + solver + ' -msh ' + msh_file + ' -solve Analysis -v2')
                """

        # ---------------------------------------- Simulation ---------------------------------------
        # create a new onelab client
        c = onelab.client(__file__)

        # get model file names with correct path
        msh_file = c.getPath('geometry.msh')
        solver = c.getPath('ind_axi_python_controlled' + '.pro')

        # Run simulations as sub clients (non blocking??)
        c.runSubClient('myGetDP', mygetdp + ' ' + solver + ' -msh ' + msh_file + ' -solve Analysis -v2')

    def visualize(self):
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

    def single_simulation(self):
        self.high_level_geo_gen()
        self.generate_mesh()
        self.excitation(f=100000, I=3)  # frequency and current
        self.file_communication()
        self.simulate()
        self.visualize()

    def freq_sweep_simulation(self, start, end, steps):
        self.high_level_geo_gen()
        self.generate_mesh()
        sweep_list = np.linspace(start, end, steps)
        for i in range(0, len(sweep_list)):
            self.excitation(f=sweep_list[i], I=3)  # frequency and current
            self.file_communication()
            self.simulate()
        self.visualize()



"""Theoretical inductivity
lz = 0.5
A = lz*(core_w-window_w)/2
mu0 = 1.256*10**(-6)
N = 10
L_theo = None
if n_air_gaps > 0:
    L_theo = N**2 * mu0 * A / airgap_h * 1000
print("Nicht mehr aktuell! Theoretical, idealized inductance: ", L_theo, "mH")
"""