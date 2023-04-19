# Python standard libraries
import os
from re import I
import numpy as np
from typing import Dict, List

# Third parry libraries
import gmsh

# Local libraries
import femmt.functions as ff
from femmt.enumerations import ComponentType, ConductorType, WindingType, CoreType
from femmt.data import FileData
from femmt.model import Conductor, Core, StrayPath, AirGaps, Insulation
from femmt.drawing import TwoDaxiSymmetric


class Mesh:
    """
    This class will create a mesh from the given model.
    """
    model: TwoDaxiSymmetric
    core: Core
    stray_path: StrayPath
    insulation: Insulation
    component_type: ComponentType
    windings: int
    air_gaps: List[AirGaps]
    correct_outer_leg: bool
    region: bool

    # File paths
    hybrid_color_png_file: str
    model_geo_fil: str
    mesh_folder_path: str
    e_m_mesh_file: str
    thermal_mesh_file: str

    # Additionaly there are all the needed lists for points, lines, curve_loops and plane_surfaces
    # See set_empty_lists()

    def __init__(self, model: TwoDaxiSymmetric, windings: List[Conductor], correct_outer_leg: bool, file_paths: FileData, region: bool = None, silent: bool = False):

        # Initialize gmsh once
        if not gmsh.isInitialized():
            gmsh.initialize()

        if silent:
            gmsh.option.set_number("General.Verbosity", 1)
        else:
            gmsh.option.set_number("General.Verbosity", 5)

        self.model = model
        self.core = model.core
        self.stray_path = model.stray_path
        self.insulation = model.insulation
        self.component_type = model.component_type
        self.windings = windings
        self.air_gaps = model.air_gaps
        self.correct_outer_leg = correct_outer_leg
        self.region = region  # Apply an outer Region or directly apply a constraint on the Core Boundary
        self.mesh_data = model.mesh_data

        # Files
        self.hybrid_color_png_file = file_paths.hybrid_color_visualize_file
        self.model_geo_file = file_paths.model_geo_file
        self.mesh_folder_path = file_paths.mesh_folder_path
        self.e_m_mesh_file = file_paths.e_m_mesh_file
        self.thermal_mesh_file = file_paths.thermal_mesh_file

    def set_empty_point_lists(self):
        # Add empty lists
        p_core = []
        p_island = []
        p_cond = []
        for num in range(len(self.windings)):
            p_cond.append([])
        p_region = []
        p_iso_core = []
        return p_core, p_island, p_cond, p_region, p_iso_core

    def set_empty_line_lists(self):

        # Curves
        l_bound_core = []
        l_bound_air = []
        l_core_air = []
        l_cond = []
        for num in range(len(self.windings)):
            l_cond.append([])
        l_region = []
        l_air_gaps_air = []
        l_iso_core = []
        return l_bound_core, l_bound_air, l_core_air, l_cond, l_region, l_air_gaps_air, l_iso_core

    def set_empty_curve_loop_lists(self):
        # Curve Loops
        curve_loop_cond = []
        for num in range(len(self.windings)):
            curve_loop_cond.append([])
        curve_loop_island = []
        curve_loop_air = []
        curve_loop_air_gaps = []
        curve_loop_iso_core = []
        # curve_loop_outer_air = []
        # curve_loop_bound = []
        return curve_loop_cond, curve_loop_island, curve_loop_air, curve_loop_air_gaps, curve_loop_iso_core

    def set_empty_plane_lists(self):
        # This is only needed for the surfaces since they are the only global variables
        # Plane Surfaces
        self.plane_surface_core = []


        self.plane_surface_cond = []
        for num in range(len(self.windings)):
            self.plane_surface_cond.append([])
        if self.core.core_type == CoreType.Single:
            self.plane_surface_air = []
        elif self.core.core_type == CoreType.Stacked:
            self.plane_surface_air_top = []
            self.plane_surface_air_bot = []
        self.plane_surface_outer_air = []
        self.plane_surface_air_gaps = []
        self.plane_surface_iso_core = []
        self.plane_surface_iso_pri_sec = []

    def single_core(self,
                    p_core: list, p_island: list,
                    l_bound_core: list, l_core_air: list, l_bound_air: list, l_air_gaps_air: list,
                    curve_loop_island: list, curve_loop_air_gaps: list):
        """sets a single core with one physical winding window
        """
        # Points
        # (index refers to sketch)

        # First point (left point of lowest air gap)
        if self.model.air_gaps.number > 0:
            p_core.append(gmsh.model.geo.addPoint(0,
                                                  self.model.p_air_gaps[0][1],
                                                  0,
                                                  self.model.p_air_gaps[0][3]))
        elif self.model.air_gaps.number == 0:
            p_core.append(None)  # dummy filled for no air gap special case
        else:
            raise Exception("Negative air gaps number?")

        # Go down and counter-clockwise
        # Four points around the core
        if self.correct_outer_leg:
            correction_of_outer_points = 0.00
        else:
            correction_of_outer_points = 0

        p_core.append(gmsh.model.geo.addPoint(0,
                                              self.model.p_outer[1][1],
                                              self.model.p_outer[1][2],
                                              self.model.p_outer[1][3]))

        p_core.append(gmsh.model.geo.addPoint(self.model.p_outer[1][0],
                                              self.model.p_outer[1][1] + correction_of_outer_points,
                                              self.model.p_outer[1][2],
                                              self.model.p_outer[1][3]))

        p_core.append(gmsh.model.geo.addPoint(self.model.p_outer[3][0],
                                              self.model.p_outer[3][1] - correction_of_outer_points,
                                              self.model.p_outer[3][2],
                                              self.model.p_outer[3][3]))

        p_core.append(gmsh.model.geo.addPoint(0,
                                              self.model.p_outer[3][1],
                                              self.model.p_outer[3][2],
                                              self.model.p_outer[3][3]))

        # Two points of highest air gap
        if self.model.air_gaps.number > 0:
            p_core.append(gmsh.model.geo.addPoint(0,
                                                  self.model.p_air_gaps[-2][1],
                                                  self.model.p_air_gaps[-2][2],
                                                  self.model.p_air_gaps[-2][3]))

            p_core.append(gmsh.model.geo.addPoint(self.model.p_air_gaps[-1][0],
                                                  self.model.p_air_gaps[-1][1],
                                                  self.model.p_air_gaps[-1][2],
                                                  self.model.p_air_gaps[-1][3]))
        else:
            p_core.append(None)  # dummy filled for no air gap special case
            p_core.append(None)  # dummy filled for no air gap special case

        # Clockwise
        # Four points of the window
        p_core.append(gmsh.model.geo.addPoint(self.model.p_window[6][0],
                                              self.model.p_window[6][1],
                                              self.model.p_window[6][2],
                                              self.model.p_window[6][3]))

        p_core.append(gmsh.model.geo.addPoint(self.model.p_window[7][0],
                                              self.model.p_window[7][1],
                                              self.model.p_window[7][2],
                                              self.model.p_window[7][3]))

        p_core.append(gmsh.model.geo.addPoint(self.model.p_window[5][0],
                                              self.model.p_window[5][1],
                                              self.model.p_window[5][2],
                                              self.model.p_window[5][3]))

        p_core.append(gmsh.model.geo.addPoint(self.model.p_window[4][0],
                                              self.model.p_window[4][1],
                                              self.model.p_window[4][2],
                                              self.model.p_window[4][3]))

        # Last point of lowest air gap
        if self.model.air_gaps.number > 0:
            p_core.append(gmsh.model.geo.addPoint(self.model.p_air_gaps[1][0],
                                                  self.model.p_air_gaps[1][1],
                                                  self.model.p_air_gaps[1][2],
                                                  self.mesh_data.c_window))
        else:
            p_core.append(None)  # dummy filled for no air gap special case

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Curves
        # (index refers to sketch)
        # To be added: Case Air Gaps directly on outer leg

        # Curves: Boundary - Core
        if self.model.air_gaps.number > 0:
            l_bound_core.append(gmsh.model.geo.addLine(p_core[0],
                                                       p_core[1]))
            l_bound_core.append(gmsh.model.geo.addLine(p_core[4],
                                                       p_core[5]))
        else:
            l_bound_core.append(gmsh.model.geo.addLine(p_core[4],
                                                       p_core[1]))
        l_bound_core.append(gmsh.model.geo.addLine(p_core[1],
                                                   p_core[2]))
        l_bound_core.append(gmsh.model.geo.addLine(p_core[2],
                                                   p_core[3]))
        l_bound_core.append(gmsh.model.geo.addLine(p_core[3],
                                                   p_core[4]))
        # Curves: Core - Air
        if self.model.air_gaps.number > 0:
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

        if self.model.air_gaps.number > 0:
            l_core_air.append(gmsh.model.geo.addLine(p_core[10],
                                                     p_core[11]))
            l_core_air.append(gmsh.model.geo.addLine(p_core[11],
                                                     p_core[0]))
        else:
            l_core_air.append(gmsh.model.geo.addLine(p_core[10],
                                                     p_core[7]))

        # Plane: Main Core --> plane_surface_core[0]
        if self.model.air_gaps.number > 0:
            curve_loop_core = gmsh.model.geo.addCurveLoop(l_bound_core + l_core_air)
            self.plane_surface_core.append(gmsh.model.geo.addPlaneSurface([-curve_loop_core]))
        else:
            curve_loop_bound_core = gmsh.model.geo.addCurveLoop(l_bound_core)
            curve_loop_core_air = gmsh.model.geo.addCurveLoop(l_core_air)
            self.plane_surface_core.append(gmsh.model.geo.addPlaneSurface([-curve_loop_bound_core, curve_loop_core_air]))

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # - Core parts between Air Gaps
        # Points of Core Islands (index refers to sketch)

        stray_path_air_gap_top_point = None
        stray_path_air_gap_bot_point = None
        l_core_air_air_gap = []

        if self.model.air_gaps.number > 1:
            for point in self.model.p_air_gaps[2:-2]:
                p_island.append(gmsh.model.geo.addPoint(*point))

            # Add two more points for closing of the air gap for a stray_path
            if self.component_type == ComponentType.IntegratedTransformer:
                stray_path_air_gap_top_point = gmsh.model.geo.addPoint(*self.model.p_close_air_gaps[0])
                stray_path_air_gap_bot_point = gmsh.model.geo.addPoint(*self.model.p_close_air_gaps[1])

            # Curves of Core Islands (index refers to sketch)
            for i in range(0, int(len(p_island) / 4)):
                if self.component_type == ComponentType.IntegratedTransformer and self.stray_path.start_index == i:
                    l_core_air_air_gap.append(gmsh.model.geo.addLine(p_island[4 * i + 0], stray_path_air_gap_bot_point))
                    l_core_air.append(gmsh.model.geo.addLine(stray_path_air_gap_bot_point, p_island[4 * i + 1]))
                    l_core_air.append(gmsh.model.geo.addLine(p_island[4 * i + 1], p_island[4 * i + 3]))
                    l_core_air.append(gmsh.model.geo.addLine(p_island[4 * i + 3], stray_path_air_gap_top_point))
                    l_core_air_air_gap.append(gmsh.model.geo.addLine(stray_path_air_gap_top_point, p_island[4 * i + 2]))
                    l_bound_core.append(gmsh.model.geo.addLine(p_island[4 * i + 2], p_island[4 * i + 0]))

                    curve_loop_island.append(gmsh.model.geo.addCurveLoop(
                        [l_core_air_air_gap[-2], l_core_air[-3], l_core_air[-2], l_core_air_air_gap[-1], l_core_air[-1], l_bound_core[-1]]))
                else:
                    # Default
                    l_core_air.append(gmsh.model.geo.addLine(p_island[4 * i + 0], p_island[4 * i + 1]))
                    l_core_air.append(gmsh.model.geo.addLine(p_island[4 * i + 1], p_island[4 * i + 3]))
                    l_core_air.append(gmsh.model.geo.addLine(p_island[4 * i + 3], p_island[4 * i + 2]))
                    l_bound_core.append(gmsh.model.geo.addLine(p_island[4 * i + 2], p_island[4 * i + 0]))

                    curve_loop_island.append(gmsh.model.geo.addCurveLoop(
                        [l_core_air[-3], l_core_air[-2], l_core_air[-1], l_bound_core[-1]]))

                # Iterative plane creation
                self.plane_surface_core.append(gmsh.model.geo.addPlaneSurface([-curve_loop_island[-1]]))

        # Curves: Boundary - Air
        if self.model.air_gaps.number == 1:
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

        # Curves: Close air gaps
        if self.model.air_gaps.number > 0:
            for i in range(self.model.air_gaps.number):
                bottom_point = p_core[11] if i == 0 else p_island[(i - 1) * 4 + 3]
                top_point = p_core[6] if i == self.model.air_gaps.number - 1 else p_island[i * 4 + 1]

                if self.component_type == ComponentType.IntegratedTransformer:
                    if self.stray_path.start_index == i:
                        # Stray path is above current air_gap
                        top_point = stray_path_air_gap_bot_point
                    elif self.stray_path.start_index + 1 == i:
                        # Stray path is below current air_gap
                        bottom_point = stray_path_air_gap_top_point

                l_air_gaps_air.append(gmsh.model.geo.addLine(bottom_point, top_point))

        for i in range(self.model.air_gaps.number):
            left = l_bound_air[i]
            right = l_air_gaps_air[i]

            if i == self.model.air_gaps.number - 1:
                top = l_core_air[0]
            elif self.component_type == ComponentType.IntegratedTransformer and self.stray_path.start_index == i:
                top = l_core_air_air_gap[0]
            else:
                top = l_core_air[7 + 3 * i]

            if i == 0:
                bottom = l_core_air[6]
            elif self.component_type == ComponentType.IntegratedTransformer and self.stray_path.start_index + 1 == i:
                bottom = l_core_air_air_gap[1]
            else:
                bottom = l_core_air[6 + 3 * i]

            curve_loop = gmsh.model.geo.addCurveLoop([left, top, bottom, right], -1, True)
            curve_loop_air_gaps.append(curve_loop)
            self.plane_surface_air_gaps.append(gmsh.model.geo.addPlaneSurface([curve_loop]))

    def stacked_core(self, p_core: list, l_bound_core: list, l_core_air: list, l_bound_air: list):
        # Points
        # (index refers to sketch)
        from operator import itemgetter
        index_bot, bot_air_gap_center_y = min(enumerate(self.model.air_gaps.midpoints[:][1]), key=itemgetter(1))
        bot_air_gap_length = self.model.air_gaps.midpoints[index_bot][2]

        index_top, top_air_gap_center_y = max(enumerate(self.model.air_gaps.midpoints[:][1]), key=itemgetter(1))
        top_air_gap_length = self.model.air_gaps.midpoints[index_top][2]

        # 0
        p_core.append(gmsh.model.geo.addPoint(0,
                                              bot_air_gap_center_y - bot_air_gap_length / 2,
                                              0,
                                              self.model.p_outer[1][3]))  # TODO accuracy of mesh
        # 1
        p_core.append(gmsh.model.geo.addPoint(0,
                                              self.model.p_outer[1][1],
                                              self.model.p_outer[1][2],
                                              self.model.p_outer[1][3]))  # TODO accuracy of mesh
        # 2 bot right
        p_core.append(gmsh.model.geo.addPoint(self.model.p_outer[1][0],
                                              self.model.p_outer[1][1],
                                              self.model.p_outer[1][2],
                                              self.model.p_outer[1][3]))  # TODO accuracy of mesh
        # 3 top right
        p_core.append(gmsh.model.geo.addPoint(self.model.p_outer[3][0],
                                              self.model.p_outer[3][1] + self.core.window_h_top + self.core.core_inner_diameter / 4,
                                              self.model.p_outer[3][2],
                                              self.model.p_outer[1][3]))  # TODO accuracy of mesh
        # 4 top center
        p_core.append(gmsh.model.geo.addPoint(0,
                                              self.model.p_outer[3][1] + self.core.window_h_top + self.core.core_inner_diameter / 4,
                                              self.model.p_outer[3][2],
                                              self.model.p_outer[1][3]))  # TODO accuracy of mesh

        # 5
        p_core.append(gmsh.model.geo.addPoint(0,
                                              top_air_gap_center_y + top_air_gap_length / 2,
                                              0,
                                              self.model.p_outer[1][3]))  # TODO accuracy of mesh

        # 6
        p_core.append(gmsh.model.geo.addPoint(self.model.p_window_top[0][0],
                                              top_air_gap_center_y + top_air_gap_length / 2,
                                              0,
                                              self.model.p_window_bot[0][3]))  # TODO accuracy of mesh

        # 7
        p_core.append(gmsh.model.geo.addPoint(self.model.p_window_top[2][0],
                                              self.model.p_window_top[2][1],
                                              self.model.p_window_top[2][2],
                                              self.model.p_window_top[2][3]))

        # 8
        p_core.append(gmsh.model.geo.addPoint(self.model.p_window_top[3][0],
                                              self.model.p_window_top[3][1],
                                              self.model.p_window_top[3][2],
                                              self.model.p_window_top[3][3]))

        # 9
        p_core.append(gmsh.model.geo.addPoint(self.model.p_window_top[1][0],
                                              self.model.p_window_top[1][1],
                                              self.model.p_window_top[1][2],
                                              self.model.p_window_top[1][3]))

        # 10
        p_core.append(gmsh.model.geo.addPoint(0,
                                              self.model.p_window_top[1][1],
                                              self.model.p_window_top[1][2],
                                              self.model.p_outer[1][3]))  # TODO air gap accuracy

        # 11 bot air gap # TODO: Code for distributed air gaps
        p_core.append(gmsh.model.geo.addPoint(0,
                                              bot_air_gap_center_y + bot_air_gap_length / 2,
                                              0,
                                              self.model.p_outer[1][3]))  # TODO: Code for distributed air gaps

        # 12 bot air gap # TODO: Code for distributed air gaps
        p_core.append(gmsh.model.geo.addPoint(self.model.p_window_bot[2][0],
                                              bot_air_gap_center_y + bot_air_gap_length / 2,
                                              0,
                                              self.model.p_window_bot[0][3]))  # TODO: Code for distributed air gaps

        # 13
        p_core.append(gmsh.model.geo.addPoint(self.model.p_window_bot[2][0],
                                              self.model.p_window_bot[2][1],
                                              self.model.p_window_bot[2][2],
                                              self.model.p_window_bot[2][3]))

        # 14
        p_core.append(gmsh.model.geo.addPoint(self.model.p_window_bot[3][0],
                                              self.model.p_window_bot[3][1],
                                              self.model.p_window_bot[3][2],
                                              self.model.p_window_bot[3][3]))

        # 15
        p_core.append(gmsh.model.geo.addPoint(self.model.p_window_bot[1][0],
                                              self.model.p_window_bot[1][1],
                                              self.model.p_window_bot[1][2],
                                              self.model.p_window_bot[1][3]))

        # 16
        p_core.append(gmsh.model.geo.addPoint(self.model.p_window_bot[0][0],
                                              self.model.p_window_bot[0][1],
                                              self.model.p_window_bot[0][2],
                                              self.model.p_window_bot[0][3]))

        # 16
        p_core.append(gmsh.model.geo.addPoint(self.model.p_window_bot[0][0],
                                              bot_air_gap_center_y - bot_air_gap_length / 2,
                                              0,
                                              self.model.p_window_bot[0][3]))  # TODO:air gap accuracy

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Curves
        # (index refers to sketch)

        # Curves: Boundary - Core
        l_bound_core.append(gmsh.model.geo.addLine(p_core[0],
                                                   p_core[1]))

        l_bound_core.append(gmsh.model.geo.addLine(p_core[1],
                                                   p_core[2]))

        l_bound_core.append(gmsh.model.geo.addLine(p_core[2],
                                                   p_core[3]))

        l_bound_core.append(gmsh.model.geo.addLine(p_core[3],
                                                   p_core[4]))

        l_bound_core.append(gmsh.model.geo.addLine(p_core[4],
                                                   p_core[5]))

        l_bound_core.append(gmsh.model.geo.addLine(p_core[10],
                                                   p_core[11]))

        # Curves: Bound - Air
        l_bound_air.append(gmsh.model.geo.addLine(p_core[5],
                                                  p_core[10]))

        l_bound_air.append(gmsh.model.geo.addLine(p_core[11],
                                                  p_core[0]))

        # Curves: Core - Air
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

        l_core_air.append(gmsh.model.geo.addLine(p_core[11],
                                                 p_core[12]))

        l_core_air.append(gmsh.model.geo.addLine(p_core[12],
                                                 p_core[13]))

        l_core_air.append(gmsh.model.geo.addLine(p_core[13],
                                                 p_core[14]))

        l_core_air.append(gmsh.model.geo.addLine(p_core[14],
                                                 p_core[15]))

        l_core_air.append(gmsh.model.geo.addLine(p_core[15],
                                                 p_core[16]))

        l_core_air.append(gmsh.model.geo.addLine(p_core[16],
                                                 p_core[17]))

        l_core_air.append(gmsh.model.geo.addLine(p_core[17],
                                                 p_core[0]))

        # Plane: Main Core --> plane_surface_core[0]
        curve_loop_core = gmsh.model.geo.addCurveLoop(l_bound_core + l_core_air)  # TODO: must be taken into account that its a kind of interrupted lines
        self.plane_surface_core.append(gmsh.model.geo.addPlaneSurface([-curve_loop_core]))

    def conductors(self, p_cond: list, l_cond: list, curve_loop_cond: list):
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Conductors
        # Points of Conductors
        for num in range(len(self.windings)):
            for i in range(self.model.p_conductor[num].shape[0]):
                p_cond[num].append(
                    gmsh.model.geo.addPoint(
                        self.model.p_conductor[num][i][0],
                        self.model.p_conductor[num][i][1],
                        0,
                        self.model.p_conductor[num][i][3]))

            # Curves of Conductors
            if self.windings[num].conductor_type in [ConductorType.RoundLitz, ConductorType.RoundSolid]:
                # Round conductor
                for i in range(int(len(p_cond[num]) / 5)):
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
                    curve_loop_cond[num].append(gmsh.model.geo.addCurveLoop([
                        l_cond[num][i * 4 + 0],
                        l_cond[num][i * 4 + 1],
                        l_cond[num][i * 4 + 2],
                        l_cond[num][i * 4 + 3]]))
                    self.plane_surface_cond[num].append(
                        gmsh.model.geo.addPlaneSurface([curve_loop_cond[num][i]]))
            else:
                # Rectangle conductor cut
                for i in range(int(len(p_cond[num]) / 4)):
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
                    self.plane_surface_cond[num].append(
                        gmsh.model.geo.addPlaneSurface([curve_loop_cond[num][i]]))

    def insulations_core_cond(self, p_iso_core: list):
        """
        Sets the rectangular electrical insulation between conductors and core.
        :return:
        """
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Insulations
        # Core to Windings
        if self.model.p_iso_core:  # Check if list is not empty
            # Points
            for iso in self.model.p_iso_core:
                p_iso = []
                for i in iso:
                    p_iso.append(gmsh.model.geo.addPoint(i[0], i[1], i[2], i[3]))
                p_iso_core.append(p_iso)
            # Lines
            l_iso_core = [[gmsh.model.geo.addLine(iso[i], iso[(i + 1) % 4]) for i in range(4)] for iso in p_iso_core]

            # Curve loop and surface
            curve_loop_iso_core = []
            self.plane_surface_iso_core = []
            for iso in l_iso_core:
                cl = gmsh.model.geo.addCurveLoop(iso)
                curve_loop_iso_core.append(cl)
                self.plane_surface_iso_core.append(gmsh.model.geo.addPlaneSurface([cl]))

    def air_single(self, l_core_air: list, l_air_gaps_air: list, curve_loop_air: list, curve_loop_cond: list, curve_loop_iso_core: list):
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
        if self.model.air_gaps.number == 0:
            l_air_tmp = l_core_air
        else:
            l_air_tmp = l_core_air[1:6] + l_air_gaps_air

            for i in range(self.model.air_gaps.number - 1):
                if self.component_type == ComponentType.IntegratedTransformer and i == self.stray_path.start_index:
                    l_air_tmp.append(l_core_air[7 + 3 * i])
                    l_air_tmp.append(l_core_air[8 + 3 * i])
                    l_air_tmp.append(l_core_air[9 + 3 * i])
                else:
                    l_air_tmp.append(l_core_air[8 + 3 * i])

        # for i in range(0, self.component.air_gaps.number):
        #    l_air_tmp.append(self.l_air_gaps_air[i])
        #    l_air_tmp.append(self.l_air_gaps_air[i+1])

        # self.curve_loop_air.append(gmsh.model.geo.addCurveLoop(l_air_tmp, -1, True))
        # for i in range(0, self.component.air_gaps.number):
        #    l_air_tmp.append(self.l_air_gaps_air[i])
        #    l_air_tmp.append(self.l_air_gaps_air[i+1])

        curve_loop_air.append(gmsh.model.geo.addCurveLoop(l_air_tmp, -1, True))

        # Need flatten list of all! conductors
        flatten_curve_loop_cond = [j for sub in curve_loop_cond for j in sub]

        # The first curve loop represents the outer bounds: self.curve_loop_air (should only contain one element)
        # The other curve loops represent holes in the surface -> For each conductor as well as each insulation
        self.plane_surface_air.append(gmsh.model.geo.addPlaneSurface(curve_loop_air + flatten_curve_loop_cond + curve_loop_iso_core))

    def air_stacked(self, l_core_air: list, l_bound_air, curve_loop_cond: list):
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
        # TODO: Split conductors in top and bottom winding window
        primary_in_top, secondary_in_top, tertiary_in_top, primary_in_bot, secondary_in_bot, tertiary_in_bot = [], [], [], [], [], []
        n_primary_in_top, n_secondary_in_top, n_tertiary_in_top, n_primary_in_bot, n_secondary_in_bot, n_tertiary_in_bot = 0, 0, 0, 0, 0, 0

        for i, ww in enumerate(self.model.winding_windows):
            for vww in ww.virtual_winding_windows:
                if i == 0:  # top window
                    try:
                        n_primary_in_top += vww.turns[0]
                    except IndexError:
                        pass
                    try:
                        n_secondary_in_top += vww.turns[1]
                    except IndexError:
                        pass
                    try:
                        n_tertiary_in_top += vww.turns[2]
                    except IndexError:
                        pass
                elif i == 1:  # bot window
                    try:
                        n_primary_in_bot += vww.turns[0]
                    except IndexError:
                        pass
                    try:
                        n_secondary_in_bot += vww.turns[1]
                    except IndexError:
                        pass
                    try:
                        n_tertiary_in_bot += vww.turns[2]
                    except IndexError:
                        pass

        try:
            primary_in_top = curve_loop_cond[0][0:n_primary_in_top]
        except IndexError:
            pass
        try:
            secondary_in_top = curve_loop_cond[1][0:n_secondary_in_top]
        except IndexError:
            pass
        try:
            tertiary_in_top = curve_loop_cond[2][0:n_tertiary_in_top]
        except IndexError:
            pass
        try:
            primary_in_bot = curve_loop_cond[0][n_primary_in_top:n_primary_in_top + n_primary_in_bot]
        except IndexError:
            pass
        try:
            secondary_in_bot = curve_loop_cond[1][n_secondary_in_top:n_secondary_in_top + n_secondary_in_bot]
        except IndexError:
            pass
        try:
            tertiary_in_bot = curve_loop_cond[2][n_tertiary_in_top:n_tertiary_in_top + n_tertiary_in_bot]
        except IndexError:
            pass

        # top window
        l_air_top = l_core_air[0:5] + [l_bound_air[0]]
        curve_loop_air_top = [gmsh.model.geo.addCurveLoop(l_air_top, -1, True)]
        flatten_curve_loop_cond_top = primary_in_top + secondary_in_top + tertiary_in_top
        curve_loop_iso_core_top = []  # TODO: insulations
        self.plane_surface_air_top.append(gmsh.model.geo.addPlaneSurface(curve_loop_air_top + flatten_curve_loop_cond_top + curve_loop_iso_core_top))

        # bot window
        l_air_bot = l_core_air[5:12] + [l_bound_air[1]]
        curve_loop_air_bot = [gmsh.model.geo.addCurveLoop(l_air_bot, -1, True)]
        flatten_curve_loop_cond_bot = primary_in_bot + secondary_in_bot + tertiary_in_bot
        curve_loop_iso_core_bot = []  # TODO: insulations
        self.plane_surface_air_bot.append(gmsh.model.geo.addPlaneSurface(curve_loop_air_bot + flatten_curve_loop_cond_bot + curve_loop_iso_core_bot))

        # TODO: How to select the conductors which are in the top and which are in the bot vww? -> Need to be cut out of the air...
        # l_air_top = l_core_air[0:5] + [l_bound_air[0]]
        # curve_loop_air_top = [gmsh.model.geo.addCurveLoop(l_air_top, -1, True)]
        # l_air_bot = l_core_air[5:12] + [l_bound_air[1]]
        # curve_loop_air_bot = [gmsh.model.geo.addCurveLoop(l_air_bot, -1, True)]
        #
        # flatten_curve_loop_cond = [j for sub in curve_loop_cond for j in sub]
        # curve_loop_iso_core_bot = []  # TODO: insulations
        # self.plane_surface_air_bot.append(
        #     gmsh.model.geo.addPlaneSurface(curve_loop_air_bot + flatten_curve_loop_cond_bot + curve_loop_iso_core_bot))

    def boundary(self, p_core: list, p_region: list, l_bound_core: list, l_bound_air: list, l_region: list):
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Boundary
        if self.region is None:
            self.l_bound_tmp = l_bound_core[:5]
            for i in range(len(l_bound_air)):
                self.l_bound_tmp.append(l_bound_air[-i - 1])
                if i != len(l_bound_air) - 1:  # last run
                    self.l_bound_tmp.append(l_bound_core[-i - 1])

        else:
            # Generate Lines of Region
            # start top left and go clockwise
            p_region.append(gmsh.model.geo.addPoint(0,
                                                    self.model.p_region_bound[2][1],
                                                    self.model.p_region_bound[2][2],
                                                    self.model.p_region_bound[2][3]))

            p_region.append(gmsh.model.geo.addPoint(self.model.p_region_bound[3][0],
                                                    self.model.p_region_bound[3][1],
                                                    self.model.p_region_bound[3][2],
                                                    self.model.p_region_bound[3][3]))

            p_region.append(gmsh.model.geo.addPoint(self.model.p_region_bound[1][0],
                                                    self.model.p_region_bound[1][1],
                                                    self.model.p_region_bound[1][2],
                                                    self.model.p_region_bound[1][3]))

            p_region.append(gmsh.model.geo.addPoint(0,
                                                    self.model.p_region_bound[0][1],
                                                    self.model.p_region_bound[0][2],
                                                    self.model.p_region_bound[0][3]))

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
            self.l_bound_tmp = [l_bound_core[4]]

            for i in range(len(l_region)):
                self.l_bound_tmp.append(l_region[i])

            self.l_bound_tmp.append(l_bound_core[0])

            for i in range(len(l_bound_air)):
                self.l_bound_tmp.append(l_bound_air[-i - 1])
                if i != len(l_bound_air) - 1:  # last run
                    self.l_bound_tmp.append(l_bound_core[-i - 1])

            # Outer Air Surface
            curve_loop_outer_air = gmsh.model.geo.addCurveLoop(l_region + l_bound_core[1:4])
            self.plane_surface_outer_air.append(gmsh.model.geo.addPlaneSurface([curve_loop_outer_air]))

    def visualize(self, visualize_before, save_png):
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Colorize model and show it if needed
        # Create mesh
        color_scheme = ff.colors_femmt_default
        colors_geometry = ff.colors_geometry_femmt_default

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
        if self.model.core.core_type == CoreType.Single:
            gmsh.model.setColor([(2, self.plane_surface_air[0])], color_scheme[colors_geometry["potting_inner"]][0],
                                color_scheme[colors_geometry["potting_inner"]][1], color_scheme[colors_geometry["potting_inner"]][2], recursive=True)
        elif self.model.core.core_type == CoreType.Stacked:
            gmsh.model.setColor([(2, self.plane_surface_air_bot[0])], color_scheme[colors_geometry["potting_inner"]][0],
                                color_scheme[colors_geometry["potting_inner"]][1], color_scheme[colors_geometry["potting_inner"]][2], recursive=True)
            gmsh.model.setColor([(2, self.plane_surface_air_top[0])], color_scheme[colors_geometry["potting_inner"]][0],
                                color_scheme[colors_geometry["potting_inner"]][1], color_scheme[colors_geometry["potting_inner"]][2], recursive=True)


        # winding colors
        for winding_number in range(len(self.windings)):
            for turn_number in range(len(self.plane_surface_cond[winding_number])):
                gmsh.model.setColor([(2, self.plane_surface_cond[winding_number][turn_number])],
                                    color_scheme[colors_geometry["winding"][winding_number]][0], color_scheme[colors_geometry["winding"][winding_number]][1],
                                    color_scheme[colors_geometry["winding"][winding_number]][2], recursive=True)

        # insulation color (inner insulation / bobbin)
        gmsh.model.setColor([(2, iso) for iso in self.plane_surface_iso_core],
                            color_scheme[colors_geometry["insulation"]][0], color_scheme[colors_geometry["insulation"]][1], color_scheme[colors_geometry["insulation"]][2], recursive=True)

        if visualize_before:
            gmsh.fltk.run()

        # Output .msh file
        # gmsh.option.setNumber("Mesh.SaveAll", 1)
        if save_png:
            gmsh.fltk.initialize()

            gmsh.write(self.hybrid_color_png_file)  # save png

    def generate_hybrid_mesh(self, color_scheme: Dict = ff.colors_femmt_default, colors_geometry: Dict = ff.colors_geometry_femmt_default,
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
        ff.femmt_print("Hybrid Mesh Generation in Gmsh")
        # Initialization
        self.set_empty_plane_lists()
        p_core, p_island, p_cond, p_region, p_iso_core = self.set_empty_point_lists()
        l_bound_core, l_bound_air, l_core_air, l_cond, l_region, l_air_gaps_air, l_iso_core = self.set_empty_line_lists()
        curve_loop_cond, curve_loop_island, curve_loop_air, curve_loop_air_gaps, curve_loop_iso_core = self.set_empty_curve_loop_lists()

        # Set path for storing the mesh file
        gmsh.model.add(os.path.join(self.e_m_mesh_file, "geometry"))

        # Define mesh for core
        if self.core.core_type == CoreType.Single:
            self.single_core(p_core, p_island,
                             l_bound_core, l_core_air, l_bound_air, l_air_gaps_air,
                             curve_loop_island, curve_loop_air_gaps)
        if self.core.core_type == CoreType.Stacked:
            self.stacked_core(p_core, l_bound_core, l_core_air, l_bound_air)

        # Define mesh for conductors
        self.conductors(p_cond, l_cond, curve_loop_cond)

        # Define mesh for conductors
        model_insulation: bool = False
        if model_insulation:
            self.insulations_core_cond(p_iso_core)

        # Define mesh for air
        if self.core.core_type == CoreType.Single:
            self.air_single(l_core_air, l_air_gaps_air, curve_loop_air, curve_loop_cond, curve_loop_iso_core)
        if self.core.core_type == CoreType.Stacked:
            self.air_stacked(l_core_air, l_bound_air, curve_loop_cond)

        # Define mesh for boundary
        self.boundary(p_core, p_region, l_bound_core, l_bound_air, l_region)

        gmsh.model.geo.synchronize()

        self.visualize(visualize_before, save_png)


        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # TODO The following algorithms try to modify the mesh in order to reduce the runtime. But maybe the synchronize() calls
        # have a high runtime. Check if thats true and when it does try to reduce the number of synchronize() calls by adding all points first and
        # embed them later together:
        # This is added here therefore the additional points are not seen in the pictures and views
        # self.forward_meshing(p_cond)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # No mesh is generated here, because generating a mesh, saving it as *.msh, loading it, appending more geometry data
        # and then mesh again can cause bugs in the mesh
        # Therefore only the model geometry is saved and the mesh will be generated later
        # -> Save file as geo: File extension must be *.geo_unrolled
        gmsh.write(self.model_geo_file)

    def generate_electro_magnetic_mesh(self, refine=0):
        ff.femmt_print("Electro Magnetic Mesh Generation in Gmsh (write physical entities)")

        gmsh.open(self.model_geo_file)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Define physical Surfaces and Curves
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Core
        self.ps_core = gmsh.model.geo.addPhysicalGroup(2, self.plane_surface_core, tag=120000)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Conductors
        self.ps_cond = []
        for num in range(len(self.windings)):
            self.ps_cond.append([])

        # Since the turns are saved for each vww all turns per winding must be collected
        flattened_turns = [0] * len(self.windings)
        for ww in self.model.winding_windows:
            for vww in ww.virtual_winding_windows:
                for index, winding in enumerate(self.windings):
                    # print(f"{index = }")
                    # print(f"{ vww.turns = }")
                    # TODO: change index_turns right no. of winding numbers, right position in list and length of list is needed
                    try:
                        flattened_turns[winding.winding_number] += vww.turns[index]
                    except:
                        pass

        for winding in self.windings:
            winding_number = winding.winding_number
            if winding.conductor_type == ConductorType.RoundLitz:
                for i in range(flattened_turns[winding_number]):
                    self.ps_cond[winding_number].append(
                        gmsh.model.geo.addPhysicalGroup(2, [self.plane_surface_cond[winding_number][i]], tag=140000 + 1000 * winding_number + i))

            else:
                if winding.parallel:
                    tags = self.plane_surface_cond[winding_number]
                    physical_group_number = gmsh.model.geo.addPhysicalGroup(2, tags, tag=130000 + 1000 * winding_number)
                    for i in range(flattened_turns[winding_number]):
                        self.ps_cond[winding_number].append(physical_group_number)
                else:
                    for i in range(flattened_turns[winding_number]):
                        self.ps_cond[winding_number].append(
                            gmsh.model.geo.addPhysicalGroup(2, [self.plane_surface_cond[winding_number][i]], tag=130000 + 1000 * winding_number + i))

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Air, air_gaps and iso (since insulation is handled as air, as well as the air gaps)
        if self.model.core.core_type == CoreType.Single:
            air_and_air_gaps = self.plane_surface_air + self.plane_surface_air_gaps + self.plane_surface_iso_core
            self.ps_air = gmsh.model.geo.addPhysicalGroup(2, air_and_air_gaps, tag=110000)
            # ps_air_ext = gmsh.model.geo.addPhysicalGroup(2, plane_surface_outer_air, tag=1001)
        elif self.model.core.core_type == CoreType.Stacked:
            air_total = self.plane_surface_air_bot + self.plane_surface_air_top
            self.ps_air = gmsh.model.geo.addPhysicalGroup(2, air_total, tag=110000)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Boundary
        self.pc_bound = gmsh.model.geo.addPhysicalGroup(1, self.l_bound_tmp, tag=111111)
        # print(f"Physical Conductor Surfaces: {self.ps_cond}")

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Set names [optional]
        gmsh.model.setPhysicalName(2, self.ps_core, "CORE")
        for num in range(len(self.windings)):
            for i in range(len(self.ps_cond[num])):
                gmsh.model.setPhysicalName(2, self.ps_cond[num][i], f"COND{num + 1}")
        gmsh.model.setPhysicalName(2, self.ps_air, "AIR")
        gmsh.model.setPhysicalName(1, self.pc_bound, "BOUND")

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
            ff.femmt_print("\n ------- \nRefined Mesh Creation ")
            # mesh the new gmsh.model using the size field
            bg_field = gmsh.model.mesh.field.add("PostView")
            # TODO: gmsh.model.mesh.field.setNumber(bg_field, "ViewTag", sf_view)
            gmsh.model.mesh.field.setAsBackgroundMesh(bg_field)
            ff.femmt_print("\nMeshing...\n")
            gmsh.model.mesh.generate(2)
        else:
            ff.femmt_print("\nMeshing...\n")
            gmsh.model.mesh.generate(2)

        if not os.path.exists(self.mesh_folder_path):
            os.mkdir(self.mesh_folder_path)

        gmsh.write(self.e_m_mesh_file)

    def generate_thermal_mesh(self, case_gap_top, case_gap_right, case_gap_bot, color_scheme, colors_geometry, visualize_before):
        ff.femmt_print("Thermal Mesh Generation in Gmsh (write physical entities)")

        gmsh.open(self.model_geo_file)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Create case around the core

        # Core point and line tags
        core_point_tags = []
        core_line_tags = []
        if self.model.air_gaps.number == 0:
            core_point_tags = [4, 3, 2, 1]
            core_line_tags = [4, 3, 2]
        else:
            core_point_tags = [5, 4, 3, 2]
            core_line_tags = [5, 4, 3]

        tl_point = core_point_tags[0]  # Top left - default 5
        tr_point = core_point_tags[1]  # Top right - default 4
        br_point = core_point_tags[2]  # Bottom right - default 3
        bl_point = core_point_tags[3]  # Bottom left - default 2

        top_line = core_line_tags[0]  # default 4
        right_line = core_line_tags[1]  # default 3
        bottom_line = core_line_tags[2]  # default 2

        # Get positions from points
        tl_point_pos = gmsh.model.getValue(0, tl_point, [])
        tr_point_pos = gmsh.model.getValue(0, tr_point, [])
        br_point_pos = gmsh.model.getValue(0, br_point, [])
        bl_point_pos = gmsh.model.getValue(0, bl_point, [])

        mesh = self.mesh_data.c_core * 4  # It typically does not need to be the same size as c_core, but it shouldn't be too big either

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
        self.ps_core = gmsh.model.geo.addPhysicalGroup(2, self.plane_surface_core, tag=120000)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Conductors
        self.ps_cond = []
        for num in range(len(self.windings)):
            self.ps_cond.append([])

        # Since the turns are saved for each vww all turns per winding must be collected
        flattened_turns = [0] * len(self.windings)
        for vww in self.model.virtual_winding_windows:
            for index, winding in enumerate(self.windings):
                flattened_turns[winding.winding_number] += vww.turns[index]

        for winding in self.windings:
            winding_number = winding.winding_number
            if winding.conductor_type == ConductorType.RoundLitz:
                for i in range(flattened_turns[winding_number]):
                    self.ps_cond[winding_number].append(
                        gmsh.model.geo.addPhysicalGroup(2, [self.plane_surface_cond[winding_number][i]], tag=140000 + 1000 * winding_number + i))
            else:
                for i in range(flattened_turns[winding_number]):
                    self.ps_cond[winding_number].append(
                        gmsh.model.geo.addPhysicalGroup(2, [self.plane_surface_cond[winding_number][i]], tag=130000 + 1000 * winding_number + i))

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Air
        self.ps_air = gmsh.model.geo.addPhysicalGroup(2, self.plane_surface_air, tag=110000)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Air gaps
        self.ps_air_gaps = gmsh.model.geo.addPhysicalGroup(2, self.plane_surface_air_gaps, tag=110001)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # insulations
        # TODO Currently insulations can only have the same material
        self.ps_insulation = gmsh.model.geo.addPhysicalGroup(2, self.plane_surface_iso_core)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Boundary
        self.thermal_boundary_region_tags = {
            "BOUNDARY_TOP": top_case_top_line,
            "BOUNDARY_TOP_RIGHT": top_right_case_top_line,
            "BOUNDARY_RIGHT_TOP": top_right_case_right_line,
            "BOUNDARY_RIGHT": right_case_right_line,
            "BOUNDARY_RIGHT_BOTTOM": bottom_right_case_right_line,
            "BOUNDARY_BOTTOM_RIGHT": bottom_right_case_bottom_line,
            "BOUNDARY_BOTTOM": bottom_case_bottom_line
        }

        for key in self.thermal_boundary_region_tags:
            self.thermal_boundary_region_tags[key] = ff.create_physical_group(1, [self.thermal_boundary_region_tags[key]], key)

        # Add surface physical groups
        # INFO: The physical groups are not created in the createRectWithPhysicalGroup because it causes a bug with the index counter when
        # 1D physical groups (lines) are added after 2D physical groups (surfaces)
        top_surface_physical_group = ff.create_physical_group(2, [top_case_surface], "TopCase")
        top_right_surface_physical_group = ff.create_physical_group(2, [top_right_case_surface], "TopRightCase")
        right_surface_physical_group = ff.create_physical_group(2, [right_case_surface], "RightCase")
        bottom_right_surface_physical_group = ff.create_physical_group(2, [bottom_right_case_surface], "BottomRightCase")
        bottom_surface_physical_group = ff.create_physical_group(2, [bottom_case_surface], "BottomCase")

        self.thermal_boundary_ps_groups = [top_surface_physical_group, top_right_surface_physical_group,
                                           right_surface_physical_group, bottom_right_surface_physical_group, bottom_surface_physical_group]

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Set names [optional]
        gmsh.model.setPhysicalName(2, self.ps_core, "CORE")
        for num in range(len(self.windings)):
            for i in range(len(self.ps_cond[num])):
                gmsh.model.setPhysicalName(2, self.ps_cond[num][i], f"COND{num + 1}")
        gmsh.model.setPhysicalName(2, self.ps_air, "AIR")
        gmsh.model.setPhysicalName(2, self.ps_air_gaps, "AIR_GAPS")
        gmsh.model.setPhysicalName(2, self.ps_insulation, "insulationS")

        # Synchronize
        gmsh.model.geo.synchronize()

        # Set colors:
        color_case = color_scheme[colors_geometry["potting_outer"]]
        color_background = color_scheme[colors_geometry["potting_inner"]]
        color_core = color_scheme[colors_geometry["core"]]
        color_air_gap = color_scheme[colors_geometry["air_gap"]]
        color_insulation = color_scheme[colors_geometry["insulation"]]

        # Set case color to core color
        for sf in [top_case_surface, top_right_case_surface, right_case_surface, bottom_right_case_surface, bottom_case_surface]:
            gmsh.model.setColor([(2, sf)], color_case[0], color_case[1], color_case[2], recursive=True)

        # core color
        for i in range(len(self.plane_surface_core)):
            gmsh.model.setColor([(2, self.plane_surface_core[i])], color_core[0], color_core[1], color_core[2], recursive=True)

        # air gap color
        if self.plane_surface_air_gaps:
            # only colorize air-gap in case of air gaps
            gmsh.model.setColor([(2, self.plane_surface_air[0]), (2, self.plane_surface_air_gaps[0])], color_air_gap[0], color_air_gap[1],
                                color_air_gap[2], recursive=True)

        # air/potting-material inside core window
        gmsh.model.setColor([(2, self.plane_surface_air[0])], color_background[0], color_background[1], color_background[2], recursive=True)

        # winding colors
        for winding_number in range(len(self.windings)):
            for turn_number in range(len(self.plane_surface_cond[winding_number])):
                gmsh.model.setColor([(2, self.plane_surface_cond[winding_number][turn_number])], color_scheme[colors_geometry["winding"][winding_number]][0],
                                    color_scheme[colors_geometry["winding"][winding_number]][1], color_scheme[colors_geometry["winding"][winding_number]][2], recursive=True)

        # insulation color (inner insulation / bobbin)
        gmsh.model.setColor([(2, iso) for iso in self.plane_surface_iso_core], color_insulation[0], color_insulation[1],
                            color_insulation[2], recursive=True)

        if visualize_before:
            gmsh.fltk.run()

        # Output .msh file
        gmsh.model.mesh.generate(2)

        if not os.path.exists(self.mesh_folder_path):
            os.mkdir(self.mesh_folder_path)

        gmsh.write(self.thermal_mesh_file)

    def inter_conductor_meshing(self, p_cond):
        p_inter = None
        for ww in self.model.winding_windows:
            for vww in ww.virtual_winding_windows:
                if vww.winding_type != WindingType.TwoInterleaved:
                    for index, winding in enumerate(vww.windings):    #TODO: proof
                        num = winding.winding_number
                        p_inter = []
                        x_inter = []
                        y_inter = []
                        j = 0

                        if winding.conductor_type == ConductorType.RoundSolid and \
                                vww.turns[index] > 1:
                            while self.model.p_conductor[num][5 * j][1] == \
                                    self.model.p_conductor[num][5 * j + 5][1]:
                                x_inter.append(
                                    0.5 * (self.model.p_conductor[num][5 * j][0] +
                                           self.model.p_conductor[num][5 * j + 5][0]))
                                j += 1
                                if j == vww.turns[index] - 1:
                                    break
                            j += 1
                            if int(vww.turns[index] / j) > 1:
                                for i in range(0, int(vww.turns[index] / j)):
                                    if 5 * j * i + 5 * j >= len(self.model.p_conductor[num][:]):
                                        break
                                    y_inter.append(0.5 * (self.model.p_conductor[num][5 * j * i][1] +
                                                          self.model.p_conductor[num][5 * j * i + 5 * j][1]))
                                for x in x_inter:
                                    for y in y_inter:
                                        p_inter.append(gmsh.model.geo.addPoint(x,
                                                                               y,
                                                                               0,
                                                                               self.mesh_data.c_center_conductor[num]))

        # TODO: Inter conductor meshing!
        if all(winding.conductor_type == ConductorType.RoundSolid for winding in self.windings):
            ff.femmt_print(f"Making use of skin based meshing\n")
            for num in range(len(self.windings)):
                for i in range(0, int(len(p_cond[num]) / 5)):
                    gmsh.model.mesh.embed(0, [p_cond[num][5 * i + 0]], 2, self.plane_surface_cond[num][i])

            # Embed points for mesh refinement
            # Inter Conductors
            for ww in self.model.winding_windows:
                for vww in ww.virtual_winding_windows:
                    if vww.winding_type != WindingType.TwoInterleaved:
                        gmsh.model.mesh.embed(0, p_inter, 2, self.plane_surface_air[0])
            # Stray path
            # mshopt gmsh.model.mesh.embed(0, stray_path_mesh_optimizer, 2, plane_surface_core[2])

    def forward_meshing(self, p_cond):
        """In this function multiple techniques in order to raise the mesh density at certain points are applied.
        :return:
        """
        def rasterize_winding_window(left_bound, right_bound, bot_bound, top_bound):

            # Winding window rasterization:
            # In order adjust the mesh density in empty parts of the winding window a grid of possible points
            # is put on the winding window. Every point that is too close to the conductors is removed.
            # Every remaining point is added to the mesh with a higher mesh density

            # min_distance = max([winding.conductor_radius for winding in self.windings]) + max(self.insulation.inner_winding_insulations)
            # min_distance = max(self.insulation.inner_winding_insulations)
            min_distance = 0  # TODO: MA Project?

            width = right_bound - left_bound
            height = top_bound - bot_bound

            number_cols = 17  # Can be changed. More points equal higher raster density
            number_rows = int(number_cols * height / width)  # Assumption: number_cols/number_rows = width/height

            cell_width = width / (number_cols + 1)
            cell_height = height / (number_rows + 1)

            # Get all possible points
            possible_points = []
            x = left_bound + cell_width / 2
            y = bot_bound + cell_height / 2
            for i in range(number_cols + 1):
                for j in range(number_rows + 1):
                    possible_points.append([x + i * cell_width, y + j * cell_height])

            fixed_points = []
            conductors = self.model.p_conductor
            for winding in range(len(self.windings)):
                for i in range(len(conductors[winding]) // 5):
                    point = conductors[winding][i * 5]
                    fixed_points.append([point[0], point[1]])

            # Because the points need to be embed into the right surface. The points now will be split between different insulations and the air in the winding window.
            # TODO Currently primary secondary insulation is not considered
            left_iso = []
            right_iso = []
            top_iso = []
            bot_iso = []
            air = []

            # insulations are currently not implemented for integrated transformers
            if self.component_type != ComponentType.IntegratedTransformer:
                iso_core_left = self.model.p_iso_core[0]
                iso_core_top = self.model.p_iso_core[1]
                iso_core_right = self.model.p_iso_core[2]
                iso_core_bot = self.model.p_iso_core[3]

            # Extract all free_points
            for i in range(len(possible_points)):
                x = possible_points[i][0]
                y = possible_points[i][1]

                # Check collision with fixed points
                valid = True
                for fixed_point in fixed_points:
                    dist = np.sqrt((fixed_point[0] - x) ** 2 + (fixed_point[1] - y) ** 2)
                    if dist < min_distance:
                        valid = False
                        break

                if not valid:
                    continue

                # Check if point is in stray_path
                if self.component_type == ComponentType.IntegratedTransformer and self.core.core_type == CoreType.Single:
                    start_index = self.stray_path.start_index
                    stray_path_top_bound = self.air_gaps.midpoints[start_index + 1][1] - self.air_gaps.midpoints[start_index + 1][2] / 2
                    stray_path_bot_bound = self.air_gaps.midpoints[start_index][1] + self.air_gaps.midpoints[start_index][2] / 2
                    stray_path_right_bound = self.stray_path.length
                    stray_path_left_bound = left_bound

                    if x > stray_path_left_bound and x < stray_path_right_bound and y > stray_path_bot_bound and y < stray_path_top_bound:
                        continue

                # Point seems to be valid. Now find out in which surface the point belongs
                point = gmsh.model.geo.addPoint(x, y, 0, 1 * self.mesh_data.c_window)

                if self.component_type != ComponentType.IntegratedTransformer:
                    if ff.point_is_in_rect(x, y, iso_core_left):
                        # Left iso
                        left_iso.append(point)
                    elif ff.point_is_in_rect(x, y, iso_core_top):
                        # Top iso
                        top_iso.append(point)
                    elif ff.point_is_in_rect(x, y, iso_core_right):
                        # Right iso
                        right_iso.append(point)
                    elif ff.point_is_in_rect(x, y, iso_core_bot):
                        # Bot iso
                        bot_iso.append(point)
                    else:
                        # Air
                        air.append(point)
                else:
                    air.append(point)
            # Call synchronize so the points will be added to the model
            gmsh.model.geo.synchronize()

            # Embed points into surfaces
            # if self.component_type != ComponentType.IntegratedTransformer:
            #     gmsh.model.mesh.embed(0, left_iso, 2, self.plane_surface_iso_core[0])
            #     gmsh.model.mesh.embed(0, top_iso, 2, self.plane_surface_iso_core[1])
            #     gmsh.model.mesh.embed(0, right_iso, 2, self.plane_surface_iso_core[2])
            #     gmsh.model.mesh.embed(0, bot_iso, 2, self.plane_surface_iso_core[3])
            return air

        if self.core.core_type == CoreType.Single:
            # Inter Conductors
            self.inter_conductor_meshing(p_cond)

            left_bound = self.core.core_inner_diameter / 2
            right_bound = self.model.r_inner
            top_bound = self.core.window_h / 2
            bot_bound = -self.core.window_h / 2

            air_tags = rasterize_winding_window(left_bound, right_bound, bot_bound, top_bound)
            gmsh.model.mesh.embed(0, air_tags, 2, self.plane_surface_air[0])

        if self.core.core_type == CoreType.Stacked:
            # Inter Conductors
            self.inter_conductor_meshing(p_cond)

            # Top Window
            top_bound = self.model.p_window_top[2][1]  # y component of top window
            bot_bound = self.model.p_window_top[0][1]  # y component of bot window
            left_bound = self.core.core_inner_diameter / 2
            right_bound = self.model.r_inner

            air_tags = rasterize_winding_window(left_bound, right_bound, bot_bound, top_bound)
            gmsh.model.mesh.embed(0, air_tags, 2, self.plane_surface_air_top[0])

            # Bot Window
            top_bound = self.model.p_window_bot[2][1]  # y component of top window
            bot_bound = self.model.p_window_bot[0][1]  # y component of bot window
            left_bound = self.core.core_inner_diameter / 2
            right_bound = self.model.r_inner

            air_tags = rasterize_winding_window(left_bound, right_bound, bot_bound, top_bound)
            gmsh.model.mesh.embed(0, air_tags, 2, self.plane_surface_air_bot[0])

        # self.visualize(visualize_before=True, save_png=False)

        # Synchronize again
        gmsh.model.geo.synchronize()
