"""Draw structures inside Onelab."""
# Python standard libraries
import numpy as np
from logging import Logger
from typing import List

# Local libraries
from femmt.enumerations import *
from femmt.data import MeshData
from femmt.model import Core, WindingWindow, AirGaps, StrayPath, Insulation


class TwoDaxiSymmetric:
    """
    Create the model of the magnetic component.

    This is done by creating lists of points which will be later added to the gmsh file.
    """

    core: Core
    winding_windows: List[WindingWindow]
    air_gaps: AirGaps
    stray_path: StrayPath
    insulation: Insulation
    component_type: ComponentType
    mesh_data: MeshData
    number_of_windings: int
    verbosity: Verbosity
    logger: Logger

    # List of points which represent the model
    # Every List is a List of 4 Points: x, y, z, mesh_factor
    p_outer: np.ndarray
    p_region_bound: np.ndarray
    p_window: np.ndarray
    p_air_gaps: np.ndarray
    # p_conductor: List[List[float]]
    p_iso_core: List[List[float]]
    p_iso_top_core: List[List[float]]
    p_iso_bot_core: List[List[float]]
    p_iso_pri_sec: List[List[float]]

    def __init__(self, core: Core, mesh_data: MeshData, air_gaps: AirGaps, winding_windows: List[WindingWindow],
                 stray_path: StrayPath, insulation: Insulation, component_type: ComponentType, number_of_windings: int, verbosity: Verbosity, logger: Logger):
        self.core = core
        self.mesh_data = mesh_data
        self.winding_windows = winding_windows
        self.air_gaps = air_gaps
        self.component_type = component_type
        self.stray_path = stray_path
        self.insulation = insulation
        self.number_of_windings = number_of_windings
        self.verbosity = verbosity
        self.logger = logger

        # -- Arrays for geometry data -- 
        # TODO Is the zero initialization necessary?
        self.p_outer = np.zeros((4, 4))
        self.p_region_bound = np.zeros((4, 4))
        if self.core.core_type == CoreType.Single:
            self.p_window = np.zeros((4 * core.number_core_windows, 4))  # TODO: why is this done for both sides?
        if self.core.core_type == CoreType.Stacked:
            self.p_window_top = np.zeros((4, 4))  # TODO: just for right side? make it the same as for single core geometry
            self.p_window_bot = np.zeros((4, 4))  # TODO: just for right side? make it the same as for single core geometry
        self.p_air_gaps = np.zeros((4 * air_gaps.number, 4))
        self.p_conductor = []
        self.p_iso_core = []
        self.p_iso_top_core = []
        self.p_iso_bot_core = []
        self.p_iso_pri_sec = []

        for i in range(number_of_windings):
            self.p_conductor.insert(i, [])

        self.r_inner = core.r_inner
        self.r_outer = core.r_outer

    def femmt_print(self, text: str):
        """
        Print text to terminal or to log-file, dependent on the current verbosity.

        :param text: text to print
        :type text: str
        """
        if not self.verbosity == Verbosity.Silent:
            self.logger.info(text)

    def draw_outer(self):
        """Draws the outer points of the main core (single core)."""
        # Outer Core
        # (A_zyl=2pi*r*h => h=0.5r=0.25core_w <=> ensure A_zyl=A_core on the tiniest point)
        self.p_outer[0][:] = [-self.r_outer, -self.core.core_h / 2, 0, self.mesh_data.c_core]

        self.p_outer[1][:] = [self.r_outer, -self.core.core_h / 2, 0, self.mesh_data.c_core]

        self.p_outer[2][:] = [-self.r_outer, self.core.core_h / 2, 0, self.mesh_data.c_core]

        self.p_outer[3][:] = [self.r_outer, self.core.core_h / 2, 0, self.mesh_data.c_core]

    def draw_single_window(self):
        """Draw a single window."""
        # At this point both windows (in a cut) are modeled
        # print(f"win: c_window: {self.component.mesh.c_window}")
        self.p_window[0] = [-self.r_inner,
                            -self.core.window_h / 2,
                            0,
                            self.mesh_data.c_window]

        self.p_window[1] = [-self.core.core_inner_diameter / 2,
                            -self.core.window_h / 2,
                            0,
                            self.mesh_data.c_window]

        self.p_window[2] = [-self.r_inner,
                            self.core.window_h / 2,
                            0,
                            self.mesh_data.c_window]

        self.p_window[3] = [-self.core.core_inner_diameter / 2,
                            self.core.window_h / 2,
                            0,
                            self.mesh_data.c_window]

        self.p_window[4] = [self.core.core_inner_diameter / 2,
                            -self.core.window_h / 2,
                            0,
                            self.mesh_data.c_window]

        self.p_window[5] = [self.r_inner,
                            -self.core.window_h / 2,
                            0,
                            self.mesh_data.c_window]

        self.p_window[6] = [self.core.core_inner_diameter / 2,
                            self.core.window_h / 2,
                            0,
                            self.mesh_data.c_window]

        self.p_window[7] = [self.r_inner,
                            self.core.window_h / 2,
                            0,
                            self.mesh_data.c_window]

    def draw_stacked_windows(self):
        """Draw stacked windows."""
        self.p_window_bot[0] = [self.core.core_inner_diameter / 2,
                                -self.core.window_h_bot / 2,
                                0,
                                self.mesh_data.c_window]

        self.p_window_bot[1] = [self.r_inner,
                                -self.core.window_h_bot / 2,
                                0,
                                self.mesh_data.c_window]

        self.p_window_bot[2] = [self.core.core_inner_diameter / 2,
                                self.core.window_h_bot / 2,
                                0,
                                self.mesh_data.c_window]

        self.p_window_bot[3] = [self.r_inner,
                                self.core.window_h_bot / 2,
                                0,
                                self.mesh_data.c_window]

        self.p_window_top[0] = [self.core.core_inner_diameter / 2,
                                self.core.window_h_bot / 2 + self.core.core_thickness,
                                0,
                                self.mesh_data.c_window]

        self.p_window_top[1] = [self.r_inner,
                                self.core.window_h_bot / 2 + self.core.core_thickness,
                                0,
                                self.mesh_data.c_window]

        self.p_window_top[2] = [self.core.core_inner_diameter / 2,
                                self.core.window_h_bot / 2 + self.core.core_thickness + self.core.window_h_top,
                                0,
                                self.mesh_data.c_window]

        self.p_window_top[3] = [self.r_inner,
                                self.core.window_h_bot / 2 + self.core.core_thickness + self.core.window_h_top,
                                0,
                                self.mesh_data.c_window]

    def draw_air_gaps(self):
        """Draw air gaps."""
        # Air gaps
        # "air_gaps" is a list with [position_tag, air_gap_position, air_gap_h, c_air_gap]
        #   - position_tag: specifies the "leg" with air gaps
        #   - air_gap_position: specifies the coordinate of the air gap's center point along the specified leg
        #   - air_gap_h: height/length of the air gap
        #   - c_air_gap: mesh accuracy factor
        # at this point the 4 corner points of each air gap are generated out of "air_gaps"

        for i in range(0, self.air_gaps.number):

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
            if self.air_gaps.midpoints[i][0] == 0:
                # The center points are transformed each into 4 corner points

                air_gap_y_position = self.air_gaps.midpoints[i][1]
                air_gap_height = self.air_gaps.midpoints[i][2]
                air_gap_length_top = self.core.core_inner_diameter / 2
                air_gap_length_bot = self.core.core_inner_diameter / 2

                # Check for stray_paths in integrated transformers
                if self.component_type == ComponentType.IntegratedTransformer:
                    if self.stray_path.start_index == i:
                        # Stray path is above current air_gap
                        air_gap_length_top = self.stray_path.length
                    elif self.stray_path.start_index + 1 == i:
                        # Stray path is below current air_gap
                        air_gap_length_bot = self.stray_path.length

                # Bottom left
                self.p_air_gaps[i * 4 + 0] = [0, air_gap_y_position - air_gap_height / 2, 0, self.mesh_data.c_air_gaps]

                # Bottom right
                self.p_air_gaps[i * 4 + 1] = [air_gap_length_bot, air_gap_y_position - air_gap_height / 2, 0, self.mesh_data.c_air_gaps]

                # Top left
                self.p_air_gaps[i * 4 + 2] = [0, air_gap_y_position + air_gap_height / 2, 0, self.mesh_data.c_air_gaps]

                # Top right
                self.p_air_gaps[i * 4 + 3] = [air_gap_length_top, air_gap_y_position + air_gap_height / 2, 0, self.mesh_data.c_air_gaps]

        # In order to close the air gap when a stray_path is added, additional points need to be added
        if self.component_type == ComponentType.IntegratedTransformer:
            top_point = [self.core.core_inner_diameter / 2,
                         self.air_gaps.midpoints[self.stray_path.start_index + 1][1] - self.air_gaps.midpoints[self.stray_path.start_index + 1][2] / 2,
                         0, self.mesh_data.c_air_gaps]
            bot_point = [self.core.core_inner_diameter / 2,
                         self.air_gaps.midpoints[self.stray_path.start_index][1] + self.air_gaps.midpoints[self.stray_path.start_index][2] / 2,
                         0, self.mesh_data.c_air_gaps]
            self.p_close_air_gaps = [top_point, bot_point]

    def draw_conductors(self, draw_top_down: bool = True):
        """
        Draw every conductor type based on the virtual_winding_window bounds.

        :param draw_top_down: True to draw from top to bottom
        :type draw_top_down: bool
        """
        for winding_window in self.winding_windows:
            for virtual_winding_window in winding_window.virtual_winding_windows:
                # Get bounds from virtual winding window
                bot_bound = virtual_winding_window.bot_bound
                top_bound = virtual_winding_window.top_bound
                left_bound = virtual_winding_window.left_bound
                right_bound = virtual_winding_window.right_bound

                # Check, if virtual winding window fits in physical window
                # print(f"{bot_bound = }\n", f"{top_bound = }\n", f"{left_bound = }\n", f"{right_bound = }\n")
                # print(f"{winding_window.max_bot_bound = }\n", f"{winding_window.max_top_bound = }\n", f"{winding_window.max_left_bound = }\n",
                # f"{winding_window.max_right_bound = }\n")
                if bot_bound < winding_window.max_bot_bound or top_bound > winding_window.max_top_bound or \
                        left_bound < winding_window.max_left_bound or right_bound > winding_window.max_right_bound:
                    # Set valid to False, so that geometry is to be neglected in geometry sweep
                    # self.valid = False
                    raise Exception("Winding does not fit into winding window!")
                else:
                    # Check the possible WindingTypes and draw accordingly
                    if virtual_winding_window.winding_type == WindingType.TwoInterleaved:
                        # Two windings in the virtual winding window
                        windings = virtual_winding_window.windings
                        winding0 = windings[0]
                        winding1 = windings[1]

                        turns = virtual_winding_window.turns
                        turns0 = turns[0]
                        turns1 = turns[1]

                        winding_numbers = [winding0.winding_number, winding1.winding_number]

                        # Now check for every winding scheme which is possible for interleaved winding type
                        if virtual_winding_window.winding_scheme == InterleavedWindingScheme.Bifilar:
                            """
                            - Bifilar interleaving means a uniform winding scheme of two conductors (prim. and sec.)
                            - Can only be used for conductors of identical radius (in terms of litz radius for 
                                stranded wires)
                            - Excess windings are placed below the bifilar ones
            
                            """
                            if winding0.conductor_radius != winding1.conductor_radius:
                                print("For bifilar winding scheme both conductors must be of the same radius!")
                            else:
                                print("Bifilar winding scheme is applied")

                            raise Exception("Bifilar winding scheme is not implemented yet.")

                        if virtual_winding_window.winding_scheme == InterleavedWindingScheme.VerticalAlternating:
                            """
                            - Vertical interleaving means a winding scheme where the two conductors are alternating 
                                in vertical (y-)direction
                            - This is practically uncommon
                            - If the turns ratio is != 1, the scheme always begins with the "higher-turns-number's" 
                                conductor
                            """

                            raise Exception("Vertical winding scheme is not implemented yet.")

                        if virtual_winding_window.winding_scheme == InterleavedWindingScheme.HorizontalAlternating:
                            """
                            - Horizontal interleaving means a winding scheme where the two conductors are alternating in 
                            horizontal  (x-)direction (Tonnenwicklung)
                            - This is practically most common
                            - If the turns ratio is != 1, the scheme always begins with the "higher-turns-number's" 
                                conductor
                                
                            """

                            # assume 2 winding transformer and dedicated stray path:
                            if draw_top_down:
                                # This will draw from top bound down

                                # Initialize the list, that counts the already placed conductors
                                N_completed = [0, 0]

                                # Initialize the starting conductor
                                if turns0 >= turns1:
                                    # Primary starts first
                                    col_cond = 0
                                else:
                                    # Secondary starts fist
                                    col_cond = 1

                                # Get winding number of current winding
                                winding_number = winding_numbers[col_cond]

                                # Initialize the x and y coordinate
                                x = left_bound + windings[col_cond].conductor_radius
                                y = top_bound - windings[col_cond].conductor_radius
                                top_window_iso_counter = 0
                                # Continue placing as long as not all conductors have been placed
                                while (turns0 - N_completed[0] != 0) or (turns1 - N_completed[1] != 0):
                                    if turns[col_cond] - N_completed[col_cond] != 0:
                                        # is this winding not already finished?
                                        if x < right_bound - windings[col_cond].conductor_radius:
                                            while y > bot_bound + windings[col_cond].conductor_radius and \
                                                    N_completed[col_cond] < turns[col_cond]:
                                                self.p_conductor[winding_number].append([
                                                    x,
                                                    y,
                                                    0,
                                                    self.mesh_data.c_center_conductor[winding_number]])

                                                self.p_conductor[winding_number].append([
                                                    x - windings[col_cond].conductor_radius,
                                                    y,
                                                    0,
                                                    self.mesh_data.c_conductor[winding_number]])

                                                self.p_conductor[winding_number].append([
                                                    x,
                                                    y + windings[col_cond].conductor_radius,
                                                    0,
                                                    self.mesh_data.c_conductor[winding_number]])

                                                self.p_conductor[winding_number].append([
                                                    x + windings[col_cond].conductor_radius,
                                                    y,
                                                    0,
                                                    self.mesh_data.c_conductor[winding_number]])

                                                self.p_conductor[winding_number].append([
                                                    x,
                                                    y - windings[col_cond].conductor_radius,
                                                    0,
                                                    self.mesh_data.c_conductor[winding_number]])

                                                N_completed[col_cond] += 1

                                                # one from top to bot
                                                y -= windings[col_cond].conductor_radius * 2 + self.insulation.cond_cond[col_cond][col_cond]

                                            # Check whether one winding is "finished" and only the other winding must be placed...
                                            if N_completed[(col_cond + 1) % 2] == turns[(col_cond + 1) % 2]:
                                                x += windings[col_cond].conductor_radius + windings[col_cond].conductor_radius + \
                                                    self.insulation.cond_cond[col_cond][(col_cond + 1) % 2]
                                            else:
                                                x += windings[col_cond].conductor_radius + windings[(col_cond + 1) % 2].conductor_radius + \
                                                    self.insulation.cond_cond[col_cond][(col_cond + 1) % 2]
                                                # TODO: only works for two windings
                                                col_cond = (col_cond + 1) % 2

                                            # Reset y
                                            winding_number = winding_numbers[col_cond]
                                            y = top_bound - windings[col_cond].conductor_radius
                                            top_window_iso_counter += 1

                                        else:
                                            break

                                    else:
                                        # is this winding already finished? - continue with the other one
                                        col_cond = (col_cond + 1) % 2

                                        # Correct the reset of y and correct x displacement
                                        x += windings[col_cond].conductor_radius - windings[(col_cond + 1) % 2].conductor_radius - \
                                            self.insulation.cond_cond[col_cond][(col_cond + 1) % 2] + self.insulation.cond_cond[col_cond][col_cond]

                                        y = top_bound - windings[col_cond].conductor_radius
                                        top_window_iso_counter -= 1
                            else:
                                # This will draw from bottom bound up.

                                # Initialize the list, that counts the already placed conductors
                                N_completed = [0, 0]

                                # Initialize the starting conductor
                                if turns0 >= turns1:
                                    col_cond = 0
                                else:
                                    col_cond = 1

                                # Get winding number of current winding
                                winding_number = winding_numbers[col_cond]

                                # Initialize the x and y coordinate
                                x = left_bound + windings[col_cond].conductor_radius
                                y = bot_bound + windings[col_cond].conductor_radius

                                # Continue placing as long as not all conductors have been placed
                                while (turns0 - N_completed[0] != 0) or \
                                        (turns1 - N_completed[1] != 0):
                                    if turns[col_cond] - N_completed[col_cond] != 0:
                                        # is this winding not already finished?
                                        if x < right_bound - windings[col_cond].conductor_radius:
                                            while y < top_bound - windings[col_cond].conductor_radius and \
                                                    N_completed[col_cond] < turns[col_cond]:
                                                self.p_conductor[winding_number].append([
                                                    x,
                                                    y,
                                                    0,
                                                    self.mesh_data.c_center_conductor[winding_number]])

                                                self.p_conductor[winding_number].append([
                                                    x - windings[col_cond].conductor_radius,
                                                    y,
                                                    0,
                                                    self.mesh_data.c_conductor[winding_number]])

                                                self.p_conductor[winding_number].append([
                                                    x,
                                                    y + windings[col_cond].conductor_radius,
                                                    0,
                                                    self.mesh_data.c_conductor[winding_number]])

                                                self.p_conductor[winding_number].append([
                                                    x + windings[col_cond].conductor_radius,
                                                    y,
                                                    0,
                                                    self.mesh_data.c_conductor[winding_number]])

                                                self.p_conductor[winding_number].append([
                                                    x,
                                                    y - windings[col_cond].conductor_radius,
                                                    0,
                                                    self.mesh_data.c_conductor[winding_number]])

                                                N_completed[col_cond] += 1

                                                # one from bot to top
                                                y += windings[col_cond].conductor_radius * 2 + self.insulation.cond_cond[col_cond][col_cond]

                                            x += windings[col_cond].conductor_radius + windings[(col_cond + 1) % 2].conductor_radius + \
                                                self.insulation.cond_cond[col_cond][(col_cond + 1) % 2]

                                            # Reset y
                                            col_cond = (col_cond + 1) % 2
                                            y = bot_bound + windings[col_cond].conductor_radius

                                        else:
                                            break

                                    else:
                                        # is this winding already finished? - continue with the other one
                                        col_cond = (col_cond + 1) % 2
                                        winding_number = winding_numbers[col_cond]

                                        # Correct the reset of y and correct x displacement
                                        x += windings[col_cond].conductor_radius - windings[(col_cond + 1) % 2].conductor_radius - \
                                            self.insulation.cond_cond[col_cond][(col_cond + 1) % 2] + self.insulation.cond_cond[col_cond][col_cond]

                                        y = bot_bound + windings[col_cond].conductor_radius

                        if virtual_winding_window.winding_scheme == InterleavedWindingScheme.VerticalStacked:
                            winding_number0, winding_number1 = winding_numbers
                            y = bot_bound + winding0.conductor_radius

                            # Initialization
                            x = left_bound + winding0.conductor_radius
                            i = 0

                            # First winding from bottom to top
                            if winding0.conductor_arrangement == ConductorArrangement.Square:

                                while y < top_bound - winding0.conductor_radius and \
                                        i < turns0:
                                    while x < right_bound - winding0.conductor_radius and \
                                            i < turns0:
                                        self.p_conductor[winding_number0].append([
                                            x,
                                            y,
                                            0,
                                            self.mesh_data.c_center_conductor[winding_number0]])
                                        self.p_conductor[winding_number0].append([
                                            x - winding0.conductor_radius,
                                            y,
                                            0,
                                            self.mesh_data.c_conductor[winding_number0]])
                                        self.p_conductor[winding_number0].append([
                                            x,
                                            y + winding0.conductor_radius,
                                            0,
                                            self.mesh_data.c_conductor[winding_number0]])
                                        self.p_conductor[winding_number0].append([
                                            x + winding0.conductor_radius,
                                            y,
                                            0,
                                            self.mesh_data.c_conductor[winding_number0]])
                                        self.p_conductor[winding_number0].append([
                                            x,
                                            y - winding0.conductor_radius,
                                            0,
                                            self.mesh_data.c_conductor[winding_number0]])
                                        i += 1
                                        x += winding0.conductor_radius * 2 + self.insulation.cond_cond[winding_number0][winding_number0]  # from left to right
                                    y += winding0.conductor_radius * 2 + self.insulation.cond_cond[winding_number0][winding_number0]  # one step from bot to top
                                    x = left_bound + winding0.conductor_radius  # always the same
                            elif winding0.conductor_arrangement == ConductorArrangement.Hexagonal:
                                base_line = True
                                while y < top_bound - winding0.conductor_radius and \
                                        i < turns0:
                                    while x < right_bound - winding0.conductor_radius and \
                                            i < turns0:
                                        self.p_conductor[winding_number0].append([
                                            x,
                                            y,
                                            0,
                                            self.mesh_data.c_center_conductor[winding_number0]])
                                        self.p_conductor[winding_number0].append([
                                            x - winding0.conductor_radius,
                                            y,
                                            0,
                                            self.mesh_data.c_conductor[winding_number0]])
                                        self.p_conductor[winding_number0].append([
                                            x,
                                            y + winding0.conductor_radius,
                                            0,
                                            self.mesh_data.c_conductor[winding_number0]])
                                        self.p_conductor[winding_number0].append([
                                            x + winding0.conductor_radius,
                                            y,
                                            0,
                                            self.mesh_data.c_conductor[winding_number0]])
                                        self.p_conductor[winding_number0].append([
                                            x,
                                            y - winding0.conductor_radius,
                                            0,
                                            self.mesh_data.c_conductor[winding_number0]])
                                        i += 1

                                        x += 2 * np.cos(np.pi / 6) * (winding0.conductor_radius + \
                                                                      self.insulation.cond_cond[winding_number0][winding_number0] / 2)

                                        # depending on what line, hexa scheme starts shifted
                                        # reset y to "new" bottom
                                        base_line = (not base_line)
                                        if base_line:
                                            y -= (winding0.conductor_radius + self.insulation.cond_cond[winding_number0][winding_number0] / 2)
                                        else:
                                            y += (winding0.conductor_radius + self.insulation.cond_cond[winding_number0][winding_number0] / 2)

                                    # Undo last base_line reset
                                    if base_line:
                                        y += (winding0.conductor_radius + self.insulation.cond_cond[winding_number0][winding_number0] / 2)
                                    else:
                                        y -= (winding0.conductor_radius + self.insulation.cond_cond[winding_number0][winding_number0] / 2)

                                    base_line = True
                                    x = left_bound + winding0.conductor_radius
                                    y += winding0.conductor_radius * 2 + self.insulation.cond_cond[winding_number0][winding_number0]
                            elif winding0.conductor_arrangement == ConductorArrangement.SquareFullWidth:
                                raise Exception("ConductorArrangement SquareFullWidth is not implemented for interleaved and vertical stacked")
                            else:
                                raise Exception(f"Unknown conductor_arrangement {winding0.conductor_arrangement}")

                            # Second winding from top to bottom

                            y = top_bound - winding1.conductor_radius
                            i = 0

                            if winding1.conductor_arrangement == ConductorArrangement.Square:
                                while y > bot_bound + winding1.conductor_radius and i < \
                                        turns1:
                                    while x < right_bound - winding1.conductor_radius and i < \
                                            turns1:
                                        self.p_conductor[winding_number1].append([
                                            x,
                                            y,
                                            0,
                                            self.mesh_data.c_center_conductor[winding_number1]])
                                        self.p_conductor[winding_number1].append([
                                            x - winding1.conductor_radius,
                                            y,
                                            0,
                                            self.mesh_data.c_conductor[winding_number1]])
                                        self.p_conductor[winding_number1].append([
                                            x,
                                            y + winding1.conductor_radius,
                                            0,
                                            self.mesh_data.c_conductor[winding_number1]])
                                        self.p_conductor[winding_number1].append([
                                            x + winding1.conductor_radius,
                                            y,
                                            0,
                                            self.mesh_data.c_conductor[winding_number1]])
                                        self.p_conductor[winding_number1].append([
                                            x,
                                            y - winding1.conductor_radius,
                                            0,
                                            self.mesh_data.c_conductor[winding_number1]])
                                        i += 1

                                        x += winding1.conductor_radius * 2 + self.insulation.cond_cond[winding_number1][winding_number1]  # from left to right
                                    # one step from bot to top
                                    y += -(winding1.conductor_radius * 2) - self.insulation.cond_cond[winding_number1][winding_number1]
                                    x = left_bound + winding1.conductor_radius  # always the same
                            elif winding1.conductor_arrangement == ConductorArrangement.Hexagonal:
                                base_line = True

                                while y > bot_bound + winding1.conductor_radius and \
                                        i < turns1:
                                    while x < right_bound - winding1.conductor_radius and \
                                            i < turns1:

                                        self.p_conductor[winding_number1].append([
                                            x,
                                            y,
                                            0,
                                            self.mesh_data.c_center_conductor[winding_number1]])
                                        self.p_conductor[winding_number1].append([
                                            x - winding1.conductor_radius,
                                            y,
                                            0,
                                            self.mesh_data.c_conductor[winding_number1]])
                                        self.p_conductor[winding_number1].append([
                                            x,
                                            y + winding1.conductor_radius,
                                            0,
                                            self.mesh_data.c_conductor[winding_number1]])
                                        self.p_conductor[winding_number1].append([
                                            x + winding1.conductor_radius,
                                            y,
                                            0,
                                            self.mesh_data.c_conductor[winding_number1]])
                                        self.p_conductor[winding_number1].append([
                                            x,
                                            y - winding1.conductor_radius,
                                            0,
                                            self.mesh_data.c_conductor[winding_number1]])

                                        i += 1
                                        x += 2 * np.cos(np.pi / 6) * \
                                            (winding1.conductor_radius + self.insulation.cond_cond[winding_number1][winding_number1] / 2)

                                        # depending on what line, hexa scheme starts shifted
                                        # reset y to "new" bottom
                                        base_line = (not base_line)
                                        if base_line:
                                            y += (winding1.conductor_radius + self.insulation.cond_cond[winding_number1][winding_number1])
                                        else:
                                            y -= (winding1.conductor_radius + self.insulation.cond_cond[winding_number1][winding_number1])

                                    # Undo last base_line reset
                                    if base_line:
                                        y -= (winding1.conductor_radius + self.insulation.cond_cond[winding_number1][winding_number1])
                                    else:
                                        y += (winding1.conductor_radius + self.insulation.cond_cond[winding_number1][winding_number1])

                                    base_line = True
                                    x = left_bound + winding1.conductor_radius
                                    # from top to bottom
                                    y += -(winding1.conductor_radius * 2) - self.insulation.cond_cond[winding_number1][winding_number1]
                            else:
                                raise Exception(f"Unknown conductor_arrangement {winding1.conductor_arrangement.name}")

                    elif virtual_winding_window.winding_type == WindingType.Single:
                        # One winding in the virtual winding window
                        winding = virtual_winding_window.windings[0]
                        turns = sum(virtual_winding_window.turns)
                        # TODO:  find another solution for this (turns = ...) (is needed in mesh.py for air_stacked) see set_winding in model
                        conductor_type = winding.conductor_type
                        winding_scheme = virtual_winding_window.winding_scheme
                        alignment = virtual_winding_window.alignment
                        placing_strategy = virtual_winding_window.placing_strategy
                        foil_vertical_placing_strategy = virtual_winding_window.foil_vertical_placing_strategy
                        foil_horizontal_placing_strategy = virtual_winding_window.foil_horizontal_placing_strategy

                        zigzag = virtual_winding_window.zigzag
                        num = winding.winding_number

                        # Check if the coil is round or rectangular
                        if conductor_type == ConductorType.RectangularSolid:
                            # Now check for each possible winding scheme
                            if winding_scheme == WindingScheme.Full:
                                # Full window conductor
                                self.p_conductor[num].append([
                                    left_bound,
                                    bot_bound,
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                self.p_conductor[num].append([
                                    right_bound,
                                    bot_bound,
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                self.p_conductor[num].append([
                                    left_bound,
                                    top_bound,
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                self.p_conductor[num].append([
                                    right_bound,
                                    top_bound,
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                center_point = self.get_center_of_rect([left_bound, bot_bound], [right_bound, bot_bound], [left_bound, top_bound],
                                                                       [right_bound, top_bound])
                                self.p_conductor[num].append([center_point[0], center_point[1], 0, self.mesh_data.c_center_conductor[num]])
                            elif winding_scheme == WindingScheme.FoilVertical:
                                # TODO Add check if turns do not fit in winding window
                                # Foil conductors where each conductor is very high and the conductors are expanding in the x-direction
                                if virtual_winding_window.wrap_para == WrapParaType.FixedThickness:
                                    # Wrap defined number of turns and chosen thickness
                                    winding.a_cell = winding.thickness * (top_bound - bot_bound)
                                    for i in range(turns):
                                        # Starting point of x depending on the distribution way if it is from left to right or vice versa.
                                        if foil_vertical_placing_strategy == FoilVerticalDistribution.HorizontalRightward:
                                            x_start = left_bound + i * winding.thickness + i * self.insulation.cond_cond[num][num]
                                            x_move = left_bound + (i + 1) * winding.thickness + i * self.insulation.cond_cond[num][num]
                                            if x_move > right_bound:
                                                break
                                        elif foil_vertical_placing_strategy == FoilVerticalDistribution.HorizontalLeftward:
                                            x_start = right_bound - i * winding.thickness - i * self.insulation.cond_cond[num][num]
                                            x_move = right_bound - (i + 1) * winding.thickness - i * self.insulation.cond_cond[num][num]
                                            if x_move < left_bound:
                                                break
                                        # Foil
                                        self.p_conductor[num].extend([
                                            [x_start, bot_bound, 0, self.mesh_data.c_conductor[num]],
                                            [x_move, bot_bound, 0, self.mesh_data.c_conductor[num]],
                                            [x_start, top_bound, 0, self.mesh_data.c_conductor[num]],
                                            [x_move, top_bound, 0, self.mesh_data.c_conductor[num]]
                                        ])
                                        # Find the center point of each turn
                                        center_point = self.get_center_of_rect(self.p_conductor[num][-4], self.p_conductor[num][-3],
                                                                               self.p_conductor[num][-2], self.p_conductor[num][-1])
                                        self.p_conductor[num].append([center_point[0], center_point[1], 0, self.mesh_data.c_center_conductor[num]])
                                # Interpolate type is where the foils will have a dynamic thickness
                                elif virtual_winding_window.wrap_para == WrapParaType.Interpolate:
                                    # Fill the allowed space in the Winding Window with a chosen number of turns
                                    # we need first to calculate the area of every turn
                                    # Find the wrap turn space
                                    winding.thickness = (right_bound - left_bound - (turns - 1) * self.insulation.cond_cond[num][num]) / turns
                                    window_height = top_bound - bot_bound
                                    winding.a_cell = winding.thickness * window_height
                                    # Update MeshData with the new thickness
                                    self.mesh_data.update_data(frequency=self.mesh_data.frequency, skin_mesh_factor=self.mesh_data.skin_mesh_factor)
                                    for i in range(turns):
                                        if foil_vertical_placing_strategy == FoilVerticalDistribution.HorizontalRightward:
                                            # Generate interpolated positions for each turn, starting from the left and moving right
                                            x_interpol = np.linspace(left_bound, right_bound + self.insulation.cond_cond[num][num], turns + 1)
                                            x_start = x_interpol[i]
                                            x_move = x_interpol[i + 1] - self.insulation.cond_cond[num][num]
                                        elif foil_vertical_placing_strategy == FoilVerticalDistribution.HorizontalLeftward:
                                            # Generate interpolated positions for each turn, starting from the right and moving left
                                            x_interpol = np.linspace(right_bound, left_bound - self.insulation.cond_cond[num][num], turns + 1)
                                            x_start = x_interpol[i]
                                            x_move = x_interpol[i + 1] + self.insulation.cond_cond[num][num]
                                        # Foil
                                        self.p_conductor[num].extend([
                                            [x_start, bot_bound, 0, self.mesh_data.c_conductor[num]],
                                            [x_move, bot_bound, 0, self.mesh_data.c_conductor[num]],
                                            [x_start, top_bound, 0, self.mesh_data.c_conductor[num]],
                                            [x_move, top_bound, 0, self.mesh_data.c_conductor[num]]
                                        ])
                                        # Find the center point of each turn
                                        center_point = self.get_center_of_rect(self.p_conductor[num][-4], self.p_conductor[num][-3],
                                                                               self.p_conductor[num][-2], self.p_conductor[num][-1])
                                        self.p_conductor[num].append([center_point[0], center_point[1], 0, self.mesh_data.c_center_conductor[num]])

                                else:
                                    raise Exception(f"Unknown wrap para type {virtual_winding_window.wrap_para}")
                                # Apply alignment
                                if alignment == Align.ToEdges:
                                    pass
                                if alignment == Align.CenterOnVerticalAxis:
                                    raise ValueError("FoilVertical Conductors can not be centered on vertical axis as they are very high")
                                if alignment == Align.CenterOnHorizontalAxis:
                                    min_x = min(point[0] for point in self.p_conductor[num])
                                    max_x = max(point[0] for point in self.p_conductor[num])
                                    occupied_width = max_x - min_x
                                    occupied_midpoint_x = min_x + (occupied_width / 2)
                                    window_midpoint_x = (left_bound + right_bound) / 2
                                    adjustment_x = window_midpoint_x - occupied_midpoint_x
                                    for i, _ in enumerate(self.p_conductor[num]):
                                        self.p_conductor[num][i][0] += adjustment_x
                            # Foil conductors where each conductor is very long and the conductors are expanding the y-direction
                            elif winding_scheme == WindingScheme.FoilHorizontal:
                                # the user can choose the thickness
                                if virtual_winding_window.wrap_para == WrapParaType.FixedThickness:
                                    # Find the turn space
                                    winding.a_cell = winding.thickness * (right_bound - left_bound)
                                    for i in range(turns):
                                        # Starting point of y depending on the distribution way if it is from left to right or vice versa.
                                        if foil_horizontal_placing_strategy == FoilHorizontalDistribution.VerticalUpward:
                                            y_start = bot_bound + i * winding.thickness + i * self.insulation.cond_cond[num][num]
                                            y_move = bot_bound + (i + 1) * winding.thickness + i * self.insulation.cond_cond[num][num]
                                            if round(y_move, 6) > round(top_bound, 6):
                                                break
                                        elif foil_horizontal_placing_strategy == FoilHorizontalDistribution.VerticalDownward:
                                            y_start = top_bound - i * winding.thickness - i * self.insulation.cond_cond[num][num]
                                            y_move = top_bound - (i + 1) * winding.thickness - i * self.insulation.cond_cond[num][num]
                                            if round(y_move) < round(bot_bound):
                                                break
                                        # Foil
                                        self.p_conductor[num].extend([
                                            [left_bound, y_start, 0, self.mesh_data.c_conductor[num]],
                                            [right_bound, y_start, 0, self.mesh_data.c_conductor[num]],
                                            [left_bound, y_move, 0, self.mesh_data.c_conductor[num]],
                                            [right_bound, y_move, 0, self.mesh_data.c_conductor[num]]
                                        ])
                                        # Find the center point of each turn
                                        center_point = self.get_center_of_rect(self.p_conductor[num][-4], self.p_conductor[num][-3],
                                                                               self.p_conductor[num][-2], self.p_conductor[num][-1])
                                        self.p_conductor[num].append([center_point[0], center_point[1], 0, self.mesh_data.c_center_conductor[num]])
                                # Interpolate type is where the foils will have a dynamic thickness
                                elif virtual_winding_window.wrap_para == WrapParaType.Interpolate:
                                    # Turn space
                                    winding.thickness = (top_bound - bot_bound - (turns - 1) * self.insulation.cond_cond[num][num]) / turns
                                    window_width = right_bound - left_bound
                                    winding.a_cell = winding.thickness * window_width
                                    # Update MeshData with the new thickness
                                    self.mesh_data.update_data(frequency=self.mesh_data.frequency, skin_mesh_factor=self.mesh_data.skin_mesh_factor)
                                    for i in range(turns):
                                        # Placing Foil horizontal rectangular conductors from bot to top
                                        if foil_horizontal_placing_strategy == FoilHorizontalDistribution.VerticalUpward:
                                            # Generate interpolated positions for each turn, starting from the bottom and moving top
                                            y_interpol = np.linspace(bot_bound, top_bound + self.insulation.cond_cond[num][num], turns + 1)
                                            y_start = y_interpol[i]
                                            y_move = y_interpol[i + 1] - self.insulation.cond_cond[num][num]
                                        elif foil_horizontal_placing_strategy == FoilHorizontalDistribution.VerticalDownward:
                                            y_interpol = np.linspace(top_bound, bot_bound - self.insulation.cond_cond[num][num], turns + 1)
                                            y_start = y_interpol[i]
                                            y_move = y_interpol[i + 1] + self.insulation.cond_cond[num][num]
                                        # Foil
                                        self.p_conductor[num].extend([
                                            [left_bound, y_start, 0, self.mesh_data.c_conductor[num]],
                                            [right_bound, y_start, 0, self.mesh_data.c_conductor[num]],
                                            [left_bound, y_move, 0, self.mesh_data.c_conductor[num]],
                                            [right_bound, y_move, 0, self.mesh_data.c_conductor[num]]
                                        ])
                                        # Find the center point of each turn
                                        center_point = self.get_center_of_rect(self.p_conductor[num][-4], self.p_conductor[num][-3],
                                                                               self.p_conductor[num][-2], self.p_conductor[num][-1])
                                        self.p_conductor[num].append([center_point[0], center_point[1], 0, self.mesh_data.c_center_conductor[num]])

                                # Apply alignment
                                if alignment == Align.ToEdges:
                                    pass
                                if alignment == Align.CenterOnVerticalAxis:
                                    min_y = min(point[1] for point in self.p_conductor[num])
                                    max_y = max(point[1] for point in self.p_conductor[num])
                                    occupied_height = max_y - min_y
                                    occupied_midpoint_y = min_y + occupied_height / 2
                                    window_midpoint_y = (bot_bound + top_bound) / 2
                                    adjustment_y = window_midpoint_y - occupied_midpoint_y
                                    for i, _ in enumerate(self.p_conductor[num]):
                                        self.p_conductor[num][i][1] += adjustment_y
                                if alignment == Align.CenterOnHorizontalAxis:
                                    raise ValueError("FoilHorizontal Conductors can not be centered on horizontal axis as they are long")

                            else:
                                raise Exception(f"Winding scheme {winding_scheme.name} is not implemented.")

                        elif conductor_type == ConductorType.RoundSolid or conductor_type == ConductorType.RoundLitz:
                            # Since round conductors have no winding scheme check for each conductor_arrangement.
                            conductor_arrangement = winding.conductor_arrangement
                            if conductor_arrangement == ConductorArrangement.Square:
                                # Define initial conditions based on the placement strategy.
                                # 16 cases will be handled here, 8 cases with consistent direction, and 8 cases with Zigzag movement.
                                vertical_first = placing_strategy in [
                                    ConductorDistribution.VerticalUpward_HorizontalRightward,
                                    ConductorDistribution.VerticalUpward_HorizontalLeftward,
                                    ConductorDistribution.VerticalDownward_HorizontalRightward,
                                    ConductorDistribution.VerticalDownward_HorizontalLeftward]

                                upward_movement = placing_strategy in [
                                    ConductorDistribution.VerticalUpward_HorizontalRightward,
                                    ConductorDistribution.VerticalUpward_HorizontalLeftward,
                                    ConductorDistribution.HorizontalRightward_VerticalUpward,
                                    ConductorDistribution.HorizontalLeftward_VerticalUpward]

                                rightward_movement = placing_strategy in [
                                    ConductorDistribution.VerticalUpward_HorizontalRightward,
                                    ConductorDistribution.VerticalDownward_HorizontalRightward,
                                    ConductorDistribution.HorizontalRightward_VerticalUpward,
                                    ConductorDistribution.HorizontalRightward_VerticalDownward]

                                # Set the starting position and step size based on initial conditions
                                if vertical_first:
                                    if upward_movement:
                                        start_y = bot_bound + winding.conductor_radius  # Start from the bottom
                                        step_y = winding.conductor_radius * 2 + self.insulation.cond_cond[num][num]
                                    else:
                                        start_y = top_bound - winding.conductor_radius  # Start from the top
                                        step_y = -(winding.conductor_radius * 2 + self.insulation.cond_cond[num][num])

                                    if rightward_movement:
                                        start_x = left_bound + winding.conductor_radius  # Moving right after completing a column
                                        step_x = winding.conductor_radius * 2 + self.insulation.cond_cond[num][num]
                                    else:
                                        start_x = right_bound - winding.conductor_radius  # Moving left after completing a column
                                        step_x = -(winding.conductor_radius * 2 + self.insulation.cond_cond[num][num])
                                # Determine if the first movement is horizontally (rightward or leftward)
                                else:
                                    if rightward_movement:
                                        start_x = left_bound + winding.conductor_radius  # start from the left
                                        step_x = winding.conductor_radius * 2 + self.insulation.cond_cond[num][num]
                                    else:
                                        start_x = right_bound - winding.conductor_radius  # Start from the right
                                        step_x = -(winding.conductor_radius * 2 + self.insulation.cond_cond[num][num])

                                    if upward_movement:
                                        start_y = bot_bound + winding.conductor_radius  # Moving top after completing a raw
                                        step_y = winding.conductor_radius * 2 + self.insulation.cond_cond[num][num]
                                    else:
                                        start_y = top_bound - winding.conductor_radius  # Moving bottom after completing a raw
                                        step_y = -(winding.conductor_radius * 2 + self.insulation.cond_cond[num][num])

                                # Loop and place conductors accordingly
                                x = start_x
                                y = start_y
                                i = 0
                                # Vertically movement
                                if vertical_first:
                                    while i < turns and left_bound + winding.conductor_radius <= x <= right_bound - winding.conductor_radius:
                                        while i < turns and bot_bound + winding.conductor_radius <= y <= top_bound - winding.conductor_radius:
                                            self.p_conductor[num].append(
                                                [x, y, 0, self.mesh_data.c_center_conductor[num]])
                                            self.p_conductor[num].append(
                                                [x - winding.conductor_radius, y, 0, self.mesh_data.c_conductor[num]])
                                            self.p_conductor[num].append(
                                                [x, y + winding.conductor_radius, 0, self.mesh_data.c_conductor[num]])
                                            self.p_conductor[num].append(
                                                [x + winding.conductor_radius, y, 0, self.mesh_data.c_conductor[num]])
                                            self.p_conductor[num].append(
                                                [x, y - winding.conductor_radius, 0, self.mesh_data.c_conductor[num]])
                                            y += step_y
                                            i += 1
                                        if not zigzag:
                                            # Start the next column with the same starting point (consistent direction)
                                            y = start_y
                                        else:
                                            # Alternating between top and bottom for the Zigzag movement
                                            step_y *= -1
                                            y += step_y
                                        # Moving one step horizontally (right or left)
                                        x += step_x
                                # Horizontally movement
                                else:
                                    while i < turns and bot_bound + winding.conductor_radius <= y <= top_bound - winding.conductor_radius:
                                        while i < turns and left_bound + winding.conductor_radius <= x <= right_bound - winding.conductor_radius:
                                            self.p_conductor[num].append(
                                                [x, y, 0, self.mesh_data.c_center_conductor[num]])
                                            self.p_conductor[num].append(
                                                [x - winding.conductor_radius, y, 0, self.mesh_data.c_conductor[num]])
                                            self.p_conductor[num].append(
                                                [x, y + winding.conductor_radius, 0, self.mesh_data.c_conductor[num]])
                                            self.p_conductor[num].append(
                                                [x + winding.conductor_radius, y, 0, self.mesh_data.c_conductor[num]])
                                            self.p_conductor[num].append(
                                                [x, y - winding.conductor_radius, 0, self.mesh_data.c_conductor[num]])
                                            x += step_x
                                            i += 1
                                        if not zigzag:
                                            # Start the next raw with the same starting point (consistent direction)
                                            x = start_x
                                        else:
                                            # Alternating between right and left for the Zigzag movement
                                            step_x *= -1
                                            x += step_x
                                        # Moving one step vertically (top or bottom)
                                        y += step_y

                                # Align to edges
                                if alignment == Align.ToEdges:
                                    # No adjustment of x or y
                                    pass

                                # Center the turns on horizontal axis
                                elif alignment == Align.CenterOnVerticalAxis:
                                    # Calculate vertical bounds of the occupied space by turns to adjust y
                                    min_y = min(point[1] for point in self.p_conductor[num])  # Find the lowest y position among all turns
                                    max_y = max(point[1] for point in self.p_conductor[num])  # Find the highest y position among all turns
                                    # The height of the space occupied by the turns is the difference between the highest and lowest y positions
                                    occupied_height = max_y - min_y
                                    # Calculate the vertical midpoint of the occupied space
                                    occupied_midpoint_y = min_y + occupied_height / 2
                                    # Calculate the midpoint of the winding window's vertical bounds as it is rectangle
                                    window_midpoint_y = (bot_bound + top_bound) / 2
                                    # Determine the adjustment needed for vertical centering
                                    adjustment_y = window_midpoint_y - occupied_midpoint_y
                                    # Apply the adjustment of y position to all points
                                    for i, _ in enumerate(self.p_conductor[num]):
                                        self.p_conductor[num][i][1] += adjustment_y

                                # Center the turns on vertical axis
                                elif alignment == Align.CenterOnHorizontalAxis:
                                    # Calculate horizontal bounds of the occupied space by turns to adjust x
                                    min_x = min(point[0] for point in self.p_conductor[num])  # Find the leftmost x position among all turns
                                    max_x = max(point[0] for point in self.p_conductor[num])  # Find the rightmost x position among all turns
                                    # The width of the space occupied by the turns is the difference between the rightmost and leftmost x positions
                                    occupied_width = max_x - min_x
                                    # Calculate the horizontal midpoint of the occupied space
                                    occupied_midpoint_x = min_x + (occupied_width / 2)
                                    # Calculate the midpoint of the winding window's horizontal bounds
                                    window_midpoint_x = (left_bound + right_bound) / 2
                                    # Determine the adjustment needed for horizontal centering
                                    adjustment_x = window_midpoint_x - occupied_midpoint_x
                                    # Apply the adjustment of x positions to all points
                                    for i, _ in enumerate(self.p_conductor[num]):
                                        self.p_conductor[num][i][0] += adjustment_x

                            elif conductor_arrangement == ConductorArrangement.Hexagonal:
                                y = bot_bound + winding.conductor_radius
                                x = left_bound + winding.conductor_radius
                                i = 0
                                base_line = True
                                # Case n_conductors higher that "allowed" is missing
                                while x < right_bound - winding.conductor_radius \
                                        and i < turns:
                                    while y < top_bound - winding.conductor_radius and \
                                            i < turns:
                                        self.p_conductor[num].append([
                                            x,
                                            y,
                                            0,
                                            self.mesh_data.c_center_conductor[num]])
                                        self.p_conductor[num].append([
                                            x - winding.conductor_radius,
                                            y,
                                            0,
                                            self.mesh_data.c_conductor[num]])
                                        self.p_conductor[num].append([
                                            x,
                                            y + winding.conductor_radius,
                                            0,
                                            self.mesh_data.c_conductor[num]])
                                        self.p_conductor[num].append([
                                            x + winding.conductor_radius,
                                            y,
                                            0,
                                            self.mesh_data.c_conductor[num]])
                                        self.p_conductor[num].append([
                                            x,
                                            y - winding.conductor_radius,
                                            0,
                                            self.mesh_data.c_conductor[num]])
                                        i += 1
                                        y += winding.conductor_radius * 2 + self.insulation.cond_cond[num][num]  # from bottom to top
                                    x += 2 * np.cos(np.pi / 6) * (winding.conductor_radius + self.insulation.cond_cond[num][num] / 2)
                                    # * np.sqrt(2 / 3 * np.pi / np.sqrt(3))  # one step from left to right
                                    # depending on what line, hexa scheme starts shifted
                                    # reset y to "new" bottom
                                    base_line = (not base_line)
                                    if base_line:
                                        y = bot_bound + winding.conductor_radius
                                    else:
                                        y = bot_bound + 2 * winding.conductor_radius + self.insulation.cond_cond[num][num] / 2

                            elif conductor_arrangement == ConductorArrangement.SquareFullWidth:
                                y = bot_bound + winding.conductor_radius
                                x = left_bound + winding.conductor_radius
                                i = 0
                                # Case n_conductors higher that "allowed" is missing
                                while round(y, 6) <= round(top_bound - winding.conductor_radius, 6) and i < turns:
                                    while x <= right_bound - winding.conductor_radius and i < turns:
                                        self.p_conductor[num].append([
                                            x,
                                            y,
                                            0,
                                            self.mesh_data.c_center_conductor[num]])
                                        self.p_conductor[num].append([
                                            x - winding.conductor_radius,
                                            y,
                                            0,
                                            self.mesh_data.c_conductor[num]])
                                        self.p_conductor[num].append([
                                            x,
                                            y + winding.conductor_radius,
                                            0,
                                            self.mesh_data.c_conductor[num]])
                                        self.p_conductor[num].append([
                                            x + winding.conductor_radius,
                                            y,
                                            0,
                                            self.mesh_data.c_conductor[num]])
                                        self.p_conductor[num].append([
                                            x,
                                            y - winding.conductor_radius,
                                            0,
                                            self.mesh_data.c_conductor[num]])
                                        i += 1
                                        # TODO: anisotrop insulation  # one step from left to right
                                        x += winding.conductor_radius * 2 + self.insulation.cond_cond[num][num]
                                    y += winding.conductor_radius * 2 + self.insulation.cond_cond[num][num]  # bot to top
                                    x = left_bound + winding.conductor_radius  # always the same

                            else:
                                raise Exception(f"Conductor arrangement {conductor_arrangement} is not implemented.")

                        else:
                            raise Exception(f"Conductor type {winding.conductor_type.name} is not implemented.")

                    elif virtual_winding_window.winding_type == WindingType.CenterTappedGroup:
                        # Primary Winding
                        winding = virtual_winding_window.windings[0]
                        turns = virtual_winding_window.turns[0]
                        num = 0

                        y = bot_bound + winding.conductor_radius
                        x = left_bound + winding.conductor_radius
                        i = 0
                        # Case n_conductors higher that "allowed" is missing
                        while y < top_bound - winding.conductor_radius and i < turns:
                            while x < right_bound - winding.conductor_radius and i < turns:
                                self.p_conductor[num].append([
                                    x,
                                    y,
                                    0,
                                    self.mesh_data.c_center_conductor[num]])
                                self.p_conductor[num].append([
                                    x - winding.conductor_radius,
                                    y,
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                self.p_conductor[num].append([
                                    x,
                                    y + winding.conductor_radius,
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                self.p_conductor[num].append([
                                    x + winding.conductor_radius,
                                    y,
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                self.p_conductor[num].append([
                                    x,
                                    y - winding.conductor_radius,
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                i += 1
                                x += winding.conductor_radius * 2 + self.insulation.cond_cond[num][num]  # from left to top
                            y += winding.conductor_radius * 2 + self.insulation.cond_cond[num][num]  # one step from left to right
                            x = left_bound + winding.conductor_radius  # always the same

                        # use y for next winding in stack
                        bot_bound = y - winding.conductor_radius - self.insulation.cond_cond[0][0] + self.insulation.cond_cond[0][1]

                        # Secondary Winding
                        winding = virtual_winding_window.windings[1]
                        turns = virtual_winding_window.turns[1]
                        num = 1
                        winding.a_cell = winding.thickness * (right_bound - left_bound)

                        for i in range(turns):
                            # CHECK if top bound is reached
                            # statement = (bot_bound + (i + 1) * winding.thickness + i * self.insulation.inner_winding_insulations[num]) <= top_bound
                            # print(statement)
                            if True:
                                low = bot_bound + i * winding.thickness + i * self.insulation.cond_cond[num][num]
                                high = bot_bound + (i + 1) * winding.thickness + i * self.insulation.cond_cond[num][num]
                                # stacking from the ground
                                self.p_conductor[num].append([
                                    left_bound,
                                    low,
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                self.p_conductor[num].append([
                                    right_bound,
                                    low,
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                self.p_conductor[num].append([
                                    left_bound,
                                    high,
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                self.p_conductor[num].append([
                                    right_bound,
                                    high,
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                center_point = self.get_center_of_rect(self.p_conductor[num][-4], self.p_conductor[num][-3],
                                                                       self.p_conductor[num][-2], self.p_conductor[num][-1])
                                self.p_conductor[num].append([center_point[0], center_point[1], 0, self.mesh_data.c_center_conductor[num]])

                        # use y for next winding in stack
                        bot_bound = high + self.insulation.cond_cond[1][1]

                        # Secondary Winding
                        winding = virtual_winding_window.windings[2]
                        turns = virtual_winding_window.turns[2]
                        num = 2
                        winding.a_cell = winding.thickness * (right_bound - left_bound)

                        for i in range(turns):
                            # CHECK if top bound is reached
                            # statement = (bot_bound + (i + 1) * winding.thickness + i * self.insulation.inner_winding_insulations[num]) <= top_bound
                            # print(statement)
                            if True:
                                low = bot_bound + i * winding.thickness + i * self.insulation.cond_cond[num][num]
                                high = bot_bound + (i + 1) * winding.thickness + i * self.insulation.cond_cond[num][num]
                                # stacking from the ground
                                self.p_conductor[num].append([
                                    left_bound,
                                    low,
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                self.p_conductor[num].append([
                                    right_bound,
                                    low,
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                self.p_conductor[num].append([
                                    left_bound,
                                    high,
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                self.p_conductor[num].append([
                                    right_bound,
                                    high,
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                center_point = self.get_center_of_rect(self.p_conductor[num][-4], self.p_conductor[num][-3],
                                                                       self.p_conductor[num][-2], self.p_conductor[num][-1])
                                self.p_conductor[num].append([center_point[0], center_point[1], 0, self.mesh_data.c_center_conductor[num]])

                    else:
                        raise Exception(f"Unknown winding type {virtual_winding_window.winding_type}")

    def check_number_of_drawn_conductors(self):
        """Check if number of conductors fits into the given winding window size."""
        needed_number_of_turns = 0
        all_windings = []

        for winding_window in self.winding_windows:
            for vww in winding_window.virtual_winding_windows:
                for _, winding in enumerate(vww.windings):
                    if winding not in all_windings:
                        all_windings.append(winding)

                needed_number_of_turns += sum(vww.turns)

        drawn_number_of_turns = 0

        for winding in all_windings:
            # Convert to numpy
            self.p_conductor[winding.winding_number] = np.asarray(self.p_conductor[winding.winding_number])

            if winding.conductor_type in [ConductorType.RoundSolid, ConductorType.RoundLitz]:
                drawn_number_of_turns += int(self.p_conductor[winding.winding_number].shape[0] / 5)  # round conductors
            else:
                drawn_number_of_turns += int(self.p_conductor[winding.winding_number].shape[0] / 5)  # rectangular conductors

        if drawn_number_of_turns != needed_number_of_turns:
            self.femmt_print(f"{drawn_number_of_turns=}")
            self.femmt_print(f"{needed_number_of_turns=}")
            raise Exception("Winding mismatch. Probably too many turns that do not fit in the winding window")

    def draw_region_single(self):
        """
        Region for Boundary Condition: Draws a rectangular outer region. Needed for air gaps in the outer flux path.

        Currently not used. Code kept for future implementations
        """
        self.p_region_bound[0][:] = [-self.r_outer * self.mesh_data.padding,
                                     -(self.core.window_h / 2 + self.core.core_thickness) * self.mesh_data.padding,
                                     0,
                                     self.mesh_data.c_core * self.mesh_data.padding]
        self.p_region_bound[1][:] = [self.r_outer * self.mesh_data.padding,
                                     -(self.core.window_h / 2 + self.core.core_thickness) * self.mesh_data.padding,
                                     0,
                                     self.mesh_data.c_core * self.mesh_data.padding]
        self.p_region_bound[2][:] = [-self.r_outer * self.mesh_data.padding,
                                     (self.core.window_h / 2 + self.core.core_thickness) * self.mesh_data.padding,
                                     0,
                                     self.mesh_data.c_core * self.mesh_data.padding]
        self.p_region_bound[3][:] = [self.r_outer * self.mesh_data.padding,
                                     (self.core.window_h / 2 + self.core.core_thickness) * self.mesh_data.padding,
                                     0,
                                     self.mesh_data.c_core * self.mesh_data.padding]

    def draw_insulations(self):
        """
        Draw insulations.

        Important:
        ---------
        Because the distance from the core to the winding is set by
        iso.core_cond, a delta, which is used in order for no overlapping lines will cause
        the "real" insulation to be slightly smaller than set by the user.
        """
        if not self.insulation.flag_insulation:
            return  # if flag_insulation is False, just return without drawing insulations

        if self.component_type == ComponentType.IntegratedTransformer:
            # if not self.stray_path:
            #
            #     # TODO: insulations implement for integrated_transformers
            #     # TODO Change back to warnings?
            #     self.femmt_print("Insulations are not set because they are not implemented for integrated transformers.")
            # else:
            iso = self.insulation
            mesh_data = self.mesh_data
            # Since there are many cases in which alternating conductors would lead to slightly different
            # mesh densities a simplification is made: Just use the lowest mesh density to be safe all the time.
            mesh_density_to_winding = min(mesh_data.c_conductor)

            mesh_density_to_core = mesh_data.c_window
            if self.insulation.max_aspect_ratio == 0:
                # If no aspect ratio is set insulations will not be drawn
                return
            else:
                insulation_delta = self.mesh_data.c_window / self.insulation.max_aspect_ratio

            # self.p_iso_core = []  # Order: Left, Top, Right, Bot
            self.p_iso_top_core = []
            self.p_iso_bot_core = []
            self.p_iso_pri_sec = []
            # Useful points for insulation creation
            left_x = self.core.core_inner_diameter / 2 + insulation_delta
            right_x = self.r_inner - insulation_delta
            if self.stray_path:
                # Useful dimensions needed for calculating some points
                window_h = self.core.window_h
                start_index = self.stray_path.start_index
                air_gap_1_position = self.air_gaps.midpoints[start_index][1]
                air_gap_1_height = self.air_gaps.midpoints[start_index][2]
                air_gap_2_position = self.air_gaps.midpoints[start_index + 1][1]
                air_gap_2_height = self.air_gaps.midpoints[start_index + 1][2]
                y_top_stray_path = air_gap_2_position - (air_gap_2_height / 2)
                y_bot_stray_path = air_gap_1_position + (air_gap_1_height / 2)
                top_section_top_y = window_h / 2 - insulation_delta
                top_section_bot_y = y_top_stray_path + insulation_delta
                bot_section_top_y = y_bot_stray_path - insulation_delta
                bot_section_bot_y = -window_h / 2 + insulation_delta
            else:
                window_h_top = self.core.window_h_top
                window_h_bot = self.core.window_h_bot
                top_section_top_y = window_h_bot / 2 + self.core.core_thickness + window_h_top - insulation_delta
                top_section_bot_y = window_h_bot / 2 + self.core.core_thickness + insulation_delta
                bot_section_top_y = window_h_bot / 2 - insulation_delta
                bot_section_bot_y = -window_h_bot / 2 + insulation_delta

            # Useful lengths
            # top core
            top_section_left_iso_width = iso.top_section_core_cond[2] - insulation_delta - insulation_delta
            top_section_top_iso_height = iso.top_section_core_cond[0] - insulation_delta - insulation_delta
            top_section_right_iso_width = iso.top_section_core_cond[3] - insulation_delta - insulation_delta
            top_section_bot_iso_height = iso.top_section_core_cond[1] - insulation_delta - insulation_delta
            # bot core
            bot_section_left_iso_width = iso.bot_section_core_cond[2] - insulation_delta - insulation_delta
            bot_section_top_iso_height = iso.bot_section_core_cond[0] - insulation_delta - insulation_delta
            bot_section_right_iso_width = iso.bot_section_core_cond[3] - insulation_delta - insulation_delta
            bot_section_bot_iso_height = iso.bot_section_core_cond[1] - insulation_delta - insulation_delta

            # Core to Pri insulation
            # top core points
            iso_top_core_left = [
                [
                    left_x,
                    top_section_top_y - top_section_top_iso_height - insulation_delta,
                    0,
                    mesh_density_to_core
                ],
                [
                    left_x + top_section_left_iso_width,
                    top_section_top_y - top_section_top_iso_height - insulation_delta,
                    0,
                    mesh_density_to_winding
                ],
                [
                    left_x + top_section_left_iso_width,
                    top_section_bot_y + top_section_bot_iso_height + insulation_delta,
                    0,
                    mesh_density_to_winding
                ],
                [
                    left_x,
                    top_section_bot_y + top_section_bot_iso_height + insulation_delta,
                    0,
                    mesh_density_to_core
                ]
            ]
            iso_top_core_top = [
                [
                    left_x,
                    top_section_top_y,
                    0,
                    mesh_density_to_core
                ],
                [
                    right_x,
                    top_section_top_y,
                    0,
                    mesh_density_to_core
                ],
                [
                    right_x,
                    top_section_top_y - top_section_top_iso_height,
                    0,
                    mesh_density_to_winding
                ],
                [
                    left_x,
                    top_section_top_y - top_section_top_iso_height,
                    0,
                    mesh_density_to_winding
                ]
            ]
            iso_top_core_right = [
                [
                    right_x - top_section_right_iso_width,
                    top_section_top_y - top_section_top_iso_height - insulation_delta,
                    0,
                    mesh_density_to_winding
                ],
                [
                    right_x,
                    top_section_top_y - top_section_top_iso_height - insulation_delta,
                    0,
                    mesh_density_to_core
                ],
                [
                    right_x,
                    top_section_bot_y + top_section_bot_iso_height + insulation_delta,
                    0,
                    mesh_density_to_core
                ],
                [
                    right_x - top_section_right_iso_width,
                    top_section_bot_y + top_section_bot_iso_height + insulation_delta,
                    0,
                    mesh_density_to_winding
                ]
            ]
            iso_top_core_bot = [
                [
                    left_x,
                    top_section_bot_y + top_section_bot_iso_height,
                    0,
                    mesh_density_to_winding
                ],
                [
                    right_x,
                    top_section_bot_y + top_section_bot_iso_height,
                    0,
                    mesh_density_to_winding
                ],
                [
                    right_x,
                    top_section_bot_y,
                    0,
                    mesh_density_to_core
                ],
                [
                    left_x,
                    top_section_bot_y,
                    0,
                    mesh_density_to_core
                ]
            ]
            # bot core points
            iso_bot_core_left = [
                [
                    left_x,
                    bot_section_top_y - bot_section_top_iso_height - insulation_delta,
                    0,
                    mesh_density_to_core
                ],
                [
                    left_x + bot_section_left_iso_width,
                    bot_section_top_y - bot_section_top_iso_height - insulation_delta,
                    0,
                    mesh_density_to_winding
                ],
                [
                    left_x + bot_section_left_iso_width,
                    bot_section_bot_y + bot_section_bot_iso_height + insulation_delta,
                    0,
                    mesh_density_to_winding
                ],
                [
                    left_x,
                    bot_section_bot_y + bot_section_bot_iso_height + insulation_delta,
                    0,
                    mesh_density_to_core
                ]
            ]
            iso_bot_core_top = [
                [
                    left_x,
                    bot_section_top_y,
                    0,
                    mesh_density_to_core
                ],
                [
                    right_x,
                    bot_section_top_y,
                    0,
                    mesh_density_to_core
                ],
                [
                    right_x,
                    bot_section_top_y - bot_section_top_iso_height,
                    0,
                    mesh_density_to_winding
                ],
                [
                    left_x,
                    bot_section_top_y - bot_section_top_iso_height,
                    0,
                    mesh_density_to_winding
                ]
            ]
            iso_bot_core_right = [
                [
                    right_x - bot_section_right_iso_width,
                    bot_section_top_y - bot_section_top_iso_height - insulation_delta,
                    0,
                    mesh_density_to_winding
                ],
                [
                    right_x,
                    bot_section_top_y - bot_section_top_iso_height - insulation_delta,
                    0,
                    mesh_density_to_core
                ],
                [
                    right_x,
                    bot_section_bot_y + bot_section_bot_iso_height + insulation_delta,
                    0,
                    mesh_density_to_core
                ],
                [
                    right_x - bot_section_right_iso_width,
                    bot_section_bot_y + bot_section_bot_iso_height + insulation_delta,
                    0,
                    mesh_density_to_winding
                ]
            ]
            iso_bot_core_bot = [
                [
                    left_x,
                    bot_section_bot_y + bot_section_bot_iso_height,
                    0,
                    mesh_density_to_winding
                ],
                [
                    right_x,
                    bot_section_bot_y + bot_section_bot_iso_height,
                    0,
                    mesh_density_to_winding
                ],
                [
                    right_x,
                    bot_section_bot_y,
                    0,
                    mesh_density_to_core
                ],
                [
                    left_x,
                    bot_section_bot_y,
                    0,
                    mesh_density_to_core
                ]
            ]
            self.p_iso_top_core = [iso_top_core_left, iso_top_core_top, iso_top_core_right, iso_top_core_bot]
            self.p_iso_bot_core = [iso_bot_core_left, iso_bot_core_top, iso_bot_core_right, iso_bot_core_bot]

        else:
            window_h = self.core.window_h
            iso = self.insulation
            mesh_data = self.mesh_data

            # Since there are many cases in which alternating conductors would lead to slightly different
            # mesh densities a simplification is made: Just use the lowest mesh density to be safe all the time.
            mesh_density_to_winding = min(mesh_data.c_conductor)

            mesh_density_to_core = mesh_data.c_window

            # Using the delta the lines and points from the insulation and the core/windings are not overlapping
            # which makes creating the mesh simpler
            # Insulation between winding and core
            # Since an aspect ratio is given the insulation delta is calculated using the length of the longest side of the triangle,
            # which is always smaller than c_window.
            if self.insulation.max_aspect_ratio == 0:
                # If no aspect ratio is set insulations will not be drawn
                return
            else:
                insulation_delta = self.mesh_data.c_window / self.insulation.max_aspect_ratio

            self.p_iso_core = []  # Order: Left, Top, Right, Bot
            self.p_iso_pri_sec = []

            # Useful points for insulation creation
            left_x = self.core.core_inner_diameter / 2 + insulation_delta
            top_y = window_h / 2 - insulation_delta
            right_x = self.r_inner - insulation_delta
            bot_y = -window_h / 2 + insulation_delta

            # Useful lengths
            left_iso_width = iso.core_cond[2] - insulation_delta - insulation_delta
            top_iso_height = iso.core_cond[0] - insulation_delta - insulation_delta
            right_iso_width = iso.core_cond[3] - insulation_delta - insulation_delta
            bot_iso_height = iso.core_cond[1] - insulation_delta - insulation_delta

            # Core to Pri insulation
            iso_core_left = [
                [
                    left_x,
                    top_y - top_iso_height - insulation_delta,
                    0,
                    mesh_density_to_core
                ],
                [
                    left_x + left_iso_width,
                    top_y - top_iso_height - insulation_delta,
                    0,
                    mesh_density_to_winding
                ],
                [
                    left_x + left_iso_width,
                    bot_y + bot_iso_height + insulation_delta,
                    0,
                    mesh_density_to_winding
                ],
                [
                    left_x,
                    bot_y + bot_iso_height + insulation_delta,
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
                    top_y - top_iso_height - insulation_delta,
                    0,
                    mesh_density_to_winding
                ],
                [
                    right_x,
                    top_y - top_iso_height - insulation_delta,
                    0,
                    mesh_density_to_core
                ],
                [
                    right_x,
                    bot_y + bot_iso_height + insulation_delta,
                    0,
                    mesh_density_to_core
                ],
                [
                    right_x - right_iso_width,
                    bot_y + bot_iso_height + insulation_delta,
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

            """
            Currently there are only core to vww insulations
            # Now create insulations for interleaved windings
            self.p_iso_interleaved = []
            for vww in self.virtual_winding_windows:
                if vww.winding_type == WindingType.Interleaved:
                    if vww.winding_scheme == InterleavedWindingScheme.Bifilar:
                        print("No insulations for bifilar winding scheme")
                    elif vww.winding_scheme == InterleavedWindingScheme.HorizontalAlternating:
                        winding_0 = vww.windings[0]
                        winding_1 = vww.windings[1]
                        turns0 = vww.turns[0]
                        turns1 = vww.turns[1]
                        current_x = vww.left_bound + 2 * winding_0.conductor_radius

                        mesh_density_left_side = 0
                        mesh_density_right_side = 0
                        if turns0 > turns1:
                            mesh_density_left_side = self.mesh_data.c_conductor[0]
                            mesh_density_right_side = self.mesh_data.c_conductor[1]
                        else:
                            mesh_density_left_side = self.mesh_data.c_conductor[1]
                            mesh_density_right_side = self.mesh_data.c_conductor[0]

                        current_iso = []
                        for index in range(self.top_window_iso_counter-1):
                            current_iso.append([
                                [
                                    current_x + insulation_delta,
                                    top_y - top_iso_height - insulation_delta,
                                    0,
                                    mesh_density_left_side
                                ],
                                [
                                    current_x + iso.vww_insulation - insulation_delta,
                                    top_y - top_iso_height - insulation_delta,
                                    0,
                                    mesh_density_right_side
                                ],
                                [
                                    current_x + iso.vww_insulation - insulation_delta,
                                    bot_y + bot_iso_height + insulation_delta,
                                    0,
                                    mesh_density_right_side
                                ],
                                [
                                    current_x + insulation_delta,
                                    bot_y + bot_iso_height + insulation_delta,
                                    0,
                                    mesh_density_left_side
                                ]
                            ])
                            # The sec winding can start first when it has a higher number of turns
                            # col_cond
                            if index % 2 == self.col_cond_start:
                                current_x += iso.vww_insulation + 2 * winding_1.conductor_radius
                            else:
                                current_x += iso.vww_insulation + 2 * winding_0.conductor_radius
                    elif vww.winding_scheme == InterleavedWindingScheme.VerticalAlternating:
                        print("No insulations for vertical alternating winding scheme")
                    elif vww.winding_scheme == InterleavedWindingScheme.VerticalStacked:
                        print("No insulations for vertical stacked winding scheme")
                    else:
                        raise Exception(f"Unknown interleaved winding scheme {vww.winding_scheme.name}")
                else:
                    # TODO Own handler for warnings?
                    print("No insulations for winding type {vww.winding_type.name}")
            """

    def draw_model(self):
        """Draw the full component.

        The full component consists of core, air gaps, insulation and the conductors.
        The insulation is only drawn if self.insulation.flag_insulation was set.
        """
        self.draw_outer()
        if self.core.core_type == CoreType.Single:
            self.draw_single_window()
            self.draw_air_gaps()
        if self.core.core_type == CoreType.Stacked:
            self.draw_stacked_windows()
        self.draw_conductors()
        self.check_number_of_drawn_conductors()

        if self.insulation.flag_insulation:  # check flag before drawing insulations
            self.draw_insulations()

        # TODO: Region
        # if self.core.core_type == CoreType.Single:
        #     self.draw_region_single()
        # if self.core.core_type == CoreType.Stacked:
        #     self.draw_region_stacked()

    @staticmethod
    def get_center_of_rect(p1: List[float], p2: List[float], p3: List[float], p4: List[float]):
        """
        Get center point of a rectangular conductor.

        :param p1: Point 1 as a x,y-List
        :type p1: List[float]
        :param p2: Point 1 as a x,y-List
        :type p2: List[float]
        :param p3: Point 1 as a x,y-List
        :type p3: List[float]
        :param p4: Point 1 as a x,y-List
        :type p4: List[float]
        """
        x_list = [p1[0], p2[0], p3[0], p4[0]]
        y_list = [p1[1], p2[1], p3[1], p4[1]]

        return np.mean(x_list), np.mean(y_list)
