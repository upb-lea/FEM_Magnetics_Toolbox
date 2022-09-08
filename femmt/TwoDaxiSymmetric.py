# Python standard libraries
import numpy as np
from typing import List

# Local libraries
from Enumerations import *
from Data import MeshData
from Model import Core, VirtualWindingWindow, AirGaps, StrayPath, Isolation

class TwoDaxiSymmetric:
    """
    This class creates the model of the magnetic component. This is done by creating lists of points which
    will be later added to the gmsh file.

    :return:

    """

    core: Core
    virtual_winding_windows: List[VirtualWindingWindow]
    air_gaps: AirGaps
    stray_path: StrayPath
    isolation: Isolation
    component_type: ComponentType
    mesh_data: MeshData
    number_of_windings: int

    # List of points which represent the model
    # Every List is a List of 4 Points: x, y, z, mesh_factor
    p_outer: List[List[float]]
    p_region_bound: List[List[float]] 
    p_window: List[List[float]] 
    p_air_gaps: List[List[float]]
    p_conductor: List[List[float]]
    p_iso_core: List[List[float]]
    p_iso_pri_sec: List[List[float]]

    def __init__(self, core: Core, mesh_data: MeshData, air_gaps: AirGaps, virtual_winding_windows: List[VirtualWindingWindow], 
                stray_path: StrayPath, isolation: Isolation, component_type: ComponentType, number_of_windings: int):
        self.core = core
        self.mesh_data = mesh_data
        self.virtual_winding_windows = virtual_winding_windows
        self.air_gaps = air_gaps
        self.component_type = component_type
        self.stray_path = stray_path
        self.isolation = isolation
        self.number_of_windings = number_of_windings

        # -- Arrays for geometry data -- 
        # TODO Is the zero initialization necessary?
        self.p_outer = np.zeros((4, 4))
        self.p_region_bound = np.zeros((4, 4))
        self.p_window = np.zeros((4 * core.number_core_windows, 4))
        self.p_air_gaps = np.zeros((4 * air_gaps.number, 4)) 
        self.p_conductor = []
        self.p_iso_core = []
        self.p_iso_pri_sec = []

        for i in range(number_of_windings):
            self.p_conductor.insert(i, [])

        self.r_inner = core.r_inner
        self.r_outer = core.r_outer

    def draw_outer(self):
        """
        Draws the outer points

        :return:
        """
        # Outer Core
        # (A_zyl=2pi*r*h => h=0.5r=0.25core_w <=> ensure A_zyl=A_core on the tiniest point)
        # TODO Case core_h is not None
        if self.core.core_h is None:
            self.p_outer[0][:] = [-self.r_outer,
                                -(self.core.window_h / 2 + self.core.core_w / 4),
                                0,
                                self.mesh_data.c_core]

            self.p_outer[1][:] = [self.r_outer,
                                -(self.core.window_h / 2 + self.core.core_w / 4),
                                0,
                                self.mesh_data.c_core]

            self.p_outer[2][:] = [-self.r_outer,
                                (self.core.window_h / 2 + self.core.core_w / 4),
                                0,
                                self.mesh_data.c_core]

            self.p_outer[3][:] = [self.r_outer,
                                (self.core.window_h / 2 + self.core.core_w / 4),
                                0,
                                self.mesh_data.c_core]
        else:
            self.p_outer[0][:] = [-self.r_outer,
                                -self.core.core_h/2,
                                0,
                                self.mesh_data.c_core]

            self.p_outer[1][:] = [self.r_outer,
                                -self.core.core_h/2,
                                0,
                                self.mesh_data.c_core]

            self.p_outer[2][:] = [-self.r_outer,
                                self.core.core_h/2,
                                0,
                                self.mesh_data.c_core]

            self.p_outer[3][:] = [self.r_outer,
                                self.core.core_h/2,
                                0,
                                self.mesh_data.c_core]

    def draw_window(self):
        # Window
        # At this point both windows (in a cut) are modeled
        # print(f"win: c_window: {self.component.mesh.c_window}")
        self.p_window[0] = [-self.r_inner,
                            -self.core.window_h / 2,
                            0,
                            self.mesh_data.c_window]

        self.p_window[1] = [-self.core.core_w / 2,
                            -self.core.window_h / 2,
                            0,
                            self.mesh_data.c_window]

        self.p_window[2] = [-self.r_inner,
                            self.core.window_h / 2,
                            0,
                            self.mesh_data.c_window]

        self.p_window[3] = [-self.core.core_w / 2,
                            self.core.window_h / 2,
                            0,
                            self.mesh_data.c_window]

        self.p_window[4] = [self.core.core_w / 2,
                            -self.core.window_h / 2,
                            0,
                            self.mesh_data.c_window]

        self.p_window[5] = [self.r_inner,
                            -self.core.window_h / 2,
                            0,
                            self.mesh_data.c_window]

        self.p_window[6] = [self.core.core_w / 2,
                            self.core.window_h / 2,
                            0,
                            self.mesh_data.c_window]

        self.p_window[7] = [self.r_inner,
                            self.core.window_h / 2,
                            0,
                            self.mesh_data.c_window]

    def draw_air_gaps(self):
        # Air gaps
        # "air_gaps" is a list with [position_tag, air_gap_position, air_gap_h, c_air_gap]
        #   - position_tag: specifies the gapped "leg"
        #   - air_gap_position: specifies the coordinate of the air gap's center point along the specified leg
        #   - air_gap_h: height/length of the air gap
        #   - c_air_gap: mesh accuracy factor
        # at this point the 4 corner points of each air gap are generated out of "air_gaps"

        mesh_accuracy = self.core.window_w / 20 * self.mesh_data.global_accuracy

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
                air_gap_length_top = self.core.core_w / 2
                air_gap_length_bot = self.core.core_w / 2

                # Check for stray_paths in integrated transformers
                if self.component_type == ComponentType.IntegratedTransformer:
                    if self.stray_path.start_index == i:
                        # Stray path is above current air_gap
                        air_gap_length_top = self.stray_path.length
                    elif self.stray_path.start_index + 1 == i:
                        # Stray path is below current air_gap
                        air_gap_length_bot = self.stray_path.length

                # Bottom left
                self.p_air_gaps[i * 4 + 0] = [0,
                                                air_gap_y_position -
                                                air_gap_height / 2,
                                                0,
                                                self.mesh_data.c_core]

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
                                                self.mesh_data.c_core]

                # Top right
                self.p_air_gaps[i * 4 + 3] = [air_gap_length_top,
                                                air_gap_y_position +
                                                air_gap_height / 2,
                                                0,
                                                mesh_accuracy]

        # In order to close the air gap when a stray_path is added, additional points need to be added
        if self.component_type == ComponentType.IntegratedTransformer:
            top_point = [self.core.core_w / 2,
                            self.air_gaps.midpoints[self.stray_path.start_index+1][1] -
                            air_gap_height / 2,
                            0,
                            mesh_accuracy]
            bot_point = [self.core.core_w / 2,
                            self.air_gaps.midpoints[self.stray_path.start_index][1] +
                            air_gap_height / 2,
                            0,
                            mesh_accuracy]
            self.p_close_air_gaps = [top_point, bot_point]

    def draw_conductors(self, draw_top_down = True):
        # Draw every conductor type based on the virtual_winding_window bounds

        for virtual_winding_window in self.virtual_winding_windows:
            # Get bounds from virtual winding window
            bot_bound = virtual_winding_window.bot_bound
            top_bound = virtual_winding_window.top_bound
            left_bound = virtual_winding_window.left_bound
            right_bound = virtual_winding_window.right_bound

            # Check the possible WindingTypes and draw accordingly
            if virtual_winding_window.winding_type == WindingType.Interleaved:
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
                        while (turns0 - N_completed[0] != 0) or \
                                (turns1 - N_completed[1] != 0):
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

                                        y -= windings[col_cond].conductor_radius * 2 + self.isolation.cond_cond[col_cond]  # one from bot to top

                                    x += windings[col_cond].conductor_radius + \
                                            windings[(col_cond + 1) % 2].conductor_radius + \
                                            self.isolation.cond_cond[2]  # from left to right

                                    # Reset y
                                    col_cond = (col_cond + 1) % 2
                                    winding_number = winding_numbers[col_cond]
                                    y = top_bound - windings[col_cond].conductor_radius
                                    top_window_iso_counter += 1
                                else:
                                    break

                            else:
                                # is this winding already finished? - continue with the other one
                                col_cond = (col_cond + 1) % 2

                                # Correct the reset of y and correct x displacement
                                x += windings[col_cond].conductor_radius - \
                                        windings[(col_cond + 1) % 2].conductor_radius \
                                        - self.isolation.cond_cond[2] + self.isolation.cond_cond[
                                            col_cond]

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

                                        y += windings[col_cond].conductor_radius * 2 + \
                                                self.isolation.cond_cond[col_cond]  # one from bot to top

                                    x += windings[col_cond].conductor_radius + \
                                            windings[(col_cond + 1) % 2].conductor_radius + \
                                            self.isolation.cond_cond[2]  # from left to right

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
                                x += windings[col_cond].conductor_radius - \
                                        windings[(col_cond + 1) % 2].conductor_radius \
                                        - self.isolation.cond_cond[2] + \
                                        self.isolation.cond_cond[
                                            col_cond]

                                y = bot_bound + windings[col_cond].conductor_radius
                
                if virtual_winding_window.winding_scheme == InterleavedWindingScheme.VerticalStacked:
                    winding_number0, winding_number1 = winding_numbers

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
                                self.p_conductor[num].append([
                                    x, 
                                    y - winding0.conductor_radius, 
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                i += 1
                                x += winding0.conductor_radius * 2 + \
                                        self.isolation.cond_cond[
                                            num]  # from left to right
                            y += winding0.conductor_radius * 2 + \
                                    self.isolation.cond_cond[
                                        num]  # one step from bot to top
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

                                    x += 2 * np.cos(np.pi / 6) * (
                                            winding0.conductor_radius +
                                            self.isolation.cond_cond[num] / 2)

                                    # depending on what line, hexa scheme starts shifted
                                    # reset y to "new" bottom
                                    base_line = (not base_line)
                                    if base_line:
                                        y -= (winding0.conductor_radius +
                                                self.isolation.cond_cond[num])
                                    else:
                                        y += (winding0.conductor_radius +
                                                self.isolation.cond_cond[num])

                                # Undo last base_line reset
                                if base_line:
                                    y += (winding0.conductor_radius +
                                            self.isolation.cond_cond[num])
                                else:
                                    y -= (winding0.conductor_radius +
                                            self.isolation.cond_cond[num])

                                base_line = True
                                x = left_bound + winding0.conductor_radius
                                y += winding0.conductor_radius + \
                                        self.component.isolation.cond_cond[num]
                    else:
                        raise Exception(f"Unknown conductor_arrangement {winding0.conductor_arrangement}")

                    # Second winding from top to bottom
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

                                x += winding1.conductor_radius * 2 + \
                                        self.isolation.cond_cond[
                                            num]  # from left to right
                            y += -(winding1.conductor_radius * 2) - \
                                    self.isolation.cond_cond[
                                        num]  # one step from bot to top
                            x = left_bound + winding1.conductor_radius  # always the same
                    elif winding1.conductor_arrangement == ConductorArrangement.Hexagonal:
                        base_line = True

                        while y > bot_bound + winding1.conductor_radius and \
                                i < turns1:
                            while x < right_bound - winding1.conductor_radius and \
                                    i < turns1:
                                print(f"i: {i} "
                                        f"x: {x} "
                                        f"y: {y} ")

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
                                x += 2 * np.cos(np.pi / 6) * (
                                        winding1.conductor_radius +
                                        self.isolation.cond_cond[num] / 2)

                                # depending on what line, hexa scheme starts shifted
                                # reset y to "new" bottom
                                base_line = (not base_line)
                                if base_line:
                                    y += (winding1.conductor_radius +
                                            self.isolation.cond_cond[num])
                                else:
                                    y -= (winding1.conductor_radius +
                                            self.isolation.cond_cond[num])

                            # Undo last base_line reset
                            if base_line:
                                y -= (winding1.conductor_radius +
                                        self.isolation.cond_cond[num])
                            else:
                                y += (winding1.conductor_radius +
                                        self.isolation.cond_cond[num])

                            base_line = True
                            x = left_bound + winding1.conductor_radius
                            # from top to bottom
                            y -= (winding1.conductor_radius +
                                    self.isolation.cond_cond[num])
                    else:
                        raise Exception(f"Unknown conductor_arrangement {winding1.conductor_arrangement}")
            elif virtual_winding_window.winding_type == WindingType.Single:
                # One winding in the virtual winding window
                winding = virtual_winding_window.windings[0]
                turns = virtual_winding_window.turns[0]
                conductor_type = winding.conductor_type
                winding_scheme = virtual_winding_window.winding_scheme

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
                    elif winding_scheme == WindingScheme.FoilVertical:
                        # TODO Add check if turns do not fit in winding window
                        # Foil conductors where each conductor is very high and the conductors are expanding in the x-direction
                        if virtual_winding_window.wrap_para == WrapParaType.FixedThickness:
                            # Wrap defined number of turns and chosen thickness
                            for i in range(turns):
                                # CHECK if right bound is reached
                                if (left_bound + (i + 1) * winding.thickness +
                                    i * self.isolation.cond_cond[num]) <= right_bound:
                                    # Foils
                                    self.p_conductor[num].append([
                                        left_bound + i * winding.thickness + i * self.isolation.cond_cond[num],
                                        bot_bound, 
                                        0, 
                                        self.mesh_data.c_conductor[num]])
                                    self.p_conductor[num].append([
                                        left_bound + (i + 1) * winding.thickness + i * self.isolation.cond_cond[num], 
                                        bot_bound, 
                                        0,
                                        self.mesh_data.c_conductor[num]])
                                    self.p_conductor[num].append([
                                        left_bound + i * winding.thickness + i * self.isolation.cond_cond[num],
                                        top_bound, 
                                        0, 
                                        self.mesh_data.c_conductor[num]])
                                    self.p_conductor[num].append([
                                        left_bound + (i + 1) * winding.thickness + i * self.isolation.cond_cond[num], 
                                        top_bound, 
                                        0,
                                        self.mesh_data.c_conductor[num]])
                        elif virtual_winding_window.wrap_para == WrapParaType.Interpolate:
                            # Fill the allowed space in the Winding Window with a chosen number of turns
                            x_interpol = np.linspace(left_bound, right_bound + self.isolation.cond_cond[num],
                                                        winding.turns + 1)
                            for i in range(turns):
                                # Foils
                                self.p_conductor[num].append([
                                    x_interpol[i], 
                                    bot_bound, 
                                    0, 
                                    self.component.mesh.c_conductor[num]])
                                self.p_conductor[num].append([
                                    x_interpol[i + 1] - self.component.isolation.cond_cond[num], 
                                    bot_bound, 
                                    0,
                                    self.component.mesh.c_conductor[num]])
                                self.p_conductor[num].append([
                                    x_interpol[i],
                                    top_bound, 
                                    0, 
                                    self.component.mesh.c_conductor[num]])
                                self.p_conductor[num].append([
                                    x_interpol[i + 1] - self.component.isolation.cond_cond[num], 
                                    top_bound, 
                                    0,
                                    self.component.mesh.c_conductor[num]])
                        else:
                            raise Exception(f"Unknown wrap para type {virtual_winding_window.wrap_para}")
                    elif winding_scheme == WindingScheme.FoilHorizontal:
                        # Foil conductors where each conductor is very long and the conductors are expanding the y-direction
                        # Stack defined number of turns and chosen thickness
                        for i in range(turns):
                            # CHECK if top bound is reached
                            if (bot_bound + (i + 1) * winding.thickness +
                                i * self.isolation.cond_cond[num]) <= top_bound:
                                # stacking from the ground
                                self.p_conductor[num].append([
                                    left_bound, 
                                    bot_bound + i * winding.thickness + i * self.isolation.cond_cond[num], 
                                    0, 
                                    self.mesh_data.c_conductor[num]])
                                self.p_conductor[num].append([
                                    right_bound,
                                    bot_bound + i * winding.thickness + i * self.isolation.cond_cond[num],
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                self.p_conductor[num].append([
                                    left_bound,
                                    bot_bound + (i + 1) * winding.thickness + i * self.isolation.cond_cond[num], 
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                self.p_conductor[num].append([
                                    right_bound,
                                    bot_bound + (i + 1) * winding.thickness + i * self.isolation.cond_cond[num],
                                    0,
                                    self.mesh_data.c_conductor[num]])      
                    else:
                        raise Exception(f"Winding scheme {winding_scheme} is not implemented.")
                elif conductor_type == ConductorType.RoundSolid or conductor_type == ConductorType.RoundLitz:
                    # Since round conductors have no winding scheme check for each conductor_arrangement
                    conductor_arrangement = winding.conductor_arrangement

                    if conductor_arrangement == ConductorArrangement.Square:
                        y = bot_bound + winding.conductor_radius
                        x = left_bound + winding.conductor_radius
                        i = 0
                        # Case n_conductors higher that "allowed" is missing
                        while x < right_bound - winding.conductor_radius and i < turns:
                            while y < top_bound - winding.conductor_radius and i < turns:
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
                                self.p_conductor[num].append(
                                    [x + winding.conductor_radius, 
                                    y, 
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                self.p_conductor[num].append([
                                    x, 
                                    y - winding.conductor_radius, 
                                    0,
                                    self.mesh_data.c_conductor[num]])
                                i += 1
                                y += winding.conductor_radius * 2 + self.isolation.cond_cond[num]  # one step from left to right
                            x += winding.conductor_radius * 2 + self.isolation.cond_cond[num]  # from left to top
                            y = bot_bound + winding.conductor_radius
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
                                y += winding.conductor_radius * 2 + \
                                        self.isolation.cond_cond[
                                            num]  # from bottom to top
                            x += 2 * np.cos(np.pi / 6) * (
                                    winding.conductor_radius +
                                    self.isolation.cond_cond[num] / 2)
                            # * np.sqrt(2 / 3 * np.pi / np.sqrt(3))  # one step from left to right
                            # depending on what line, hexa scheme starts shifted
                            # reset y to "new" bottom
                            base_line = (not base_line)
                            if base_line:
                                y = bot_bound + winding.conductor_radius
                            else:
                                y = bot_bound + 2 * winding.conductor_radius + \
                                    self.isolation.cond_cond[
                                        num] / 2
                    elif conductor_arrangement == ConductorArrangement.SquareFullWidth:
                        y = bot_bound + winding.conductor_radius
                        x = left_bound + winding.conductor_radius
                        i = 0
                        # Case n_conductors higher that "allowed" is missing
                        while y < top_bound - winding.conductor_radius \
                                and i < turns:
                            while x < right_bound - winding.conductor_radius \
                                    and i < turns:
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
                                x += winding.conductor_radius * 2 + \
                                        self.isolation.cond_cond[
                                            num]  # from left to top
                            y += winding.conductor_radius * 2 + \
                                    self.isolation.cond_cond[
                                        num]  # one step from left to right
                            x = left_bound + winding.conductor_radius  # always the same
                    else:
                        raise Exception(f"Conductor arrangement {conductor_arrangement} is not implemented.")
                else:
                    raise Exception(f"Conductor shape {winding.conductor_shape} is not implemented.")
            else:
                raise Exception(f"Unknown winding type {virtual_winding_window.winding_type}")

        # Checking the Conductors
        for index, winding in enumerate(virtual_winding_window.windings):
            num = winding.winding_number

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
            if winding.conductor_type in [ConductorType.RoundSolid, ConductorType.RoundLitz]:
                if int(self.p_conductor[num].shape[0] / 5) < virtual_winding_window.turns[index]:
                    # Warning: warnings.warn("Too many turns that do not fit in the winding window.")
                    # Correct: self.component.windings[num].turns = int(self.p_conductor[num].shape[0]/5)
                    # TODO: break, but remove warning. valid bit should be set to False
                    #  Code must go to the next parameter-iteration step for geometric sweep
                    self.valid = False
                    # TODO Tell the user which winding window
                    raise Exception(f"Too many turns that do not fit in the winding window {str(virtual_winding_window)}")

            # Region for Boundary Condition
            self.p_region_bound[0][:] = [-self.r_outer * self.mesh_data.padding,
                                        -(self.core.window_h / 2 + self.core.core_w / 4)
                                        * self.mesh_data.padding,
                                        0,
                                        self.mesh_data.c_core * self.mesh_data.padding]
            self.p_region_bound[1][:] = [self.r_outer * self.mesh_data.padding,
                                        -(self.core.window_h / 2 + self.core.core_w / 4)
                                        * self.mesh_data.padding,
                                        0,
                                        self.mesh_data.c_core * self.mesh_data.padding]
            self.p_region_bound[2][:] = [-self.r_outer * self.mesh_data.padding,
                                        (self.core.window_h / 2 + self.core.core_w / 4)
                                        * self.mesh_data.padding,
                                        0,
                                        self.mesh_data.c_core * self.mesh_data.padding]
            self.p_region_bound[3][:] = [self.r_outer * self.mesh_data.padding,
                                        (self.core.window_h / 2 + self.core.core_w / 4)
                                        * self.mesh_data.padding,
                                        0,
                                        self.mesh_data.c_core * self.mesh_data.padding]

    def draw_isolations(self):
        """
        IMPORTANT
        Because the distance from the core to the winding is set by
        iso.core_cond, a delta, which is used in order for no overlapping lines will cause
        the "real" isolation to be slightly smaller than set by the user.
        """

        window_h = self.core.window_h
        iso = self.isolation
        mesh_data = self.mesh_data
        isolation_delta = self.isolation.isolation_delta

        # Since there are many cases in which alternating conductors would lead to slightly different
        # mesh densities a simplification is made: Just use the lowest mesh density to be safe all the time.
        mesh_density_to_winding = min(mesh_data.c_conductor)

        mesh_density_to_core = mesh_data.c_window

        # Using the delta the lines and points from the isolation and the core/windings are not overlapping
        # which makes creating the mesh more simpler
        # Isolation between winding and core
        iso_core_delta_left = isolation_delta["core_left"]
        iso_core_delta_top = isolation_delta["core_top"]
        iso_core_delta_right = isolation_delta["core_right"]
        iso_core_delta_bot = isolation_delta["core_bot"]
        iso_iso_delta = isolation_delta["iso_iso"]
        iso_winding_delta_left = isolation_delta["winding_left"]
        iso_winding_delta_top = isolation_delta["winding_top"]
        iso_winding_delta_right = isolation_delta["winding_right"]
        iso_winding_delta_bot = isolation_delta["winding_bot"]

        self.p_iso_core = [] # Order: Left, Top, Right, Bot
        self.p_iso_pri_sec = []

        if self.component_type == ComponentType.IntegratedTransformer:
            # TODO implement for integrated_transformers
            # TODO Change back to warnings?
            print("Isolations are not set because they are not implemented for integrated transformers.")
        else:
            # Useful points for isolation creation
            left_x = self.core.core_w / 2 + iso_core_delta_left
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

            """
            Since the way virtual winding windows are handled is changed this will be changed too

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
            """

    def draw_model(self):

        self.draw_outer()

        self.draw_window()

        self.draw_air_gaps()

        self.draw_conductors()

        self.draw_isolations()