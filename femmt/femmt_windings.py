from enum import Enum
from typing import List
#from femmt_model_classes import *


class VirtualWindingWindow:
    """
    A virtual winding window is the area, where either some kind of interleaved conductors or a one winding
    (primary, secondary,...) is placed in a certain way.

    An instance of this class will be automatically created when the Winding is added to the MagneticComponent
    """

    # Rectangular frame:
    bot_bound: float
    top_bound: float
    left_bound: float
    right_bound: float

    def __init__(self, bot_bound: float, top_bound: float, left_bound: float, right_bound: float):
        self.bot_bound = bot_bound
        self.top_bound = top_bound
        self.left_bound = left_bound
        self.right_bound = right_bound

    def __repr__(self):
        return f"bot: {self.bot_bound}, top: {self.top_bound}, left: {self.left_bound}, right: {self.right_bound}"

class WindingWindow():
    max_bot_bound: float
    max_top_bound: float
    max_left_bound: float
    max_right_bound: float

    virtual_winding_windows: List[VirtualWindingWindow]

    def __init__(self, core, isolation, horizontal_split_factor: float, vertical_split_factor: float):
        
        #self.max_bot_bound = -core.window_h / 2 + isolation.core_cond[0]
        #self.max_top_bound = core.window_h / 2 - isolation.core_cond[1]
        #self.max_left_bound = core.core_w / 2 + isolation.core_cond[2]
        #self.max_right_bound = core.r_inner - isolation.core_cond[3]

        
        self.max_bot_bound = 0
        self.max_top_bound = 10
        self.max_left_bound = 0
        self.max_right_bound = 10

        # Isolation between vwws
        # isolation_vww = isolation.cond_cond[2]
        isolation_vww = 1

        # Splits
        horizontal_split = self.max_left_bound + (self.max_right_bound - self.max_left_bound) * horizontal_split_factor
        vertical_split = self.max_top_bound - abs(self.max_bot_bound - self.max_top_bound) * vertical_split_factor

        # Create a grid for the virtual winding windows
        top_left_vww = VirtualWindingWindow(bot_bound = horizontal_split + isolation_vww / 2,
                                            top_bound = self.max_top_bound,
                                            left_bound = self.max_left_bound,
                                            right_bound = vertical_split - isolation_vww / 2)

        top_right_vww = VirtualWindingWindow(bot_bound = horizontal_split + isolation_vww / 2,
                                            top_bound = self.max_top_bound,
                                            left_bound = vertical_split + isolation_vww / 2,
                                            right_bound = self.max_right_bound)

        bot_left_vww = VirtualWindingWindow(bot_bound = self.max_bot_bound,
                                            top_bound = horizontal_split - isolation_vww / 2,
                                            left_bound = self.max_left_bound,
                                            right_bound = vertical_split - isolation_vww / 2)

        bot_right_vww = VirtualWindingWindow(bot_bound = self.max_bot_bound,
                                            top_bound = horizontal_split - isolation_vww / 2,
                                            left_bound = vertical_split + isolation_vww / 2,
                                            right_bound = self.max_right_bound)

        self.virtual_winding_windows = [top_left_vww, top_right_vww, bot_left_vww, bot_right_vww]

    def combine_vww(self, vww1, vww2):
        index1 = self.virtual_winding_windows.index(vww1)
        index2 = self.virtual_winding_windows.index(vww2)

        if index2-index1 == 3:
            raise Exception("Cannot combine top left and bottom right.")

        self.virtual_winding_windows.remove(vww1)
        self.virtual_winding_windows.remove(vww2)

        new_vww = VirtualWindingWindow(bot_bound = min(vww1.bot_bound, vww2.bot_bound), 
                                    top_bound = max(vww1.top_bound, vww2.top_bound), 
                                    left_bound = min(vww1.left_bound, vww2.left_bound), 
                                    right_bound = max(vww1.right_bound, vww2.right_bound))

        self.virtual_winding_windows.append(new_vww)

        return new_vww


if __name__ == "__main__":
    winding_window = WindingWindow(None, None, 0.5, 0.5)
    top_left, top_right, bot_left, bot_right = winding_window.virtual_winding_windows


    winding_window.combine_vww(top_left, top_right)

    for window in winding_window.virtual_winding_windows:
        print(window)
