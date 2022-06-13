from enum import Enum
from typing import List

class AirGapMethod(Enum):
    Center = "center"
    Percent = "percent"
    Manually = "manually"

class AirGapLegPosition(Enum):
    LeftLeg = -1
    CenterLeg = 0
    RightLeg = 1

class AirGaps():
    method: AirGapMethod
    leg_positions: List[AirGapLegPosition] = []
    positions: List[float] = []
    heights: List[float] = []

    def __init__(self, method: AirGapMethod):
        self.method = method

    def add_air_gap(self, leg_position: AirGapLegPosition, position_value: float, height: float):
        self.leg_positions.append(leg_position)
        self.positions.append(position_value)
        self.heights.append(height)

    def get_all_parameters(self):
        position_tag = []
        air_gap_position = []
        air_gap_h = []

        for leg_position, position, height in zip(self.leg_positions, self.positions, self.heights):
            position_tag.append(leg_position.value)
            air_gap_position.append(position)
            air_gap_h.append(height)

        return {
            "method": self.method.value,
            "n_air_gaps": len(self.leg_positions),
            "position_tag": position_tag,
            "air_gap_position": air_gap_position,
            "air_gap_h": air_gap_h 
        }

    