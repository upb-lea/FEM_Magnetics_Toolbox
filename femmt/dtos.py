from dataclasses import dataclass


@dataclass
class SingleCoreDimensions:
    core_inner_diameter: float
    window_w: float
    window_h: float


@dataclass
class StackedCoreDimensions:
    core_inner_diameter: float
    window_w: float
    window_h_top: float
    window_h_bot: float
