from dataclasses import dataclass
from femmt.enumerations import *


@dataclass
class SingleCoreDimensions:
    core_inner_diameter: float
    window_w: float
    window_h: float
    core_h: float


@dataclass
class StackedCoreDimensions:
    core_inner_diameter: float
    window_w: float
    window_h_top: float
    window_h_bot: float


@dataclass
class ConductorRow:
    number_of_conds_per_winding: int
    number_of_conds_per_row: int
    row_height: float
    winding_tag: WindingTag
    number_of_rows: int
    additional_bobbin: float


@dataclass
class ThreeWindingIsolation:
    primary_to_primary: float
    primary_to_secondary: float
    primary_to_tertiary: float
    secondary_to_primary: float
    secondary_to_secondary: float
    secondary_to_tertiary: float
    tertiary_to_primary: float
    tertiary_to_secondary: float
    tertiary_to_tertiary: float


@dataclass
class CenterTappedGroup:
    primary_number_of_rows: int
    secondary_number_of_rows: int
    primary_rest: int
    secondary_rest: int
    stack: [WindingTag]


@dataclass
class ConductorStack:
    number_of_groups: int
    number_of_single_rows: int
    order: [int]


@dataclass
class StackIsolation:
    thickness: float


@dataclass
class WireMaterial:
    name: str
    sigma: float
    temperature: float
    temperature_coefficient: float
    thermal_conductivity: float
    volumetric_mass_density: float

@dataclass
class TransformerInductance:
    l_h_conc: float
    l_s_conc: float
    n_conc: float
    M: float
    L_1_1: float
    L_2_2: float

@dataclass
class ThreeWindingTransformerInductance:
    M_12: float
    M_13: float
    M_23: float
    L_s1: float
    L_s2: float
    L_s3: float
    L_h: float
    n_12: float
    n_13: float
    n_23: float
    L_s12: float
    L_s13: float
    L_s23: float