"""Data transfer objects (DTOs) used by this toolbox."""
from dataclasses import dataclass
from femmt.enumerations import WindingTag
from typing import Optional, List


@dataclass
class SingleCoreDimensions:
    """Defines the dimensions of a default core."""

    core_inner_diameter: float
    window_w: float
    window_h: float
    core_h: float


@dataclass
class StackedCoreDimensions:
    """Defines the dimensions of a stacked core. A stacked core is made of a transformer section and an inductor."""

    core_inner_diameter: float
    window_w: float
    window_h_top: float
    window_h_bot: float


@dataclass
class ConductorRow:
    """Defines the conductors in one row."""

    number_of_conds_per_winding: int
    number_of_conds_per_row: Optional[int]
    row_height: Optional[float]
    winding_tag: WindingTag
    number_of_rows: Optional[int]
    additional_bobbin: float


@dataclass
class ThreeWindingIsolation:
    """Defines the insulation of a three-winding transformer."""

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
    """Definitions for the center tapped group. A group is made of several primary and secondary rows."""

    primary_number_of_rows: Optional[int]
    secondary_number_of_rows: Optional[int]
    primary_rest: Optional[int]
    secondary_rest: Optional[int]
    stack: List[WindingTag]


@dataclass
class ConductorStack:
    """Definitions for the conductor stack."""

    number_of_groups: int
    number_of_single_rows: Optional[int]
    order: List[int]


@dataclass
class StackIsolation:
    """Definition for the stack insulation."""

    thickness: float


@dataclass
class WireMaterial:
    """Definitions for the wire material."""

    name: str
    sigma: float
    temperature: float
    temperature_coefficient: float
    thermal_conductivity: float
    volumetric_mass_density: float


@dataclass
class TransformerInductance:
    """Inductance definitions for a two-winding transformer."""

    l_h_conc: float
    l_s_conc: float
    n_conc: float
    M: float
    L_1_1: float
    L_2_2: float


@dataclass
class ThreeWindingTransformerInductance:
    """Inductance definitions for a three-winding transformer."""

    M_12: Optional[float]
    M_13: Optional[float]
    M_23: Optional[float]
    L_s1: Optional[float]
    L_s2: Optional[float]
    L_s3: Optional[float]
    L_h: Optional[float]
    n_12: Optional[float]
    n_13: Optional[float]
    n_23: Optional[float]
    L_s12: Optional[float]
    L_s13: Optional[float]
    L_s23: Optional[float]
