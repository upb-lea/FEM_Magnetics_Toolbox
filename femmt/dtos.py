"""Data transfer objects (DTOs) used by this toolbox."""
from dataclasses import dataclass
from femmt.enumerations import WindingTag


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
    number_of_conds_per_row: int | None
    row_height: float | None
    winding_tag: WindingTag
    number_of_rows: int | None
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

    primary_number_of_rows: int | None
    secondary_number_of_rows: int | None
    primary_rest: int | None
    secondary_rest: int | None
    stack: list[WindingTag]


@dataclass
class ConductorStack:
    """Definitions for the conductor stack."""

    number_of_groups: int
    number_of_single_rows: int | None
    order: list[int]


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

    M_12: float | None
    M_13: float | None
    M_23: float | None
    L_s1: float | None
    L_s2: float | None
    L_s3: float | None
    L_h: float | None
    n_12: float | None
    n_13: float | None
    n_23: float | None
    L_s12: float | None
    L_s13: float | None
    L_s23: float | None
