"""Data transfer objects (DTOs) for the inductor optimization."""
# python libraries
import dataclasses
from typing import List, Optional

# 3rd party libraries

# own libraries
from materialdatabase.dtos import MaterialCurve

@dataclasses.dataclass
class WorkingDirectories:
    """Working directories."""

    pass

@dataclasses.dataclass
class InductorInsulationDTO:
    """Insulation distances for the inductor geometry."""

    core_left: float
    core_right: float
    core_top: float
    core_bot: float
    primary_to_primary: float

@dataclasses.dataclass
class InductorOptimizationDTO:
    """Contains boundary parameters for the inductor optimization."""

    material_name_list: List[str]
    core_name_list: Optional[List[str]]
    core_inner_diameter_list: Optional[List[float]]
    window_w_list: Optional[List[float]]
    window_h_list: Optional[List[float]]
    litz_wire_list: List[str]
    insulations: InductorInsulationDTO
    target_inductance: float
    temperature: float
    time_current_vec: List
    working_directory: str

@dataclasses.dataclass
class InductorOptimizationTargetAndFixedParameters:
    """
    Inductor optimization target and fixed parameters.

    These parameters are calculated from the inductor input configuration (InductorOptimizationDTO).
    """

    i_rms: float
    i_peak: float
    material_dto_curve_list: List[MaterialCurve]
    time_extracted_vec: List
    current_extracted_vec: List
    fundamental_frequency: float
    working_directories: WorkingDirectories
