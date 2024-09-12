"""Data transfer objects (DTOs) for the inductor optimization."""
# python libraries
import dataclasses
from typing import List, Optional

# 3rd party libraries

# own libraries
from materialdatabase.dtos import MaterialCurve
from femmt.enumerations import *
from femmt.optimization.ito_dtos import WorkingDirectories

@dataclasses.dataclass
class InductorMaterialDataSources:
    """Data sources for the FEM simulation."""

    permeability_datasource: MaterialDataSource
    permeability_datatype: MeasurementDataType
    permeability_measurement_setup: MeasurementSetup
    permittivity_datasource: MaterialDataSource
    permittivity_datatype: MeasurementDataType
    permittivity_measurement_setup: MeasurementSetup

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

    # general parameters
    inductor_study_name: str
    inductor_optimization_directory: str

    # target parameter
    target_inductance: float

    # fixed parameters
    insulations: InductorInsulationDTO
    temperature: float
    time_current_vec: List

    # optimization parameters
    material_name_list: List[str]
    core_name_list: Optional[List[str]]
    core_inner_diameter_list: Optional[List[float]]
    window_w_list: Optional[List[float]]
    window_h_list: Optional[List[float]]
    litz_wire_list: List[str]

    # FEM simulation
    material_data_sources: InductorMaterialDataSources

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
    fft_frequency_list: List[float]
    fft_amplitude_list: List[float]
    fft_phases_list: List[float]

@dataclasses.dataclass
class FemInput:
    """Input DTO for a FEM simulation within the inductor optimization."""

    # general parameters
    working_directory: str
    simulation_name: str

    # material and geometry parameters
    material_name: str
    litz_wire_name: str
    core_inner_diameter: float
    window_w: float
    window_h: float
    air_gap_length: float
    turns: float
    insulations: InductorInsulationDTO

    # data sources
    material_data_sources: InductorMaterialDataSources

    # operating point conditions
    temperature: float
    fundamental_frequency: float
    fft_frequency_list: List[float]
    fft_amplitude_list: List[float]
    fft_phases_list: List[float]

@dataclasses.dataclass
class FemOutput:
    """Output DTO for a FEM simulation within the inductor optimization."""

    fem_inductance: float
    fem_p_loss_winding: float
    fem_eddy_core: float
    fem_core: float
