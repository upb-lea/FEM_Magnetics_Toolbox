"""DTOs for the stacked transformer optimization."""
# python libraries
from dataclasses import dataclass
from typing import List, Optional

# 3rd party libraries
import numpy as np
from materialdatabase.dtos import MaterialCurve
from femmt.enumerations import *

@dataclass
class WorkingDirectories:
    """Working directories for an integrated transformer optimization."""

    fem_working_directory: str
    reluctance_model_results_directory: str
    fem_simulation_results_directory: str
    fem_simulation_filtered_results_directory: str
    fem_thermal_simulation_results_directory: str
    fem_thermal_filtered_simulation_results_directory: str

@dataclass
class StoTargetAndFixedParameters:
    """
    Stacked-transformer optimization target and fixed parameters.

    These parameters are calculated from the stacked-transformer input configuration (StoSingleInputConfig).
    """

    i_rms_1: float
    i_rms_2: float
    i_peak_1: float
    i_peak_2: float
    i_phase_deg_1: float
    i_phase_deg_2: float
    material_dto_curve_list: List[MaterialCurve]
    time_extracted_vec: List
    current_extracted_1_vec: List
    current_extracted_2_vec: List
    fundamental_frequency: float
    target_inductance_matrix: np.ndarray
    working_directories: WorkingDirectories

@dataclass
class StoInsulation:
    """Insulation definition for stacked transformer optimization."""

    # insulation for top core window
    iso_window_top_core_top: float
    iso_window_top_core_bot: float
    iso_window_top_core_left: float
    iso_window_top_core_right: float
    # insulation for bottom core window
    iso_window_bot_core_top: float
    iso_window_bot_core_bot: float
    iso_window_bot_core_left: float
    iso_window_bot_core_right: float
    # winding-to-winding insulation
    iso_primary_to_primary: float
    iso_secondary_to_secondary: float
    iso_primary_to_secondary: float

@dataclass
class StackedTransformerMaterialDataSources:
    """Data sources for the FEM simulation."""

    permeability_datasource: MaterialDataSource
    permeability_datatype: MeasurementDataType
    permeability_measurement_setup: MeasurementSetup
    permittivity_datasource: MaterialDataSource
    permittivity_datatype: MeasurementDataType
    permittivity_measurement_setup: MeasurementSetup

@dataclass
class StoSingleInputConfig:
    """
    Configuration to simulate a stacked transformer.

    Input parameters are the target parameters, current vectors and the parameters to sweep.
    Also specifies the working directory where to store the results.
    """

    stacked_transformer_study_name: str

    # target parameters
    l_s12_target: float
    l_h_target: float
    n_target: float

    # operating point: current waveforms and temperature
    time_current_1_vec: np.ndarray
    time_current_2_vec: np.ndarray
    temperature: float

    # sweep parameters: geometry and materials
    n_p_top_min_max_list: list
    n_p_bot_min_max_list: list
    material_list: list
    core_name_list: Optional[List[str]]
    core_inner_diameter_min_max_list: Optional[List[float]]
    window_w_min_max_list: Optional[List[float]]
    window_h_bot_min_max_list: list
    primary_litz_wire_list: list
    secondary_litz_wire_list: list

    # maximum limitation for transformer total height and core volume
    max_transformer_total_height: float
    max_core_volume: float

    # fix parameters: insulations
    insulations: StoInsulation

    # misc
    stacked_transformer_optimization_directory: str
    fft_filter_value_factor: float
    mesh_accuracy: float

    # data sources
    material_data_sources: StackedTransformerMaterialDataSources

@dataclass
class FemInput:
    """Input DTO for a FEM simulation within the stacked transformer optimization."""

    # general parameters
    working_directory: str
    simulation_name: str

    # material and geometry parameters
    material_name: str
    primary_litz_wire_name: str
    secondary_litz_wire_name: str
    core_inner_diameter: float
    window_w: float
    window_h_top: float
    window_h_bot: float
    air_gap_length_top: float
    air_gap_length_bot: float
    turns_primary_top: float
    turns_primary_bot: float
    turns_secondary_bot: float
    insulations: StoInsulation

    # data sources
    material_data_sources: StackedTransformerMaterialDataSources

    # operating point conditions
    temperature: float
    fundamental_frequency: float
    time_current_1_vec: np.array
    time_current_2_vec: np.array

@dataclass
class FemOutput:
    """Output DTO for a FEM simulation within the stacked transformer optimization."""

    n_conc: float
    l_s_conc: float
    l_h_conc: float
    p_loss_winding_1: float
    p_loss_winding_2: float
    eddy_core: float
    core: float
