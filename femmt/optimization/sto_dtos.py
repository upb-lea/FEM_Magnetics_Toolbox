"""DTOs for the stacked transformer optimization."""
# python libraries
from dataclasses import dataclass
from typing import List, Optional

# 3rd party libraries
import numpy as np
from materialdatabase.dtos import MaterialCurve
from femmt.enumerations import *
from magnethub.loss import LossModel

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
    magnet_hub_model_list: List[LossModel]
    time_extracted_vec: List
    current_extracted_1_vec: List
    current_extracted_2_vec: List
    fundamental_frequency: float
    target_inductance_matrix: np.ndarray
    working_directories: WorkingDirectories
    # winding 1
    fft_frequency_list_1: List[float]
    fft_amplitude_list_1: List[float]
    fft_phases_list_1: List[float]

    # winding 2
    fft_frequency_list_2: List[float]
    fft_amplitude_list_2: List[float]
    fft_phases_list_2: List[float]

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
    volume: float

@dataclass
class ReluctanceModelInput:
    """Input DTO for reluctance model simulation within the inductor optimization."""

    target_inductance_matrix: np.array
    core_inner_diameter: float
    window_w: float
    window_h_bot: float
    window_h_top: float
    turns_1_top: int
    turns_1_bot: int
    turns_2_bot: int
    litz_wire_name_1: str
    litz_wire_diameter_1: float
    litz_wire_name_2: str
    litz_wire_diameter_2: float

    insulations: StoInsulation
    material_dto: MaterialCurve
    magnet_material_model: LossModel

    temperature: float
    time_extracted_vec: List
    current_extracted_vec_1: List
    current_extracted_vec_2: List
    fundamental_frequency: float

    i_rms_1: float
    i_rms_2: float

    primary_litz_dict: dict
    secondary_litz_dict: dict

    # # winding 1
    fft_frequency_list_1: List[float]
    fft_amplitude_list_1: List[float]
    fft_phases_list_1: List[float]
    #
    # # winding 2
    fft_frequency_list_2: List[float]
    fft_amplitude_list_2: List[float]
    fft_phases_list_2: List[float]

@dataclass
class ReluctanceModelOutput:
    """output DTO for reluctance model simulation within the inductor optimization."""

    # set additional attributes
    p_hyst: float
    p_hyst_top: float
    p_hyst_bot: float
    p_hyst_middle: float
    b_max_top: float
    b_max_bot: float
    b_max_middle: float
    winding_1_loss: float
    winding_2_loss: float
    l_top_air_gap: float
    l_bot_air_gap: float
    volume: float
    area_to_heat_sink: float
    p_loss: float
