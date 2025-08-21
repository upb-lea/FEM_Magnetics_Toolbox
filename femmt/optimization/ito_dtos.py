"""DTOs for the integrated transformer optimization."""
# python libraries
from dataclasses import dataclass

# 3rd party libraries
import numpy as np
from materialdatabase.meta.data_classes import MaterialCurve
from materialdatabase.meta.data_enums import DataSource
from magnethub.loss import LossModel
from femmt.enumerations import *

@dataclass
class ItoInsulation:
    """Insulation definition for the integrated transformer optimization."""

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
class IntegratedTransformerMaterialDataSources:
    """Data sources for the FEM simulation."""

    permeability_datasource: MaterialDataSource
    permeability_datatype: DataSource
    permittivity_datasource: MaterialDataSource
    permittivity_datatype: DataSource

@dataclass
class ItoSingleInputConfig:
    """
    Configuration to simulate an integrated transformer.

    Input parameters are the target parameters, current vectors and the parameters to sweep.
    Also specifies the working directory where to store the results.
    """

    integrated_transformer_study_name: str
    integrated_transformer_optimization_directory: str

    # target parameters
    l_s12_target: float
    l_h_target: float
    n_target: float

    # fix input parameters
    time_current_1_vec: np.ndarray
    time_current_2_vec: np.ndarray

    # parameters to optimize
    material_list: list
    core_name_list: list | None
    core_inner_diameter_min_max_list: list | None
    window_w_min_max_list: list | None
    window_h_top_min_max_list: list | None
    window_h_bot_min_max_list: list | None
    n_1_top_min_max_list: list
    n_1_bot_min_max_list: list
    n_2_top_min_max_list: list
    n_2_bot_min_max_list: list
    factor_max_flux_density: float
    litz_wire_list_1: list[str]
    litz_wire_list_2: list[str]
    temperature: float

    # fix parameters: insulations
    insulations: ItoInsulation

    # data sources
    material_data_sources: IntegratedTransformerMaterialDataSources

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
class ItoTargetAndFixedParameters:
    """
    Integrated-transformer optimization target and fixed parameters.

    These parameters are calculated from the integrated-transformer input configuration (ItoSingleInputConfig).
    """

    i_rms_1: float
    i_rms_2: float
    i_peak_1: float
    i_peak_2: float
    i_phase_deg_1: float
    i_phase_deg_2: float
    material_dto_curve_list: list[MaterialCurve]
    magnet_hub_model_list: list[LossModel]
    time_extracted_vec: list
    current_extracted_1_vec: list
    current_extracted_2_vec: list
    fundamental_frequency: float
    target_inductance_matrix: np.ndarray
    working_directories: WorkingDirectories

    # winding 1
    fft_frequency_list_1: list[float]
    fft_amplitude_list_1: list[float]
    fft_phases_list_1: list[float]

    # winding 2
    fft_frequency_list_2: list[float]
    fft_amplitude_list_2: list[float]
    fft_phases_list_2: list[float]

@dataclass
class ItoSingleResultFile:
    """
    Dataclass to store the reluctance model simulation results.

    Contains concrete geometry parameters as well as the calculated results.
    """

    case: int
    # geometry parameters
    air_gap_top: float
    air_gap_bot: float
    air_gap_middle: float
    n_p_top: int
    n_p_bot: int
    n_s_top: int
    n_s_bot: int
    window_h_top: float
    window_h_bot: float
    window_w: float
    core_material: str
    core_inner_diameter: float
    litz_wire_name_1: str
    litz_wire_name_2: str

    # reluctance model results
    flux_top_max: float
    flux_bot_max: float
    flux_stray_max: float
    flux_density_top_max: float
    flux_density_bot_max: float
    flux_density_stray_max: float
    p_hyst: float
    litz_wire_loss_1: float
    litz_wire_loss_2: float
    core_2daxi_total_volume: float
    total_loss: float

@dataclass
class ItoReluctanceModelInput:
    """Input DTO for reluctance model simulation within the inductor optimization."""

    target_inductance_matrix: np.array
    core_inner_diameter: float
    window_w: float
    window_h_bot: float
    window_h_top: float
    turns_1_top: int
    turns_2_top: int
    turns_1_bot: int
    turns_2_bot: int
    litz_wire_name_1: str
    litz_wire_diameter_1: float
    litz_wire_name_2: str
    litz_wire_diameter_2: float

    insulations: ItoInsulation
    material_dto: MaterialCurve
    magnet_material_model: LossModel

    temperature: float
    time_extracted_vec: list
    current_extracted_vec_1: list
    current_extracted_vec_2: list
    fundamental_frequency: float

    i_rms_1: float
    i_rms_2: float

    litz_dict_1: dict
    litz_dict_2: dict

    # winding 1
    fft_frequency_list_1: list[float]
    fft_amplitude_list_1: list[float]
    fft_phases_list_1: list[float]

    # winding 2
    fft_frequency_list_2: list[float]
    fft_amplitude_list_2: list[float]
    fft_phases_list_2: list[float]

@dataclass
class ItoReluctanceModelOutput:
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
    l_middle_air_gap: float
    volume: float
    area_to_heat_sink: float
    p_loss: float
