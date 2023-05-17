# python libraries
from dataclasses import dataclass
from typing import List

# 3rd party libraries
import numpy as np
from materialdatabase.dtos import MaterialCurve



@dataclass
class StoInsulation:
    iso_top_core: float
    iso_bot_core: float
    iso_left_core: float
    iso_right_core: float
    iso_primary_to_primary: float
    iso_secondary_to_secondary: float
    iso_primary_to_secondary:float

@dataclass
class StoSingleInputConfig:
    """
    Configuration to simulate a stacked transformer.

    Input parameters are the target parameters, current vectors and the parameters to sweep.
    Also specifies the working directory where to store the results.
    """
    # target parameters
    l_s_target: float
    l_h_target: float
    n_target: float

    # operating point: current waveforms and temperature
    time_current_1_vec: np.ndarray
    time_current_2_vec: np.ndarray
    temperature: float

    # sweep parameters: geometry and materials
    material_list: list
    core_inner_diameter_min_max_list: list
    window_w_min_max_list: list
    window_h_top_min_max_list: list
    window_h_bot_min_max_list: list
    factor_max_flux_density: float
    primary_litz_wire_list: list
    metal_sheet_thickness: list

    # fix parameters: insulations
    insulations: StoInsulation

    # misc
    working_directory: str



@dataclass
class WorkingDirectories:
    """
    Working directories for an integrated transformer optimization
    """
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
