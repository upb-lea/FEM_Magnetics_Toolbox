# python libraries
from dataclasses import dataclass
from typing import List

# 3rd party libraries
import numpy as np


@dataclass
class ItoSingleInputConfig:
    l_s_target: float
    l_h_target: float
    n_target: float
    time_current_1_vec: np.ndarray
    time_current_2_vec: np.ndarray
    material_list: list
    core_inner_diameter_min_max_list: list
    window_w_min_max_list: list
    window_h_top_min_max_list: list
    window_h_bot_min_max_list: list
    factor_max_flux_density: float
    n_p_top_min_max_list: list
    n_p_bot_min_max_list: list
    n_s_top_min_max_list: list
    n_s_bot_min_max_list: list
    primary_litz_wire_list: list
    secondary_litz_wire_list: list
    temperature: float
    working_directory: str

@dataclass
class MaterialCurve:
    material_name: str
    material_mu_r_abs: float
    material_flux_density_vec: np.ndarray
    material_mu_r_imag_vec: np.ndarray
    material_mu_r_real_vec: np.ndarray
    saturation_flux_density: float


@dataclass
class ItoTargetAndFixedParameters:
    i_rms_1: float
    i_rms_2: float
    material_dto_curve_list: List[MaterialCurve]
    time_extracted_vec: List
    current_extracted_1_vec: List
    current_extracted_2_vec: List
    fundamental_frequency: float
    target_inductance_matrix: np.ndarray
    fem_working_directory: str
    fem_simulation_results_directory: str
    reluctance_model_results_directory: str

@dataclass
class SweepTensor:
    t1_n_p_top: np.ndarray
    t1_n_p_bot: np.ndarray
    t1_n_s_top: np.ndarray
    t1_n_s_bot: np.ndarray
    t1_window_h_top: np.ndarray
    t1_window_h_bot: np.ndarray
    t1_window_w: np.ndarray
    t1_core_material: list
    t1_core_inner_diameter: np.ndarray
    t1_primary_litz_wire_list: list
    t1_secondary_litz_wire_list: list
    time_current_1_vec: np.ndarray
    time_current_2_vec: np.ndarray
    l_s_target_value: float
    l_h_target_value: float
    n_target_value: float
    factor_max_flux_density: float
@dataclass
class ItoSingleResultFile:
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
    primary_litz_wire: str
    secondary_litz_wire: str

    # reluctance model results
    flux_top_max: float
    flux_bot_max: float
    flux_stray_max: float
    flux_density_top_max: float
    flux_density_bot_max: float
    flux_density_stray_max: float
    p_hyst: float
    primary_litz_wire_loss: float
    secondary_litz_wire_loss: float
    core_2daxi_total_volume: float
    total_loss: float

    # fem simulation results
