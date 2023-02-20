# python libraries
from dataclasses import dataclass

# femmt libraries


# 3rd party libraries
import numpy as np


@dataclass
class InputConfig:
    l_s_target: float
    l_h_target: float
    n_target: float
    time_current_1_vec: np.ndarray
    time_current_2_vec: np.ndarray
    material_list: list
    core_inner_diameter_min_max_count_list: list
    window_w_min_max_count_list: list
    window_h_top_min_max_count_list: list
    window_h_bot_min_max_count_list: list
    factor_max_flux_density: float
    n_p_top_min_max_list: list
    n_p_bot_min_max_list: list
    n_s_top_min_max_list: list
    n_s_bot_min_max_list: list


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
    time_current_1_vec: np.ndarray
    time_current_2_vec: np.ndarray
    l_s_target_value: np.ndarray
    l_h_target_value: np.ndarray
    n_target_value: np.ndarray
    factor_max_flux_density: np.ndarray
@dataclass

class ResultFile:
    air_gap_top: float
    air_gap_bot: float
    n_p_top: int
    n_p_bot: int
    n_s_top: int
    n_s_bot: int
    window_h_top: float
    window_h_bot: float
    window_w: float
    mu_r_abs: float
    core_inner_diameter: float
    flux_top_max: float
    flux_bot_max: float
    flux_stray_max: float
    p_hyst: float
    core_2daxi_total_volume: float