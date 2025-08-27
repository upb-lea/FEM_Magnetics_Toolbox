"""Includes the DTOs to perform a transformer optimization."""
# python libraries
from dataclasses import dataclass

# 3rd party libraries
import numpy as np
from femmt.optimization.optimization_dtos import MaterialDataSources, WorkingDirectories

@dataclass
class ToInsulation:
    """Insulation settings for the transformer optimization."""

    iso_top_core: float
    iso_bot_core: float
    iso_left_core_min: float
    iso_right_core: float
    iso_primary_to_primary: float
    iso_secondary_to_secondary: float
    iso_primary_to_secondary: float
    iso_primary_inner_bobbin: float


@dataclass
class ToSingleInputConfig:
    """
    Configuration to simulate a stacked transformer.

    Input parameters are the target parameters, current vectors and the parameters to sweep.
    Also specifies the working directory where to store the results.
    """

    # target parameters
    l_s12_target: float
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
    window_h_min_max_list: list
    primary_litz_wire_list: list
    metal_sheet_thickness_list: list
    interleaving_scheme_list: list
    interleaving_type_list: list

    # maximum limitation for transformer total height and core volume
    max_core_volume: float

    # fix parameters: insulations
    insulations: ToInsulation

    # misc
    working_directory: str
    fft_filter_value_factor: float
    mesh_accuracy: float

    # data sources
    material_datasources: MaterialDataSources


@dataclass
class ThermalConfig:
    """Thermal configuration for a transformer optimization."""

    thermal_conductivity_dict: dict
    case_gap_top: float
    case_gap_right: float
    case_gap_bot: float
    boundary_temperatures: dict
    boundary_flags: dict

@dataclass
class ToTargetAndFixedParameters:
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
    material_name_list: list[str]
    material_complex_mu_r_list: list[float]
    time_extracted_vec: list
    current_extracted_1_vec: list
    current_extracted_2_vec: list
    fundamental_frequency: float
    target_inductance_matrix: np.ndarray
    working_directories: WorkingDirectories


@dataclass
class CurrentWorkingPoint:
    """Stores the working point of currents together with a human-readable name."""

    name: str
    time_current_1_vec: np.ndarray | list
    time_current_2_vec: np.ndarray | list
