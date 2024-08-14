"""Contains different functions, used by the whole FEMMT functions."""
# Python standard libraries
import json
import pkg_resources
import subprocess
import sys
import os
import warnings
from typing import Union, List, Tuple, Dict
from scipy.integrate import quadrature


# Third parry libraries
import gmsh
from matplotlib import pyplot as plt
import numpy.typing as npt
import numpy as np

# Local libraries
from femmt.constants import *
from femmt.enumerations import ConductorType
from femmt.dtos import *

colors_femmt_default = {"blue": (28, 113, 216),
                        'red': (192, 28, 40),
                        "green": (46, 194, 126),
                        "orange": (230, 97, 0),
                        "purple": (129, 61, 156),
                        "brown": (134, 94, 60),
                        "grey": (119, 118, 123),
                        "yellow": (245, 194, 17),
                        "black": (0, 0, 0),
                        "white": (255, 255, 255)
                        }

colors_geometry_femmt_default = {
    "core": "black",
    "air_gap": "yellow",
    "winding": ["orange", "brown", "yellow", "green", "red", "black", "grey", "blue", "orange", "purple", "grey",
                "blue", "orange", "purple"],
    "insulation": "blue",
    "potting_inner": "grey",
    "potting_outer": "grey",
}

colors_ba_jonas = {"blue": (28, 113, 216),
                   'red': (213, 6, 6),
                   "green": (6, 213, 6),
                   "orange": (230, 97, 0),
                   "purple": (129, 61, 156),
                   "brown": (134, 94, 60),
                   "grey": (193, 193, 193),
                   "yellow": (255, 171, 6),
                   "black": (58, 58, 58),
                   "white": (255, 255, 255),
                   "grey_dark": (109, 109, 109),
                   "grey_dark_dark": (50, 50, 50)
                   }

colors_geometry_ba_jonas = {
    "core": "black",
    "air_gap": "yellow",
    "winding": ["green", "red", "yellow"],
    "insulation": "grey_dark",
    "potting_inner": "grey",
    "potting_outer": "grey_dark_dark",
}

colors_geometry_draw_only_lines = {
    "core": "grey_dark",
    "air_gap": "grey_dark",
    "winding": ["green", "green", "green"],
    "insulation": "grey_dark",
    "potting_inner": "grey_dark",
    "potting_outer": "grey_dark",
}


def core_database() -> Dict:
    """
    Return a core geometry for defined core structure.

    All dimensions are nominal dimensions without consideration of tolerances.

    For PQ core sizes (e.g. PQ 40/40), it has been found out that
    core_dimension_x / core_dimension_y = 1.45, the error over all available shapes is maximum 7%
    (compared to datasheet value)
    Derivation:
    core_list: ['PQ 20/20', 'PQ 26/25', 'PQ 32/30', 'PQ 35/35', 'PQ 40/40', 'PQ 50/50']
    factor core_dimension_x / core_dimension_y = [1.46, 1.39, 1.45, 1.35, 1.44, 1.56]
    mean over the factors = 1.45
    max derivation / mean = 1.07 (< 7% accuracy)

    :return: Dict including core_h, core_inner_diameter, window_h, window_w
    :rtype: Dict
    """
    core_dict = {}

    # -----------------------
    # PQ Cores
    # -----------------------

    core_dict["PQ 16/11.6"] = {
        "core_h": 11.6e-3,
        "core_inner_diameter": 7e-3,
        "window_h": 6.7e-3,
        "window_w": 3.7e-3,
        "core_dimension_x": 16.4e-3,
        "core_dimension_y": 11.2e-3,
    }
    core_dict["PQ 20/16"] = {
        "core_h": 16.2e-3,
        "core_inner_diameter": 8.8e-3,
        "window_h": 10.3e-3,
        "window_w": (18 - 8.8) / 2 * 1e-3,
        "core_dimension_x": 20.5e-3,
        "core_dimension_y": 14.0e-3,
    }
    core_dict["PQ 20/20"] = {
        "core_h": 20.2e-3,
        "core_inner_diameter": 8.8e-3,
        "window_h": 14.3e-3,
        "window_w": 4.6e-3,
        "core_dimension_x": 20.5e-3,
        "core_dimension_y": 14.0e-3,
    }
    core_dict["PQ 26/20"] = {
        "core_h": 20.16e-3,
        "core_inner_diameter": 12e-3,
        "window_h": 11.5e-3,
        "window_w": 5.25e-3,
        "core_dimension_x": 26.5e-3,
        "core_dimension_y": 19.0e-3,
    }
    core_dict["PQ 26/25"] = {
        "core_h": 24.76e-3,
        "core_inner_diameter": 12e-3,
        "window_h": 16.1e-3,
        "window_w": (22.5 - 12) / 2 * 1e-3,
        "core_dimension_x": 26.5e-3,
        "core_dimension_y": 19.0e-3,
    }
    core_dict["PQ 32/20"] = {
        "core_h": 20.5e-3,
        "core_inner_diameter": 13.45e-3,
        "window_h": 11.5e-3,
        "window_w": (27.5 - 13.45) / 2 * 1e-3,
        "core_dimension_x": 32.0e-3,
        "core_dimension_y": 22.0e-3,
    }
    core_dict["PQ 32/30"] = {
        "core_h": 30.35e-3,
        "core_inner_diameter": 13.45e-3,
        "window_h": 21.3e-3,
        "window_w": (27.5 - 13.45) / 2 * 1e-3,
        "core_dimension_x": 32.0e-3,
        "core_dimension_y": 22.0e-3,
    }
    core_dict["PQ 35/35"] = {
        "core_h": 34.8e-3,
        "core_inner_diameter": 14.35e-3,
        "window_h": 25e-3,
        "window_w": (32 - 14.35) / 2 * 1e-3,
        "core_dimension_x": 35.1e-3,
        "core_dimension_y": 26.0e-3,
    }
    core_dict["PQ 40/30"] = {
        "core_h": 30.3e-3,
        "core_inner_diameter": 14.9e-3,
        "window_h": 20e-3,
        "window_w": (37 - 14.9) / 2 * 1e-3,
        "core_dimension_x": 40.5e-3,
        "core_dimension_y": 28.0e-3,
    }
    core_dict["PQ 40/40"] = {
        "core_h": 39.8e-3,
        "core_inner_diameter": 14.9e-3,
        "window_h": 29.5e-3,
        "window_w": (37 - 14.9) / 2 * 1e-3,
        "core_dimension_x": 40.5e-3,
        "core_dimension_y": 28.0e-3,
    }
    core_dict["PQ 50/40"] = {
        "core_h": 40e-3,
        "core_inner_diameter": 20e-3,
        "window_h": 26.1e-3,
        "window_w": (44 - 20) / 2 * 1e-3,
        "core_dimension_x": 50.0e-3,
        "core_dimension_y": 32.0e-3,
    }
    core_dict["PQ 50/50"] = {
        "core_h": 50e-3,
        "core_inner_diameter": 20e-3,
        "window_h": 36.1e-3,
        "window_w": (44 - 20) / 2 * 1e-3,
        "core_dimension_x": 50.0e-3,
        "core_dimension_y": 32.0e-3,
    }
    core_dict["PQ 65/60"] = {
        "core_h": 60e-3,
        "core_inner_diameter": 26e-3,
        "window_h": 42e-3,
        "window_w": (65 - 26) / 2 * 1e-3,
        "core_dimension_x": 65.0e-3,
        "core_dimension_y": 45.0e-3,
    }
    # -----------------------
    # PM Cores
    # -----------------------

    core_dict["PM 114/93"] = {
        "core_h": 93 * 1e-3,
        "core_inner_diameter": pm_core_inner_diameter_calculator(43e-3, 5.4e-3),
        "window_h": 63 * 1e-3,
        "window_w": (88 - 43) / 2 * 1e-3,
    }
    core_dict["PM 50/39"] = {
        "core_h": 39 * 1e-3,
        "core_inner_diameter": pm_core_inner_diameter_calculator(20e-3, 5.4e-3),
        "window_h": 26.4 * 1e-3,
        "window_w": (39 - 20) / 2 * 1e-3,
    }
    core_dict["PM 62/49"] = {
        "core_h": 49 * 1e-3,
        "core_inner_diameter": pm_core_inner_diameter_calculator(25.5e-3, 5.4e-3),
        "window_h": 33.4 * 1e-3,
        "window_w": (48.8 - 25.5) / 2 * 1e-3,
    }
    core_dict["PM 74/59"] = {
        "core_h": 59 * 1e-3,
        "core_inner_diameter": pm_core_inner_diameter_calculator(29.5e-3, 5.4e-3),
        "window_h": 40.7e-3,
        "window_w": (57.5 - 29.5) / 2 * 1e-3,
    }
    core_dict["PM 87/70"] = {
        "core_h": 70 * 1e-3,
        "core_inner_diameter": pm_core_inner_diameter_calculator(31.7e-3, 8.5e-3),
        "window_h": 48 * 1e-3,
        "window_w": (67.1 - 31.7) / 2 * 1e-3,
    }
    return core_dict


def litz_database() -> Dict:
    """
    Return litz parameters for defined litz wires.

    :return: Dict including litz parameters like strand_numbers, strand_radii and conductor_radii
    :rtype: Dict
    """
    litz_dict = {}

    litz_dict["1.5x105x0.1"] = {"strands_numbers": 105,
                                "strand_radii": 0.1e-3 / 2,
                                "conductor_radii": 1.5e-3 / 2,
                                "ff": "",
                                "manufacturer": "PACK",
                                "material_number": "",
                                "litz": "RUPALIT V155",
                                "insulation": "textile"}
    litz_dict["1.4x200x0.071"] = {"strands_numbers": 200,
                                  "strand_radii": 0.071e-3 / 2,
                                  "conductor_radii": 1.4e-3 / 2,
                                  "ff": "",
                                  "manufacturer": "PACK",
                                  "material_number": "",
                                  "litz": "RUPALIT V155",
                                  "insulation": "textile"}
    litz_dict["2.0x405x0.071"] = {"strands_numbers": 405,
                                  "strand_radii": 0.071e-3 / 2,
                                  "conductor_radii": 2.0e-3 / 2,
                                  "ff": "",
                                  "manufacturer": "",
                                  "material_number": "",
                                  "litz": "",
                                  "insulation": "unknown blue plastic"}
    litz_dict["2.0x800x0.05"] = {"strands_numbers": 800,
                                 "strand_radii": 0.05e-3 / 2,
                                 "conductor_radii": 2e-3 / 2,
                                 "ff": "",
                                 "manufacturer": "Elektrisola",
                                 "material_number": "12104184",
                                 "litz": "",
                                 "insulation": ""
                                 }
    litz_dict["1.1x60x0.1"] = {"strands_numbers": 60,
                               "strand_radii": 0.1e-3 / 2,
                               "conductor_radii": 1.1e-3 / 2,
                               "ff": "",
                               "manufacturer": "PACK",
                               "material_number": "",
                               "litz": "RUPALIT V155",
                               "insulation": "textile"
                               }
    litz_dict["1.35x200x0.071"] = {"strands_numbers": 200,
                                   "strand_radii": 0.071e-3 / 2,
                                   "conductor_radii": 1.35e-3 / 2,
                                   "ff": "",
                                   "manufacturer": "PACK",
                                   "material_number": "",
                                   "litz": "RUPALIT V155",
                                   "insulation": "textile"}

    litz_dict["3.2x2100x0.05"] = {"strands_numbers": 2100,
                                  "strand_radii": 0.05e-3 / 2,
                                  "conductor_radii": 3.2e-3 / 2,
                                  "ff": "",
                                  "manufacturer": "PACK",
                                  "material_number": "AB21220373",
                                  "litz": "RUPALIT V155",
                                  "insulation": "textile"
                                  }

    litz_dict["4.6x2160x0.071"] = {"strands_numbers": 2160,
                                   "strand_radii": 0.071e-3 / 2,
                                   "conductor_radii": 4.6e-3 / 2,
                                   "ff": "",
                                   "manufacturer": "PACK",
                                   "material_number": "AB21225497",
                                   "litz": "RUPALIT V155",
                                   "insulation": "textile"
                                   }

    litz_dict["2.9x1200x0.06"] = {"strands_numbers": 1200,
                                  "strand_radii": 0.06e-3 / 2,
                                  "conductor_radii": 2.9e-3 / 2,
                                  "ff": "",
                                  "manufacturer": "Elektrisola",
                                  "material_number": "",
                                  "litz": "",
                                  "insulation": "unknown plastic"}

    litz_dict["2.6x1000x0.06"] = {"strands_numbers": 1000,
                                  "strand_radii": 0.06e-3 / 2,
                                  "conductor_radii": 2.6e-3 / 2,
                                  "ff": "",
                                  "manufacturer": "Elektrisola",
                                  "material_number": "",
                                  "litz": "",
                                  "insulation": "unknown plastic"}

    litz_dict["1.8x512x0.05"] = {"strands_numbers": 512,
                                 "strand_radii": 0.05e-3 / 2,
                                 "conductor_radii": 1.8e-3 / 2,
                                 "ff": "",
                                 "manufacturer": "PACK",
                                 "material_number": "AB21217207",
                                 "litz": "RUPALIT Safety VB155",
                                 "insulation": "3 layers Mylar"}

    litz_dict["2.3x600x0.071"] = {"strands_numbers": 600,
                                  "strand_radii": 0.071e-3 / 2,
                                  "conductor_radii": 2.3e-3 / 2,
                                  "ff": "",
                                  "manufacturer": "PACK",
                                  "material_number": "AB21220522",
                                  "litz": "RUPALIT Safety Profil V155",
                                  "insulation": "3 layers Mylar"}

    litz_dict["2.8x400x0.1"] = {"strands_numbers": 400,
                                "strand_radii": 0.1e-3 / 2,
                                "conductor_radii": 2.8e-3 / 2,
                                "ff": "",
                                "manufacturer": "PACK",
                                "material_number": "AB21222210",
                                "litz": "RUPALIT Safety V155",
                                "insulation": "3 layers Mylar"}

    litz_dict["1.71x140x0.1"] = {"strands_numbers": 140,
                                 "strand_radii": 0.1e-3 / 2,
                                 "conductor_radii": 1.71e-3 / 2,
                                 "ff": "",
                                 "manufacturer": "",
                                 "material_number": "",
                                 "litz": "",
                                 "insulation": ""}

    litz_dict["1.7x500x0.06"] = {"strands_numbers": 500,
                                 "strand_radii": 0.06e-3 / 2,
                                 "conductor_radii": 1.7e-3 / 2,
                                 "ff": "",
                                 "manufacturer": "",
                                 "material_number": "",
                                 "litz": "",
                                 "insulation": ""}

    return litz_dict


def wire_material_database() -> Dict[str, WireMaterial]:
    """
    Return wire materials e.g. copper, aluminium in a dictionary.

    :return: Dict with materials and conductivity
    :rtype: Dict
    """
    wire_material = {}

    wire_material["Copper"] = WireMaterial(
        name="copper",
        sigma=5.8e7,
        temperature=25,
        temperature_coefficient=3.9e-3,
        thermal_conductivity=400,
        volumetric_mass_density=8920,
    )

    wire_material["Aluminium"] = WireMaterial(
        name="aluminium",
        sigma=3.7e7,
        temperature=25,
        temperature_coefficient=3.9e-3,
        thermal_conductivity=235,
        volumetric_mass_density=2699,
    )

    return wire_material


def conductivity_temperature(material: str, temperature: float) -> float:
    """
    Calculate the conductivity for a certain temperature of the material.

    :param material: material name, e.g. "copper"
    :type material: str
    :param temperature: temperature in °C
    :type temperature: float
    :return: conductivity of material at given temperature
    :rtype: float
    """
    material_from_database = wire_material_database()[material]

    sigma_database = material_from_database.sigma
    temperature_database = material_from_database.temperature
    temperature_coefficient_database = material_from_database.temperature_coefficient

    resistance_temperature = 1 / sigma_database * (
        1 + temperature_coefficient_database * (temperature - temperature_database))
    sigma_temperature = 1 / resistance_temperature

    return sigma_temperature


def create_folders(*args) -> None:
    """Create folders for every given folder path (if it does not exist)."""
    for folder in list(args):
        if not os.path.exists(folder):
            os.mkdir(folder)


def cost_material_database() -> Dict:
    """
    Return costs for core and winding. This is split in material and fabrication costs.

    Both, material and fabrication costs have a euro_per_kilogram and a euro_per_unit (fix costs) price.

    Source: R. Burkart and J. Kolar 'Component Cost Models for Multi-Objective Optimizations of #
    Switched-Mode Power Converter' 2013.

    These are outdated prices (year 2013). Update needed in the future.
    """
    cost_database = {}
    cost_database["core_euro_per_kilogram"] = {"ferrite": 5.5,
                                               "amorphous": 16,
                                               "nanocristalline": 23,
                                               "high_si_steel": 12,
                                               "goes": 2.5}
    cost_database["winding_material_euro_per_kilogram"] = {ConductorType.RoundSolid.name: 10,
                                                           "flat": 10,
                                                           ConductorType.RectangularSolid.name: 20,
                                                           ConductorType.RoundLitz.name: -1}

    cost_database["winding_material_euro_per_unit"] = {ConductorType.RoundSolid.name: 1,
                                                       "flat": 2,
                                                       ConductorType.RectangularSolid.name: 2,
                                                       ConductorType.RoundLitz.name: 1}

    cost_database["winding_fabrication_euro_per_kilogram"] = {ConductorType.RoundSolid.name: 7,
                                                              "flat": 21,
                                                              ConductorType.RectangularSolid.name: 14,
                                                              ConductorType.RoundLitz.name: 7}
    cost_database["winding_fabrication_euro_per_unit"] = {ConductorType.RoundSolid.name: 2,
                                                          "flat": 4,
                                                          ConductorType.RectangularSolid.name: 2.5,
                                                          ConductorType.RoundLitz.name: 2}

    cost_database["winding_material_euro_per_kilogram_for_litz"] = {"sigma_numerator": 15,
                                                                    "sigma_denominator": 0.45}

    cost_database["gross_margin"] = 0.25

    return cost_database


def pm_core_inner_diameter_calculator(inner_core_diameter: float, hole_diameter: float) -> np.array:
    """
    Calculate the effective inner core diameter without the hole often used in PM-cores.

    :param inner_core_diameter: inner core diameter
    :type inner_core_diameter: float
    :param hole_diameter: hole diameter
    :type hole_diameter: float
    :return: effective inner core diameter without hole
    :rtype: np.array
    """
    area_inner_core_inner_diameterithout_hole = (inner_core_diameter / 2) ** 2 * np.pi
    area_hole = (hole_diameter / 2) ** 2 * np.pi
    area_total = area_inner_core_inner_diameterithout_hole - area_hole

    return np.around(2 * np.sqrt(area_total / np.pi), decimals=4)


def install_pyfemm_if_missing() -> None:
    """Installs femm-software pip package in case of running on Windows machine. Windows users only."""
    required = {'pyfemm'}
    installed = {pkg.key for pkg in pkg_resources.working_set}
    missing = required - installed

    if missing:
        print("Missing 'pyfemm' installation.")
        print("Installing 'pyfemm' ...")
        python = sys.executable
        subprocess.check_call([python, '-m', 'pip', 'install', *missing], stdout=subprocess.DEVNULL)
        print("'pyfemm' is now installed!")

def litz_calculate_number_strands(n_layers: int) -> int:
    """
    Return the number of strands in a hexagonal litz winding with a specified number of layers (n_layers).

    CAUTION: Zero number of layers corresponds to a single strand.

    :param n_layers: number of litz_layers
    :type n_layers: int

    :return: number of strands in a litz wire
    :rtype: int

    """
    return 3 * (n_layers + 1) ** 2 - 3 * (n_layers + 1) + 1


def litz_calculate_number_layers(n_strands: int) -> int:
    """
    Return the number of layers in a hexagonal litz winding with a specified number of strands (n_strands).

    .. note:: Zero number of layers corresponds to a single strand.

    :param n_strands: Number of strands in a litz
    :type n_strands: int

    :return: number of layers for a litz
    :rtype: int
    """
    return np.sqrt(0.25 + (n_strands - 1) / 3) - 0.5


def fft(period_vector_t_i: npt.ArrayLike, sample_factor: int = 1000, plot: str = 'no', mode: str = 'rad',
        f0: Union[float, None] = None, title: str = 'ffT', filter_type: str = 'factor',
        filter_value_factor: float = 0.01, filter_value_harmonic: int = 100,
        figure_size: Tuple = None, figure_directory: str = None) -> npt.NDArray[list]:
    """
    Calculate the FFT for a given input signal. Input signal is in vector format and should include one period.

    Output vector includes only frequencies with amplitudes > 1% of input signal

    :Minimal Example:

    >>> import femmt as fmt
    >>> import numpy as np
    >>> example_waveform = np.array([[0, 1.34, 3.14, 4.48, 6.28],[-175.69, 103.47, 175.69, -103.47,-175.69]])
    >>> out = fmt.fft(example_waveform, plot='yes', mode='rad', f0=25000, title='ffT input current')

    :param period_vector_t_i: numpy-array [[time-vector[,[current-vector]]. One period only
    :type period_vector_t_i: np.array
    :param sample_factor: f_sampling/f_period, defaults to 1000
    :type sample_factor: int
    :param plot: insert anything else than "no" or 'False' to show a plot to visualize input and output
    :type plot: str
    :param mode: 'rad'[default]: full period is 2*pi, 'deg': full period is 360°, 'time': time domain.
    :type mode: str
    :param f0: fundamental frequency. Needs to be set in 'rad'- or 'deg'-mode
    :type f0: float
    :param title: plot window title, defaults to 'ffT'
    :type title: str
    :param filter_type: 'factor'[default] or 'harmonic' or 'disabled'.
    :type filter_type: str
    :param filter_value_factor: filters out amplitude-values below a certain factor of max. input amplitude.
        Should be 0...1, default to 0.01 (1%)
    :type filter_value_factor: float
    :param filter_value_harmonic: filters out harmonics up to a certain number. Default value is 100.
        Note: count 1 is DC component, count 2 is the fundamental frequency
    :type filter_value_harmonic: int
    :param figure_directory: full path with file extension
    :type figure_directory: Tuple
    :param figure_size: None for auto-fit; fig_size for matplotlib (width, length)
    :type figure_size: Tuple

    :return: numpy-array [[frequency-vector],[amplitude-vector],[phase-vector]]
    :rtype: npt.NDArray[list]
    """
    # check for correct input parameter
    if (mode == 'rad' or mode == 'deg') and f0 is None:
        raise ValueError("if mode is 'rad' or 'deg', a fundamental frequency f0 must be set")
    # check for input is list. Convert to numpy-array
    if isinstance(period_vector_t_i, list):
        if plot != 'no' and plot is not False:
            print("Input is list, convert to np.array()")
        period_vector_t_i = np.array(period_vector_t_i)

    # first value of time vector must be zero
    if period_vector_t_i[0][0] != 0:
        raise ValueError("Period vector must start with 0 seconds!")

    # mode pre-calculation
    if mode == 'rad':
        period_vector_t_i[0] = period_vector_t_i[0] / (2 * np.pi * f0)
    elif mode == 'deg':
        period_vector_t_i[0] = period_vector_t_i[0] / (360 * f0)
    elif mode != 'time':
        raise ValueError("Mode not available. Choose: 'rad', 'deg', 'time'")

    t = period_vector_t_i[0]
    i = period_vector_t_i[1]

    # fft-function works per default in time domain
    t_interp = np.linspace(0, t[-1], sample_factor)
    i_interp = np.interp(t_interp, t, i)

    f0 = round(1 / t[-1])
    Fs = f0 * sample_factor

    # frequency domain
    f = np.linspace(0, (sample_factor - 1) * f0, sample_factor)
    x = np.fft.fft(i_interp)
    x_mag = np.abs(x) / sample_factor
    phi_rad = np.angle(x)

    f_corrected = f[0:int(sample_factor / 2 + 1)]
    x_mag_corrected = 2 * x_mag[0:int(sample_factor / 2 + 1)]
    x_mag_corrected[0] = x_mag_corrected[0] / 2
    phi_rad_corrected = phi_rad[0:int(sample_factor / 2 + 1)]

    f_out = []
    x_out = []
    phi_rad_out = []
    if filter_type.lower() == 'factor':
        for count, _ in enumerate(x_mag_corrected):
            if x_mag_corrected[count] > filter_value_factor * max(abs(i)):
                f_out.append(f_corrected[count])
                x_out.append(x_mag_corrected[count])
                phi_rad_out.append(phi_rad_corrected[count])
    elif filter_type.lower() == 'harmonic':
        for count, _ in enumerate(x_mag_corrected):
            if count < filter_value_harmonic:
                f_out.append(f_corrected[count])
                x_out.append(x_mag_corrected[count])
                phi_rad_out.append(phi_rad_corrected[count])
    elif filter_type.lower() == 'disabled':
        f_out = f_corrected
        x_out = x_mag_corrected
        phi_rad_out = phi_rad_corrected
    else:
        raise ValueError(
            f"filter_type '{filter_value_harmonic}' not available: Must be 'factor','harmonic' or 'disabled ")

    if plot != 'no' and plot is not False:
        print(f"{title=}")
        print(f"{t[-1]=}")
        print(f"{f0=}")
        print(f"{Fs=}")
        print(f"{sample_factor=}")
        print(f"f_out = {np.around(f_out, decimals=0)}")
        print(f"x_out = {np.around(x_out, decimals=3)}")
        print(f"phi_rad_out = {np.around(phi_rad_out, decimals=3)}")

        reconstructed_signal = 0
        for i_range in range(len(f_out)):
            reconstructed_signal += x_out[i_range] * np.cos(
                2 * np.pi * f_out[i_range] * t_interp + phi_rad_out[i_range])

        fig, [ax1, ax2, ax3] = plt.subplots(num=title, nrows=3, ncols=1, figsize=figure_size)
        ax1.plot(t, i, label='original signal')
        ax1.plot(t_interp, reconstructed_signal, label='reconstructed signal')
        ax1.grid()
        ax1.set_title('Signal')
        ax1.set_xlabel('time in s')
        ax1.set_ylabel('Amplitude')
        ax1.legend()
        ax2.stem(f_out, x_out)
        ax2.grid()
        ax2.set_title('ffT')
        ax2.set_xlabel('Frequency in Hz')
        ax2.set_ylabel('Amplitude')

        ax3.stem(f_out, phi_rad_out)
        ax3.grid()
        ax3.set_title('ffT')
        ax3.set_xlabel('Frequency in Hz')
        ax3.set_ylabel('Phase in rad')

        plt.tight_layout()
        if figure_directory is not None:
            plt.savefig(figure_directory, bbox_inches="tight")
        plt.show()

    return np.array([f_out, x_out, phi_rad_out])


def plot_fourier_coefficients(frequency_list, amplitude_list, phi_rad_list, sample_factor: int = 1000,
                              figure_directory: str = None):
    """Plot fourier coefficients in a visual figure."""
    # dc and ac handling
    nonzero_frequencies = [f for f in frequency_list if f != 0]
    if nonzero_frequencies:
        time_period = 1 / min(nonzero_frequencies)
    else:
        time_period = 1
    # time_period = 1 / min(frequency_list)

    t_interp = np.linspace(0, time_period, sample_factor)
    reconstructed_signal = 0
    for i_range, _ in enumerate(frequency_list):
        reconstructed_signal += amplitude_list[i_range] * np.cos(
            2 * np.pi * frequency_list[i_range] * t_interp + phi_rad_list[i_range])

    fig, [ax1, ax2, ax3] = plt.subplots(nrows=3, ncols=1)
    ax1.plot(t_interp, reconstructed_signal, label='reconstructed signal')
    ax1.grid()
    ax1.set_title('Signal')
    ax1.set_xlabel('time in s')
    ax1.set_ylabel('Amplitude')
    ax1.legend()
    ax2.stem(frequency_list, amplitude_list)
    ax2.grid()
    ax2.set_title('ffT')
    ax2.set_xlabel('Frequency in Hz')
    ax2.set_ylabel('Amplitude')

    ax3.stem(frequency_list, phi_rad_list)
    ax3.grid()
    ax3.set_title('ffT')
    ax3.set_xlabel('Frequency in Hz')
    ax3.set_ylabel('Phase in rad')

    plt.tight_layout()
    if figure_directory is not None:
        plt.savefig(figure_directory, bbox_inches="tight")
    plt.close('all')  # close the figures to remove the warning when you run many figures

    # plt.show()


def compare_fft_list(input_data_list: list, sample_factor: int = 1000, mode: str = 'rad',
                     f0: Union[float, None] = None) -> None:
    """
    Generate fft curves from input curves and compare them to each other.

    :Minimal Example:

    >>> example_waveform = np.array([[0, 1.34, 3.14, 4.48, 6.28],[-175.69, 103.47, 175.69, -103.47,-175.69]])
    >>> example_waveform_2 = np.array([[0, 0.55, 3.14, 3.69, 6.28],[-138.37, 257.58, 138.37, -257.58, -138.37]])
    >>> compare_fft_list([example_waveform, example_waveform_2], mode='rad', f0=25000)

    :param input_data_list: list of fft-compatible numpy-arrays [element, element, ... ], each element format
        like [[time-vector[,[current-vector]]. One period only
    :param mode: 'rad'[default]: full period is 2*pi, 'deg': full period is 360°, 'time': time domain.
    :type mode: str
    :param f0: fundamental frequency. Needs to be set in 'rad'- or 'deg'-mode
    :type f0: float
    :param sample_factor: samle factor, defaults to 1000
    :type sample_factor: int
    """
    out = []
    for count, _ in enumerate(input_data_list):
        out.append([fft(input_data_list[count], sample_factor=sample_factor, plot='no', mode=mode, f0=f0)])

    fig, axs = plt.subplots(2, len(input_data_list), sharey=True)
    for count, _ in enumerate(input_data_list):
        axs[0, count].plot(input_data_list[count][0], input_data_list[count][1], label='original signal')
        axs[0, count].grid()
        axs[0, count].set_xlabel('time in s')
        axs[0, count].set_ylabel('Amplitude')
        axs[1, count].stem(out[count][0][0], out[count][0][1])
        axs[1, count].grid()
        axs[1, count].set_xlabel('frequency in Hz')
        axs[1, count].set_ylabel('Amplitude')

        # ax1.plot(t_interp, reconstructed_signal, label='reconstructed signal')

    plt.tight_layout()
    plt.show()


def store_as_npy_in_directory(dir_path: str, file_name: str, numpy_data) -> None:
    """
    Store a numpy array in a given directory.

    :param dir_path: directory path
    :type dir_path: str
    :param file_name: file name
    :type file_name: str
    :param numpy_data: numpy array
    :type numpy_data:
    """
    if not os.path.isdir(dir_path):
        os.mkdir(dir_path)
    np.save(dir_path + "/" + file_name, numpy_data)


def get_dicts_with_keys_and_values(data, **kwargs) -> Dict:
    """
    Return a list of dictionaries out of a list of dictionaries which contains pairs of the given key(s) and value(s).

    :param data: list of dicts
    :type data: List
    :param kwargs: keys and values in dicts
    """
    invalid_index = []
    for n, dictionary in enumerate(data):
        for key, value in kwargs.items():
            if not (key in dictionary and value == dictionary[key]):
                invalid_index.append(n)
                break
    valid_data = np.delete(data, invalid_index)
    return valid_data


def get_dict_with_unique_keys(data: list[dict], *keys) -> Dict:
    """
    Return a dictionary out of a list of dictionaries which contains the given key(s).

    :param data: list of dicts
    :type data: list[dict]
    :param keys: keys in dicts
    :return:
    """
    invalid_index = []
    for n, dictionary in enumerate(data):
        for key in keys:
            if key not in dictionary:
                invalid_index.append(n)
                break
    valid_data = np.delete(data, invalid_index)
    if len(valid_data) != 1:
        warnings.warn("Keyword(s) not unique!", stacklevel=2)
    # Only one dictionary shall survive --> choose element 0
    return valid_data[0]


def find_common_frequencies(frequency_list_1: List, amplitude_list_1: List, phase_list_1_rad_or_deg: List,
                            frequency_list_2: List, amplitude_list_2: List, phase_list_2_rad_or_deg: List) -> List:
    """
    Find common frequencies and returns a list of intersections.

    :param amplitude_list_1: Amplitudes signal 1
    :type amplitude_list_1: List
    :param phase_list_1_rad_or_deg: Phases signal 1, can be degree or rad. return is same as input.
    :type phase_list_1_rad_or_deg: List
    :param frequency_list_1: Frequencies signal 1
    :type frequency_list_1: List
    :param amplitude_list_2: Amplitudes signal 2
    :type amplitude_list_2: List
    :param phase_list_2_rad_or_deg: Phases signal 2, can be degree or rad. return is same as input
    :type phase_list_2_rad_or_deg: List
    :param frequency_list_2: Frequencies signal 2
    :type frequency_list_2: List
    :return: [current_pair_list, phase_pair_list, common_frequency_list]
    :rtype: Tuple

    :Example:

    >>> import femmt as fmt
    >>> frequency_1 = [50, 100, 150, 200]
    >>> frequency_2 = [50, 100, 150, 170, 200]
    >>> amplitude_1 = [1, 2, 3, 4]
    >>> amplitude_2 = [5, 6, 7, 8, 9]
    >>> phase_1 = [10, 20, 30, 40]
    >>> phase_2 = [50, 60, 70, 80, 90]
    >>> common_f, common_a, common_phase = fmt.find_common_frequencies(frequency_1, amplitude_1, phase_1,
    >>>     frequency_2, amplitude_2, phase_2)
    :Returns:
    >>> common_f = [200, 50, 100, 150]
    >>> common_a = [[4, 9], [1, 5], [2, 6], [3, 7]]
    >>> common_phase = [[40, 90], [10, 50], [20, 60], [30, 70]]
    """
    common_amplitude_list_1 = []
    common_phase_list_1 = []
    common_amplitude_list_2 = []
    common_phases_list_2 = []

    common_frequency_list = list(set(frequency_list_1).intersection(frequency_list_2))
    print(f"{common_frequency_list=}")
    common_frequency_list.sort()
    print(f"{common_frequency_list=}")

    # Delete the corresponding phases and amplitudes
    if isinstance(amplitude_list_1, List):
        for frequency in common_frequency_list:
            common_amplitude_list_1.append(amplitude_list_1[frequency_list_1.index(frequency)])
            common_phase_list_1.append(phase_list_1_rad_or_deg[frequency_list_1.index(frequency)])
            common_amplitude_list_2.append(amplitude_list_2[frequency_list_2.index(frequency)])
            common_phases_list_2.append(phase_list_2_rad_or_deg[frequency_list_2.index(frequency)])
    elif isinstance(amplitude_list_1, np.ndarray):
        for frequency in common_frequency_list:
            common_amplitude_list_1.append(amplitude_list_1[np.where(frequency_list_1 == frequency)][0])
            common_phase_list_1.append(phase_list_1_rad_or_deg[np.where(frequency_list_1 == frequency)][0])
            common_amplitude_list_2.append(amplitude_list_2[np.where(frequency_list_2 == frequency)][0])
            common_phases_list_2.append(phase_list_2_rad_or_deg[np.where(frequency_list_2 == frequency)][0])
    else:
        warnings.warn("Either a list or a np.ndarray must be provided!", stacklevel=2)

    current_pair_list = list(map(list, zip(common_amplitude_list_1, common_amplitude_list_2)))
    phase_pair_list = list(map(list, zip(common_phase_list_1, common_phases_list_2)))

    return [common_frequency_list, current_pair_list, phase_pair_list]


def sort_out_small_harmonics(frequency_list: List, amplitude_pair_list: List,
                             phase_pair_list_rad_or_deg: List, sort_out_factor: float) -> List:
    """
    Sort out small harmonics from a given fft-output of a signal.

    :param frequency_list: List of input frequencies
    :type frequency_list: List
    :param amplitude_pair_list: list of amplitude pairs
    :type amplitude_pair_list: List
    :param phase_pair_list_rad_or_deg: list of phase pairs (can be rad or degree)
    :type phase_pair_list_rad_or_deg: List
    :param sort_out_factor: sort out factor [0...1]
    :type sort_out_factor: float
    :return: [frequency_list, amplitude_pair_list, phase_pair_list_rad_or_deg]
    :rtype: List
    """
    # Calculate geometric lengths
    amp_tot = np.sqrt(np.sum(np.array(amplitude_pair_list) ** 2, axis=0))
    # amp_tot = np.max(amplitude_pairs, axis=0)

    invalid_index = []
    for n, amplitude_pair in enumerate(amplitude_pair_list):
        if all(amplitude / amp_tot[i] < sort_out_factor for i, amplitude in enumerate(amplitude_pair)):
            invalid_index.append(n)

    phase_pair_list_rad_or_deg = np.delete(phase_pair_list_rad_or_deg, invalid_index, axis=0)
    amplitude_pair_list = np.delete(amplitude_pair_list, invalid_index, axis=0)
    frequency_list = np.delete(frequency_list, invalid_index)

    return [frequency_list, amplitude_pair_list, phase_pair_list_rad_or_deg]


# Reluctance Model [with calculation]

def calculate_cylinder_volume(cylinder_diameter: float, cylinder_height: float):
    """
    Calculate the volume of an ideal cylinder.

    This function is uses e.g. to calculate the volume of the inner core part.
    :param cylinder_height: height of cylinder
    :type cylinder_height: float
    :param cylinder_diameter: diameter of cylinder
    :type cylinder_diameter: float
    :returns: volume
    :rtype: float
    """
    return (cylinder_diameter / 2) ** 2 * np.pi * cylinder_height


def create_physical_group(dim, entities, name):
    """Greate a physical group, what is used inside ONELAB."""
    tag = gmsh.model.addPhysicalGroup(dim, entities)
    gmsh.model.setPhysicalName(dim, tag, name)

    return tag


def visualize_simulation_results(simulation_result_file_path: str, store_figure_file_path: str, show_plot=True) -> None:
    """Visualize the simulation results by a figure."""
    with open(simulation_result_file_path, "r") as fd:
        loaded_results_dict = json.loads(fd.read())

    # Initialize accumulators for cumulative losses and inductances
    cumulative_core_hysteresis = 0
    cumulative_core_eddy = 0
    cumulative_losses = []
    cumulative_inductances = []
    windings_labels = []

    for index, sweep in enumerate(loaded_results_dict["single_sweeps"]):
        freq = sweep['f']
        loss_core_eddy_current = sweep.get("core_eddy_losses", 0)
        loss_core_hysteresis = sweep.get("core_hyst_losses", 0)

        # Accumulate core losses
        cumulative_core_hysteresis += loss_core_hysteresis
        cumulative_core_eddy += loss_core_eddy_current

        # Plotting for each frequency
        fig, ax = plt.subplots()
        ax.bar(0, loss_core_hysteresis, width=0.35, label='Core Hysteresis Loss')
        ax.bar(0, loss_core_eddy_current, bottom=loss_core_hysteresis, width=0.35, label='Core Eddy Current Loss')

        for i in range(1, 4):
            winding_key = f"winding{i}"
            if winding_key in sweep:
                inductance = sweep[winding_key].get("flux_over_current", [0])[0]
                loss = sweep[winding_key].get("winding_losses", 0)

                if len(cumulative_losses) < i:
                    cumulative_losses.append(loss)
                    cumulative_inductances.append(inductance)
                    windings_labels.append(f"Winding {i}")
                else:
                    cumulative_losses[i - 1] += loss
                    cumulative_inductances[i - 1] += inductance

                # Plot for current frequency
                ax.bar(i, loss, width=0.35, label=f'{windings_labels[i - 1]} Loss at {freq} Hz')

        ax.set_ylabel('Losses in W')
        ax.set_title(f'Loss Distribution at {freq} Hz')
        ax.legend()
        plt.grid(True)

        # Save plot for the current frequency
        base_path, ext = os.path.splitext(store_figure_file_path)
        filename = f"{base_path}_{index}{ext}"
        plt.savefig(filename, bbox_inches="tight")
        plt.close(fig)

    # Plot cumulative results for core and windings
    fig, ax = plt.subplots()
    ax.bar(0, cumulative_core_hysteresis, width=0.35, label='Cumulative Core Hysteresis Loss')
    ax.bar(0, cumulative_core_eddy, bottom=cumulative_core_hysteresis, width=0.35, label='Cumulative Core Eddy Current Loss')

    for index, loss in enumerate(cumulative_losses):
        ax.bar(index + 1, loss, width=0.35, label=f'{windings_labels[index]} Cumulative Loss')

    ax.set_ylabel('total Losses in W')
    ax.set_title('Loss Distribution in Magnetic Components for all frequencies')
    ax.set_xticks(range(len(windings_labels) + 1))
    ax.set_xticklabels(['Core'] + windings_labels)
    ax.legend()
    plt.grid(True)
    base_path, ext = os.path.splitext(store_figure_file_path)
    cumulative_filename = f"{base_path}_total_freq{ext}"
    plt.savefig(cumulative_filename, bbox_inches="tight")

    if show_plot:
        plt.show()

    return loaded_results_dict

def point_is_in_rect(x, y, rect):
    """Check if a given x-y point is inside a rectangular field (e.g. inside a conductor)."""
    # x, y of the point
    # List of 4 points given as tuples with (x, y) in the order top-right, top-left, bottom-right, bottom-left

    # Return true if point is in rect
    if y < rect[0][1] and y > rect[3][1] and x > rect[0][0] and x < rect[1][0]:
        return True
    return False


def get_number_of_turns_of_winding(winding_windows, windings: List, winding_number: int):
    """Get the number of turns of a winding."""
    turns = 0
    for ww in winding_windows:
        for vww in ww.virtual_winding_windows:
            for index, winding in enumerate(windings):
                if winding.winding_number == winding_number:
                    # TODO: change index_turns right no. of winding numbers, right position in list and
                    #  length of list is needed
                    try:
                        turns += vww.turns[index]
                    except:
                        pass
    return turns


def cost_function_core(core_weight: float, core_type: str = "ferrite") -> float:
    """
    Calculate core material costs depending on material and weight.

    :param core_weight: core weight in kg
    :type core_weight: float
    :param core_type: core type. Can be "ferrite", "amorphous", "nanocristalline", "high_si_steel", "goes"
    :type core_type: str
    :return: costs of core in euro
    :rtype: float
    """
    cost_database = cost_material_database()
    sigma_core = cost_database["core_euro_per_kilogram"][core_type]

    return sigma_core * core_weight


def cost_function_winding(wire_weight_list: List[float], wire_type_list: List[str],
                          single_strand_cross_section_list: Union[List[float], None] = None):
    """
    Calculate single winding material and fabrication costs depending on winding-type and weight.

    Reference: Ralph Burkart and Johann W. Kolar: "Component Cost Models for Multi-Objective Optimizations
    of Switched-Mode Power Converters"

    :param wire_weight_list: winding weight in kg in list-form
    :type wire_weight_list: List[float]
    :param wire_type_list: winding type. Must fit to enum-names in ConductorType-Enum
    :type wire_type_list: List[str]
    :param single_strand_cross_section_list: single strand cross-section in list-form
    :type single_strand_cross_section_list: List[float]
    :return: winding cost of single winding
    :rtype: float
    """
    if single_strand_cross_section_list is None:
        single_strand_cross_section_list = []

    cost_database = cost_material_database()
    winding_cost_list = []

    for winding_count, winding_weight in enumerate(wire_weight_list):
        # material cost (per kilogram and per unit)
        sigma_material_winding_euro_per_kilogram = cost_database["winding_material_euro_per_kilogram"][
            wire_type_list[winding_count]]
        if sigma_material_winding_euro_per_kilogram == -1:
            # case for special litz wire calculation. Additional data is loaded from cost_database.
            sigma_material_winding_euro_per_kilogram = cost_database["winding_material_euro_per_kilogram_for_litz"][
                "sigma_numerator"] / (single_strand_cross_section_list[winding_count] * 1e6 + cost_database[
                    "winding_material_euro_per_kilogram_for_litz"]["sigma_denominator"])

        winding_material_euro_per_unit = cost_database["winding_material_euro_per_unit"][wire_type_list[winding_count]]

        winding_material_cost = sigma_material_winding_euro_per_kilogram * winding_weight + \
            winding_material_euro_per_unit

        # fabrication cost (per kilogram and per unit)
        sigma_fabrication_euro_per_kilogram = cost_database["winding_fabrication_euro_per_kilogram"][
            wire_type_list[winding_count]]
        fabrication_material_euro_per_unit = cost_database["winding_fabrication_euro_per_unit"][
            wire_type_list[winding_count]]

        winding_fabrication_cost = sigma_fabrication_euro_per_kilogram * winding_weight + fabrication_material_euro_per_unit

        winding_cost_list.append(winding_material_cost + winding_fabrication_cost)

    return winding_cost_list


def cost_function_total(core_weight: float, core_type: str, wire_weight_list: List[float], wire_type_list: List[str],
                        single_strand_cross_section_list: Union[None, List[float]] = None) -> float:
    """
    Calculate the total costs for an inductive element.

    This includes material costs for core and winding, fabrication costs for core and winding and manufacturer margin

    Reference: Ralph Burkart and Johann W. Kolar: "Component Cost Models for Multi-Objective Optimizations of
    Switched-Mode Power Converters"

    :param core_weight: core weight in kg
    :type core_weight: float
    :param core_type: core type. Can be "ferrite", "amorphous", "nanocristalline", "high_si_steel", "goes"
    :type core_type: str
    :param wire_weight_list: winding weight in kg
    :type wire_weight_list: float
    :param wire_type_list: winding type in list-form. Must fit to enum-names in ConductorType-Enum
    :type wire_type_list: List[str]
    :param single_strand_cross_section_list: single strand cross-section in list-form
    :type single_strand_cross_section_list: List[float]
    :return: total costs for inductive element
    :rtype: float
    """
    if single_strand_cross_section_list is None:
        single_strand_cross_section_list = []

    cost_database = cost_material_database()

    cost_core = cost_function_core(core_weight, core_type)

    cost_winding_list = cost_function_winding(wire_weight_list, wire_type_list, single_strand_cross_section_list)
    cost_winding = sum(cost_winding_list)

    total_cost_including_margin = 1 / (1 - cost_database["gross_margin"]) * (cost_core + cost_winding)

    return total_cost_including_margin


def find_result_log_file(result_log_folder: str, keyword_list: list, value_min_max: list):
    """
    Find a result log-file in a folder with many result-log files.

    Check a dictornary keyword list for matching a certain value (equel, greater equal, smaller equal).

    :param result_log_folder: filepath to result-log folder
    :type result_log_folder: str
    :param keyword_list: list with hirarchical keywords for dictionary structure, e.g. ["simulation_settings", "core", "core_inner_diameter"]
    :type keyword_list: list
    :param value_min_max: value to check for
    :type value_min_max: list

    :Example:

    Check for files with a core inner diameter smaller equal than 0.02 m.
    >>> import femmt as fmt
    >>> fmt.find_result_log_file("/home/filepath/fem_simulation_data", ["simulation_settings", "core",
    >>>     "core_inner_diameter"],[0.015, 0.02])
    """
    files_list = os.listdir(result_log_folder)

    value_min = value_min_max[0]
    value_max = value_min_max[1]

    for file in files_list:
        file_path = os.path.join(result_log_folder, file)
        with open(file_path, "r") as fd:
            full_data = json.loads(fd.read())

        if len(keyword_list) == 2:
            data_to_compare = full_data[keyword_list[0]][keyword_list[1]]
        elif len(keyword_list) == 3:
            data_to_compare = full_data[keyword_list[0]][keyword_list[1]][keyword_list[2]]
        elif len(keyword_list) == 4:
            data_to_compare = full_data[keyword_list[0]][keyword_list[1]][keyword_list[2]][keyword_list[3]]
        elif len(keyword_list) == 5:
            data_to_compare = full_data[keyword_list[0]][keyword_list[1]][keyword_list[2]][keyword_list[3]][
                keyword_list[4]]

        if value_min <= data_to_compare <= value_max:
            print(f"{value_min} <= {data_to_compare} <= {value_max} for file named {file}")


def wave_vector(f: float, complex_permeability: complex, complex_permittivity: complex, conductivity: float):
    """
    Calculate the wave-vector of a signal inside the core material with its material parameters.

    :param f: frequency
    :type f: float
    :param complex_permeability: complex permeability
    :type complex_permeability: complex
    :param complex_permittivity: complex permittivity
    :type complex_permittivity: complex
    :param conductivity: conductivity of the core material
    :type conductivity: float
    """
    omega = 2 * np.pi * f
    j = complex(0, 1)
    complex_equivalent_permittivity = complex_permittivity - j * conductivity / omega
    return omega * np.sqrt(complex_permeability * complex_equivalent_permittivity)


def axial_wavelength(f, complex_permeability, complex_permittivity, conductivity):
    """Calculate the axial wavelength for a given frequency."""
    k = wave_vector(f, complex_permeability, complex_permittivity, conductivity)
    return 2 * np.pi / k.real


def check_mqs_condition(radius: float, frequency: float, complex_permeability: float, complex_permittivity: float,
                        conductivity: float, relative_margin_to_first_resonance: float = 0.5, silent: bool = False):
    """
    Check if the condition for a magnetoquasistatic simulation is fulfilled.

    Calculates the ratio (core-diameter / wavelength) and includes a safety margin factor of 0.5.
    In case of ratio > 1, the simulated frequency is too high. A magnetoquasistatic simulation will not lead to good
    results. It is recommended to reduce the frequency or use a full-wave solver (not supported by FEMMT).

    :param radius: core radius
    :type radius: float
    :param frequency: frequency in Hz
    :type frequency: float
    :param complex_permeability: complex permeability
    :type complex_permeability: float
    :param complex_permittivity: complex permittivity
    :type complex_permittivity: float
    :param conductivity: core conductivity
    :type conductivity: float
    :param relative_margin_to_first_resonance: relative margin to the first resonance. Defaults to 0.5.
    :type relative_margin_to_first_resonance: float
    :param silent: True for no terminal output
    :type silent: bool
    """
    if frequency == 0:
        raise ValueError("check_mqs_condition() only works for frequencies != 0")

    axial_lambda = axial_wavelength(frequency, complex_permeability, complex_permittivity, conductivity)
    diameter_to_wavelength_ratio_of_first_resonance = 0.7655
    diameter_to_wavelength_ratio = 2 * radius / axial_lambda
    if diameter_to_wavelength_ratio > diameter_to_wavelength_ratio_of_first_resonance * relative_margin_to_first_resonance:
        # raise Warning(f"Resonance Ratio: {diameter_to_wavelength_ratio / diameter_to_wavelength_ratio_of_first_resonance} - "
        #               f"1 means 1st resonance - should be kept well below 1 to ensure MQS approach to be correct! ")
        if not silent:
            print(f"Resonance Ratio: {diameter_to_wavelength_ratio / diameter_to_wavelength_ratio_of_first_resonance}")


def create_open_circuit_excitation_sweep(I0, n, frequency):
    """Create a circuit excitation sweep with the other windings unloaded."""
    frequencies = [frequency] * n
    currents = [[0] * n for _ in range(n)]
    phases = [[180] * n for _ in range(n)]

    for x in range(0, n):
        for y in range(0, n):
            if x == y:
                currents[x][y] = I0
                phases[x][y] = 0

    return frequencies, currents, phases


def list_to_complex(complex_list: list):
    """
    Brings a list of two numbers (where first is real part, second is imaginary part) into a python specific complex number.

    :param complex_list:
    :type complex_list: list
    :return: complex number
    :rtype: complex
    """
    return complex(complex_list[0], complex_list[1])


def get_self_inductances_from_log(log: Dict) -> List:
    """
    Read the self-inductances from the result log file (dictionary).

    :param log: Result log dictionary
    :type log: Dict
    :return: self-inductances in a list
    :rtype: List
    """
    self_inductances = []
    for ol_index, open_loop_result in enumerate(log["single_sweeps"]):
        active_winding_name = f"winding{ol_index + 1}"
        self_inductances.append(list_to_complex(open_loop_result[active_winding_name]["flux_over_current"]))
    return self_inductances


def get_flux_linkages_from_log(log: Dict) -> List:
    """
    Read the flux-linkages from the result log file (dictionary).

    :param log: Result log dictionary
    :type log: Dict
    :return: flux-linkages in a list
    :rtype: List
    """
    flux_linkages = []
    for ol_index, open_loop_result in enumerate(log["single_sweeps"]):
        flux_linkages.append([])
        for winding_index in range(0, len(log["single_sweeps"])):
            flux_linkages[ol_index].append(list_to_complex(open_loop_result[f"winding{winding_index + 1}"]["flux"]))
    return flux_linkages


def get_coupling_matrix(flux_linkages: List) -> np.array:
    """
    Calculate the coupling factors from the given flux linkages.

    :param flux_linkages: flux-linkages
    :type flux_linkages: List
    :return: coupling-matrix in a matrix (np.array)
    :rtype: np.array
    """
    coupling_matrix = [[None] * len(flux_linkages) for _ in range(len(flux_linkages))]
    for self_index in range(0, len(flux_linkages)):
        for cross_index in range(0, len(flux_linkages)):
            coupling_matrix[cross_index][self_index] = flux_linkages[cross_index][self_index].real / \
                flux_linkages[self_index][self_index].real
    return coupling_matrix


def get_mean_coupling_factors(coupling_matrix: np.array):
    """
    Calculate the mean coupling factors from the coupling matrix.

    :param coupling_matrix: matrix with coupling factors between windings
    :type coupling_matrix: np.array
    """
    mean_coupling_factors = [[None] * len(coupling_matrix) for _ in range(len(coupling_matrix))]
    for self_index in range(0, len(coupling_matrix)):
        for cross_index in range(0, len(coupling_matrix)):
            mean_coupling_factors[cross_index][self_index] = (coupling_matrix[cross_index][self_index] * \
                                                              coupling_matrix[self_index][cross_index]) ** 0.5
    return mean_coupling_factors


def get_inductance_matrix(self_inductances, mean_coupling_factors, coupling_matrix):
    """Get the inductance matrix from self_inductances, mean_coupling_factors and the coupling_matrix."""
    inductance_matrix = [[None] * len(mean_coupling_factors) for _ in range(len(mean_coupling_factors))]
    for x in range(0, len(coupling_matrix)):
        for y in range(0, len(coupling_matrix)):
            inductance_matrix[x][y] = mean_coupling_factors[x][y] * (self_inductances[x] * self_inductances[y]) ** 0.5
    return inductance_matrix


def visualize_flux_linkages(flux_linkages: List, silent: bool) -> None:
    """
    Print the flux linkages to the terminal (or file-) output.

    :param flux_linkages: flux-linkages in a list
    :type flux_linkages: List
    :param silent: True for no output
    :type silent: bool
    """
    string_to_print = ""
    for x in range(0, len(flux_linkages)):
        for y in range(0, len(flux_linkages)):
            string_to_print += f"Phi_{x+1}{y+1} = {flux_linkages[x][y]}     Induced by I_{y+1} in Winding{x+1}\n"
    if not silent:
        print("\nFluxes: ")
        print(string_to_print)


def visualize_self_inductances(self_inductances: List | np.array, flux_linkages: list | np.array, silent: bool) -> None:
    """
    Print the self-inductances to the terminal (or file-) output.

    :param self_inductances: self-inductances in H in a list or numpy array
    :type self_inductances: List | np.array
    :param flux_linkages: flux linkages
    :type flux_linkages: list | np.array
    :param silent: True for no output
    :type silent: bool
    """
    string_to_print = ""
    for x in range(0, len(flux_linkages)):
        string_to_print += f"L_{x+1}_{x+1} = {self_inductances[x]}\n"
    if not silent:
        print("\n"
              "Self Inductances: ")
        print(string_to_print)


def visualize_self_resistances(self_inductances: List, flux_linkages: List, frequency: float, silent: bool) -> None:
    """
    Calculate and print the self resistances to the terminal (or file-) output.

    :param self_inductances: self-inductances in a list
    :type self_inductances: List
    :param flux_linkages: flux-linkage
    :type flux_linkages: List
    :param frequency: Frequency in Hz
    :type frequency: float
    :param silent: True for no output
    :type silent: bool
    """
    string_to_print = ""
    for x in range(0, len(flux_linkages)):
        string_to_print += f"Z_{x+1}_{x+1} = {self_inductances[x].imag*2*np.pi*frequency}\n"
    if not silent:
        print("\n"
              "Self Resistances: ")
        print(string_to_print)


def visualize_coupling_factors(coupling_matrix: np.array, flux_linkages: List, silent: bool):
    """
    Print the coupling factors to the terminal (or file-) output.

    :param coupling_matrix: matrix with coupling factors between the windings
    :type coupling_matrix: np.array
    :param flux_linkages: flux-linkages in a list
    :type flux_linkages: List
    :param silent: True for no output
    :type silent: bool
    """
    string_to_print = ""
    for x in range(0, len(flux_linkages)):
        for y in range(0, len(coupling_matrix)):
            string_to_print += f"K_{x + 1}{y + 1} = Phi_{x + 1}{y + 1} / Phi_{y + 1}{y + 1} = {coupling_matrix[x][y]}\n"
    if not silent:
        print("\n"
              "Coupling Factors: ")
        print(string_to_print)


def visualize_mean_coupling_factors(mean_coupling_factors: List, silent: bool):
    """
    Print the mean coupling factors to the terminal (or file-) output.

    :param mean_coupling_factors: mean_coupling_factors in a list
    :type mean_coupling_factors: List
    :param silent: True for no output
    :type silent: bool
    """
    string_to_print = ""
    for x in range(0, len(mean_coupling_factors)):
        for y in range(0, len(mean_coupling_factors)):
            string_to_print += f"k_{x + 1}{y + 1} = Sqrt(K_{x + 1}{y + 1} * K_{y + 1}{x + 1}) = M_{x + 1}{y + 1} " \
                               f"/ Sqrt(L_{x + 1}_{x + 1} * L_{y + 1}_{y + 1}) = {mean_coupling_factors[x][y]}\n"
        if not silent:
            print("\n"
                  "Mean Coupling Factors: ")
            print(string_to_print)


def visualize_mean_mutual_inductances(inductance_matrix: np.array, silent: bool):
    """
    Print the mean mutal inductances to the terminal (or file-) output.

    :param inductance_matrix: inductance matrix
    :type inductance_matrix: np.array
    :param silent: True for no output
    :type silent: bool

    e.g.  M_12 = M_21 = k_12 * (L_11 * L_22) ** 0.5
    """
    string_to_print = ""
    for x in range(0, len(inductance_matrix)):
        for y in range(0, len(inductance_matrix)):
            if x == y:
                pass
            else:
                string_to_print += f"M_{x + 1}{y + 1} = {inductance_matrix[x][y].real}\n"
    if not silent:
        print("\n"
              "Mean Mutual Inductances: ")
        print(string_to_print)


def visualize_mutual_inductances(self_inductances: List, coupling_factors: List, silent: bool):
    """
    Print the mutal inductances to the terminal (or file-) output.

    :param self_inductances: Matrix with self inductances
    :type self_inductances: List
    :param coupling_factors: Matrix with coupling factors
    :type coupling_factors: List
    :param silent: True for no output
    :type silent: bool

    e.g. M_12 = L_11 * K_21  !=   M_21 = L_22 * K_12   (ideally, they are the same)
    """
    string_to_print = ""
    for x in range(0, len(coupling_factors)):
        for y in range(0, len(coupling_factors)):
            if x == y:
                pass
            else:
                string_to_print += f"M_{x + 1}{y + 1} = {self_inductances[y].real * coupling_factors[x][y]}\n"
    if not silent:
        print("\n"
              "Mutual Inductances: ")
        print(string_to_print)


def visualize_inductance_matrix_coefficients(inductance_matrix: np.array, silent: bool):
    """Visualize the inductance matrix coefficients in the terminal.

    e.g. M_12 = L_11 * K_21  !=   M_21 = L_22 * K_12   (ideally, they are the same)

    :param inductance_matrix: inductance matrix of transformer
    :type inductance_matrix: np.array
    :param silent: False to show the terminal output
    :type silent: bool
    """
    string_to_print = ""
    for x in range(0, len(inductance_matrix)):
        for y in range(0, len(inductance_matrix)):
            if x == y:
                string_to_print += f"L_{x + 1}{y + 1} = {inductance_matrix[x][y].real}\n"
            else:
                string_to_print += f"M_{x + 1}{y + 1} = {inductance_matrix[x][y].real}\n"
    if not silent:
        print("\n"
              "Inductance Matrix Coefficients: ")
        print(string_to_print)


def visualize_inductance_matrix(inductance_matrix, silent: bool):
    """Visualize the inductance matrix in the terminal.

    e.g. M_12 = L_11 * K_21  !=   M_21 = L_22 * K_12   (ideally, they are the same)
    """
    string_to_print = ""
    for x in range(0, len(inductance_matrix)):
        for y in range(0, len(inductance_matrix)):
            string_to_print += f"{np.round(inductance_matrix[x][y].real, 12)} "
        string_to_print += "\n"

    if not silent:
        print("\n"
              "Inductance Matrix: ")
        print(string_to_print)

def calculate_quadrature_integral(time_steps: List[float], data: List[float]) -> float:
    """
    Calculate the integral of given data over specific time steps using the quadrature method.

    :param time_steps: List of time steps.
    :type time_steps: List[float]
    :param data: List of data corresponding to each timestep.
    :type data: List[float]
    :return: The calculated integral.
    :rtype: float
    """
    func = lambda x: np.interp(x, time_steps, data)
    return quadrature(func, time_steps[0], time_steps[-1])[0]

def calculate_squared_quadrature_integral(time_steps: List[float], data: List[float]) -> float:
    """
    Calculate the integral of squared given data over specific time steps using the quadrature method..

    :param time_steps: List of time steps.
    :type time_steps: List[float]
    :param data: List of data corresponding to each timestep.
    :type data: List[float]
    :return: The calculated integral.
    :rtype: float
    """
    func = lambda x: np.interp(x, time_steps, data) ** 2
    return quadrature(func, time_steps[0], time_steps[-1])[0]

def calculate_average(integral: float, time_steps: List[float]) -> float:
    """
    Compute the average in general.

    :param integral: The integral value.
    :type integral: float
    :param time_steps: List of time steps.
    :type time_steps: List[float]

    Returns:
    :return: The calculated average.
    :rtype: float.
    """
    total_time = time_steps[-1] - time_steps[0]
    if total_time == 0:
        raise ValueError("Total time cannot be zero.")
    return integral / total_time

def calculate_rms(squared_integral: float, time_steps: List[float]) -> float:
    """
    Compute the RMS.

    :param squared_integral: The integral value.
    :type squared_integral: float
    :param time_steps: List of time steps.
    :type time_steps: List[float]

    Returns:
    :return: The calculated average.
    :rtype: float.
    """
    total_time = time_steps[-1] - time_steps[0]
    if total_time == 0:
        raise ValueError("Total time cannot be zero.")
    mean_square = squared_integral / total_time  # Calculate the mean of the square of the data
    return np.sqrt(mean_square)  # Take the square root to get RMS value


def convert_air_gap_corner_points_to_center_and_distance(corner_points: list) -> list:
    """
    Convert the list-defined air_gap_corner_points from a "two_d_axi" object to center points and lengths as to separate lists.

    :param corner_points: in usage of magnetic component -> "self.two_d_axi.p_air_gaps.tolist()"
    :type corner_points: list
    :returns: centers and heights of the air gaps
    :rtype: list
    """
    centers = []
    heights = []
    # 4 corner points make up one air gap
    for n_air_gap in range(0, int(len(corner_points)/4)):
        width = corner_points[n_air_gap*4 + 1][0] - corner_points[n_air_gap*4 + 0][0]
        height = corner_points[n_air_gap * 4 + 2][1] - corner_points[n_air_gap * 4 + 0][1]
        heights.append(height)
        centers.append(
            [
                corner_points[n_air_gap*4 + 0][0] + width/2,  # x-coordinate of air gap's center
                corner_points[n_air_gap*4 + 2][1] - height/2  # y-coordinate of air gap's center
            ]
        )
    return centers, heights


if __name__ == '__main__':
    pass
