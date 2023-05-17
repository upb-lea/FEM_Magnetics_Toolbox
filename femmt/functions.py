# Python standard libraries
import json
import random
import string
import pkg_resources
import subprocess
import sys
import os
import time
import warnings
from typing import Union, List, Tuple, Dict

# Third parry libraries
import gmsh
import pandas as pd
from matplotlib import pyplot as plt
import numpy.typing as npt
import numpy as np

# Local libraries
from femmt.constants import *
from femmt.enumerations import ConductorType

# Needed for femmt_print
silent = False

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
                    "winding": ["orange", "brown", "yellow", "green", "red", "black", "grey", "blue", "orange", "purple"],
                    "insulation": "blue",
                    "potting_inner": "grey",
                    "potting_outer": "grey",
                }



colors_ba_jonas = {"blue": (28, 113, 216),
                        'red': (213, 6, 6),
                        "green":  (6, 213, 6),
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
    Returns core geometry for defined core structures

    All dimensions are nominal dimensions without consideration of tolerances.

    For PQ core sizes (e.g. PQ 40/40), it has been found out that
    core_dimension_x / core_dimension_y = 1.45, the error over all available shapes is maximum 7% (compared to datasheet value)
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
        "window_w": 5.85e-3,
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
        "window_w": (22.5-12)/2*1e-3,
        "core_dimension_x": 26.5e-3,
        "core_dimension_y": 19.0e-3,
    }
    core_dict["PQ 32/20"] = {
        "core_h": 20.5e-3,
        "core_inner_diameter": 13.45e-3,
        "window_h": 11.5e-3,
        "window_w": (27.5-13.45)/2 * 1e-3,
        "core_dimension_x": 32.0e-3,
        "core_dimension_y": 22.0e-3,
    }
    core_dict["PQ 32/30"] = {
        "core_h": 30.35e-3,
        "core_inner_diameter": 13.45e-3,
        "window_h": 21.3e-3,
        "window_w": (27.5-13.45)/2 * 1e-3,
        "core_dimension_x": 32.0e-3,
        "core_dimension_y": 22.0e-3,
    }
    core_dict["PQ 35/35"] = {
        "core_h": 34.8e-3,
        "core_inner_diameter": 14.35e-3,
        "window_h": 25e-3,
        "window_w": (32-14.35)/2 * 1e-3,
        "core_dimension_x": 35.1e-3,
        "core_dimension_y": 26.0e-3,
    }
    core_dict["PQ 40/30"] = {
        "core_h": 30.3e-3,
        "core_inner_diameter": 14.9e-3,
        "window_h": 20e-3,
        "window_w": (37-14.9)/2 * 1e-3,
        "core_dimension_x": 40.5e-3,
        "core_dimension_y": 28.0e-3,
    }
    core_dict["PQ 40/40"] = {
        "core_h": 39.8e-3,
        "core_inner_diameter": 14.9e-3,
        "window_h": 29.5e-3,
        "window_w": (37-14.9)/2 * 1e-3,
        "core_dimension_x": 40.5e-3,
        "core_dimension_y": 28.0e-3,
    }
    core_dict["PQ 50/40"] = {
        "core_h": 40e-3,
        "core_inner_diameter": 20e-3,
        "window_h": 26.1e-3,
        "window_w": (44-20)/2 * 1e-3,
        "core_dimension_x": 50.0e-3,
        "core_dimension_y": 32.0e-3,
    }
    core_dict["PQ 50/50"] = {
        "core_h": 50e-3,
        "core_inner_diameter": 20e-3,
        "window_h": 36.1e-3,
        "window_w": (44-20)/2 * 1e-3,
        "core_dimension_x": 50.0e-3,
        "core_dimension_y": 32.0e-3,
    }
    core_dict["PQ 65/60"] = {
        "core_h": 60e-3,
        "core_inner_diameter": 26e-3,
        "window_h": 42e-3,
        "window_w": (65-26)/2 * 1e-3,
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
        "window_w":  (88-43)/2 * 1e-3,
    }
    core_dict["PM 50/39"] = {
        "core_h": 39 * 1e-3,
        "core_inner_diameter": pm_core_inner_diameter_calculator(20e-3, 5.4e-3),
        "window_h": 26.4 * 1e-3,
        "window_w": (39-20)/2 * 1e-3,
    }
    core_dict["PM 62/49"] = {
        "core_h": 49 * 1e-3,
        "core_inner_diameter": pm_core_inner_diameter_calculator(25.5e-3, 5.4e-3),
        "window_h": 33.4 * 1e-3,
        "window_w": (48.8-25.5)/2 * 1e-3,
    }
    core_dict["PM 74/59"] = {
        "core_h": 59 * 1e-3,
        "core_inner_diameter": pm_core_inner_diameter_calculator(29.5e-3, 5.4e-3),
        "window_h": 40.7e-3,
        "window_w": (57.5-29.5)/2 * 1e-3,
    }
    core_dict["PM 87/70"] = {
        "core_h": 70 * 1e-3,
        "core_inner_diameter": pm_core_inner_diameter_calculator(31.7e-3, 8.5e-3),
        "window_h": 48 * 1e-3,
        "window_w": (67.1-31.7)/2 * 1e-3,
    }
    return core_dict


def litz_database() -> Dict:
    """
    Returns litz parameters for defined litz wires.

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


    return litz_dict


def wire_material_database() -> Dict:
    """
    Returns wire materials e.g. copper, aluminium in a dictionary

    :return: Dict with materials and conductivity
    :rtype: Dict
    """
    wire_material = {}

    wire_material["Copper"] = {
        "sigma": 5.8e7,
        "thermal_conductivity": 400,
        "volumetric_mass_density": 8920,
    }

    wire_material["Aluminium"] = {
        "sigma": 3.7e7,
        "thermal_conductivity": 235,
        "volumetric_mass_density": 2699,
    }

    return wire_material

def create_folders(*args) -> None:
    """Creates folder for every given folder path (if it does not exist).
    """
    for folder in list(args):
        if not os.path.exists(folder):
            os.mkdir(folder)

def cost_material_database() -> Dict:
    """
    Returns costs for core and winding.
    This is split in material and fabrication costs.
    Both, material and fabrication costs have a euro_per_kilogram and a euro_per_unit (fix costs) price.

    Source: R. Burkart and J. Kolar 'Component Cost Models for Multi-Objective Optimizations of Switched-Mode Power Converter'
    2013.

    These are outdated prices (year 2013). Update needed in future.
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
    Calculates the effective inner core diameter without the hole
    Often used in PM-cores

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
    """
    Windows users only.
    Installs femm-software pip package in case of running on Windows machine

    :return: None

    """
    required = {'pyfemm'}
    installed = {pkg.key for pkg in pkg_resources.working_set}
    missing = required - installed

    if missing:
        print("Missing 'pyfemm' installation.")
        print("Installing 'pyfemm' ...")
        python = sys.executable
        subprocess.check_call([python, '-m', 'pip', 'install', *missing], stdout=subprocess.DEVNULL)
        print("'pyfemm' is now installed!")

def inner_points(a, b, input_points):
    """
    Returns the input points that have a common coordinate as the two
    interval borders a and b

    Use case: Air gap generation. Helps to generate the core structure form class AirGap.

    :param a: first point
    :param b: second point
    :param input_points:

    :return:

    """
    [min_point, max_point] = [None, None]
    output = input_points
    dim = None
    # Find equal dimension
    for i in range(0, 3):
        if a[i] == b[i] and a[i] != 0:
            dim = i
    if dim is None:
        raise Exception("Given points do not have a common dimension")
    n = 0
    while n < output.shape[0]:
        if a[dim] != output[n, dim]:
            output = np.delete(output, n, 0)
        else:
            n += 1
    if output.shape[0] == 0:
        raise Exception("Not implemented Error: No air gaps between interval borders")
    if output.shape[0] % 2 == 1:
        raise Exception("Odd number of input points")
    if dim == 2:
        raise Exception("Not implemented Error: Only 2D is implemented")
    dim2 = (dim+1) % 2
    if output.shape[0] >= 2:
        argmax = np.argmax(output[:, dim2])
        output = np.delete(output, argmax, 0)
        argmin = np.argmin(output[:, dim2])
        output = np.delete(output, argmin, 0)
        #if output.shape[0] == 0:
            #print("Only one air gap in this leg. No island needed.")
    return output


def min_max_inner_points(a, b, input_points):
    """
    Returns the input points that have a common coordinate and
    the minimum distance from the interval borders.

    :param a: first point
    :param b: second point
    :param input_points:

    :return:

    """

    [min_point, max_point] = [None, None]
    buffer = input_points
    dim = None
    # Find equal dimension
    for i in range(0, 3):
        if a[i] == b[i] and a[i] != 0:
            dim = i
    if dim is None:
        raise Exception("Given points do not have a common dimension")
    n = 0
    while n < buffer.shape[0]:
        if a[dim] != buffer[n, dim]:
            buffer = np.delete(buffer, n, 0)
        else:
            n += 1
    if buffer.shape[0] == 0:
        femmt_print("No air gaps between interval borders")
    if buffer.shape[0] % 2 == 1:
        raise Exception("Odd number of input points")
    if dim == 2:
        raise Exception("Not implemented Error: Only 2D is implemented")
    dim2 = (dim+1) % 2
    if buffer.shape[0] >= 2:
        argmax = np.argmax(buffer[:, 1])
        max_point = buffer[argmax]
        argmin = np.argmin(buffer[:, 1])
        min_point = buffer[argmin]
    return [min_point, max_point]
    

def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


def litz_calculate_number_strands(n_layers: int) -> int:
    """
    Returns the number of strands in a hexagonal litz winding with a
    specified number of layers (n_layers). CAUTION: Zero number of
    layers corresponds to a single strand.

    :param n_layers: number of litz_layers
    :type n_layers: int

    :return: number of strands in a litz wire
    :rtype: int

    """
    return 3 * (n_layers + 1) ** 2 - 3 * (n_layers + 1) + 1


def litz_calculate_number_layers(n_strands: int) -> int:
    """
    Returns the number of layers in a hexagonal litz winding with a
    specified number of strands (n_strands).

    .. note:: Zero number of layers corresponds to a single strand.

    :param n_strands: Number of strands in a litz
    :type n_strands: int

    :return: number of layers for a litz
    :rtype: int

    """
    return np.sqrt(0.25+(n_strands-1)/3)-0.5


def fft(period_vector_t_i: npt.ArrayLike, sample_factor: int = 1000, plot: str = 'no', mode: str = 'rad',
        f0: Union[float, None] = None, title: str = 'ffT', filter_type: str = 'factor',
        filter_value_factor: float = 0.01, filter_value_harmonic: int = 100,
        figure_size: Tuple=None, figure_directory: str=None) -> npt.NDArray[list]:
    """
    A fft for an input signal. Input signal is in vector format and should include one period.

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
        if plot != 'no' and plot != False:
            print("Input is list, convert to np.array()")
        period_vector_t_i = np.array(period_vector_t_i)

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
        for count, value in enumerate(x_mag_corrected):
            if x_mag_corrected[count] > filter_value_factor * max(abs(i)):
                f_out.append(f_corrected[count])
                x_out.append(x_mag_corrected[count])
                phi_rad_out.append(phi_rad_corrected[count])
    elif filter_type.lower() == 'harmonic':
        for count, value in enumerate(x_mag_corrected):
            if count < filter_value_harmonic:
                f_out.append(f_corrected[count])
                x_out.append(x_mag_corrected[count])
                phi_rad_out.append(phi_rad_corrected[count])
    elif filter_type.lower() == 'disabled':
        f_out = f_corrected
        x_out = x_mag_corrected
        phi_rad_out = phi_rad_corrected
    else:
        raise ValueError(f"filter_type '{filter_value_harmonic}' not available: Must be 'factor','harmonic' or 'disabled ")


    if plot != 'no' and plot != False:
        print(f"{title = }")
        print(f"{t[-1] = }")
        print(f"{f0 = }")
        print(f"{Fs = }")
        print(f"{sample_factor = }")
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


def plot_fourier_coefficients(frequency_list, amplitude_list, phi_rad_list, sample_factor: int = 1000, figure_directory: str = None):
    time_period = 1 / min(frequency_list)

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

    #plt.show()


def compare_fft_list(input_data_list: list, sample_factor: float = 1000,  mode: str = 'rad', f0: Union[float,None] = None) -> None:
    """
    generate fft curves from input curves and compare them to each other

    :Minimal Example:

    >>> example_waveform = np.array([[0, 1.34, 3.14, 4.48, 6.28],[-175.69, 103.47, 175.69, -103.47,-175.69]])
    >>> example_waveform_2 = np.array([[0, 0.55, 3.14, 3.69, 6.28],[-138.37, 257.58, 138.37, -257.58, -138.37]])
    >>> compare_fft_list([example_waveform, example_waveform_2], mode='rad', f0=25000)

    :param input_data_list: list of fft-compatible numpy-arrays [element, element, ... ], each element format like [[time-vector[,[current-vector]]. One period only
    :param mode: 'rad'[default]: full period is 2*pi, 'deg': full period is 360°, 'time': time domain.
    :type mode: str
    :param f0: fundamental frequency. Needs to be set in 'rad'- or 'deg'-mode
    :type f0: float
    :param sample_factor: samle factor, defaults to 1000
    :type sample_factor: float

    :return: plot

    """
    out = []
    for count, value in enumerate(input_data_list):
        out.append([fft(input_data_list[count], sample_factor=sample_factor, plot='no', mode=mode, f0=f0)])

    fig, axs = plt.subplots(2, len(input_data_list), sharey=True)
    for count, value in enumerate(input_data_list):
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
    Stores a numpy array in a given directory

    :param dir_path: directory path
    :type dir_path: str
    :param file_name: file name
    :type file_name: str
    :param numpy_data: numpy array
    :type numpy_data:
    :return: None
    :rtype: None
    """
    if not os.path.isdir(dir_path):
        os.mkdir(dir_path)
    np.save(dir_path + "/" + file_name, numpy_data)


def data_logging(sim_choice):
    """
    !!! not finally implemented !!!

    This method shall do the saving and loading of results! with date and time

    :return:

    """
    frequencies = None
    # Data Logging with date and time
    datetime = time.strftime("%Y%m%d_%H%M%S")

    target_femm = 'Pv_FEMM_' + datetime + '.json'
    target_femmt = 'Pv_FEMMT_' + datetime + '.json'

    # Either read old data or create
    if sim_choice != 'show':
      # pseudo 2D dataframe: ['strand number, strand radius'], 'frequency'  | const. litz radius --> ff
      if not os.path.isfile(target_femm):
          df_pv_femm = pd.DataFrame([], index=frequencies, columns=[])
          df_pv_femmt = pd.DataFrame([], index=frequencies, columns=[])
      else:
          # Read Loss Data
          df_pv_femm = pd.read_json(target_femm)
          df_pv_femmt = pd.read_json(target_femmt)
      # print(df_pv_femm, df_pv_femmt)


def get_dicts_with_keys_and_values(data, **kwargs):
    """
    Returns a list of dictionaries out of a list of dictionaries which contains pairs of the given key(s) and value(s).
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


def get_dict_with_unique_keys(data, *keys):
    """
    Returns a dictionary out of a list of dictionaries which contains the given key(s).
    :param data: list of dicts
    :param keys: keys in dicts
    :return:
    """
    invalid_index = []
    for n, dictionary in enumerate(data):
        for key in keys:
            if not key in dictionary:
                invalid_index.append(n)
                break
    valid_data = np.delete(data, invalid_index)
    if len(valid_data) != 1:
        warnings.warn("Keyword(s) not unique!")
    # Only one dictionary shall survive --> choose element 0
    return valid_data[0]


def find_common_frequencies(frequency_list_1: List, amplitude_list_1: List, phase_list_1_rad_or_deg: List,
                            frequency_list_2: List, amplitude_list_2: List, phase_list_2_rad_or_deg: List) -> List:
    """
    Finds common frequencies and returns a list of intersections.

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
    >>> common_f, common_a, common_phase = fmt.find_common_frequencies(frequency_1, amplitude_1, phase_1, frequency_2, amplitude_2, phase_2)
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
    print(f"{common_frequency_list = }")
    common_frequency_list.sort()
    print(f"{common_frequency_list = }")

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
        warnings.warn("Either a list or a np.ndarray must be provided!")

    current_pair_list = list(map(list, zip(common_amplitude_list_1, common_amplitude_list_2)))
    phase_pair_list = list(map(list, zip(common_phase_list_1, common_phases_list_2)))

    return [common_frequency_list, current_pair_list, phase_pair_list]


def sort_out_small_harmonics(frequency_list: List, amplitude_pair_list: List,
                             phase_pair_list_rad_or_deg: List, sort_out_factor: float) -> List:
    """
    Sout out small harmonics.
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
    Calculates the volume of an ideal cylinder. This function is uses e.g. to calculate the volume
    of the inner core part.
    :param cylinder_height: height of cylinder
    :type cylinder_height: float
    :param cylinder_diameter: diameter of cylinder
    :type cylinder_diameter: float
    :returns: volume
    :rtype: float
    """
    return (cylinder_diameter / 2) ** 2 * np.pi * cylinder_height



def create_physical_group(dim, entities, name):
    tag = gmsh.model.addPhysicalGroup(dim, entities)
    gmsh.model.setPhysicalName(dim, tag, name)

    return tag


def visualize_simulation_results(simulation_result_file_path: str, store_figure_file_path: str, show_plot = True) -> None:
    with open(simulation_result_file_path, "r") as fd:
        loaded_results_dict = json.loads(fd.read())

    inductance = loaded_results_dict["single_sweeps"][0]["winding1"]["self_inductance"][0]
    loss_core_eddy_current = loaded_results_dict["total_losses"]["eddy_core"]
    loss_core_hysteresis = loaded_results_dict["total_losses"]["hyst_core_fundamental_freq"]
    loss_winding_1 = loaded_results_dict["total_losses"]["winding1"]["total"]

    femmt_print(inductance)
    femmt_print(loss_core_eddy_current)
    femmt_print(loss_core_hysteresis)
    femmt_print(loss_winding_1)

    bar_width = 0.35
    plt.bar(0, loss_core_hysteresis, width=bar_width)
    plt.bar(0, loss_core_eddy_current, bottom= loss_core_hysteresis, width=bar_width)
    plt.bar(1, loss_winding_1, width=bar_width)
    plt.legend(['Hysteresis loss', 'Eddy current loss', 'Winding loss'])
    plt.ylabel('Losses in W')
    plt.xticks([0, 1, 2], ['Core', 'Winding 1', 'Winding 2'])
    plt.grid()
    plt.savefig(store_figure_file_path, bbox_inches="tight")

    if show_plot:
        plt.show()

    return loaded_results_dict

def point_is_in_rect(x, y, rect):
    # x, y of the point
    # List of 4 points given as tuples with (x, y) in the order top-right, top-left, bottom-right, bottom-left

    # Return true if point is in rect
    if y < rect[0][1] and y > rect[3][1] and x > rect[0][0] and x < rect[1][0]:
        return True
    
    return False

def set_silent_status(s: bool) -> None:
    """
    Variable to store the silent status to show terminal outputs or not.
    Using silent-mode increases speed of simulation significantly.
    :param s: status for silent mode True/False
    :type s: bool
    :return: None
    :rtype: None
    """
    global silent
    silent = s

def femmt_print(text: str, end:str='\n') -> None:
    """
    Print-function for femmt. Uses a silent-Flag as global variable, to show terminal outputs or not.
    Using silent-mode increases speed of simulation significantly.

    :param text: text to print
    :type text: str
    :param end: command text ends with, defaults to '\n'
    :type: end: str
    :return: None
    :rtype: None
    """
    if not silent:
        print(text, end)

def cost_function_core(core_weight: float, core_type: str = "ferrite") -> float:
    """
    Calculates core material costs depending on material and weight.

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


def cost_function_winding(wire_weight_list: List[float], wire_type_list: List[str], single_strand_cross_section_list: Union[List[float], None] = None):
    """
    Calculates single winding material and fabrication costs depending on winding-type and weight

    Reference: Ralph Burkart and Johann W. Kolar: "Component Cost Models for Multi-Objective Optimizations of Switched-Mode Power Converters"

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
        sigma_material_winding_euro_per_kilogram = cost_database["winding_material_euro_per_kilogram"][wire_type_list[winding_count]]
        if sigma_material_winding_euro_per_kilogram == -1:
            # case for special litz wire calculation. Additional data is loaded from cost_database.
            sigma_material_winding_euro_per_kilogram = cost_database["winding_material_euro_per_kilogram_for_litz"]["sigma_numerator"] / (single_strand_cross_section_list[winding_count] * 1e6 + cost_database["winding_material_euro_per_kilogram_for_litz"]["sigma_denominator"])

        winding_material_euro_per_unit = cost_database["winding_material_euro_per_unit"][wire_type_list[winding_count]]

        winding_material_cost = sigma_material_winding_euro_per_kilogram * winding_weight + winding_material_euro_per_unit

        # fabrication cost (per kilogram and per unit)
        sigma_fabrication_euro_per_kilogram = cost_database["winding_fabrication_euro_per_kilogram"][wire_type_list[winding_count]]
        fabrication_material_euro_per_unit = cost_database["winding_fabrication_euro_per_unit"][wire_type_list[winding_count]]

        winding_fabrication_cost = sigma_fabrication_euro_per_kilogram * winding_weight + fabrication_material_euro_per_unit

        winding_cost_list.append(winding_material_cost + winding_fabrication_cost)

    return winding_cost_list



def cost_function_total(core_weight: float, core_type: str, wire_weight_list: List[float], wire_type_list: List[str],
                        single_strand_cross_section_list: Union[None, List[float]] = None) -> float:
    """
    Calculates the total costs for an inductive element.
    This includes material costs for core and winding, fabrication costs for core and winding and manufacturer margin

    Reference: Ralph Burkart and Johann W. Kolar: "Component Cost Models for Multi-Objective Optimizations of Switched-Mode Power Converters"

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
    find a result log-file in a folder with many result-log files.
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
    >>> fmt.find_result_log_file("/home/filepath/fem_simulation_data", ["simulation_settings", "core", "core_inner_diameter"],[0.015, 0.02])

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
            data_to_compare = full_data[keyword_list[0]][keyword_list[1]][keyword_list[2]][keyword_list[3]][keyword_list[4]]

        if value_min <= data_to_compare <= value_max:
            print(f"{value_min} <= {data_to_compare} <= {value_max} for file named {file}")


def wave_vector(f, complex_permeability, complex_permittivity, conductivity):
    omega = 2*np.pi*f
    j = complex(0, 1)
    complex_equivalent_permittivity = complex_permittivity - j * conductivity / omega
    return omega * np.sqrt(complex_permeability * complex_equivalent_permittivity)


def axial_wavelength(f, complex_permeability, complex_permittivity, conductivity):
    k = wave_vector(f, complex_permeability, complex_permittivity, conductivity)
    return 2 * np.pi / k.real


def check_mqs_condition(radius, f, complex_permeability, complex_permittivity, conductivity, relative_margin_to_first_resonance=0.5):
    axial_lambda = axial_wavelength(f, complex_permeability, complex_permittivity, conductivity)
    diameter_to_wavelength_ratio_of_first_resonance = 0.7655
    diameter_to_wavelength_ratio = 2 * radius / axial_lambda
    if diameter_to_wavelength_ratio > diameter_to_wavelength_ratio_of_first_resonance * relative_margin_to_first_resonance:
        # raise Warning(f"Resonance Ratio: {diameter_to_wavelength_ratio / diameter_to_wavelength_ratio_of_first_resonance} - "
        #               f"1 means 1st resonance - should be kept well below 1 to ensure MQS approach to be correct! ")
        femmt_print(f"Resonance Ratio: {diameter_to_wavelength_ratio / diameter_to_wavelength_ratio_of_first_resonance}")

if __name__ == '__main__':
    pass
