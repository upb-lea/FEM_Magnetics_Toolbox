# Python standard libraries
import numpy.typing as npt
import numpy as np
import json
import random
import string
import pkg_resources
import subprocess
import sys
import os
import pandas as pd
import time
import warnings
from typing import Union, List, Tuple, Dict
from matplotlib import pyplot as plt
from femmt.Enumerations import ConductorType

# Third parry libraries
import gmsh

# Local libraries
from femmt.Enumerations import *

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
                    "winding": ["orange", "brown", "yellow"],
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

    :return: Dict including core_h, core_inner_diameter, window_h, window_w
    :rtype: Dict
    """
    core_dict = {}

    # -----------------------
    # PQ Cores
    # -----------------------

    core_dict["PQ 16/11.5"] = {
        "core_h": 11.6e-3,
        "core_inner_diameter": 7e-3,
        "window_h": 6.7e-3,
        "window_w": 3.7e-3,
    }
    core_dict["PQ 20/16"] = {
        "core_h": 16.2e-3,
        "core_inner_diameter": 8.8e-3,
        "window_h": 10.3e-3,
        "window_w": 5.85e-3,
    }
    core_dict["PQ 20/20"] = {
        "core_h": 20.2e-3,
        "core_inner_diameter": 8.8e-3,
        "window_h": 14.3e-3,
        "window_w": 4.6e-3,
    }
    core_dict["PQ 26/20"] = {
        "core_h": 20.16e-3,
        "core_inner_diameter": 12e-3,
        "window_h": 11.5e-3,
        "window_w": 5.25e-3,
    }
    core_dict["PQ 26/25"] = {
        "core_h": 24.76e-3,
        "core_inner_diameter": 12e-3,
        "window_h": 16.1e-3,
        "window_w": (22.5-12)/2*1e-3,
    }
    core_dict["PQ 32/20"] = {
        "core_h": 20.5e-3,
        "core_inner_diameter": 13.45e-3,
        "window_h": 11.5e-3,
        "window_w": (27.5-13.45)/2 * 1e-3,
    }
    core_dict["PQ 32/30"] = {
        "core_h": 30.35e-3,
        "core_inner_diameter": 13.45e-3,
        "window_h": 21.3e-3,
        "window_w": (27.5-13.45)/2 * 1e-3,
    }
    core_dict["PQ 35/35"] = {
        "core_h": 34.8e-3,
        "core_inner_diameter": 14.35e-3,
        "window_h": 25e-3,
        "window_w": (32-14.35)/2 * 1e-3,
    }
    core_dict["PQ 40/30"] = {
        "core_h": 30.3e-3,
        "core_inner_diameter": 14.9e-3,
        "window_h": 20e-3,
        "window_w": (37-14.9)/2 * 1e-3,
    }
    core_dict["PQ 40/40"] = {
        "core_h": 39.8e-3,
        "core_inner_diameter": 14.9e-3,
        "window_h": 29.5e-3,
        "window_w": (37-14.9)/2 * 1e-3,
    }
    core_dict["PQ 50/40"] = {
        "core_h": 40e-3,
        "core_inner_diameter": 20e-3,
        "window_h": 26.1e-3,
        "window_w": (44-20)/2 * 1e-3,
    }
    core_dict["PQ 50/50"] = {
        "core_h": 50e-3,
        "core_inner_diameter": 20e-3,
        "window_h": 36.1e-3,
        "window_w": (44-20)/2 * 1e-3,
    }
    core_dict["PQ 65/60"] = {
        "core_h": 60e-3,
        "core_inner_diameter": 26e-3,
        "window_h": 21e-3,
        "window_w": (65-26)/2 * 1e-3,
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

    litz_dict["1.5x105x0.1"] = {"implicit": "implicit_ff",
                                "strands_numbers": 105,
                                "strand_radii": 0.1e-3 / 2,
                                "conductor_radii": 1.5e-3 / 2,
                                "ff": ""}
    litz_dict["1.4x200x0.071"] = {"implicit": "implicit_ff",
                                  "strands_numbers": 200,
                                  "strand_radii": 0.071e-3 / 2,
                                  "conductor_radii": 1.4e-3 / 2,
                                  "ff": ""}
    litz_dict["2.0x405x0.071"] = {"implicit": "implicit_ff",
                                  "strands_numbers": 405,
                                  "strand_radii": 0.071e-3 / 2,
                                  "conductor_radii": 2.0e-3 / 2,
                                  "ff": ""}
    litz_dict["2.0x800x0.05"] = {"implicit": "implicit_ff",
                                 "strands_numbers": 800,
                                 "strand_radii": 0.05e-3 / 2,
                                 "conductor_radii": 2e-3 / 2,
                                 "ff": ""}

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

def cost_material_database() -> Dict:
    """
    Returns costs for core and winding.
    This is splitted in material and fabrication costs.
    Both, material and fabrication costs have a euro_per_kilogram and a euro_per_unit (fixcosts) price.

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
                                                                    "sigma_denumerator": 0.45}

    cost_database["gross_margin"] = 0.25


    return cost_database

def pm_core_inner_diameter_calculator(inner_core_diameter: float, hole_diameter: float) -> float:
    """
    Calculates the effective inner core diameter without the hole
    Often used in PM-cores

    :param inner_core_diameter: inner core diameter
    :type inner_core_diameter: float
    :param hole_diameter: hole diameter
    :type hole_diameter: float
    :return: effective inner core diameter without hole
    :rtype: float
    """
    area_inner_core_inner_diameterithout_hole = (inner_core_diameter / 2) ** 2 * np.pi
    area_hole = (hole_diameter / 2) ** 2 * np.pi
    area_total = area_inner_core_inner_diameterithout_hole - area_hole

    return np.around(2 * np.sqrt(area_total / np.pi), decimals=4)


def install_pyfemm_if_missing() -> None:
    """
    Windows users only.
    Installs femm-software pip package in case of running on windows machine

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
    [min, max] = [None, None]
    output = input_points
    dim = None
    # Find equal dimension
    for i in range(0, 3):
        if a[i] == b[i] and a[i] != 0:
            dim = i
    if dim == None:
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
        raise Exception("Not implemented Error: Only 2D is implemeted")
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

    [min, max] = [None, None]
    buffer = input_points
    dim = None
    # Find equal dimension
    for i in range(0, 3):
        if a[i] == b[i] and a[i] != 0:
            dim = i
    if dim == None:
        raise Exception("Given points do not have a common dimension")
    n = 0
    while n < buffer.shape[0]:
        if a[dim] != buffer[n, dim]:
            buffer = np.delete(buffer, n, 0)
        else:
            n += 1
    if buffer.shape[0] == 0:
        print("No air gaps between interval borders")
    if buffer.shape[0] % 2 == 1:
        raise Exception("Odd number of input points")
    if dim == 2:
        raise Exception("Not implemented Error: Only 2D is implemeted")
    dim2 = (dim+1) % 2
    if buffer.shape[0] >= 2:
        argmax = np.argmax(buffer[:, 1])
        max = buffer[argmax]
        argmin = np.argmin(buffer[:, 1])
        min = buffer[argmin]
    return [min, max]
    

def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


def NbrStrands(n_layers: int) -> int:
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


def NbrLayers(n_strands: int) -> int:
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


def fft(period_vector_t_i: npt.ArrayLike, sample_factor: float = 1000, plot: str = 'no', mode: str = 'rad',
        f0: Union[float, None] = None, title: str = 'ffT', filter_type: str = 'factor',
        filter_value_factor: float = 0.01, filter_value_harmonic: int = 100,
        figure_size: Tuple=None, figure_directory: str=None) -> npt.NDArray[list]:
    """
    A fft for a input signal. Input signal is in vector format and should include one period.

    Output vector includes only frequencies with amplitudes > 1% of input signal

    :Minimal Example:

    >>> import femmt as fmt
    >>> import numpy as np
    >>> example_waveform = np.array([[0, 1.34, 3.14, 4.48, 6.28],[-175.69, 103.47, 175.69, -103.47,-175.69]])
    >>> out = fmt.fft(example_waveform, plot='yes', mode='rad', f0=25000, title='ffT input current')

    :param period_vector_t_i: numpy-array [[time-vector[,[current-vector]]. One period only
    :type period_vector_t_i: np.array
    :param sample_factor: f_sampling/f_period, defaults to 1000
    :type sample_factor: float
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
    :param figure_size: None for auto fit; fig_size for matplotlib (width, length)
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
        raise ValueError("Mode not availabe. Choose: 'rad', 'deg', 'time'")

    t = period_vector_t_i[0]
    i = period_vector_t_i[1]

    # fft-function works per default in time domain
    # time domain
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
    for n, dict in enumerate(data):
        for key in keys:
            if not key in dict:
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
mu0 = 4e-7*np.pi

def r_basic_round_inf(air_gap_radius, air_gap_basic_hight, core_hight):
    """
    Do not use this function directly!
    Use it indirectly by using
     - r_air_gap_round_round
     - r_air_gap_round_inf
    instead!

    This function calculates the r_basic for a round to infinite structure according to the following paper:
    [according to "A Novel Approach for 3D Air Gap Reluctance Calculations" - J. Mühlethaler, J.W. Kolar, A. Ecklebe]

    :param air_gap_radius: air gap radius
    :param air_gap_basic_hight: air gap hight for the BASIC-AIR-GAP (e.g. if you use a round-round structure, this is half of the total air gap).
    :param core_hight: core hight
    :return: basic reluctance for round - infinite structure
    """
    conductance_basic = mu0 * (air_gap_radius * 2 / 2 / air_gap_basic_hight + 2 / np.pi * (1 + np.log(np.pi * core_hight / 4 / air_gap_basic_hight)))

    return 1 / conductance_basic

def sigma_round(r_equivalent, air_gap_radius, air_gap_total_hight):
    """
    Do not use this function directly!
    Use it indirectly by using
     - r_air_gap_round_round
     - r_air_gap_round_inf
    instead!

    :param r_equivalent: this is a series/parallel connection of r_basic, depending on the air gap structure
    :param air_gap_radius: air gap radius
    :param air_gap_total_hight: air gap total hight (for the total air gap, also for round-round structures)
    :return: fringing factor 'sigma'
    """
    return r_equivalent * mu0 * air_gap_radius / air_gap_total_hight

def r_air_gap_round_round(air_gap_total_hight, core_inner_diameter , core_hight_upper, core_hight_lower):
    """
    Returns the reluctance of a round-round air gap structure and includes finging effects.

    :param air_gap_total_hight: total air gap hight of the air gap
    :param core_inner_diameter: core inner diameter
    :param core_hight_upper: core hight upper (needed for better calculating fringing effects)
    :param core_hight_lower: core hight lower (needed for better calculating fringing effects)
    :return: air gap reluctance for round-round structure including fringing effects
    """
    core_inner_diameter = np.array(core_inner_diameter)
    air_gap_radius = core_inner_diameter / 2

    air_gap_basic_hight = air_gap_total_hight / 2
    r_basic_upper = r_basic_round_inf(air_gap_radius, air_gap_basic_hight, core_hight_upper)
    r_basic_lower = r_basic_round_inf(air_gap_radius, air_gap_basic_hight, core_hight_lower)
    print(f"{r_basic_upper = }")

    r_equivalent_round_round = r_basic_upper + r_basic_lower

    sigma = sigma_round(r_equivalent_round_round, air_gap_radius, air_gap_total_hight)
    if np.any(sigma) > 1:
        raise Exception("Failure in calculting reluctance. Sigma was calculated to >1. Check input parameters!")

    r_air_gap_ideal = air_gap_total_hight / mu0 / np.pi / (air_gap_radius ** 2)
    r_air_gap = sigma ** 2 * r_air_gap_ideal

    return r_air_gap

def r_air_gap_round_inf(air_gap_total_hight, core_inner_diameter, core_hight):
    """
    Returns the reluctance of a round-infinite air gap structure and includes fringing effects

    :param air_gap_total_hight: total air gap hight of the air gap
    :param core_inner_diameter: core inner diameter
    :param core_hight: core hight (needed for better calculating fringing effects)
    :return: air gap reluctance for round-inf structure including fringing effects
    """

    air_gap_total_hight = np.array(air_gap_total_hight)
    core_inner_diameter = np.array(core_inner_diameter)
    core_hight = np.array(core_hight)

    air_gap_radius = core_inner_diameter / 2
    r_basic = r_basic_round_inf(air_gap_radius, air_gap_total_hight, core_hight)

    r_equivalent_round_inf = r_basic
    sigma = sigma_round(r_equivalent_round_inf, air_gap_radius, air_gap_total_hight)

    r_air_gap_ideal = air_gap_total_hight / mu0 / np.pi / (air_gap_radius ** 2)
    r_air_gap = sigma ** 2 * r_air_gap_ideal

    return r_air_gap




def r_basic_tablet_cyl(tablet_hight, air_gap_basic_hight, tablet_radius):
    """
    Do not use this function directly!
    Use it indirectly by using
     - r_air_gap_tablet_cyl
    instead!

    This function calculates the r_basic for a round to infinite structure according to the following paper:
    [according to "A Novel Approach for 3D Air Gap Reluctance Calculations" - J. Mühlethaler, J.W. Kolar, A. Ecklebe]

    Note: this is the same function as r_basic_round_inf, but with clear variable names for tablet-cylinder structure

    :param tablet_hight: tablet hight = air gap width for tablet-cylinder structure
    :param air_gap_basic_hight: air gap hight for the BASIC-AIR-GAP (e.g. if you use a round-round structure, this is half of the total air gap).
    :param tablet_radius: tablet radius
    :return: basic reluctance for tablet - cylinder structure
    """
    conductance_basic = mu0 * (tablet_hight / 2 / air_gap_basic_hight + 2 / np.pi * (1 + np.log(np.pi * tablet_radius / 4 / air_gap_basic_hight)))

    return 1 / conductance_basic

def sigma_tablet_cyl(r_equivalent, tablet_hight, air_gap_total_hight):
    """
    Do not use this function directly!
    Use it indirectly by using
     - r_air_gap_tablet_cyl
    instead!

    Note: this is the same function as sigma_round, but with clear variable names for tablet-cylinder structure

    :param r_equivalent: this is a series/parallel connection of r_basic, depending on the air gap structure
    :param tablet_hight: tablet hight
    :param air_gap_total_hight: air gap total hight (for the total air gap)
    :return: fringing factor 'sigma' for tablet - cylinder structure
    """
    return r_equivalent * mu0 * tablet_hight / air_gap_total_hight


def r_air_gap_tablet_cyl(tablet_hight, air_gap_total_hight, r_outer):
    """
    Returns the reluctance of a cylinder-tablet air gap structure and includes fringing effects
    This function calculates the air gap reluctance for a 2D-axisymmetric core.

    :param tablet_hight: tablet hight
    :param air_gap_total_hight: total air gap hight
    :param r_outer: radius of outer core window
    :return: air gap reluctance for tablet - cylinder structure including air gap fringing
    """

    # translate practical core dimensions to non-practial air-gap dimensions
    tablet_radius = r_outer - air_gap_total_hight

    air_gap_basic_hight = air_gap_total_hight
    r_basic = r_basic_tablet_cyl(tablet_hight, air_gap_basic_hight, tablet_radius)

    r_equivalent = r_basic / 2
    sigma = sigma_tablet_cyl(r_equivalent, tablet_hight, air_gap_total_hight)
    if sigma > 1:
        raise Exception("Failure in calculting reluctance. Sigma was calculated to >1. Check input parameters!")

    r_air_gap_ideal = np.log(r_outer / (r_outer - air_gap_total_hight)) / 2 / mu0 / np.pi / tablet_hight

    r_air_gap = sigma * r_air_gap_ideal

    return r_air_gap


def r_air_gap_tablet_cyl_no_2d_axi(tablet_hight, tablet_diameter, r_outer, real_core_width_no_2d_axi):
    """
    Returns the reluctance of a cylinder-tablet air gap structure and includes fringing effects
    Note:
    This function differes from r_air_gap_tablet_cyl (ideal 2D axisymmetric core). Here, the air gap reluctance for
    a non-2D-axisymmetric core is taken into account, as a real PQ core is open at the side. So, there is no air gap
    taken into account for the side-sections. The new real_core_width_no_2d_axi parameter describes the width of the
    core when you are in a xy-coordinate system.

    :param tablet_hight: tablet hight
    :param tablet_diameter: tablet diameter
    :param r_outer: radius of outer core window
    :param real_core_width_no_2d_axi: core width for a real core (e.g. PQ-core) in xy-coordinate system.
    :return: air gap reluctance for tablet - cylinder structure including air gap fringing
    """

    if tablet_diameter / 2 >= r_outer:
        raise Exception("tablet radius is greater than r_outer")

    # translate practical core dimensions to non-practial air-gap dimensions
    air_gap_total_hight = r_outer - tablet_diameter / 2

    air_gap_basic_hight = air_gap_total_hight
    r_basic = r_basic_tablet_cyl(tablet_hight, air_gap_basic_hight, tablet_diameter / 2)

    r_equivalent = r_basic / 2
    sigma = sigma_tablet_cyl(r_equivalent, tablet_hight, air_gap_total_hight)
    if sigma > 1:
        raise Exception("Failure in calculting reluctance. Sigma was calculated to >1. Check input parameters!")

    r_air_gap_ideal = np.log(r_outer / (r_outer - air_gap_total_hight)) / mu0 / (2 * np.pi - 4 * np.arccos(real_core_width_no_2d_axi / 2 / r_outer)) / tablet_hight

    r_air_gap = sigma * r_air_gap_ideal

    return r_air_gap

def r_core_tablet(tablet_hight, tablet_radius, mu_r, core_inner_diameter):
    """
    Calculates the magentic resistance of the core tablet

    :param tablet_hight: tablet hight
    :param tablet_radius: tablet radius
    :param mu_r: relative permeability (mu_r) of the core material from datasheet
    :param core_inner_diameter: core inner diameter. For idealized core material, this value can be 0.001.
    """

    return np.log(tablet_radius / (core_inner_diameter / 2)) / ( 2 * np.pi * mu0 * mu_r * tablet_hight)


def r_core_top_bot_radiant(core_inner_diameter, window_w, mu_r, core_top_bot_hight):
    """
    Calculates the top or bottom core material part

    :param core_inner_diameter: core inner diameter
    :param window_w: width of winding window
    :param mu_r: relative permeability (mu_r) of the core material from datasheet
    :param core_top_bot_hight: hight of the core material top / bottom of the winding window
    """

    return np.log( (core_inner_diameter + 2 * window_w) / core_inner_diameter) / ( 2 * np.pi * mu0 * mu_r * core_top_bot_hight)

def r_core_round(core_inner_diameter, core_round_hight, mu_r):
    """
    Calculates the core reluctance for a round structure

    :param core_round_hight: hight of the round core part section
    :param core_inner_diameter: core inner diameter
    :param mu_r: relative permeability (mu_r) of the core material from datasheet
    """
    return core_round_hight / ( mu0 * mu_r * (core_inner_diameter / 2) ** 2 * np.pi)


def r_top_bot_stray(core_inner_diameter, air_gap_middle_leg_list, window_w, window_h, stray_path_air_gap_length, mu_r, start_index, position_air_gap_percent_list):

    # core geometry calculations

    middle_leg_top_total_hight = (100 - position_air_gap_percent_list[start_index + 1]) * 0.01 * window_h + air_gap_middle_leg_list[start_index + 1] / 2
    middle_leg_bot_total_hight = position_air_gap_percent_list[start_index] * 0.01 * window_h + air_gap_middle_leg_list[start_index] / 2
    tablet_hight = window_h - middle_leg_top_total_hight - middle_leg_bot_total_hight

    start_index_air_gap_top = start_index + 1
    stop_index_air_gap_top = len(air_gap_middle_leg_list) - 1
    #start_index_air_gap_bot = 0
    #stop_index_air_gap_bot = start_index

    air_gap_middle_leg_bot_list = []
    position_air_gap_percent_bot_list = []
    for count, air_gap in enumerate(air_gap_middle_leg_list):
        if count <= start_index:
            air_gap_middle_leg_bot_list.append(air_gap)
            position_air_gap_percent_bot_list.append(position_air_gap_percent_list[count])


    # air gap geometry calculations
    air_gap_total_length_top = 0.0
    air_gap_middle_leg_top_list = []
    air_gap_middle_leg_top_percent_list = []
    for count in range(start_index + 1, len(air_gap_middle_leg_list)):
        #print(f"{air_gap_middle_leg_list[count] = }")
        air_gap_total_length_top += air_gap_middle_leg_list[count]

    air_gap_total_length_bot = 0.0
    for count in range(0, start_index + 1):
        air_gap_total_length_bot += air_gap_middle_leg_list[count]
    #print(f"{air_gap_total_length_top = }")
    #print(f"{air_gap_total_length_bot = }")

    # calculate r_top
    core_round_hight_top = middle_leg_top_total_hight - air_gap_total_length_top

    r_core_round_top = r_core_round(core_inner_diameter, core_round_hight_top, mu_r)
    #ToDo: radiant calculation core hight is very simplified using core_inner_diameter/2
    r_core_top_radiant = r_core_top_bot_radiant(core_inner_diameter, window_w, mu_r, core_inner_diameter/2)
    r_core_outer_top = r_core_round_top
    r_core_top = r_core_round_top + r_core_top_radiant + r_core_outer_top

    # calculate r_air_gap_top
    # sweep trough top air gap list
    r_air_gap_top = 0
    for count in range(start_index_air_gap_top, stop_index_air_gap_top+1):
        # this routine sums up all air gaps inside the middle leg above the stray path
        # there are three different types of caluclations
        #  - first air gap is definitely a round_inf structure
        #  - in cases of last air gap is >95% of window_h, last air gap is treated as round_inf structure
        #  - other air gaps are treated as round_round structure
        if count == start_index + 1:
            # this is for the first air gap, what is definetly a round_inf air gap
            if start_index_air_gap_top == stop_index_air_gap_top:
                core_hight = (100 - position_air_gap_percent_list[count]) * 0.01 * window_h
            else:
                core_hight = (position_air_gap_percent_list[count + 1] - position_air_gap_percent_list[count]) * 0.01 *window_h / 2
            r_air_gap_top += r_air_gap_round_inf(air_gap_middle_leg_list[count], core_inner_diameter, core_hight)
            #print('### Case 1: first air gap for top')
            #print(f"{core_hight = }")
            #print(f"{r_air_gap_top = }")

        elif position_air_gap_percent_list[count] > 95 and count != start_index_air_gap_top:
            # this is for the last air gap in case of very close to the top core (95%), so there will be the assumption for a round-inf air gap
            core_hight = (position_air_gap_percent_list[stop_index_air_gap_top] - position_air_gap_percent_list[stop_index_air_gap_top - 1]) * 0.01 * window_h / 2
            r_air_gap_top += r_air_gap_round_inf(air_gap_middle_leg_list[count], core_inner_diameter, core_hight)
            #print('### Case 2: last air gap for top')
            #print(f"{core_hight = }")
            #print(f"{r_air_gap_top = }")
        else:
            # air gap in the middle between tablet and top air gap
            if count + 1 < stop_index_air_gap_top:
                # this is for multiple air gaps in the top-section. Calculation of core hight is only to the next air gap
                core_hight_upper = (position_air_gap_percent_list[count + 1] - position_air_gap_percent_list[count]) * 0.01 * window_h
            else:
                # this is for the last (upper) air gap in the top-section. Calculation of core_hight is until the end of the window
                core_hight_upper = (100 - position_air_gap_percent_list[count]) * 0.01 * window_h
            core_hight_lower = (position_air_gap_percent_list[count] - position_air_gap_percent_list[count - 1]) * 0.01 * window_h
            r_air_gap_top += r_air_gap_round_round(air_gap_middle_leg_list[count], core_inner_diameter, core_hight_upper, core_hight_lower)
            #print('### Case 3: middle air gap for top')
            #print(f"{core_hight_upper = }")
            #print(f"{core_hight_lower = }")
            #print(f"{r_air_gap_top = }")


    r_top = r_core_top + r_air_gap_top

    # calculate r_bot
    core_round_hight_bot = middle_leg_bot_total_hight - air_gap_total_length_bot

    r_core_round_bot = r_core_round(core_inner_diameter, core_round_hight_bot, mu_r)
    #ToDo: radiant calculation core hight is very simplified using core_inner_diameter/2
    r_core_bot_radiant = r_core_top_bot_radiant(core_inner_diameter, window_w, mu_r, core_inner_diameter/2)
    r_core_outer_bot = r_core_round_bot
    r_core_bot = r_core_round_bot + r_core_bot_radiant + r_core_outer_bot

    #print('##########################################')
    #print("#### bottom air gap  ####")
    #print('##########################################')

    # calculate r_air_gap_bot
    # sweep trough bot air gap list
    r_air_gap_bot = 0
    for count, air_gap in enumerate(air_gap_middle_leg_bot_list):
        #print(f"{air_gap = }")
        # this routine sums up all air gaps inside the middle leg below the stray path
        # there are three different types of caluclations
        #  - last air gap is definitely a round_inf structure
        #  - in cases of first air gap is <5% of window_h, first air gap is treated as round_inf structure
        #  - other air gaps are treated as round_round structure
        if count == position_air_gap_percent_bot_list.index(max(position_air_gap_percent_bot_list)):
            # this is for the last air gap, what is definetly a round_inf air gap
            # the check checks not for the last position in the list, but for the real highest percent value
            if len(position_air_gap_percent_bot_list) == 1:
                core_hight = (position_air_gap_percent_bot_list[count]) * 0.01 * window_h
            else:
                core_hight = (position_air_gap_percent_bot_list[count] - position_air_gap_percent_bot_list[count - 1]) * 0.01 * window_h / 2
            r_air_gap_bot += r_air_gap_round_inf(air_gap_middle_leg_list[count], core_inner_diameter, core_hight)
            #print('### Case 1: first air gap for bot')
            #print(f"{core_hight = }")
            #print(f"{r_air_gap_bot = }")

        elif position_air_gap_percent_list[count] < 5  and count == 0:
            # this is for the first air gap in case of very close to the bot core (<5%), so there will be the assumption for a round-inf air gap
            core_hight = (position_air_gap_percent_list[1] - position_air_gap_percent_list[0]) * 0.01 * window_h / 2
            r_air_gap_bot += r_air_gap_round_inf(air_gap_middle_leg_list[count], core_inner_diameter, core_hight)
            #print('### Case 2: last air gap for bot')
            #print(f"{core_hight = }")
            #print(f"{r_air_gap_bot = }")
        else:
            # this is for multiple air gaps in the bot-section. Calculation of core hight is only to the next air gap
            if count + 1 < len(position_air_gap_percent_bot_list):
                core_hight_lower = (position_air_gap_percent_bot_list[count + 1] - position_air_gap_percent_bot_list[count]) * 0.01 * window_h / 2
            else:
                # this is for the first (lower) air gap in the bot-section. Calculation of core_hight is until the end of the window
                core_hight_lower = (position_air_gap_percent_bot_list[count]) * 0.01 * window_h / 2

            core_hight_upper  = (position_air_gap_percent_bot_list[count + 1] - position_air_gap_percent_bot_list[count]) * 0.01 * window_h / 2
            r_air_gap_bot += r_air_gap_round_round(air_gap_middle_leg_list[count], core_inner_diameter, core_hight_upper, core_hight_lower)
            #print('### Case 3: middle air gap for bot')
            #print(f"{core_hight_upper = }")
            #print(f"{core_hight_lower = }")
            #print(f"{r_air_gap_bot = }")

    r_bot = r_core_bot + r_air_gap_bot


    # calculate r_stray
    r_stray_air_gap = r_air_gap_tablet_cyl(tablet_hight, stray_path_air_gap_length, core_inner_diameter / 2 + window_w)
    r_stray_core = r_core_tablet(tablet_hight, core_inner_diameter / 2 + window_w - stray_path_air_gap_length, mu_r, core_inner_diameter)
    r_stray = r_stray_air_gap + r_stray_core

    return r_top, r_bot, r_stray

def calculate_reluctances(winding_matrix, inductance_matrix):
    """
    Calculates the Reluctance Matrix.
    Everything must be numpy!

    L. Keuck, "Entwurf eines einstufigen Ladewandlers auf Basis eines LLC-Resonanzwandlers", dissertation 2023

    :param winding_matrix: winding matrix
    :param inductance_matrix: inductance matrix
    :return: reluctance[-matrix]

    inductance matrix e.g.
    L = [ [L_11, M], [M, L_22] ]

    winding matrix e.g.
    N = [ [N_1a, N_2b], [N_1b, N_2b] ]

    """

    # Reluctance Matrix
    if np.ndim(winding_matrix) == 0:
        L_invert = 1 / inductance_matrix
    else:
        L_invert = np.linalg.inv(inductance_matrix)

    return np.matmul(np.matmul(winding_matrix, L_invert), np.transpose(winding_matrix))


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

    print(inductance)
    print(loss_core_eddy_current)
    print(loss_core_hysteresis)
    print(loss_winding_1)

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

def set_silent_status(s: bool):
    global silent
    silent = s

def femmt_print(text, end='\n'):
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


def cost_function_winding(wire_weight_list: List[float], wire_type_list: List[str], single_strand_cross_section_list: List[float] = []):
    """
    Calculates single winding material and fabrication costs depending on winding-type and weight

    :param wire_weight_list: winding weight in kg in list-form
    :type wire_weight_list: List[float]
    :param wire_type_list: winding type. Must fit to enum-names in ConductorType-Enum
    :type wire_type_list: List[str]
    :param single_strand_cross_section_list: single strand cross section in list-form
    :type single_strand_cross_section_list: List[float]
    :return: winding cost of single winding
    :rtype: float
    """
    cost_database = cost_material_database()
    winding_cost_list = []

    for winding_count, winding_weight in enumerate(wire_weight_list):
        # material cost (per kilogram and per unit)
        sigma_material_winding_euro_per_kilogram = cost_database["winding_material_euro_per_kilogram"][wire_type_list[winding_count]]
        if sigma_material_winding_euro_per_kilogram == -1:
            # case for special litz wire calculation. Additional data is loaded from cost_database.
            sigma_material_winding_euro_per_kilogram = cost_database["winding_material_euro_per_kilogram_for_litz"]["sigma_numerator"] / (single_strand_cross_section_list[winding_count] * 1e6 + cost_database["winding_material_euro_per_kilogram_for_litz"]["sigma_denumerator"])

        winding_material_euro_per_unit = cost_database["winding_material_euro_per_unit"][wire_type_list[winding_count]]

        winding_material_cost = sigma_material_winding_euro_per_kilogram * winding_weight + winding_material_euro_per_unit

        # fabrication cost (per kilogram and per unit)
        sigma_fabrication_euro_per_kilogram = cost_database["winding_fabrication_euro_per_kilogram"][wire_type_list[winding_count]]
        fabrication_material_euro_per_unit = cost_database["winding_fabrication_euro_per_unit"][wire_type_list[winding_count]]

        winding_fabrication_cost = sigma_fabrication_euro_per_kilogram * winding_weight + fabrication_material_euro_per_unit

        winding_cost_list.append(winding_material_cost + winding_fabrication_cost)

    return winding_cost_list



def cost_function_total(core_weight: float, core_type: str, wire_weight_list: List[float], wire_type_list: List[str],
                        single_strand_cross_section_list: List[float] = []) -> float:
    """
    Calculates the total costs for a inductive element.
    This includes material costs for core and winding, fabrication costs for core and winding and manufacturer margin

    :param core_weight: core weight in kg
    :type core_weight: float
    :param core_type: core type. Can be "ferrite", "amorphous", "nanocristalline", "high_si_steel", "goes"
    :type core_type: str
    :param wire_weight_list: winding weight in kg
    :type wire_weight_list: float
    :param wire_type_list: winding type in list-form. Must fit to enum-names in ConductorType-Enum
    :type wire_type_list: List[str]
    :param single_strand_cross_section_list: single strand cross section in list-form
    :type single_strand_cross_section_list: List[float]
    :return: total costs for inductive element
    :rtype: float
    """
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

        if value_min <= data_to_compare and data_to_compare <= value_max:
            print(f"{value_min} <= {data_to_compare} <= {value_max} for file named {file}")


if __name__ == '__main__':
    pass
