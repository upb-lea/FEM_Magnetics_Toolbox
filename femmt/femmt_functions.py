# Usual Python libraries
import numpy.typing as npt
import numpy as np
from matplotlib import pyplot as plt
import json
import random
import string
import pkg_resources
import subprocess
import sys
import os
import pandas as pd
import time
import gmsh
import warnings
from typing import Union, List, Tuple, Dict


def core_database(core_type: str) -> Dict:
    """
    Returns core geometry for defined core structures

    :param core_type: core type, e.g. "PQ4040"
    :type core_type: str
    :return: Dict including core_h, core_w, window_h, window_w
    :rtype: Dict
    """
    core_dict = {}
    PQ4040_dict = {
        "core_h": 39.8e-3,
        "core_w": 14.5e-3,
        "window_h": 29.5e-3 ,
        "window_w": 11.05e-3,
    }
    core_dict["PQ4040"] = PQ4040_dict

    return core_dict[core_type]


def litz_database(litz_type: str) -> Dict:
    """
    Returns litz parameters for defined litz wires.

    :param litz_type: str
    :type litz_type: str
    :return: Dict including litz parameters like strand_numbers, strand_radii and conductor_radii
    :rtype: Dict
    """
    litz_dict = {}

    litz_dict["105x0.1"] = {"strands_numbers": 105,
                            "strand_radii": 0.1e-3 / 2,
                            "conductor_radii": 1.5e-3 / 2}
    litz_dict["200x0.071"] = {"strands_numbers": 200,
                              "strand_radii": 0.071e-3 / 2,
                              "conductor_radii": 1.4e-3 / 2}
    litz_dict["405x0.071"] = {"strands_numbers": 405,
                              "strand_radii": 0.071e-3 / 2,
                              "conductor_radii": 2.0e-3 / 2}
    litz_dict["800x0.05"] = {"strands_numbers": 800,
                             "strand_radii": 0.05e-3 / 2,
                             "conductor_radii": 2e-3 / 2}

    return litz_dict[litz_type]


# Static Functions
#  Used somewhere in the Code of femmt.py
def install_femm_if_missing() -> None:
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
        print(f"f_out = {np.around(f_out, 0)}")
        print(f"x_out = {np.around(x_out, 1)}")
        print(f"phi_rad_out = {np.around(phi_rad_out, 1)}")

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


def r_basis(basis_air_gap_length: float, basis_air_gap_width: float, basis_air_gap_height_core_material: float) -> None:
    """
    1-dim reluctance per-unit-of-length
    [according to "A Novel Approach for 3D Air Gap Reluctance Calculations" - J. Mühlethaler, J.W. Kolar, A. Ecklebe]

    :param basis_air_gap_width: width of basis air gap
    :type basis_air_gap_width: float
    :param basis_air_gap_length: length of basis air gap
    :type basis_air_gap_length: float
    :param basis_air_gap_height_core_material: height of core material direct on the air gap
    :type basis_air_gap_height_core_material: float

    :return: basis air gap
    :rtype: float

    """
    if basis_air_gap_length <= 0:
        basis_air_gap_length = 0.0000001
    return 1 / (mu0 * (basis_air_gap_width / 2 / basis_air_gap_length + 2 / np.pi * (1 + np.log(np.pi * basis_air_gap_height_core_material / 4 / basis_air_gap_length))))


def sigma(l, w, r_equivalent):
    """
    1-dim fringing factor
    [according to "A Novel Approach for 3D Air Gap Reluctance Calculations" - J. Mühlethaler, J.W. Kolar, A. Ecklebe]

    :param w:
    :type w: float
    :param l:
    :type l: float
    :param r_equivalent:
    :type r_equivalent: float

    :return:
    :rtype: float

    """
    return r_equivalent / (l / mu0 / w)


def r_round_inf(air_gap_length: float, sigma: float, air_gap_radius: float) -> float:
    """
    3-dim reluctance for 2-dim axial-symmetric air gaps
    Round to infinity structure

    :param sigma: fringing factor
    :type sigma: float
    :param air_gap_radius: air gap radius
    :type air_gap_radius: float
    :param air_gap_length: air gap length
    :type air_gap_length: float

    :return: air gap reluctance
    :rtype: float

    """
    return sigma ** 2 * air_gap_length / mu0 / air_gap_radius ** 2 / np.pi


def r_round_round(air_gap_length: float, sigma: float, air_gap_radius: float) -> float:
    """
    3-dim reluctance for 2-dim axial-symmetric air gaps
    Round to round structure

    :param sigma: fringing factor
    :type sigma: float
    :param air_gap_radius: air gap radius
    :type air_gap_radius: float
    :param air_gap_length: air gap length
    :type air_gap_length: float

    :return: air gap reluctance
    :rtype: float

    """
    return sigma ** 2 * air_gap_length / mu0 / air_gap_radius ** 2 / np.pi


def r_cyl_cyl(air_gap_length: float, sigma: float, air_gap_width: float, radius_outer) -> float:
    """
    Cylinder-cylinder air gap reluctance

    :param air_gap_length: air gap length
    :type air_gap_length: float
    :param sigma: fringing factor
    :type sigma: float
    :param air_gap_width: air gap width
    :type air_gap_width: float
    :param radius_outer: outer cylinder radius
    :type radius_outer: float

    :return: air gap reluctance
    :rtype: float

    """
    return sigma * np.log(radius_outer / (radius_outer - air_gap_length)) / 2 / mu0 / np.pi / air_gap_width


def r_cyl_cyl_real(air_gap_length: float, sigma: float, air_gap_width: float, radius_outer: float,
                   real_core_heigth: float) -> float:
    """

    :param air_gap_length: air gap length
    :type air_gap_length: float
    :param sigma: fringing factor
    :type sigma: float
    :param air_gap_width: air gap width
    :type air_gap_width: float
    :param radius_outer: outer cylinder radius
    :type radius_outer: float

    :return: air gap reluctance
    :rtype: float

    """
    return sigma * np.log(radius_outer / (radius_outer - air_gap_length)) / mu0 / (2 * np.pi - 4 * np.arccos(real_core_heigth / 2 / radius_outer)) / air_gap_width


def r_cheap_cyl_cyl(radius_outer: float, air_gap_length: float, air_gap_width: float) -> float:
    """
    Simplified method to calculate cylinder-cylinder air gap reluctance

    :param radius_outer: outer cylinder radius
    :type radius_outer: float
    :param air_gap_length: air gap length
    :type air_gap_length: float
    :param air_gap_width: air gap width
    :type air_gap_width: float

    :return: air gap reluctance
    :rtype: float

    """
    r_i = radius_outer - air_gap_length
    return (radius_outer - r_i) / mu0 / air_gap_width / np.pi / (radius_outer + r_i)


def calculate_reluctances(N, L):
    """
    Calculates the Reluctance Matrix.
    Everything must be numpy!

    :return: reluctance[-matrix]

    """

    # Reluctance Matrix
    if np.ndim(N) == 0:
        L_invert = 1 / L
    else:
        L_invert = np.linalg.inv(L)

    return np.matmul(np.matmul(N, L_invert), np.transpose(N))

def create_physical_group(dim, entities, name):
    tag = gmsh.model.addPhysicalGroup(dim, entities)
    gmsh.model.setPhysicalName(dim, tag, name)

    return tag


if __name__ == '__main__':
    pass
