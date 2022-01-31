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
import warnings
from typing import Union, List, Tuple


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


def NbrStrands(n_layers):
    """
    Returns the number of strands in a hexagonal litz winding with a
    specified number of layers (n_layers). CAUTION: Zero number of
    layers corresponds to a single strand.

    :param n_layers:

    :return:

    """
    return 3 * (n_layers + 1) ** 2 - 3 * (n_layers + 1) + 1


def NbrLayers(n_strands):
    """
    Returns the number of layers in a hexagonal litz winding with a
    specified number of strands (n_strands).

    .. note:: Zero number of layers corresponds to a single strand.

    :param n_strands:

    :return:

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

    >>> example_waveform = np.array([[0, 1.34, 3.14, 4.48, 6.28],[-175.69, 103.47, 175.69, -103.47,-175.69]])
    >>> out = fft(example_waveform, plot=True, mode='rad', f0=25000, title='ffT input current')

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


def store_as_npy_in_directory(dir_path, file_name, data):
    if not os.path.isdir(dir_path):
        os.mkdir(dir_path)
    np.save(dir_path + "/" + file_name, data)


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
    :param kwargs: keys and values in dicts
    :return:
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


def find_common_frequencies(amplitudes1, phases1, frequencies1, amplitudes2, phases2, frequencies2):
    """

    :param amplitudes1:
    :param phases1:
    :param frequencies1:
    :param amplitudes2:
    :param phases2:
    :param frequencies2:
    :return:
    """
    common_amplitudes1 = []
    common_phases1 = []
    common_amplitudes2 = []
    common_phases2 = []

    # Find all common frequencies
    # print(f"{frequencies1=}")
    # print(f"{frequencies2=}")
    # print(f"{phases1=}")
    # print(f"{phases2=}")
    # print(f"{amplitudes1=}")
    # print(f"{amplitudes2=}")

    common_frequencies = list(set(frequencies1).intersection(frequencies2))
    # print(f"{common_frequencies=}")

    # Delete the corresponding phases and amplitudes
    if isinstance(amplitudes1, list):
        for frequency in common_frequencies:
            common_amplitudes1.append(amplitudes1[frequencies1.index(frequency)])
            common_phases1.append(phases1[frequencies1.index(frequency)])
            common_amplitudes2.append(amplitudes2[frequencies2.index(frequency)])
            common_phases2.append(phases2[frequencies2.index(frequency)])
    if isinstance(amplitudes1, np.ndarray):
        for frequency in common_frequencies:
            common_amplitudes1.append(amplitudes1[np.where(frequencies1==frequency)][0])
            common_phases1.append(phases1[np.where(frequencies1==frequency)][0])
            common_amplitudes2.append(amplitudes2[np.where(frequencies2==frequency)][0])
            common_phases2.append(phases2[np.where(frequencies2==frequency)][0])
    else:
        warnings.warn("Either a list or a np.ndarray must be provided!")

    current_pairs = list(map(list, zip(common_amplitudes1, common_amplitudes2)))
    phase_pairs = list(map(list, zip(common_phases1, common_phases2)))

    return current_pairs, phase_pairs, common_frequencies


def sort_out_small_harmonics(phase_pairs, amplitude_pairs, frequencies, limes):
    # Calculate geometric lengths
    amp_tot = np.sqrt(np.sum(np.array(amplitude_pairs)**2, axis=0))
    # amp_tot = np.max(amplitude_pairs, axis=0)

    invalid_index = []
    for n, amplitude_pair in enumerate(amplitude_pairs):
        if all(amplitude/amp_tot[i] < limes for i, amplitude in enumerate(amplitude_pair)):
            invalid_index.append(n)

    phase_pairs = np.delete(phase_pairs, invalid_index, axis=0)
    amplitude_pairs = np.delete(amplitude_pairs, invalid_index, axis=0)
    frequencies = np.delete(frequencies, invalid_index)

    return phase_pairs, amplitude_pairs, frequencies



# Reluctance Model [with calculation]
mu0 = 4e-7*np.pi


def r_basis(l, w, h):
    """
    1-dim reluctance per-unit-of-length
    [according to "A Novel Approach for 3D Air Gap Reluctance Calculations" - J. Mühlethaler, J.W. Kolar, A. Ecklebe]

    :param w:
    :param l:
    :param h:

    :return:

    """
    if l <= 0:
        l = 0.0000001
    return 1 / (mu0 * (w/2/l + 2/np.pi * (1+np.log(np.pi*h/4/l))))


def sigma(l, w, R_equivalent):
    """
    1-dim fringing factor
    [according to "A Novel Approach for 3D Air Gap Reluctance Calculations" - J. Mühlethaler, J.W. Kolar, A. Ecklebe]

    :param w:
    :param l:
    :param R_equivalent:

    :return:

    """
    return R_equivalent / (l/mu0/w)


def r_round_inf(l, sigma, r):
    """
    3-dim reluctance for 2-dim axial-symmetric air gaps

    :param sigma:
    :param r:
    :param l:

    :return:

    """
    return sigma**2 * l/mu0/r**2/np.pi


def r_round_round(l, sigma, r):
    """
    3-dim reluctance for 2-dim axial-symmetric air gaps

    :param sigma:
    :param r:
    :param l:

    :return:

    """
    return sigma**2 * l/mu0/r**2/np.pi


def r_cyl_cyl(l, sigma, w, r_o):
    """

    :param l:
    :param sigma:
    :param w:
    :param r_o:

    :return:

    """
    return sigma * np.log(r_o/(r_o-l)) / 2/mu0/np.pi/w


def r_cyl_cyl_real(l, sigma, w, r_o, h_real_core):
    """

    :param l:
    :param sigma:
    :param w:
    :param r_o:

    :return:

    """
    return sigma * np.log(r_o/(r_o-l)) / mu0 / (2*np.pi - 4*np.arccos(h_real_core/2/r_o)) / w


def r_cheap_cyl_cyl(r_o, l, w):
    """

    :param r_o:
    :param l:
    :param w:

    :return:

    """
    r_i = r_o - l
    return (r_o-r_i) / mu0/w/np.pi/(r_o+r_i)


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

if __name__ == '__main__':
    pass
