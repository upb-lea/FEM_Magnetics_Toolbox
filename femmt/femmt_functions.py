# Usual Python libraries
import numpy.typing as npt
import numpy as np
import pathlib
from matplotlib import pyplot as plt
import json
import random
import string
from typing import Union
import pkg_resources
import subprocess
import sys
import os
import pandas as pd
import pathlib
import time


# Static Functions
#  Used somewhere in the Code of femmt.py
def install_femm_if_missing():
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
    :param a:
    :param b:
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
    :param a:
    :param b:
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


def call_for_path(destination, config_file="config.json"):
    """
    asks the user to enter the filepath of a destinated file WITHOUT the suffix
    stores a the filepath as a python string declaration in a config file
    returns the filepath
    :param destination:
    :param config_file:
    :return:
    """
    # pickle_file = open(config_file, "w")
    # path = input(f"Please enter the parent folder path of {destination} in ways of 'C:.../onelab-Windows64/': ")
    # pickle.dumps(path, pickle_file) # f"{destination} = '{path}'\n")
    # pickle_file.close()

    # Find out the path of installed module, or in case of running directly from git, find the path of git repository
    module_file_path = pathlib.Path(__file__).parent.absolute()

    path = input(f"Please enter the parent folder path of {destination} in ways of 'C:.../onelab-Windows64/': ")
    dict = {"onelab": path}
    file = open(module_file_path / config_file, 'w', encoding='utf-8')
    json.dump(dict, file, ensure_ascii=False)
    file.close()

    return path


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
    specified number of strands (n_strands). CAUTION: Zero number of
    layers corresponds to a single strand.
    :param n_strands:
    :return:
    """
    return np.sqrt(0.25+(n_strands-1)/3)-0.5


def fft(period_vector_t_i: npt.ArrayLike, sample_factor: float = 1000, plot: str = 'no', rad: str ='no',
        f0: Union[float, None]=None, title: str='ffT') -> npt.NDArray[list]:
    """
    A fft for a input signal. Input signal is in vector format and should include one period.

    Output vector includes only frequencies with amplitudes > 1% of input signal

    Minimal example:
    example_waveform = np.array([[0, 1.34, 3.14, 4.48, 6.28],[-175.69, 103.47, 175.69, -103.47,-175.69]])
    out = fft(example_waveform, plot=True, rad='yes', f0=25000, title='ffT input current')

    :param period_vector_t_i: numpy-array [[time-vector[,[current-vector]]. One period only
    :param sample_factor: f_sampling/f_period, defaults to 1000
    :param plot: insert anything else than "no" or 'False' to show a plot to visualize input and output
    :param rad: 'no' for time domain input vector, anything else than 'no' for 2pi-time domain
    :param f0: set when rad != 'no' and rad != False
    :param title: plot window title, defaults to 'ffT'
    :return: numpy-array [[frequency-vector],[amplitude-vector],[phase-vector]]
    """

    t = period_vector_t_i[0]
    i = period_vector_t_i[1]

    if rad != 'no' and rad!=False:
        if f0 is None:
            raise ValueError("if rad!='no', a fundamental frequency f0 must be set")
        else:
            period_vector_t_i[0] = period_vector_t_i[0] / (2 * np.pi * f0)

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
    for count, value in enumerate(x_mag_corrected):
        if x_mag_corrected[count] > 0.01 * max(i):
            f_out.append(f_corrected[count])
            x_out.append(x_mag_corrected[count])
            phi_rad_out.append(phi_rad_corrected[count])

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

        fig, [ax1, ax2] = plt.subplots(num=title, nrows=2, ncols=1)
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
        plt.tight_layout()
        plt.show()

    return np.array([f_out, x_out, phi_rad_out])


def compare_fft_list(list: list, rad: float = 'no', f0: Union[float,None] = None) -> None:
    """
    generate fft curves from input curves and compare them to each other

    minimal example:
    example_waveform = np.array([[0, 1.34, 3.14, 4.48, 6.28],[-175.69, 103.47, 175.69, -103.47,-175.69]])
    example_waveform_2 = np.array([[0, 0.55, 3.14, 3.69, 6.28],[-138.37, 257.58, 138.37, -257.58, -138.37]])
    compare_fft_list([example_waveform, example_waveform_2], rad='yes', f0=25000)

    :param list: list of fft-compatible numpy-arrays [element, element, ... ], each element format like [[time-vector[,[current-vector]]. One period only
    :param rad: 'no' for time domain input vector, anything else than 'no' for 2pi-time domain
    :param f0: set when rad != 'no'
    :return: plot
    """

    out = []
    for count, value in enumerate(list):
        out.append([fft(list[count], sample_factor=1000, plot='no', rad=rad, f0=f0)])

    fig, axs = plt.subplots(2, len(list), sharey=True)
    for count, value in enumerate(list):
        axs[0, count].plot(list[count][0], list[count][1], label='original signal')
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


def r_cheap_cyl_cyl(r_o, l, w):
    """

    :param r_o:
    :param l:
    :param w:
    :return:
    """
    r_i = r_o - l
    return (r_o-r_i) / mu0/w/np.pi/(r_o+r_i)

