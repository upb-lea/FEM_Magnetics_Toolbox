import femmt
from femmt import MagneticComponent
import numpy as np
import itertools
import re
import os
import csv
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d


def get_data(directory: str, file_name: str):
    with open(directory + "/" + file_name, newline='\n') as file:
        raw_data = [re.split(r'\t+', line[0]) for line in list(csv.reader(file))]

    return np.array([[float(i) for i in j] for j in raw_data])


def crop_data(raw_data, lower_bound: float, upper_bound: float):
    data = []
    for dat in raw_data:
        if (lower_bound < dat[0] < upper_bound):
            data.append(dat)
    return np.array(data)


def load_and_extract(parameter: str, read_directory: str, write_directory: str, temperatures: list[int], frequencies: list[int], material: str, save_data: bool = False, do_plot: bool = False, save_plot: bool = False):

    if do_plot:
        plt.figure(figsize=(5, 2.5))

    for temperature in temperatures:
        for frequency in frequencies:
            file_name = f"{parameter}_{int(frequency/1000)}kHz_{material}_{temperature}C.txt"

            print(temperature, frequency)

            raw_data = get_data(directory=read_directory, file_name=file_name)
            data = crop_data(raw_data=raw_data, lower_bound=0.03, upper_bound=0.25)

            b = data[:, 0]
            mur = data[:, 1]

            if not os.path.isdir(write_directory):
                os.mkdir(write_directory)
            np.savetxt(write_directory + f"b_{int(frequency/1000)}kHz_{material}_{temperature}C.csv", b, newline=', ', fmt='%f')
            np.savetxt(write_directory + f"{parameter}_{int(frequency/1000)}kHz_{material}_{temperature}C.csv", mur, newline=', ', fmt='%f')

            if do_plot:
                plt.plot(1000 * data[:, 0], data[:, 1], "--", label=f"{frequency}, {temperature} °C")

    if do_plot:
        if parameter=="mu_r":
            plt.ylabel(r"$\mu_\mathrm{r}  /  \mu_0$")
            # plt.axhline(y=3000, color='k', linestyle='-', label=r"inital $\mu_\mathrm{r}$")
        if parameter=="mu_phi":
            plt.ylabel(r"$\zeta_\mathrm{\mu}$")
        plt.xlabel("$B$ / mT")
        # plt.ylim(2600, 3600)
        plt.legend(ncol=2)
        plt.grid()
        if save_plot:
            plt.savefig("C:/Users/tillp/sciebo/Exchange_FEMMT/04_Documentation/mu_r.pdf", bbox_inches="tight")
        plt.show()


def to_arithmetic_form(amplitude: np.ndarray, phase: np.ndarray):
    real_part = amplitude * np.cos(phase)
    imaginary_part = amplitude * np.sin(phase)
    return real_part, imaginary_part


def PolyCoefficients(x, coeffs):
    """ Returns a polynomial for ``x`` values for the ``coeffs`` provided.

    The coefficients must be in ascending order (``x**0`` to ``x**o``).
    """
    o = len(coeffs)
    print(f'# This is a polynomial of order {ord}.')
    y = 0
    for i in range(o):
        y += coeffs[i]*x**i
    return y


def create_arithmetic_form(temperatures: list[int], frequencies: list[int], material: str, mu_r_frequency: int = None, poly_degree: int = 1):
    # static read-directory
    read_directory=f"materials/{material}/"
    write_directory = f"materials/{material}/"

    for temperature in temperatures:
        for frequency in frequencies:
            if mu_r_frequency is not None:
                b_mu_r = np.genfromtxt(read_directory+f"mu_r/b_{int(mu_r_frequency / 1000)}kHz_{material}_{temperature}C.csv", delimiter=',', dtype=None)
                mu_r = np.genfromtxt(read_directory+f"mu_r/mu_r_{int(mu_r_frequency / 1000)}kHz_{material}_{temperature}C.csv", delimiter=',', dtype=None)
            else:
                b_mu_r = np.genfromtxt(read_directory+f"mu_r/b_{int(frequency / 1000)}kHz_{material}_{temperature}C.csv", delimiter=',', dtype=None)
                mu_r = np.genfromtxt(read_directory+f"mu_r/mu_r_{int(frequency / 1000)}kHz_{material}_{temperature}C.csv", delimiter=',', dtype=None)

            # Sort amplitude permeability values
            arr1inds = b_mu_r.argsort()
            b_mu_r_sorted = np.flip(b_mu_r[arr1inds[::-1]])
            mu_r_sorted = np.flip(mu_r[arr1inds[::-1]])

            # Clean amplitude permeability data for way to small values (typically at b=0)
            to_delete = []
            for i, val in enumerate(mu_r_sorted):
                # if val <= 0.5*np.mean(mu_r_sorted):
                if val == 0:
                    to_delete.append(i)
            b_mu_r_sorted = np.delete(b_mu_r_sorted, to_delete)
            mu_r_sorted = np.delete(mu_r_sorted, to_delete)

            # Fit polynomial factors
            z_mur = np.polyfit(b_mu_r_sorted, mu_r_sorted, poly_degree)
            z_mur = np.flip(z_mur, 0)
            mu_r_fitted = PolyCoefficients(x=b_mu_r_sorted, coeffs=z_mur)


            plt.plot(b_mu_r_sorted, mu_r_sorted, label=f"mu_r at {mu_r_frequency}Hz")
            plt.plot(b_mu_r_sorted, mu_r_fitted, label=f"mu_r fitted at {mu_r_frequency}Hz")

            # Get permeability loss angle data
            b = np.genfromtxt(read_directory+f"mu_phi/b_{int(frequency / 1000)}kHz_{material}_{temperature}C.csv", delimiter=',', dtype=None)
            mu_phi = np.genfromtxt(read_directory+f"mu_phi/mu_phi_{int(frequency / 1000)}kHz_{material}_{temperature}C.csv", delimiter=',', dtype=None)

            # Fit the permeability loss angle data
            # Extrapolate the permeability loss angle data
            # print(temperature, frequency)
            # kernel_size = 5
            # kernel = np.ones(kernel_size) / kernel_size
            # print(np.concatenate((np.sort(mu_phi), np.sort(mu_phi)[-kernel_size:])))
            # print(np.sort(mu_phi))
            # print(np.sort(mu_phi)[-kernel_size:])
            # mu_phi_smoothed = np.convolve(np.concatenate((sorted(mu_phi), sorted(mu_phi)[-kernel_size:-1])), kernel, mode='valid')
            # plt.plot(sorted(b), mu_phi_smoothed)
            # f_mu_r = interp1d(b_mu_r, mu_r)

            # print(b, mu_phi)

            # Sort permeability loss angle data
            arr1inds = b.argsort()
            b_sorted = np.flip(b[arr1inds[::-1]])
            mu_phi_sorted = np.flip(mu_phi[arr1inds[::-1]])
            # plt.plot(b, mu_phi)
            # plt.plot(b_sorted, mu_phi_sorted)

            # Extend b
            b_full_range = np.concatenate((b_sorted, max(b_sorted)+b_sorted))

            # Fit polynomial factors
            z = np.polyfit(b_sorted, mu_phi_sorted, 1)
            # z = np.polyfit(b_sorted, mu_phi_sorted, poly_degree-1)
            z = np.flip(z, 0)
            phi_full_range = PolyCoefficients(x=b_full_range, coeffs=z)
            if not os.path.isdir(write_directory+"/phi_full_range"):
                os.mkdir(write_directory+"/phi_full_range")
            np.savetxt(write_directory + f"phi_full_range/phi_full_range_{int(frequency/1000)}kHz_{material}_{temperature}C.csv", phi_full_range, newline=', ', fmt='%f')


            # mu_r with b from current mu_phi data
            # plt.plot(b_sorted, PolyCoefficients(x=b_sorted, coeffs=z_mur))

            # plt.plot(b_sorted, phi)
            # plt.plot(b_sorted, mu_phi_sorted)

            # Calculate real and imaginary part of complex permeability
            # mu_imag_1 = np.sin(np.deg2rad(mu_phi_sorted)) * PolyCoefficients(x=b_sorted, coeffs=z_mur)
            # plt.plot(b_sorted, mu_imag_1)
            # plt.plot(b_full_range, PolyCoefficients(x=b_full_range, coeffs=z_mur), label="amplitude")

            b_max = 0.4

            mu_imag = np.sin(np.deg2rad(phi_full_range)) * PolyCoefficients(x=b_full_range, coeffs=z_mur)  # * (1+np.sign(b_max-b_full_range))/2
            plt.plot(b_full_range, mu_imag, label="imaginary part")

            mu_real

            mu_real = np.cos(np.deg2rad(phi_full_range)) * PolyCoefficients(x=b_full_range, coeffs=z_mur)  # * (1+np.sign(b_max-b_full_range))/2 + (1-np.sign(b_max-b_full_range))/2
            plt.plot(b_full_range, mu_real, label="real part")

            plt.legend()
            plt.xlabel("$B$ / T")
            plt.grid()
            plt.show()

            if not os.path.isdir(write_directory+"/mu_real"):
                os.mkdir(write_directory+"/mu_real")
            if not os.path.isdir(write_directory+"/mu_imag"):
                os.mkdir(write_directory+"/mu_imag")

            np.savetxt(write_directory + f"mu_real/mu_real_{int(frequency/1000)}kHz_{material}_{temperature}C.csv", mu_real, newline=', ', fmt='%f')
            np.savetxt(write_directory + f"mu_real/b_{int(frequency/1000)}kHz_{material}_{temperature}C.csv", b_full_range, newline=', ', fmt='%f')
            np.savetxt(write_directory + f"mu_imag/mu_imag_{int(frequency/1000)}kHz_{material}_{temperature}C.csv", mu_imag, newline=', ', fmt='%f')
            np.savetxt(write_directory + f"mu_imag/b_{int(frequency/1000)}kHz_{material}_{temperature}C.csv", b_full_range, newline=', ', fmt='%f')
