import femmt
from femmt import MagneticComponent
import numpy as np
import itertools
import re
import csv
from matplotlib import pyplot as plt

plt.figure(figsize=(6, 3))


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


# directory = "C:/Users/tillp/OneDrive/Documents/GitHub/FEM_Magnetics_Toolbox/data/materials/N95/N95/100C"
directory = "C:/Users/Aniket/sciebo/Exchange_FEMMT/05_Materials/data/Loss_Data_PHD_Keuck/mu_r_Plot"
frequencies_in_kHz = [300]  # [100, 200, 300]
temperatures = [30]

for j, temperature in enumerate(temperatures):
    for i, frequency_in_kHz in enumerate(frequencies_in_kHz):
        # file_name = f"mu_phi_{frequency_in_kHz}kHz_N95_{temperature}C.txt"
        file_name = f"mu_r_{frequency_in_kHz}kHz_N95_{temperature}C.txt"
        # Get raw data from measurements
        with open(directory + "/" + file_name, newline='\n') as file:
            phi = [re.split(r'\t+', line[0]) for line in list(csv.reader(file))]

        phi_raw = np.array([[float(i) for i in j] for j in phi])

        print(phi_raw[:, 1])

        # Fit polynomial factors
        poly_degree = 1
        z = np.polyfit(phi_raw[:, 0], phi_raw[:, 1], poly_degree)
        z = np.flip(z, 0)
        print(z)

        # Construct fitted function
        b = phi_raw[:, 0]
        phi = PolyCoefficients(x=b, coeffs=z)
        mu_imag = np.sin(phi*np.pi/180) * 3000

        # Extrapolate
        b_extrapolated = np.concatenate((b, b+b.max(0)), axis=0)
        # Add zero value pair
        b_extrapolated = np.concatenate((np.array([0]), b_extrapolated), axis=0)
        # print(f"{b_extrapolated=}")

        phi_extrapolated = PolyCoefficients(x=b_extrapolated, coeffs=z)
        np.insert(phi_extrapolated, 0, 0.99*phi_extrapolated[0], axis=0)

        # print(f"{phi_extrapolated=}")
        mu_imag_extrapolated = np.sin(phi_extrapolated*np.pi/180) * 3000

        # Saving
        np.savetxt(directory + "/ONELAB_READY_B_" + file_name, b_extrapolated, newline=', ')
        np.savetxt(directory + "/ONELAB_READY_MU_IMAG_" + file_name, mu_imag_extrapolated, newline=', ')

        if j == 0:
            col = "r"
        else:
            col = "b"

        # Visualize
        plt.plot(1000*np.array(phi_raw[:, 0]), np.sin(phi_raw[:, 1]*np.pi/180) * 3000, color=col, label=f"{frequency_in_kHz} kHz, {temperature} C,raw data")
        # plt.plot(b, mu_imag, label=f"{frequency_in_kHz}")
        plt.plot(1000*np.array(b_extrapolated), mu_imag_extrapolated, "--", color=col, label=f"{frequency_in_kHz} kHz, {temperature} C, curve fitted")

plt.ylabel("$\mu''/\mu_0$")
plt.xlabel("$B$ / mT")
plt.legend()
plt.grid()
plt.savefig("C:/Users/Aniket/sciebo/Exchange_FEMMT/04_Documentation//mu_imag.pdf", bbox_inches="tight")
plt.show()
