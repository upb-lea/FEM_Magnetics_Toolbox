"""Functions to calculate reluctance models."""
# python libraries
import logging

# femmt libraries
from femmt.constants import *
import femmt.functions as ff

# 3rd party libraries
import numpy as np
import scipy
from matplotlib import pyplot as plt

logger = logging.getLogger(__name__)

def calculate_ls_lh_n_from_inductance_matrix(inductance_matrix: float | np.ndarray):
    """
    Calculate the transformer primary concentrated circuit parameters from matrix.

    :param inductance_matrix: input reluctance matrix in form of [[L_11, M], [M, L_22]] in H
    :type inductance_matrix: float | np.ndarray
    :return l_s: primary concentrated stray inductance in H
    :return l_h: primary concentrated main inductance in H
    :return n: ratio
    """
    l_11 = inductance_matrix[0, 0]
    l_22 = inductance_matrix[1, 1]

    mutal_inductance = inductance_matrix[1, 0]

    coupling_factor = mutal_inductance / (np.sqrt(l_11 * l_22))
    n = mutal_inductance / l_22
    l_s = (1 - coupling_factor ** 2) * l_11
    l_h = coupling_factor ** 2 * l_11

    return l_s, l_h, n


def calculate_inductance_matrix_from_ls_lh_n(l_s_target_value: float | np.ndarray, l_h_target_value: float | np.ndarray,
                                             n_target_value: float | np.ndarray) -> float | np.ndarray:
    """
    Calculate the inductance matrix from ls, lh, n parameters.

    :param l_s_target_value: serial inductance in H
    :type l_s_target_value: float | np.ndarray
    :param l_h_target_value: mutual inductance in H
    :type l_h_target_value: float | np.ndarray
    :param n_target_value: transfer ratio
    :type n_target_value: float | np.ndarray
    :return: inductance matrix in H
    :rtype: float | np.ndarray
    """
    inductance_matrix = [
        [l_s_target_value + l_h_target_value, l_h_target_value / n_target_value],
        [l_h_target_value / n_target_value, l_h_target_value / (n_target_value ** 2)]]
    return inductance_matrix


def power_losses_hysteresis_cylinder_radial_direction_mu_r_imag(
        flux, cylinder_height, cylinder_inner_radius, cylinder_outer_radius,
        fundamental_frequency, mu_r_abs, flux_density_data_vec, mu_r_imag_data_vec):
    """
    Calculate the hysteresis losses inside a cylinder, where the flux flows in radial direction.

    :param flux: flux in Wb
    :type flux: float | np.ndarray
    :param cylinder_height: cylinder height in m
    :type cylinder_height: float | np.ndarray
    :param cylinder_inner_radius: cylinder inner radius in m
    :type cylinder_inner_radius: float | np.ndarray
    :param cylinder_outer_radius: cylinder outer radius in m
    :type cylinder_outer_radius: float | np.ndarray
    :param fundamental_frequency: fundamental frequency in Hz
    :type fundamental_frequency: float | np.ndarray
    :param mu_r_imag_data_vec: imaginary part of u_r as vector
    :type mu_r_imag_data_vec: list | np.ndarray
    :param mu_r_abs: absolute value of mu_r: abs(mu_r)
    :type mu_r_abs: float | np.ndarray
    :param flux_density_data_vec: flux-density data vector in T
    :type flux_density_data_vec: list | np.ndarray
    """

    def flux_density_cylinder_envelope(cylinder_radius: float | np.ndarray, flux_in_cylinder: float | np.ndarray,
                                       height_of_cylinder: float | np.ndarray):
        """
        Helper-function, what is used as a function to integrate by scipy.integrate.quad.

        It calculates the flux density in a cylinder envelope. By using the integration function, the flux density
        in a volume can be calculated, as done in the superordinate function.

        Note: function parameter names differ from outer parameters to avoid 'shadows name from outer scope'.

        :param cylinder_radius: cylinder radius
        :type cylinder_radius: float | np.ndarray
        :param flux_in_cylinder: flux trough cylinder envelope depending on its radius
        :type flux_in_cylinder: float | np.ndarray
        :param height_of_cylinder: cylinder height
        :type height_of_cylinder: float | np.ndarray
        """
        return flux_in_cylinder / (2 * np.pi * cylinder_radius * height_of_cylinder)

    def power_loss_density_cylinder_envelope(cylinder_radius: float | np.ndarray, flux_in_cylinder: float | np.ndarray,
                                             height_of_cylinder: float | np.ndarray) -> float | np.ndarray:
        """
        Helper-function, what is used as a function to integrate by scipy.integrate.quad.

        It calculates the power losses in a cylinder envelope. Together with the integration function, the hysteresis
        losses are calculated in a volumetric cylinder where the flux is orientated in radiant direction.

        volume = cylinder_diameter * pi * cylinder_height

        Basic formula for hysteresis losses:
        p_hyst = volume * 0.5 * omega * u_imag * (max(flux_density) / mu_0 / mu_r) ** 2
        -> p_hyst = volume * pi * frequency * u_imag * (max(flux_density) / mu_0 / mu_r) ** 2

        Note: The [0] in the return for only handling the result itself to the output, not the error.

        Note: function parameter names differ from outer parameters to avoid 'shadows name from outer scope'.

        :param cylinder_radius: cylinder radius in m
        :type cylinder_radius: float | np.ndarray
        :param flux_in_cylinder: flux in Wb inside the cylinder
        :type flux_in_cylinder: float | np.ndarray
        :param height_of_cylinder: height of the cylinder in m
        :type height_of_cylinder: float | np.ndarray
        """
        mu_r_imag = np.interp(flux_density_cylinder_envelope(cylinder_radius, flux_in_cylinder, height_of_cylinder),
                              flux_density_data_vec, mu_r_imag_data_vec)

        return 2 * np.pi * cylinder_radius * height_of_cylinder * np.pi * fundamental_frequency * mu_0 * mu_r_imag * \
            (flux_density_cylinder_envelope(cylinder_radius,
                                            flux_in_cylinder, height_of_cylinder) / mu_r_abs / mu_0) ** 2

    return scipy.integrate.quad(power_loss_density_cylinder_envelope, cylinder_inner_radius, cylinder_outer_radius,
                                args=(flux, cylinder_height), epsabs=1e-4)[0]


def hyst_losses_core_half_mu_r_imag(core_inner_diameter: float | np.ndarray, window_h_half: float | np.ndarray, window_w: float | np.ndarray,
                                    mu_r_abs: float | np.ndarray, flux_max: float | np.ndarray, fundamental_frequency: float | np.ndarray,
                                    flux_density_data_vec: list | np.ndarray, mu_r_imag_data_vec: list | np.ndarray) -> float | np.ndarray:
    """
    Calculate the losses of a core cylinder half.

    means: losses in inner cylinder + losses in outer cylinder + losses in ONE cylinder_radial e.g. for top or bot

    Note: To calculate the hysteresis losses of an inductor, you need to run this function twice
    with each the half window_h

    :param core_inner_diameter: core inner diameter in m
    :type core_inner_diameter: float | np.ndarray
    :param window_h_half: window height of the core-half to consider in m
    :type window_h_half: float | np.ndarray
    :param window_w: window width in m
    :type window_w: float | np.ndarray
    :param mu_r_abs: (absolute) mu_r value from core material datasheet
    :type mu_r_abs: float | np.ndarray
    :param mu_r_imag_data_vec: imaginary part of mu_r as data vector
    :type mu_r_imag_data_vec: list | np.ndarray
    :param flux_max: maximum flux in Wb what appears in the core-half
    :type flux_max: float | np.ndarray
    :param fundamental_frequency: fundamental frequency of flux in Hz
    :type fundamental_frequency: float | np.ndarray
    :param flux_density_data_vec: flux density as data vector
    :type flux_density_data_vec: float | np.ndarray
    :return: hysteresis losses of the core half in W
    :rtype: float | np.ndarray
    """
    core_cross_section = (core_inner_diameter / 2) ** 2 * np.pi
    flux_density_max = flux_max / core_cross_section
    volume_cylinder_inner = core_cross_section * window_h_half

    losses_cylinder_inner = power_loss_hysteresis_simple_volume_mu_r_imag(fundamental_frequency, flux_density_max,
                                                                          mu_r_abs, volume_cylinder_inner,
                                                                          flux_density_data_vec, mu_r_imag_data_vec)

    cylinder_inner_radius = core_inner_diameter / 2
    cylinder_outer_radius = core_inner_diameter / 2 + window_w
    losses_cylinder_radial = power_losses_hysteresis_cylinder_radial_direction_mu_r_imag(
        flux_max, core_inner_diameter / 4, cylinder_inner_radius, cylinder_outer_radius,
        fundamental_frequency, mu_r_abs, flux_density_data_vec, mu_r_imag_data_vec)
    return 2 * losses_cylinder_inner + losses_cylinder_radial


def calculate_core_2daxi_total_volume(core_inner_diameter: float | np.ndarray, window_h: float | np.ndarray, window_w: float | np.ndarray):
    """
    Calculate the full volume of a rotationally symmetric core.

    Inside material (windings, air) also belong to this volume.
    This is the total volume used by the magnetic itself.

    :param core_inner_diameter: core inner diameter
    :type core_inner_diameter: float | np.ndarray
    :param window_h: winding window height
    :type window_h: float | np.ndarray
    :param window_w: winding window width
    :type window_w: float | np.ndarray
    """
    outer_core_radius = calculate_r_outer(core_inner_diameter, window_w)

    core_2daxi_total_volume = outer_core_radius ** 2 * np.pi * (window_h + core_inner_diameter / 2)

    return core_2daxi_total_volume


def calculate_r_outer(core_inner_diameter: float | np.ndarray, window_w: float | np.ndarray,
                      outer_core_cross_section_scale: float | np.ndarray = 1.0) -> float | np.ndarray:
    """
    Calculate outer core radius.

    Default assumption: outer core cross-section is same as inner core cross-section.

    :param outer_core_cross_section_scale: scales the outer legs cross-section relative to the center leg cross-section
    :type outer_core_cross_section_scale: float | np.ndarray
    :param core_inner_diameter: core inner diameter in m
    :type core_inner_diameter: float | np.ndarray
    :param window_w: width of core window in m
    :type window_w: float | np.ndarray
    """
    outer_core_radius = np.sqrt(outer_core_cross_section_scale * (core_inner_diameter/2)**2 + (core_inner_diameter / 2 + window_w)**2)
    return outer_core_radius


def power_losses_hysteresis_cylinder_radial_direction(
        flux: float | np.ndarray, cylinder_height: float | np.ndarray, cylinder_inner_radius: float | np.ndarray,
        cylinder_outer_radius: float | np.ndarray, fundamental_frequency: float | np.ndarray, mu_r_imag: float | np.ndarray,
        mu_r_abs: float | np.ndarray) -> float | np.ndarray:
    """
    Calculate the hysteresis losses inside a cylinder, where the flux flows in radial direction.

    :param flux: flux in Wb
    :type flux: float | np.ndarray
    :param cylinder_height: cylinder height in m
    :type cylinder_height: float | np.ndarray
    :param cylinder_inner_radius: cylinder inner radius in m
    :type cylinder_inner_radius: float | np.ndarray
    :param cylinder_outer_radius: cylinder outer radius in m
    :type cylinder_outer_radius: float | np.ndarray
    :param fundamental_frequency: fundamental frequency in Hz
    :type fundamental_frequency: float | np.ndarray
    :param mu_r_imag: imaginary part of u_r
    :type mu_r_imag: float | np.ndarray
    :param mu_r_abs: absolute value of mu_r: abs(mu_r)
    :type mu_r_abs: float | np.ndarray
    :return: Hysteresis losses in W of cylinder with flux in radial direction
    :rtype: float | np.ndarray
    """
    def flux_density_cylinder_envelope(cylinder_radius: float | np.ndarray, flux_in_cylinder: float | np.ndarray,
                                       height_of_cylinder: float | np.ndarray) -> float | np.ndarray:
        """
        Helper-function, what is used as a function to integrate by scipy.integrate.quad.

        It calculates the flux density in a cylinder envelope. By using the integration function, the flux density
        in a volume can be calculated, as done in the superordinate function.

        Note: function parameter names differ from outer parameters to avoid 'shadows name from outer scope'.

        :param cylinder_radius: cylinder radius in m
        :type cylinder_radius: float | np.ndarray
        :param flux_in_cylinder: flux in Wb trough cylinder envelope depending on its radius
        :type flux_in_cylinder: float | np.ndarray
        :param height_of_cylinder: cylinder height in m
        :type height_of_cylinder: float | np.ndarray
        :return: Flux density in T
        :rtype: float | np.ndarray
        """
        return flux_in_cylinder / (2 * np.pi * cylinder_radius * height_of_cylinder)

    def power_loss_density_cylinder_envelope(cylinder_radius: float | np.ndarray, flux_in_cylinder: float | np.ndarray,
                                             height_of_cylinder: float | np.ndarray) -> float | np.ndarray:
        """
        Helper-function, what is used as a function to integrate by scipy.integrate.quad.

        It calculates the power losses in a cylinder envelope. Together with the integration function, the hysteresis
        losses are calculated in a volumetric cylinder where the flux is orientated in radiant direction.

        volume = cylinder_diameter * pi * cylinder_height

        Basic formula for hysteresis losses:
        p_hyst = volume * 0.5 * omega * u_imag * (max(flux_density) / mu_0 / mu_r) ** 2
        -> p_hyst = volume * pi * frequency * u_imag * (max(flux_density) / mu_0 / mu_r) ** 2

        Note: this is the simplified version to calculate the losses using a fixed mu_r_imag.
        In case of using the u_r_imag in dependence of the different flux density, what varies inside the cylinder,
        use this as a hint to modify the formula.
        AnalyticalCoreData.f_N95_mu_imag(fundamental_frequency, flux_density_cylinder_envelope(cylinder_radius,
         flux, cylinder_height))

        The [0] in the return for only handling the result itself to the output, not the error.

        Note: function parameter names differ from outer parameters to avoid 'shadows name from outer scope'.

        :param cylinder_radius: cylinder radius in m
        :type cylinder_radius: float | np.ndarray
        :param flux_in_cylinder: Flux in cylinder in Wb
        :type flux_in_cylinder: float | np.ndarray
        :param height_of_cylinder: height of cylinder in m
        :type height_of_cylinder: float | np.ndarray
        :return: Power loss density in W/m³
        :rtype: float | np.ndarray
        """
        return 2 * np.pi * cylinder_radius * height_of_cylinder * np.pi * fundamental_frequency * mu_0 * mu_r_imag * \
            (flux_density_cylinder_envelope(cylinder_radius,
                                            flux_in_cylinder, height_of_cylinder) / mu_r_abs / mu_0) ** 2

    return scipy.integrate.quad(power_loss_density_cylinder_envelope, cylinder_inner_radius, cylinder_outer_radius,
                                args=(flux, cylinder_height), epsabs=1e-4)[0]


def hyst_losses_core_half(core_inner_diameter: float | np.ndarray, window_h_half: float | np.ndarray, window_w: float | np.ndarray,
                          mu_r_imag: float | np.ndarray, mu_r_abs: float | np.ndarray, flux_max: float | np.ndarray,
                          fundamental_frequency: float | np.ndarray) -> float | np.ndarray:
    """
    Calculate the losses of a core cylinder half.

    means: losses in inner cylinder + losses in outer cylinder + losses in ONE cylinder_radial e.g. for top or bot

    Note: To calculate the hysteresis losses of an inductor, you need to run this
    function twice with each the half window_h

    :param core_inner_diameter: core inner diameter in m
    :type core_inner_diameter: float | np.ndarray
    :param window_h_half: window height in m of the core-half to consider
    :type window_h_half: float | np.ndarray
    :param window_w: window width in m
    :type window_w: float | np.ndarray
    :param mu_r_abs: (absolute) mu_r value from core material datasheet
    :type mu_r_abs: float | np.ndarray
    :param mu_r_imag: imaginary part of mu_r
    :type mu_r_imag: float | np.ndarray
    :param flux_max: maximum flux in Wb what appears in the core-half
    :type flux_max: float | np.ndarray
    :param fundamental_frequency: fundamental frequency of flux in Hz
    :type fundamental_frequency: float | np.ndarray
    :return: Hysteresis losses in W of the core half
    :rtype: float | np.ndarray
    """
    core_cross_section = (core_inner_diameter / 2) ** 2 * np.pi
    flux_density_max = flux_max / core_cross_section
    volume_cylinder_inner = core_cross_section * window_h_half
    losses_cylinder_inner = power_loss_hysteresis_simple_volume(fundamental_frequency, mu_r_imag,
                                                                flux_density_max, mu_r_abs, volume_cylinder_inner)

    cylinder_inner_radius = core_inner_diameter / 2
    cylinder_outer_radius = core_inner_diameter / 2 + window_w
    losses_cylinder_radial = power_losses_hysteresis_cylinder_radial_direction(
        flux_max, core_inner_diameter / 4, cylinder_inner_radius, cylinder_outer_radius,
        fundamental_frequency, mu_r_imag, mu_r_abs)
    return 2 * losses_cylinder_inner + losses_cylinder_radial


def calculate_reluctance_matrix(winding_matrix: np.array, inductance_matrix: np.array) -> np.array:
    """
    Calculate the reluctance matrix. Everything must be numpy.

    L. Keuck, "Entwurf eines einstufigen Ladewandlers auf Basis eines LLC-Resonanzwandlers", dissertation 2023

    :param winding_matrix: winding matrix
    :type winding_matrix: np.array
    :param inductance_matrix: inductance matrix
    :type inductance_matrix: np.array
    :return: reluctance[-matrix]
    :rtype: np.array

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


def calculate_inductance_matrix(reluctance_matrix: float | np.ndarray, winding_matrix: float | np.ndarray) -> float | np.ndarray:
    """
    Calculate the inductance matrix out of reluctance matrix and winding matrix.

    :param reluctance_matrix: matrix of transformer reluctance
    :type reluctance_matrix: float | np.ndarray
    :param winding_matrix: matrix of transformer windings
    :type winding_matrix: float | np.ndarray
    :return: inductance matrix in H
    :rtype: float | np.ndarray


    reluctance matrix e.g.
    r = [ [], [] ]

    winding matrix e.g.
    N = [ [N_1a, N_2b], [N_1b, N_2b] ]

    returns inductance matrix e.g.
    L = [ [L_11, M], [M, L_22] ]
    """
    if np.ndim(reluctance_matrix) == 0:
        reluctance_matrix_invert = 1 / reluctance_matrix
    else:
        reluctance_matrix_invert = np.linalg.inv(reluctance_matrix)

    return np.matmul(np.matmul(np.transpose(winding_matrix), reluctance_matrix_invert), winding_matrix)


def calculate_flux_matrix(reluctance_matrix: float | np.ndarray, winding_matrix: float | np.ndarray,
                          current_matrix: float | np.ndarray) -> float | np.ndarray:
    """
    Calculate the flux for e.g. an integrated transformer.

    reluctance matrix e.g.
    reluctance_matrix = [ [r1, r2], [r3, r4] ]

    winding matrix e.g.
    winding_matrix = [ [N_1a, N_2b], [N_1b, N_2b] ]

    current matrix e.g.
    current_matrix = [ [current_1], [current_2] ]

    returns flux matrix e.g.
    flux_matrix = [ [flux_1], [flux_2] ]

    :param reluctance_matrix: matrix of transformer reluctance
    :type reluctance_matrix: float | np.ndarray
    :param winding_matrix: matrix of transformer windings
    :type winding_matrix: float | np.ndarray
    :return: inductance matrix in H
    :rtype: float | np.ndarray
    :param current_matrix: matrix of currents in A
    :type current_matrix: float | np.ndarray
    """
    if np.ndim(reluctance_matrix) == 0:
        reluctance_matrix_invert = 1 / reluctance_matrix
    else:
        reluctance_matrix_invert = np.linalg.inv(reluctance_matrix)

    return np.matmul(np.matmul(reluctance_matrix_invert, winding_matrix), current_matrix)


def time_vec_current_vec_from_time_current_vec(time_current_vec: list[list[float]] | np.ndarray):
    """
    Split a time-current vector into time and current vector.

    :param time_current_vec: [time_vector, current_vector]
    :type time_current_vec: list[list[float]]
    """
    return time_current_vec[0], time_current_vec[1]


def flux_vec_from_current_vec(current_vec_1, current_vec_2, winding_matrix, inductance_matrix):
    """
    Calculate the integrated transformer flux from current vectors.

    :param current_vec_2: current vector e.g. [[time1, time2, ...], [current1, current2, ...]]
    :type current_vec_2: np.array
    :param current_vec_1: current vector e.g. [[time1, time2, ...], [current1, current2, ...]]
    :type current_vec_1: np.array
    :param winding_matrix: winding matrix e.g. [[... ], [...]] in shape (2,2)
    :type winding_matrix: np.array
    :param inductance_matrix: inductance matrix e.g. [[... ], [...]] in shape (2,2)
    :type inductance_matrix: np.array
    """
    flux_top_vec = []
    flux_bot_vec = []
    flux_stray_vec = []

    for count, _ in enumerate(current_vec_1):
        current_value_timestep = [current_vec_1[count], current_vec_2[count]]

        # simplified formula: flux = L * I / N
        [flux_top_timestep, flux_bot_timestep] = np.matmul(np.matmul(np.linalg.inv(np.transpose(winding_matrix)),
                                                                     inductance_matrix),
                                                           np.transpose(current_value_timestep))
        flux_stray_timestep = flux_bot_timestep - flux_top_timestep

        flux_top_vec.append(flux_top_timestep)
        flux_bot_vec.append(flux_bot_timestep)
        flux_stray_vec.append(flux_stray_timestep)

    return flux_top_vec, flux_bot_vec, flux_stray_vec


def visualize_current_and_flux(time, flux_top_vec, flux_bot_vec, flux_stray_vec, current_1_vec, current_2_vec):
    """
    Visualize current and flux over time for the integrated transformer.

    :param time: time vector
    :type time: list | np.ndarray
    :param flux_top_vec: flux vector in top core part
    :type flux_top_vec: list | np.ndarray
    :param flux_bot_vec: flux vector in bottom core part
    :type flux_bot_vec: list | np.ndarray
    :param flux_stray_vec: flux vector in stray core part
    :type flux_stray_vec: list | np.ndarray
    :param current_1_vec: current vector of winding 1
    :type current_1_vec: list | np.ndarray
    :param current_2_vec: current vector of winding 2
    :type current_2_vec: list | np.ndarray

    """
    figure, axis = plt.subplots(2, figsize=(4, 4))

    axis[0].plot(time, current_1_vec, label=r"$I_{\mathrm{in}}$")
    axis[0].plot(time, -np.array(current_2_vec), label=r"$I_{\mathrm{out}}$")

    axis[1].plot(time, 1000 * np.array(flux_top_vec), label=r"$\mathcal{\phi}_{\mathrm{top}}$")
    axis[1].plot(time, 1000 * np.array(flux_bot_vec), label=r"$\mathcal{\phi}_{\mathrm{bot}}$")
    axis[1].plot(time, 1000 * np.array(flux_stray_vec), label=r"$\mathcal{\phi}_{\mathrm{stray}}$")

    axis[1].set_ylabel("Magnetic fluxes / mWb")
    axis[1].set_xlabel(r"$t$ / $\mathrm{\mu s}$")
    axis[0].set_ylabel("Currents / A")
    axis[0].set_xlabel(r"$t$ / $\mathrm{\mu s}$")

    axis[0].legend()
    axis[1].legend()
    axis[0].grid()
    axis[1].grid()
    plt.show()


def max_value_from_value_vec(*args):
    """
    Return the peak values from the vectors.

    :param args: value_vector
    :return: peak_value_from_vector
    """
    peak_list = []
    for value_vector in args:
        peak = max([abs(value) for value in value_vector])
        peak_list.append(peak)

    return tuple(peak_list)


def phases_deg_from_time_current(time_vec: list, *args):
    """
    Return the phases_deg of the peaks.

    :param time_vec: time vector with time steps
    :type time_vec: list
    :param args: vectors of current
    """
    period = time_vec[-1]
    phases_deg = []

    for current_vec in args:
        time_max = time_vec[np.array(current_vec).argmax(axis=0)]
        phase = time_max / period * 360
        phases_deg.append(phase)

    return tuple(phases_deg)


def power_loss_hysteresis_simple_volume(fundamental_frequency: float | np.ndarray, mu_r_imag: float | np.ndarray,
                                        flux_density_max: float | np.ndarray,
                                        mu_r_abs: float | np.ndarray, core_volume: float | np.ndarray) -> float | np.ndarray:
    """
    Calculate the hysteresis losses depending on the input parameters.

    The output are the losses for a certain volume of core.

    :param fundamental_frequency: fundamental frequency in Hz
    :type fundamental_frequency: float | np.ndarray
    :param mu_r_imag: imaginary part of u_r
    :type mu_r_imag: float | np.ndarray
    :param flux_density_max: maximum flux density
    :type flux_density_max: float | np.ndarray
    :param mu_r_abs: abs(mu_r)
    :type mu_r_abs: float | np.ndarray
    :param core_volume: core volume in m³
    :type core_volume: float | np.ndarray
    :return: Hysteresis losses in W for a given volume
    :rtype: float | np.ndarray
    """
    return core_volume * np.pi * fundamental_frequency * mu_r_imag * mu_0 * (flux_density_max / mu_0 / mu_r_abs) ** 2


def power_loss_hysteresis_simple_volume_mu_r_imag(fundamental_frequency: float | np.ndarray, flux_density_max: float | np.ndarray,
                                                  mu_r_abs: float | np.ndarray, core_volume: float | np.ndarray,
                                                  flux_density_data_vec: float | np.ndarray, mu_r_imag_data_vec: float | np.ndarray):
    """
    Calculate the hysteresis losses depending on the input parameters.

    The output are the losses for a certain volume of core.

    :param fundamental_frequency: fundamental frequency in Hz
    :type fundamental_frequency: float | np.ndarray
    :param mu_r_imag_data_vec: imaginary part of u_r as data vector
    :type mu_r_imag_data_vec: float | np.ndarray
    :param flux_density_max: maximum flux density
    :type flux_density_max: float | np.ndarray
    :param mu_r_abs: abs(mu_r)
    :type mu_r_abs: float | np.ndarray
    :param core_volume: core volume
    :type core_volume: float | np.ndarray
    :param flux_density_data_vec: flux density as data input vector
    :type flux_density_data_vec: float | np.ndarray
    """
    mu_r_imag = np.interp(flux_density_max, flux_density_data_vec, mu_r_imag_data_vec)

    return core_volume * np.pi * fundamental_frequency * mu_r_imag * mu_0 * (flux_density_max / mu_0 / mu_r_abs) ** 2


def r_basic_round_inf(air_gap_radius: float | np.ndarray, air_gap_basic_height: float | np.ndarray,
                      core_height: float | np.ndarray) -> float | np.ndarray:
    """
    Calculate the r_basic for a round to infinite structure. Do not use this function directly.

    Use it indirectly by using
     - r_air_gap_round_round
     - r_air_gap_round_inf
    instead!

    This function calculates the r_basic for a round to infinite structure according to the following paper:
    [according to "A Novel Approach for 3D Air Gap Reluctance Calculations" - J. Mühlethaler, J.W. Kolar, A. Ecklebe]

    :param air_gap_radius: air gap radius in m
    :type air_gap_radius: float | np.ndarray
    :param air_gap_basic_height: air gap height in m for the BASIC-AIR-GAP (e.g. if you use a round-round structure,
        this is half of the total air gap).
    :type air_gap_basic_height: float | np.ndarray
    :param core_height: core height in m
    :type core_height: float | np.ndarray
    :return: basic reluctance for round - infinite structure
    :rtype: float | np.ndarray
    """
    conductance_basic = mu_0 * (air_gap_radius * 2 / 2 / air_gap_basic_height + 2 / np.pi * \
                                (1 + np.log(np.pi * core_height / 4 / air_gap_basic_height)))

    return 1 / conductance_basic


def sigma_round(r_equivalent: float | np.ndarray, air_gap_radius: float | np.ndarray, air_gap_total_height: float | np.ndarray)\
        -> float | np.ndarray:
    """
    Calculate sigma for a round structure. Do not use this function directly.

    Use it indirectly by using
     - r_air_gap_round_round
     - r_air_gap_round_inf
    instead!

    :param r_equivalent: this is a series/parallel connection of r_basic, depending on the air gap structure
    :type r_equivalent: float | np.ndarray
    :param air_gap_radius: air gap radius in m
    :type air_gap_radius: float | np.ndarray
    :param air_gap_total_height: air gap total height in m (for the total air gap, also for round-round structures)
    :type air_gap_total_height: float | np.ndarray
    :return: fringing factor 'sigma'
    :rtype: float | np.ndarray
    """
    return r_equivalent * mu_0 * air_gap_radius / air_gap_total_height


def r_air_gap_round_round(air_gap_total_height: float | np.ndarray, core_inner_diameter: float | np.ndarray, core_height_upper: float | np.ndarray,
                          core_height_lower: float | np.ndarray) -> float | np.ndarray:
    """
    Return the reluctance of a round-round air gap structure and includes fringing effects.

    :param air_gap_total_height: total air gap height of the air gap
    :type air_gap_total_height: float | np.ndarray
    :param core_inner_diameter: core inner diameter
    :type core_inner_diameter: float | np.ndarray
    :param core_height_upper: core height upper (needed for better calculating fringing effects)
    :type core_height_upper: float | np.ndarray
    :param core_height_lower: core height lower (needed for better calculating fringing effects)
    :type core_height_lower: float | np.ndarray
    :return: air gap reluctance for round-round structure including fringing effects
    :rtype: float | np.ndarray
    """
    if np.any(air_gap_total_height == 0):
        raise ValueError("'air_gap_total_height' can not be Zero!")

    air_gap_total_height = np.array(air_gap_total_height)
    core_inner_diameter = np.array(core_inner_diameter)
    core_height_upper = np.array(core_height_upper)
    core_height_lower = np.array(core_height_lower)
    air_gap_radius = core_inner_diameter / 2

    air_gap_basic_height = air_gap_total_height / 2
    r_basic_upper = r_basic_round_inf(air_gap_radius, air_gap_basic_height, core_height_upper)
    r_basic_lower = r_basic_round_inf(air_gap_radius, air_gap_basic_height, core_height_lower)

    r_equivalent_round_round = r_basic_upper + r_basic_lower

    sigma = sigma_round(r_equivalent_round_round, air_gap_radius, air_gap_total_height)
    if np.any(sigma > 1):
        raise Exception(f"Failure in calculating reluctance. Sigma was calculated to >1. Check input parameters!\n"
                        f"\t{air_gap_total_height=}\n"
                        f"\t{core_inner_diameter=}\n"
                        f"\t{core_height_upper=}\n"
                        f"\t{core_height_lower=}\n"
                        f"\t{r_equivalent_round_round=}\n"
                        f"\t{air_gap_radius=}")

    r_air_gap_ideal = air_gap_total_height / mu_0 / np.pi / (air_gap_radius ** 2)
    r_air_gap = sigma ** 2 * r_air_gap_ideal

    return r_air_gap


def r_air_gap_round_round_sct(air_gap_total_height: float | np.ndarray, core_inner_diameter: float | np.ndarray,
                              core_height_upper: float | np.ndarray, core_height_lower: float | np.ndarray,
                              target_reluctance: float | np.ndarray) -> float | np.ndarray:
    """
    Calculate the air gap length of a round-round structure by a given target reluctance.

    :param air_gap_total_height: total height in m of all air gaps
    :type air_gap_total_height: float | np.ndarray
    :param core_inner_diameter: core inner diameter in m
    :type core_inner_diameter: float | np.ndarray
    :param target_reluctance: target reluctance
    :type target_reluctance: float | np.ndarray
    :param core_height_upper: height of upper core part in m
    :type core_height_upper: float | np.ndarray
    :param core_height_lower: height of lower core part in m
    :type core_height_lower: float | np.ndarray
    """
    return r_air_gap_round_round(air_gap_total_height, core_inner_diameter, core_height_upper, core_height_lower) - \
        target_reluctance


def r_air_gap_round_inf(air_gap_total_height: float | np.ndarray, core_inner_diameter: float | np.ndarray, core_height: float | np.ndarray):
    """
    Return the reluctance of a round-infinite air gap structure and includes fringing effects.

    :param air_gap_total_height: total air gap height of the air gap in m
    :type air_gap_total_height: float | np.ndarray
    :param core_inner_diameter: core inner diameter in m
    :type core_inner_diameter: float | np.ndarray
    :param core_height: core height in m (needed for better calculating fringing effects)
    :type core_height: float | np.ndarray
    :return: air gap reluctance for round-inf structure including fringing effects
    :rtype: float | np.ndarray
    """
    air_gap_total_height = np.array(air_gap_total_height)
    core_inner_diameter = np.array(core_inner_diameter)
    core_height = np.array(core_height)

    air_gap_radius = core_inner_diameter / 2
    r_basic = r_basic_round_inf(air_gap_radius, air_gap_total_height, core_height)

    r_equivalent_round_inf = r_basic
    sigma = sigma_round(r_equivalent_round_inf, air_gap_radius, air_gap_total_height)

    r_air_gap_ideal = air_gap_total_height / mu_0 / np.pi / (air_gap_radius ** 2)
    r_air_gap = sigma ** 2 * r_air_gap_ideal

    return r_air_gap


def r_air_gap_round_inf_sct(air_gap_total_height: float | np.ndarray, core_inner_diameter: float | np.ndarray,
                            core_height: float | np.ndarray, target_reluctance: float | np.ndarray) -> float | np.ndarray:
    """
    Calculate the air gap length of a round infinite structure by a given target reluctance.

    :param air_gap_total_height: total height in m of all air gaps
    :type air_gap_total_height: float | np.ndarray
    :param core_inner_diameter: core inner diameter in m
    :type core_inner_diameter: float | np.ndarray
    :param target_reluctance: target reluctance
    :type target_reluctance: float | np.ndarray
    :param core_height: height of core in m
    :type core_height: float | np.ndarray
    """
    return r_air_gap_round_inf(air_gap_total_height, core_inner_diameter, core_height) - target_reluctance


def r_basic_tablet_cyl(tablet_height: float | np.ndarray, air_gap_basic_height: float | np.ndarray,
                       tablet_radius: float | np.ndarray) -> float | np.ndarray:
    """
    Calculate the r_basic for round to infinite structure. Do not use this function directly.

    Use it indirectly by using
     - r_air_gap_tablet_cyl
    instead!

    This function calculates the r_basic for a round to infinite structure according to the following paper:
    [according to "A Novel Approach for 3D Air Gap Reluctance Calculations" - J. Mühlethaler, J.W. Kolar, A. Ecklebe]

    Note: this is the same function as r_basic_round_inf, but with clear variable names for tablet-cylinder structure

    :param tablet_height: tablet height in m = air gap width for tablet-cylinder structure
    :type tablet_height: float | np.ndarray
    :param air_gap_basic_height: air gap height in m for the BASIC-AIR-GAP (e.g. if you use a round-round structure, this
        is half of the total air gap).
    :type air_gap_basic_height: float | np.ndarray
    :param tablet_radius: tablet radius in m
    :type tablet_radius: float | np.ndarray
    :return: basic reluctance for tablet in m - cylinder structure
    :rtype: float | np.ndarray
    """
    if air_gap_basic_height == 0:
        raise ZeroDivisionError(f"Division by zero: {air_gap_basic_height=}")

    conductance_basic = mu_0 * (tablet_height / 2 / air_gap_basic_height + 2 / \
                                np.pi * (1 + np.log(np.pi * tablet_radius / 4 / air_gap_basic_height)))

    return 1 / conductance_basic


def sigma_tablet_cyl(r_equivalent: float | np.ndarray, tablet_height: float | np.ndarray, air_gap_total_height: float | np.ndarray):
    """
    Do not use this function directly! Calculate sigma for a tablet-cylinder structure.

    Use it indirectly by using
     - r_air_gap_tablet_cyl
    instead!

    Note: this is the same function as sigma_round, but with clear variable names for tablet-cylinder structure

    :param r_equivalent: this is a series/parallel connection of r_basic, depending on the air gap structure
    :type r_equivalent: float | np.ndarray
    :param tablet_height: tablet height
    :type tablet_height: float | np.ndarray
    :param air_gap_total_height: air gap total height (for the total air gap)
    :type air_gap_total_height: float | np.ndarray
    :return: fringing factor 'sigma' for tablet - cylinder structure
    :rtype: float | np.ndarray
    """
    return r_equivalent * mu_0 * tablet_height / air_gap_total_height


def r_air_gap_tablet_cyl(tablet_height: float | np.ndarray, air_gap_total_height: float | np.ndarray, core_inner_diameter: float | np.ndarray,
                         window_w: float | np.ndarray) -> float | np.ndarray:
    """
    Return the reluctance of a cylinder-tablet air gap structure and includes fringing effects.

    This function calculates the air gap reluctance for a 2D-axisymmetric core.

    :param tablet_height: tablet height in m
    :type tablet_height: float | np.ndarray
    :param air_gap_total_height: total air gap height in m
    :type air_gap_total_height: float | np.ndarray
    :param core_inner_diameter: core inner diameter in m
    :type core_inner_diameter: float | np.ndarray
    :param window_w: core window width in m
    :type window_w: float | np.ndarray
    :return: air gap reluctance for tablet - cylinder structure including air gap fringing
    :rtype: float | np.ndarray
    """
    r_inner = core_inner_diameter / 2 + window_w

    # translate practical core dimensions to non-practical air-gap dimensions
    tablet_radius = r_inner - air_gap_total_height

    air_gap_basic_height = air_gap_total_height
    r_basic = r_basic_tablet_cyl(tablet_height, air_gap_basic_height, tablet_radius)

    r_equivalent = r_basic / 2
    sigma = sigma_tablet_cyl(r_equivalent, tablet_height, air_gap_total_height)
    if np.any(sigma > 1):
        raise Exception("Failure in calculating reluctance. Sigma was calculated to >1. Check input parameters!")

    r_air_gap_ideal = np.log(r_inner / (r_inner - air_gap_total_height)) / 2 / mu_0 / np.pi / tablet_height

    r_air_gap = sigma * r_air_gap_ideal

    return r_air_gap


def r_air_gap_tablet_cylinder_sct(air_gap_total_height: float | np.ndarray, core_inner_diameter: float | np.ndarray,
                                  tablet_height: float | np.ndarray, window_w: float | np.ndarray,
                                  target_reluctance: float | np.ndarray) -> float | np.ndarray:
    """
    Calculate the air gap length of a table cylinder by a given target reluctance.

    :param air_gap_total_height: total height in m of all air gaps
    :type air_gap_total_height: float | np.ndarray
    :param core_inner_diameter: core inner diameter in m
    :type core_inner_diameter: float | np.ndarray
    :param tablet_height: height of tablet in m
    :type tablet_height: float | np.ndarray
    :param window_w: window width in m
    :type window_w: float | np.ndarray
    :param target_reluctance: target reluctance
    :type target_reluctance: float | np.ndarray
    """
    return r_air_gap_tablet_cyl(tablet_height, air_gap_total_height, core_inner_diameter, window_w) - target_reluctance


def r_air_gap_tablet_cyl_no_2d_axi(tablet_height, air_gap_total_length, core_inner_diameter, window_w):
    """
    Return the reluctance of a cylinder-tablet air gap structure and includes fringing effects.

    Note:
    ----
    This function differs from r_air_gap_tablet_cyl (ideal 2D axisymmetric core). Here, the air gap reluctance for
    a non-2D-axisymmetric core is taken into account, as a real PQ core is open at the side. So, there is no air gap
    taken into account for the side-sections. The new core_dimension_y parameter describes the width of the
    core when you are in a xy-coordinate system.

    :param tablet_height: tablet height in m
    :type tablet_height: float | np.ndarray
    :param air_gap_total_length: air gap total length
    :type air_gap_total_length: float | np.ndarray
    :param core_inner_diameter: core inner diameter in m
    :type core_inner_diameter: float | np.ndarray
    :param window_w: core window width in m
    :type window_w: float | np.ndarray
    :return: air gap reluctance for tablet - cylinder structure including air gap fringing
    :rtype: float | np.ndarray
    """
    r_inner = core_inner_diameter / 2 + window_w

    if np.any(air_gap_total_length >= window_w):
        raise Exception("air_gap_total_height is greater than window_w")

    air_gap_basic_height = air_gap_total_length
    r_basic = r_basic_tablet_cyl(tablet_height, air_gap_basic_height, (core_inner_diameter + 2 * window_w - 2 * \
                                                                       air_gap_total_length) / 2)

    r_equivalent = r_basic / 2
    sigma = sigma_tablet_cyl(r_equivalent, tablet_height, air_gap_total_length)
    if np.any(sigma > 1):
        raise Exception("Failure in calculating reluctance. Sigma was calculated to >1. Check input parameters!")

    # Note:
    # the circumference of the air gap differs for open cores (e.g. PQ40/40) to closed ones
    # (ideal rotationally symmetric)
    # The circumference is no more (diameter * pi), but (2 * pi - 4 * alpha) * (core_inner_diameter/2 + window_w),
    # with alpha = arccos(core_dimension_x / (core_inner_diameter + 2 * window_w))
    # See: Dissertation Lukas Keuck
    # For equal pq core sizes (e.g. PQ 40/40), it has been found out that
    #     core_dimension_x / core_dimension_y = 1.45, the error over all available shapes is maximum 7%
    #     (compared to datasheet value)
    # Now, the new and partly circumference of the stray-path air gap can be calculated
    # First, the core dimension_y needs to be calculated.
    # Formular 1: core_dimension_x / core_dimension_y = 1.45
    # Formular 2: core_dimension_x * core_dimension_y - (core_inner_diameter / 2 + window_w) ** 2 * np.pi
    #     = (core_inner_diameter / 2 ) ** 2 * np.pi
    # Formular 2 assumes that the outer core cross-section of the core is equal to the inner core cross-section
    # Formular 1 & 2 needs to be solved to get core_dimension_y:

    core_dimension_y = np.sqrt((core_inner_diameter ** 2 / 4 + \
                                (core_inner_diameter / 2 + window_w) ** 2) * np.pi / 1.45)
    r_air_gap_ideal_partly = np.log(r_inner / (r_inner - air_gap_total_length)) / \
        mu_0 / (2 * np.pi - 4 * np.arccos(core_dimension_y / 2 / r_inner)) / tablet_height

    r_air_gap = sigma * r_air_gap_ideal_partly

    return r_air_gap


def r_core_tablet(tablet_height: float | np.ndarray, tablet_radius: float | np.ndarray, mu_r_abs: float | np.ndarray,
                  core_inner_diameter: float | np.ndarray):
    """
    Calculate the magnetic resistance of the core tablet.

    :param tablet_height: tablet height in m
    :type tablet_height: float | np.ndarray
    :param tablet_radius: tablet radius in m
    :type tablet_radius: float | np.ndarray
    :param mu_r_abs: relative permeability (mu_r) of the core material from datasheet
    :type mu_r_abs: float | np.ndarray
    :param core_inner_diameter: core inner diameter. For idealized core material, this value can be 0.001.
    :type core_inner_diameter: float | np.ndarray
    """
    return np.log(tablet_radius / (core_inner_diameter / 2)) / (2 * np.pi * mu_0 * mu_r_abs * tablet_height)

def r_core_tablet_2(tablet_height, tablet_radius, mu_r_abs):
    """
     Calculate the magnetic resistance (reluctance) of the core tablet.

    :param tablet_height: The height of the core tablet.
    :type tablet_height: float
    :param tablet_radius: The radius of the core tablet.
    :type tablet_radius: float
    :param mu_r_abs: The absolute relative permeability (mu_r) of the core material.
    :type mu_r_abs: float
    :return: The magnetic reluctance of the core tablet.
    :rtype: float
    """
    return tablet_height / (np.pi * mu_0 * mu_r_abs * tablet_radius ** 2)


def r_core_top_bot_radiant(core_inner_diameter: float | np.ndarray, window_w: float | np.ndarray, mu_r_abs: float | np.ndarray,
                           core_top_bot_height: float | np.ndarray) -> float | np.ndarray:
    """
    Calculate the top or bottom core material part.

    :param core_inner_diameter: core inner diameter
    :type core_inner_diameter: float | np.ndarray
    :param window_w: width of winding window
    :type window_w: float | np.ndarray
    :param mu_r_abs: relative permeability (mu_r) of the core material from datasheet
    :type mu_r_abs: float | np.ndarray
    :param core_top_bot_height: height of the core material top / bottom of the winding window
    :type core_top_bot_height: float | np.ndarray
    """
    return np.log((core_inner_diameter + 2 * window_w) / core_inner_diameter) / \
        (2 * np.pi * mu_0 * mu_r_abs * core_top_bot_height)


def r_core_round(core_inner_diameter: float | np.ndarray, core_round_height: float | np.ndarray, mu_r_abs: float | np.ndarray):
    """
    Calculate the core reluctance for a round structure.

    :param core_round_height: height of the round core part section
    :type core_round_height: float | np.ndarray
    :param core_inner_diameter: core inner diameter
    :type core_inner_diameter: float | np.ndarray
    :param mu_r_abs: relative permeability (mu_r) of the core material from datasheet
    :type mu_r_abs: float | np.ndarray
    """
    return core_round_height / (mu_0 * mu_r_abs * (core_inner_diameter / 2) ** 2 * np.pi)


def resistance_solid_wire(core_inner_diameter: float, window_w: float, turns_count: int, conductor_radius: float,
                          material: str = 'Copper', temperature: float = 100) -> float:
    """
    Calculate the resistance of a solid wire, considering the winding temperature and a middle turn length.

    :param core_inner_diameter: core inner diameter
    :type core_inner_diameter: float
    :param window_w: winding window width
    :type window_w: float
    :param turns_count: number of turns
    :type turns_count: int
    :param conductor_radius: conductor radius
    :type conductor_radius: float
    :param material: Material, e.g. "Copper" or "Aluminum"
    :type material: str
    :param temperature: temperature in °C
    :type temperature: float
    :return: total resistance of wire
    :rtype: float
    """
    # simplification: middle turn length
    # figure out middle length of one turn for given geometry
    turn_radius = core_inner_diameter / 2 + conductor_radius
    turn_count = 1
    total_turn_length = 0.0
    while turn_radius <= (core_inner_diameter / 2 + window_w - conductor_radius):
        total_turn_length += turn_radius * 2 * np.pi
        turn_radius += 2 * conductor_radius
        turn_count += 1
    middle_turn_length = total_turn_length / turn_count

    total_turn_length = turns_count * middle_turn_length

    sigma_copper = ff.conductivity_temperature(material, temperature)

    # return R = rho * l / A
    return total_turn_length / (conductor_radius ** 2 * np.pi) / sigma_copper

def resistance_litz_wire(core_inner_diameter: float, window_w: float, window_h: float, turns_count: int,
                         iso_core_top: float, iso_core_bot: float, iso_core_left: float, iso_core_right: float,
                         iso_primary_to_primary: float, litz_wire_name: str, material: str = "Copper",
                         scheme: str = "vertical_first", temperature: float = 100) -> float:
    """
    Detailed DC resistance calculation for a given number of litz wire turns in a given core window.

    :param core_inner_diameter: core inner diameter in meter
    :type core_inner_diameter: float
    :param window_w: winding window width in meter
    :type window_w: float
    :param turns_count: number of turns
    :type turns_count: int
    :param material: Material, e.g. "Copper" or "Aluminum"
    :type material: str
    :param temperature: temperature in °C
    :type temperature: float
    :param window_h: window height in meter
    :type window_h: float
    :param iso_core_top: insulation (bobbin) from winding window to top core in meter
    :type iso_core_top: float
    :param iso_core_bot: insulation (bobbin) from winding window to bot core in meter
    :type iso_core_bot: float
    :param iso_core_left: insulation (bobbin) from winding window to left core in meter
    :type iso_core_left: float
    :param iso_core_right: insulation (bobbin) from winding window to right core in meter
    :type iso_core_right: float
    :param iso_primary_to_primary: wire-to-wire insulation in meter
    :type iso_primary_to_primary: float
    :param litz_wire_name: litz wire name, e.g. "1.5x105x0.1"
    :type litz_wire_name: str
    :param scheme: scheme, e.g. "vertical_first" (finish vertical wires before going to next layer) or "horizontal_first"
    :type scheme: str
    :return: total resistance of wire in Ohm
    :rtype: float
    """
    litz_wire = ff.litz_database()[litz_wire_name]
    litz_wire_diameter = 2 * litz_wire["conductor_radii"]
    litz_wire_effective_area = litz_wire["strands_numbers"] * litz_wire["strand_radii"] ** 2 * np.pi

    # calculate total 2D-axi symmetric volume of the core:
    # formula: number_turns_per_row = (available_width + primary_to_primary) / (wire_diameter + primary_to_primary)
    available_width = window_w - iso_core_left - iso_core_right
    possible_number_turns_per_row = int(
        (available_width + iso_primary_to_primary) / (litz_wire_diameter + iso_primary_to_primary))
    available_height = window_h - iso_core_top - iso_core_bot
    possible_number_turns_per_column = int(
        (available_height + iso_primary_to_primary) / (litz_wire_diameter + iso_primary_to_primary))

    if scheme == "vertical_first":
        number_of_columns = np.ceil(turns_count / possible_number_turns_per_column)
        length_row_per_turn = np.zeros(int(number_of_columns))

        # figure out the turn length per row
        r_inner = np.array([core_inner_diameter / 2 + iso_core_left])
        middle_radius_per_column = np.array([])
        for count, _ in enumerate(length_row_per_turn):
            middle_radius_per_column = np.append(middle_radius_per_column, r_inner + count * iso_primary_to_primary + (count + 0.5) * litz_wire_diameter)
        turn_length_per_column = 2 * middle_radius_per_column * np.pi  # diameter * pi

        # figure out the windings per column
        windings_per_column = possible_number_turns_per_column * np.ones(int(number_of_columns))
        last_column_turns = np.mod(turns_count, possible_number_turns_per_column)
        windings_per_column[-1] = last_column_turns

        # get the total turn length
        total_turn_length = np.sum(np.multiply(windings_per_column, turn_length_per_column))
    elif scheme == "horizontal_first":
        number_of_columns = (window_w - iso_core_left - iso_core_right + iso_primary_to_primary) / (litz_wire_diameter + iso_primary_to_primary)
        length_row_per_turn = np.zeros(int(number_of_columns))

        # figure out the turn length per row
        r_inner = np.array([core_inner_diameter / 2 + iso_core_left])
        middle_radius_per_column = np.array([])
        for count, _ in enumerate(length_row_per_turn):
            middle_radius_per_column = np.append(middle_radius_per_column, r_inner + count * iso_primary_to_primary + (count + 0.5) * litz_wire_diameter)
        turn_length_per_column = 2 * middle_radius_per_column * np.pi  # diameter * pi

        # figure out the windings per row
        number_of_rows = np.ceil(turns_count / possible_number_turns_per_row)
        windings_per_row = possible_number_turns_per_row * np.ones(int(number_of_rows))
        last_row_turns = np.mod(turns_count, possible_number_turns_per_row)
        windings_per_row[-1] = last_row_turns

        # get the total turn length
        total_turn_length = 0
        for windings_in_row in windings_per_row:
            logger.info(f"{windings_in_row=}")
            total_turn_length += np.sum(turn_length_per_column[:int(windings_in_row)])
    else:
        raise ValueError(f"{scheme} not defined. Must be 'horizontal_first' or 'vertical_first'.")

    sigma_copper = ff.conductivity_temperature(material, temperature)
    # return R = rho * l / A
    return total_turn_length / litz_wire_effective_area / sigma_copper

def i_rms(time_current_matrix: np.array) -> float:
    """
    RMS calculation from a time-current-vector.

    :param time_current_matrix: time and current in format [[0, 0.5e-6, 2.5e-6, 3e-6], [16.55, -10.55, -16.55, 10.55]]
    :type time_current_matrix: np.array
    :return: rms current
    :rtype: float
    """
    time = time_current_matrix[0]
    current = time_current_matrix[1]

    square_integral_sum = 0

    # set up function
    for count in np.arange(np.shape(time_current_matrix)[1] - 1):
        # figure out linear equation between points
        y_axis = current[count]
        delta_time = time[count + 1] - time[count]
        delta_current = current[count + 1] - current[count]
        gradient = delta_current / delta_time

        # calculate solution of (partly) square integral
        square_integral = gradient ** 2 / 3 * delta_time ** 3 + gradient * delta_time ** 2 * \
            y_axis + y_axis ** 2 * delta_time

        # add (partly) square integral to total integration sum
        square_integral_sum += square_integral

    # return "mean" and "root" to finalize rms calculation
    return np.sqrt(square_integral_sum / time[-1])

def calc_skin_depth(frequency: float, material_name: str = "Copper", temperature: float = 100) -> float:
    """
    Calculate the skin depth for the given material, frequency and temperature.

    :param frequency: Frequency in Hz
    :type frequency: float
    :param material_name: material name, e.g. 'Copper' or 'Aluminum'
    :type material_name: str
    :param temperature: Temperature in °C
    :type temperature: float
    :return: Skin depth in meter
    :rtype: float
    """
    rho = 1 / ff.conductivity_temperature(material_name, temperature)
    skin_depth = np.sqrt(rho / np.pi / frequency / mu_0)
    return skin_depth


def calc_proximity_factor(litz_wire_name: str, number_turns: int, window_h: float, iso_core_top: float, iso_core_bot: float,
                          frequency: float, litz_wire_material_name: str = 'Copper', temperature: float = 100) -> float:
    """
    Calculate the proximity factor for cores without air gaps, according to Charles R. Sullivan: "Simpliﬁed Design Method for Litz Wire", 2014 APEC.

    :param litz_wire_name: litz wire name from the litz_database(), e.g. "1.5x105x0.1", "1.4x200x0.071"
    :type litz_wire_name: str
    :param number_turns: number of turns
    :type number_turns: int
    :param window_h: window height in meter
    :type window_h: float
    :param iso_core_top: insulation (bobbin) from winding to top core
    :type iso_core_top: float
    :param iso_core_bot: insulation (bobbin) from winding to bot core
    :type iso_core_bot: float
    :param frequency: Frequency in Hz
    :type frequency: float
    :param litz_wire_material_name: material name, e.g. 'Copper' or 'Aluminum'
    :type litz_wire_material_name: str
    :param temperature: Temperature in °C
    :type temperature: float
    """
    litz_wire = ff.litz_database()[litz_wire_name]
    nominator = (np.pi * litz_wire["strands_numbers"] * number_turns) ** 2 * (litz_wire["strand_radii"] * 2) ** 6
    b = window_h - iso_core_top - iso_core_bot
    skin_depth = calc_skin_depth(frequency, litz_wire_material_name, temperature)
    denominator = 192 * skin_depth ** 4 * b ** 2

    proximity_factor = 1 + nominator / denominator
    return proximity_factor

def calc_proximity_factor_air_gap(litz_wire_name: str, number_turns: int, r_1: float, frequency: float, winding_area: float,
                                  litz_wire_material_name: str = 'Copper', temperature: float = 100) -> float:
    """
    Calculate the proximity factor for cores with center air gap, according to Charles R. Sullivan: "Simpliﬁed Design Method for Litz Wire", 2014 APEC.

    :param litz_wire_name: litz wire name from the litz_database(), e.g. "1.5x105x0.1", "1.4x200x0.071"
    :type litz_wire_name: str
    :param number_turns: number of turns
    :type number_turns: int
    :param r_1: Radius from air gap to winding in meter (see mentioned paper for a clear drawing)
    :type r_1: float
    :param winding_area: winding area in square meters (m²). Used to calculate r_2 from the mentioned paper.
    :type winding_area: float
    :param frequency: Frequency in Hz
    :type frequency: float
    :param litz_wire_material_name: material name, e.g. 'Copper' or 'Aluminum'
    :type litz_wire_material_name: str
    :param temperature: Temperature in °C
    :type temperature: float
    """
    litz_wire = ff.litz_database()[litz_wire_name]

    nominator = (np.pi * litz_wire["strands_numbers"] * number_turns) ** 2 * (litz_wire["strand_radii"] * 2) ** 6

    # single air gap
    r_2 = np.sqrt(2 * winding_area / np.pi + r_1 ** 2)

    b_eff = np.pi * (0.693 * r_1 + 0.307 * r_2 ** 0.91 * r_1 ** 0.09)

    skin_depth = calc_skin_depth(frequency, litz_wire_material_name, temperature)

    denominator = 192 * skin_depth ** 4 * b_eff ** 2

    proximity_factor = 1 + nominator / denominator
    return proximity_factor
