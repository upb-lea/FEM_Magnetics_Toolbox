# python libraries


# femmt libraries
from femmt.constants import *
import femmt.functions as ff

# 3rd party libraries
import numpy as np
import scipy
from matplotlib import pyplot as plt



def calculate_ls_lh_n_from_inductance_matrix(inductance_matrix):
    """
    Calculates the transformer primary concentrated circuit parameters from matrix
    :param inductance_matrix: input reluctance matrix in form of [[L_11, M], [M, L_22]]
    :return l_s: primary concentrated stray inductance
    :return l_h: primary concentrated main inductance
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



def power_losses_hysteresis_cylinder_radial_direction_mu_r_imag(flux, cylinder_height, cylinder_inner_radius, cylinder_outer_radius,
                                                                fundamental_frequency, mu_r_abs, flux_density_data_vec, mu_r_imag_data_vec):
    """
    This function calculates the hysteresis losses inside a cylinder, where the flux flows in radial direction.

    :param flux: flux
    :param cylinder_height: cylinder height
    :param cylinder_inner_radius: cylinder inner radius
    :param cylinder_outer_radius: cylinder outer radius
    :param fundamental_frequency: fundamental frequency
    :param mu_r_imag_data_vec: imaginary part of u_r as vector
    :param mu_r_abs: absolute value of mu_r: abs(mu_r)
    :param flux_density_data_vec: flux-density data vector

    """

    def flux_density_cylinder_envelope(cylinder_radius, flux_in_cylinder, height_of_cylinder):
        """
        This is a helper-function, what is used as a function to integrate by scipy.integrate.quad.
        It calculates the flux density in a cylinder envelope. By using the integration function, the flux density
        in a volume can be calculated, as done in the superordinate function.

        Note: function parameter names differ from outer parameters to avoid 'shadows name from outer scope'.

        :param cylinder_radius: cylinder radius
        :param flux_in_cylinder: flux trough cylinder envelope depending on its radius
        :param height_of_cylinder: cylinder height
        """
        return flux_in_cylinder / (2 * np.pi * cylinder_radius * height_of_cylinder)

    def power_loss_density_cylinder_envelope(cylinder_radius, flux_in_cylinder, height_of_cylinder):
        """
        This is a helper-function, what is used as a function to integrate by scipy.integrate.quad
        It calculates the power losses in a cylinder envelope. Together with the integration function, the hysteresis
        losses are calculated in a volumetric cylinder where the flux is orientated in radiant direction.

        volume = cylinder_diameter * pi * cylinder_height

        Basic formula for hysteresis losses:
        p_hyst = volume * 0.5 * omega * u_imag * (max(flux_density) / mu_0 / mu_r) ** 2
        -> p_hyst = volume * pi * frequency * u_imag * (max(flux_density) / mu_0 / mu_r) ** 2

        Note: The [0] in the return for only handling the result itself to the output, not the error.

        Note: function parameter names differ from outer parameters to avoid 'shadows name from outer scope'.
        """
        mu_r_imag = np.interp(flux_density_cylinder_envelope(cylinder_radius, flux_in_cylinder, height_of_cylinder), flux_density_data_vec, mu_r_imag_data_vec)

        return 2 * np.pi * cylinder_radius * height_of_cylinder * np.pi * fundamental_frequency * mu_0 * mu_r_imag * (flux_density_cylinder_envelope(cylinder_radius, flux_in_cylinder, height_of_cylinder) / mu_r_abs / mu_0) ** 2

    return scipy.integrate.quad(power_loss_density_cylinder_envelope, cylinder_inner_radius, cylinder_outer_radius, args=(flux, cylinder_height),
                                epsabs=1e-4)[0]

def hyst_losses_core_half_mu_r_imag(core_inner_diameter, window_h_half, window_w, mu_r_abs, flux_max, fundamental_frequency, flux_density_data_vec, mu_r_imag_data_vec):
    """
    Calculates the losses of a core cylinder half.
    means: losses in inner cylinder + losses in outer cylinder + losses in ONE cylinder_radial e.g. for top or bot

    Note: To calculate the hysteresis losses of an inductor, you need to run this function twice with each the half window_h

    :param core_inner_diameter: core inner diameter
    :param window_h_half: window height of the core-half to consider
    :param window_w: window width
    :param mu_r_abs: (absolute) mu_r value from core material datasheet
    :param mu_r_imag_data_vec: imaginary part of mu_r as data vector
    :param flux_max: maximum flux what appears in the core-half
    :param fundamental_frequency: fundamental frequency of flux
    :param flux_density_data_vec: flux density as data vector
    """
    core_cross_section = (core_inner_diameter / 2) ** 2 * np.pi
    flux_density_max = flux_max / core_cross_section
    volume_cylinder_inner = core_cross_section * window_h_half

    losses_cylinder_inner = power_loss_hysteresis_simple_volume_mu_r_imag(fundamental_frequency, flux_density_max, mu_r_abs, volume_cylinder_inner,
                                                      flux_density_data_vec, mu_r_imag_data_vec)

    cylinder_inner_radius = core_inner_diameter / 2
    cylinder_outer_radius = core_inner_diameter / 2 + window_w
    losses_cylinder_radial = power_losses_hysteresis_cylinder_radial_direction_mu_r_imag(flux_max, core_inner_diameter / 4,
                                                                               cylinder_inner_radius, cylinder_outer_radius,
                                                                               fundamental_frequency, mu_r_abs, flux_density_data_vec, mu_r_imag_data_vec)
    return 2 * losses_cylinder_inner + losses_cylinder_radial



def calculate_core_2daxi_total_volume(core_inner_diameter, window_h, window_w):
    """
    Calculates the full volume of a rotationally symmetric core.
    Inside material (windings, air) also belong to this volume.
    This is the total volume used by the magnetic itself.

    :param core_inner_diameter: core inner diameter
    :param window_h: winding window height
    :param window_w: winding window width
    """
    outer_core_radius = np.sqrt( core_inner_diameter ** 2 / 2 + core_inner_diameter * window_w + window_w ** 2)

    return outer_core_radius ** 2 * np.pi * (window_h + core_inner_diameter / 2)

def power_losses_hysteresis_cylinder_radial_direction(flux, cylinder_height, cylinder_inner_radius, cylinder_outer_radius,
                                                      fundamental_frequency, mu_r_imag, mu_r_abs):
    """
    This function calculates the hysteresis losses inside a cylinder, where the flux flows in radial direction.

    :param flux: flux
    :param cylinder_height: cylinder height
    :param cylinder_inner_radius: cylinder inner radius
    :param cylinder_outer_radius: cylinder outer radius
    :param fundamental_frequency: fundamental frequency
    :param mu_r_imag: imaginary part of u_r
    :param mu_r_abs: absolute value of mu_r: abs(mu_r)

    """

    def flux_density_cylinder_envelope(cylinder_radius, flux_in_cylinder, height_of_cylinder):
        """
        This is a helper-function, what is used as a function to integrate by scipy.integrate.quad.
        It calculates the flux density in a cylinder envelope. By using the integration function, the flux density
        in a volume can be calculated, as done in the superordinate function.

        Note: function parameter names differ from outer parameters to avoid 'shadows name from outer scope'.

        :param cylinder_radius: cylinder radius
        :param flux_in_cylinder: flux trough cylinder envelope depending on its radius
        :param height_of_cylinder: cylinder height
        """
        return flux_in_cylinder / (2 * np.pi * cylinder_radius * height_of_cylinder)

    def power_loss_density_cylinder_envelope(cylinder_radius, flux_in_cylinder, height_of_cylinder):
        """
        This is a helper-function, what is used as a function to integrate by scipy.integrate.quad
        It calculates the power losses in a cylinder envelope. Together with the integration function, the hysteresis
        losses are calculated in a volumetric cylinder where the flux is orientated in radiant direction.

        volume = cylinder_diameter * pi * cylinder_height

        Basic formula for hysteresis losses:
        p_hyst = volume * 0.5 * omega * u_imag * (max(flux_density) / mu_0 / mu_r) ** 2
        -> p_hyst = volume * pi * frequency * u_imag * (max(flux_density) / mu_0 / mu_r) ** 2

        Note: this is the simplified version to calculate the losses using a fixed mu_r_imag.
        In case of using the u_r_imag in dependence of the different flux density, what varies inside the cylinder,
        use this as a hint to modify the formula.
        AnalyticalCoreData.f_N95_mu_imag(fundamental_frequency, flux_density_cylinder_envelope(cylinder_radius, flux, cylinder_height))

        The [0] in the return for only handling the result itself to the output, not the error.

        Note: function parameter names differ from outer parameters to avoid 'shadows name from outer scope'.
        """

        return 2 * np.pi * cylinder_radius * height_of_cylinder * np.pi * fundamental_frequency * mu_0 * mu_r_imag * (flux_density_cylinder_envelope(cylinder_radius, flux_in_cylinder, height_of_cylinder) / mu_r_abs / mu_0) ** 2

    return scipy.integrate.quad(power_loss_density_cylinder_envelope, cylinder_inner_radius, cylinder_outer_radius, args=(flux, cylinder_height),
                                epsabs=1e-4)[0]


def hyst_losses_core_half(core_inner_diameter, window_h_half, window_w, mu_r_imag, mu_r_abs, flux_max, fundamental_frequency):
    """
    Calculates the losses of a core cylinder half.
    means: losses in inner cylinder + losses in outer cylinder + losses in ONE cylinder_radial e.g. for top or bot

    Note: To calculate the hysteresis losses of an inductor, you need to run this function twice with each the half window_h

    :param core_inner_diameter: core inner diameter
    :param window_h_half: window height of the core-half to consider
    :param window_w: window width
    :param mu_r_abs: (absolute) mu_r value from core material datasheet
    :param mu_r_imag: imaginary part of mu_r
    :param flux_max: maximum flux what appears in the core-half
    :param fundamental_frequency: fundamental frequency of flux
    """
    core_cross_section = (core_inner_diameter / 2) ** 2 * np.pi
    flux_density_max = flux_max / core_cross_section
    volume_cylinder_inner = core_cross_section * window_h_half
    losses_cylinder_inner = power_loss_hysteresis_simple_volume(fundamental_frequency, mu_r_imag, flux_density_max, mu_r_abs, volume_cylinder_inner)

    cylinder_inner_radius = core_inner_diameter / 2
    cylinder_outer_radius = core_inner_diameter / 2 + window_w
    losses_cylinder_radial = power_losses_hysteresis_cylinder_radial_direction(flux_max, core_inner_diameter / 4,
                                                                               cylinder_inner_radius, cylinder_outer_radius,
                                                                               fundamental_frequency, mu_r_imag, mu_r_abs)
    return 2 * losses_cylinder_inner + losses_cylinder_radial

def calculate_reluctance_matrix(winding_matrix, inductance_matrix):
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

def calculate_inductance_matrix(reluctance_matrix, winding_matrix):
    """
    Calculates the inductance matrix out of reluctance matrix and winding matrix.

    :param reluctance_matrix: matrix of transformer reluctance
    :param winding_matrix: matrix of transformer windings
    :return: inductance matrix


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

def calculate_flux_matrix(reluctance_matrix, winding_matrix, current_matrix):
    """
    calculates the flux for e.g. an integrated transformer

    reluctance matrix e.g.
    reluctance_matrix = [ [r1, r2], [r3, r4] ]

    winding matrix e.g.
    winding_matrix = [ [N_1a, N_2b], [N_1b, N_2b] ]

    current matrix e.g.
    current_matrix = [ [current_1], [current_2] ]

    returns flux matrix e.g.
    flux_matrix = [ [flux_1], [flux_2] ]
    """


    if np.ndim(reluctance_matrix) == 0:
        reluctance_matrix_invert = 1 / reluctance_matrix
    else:
        reluctance_matrix_invert = np.linalg.inv(reluctance_matrix)

    return np.matmul(np.matmul(reluctance_matrix_invert, winding_matrix), current_matrix)


def time_vec_current_vec_from_time_current_vec(time_current_vec):
    return time_current_vec[0], time_current_vec[1]


def flux_vec_from_current_vec(current_vec_1, current_vec_2, winding_matrix, inductance_matrix):
    """
    calculates the integrated transformer flux from current vectors

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

    for count, value in enumerate(current_vec_1):
        # Negative sign is placed here
        current_value_timestep = [current_vec_1[count], -current_vec_2[count]]

        # simplified formula: flux = L * I / N
        [flux_top_timestep, flux_bot_timestep] = np.matmul(np.matmul(np.linalg.inv(np.transpose(winding_matrix)), inductance_matrix),
                                                           np.transpose(current_value_timestep))
        flux_stray_timestep = flux_bot_timestep - flux_top_timestep

        flux_top_vec.append(flux_top_timestep)
        flux_bot_vec.append(flux_bot_timestep)
        flux_stray_vec.append(flux_stray_timestep)

    return flux_top_vec, flux_bot_vec, flux_stray_vec


def visualize_current_and_flux(time, flux_top_vec, flux_bot_vec, flux_stray_vec, current_1_vec, current_2_vec):
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
    Returns the peak values from the vectors

    :param args: vector
    :return: peak_value_from_vector
    """
    peak_list = []
    for vector in args:
        peak = max([abs(timestep) for timestep in vector])
        peak_list.append(peak)

    return tuple(peak_list)


def phases_deg_from_time_current(time_vec, *args):
    """
    Returns the phases_deg of the peaks

    :param time_vec: time vector with time steps
    :param args: vectors of current
    """
    period = time_vec[-1]
    phases_deg = []

    for current_vec in args:
        time_max = time_vec[np.array(current_vec).argmax(axis=0)]
        phase = time_max / period * 360
        phases_deg.append(phase)

    return tuple(phases_deg)


def power_loss_hysteresis_simple_volume(fundamental_frequency, mu_r_imag, flux_density_max, mu_r_abs, core_volume):
    """
    Calculates the hysteresis losses depending on the input parameters.
    The output are the losses for a certain volume of core.

    :param fundamental_frequency: fundamental frequency in Hz
    :param mu_r_imag: imaginary part of u_r
    :param flux_density_max: maximum flux density
    :param mu_r_abs: abs(mu_r)
    :param core_volume: core volume
    """

    return core_volume * np.pi * fundamental_frequency * mu_r_imag * mu_0 * (flux_density_max / mu_0 / mu_r_abs) ** 2



def power_loss_hysteresis_simple_volume_mu_r_imag(fundamental_frequency, flux_density_max, mu_r_abs, core_volume, flux_density_data_vec, mu_r_imag_data_vec):
    """
    Calculates the hysteresis losses depending on the input parameters.
    The output are the losses for a certain volume of core.

    :param fundamental_frequency: fundamental frequency in Hz
    :param mu_r_imag_data_vec: imaginary part of u_r as data vector
    :param flux_density_max: maximum flux density
    :param mu_r_abs: abs(mu_r)
    :param core_volume: core volume
    :param flux_density_data_vec: flux density as data input vector
    """
    mu_r_imag = np.interp(flux_density_max, flux_density_data_vec, mu_r_imag_data_vec)

    return core_volume * np.pi * fundamental_frequency * mu_r_imag * mu_0 * (flux_density_max / mu_0 / mu_r_abs) ** 2


def r_basic_round_inf(air_gap_radius, air_gap_basic_height, core_height):
    """
    Do not use this function directly!
    Use it indirectly by using
     - r_air_gap_round_round
     - r_air_gap_round_inf
    instead!

    This function calculates the r_basic for a round to infinite structure according to the following paper:
    [according to "A Novel Approach for 3D Air Gap Reluctance Calculations" - J. Mühlethaler, J.W. Kolar, A. Ecklebe]

    :param air_gap_radius: air gap radius
    :param air_gap_basic_height: air gap height for the BASIC-AIR-GAP (e.g. if you use a round-round structure, this is half of the total air gap).
    :param core_height: core height
    :return: basic reluctance for round - infinite structure
    """
    conductance_basic = mu_0 * (air_gap_radius * 2 / 2 / air_gap_basic_height + 2 / np.pi * (1 + np.log(np.pi * core_height / 4 / air_gap_basic_height)))

    return 1 / conductance_basic

def sigma_round(r_equivalent, air_gap_radius, air_gap_total_height):
    """
    Do not use this function directly!
    Use it indirectly by using
     - r_air_gap_round_round
     - r_air_gap_round_inf
    instead!

    :param r_equivalent: this is a series/parallel connection of r_basic, depending on the air gap structure
    :param air_gap_radius: air gap radius
    :param air_gap_total_height: air gap total height (for the total air gap, also for round-round structures)
    :return: fringing factor 'sigma'
    """
    return r_equivalent * mu_0 * air_gap_radius / air_gap_total_height

def r_air_gap_round_round(air_gap_total_height, core_inner_diameter, core_height_upper, core_height_lower):
    """
    Returns the reluctance of a round-round air gap structure and includes finging effects.

    :param air_gap_total_height: total air gap height of the air gap
    :param core_inner_diameter: core inner diameter
    :param core_height_upper: core height upper (needed for better calculating fringing effects)
    :param core_height_lower: core height lower (needed for better calculating fringing effects)
    :return: air gap reluctance for round-round structure including fringing effects
    """
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
        raise Exception("Failure in calculating reluctance. Sigma was calculated to >1. Check input parameters!")

    r_air_gap_ideal = air_gap_total_height / mu_0 / np.pi / (air_gap_radius ** 2)
    r_air_gap = sigma ** 2 * r_air_gap_ideal

    return r_air_gap

def r_air_gap_round_round_sct(air_gap_total_height, core_inner_diameter, core_height_upper, core_height_lower, target_reluctance):
    return r_air_gap_round_round(air_gap_total_height, core_inner_diameter, core_height_upper, core_height_lower) - target_reluctance

def r_air_gap_round_inf(air_gap_total_height, core_inner_diameter, core_height):
    """
    Returns the reluctance of a round-infinite air gap structure and includes fringing effects

    :param air_gap_total_height: total air gap height of the air gap
    :param core_inner_diameter: core inner diameter
    :param core_height: core height (needed for better calculating fringing effects)
    :return: air gap reluctance for round-inf structure including fringing effects
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

def r_air_gap_round_inf_sct(air_gap_total_height, core_inner_diameter, core_height, target_reluctance):
    return r_air_gap_round_inf(air_gap_total_height, core_inner_diameter, core_height) - target_reluctance



def r_basic_tablet_cyl(tablet_height, air_gap_basic_height, tablet_radius):
    """
    Do not use this function directly!
    Use it indirectly by using
     - r_air_gap_tablet_cyl
    instead!

    This function calculates the r_basic for a round to infinite structure according to the following paper:
    [according to "A Novel Approach for 3D Air Gap Reluctance Calculations" - J. Mühlethaler, J.W. Kolar, A. Ecklebe]

    Note: this is the same function as r_basic_round_inf, but with clear variable names for tablet-cylinder structure

    :param tablet_height: tablet height = air gap width for tablet-cylinder structure
    :param air_gap_basic_height: air gap height for the BASIC-AIR-GAP (e.g. if you use a round-round structure, this is half of the total air gap).
    :param tablet_radius: tablet radius
    :return: basic reluctance for tablet - cylinder structure
    """
    conductance_basic = mu_0 * (tablet_height / 2 / air_gap_basic_height + 2 / np.pi * (1 + np.log(np.pi * tablet_radius / 4 / air_gap_basic_height)))

    return 1 / conductance_basic

def sigma_tablet_cyl(r_equivalent, tablet_height, air_gap_total_height):
    """
    Do not use this function directly!
    Use it indirectly by using
     - r_air_gap_tablet_cyl
    instead!

    Note: this is the same function as sigma_round, but with clear variable names for tablet-cylinder structure

    :param r_equivalent: this is a series/parallel connection of r_basic, depending on the air gap structure
    :param tablet_height: tablet height
    :param air_gap_total_height: air gap total height (for the total air gap)
    :return: fringing factor 'sigma' for tablet - cylinder structure
    """
    return r_equivalent * mu_0 * tablet_height / air_gap_total_height


def r_air_gap_tablet_cyl(tablet_height, air_gap_total_height, core_inner_diameter, window_w):
    """
    Returns the reluctance of a cylinder-tablet air gap structure and includes fringing effects
    This function calculates the air gap reluctance for a 2D-axisymmetric core.

    :param tablet_height: tablet height in m
    :param air_gap_total_height: total air gap height in m
    :param core_inner_diameter: core inner diameter in m
    :param window_w: core window width in m
    :return: air gap reluctance for tablet - cylinder structure including air gap fringing
    """

    r_outer = core_inner_diameter / 2 + window_w

    # translate practical core dimensions to non-practical air-gap dimensions
    tablet_radius = r_outer - air_gap_total_height

    air_gap_basic_height = air_gap_total_height
    r_basic = r_basic_tablet_cyl(tablet_height, air_gap_basic_height, tablet_radius)

    r_equivalent = r_basic / 2
    sigma = sigma_tablet_cyl(r_equivalent, tablet_height, air_gap_total_height)
    if np.any(sigma > 1):
        raise Exception("Failure in calculating reluctance. Sigma was calculated to >1. Check input parameters!")

    r_air_gap_ideal = np.log(r_outer / (r_outer - air_gap_total_height)) / 2 / mu_0 / np.pi / tablet_height

    r_air_gap = sigma * r_air_gap_ideal

    return r_air_gap

def r_air_gap_tablet_cylinder_sct(air_gap_total_height, core_inner_diameter, tablet_height, window_w, target_reluctance):
    return r_air_gap_tablet_cyl(tablet_height, air_gap_total_height, core_inner_diameter, window_w) - target_reluctance

def r_air_gap_tablet_cyl_no_2d_axi(tablet_height, air_gap_total_length, core_inner_diameter, window_w):
    """
    Returns the reluctance of a cylinder-tablet air gap structure and includes fringing effects
    Note:
    This function differs from r_air_gap_tablet_cyl (ideal 2D axisymmetric core). Here, the air gap reluctance for
    a non-2D-axisymmetric core is taken into account, as a real PQ core is open at the side. So, there is no air gap
    taken into account for the side-sections. The new core_dimension_y parameter describes the width of the
    core when you are in a xy-coordinate system.

    :param tablet_height: tablet height in m
    :param air_gap_total_length: air gap total length
    :param core_inner_diameter: core inner diameter in m
    :param window_w: core window width in m
    :return: air gap reluctance for tablet - cylinder structure including air gap fringing
    """

    r_outer = core_inner_diameter / 2 + window_w

    if np.any(air_gap_total_length >= window_w):
        raise Exception("air_gap_total_height is greater than window_w")

    air_gap_basic_height = air_gap_total_length
    r_basic = r_basic_tablet_cyl(tablet_height, air_gap_basic_height, (core_inner_diameter + 2 * window_w - 2 * air_gap_total_length) / 2)

    r_equivalent = r_basic / 2
    sigma = sigma_tablet_cyl(r_equivalent, tablet_height, air_gap_total_length)
    if np.any(sigma > 1):
        raise Exception("Failure in calculating reluctance. Sigma was calculated to >1. Check input parameters!")

    # Note:
    # the circumference of the air gap differs for open cores (e.g. PQ40/40) to closed ones (ideal rotationally symmetric)
    # The circumference is no more (diameter * pi), but (2 * pi - 4 * alpha) * (core_inner_diameter/2 + window_w), with alpha = arccos(core_dimension_x / (core_inner_diameter + 2 * window_w))
    # See: Dissertation Lukas Keuck
    # For equal pq core sizes (e.g. PQ 40/40), it has been found out that
    #     core_dimension_x / core_dimension_y = 1.45, the error over all available shapes is maximum 7% (compared to datasheet value)
    # Now, the new and partly circumference of the stray-path air gap can be calculated
    # First, the core dimension_y needs to be calculated.
    # Formular 1: core_dimension_x / core_dimension_y = 1.45
    # Formular 2: core_dimension_x * core_dimension_y - (core_inner_diameter / 2 + window_w) ** 2 * np.pi = (core_inner_diameter / 2 ) ** 2 * np.pi
    # Formular 2 assumes that the outer core cross-section of the core is equal to the inner core cross-section
    # Formular 1 & 2 needs to be solved to get core_dimension_y:

    core_dimension_y = np.sqrt( ( core_inner_diameter ** 2 / 4 + (core_inner_diameter / 2 + window_w) ** 2  ) * np.pi / 1.45)
    r_air_gap_ideal_partly = np.log(r_outer / (r_outer - air_gap_total_length)) / mu_0 / (2 * np.pi - 4 * np.arccos(core_dimension_y / 2 / r_outer)) / tablet_height

    r_air_gap = sigma * r_air_gap_ideal_partly

    return r_air_gap

def r_core_tablet(tablet_height, tablet_radius, mu_r_abs, core_inner_diameter):
    """
    Calculates the magentic resistance of the core tablet

    :param tablet_height: tablet height
    :param tablet_radius: tablet radius
    :param mu_r_abs: relative permeability (mu_r) of the core material from datasheet
    :param core_inner_diameter: core inner diameter. For idealized core material, this value can be 0.001.
    """

    return np.log(tablet_radius / (core_inner_diameter / 2)) / (2 * np.pi * mu_0 * mu_r_abs * tablet_height)

def r_core_top_bot_radiant(core_inner_diameter, window_w, mu_r_abs, core_top_bot_height):
    """
    Calculates the top or bottom core material part

    :param core_inner_diameter: core inner diameter
    :param window_w: width of winding window
    :param mu_r_abs: relative permeability (mu_r) of the core material from datasheet
    :param core_top_bot_height: height of the core material top / bottom of the winding window
    """

    return np.log( (core_inner_diameter + 2 * window_w) / core_inner_diameter) / (2 * np.pi * mu_0 * mu_r_abs * core_top_bot_height)

def r_core_round(core_inner_diameter, core_round_height, mu_r_abs):
    """
    Calculates the core reluctance for a round structure

    :param core_round_height: height of the round core part section
    :param core_inner_diameter: core inner diameter
    :param mu_r_abs: relative permeability (mu_r) of the core material from datasheet
    """
    return core_round_height / (mu_0 * mu_r_abs * (core_inner_diameter / 2) ** 2 * np.pi)


def resistance_solid_wire(core_inner_diameter: float, window_w: float, turns_count: int, conductor_radius: float,
                          material: str='Copper') -> float:
    """
    Calculates the resistance of a solid wire.

    :param core_inner_diameter: core inner diameter
    :type core_inner_diameter: float
    :param window_w: winding window width
    :type window_w: float
    :param turns_count: number of turns
    :type turns_count: int
    :param conductor_radius: conductor radius
    :type conductor_radius: float
    :param material: Material, e.g. "Copper" or "Aluminium"
    :type material: str
    :return: total resistance of wire
    :rtype: float
    """

    # simplification: middle turn length
    # figure out middle length of one turn for given geometry
    turn_radius = core_inner_diameter / 2 + conductor_radius
    turn_count = 1
    total_turn_length = 0
    while turn_radius <= (core_inner_diameter / 2 + window_w - conductor_radius):
        total_turn_length += turn_radius * 2 * np.pi
        turn_radius += 2 * conductor_radius
        turn_count += 1
    middle_turn_length = total_turn_length / turn_count

    total_turn_length = turns_count * middle_turn_length

    sigma_copper = ff.wire_material_database()[material]["sigma"]

    # return R = rho * l / A
    return total_turn_length / (conductor_radius ** 2 * np.pi) / sigma_copper

def i_rms(time_current_matrix: np.array) -> float:
    """
    RMS calculation from a time-current-vector

    :param time_current_matrix: time and current in format [[0, 0.5e-6, 2.5e-6, 3e-6, 5e-6], [16.55, -10.55, -16.55, 10.55, 16.55]]
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
        square_integral = gradient ** 2 / 3 * delta_time ** 3 + gradient * delta_time ** 2 * y_axis + y_axis ** 2 * delta_time

        # add (partly) square integral to total integration sum
        square_integral_sum += square_integral

    # return "mean" and "root" to finalize rms calculation
    return np.sqrt(square_integral_sum / time[-1])

