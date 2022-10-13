import numpy as np
import itertools
from Functions import get_dicts_with_keys_and_values

# ----------------------------------------------------------------------------------------------------------------------
# Goal Parameters
wp_f_switch = [200000]
wp_lambda = 17
n = 3.2
ratio_lm_ls = 10
power_nom = 2000
power_max = 2500
voltage_nom = 235  # secondary
voltage_min = 175  # secondary
wp_power_voltage_out = [[power_max, voltage_min], [power_nom, voltage_nom]]

# Load Nikolas' Data
path = 'C:/Users/tillp/sciebo/Exchange_DAB/current_shapes.npy'
data = np.load(path, allow_pickle=True)


working_points = []  # len = frequencies;
L_goal = {}
for f_switch in wp_f_switch:
    # Goal Inductances
    L_s1 = wp_lambda/f_switch
    L_h = ratio_lm_ls * L_s1

    # Inductance Matrix
    L_11 = L_s1 + L_h
    L_22 = L_h / n ** 2
    M = L_h / n
    L = np.array([[L_11, M], [M, L_22]])

    L_goal[f"{f_switch}"] = L

    wp_data = []
    for power, voltage in wp_power_voltage_out:
        wp_data.append(get_dicts_with_keys_and_values(data, power=power, voltage=voltage, frequency=f_switch,
                                                      wp_lambda=wp_lambda, wp_n=n, ratio_lm_ls=ratio_lm_ls))


    working_points.append(np.array(wp_data).flatten())

# ----------------------------------------------------------------------------------------------------------------------
# TODO: Dicts schon hier packen
# Reluctance Model Parameters
window_h = [0.029]
window_w = [0.01105]
core_w = [0.0145]
real_core_width = [28e-3]
midpoint = [29, 30, 31]  # [25, 30, 35, 40, 45]  # stray_path.midpoint
b_stray_rel_overshoot = [1.4, 1.5, 1.6]

width = []  # [0.004]  # stray_path.width


N1 = np.arange(28, 31)  # Turns in main window
N2 = np.arange(6, 9)  #  Turns in main window
Ns1 = np.arange(0, 1)  # Turns in stray window
Ns2 = np.arange(5, 8)  #  Turns in stray window

"""
N1 = np.arange(26, 27)  # [27, 20]  # Turns in main window
N2 = np.arange(6, 8)  # [7]  # Turns in main window
Ns1 = np.arange(0, 2)  # [5]  # Turns in stray window
Ns2 = np.arange(4, 5)  # [6]  # Turns in stray window
"""

N_flat = list(itertools.product(N1, N2, Ns1, Ns2))
N = [np.reshape(N_flat_single, (2, 2)) for N_flat_single in N_flat]

# Create List of Dictionaries for Reluctance Model
reluctance_parameter_categories = ["window_w", "window_h", "core_w", "real_core_width", "midpoint", "N", "b_stray_rel_overshoot"]
reluctance_parameter_values = list(itertools.product(window_w, window_h, core_w, real_core_width, midpoint, N, b_stray_rel_overshoot))

reluctance_parameters = []
for objects in reluctance_parameter_values:
    reluctance_parameters.append({key: value for key, value in zip(reluctance_parameter_categories, objects)})

# ----------------------------------------------------------------------------------------------------------------------
# Litz Definitions
litzes = {
    "prim": {"strand_radius": 0.0355e-3, "N_strands": 405, "conductor_radius": 0.90e-3, "isolation": 50e-6},
    "sec": {"strand_radius": 0.0355e-3, "N_strands": 600, "conductor_radius": 1.1e-3, "isolation": 100e-6}}


# Create List of Dictionaries for FEM simulations
non_reluctance_categories = ["litzes"]
non_reluctance_values = [litzes]
# print(non_reluctance_values)


# non_reluctance_parameters = []
# for objects in non_reluctance_values:
#     non_reluctance_parameters.append({key: value for key, value in zip(non_reluctance_categories, objects)})
non_reluctance_parameters = [{"litzes": litzes}]

# print(non_reluctance_parameters)
# print(len(non_reluctance_parameters))


# Unused:
ki = 0.53
alpha = 1.50
beta = 2.38
