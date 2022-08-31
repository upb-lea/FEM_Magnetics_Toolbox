import femmt as fmt
import numpy as np
import os
from femmt.femmt_reluctance import MagneticCircuit

# Component type
component = "inductor"

# Working directory can be set arbitrarily
working_directory = os.path.join(os.path.dirname(__file__), '..')

# MagneticComponent class object
geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=working_directory)

# Formation of core-geometry list from core database
core_db1 = fmt.core_database()["PQ 20/16"]
core_db2 = fmt.core_database()["PQ 20/20"]
core_db3 = fmt.core_database()["PQ 26/20"]
core_db4 = fmt.core_database()["PQ 26/25"]

core_width_list = [core_db1["core_w"], core_db2["core_w"], core_db3["core_w"], core_db4["core_w"]]
core_window_w_list = [core_db1["window_w"], core_db2["window_w"], core_db3["window_w"], core_db4["window_w"]]
core_window_h_list = [core_db1["window_h"], core_db2["window_h"], core_db3["window_h"], core_db4["window_h"]]

# Air-gap and core parameters list creation
no_of_turns = list(np.linspace(5, 19, 15))
n_air_gaps = list(np.linspace(1, 5, 5))
air_gap_length = list(np.linspace(0.000001, 0.0005, 200))  # 200
air_gap_position = list(np.linspace(0, 100, 101))  # 101
mu_rel = [2700, 3000, 3100, 3200]
mult_air_gap_type = [1, 2]  # 'Type1 = with corner air-gaps; 'Type2' = without air-gaps; 'Type0' = single air-gap

goal_inductance = 200 * 1e-6
conductor_radius = 0.0013
winding_factor = 0.91

# Reluctance model simulation
mc1 = MagneticCircuit(core_width_list, core_window_h_list, core_window_w_list, no_of_turns, n_air_gaps, air_gap_length, air_gap_position, mu_rel, mult_air_gap_type)
mc1.core_reluctance()
mc1.air_gap_reluctance()
mc1.data_matrix = np.hstack((mc1.data_matrix, np.reshape(mc1.core_h_middle, (mc1.data_matrix_len, 1))))  # position: 10
mc1.data_matrix = np.hstack((mc1.data_matrix, np.reshape(mc1.r_inner, (mc1.data_matrix_len, 1))))  # position: 11
mc1.data_matrix = np.hstack((mc1.data_matrix, np.reshape(mc1.r_outer, (mc1.data_matrix_len, 1))))  # position: 12
mc1.data_matrix = np.hstack((mc1.data_matrix, np.reshape(mc1.area[:, 0], (mc1.data_matrix_len, 1))))  # position: 13
mc1.data_matrix = np.hstack((mc1.data_matrix, np.reshape(mc1.area[:, 4], (mc1.data_matrix_len, 1))))  # position: 14

data_matrix_1 = mc1.data_matrix[np.where((mc1.data_matrix[:, 4] * np.pi * conductor_radius * conductor_radius) < (winding_factor * mc1.data_matrix[:, 1] * mc1.data_matrix[:, 2]))]

# First filter: Based on +-10% inductance tolerance band
FEM_data_matrix = data_matrix_1[np.where((data_matrix_1[:, 9] > 0.9 * goal_inductance) *
                                           (data_matrix_1[:, 9] < 1.1 * goal_inductance))]

max_current = 3 # With an assumption of sinusoidal current waveform
freq = 100 * 1e3 # magnetic flux frequency
mu_imag = 100
Cu_sigma = 5.96 * 1e7


volume_center = (np.pi * FEM_data_matrix[:, 0] ** 2) * (FEM_data_matrix[:, 1] + FEM_data_matrix[:, 10] - (FEM_data_matrix[:, 5] * FEM_data_matrix[:, 6]))
volume_outer = (np.pi * ((FEM_data_matrix[:, 12] ** 2) - (FEM_data_matrix[:, 11] ** 2))) * (FEM_data_matrix[:, 1] + FEM_data_matrix[:, 10])

max_flux = (FEM_data_matrix[:, 9] * max_current) / FEM_data_matrix[:, 4]    # max_flux = L * i_max / N

B_max_center = max_flux / FEM_data_matrix[:, 13]
P_hyst_center = - 0.5 * (2 * np.pi * freq) * mu_imag * fmt.mu0 * (B_max_center / (fmt.mu0 * FEM_data_matrix[:, 3])) ** 2
P_hyst_density_center = P_hyst_center * volume_center

B_max_outer = max_flux / FEM_data_matrix[:, 14]
P_hyst_outer = - 0.5 * (2 * np.pi * freq) * mu_imag * fmt.mu0 * (B_max_outer / (fmt.mu0 * FEM_data_matrix[:, 3])) ** 2
P_hyst_density_outer = P_hyst_outer * volume_outer

P_hyst_density_middle = - 0.5 * (2 * np.pi * freq) * mu_imag * fmt.mu0 * ((max_flux / (fmt.mu0 * FEM_data_matrix[:, 3])) ** 2) * (1/(2 * np.pi * FEM_data_matrix[:, 10])) * np.log((FEM_data_matrix[:, 11] * 2) / FEM_data_matrix[:, 0])

P_hyst_density_total = P_hyst_density_center + (2 * P_hyst_density_middle) + P_hyst_density_outer
FEM_data_matrix = np.hstack((FEM_data_matrix, np.reshape(P_hyst_density_total, (len(P_hyst_density_total), 1)))) # position: 15

# Winding loss (only DC loss)
Resistance = (FEM_data_matrix[:, 4] * 2 * np.pi * (FEM_data_matrix[:, 0] / 2 + conductor_radius)) / Cu_sigma * (np.pi * conductor_radius * conductor_radius)
DC_loss = ((max_current ** 2) / 2) * Resistance
FEM_data_matrix = np.hstack((FEM_data_matrix, np.reshape(DC_loss, (len(DC_loss), 1)))) # position: 16

sorted_FEM_data_matrix = FEM_data_matrix[FEM_data_matrix[:, 15].argsort()]

qwerty = 1
# def aut_action_run_simulation(sim_value):
#     geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor)
#     core = fmt.Core(core_w=sim_value[0], window_h=sim_value[1], window_w=sim_value[2],
#                     mu_rel=3100, phi_mu_deg=12,
#                     sigma=0.6)
#     geo.set_core(core)
#
#     # 3. set air gap parameters
#     air_gaps = fmt.AirGaps(fmt.AirGapMethod.Manually, core)
#     air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.001, sim_value[5])
#     geo.set_air_gaps(air_gaps)
#
#     # 4. set conductor parameters: use solid wires
#     winding = fmt.Winding(sim_value[3], 0, fmt.Conductivity.Copper,
#                           fmt.WindingType.Primary,
#                           fmt.WindingScheme.Square)
#     winding.set_litz_conductor(None, 600, 35.5e-6, 0.6)
#     # winding.set_solid_conductor(0.0015)
#     geo.set_windings([winding])
#
#     # 5. set isolations
#     isolation = fmt.Isolation()
#     isolation.add_core_isolations(0.001, 0.001, 0.002, 0.001)
#     isolation.add_winding_isolations(0.0001)
#     geo.set_isolation(isolation)
#
#     # 5. create the model
#     geo.create_model(freq=100000, visualize_before=False, save_png=False)
#
#     # 6.a. start simulation
#     geo.single_simulation(freq=100000, current=[3], show_results=False)
#
#
# # FEM_data_matrix = [core_w, window_h, window_w, no_of_turns, n_air_gaps, air_gap_h, inductance]
# # FEM_data_matrix = [[0.0149, 0.0295, 0.01105, 8, 1, 0.0001, 0.0002],
# #                    [0.0149, 0.0295, 0.01105, 8, 1, 0.0002, 0.0002],
# #                    [0.0149, 0.0295, 0.01105, 8, 1, 0.0003, 0.0002],
# #                    [0.0149, 0.0295, 0.01105, 8, 1, 0.0004, 0.0002]]
# num_r = len(FEM_data_matrix)
# m = []
# for i in range(num_r):
#     m.append([float(x) for x in FEM_data_matrix[i, :]])
#     for j in m:
#         aut_action_run_simulation(j)
