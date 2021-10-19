from femmt import MagneticComponent
import numpy as np
import itertools

#                                               -- Definitions --
# ----------------------------------------------------------------------------------------------------------------------
mu0 = 4e-7*np.pi
mur = 3000
inductance_calculation = 1
skin_accuracy = 0.5
frequency = 250000
I0 = 6 * 1.25  # 6 is peak of current wave_form

# Goal Parameters
L_s1 = 6.8e-5
L_h = 10 * L_s1
n = 3.2
# Inductance Matrix
L_11 = L_s1 + L_h
L_22 = L_h / n**2
M = L_h / n
L = np.array([[L_11, M], [M, L_22]])

# TODO: Dicts schon hier packen
# Reluctance Model Parameters
window_h = [0.0295]
window_w = [0.012]
core_w = [0.015]
midpoint = [35]
b_stray_rel_overshoot = [2]
width = []  # [0.004]
N1 = [27, 20]  # Turns in main window
N2 = [7]  # Turns in main window
Ns1 = [5]  # Turns in stray window
Ns2 = [6]  # Turns in stray window

# N1 = np.arange(10, 40)  # [27, 20]  # Turns in main window
# N2 = np.arange(10, 40)  # [7]  # Turns in main window
# Ns1 = np.arange(10, 20)  # [5]  # Turns in stray window
# Ns2 = np.arange(10, 20)  # [6]  # Turns in stray window

N_flat = list(itertools.product(N1, N2, Ns1, Ns2))
N = [np.reshape(N_flat_single, (2, 2)) for N_flat_single in N_flat]

# Create List of Dictionaries for Reluctance Model
reluctance_parameter_categories = ["window_w", "window_h", "core_w", "midpoint", "N", "b_stray_rel_overshoot"]
reluctance_parameter_values = list(itertools.product(window_w, window_h, core_w, midpoint, N, b_stray_rel_overshoot))

reluctance_parameters = []
for objects in reluctance_parameter_values:
    reluctance_parameters.append({key: value for key, value in zip(reluctance_parameter_categories, objects)})


#                                           -- Reluctance Model --
# ----------------------------------------------------------------------------------------------------------------------
geo = MagneticComponent(component_type="integrated_transformer")


valid_reluctance_parameters = geo.reluctance_model.air_gap_design(L_goal=L,
                                                                  parameters_init=reluctance_parameters,
                                                                  current=[1.25*I0, -1.25*I0*3.2],
                                                                  b_max=0.3, b_stray=0.25,
                                                                  stray_path_parametrization="max_flux")


#                                         -- FEM Simulation --
# --------------------------------------------------------------------------------------------------------------
# Strand Parameters
strand_radius = [0.025e-3]
N_strands_prim = [300, 450]
N_strands_sec = [300, 450]

# Create List of Dictionaries for FEM simulations
non_reluctance_categories = ["strand_radius", "N_strands_prim", "N_strands_sec"]
non_reluctance_values = list(itertools.product(strand_radius, N_strands_prim, N_strands_sec))
# print(non_reluctance_values)


non_reluctance_parameters = []
for objects in non_reluctance_values:
    non_reluctance_parameters.append({key: value for key, value in zip(non_reluctance_categories, objects)})
# print(non_reluctance_parameters)
# print(len(non_reluctance_parameters))


# Bring together the
# print(f"{valid_reluctance_parameters=}")
# print(len(valid_reluctance_parameters))
FEM_parameters = []
for reluctance_parameters in valid_reluctance_parameters:
    for non_reluctance_parameter in non_reluctance_parameters:
        # print(non_reluctance_parameter)
        # print(reluctance_parameters)
        FEM_parameters.append(dict(reluctance_parameters, **non_reluctance_parameter))

# print(FEM_parameters)
# print(len(FEM_parameters))


do_simulate = True

geo.ki = 0.53
geo.alpha = 1.50
geo.beta = 2.38

for n_par, parameters in enumerate(FEM_parameters):
    # print(f"{parameters=}")

    if do_simulate:

        # -- FFT --
        print(f"{n_par=}")

        # -- Component Design --
        geo.core.update(type="EI",
                        window_h=parameters["window_h"], window_w=parameters["window_w"], core_w=parameters["core_w"],
                        non_linear=False, material=95_100, re_mu_rel=mur)

        geo.air_gaps.update(method="percent",
                            n_air_gaps=2,
                            position_tag=[0, 0],
                            air_gap_h=[parameters["R_bot"], parameters["R_top"]],
                            air_gap_position=[parameters["midpoint"] -
                                              (parameters["stray_path_width"] + parameters["R_bot"]) / 2
                                              / parameters["window_h"] * 100,
                                              parameters["midpoint"] +
                                              (parameters["stray_path_width"] + parameters["R_top"]) / 2
                                              / parameters["window_h"] * 100])

        geo.stray_path.update(start_index=0,
                              radius=geo.core.core_w/2+geo.core.window_w-parameters["R_stray"])

        geo.update_conductors(n_turns=[[parameters["N"][0, 0], parameters["N"][1, 0]],
                                        [parameters["N"][0, 1], parameters["N"][1, 1]]],
                              conductor_type=["litz", "litz"],
                              winding=["interleaved", "interleaved"],
                              scheme=["horizontal", "horizontal"],  # ["square", "square"]  "horizontal"
                              conductor_radii=[0.0005, 0.0009],
                              litz_para_type=['implicit_litz_radius', 'implicit_litz_radius'],
                              strands_numbers=[parameters["N_strands_prim"], parameters["N_strands_sec"]],
                              ff=[0.5, 0.5],
                              strand_radii=[parameters["strand_radius"], parameters["strand_radius"]],
                              cond_cond_isolation=[0.0002, 0.0002, 0.0004],
                              core_cond_isolation=[0.002])

        # -- Simulation --
        if inductance_calculation:
            geo.get_inductances(I0=I0, op_frequency=frequency, skin_mesh_factor=skin_accuracy)
            # print(f"\n-------------------------------\n"
            #       f"Reluctance Model Values:\n"
            #       f"L_11 = {L_11}\n"
            #       f"L_22 = {L_22}\n"
            #       f"M = {M}\n"
            #
            #       f"R = {R}\n"
            #       f"R_stray_a = {R_stray}\n"
            #       f"R_top_a = {R_top}\n"
            #       f"R_bot_a = {R_bot}\n"
            #       # f"R_stray_cheap = {R_stray_cheap}\n"
            #       f"\n"
            #       # f"air_gap_h_stray = {air_gap_h_stray}\n"
            #       # f"air_gap_h_top = {air_gap_h_top}\n"
            #       # f"air_gap_h_bot = {air_gap_h_bot}\n"
            #       )
            #
            # L_FEM_11 = geo.L_h_conc + geo.L_s_conc
            # L_FEM_22 = geo.L_h_conc / geo.n_conc ** 2
            # M_FEM = geo.L_h_conc / geo.n_conc
            #
            # L_FEM = np.array([[L_FEM_11, -M_FEM], [-M_FEM, L_FEM_22]])
            # """
            # print(f"\n"
            #       f"L_FEM = {L_FEM}\n"
            #       f"L_FEM-1 = {np.linalg.inv(L_FEM)}\n"
            #       f"det(L_FEM-1) = {np.linalg.det(np.linalg.inv(L_FEM))}\n"
            #       f"L_a = {L_a}\n"
            #       f"L_a-1 = {np.linalg.inv(L_a)}\n"
            #       f"det(L_a-1) = {np.linalg.det(np.linalg.inv(L_a))}\n"
            #       f"L_a-1/4 = {np.linalg.inv(L_a)/4}")
            # """
            #
            # R_FEM = np.matmul(np.matmul(N, np.linalg.inv(L_FEM)), np.transpose(N))
            # print(f"TEST:{np.matmul(np.linalg.inv(L_FEM), L_FEM)}")
            # # R_something = np.matmul(np.matmul(N, np.linalg.inv(L_a)), np.transpose(N))
            # # print(L_FEM)
            # R_FEM_stray = -R_FEM[0, 1]
            # R_FEM_top = R_FEM[0, 0] - R_FEM_stray
            # R_FEM_bot = R_FEM[1, 1] - R_FEM_stray
            # print(f"\n-------------------------------\n"
            #       f"FEM Values:\n"
            #       f"L_FEM_11 = {L_FEM_11}\n"
            #       f"L_FEM_22 = {L_FEM_22}\n"
            #       f"M_FEM = {M_FEM}\n"
            #       # f"R_something = {R_something}\n"
            #       f"R_FEM = {R_FEM}\n"
            #       f"R_FEM_stray = {R_FEM_stray}\n"
            #       f"R_FEM_top = {R_FEM_top}\n"
            #       f"R_FEM_bot = {R_FEM_bot}\n")

        else:
            # geo.single_simulation(freq=frequency, current=[I0, N1/N2*I0], phi_deg=[0, 180],
            # skin_mesh_factor=skin_accuracy)
            geo.single_simulation(freq=frequency,
                                  current=[I0, 3.2*I0],
                                  phi_deg=[0, 180],
                                  skin_mesh_factor=skin_accuracy)

"""
# -- Sweep Parameters --
# - Windings
N1 = 27 # Turns in main window
N2 = 7  # Turns in main window
Ns1 = 5  # Turns in stray window
Ns2 = 6  # Turns in stray window
# - Air gaps
air_gap_h_bot = 0.0012
air_gap_h_top = 0.0002
air_gap_h_stray = 0.0002
# - Core
stray_path_pos = 30
# - Conductors
strand_radius = 0.025e-3
N_strands_prim = 400
N_strands_sec = 800
"""