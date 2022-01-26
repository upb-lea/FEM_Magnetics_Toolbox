from femmt_functions import store_as_npy_in_directory
import numpy as np
from matplotlib import pyplot as plt
from operator import itemgetter
from femmt import MagneticComponent
from femmt_functions import get_dict_with_unique_keys, get_dicts_with_keys_and_values, find_common_frequencies, \
    sort_out_small_harmonics, store_as_npy_in_directory
from DAB_Input_Data import working_points, non_reluctance_parameters, L_goal, power_nom, power_max



pathC = f"C:/Users/tillp/OneDrive/Documents/GitHub/FEM_Magnetics_Toolbox/femmt/MA/final/local_4/Optimal_Result.npy"
pathD = f"C:/Users/tillp/OneDrive/Documents/GitHub/FEM_Magnetics_Toolbox/femmt/MA/final/local_3/Result_FEM_parameters_200000.npy"

print(np.load(pathC, allow_pickle=True))

result = np.load(pathD, allow_pickle=True)
parameters = dict(result[0])
print(parameters)
print(type(parameters))
# parameters["N"] = np.array([[2, 2], [1, 0]])
# parameters["R_top"] = 1e-6
parameters["R_bot"] = np.round(parameters["R_bot"] * 0.95, 4)
parameters["R_stray"] = np.round(parameters["R_stray"] * 1.0, 4)
parameters["stray_path_width"] = 0.0045 + 0.0045  # np.round(parameters["stray_path_width"] * 1.15, 5)
parameters["midpoint"] = 31
parameters["real_core_width"] = 0.024
parameters["window_h"] = 0.029 + 0.0045

print(f"{parameters['R_stray']=}")
print(f"{parameters['R_bot']=}")
print(f"{parameters['stray_path_width']=}")

f_opt = parameters["frequency"]




wp_data = None
for working_point in working_points:
    dicts_nom = get_dicts_with_keys_and_values(working_point, power=power_nom)
    dict_nom_in = get_dict_with_unique_keys(dicts_nom, 'wp_ib_il_phase_rad_vec')
    wp_frequency = dict_nom_in["frequency"]

    if wp_frequency == f_opt:
        wp_data = working_point


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Nominal Power Data
dicts_nom = get_dicts_with_keys_and_values(wp_data, power=power_nom)
dict_nom_in = get_dict_with_unique_keys(dicts_nom, 'wp_ib_il_phase_rad_vec')
dict_nom_out = get_dict_with_unique_keys(dicts_nom, 'wp_ob_il_phase_rad_vec')

# Working point switching frequency
frequency = dict_nom_in["frequency"]

# Working point Input current data
frequencies_in = dict_nom_in['wp_ib_il_frequency_vec']
currents_in = dict_nom_in['wp_ib_il_amplitude_vec']
phases_in = dict_nom_in['wp_ib_il_phase_rad_vec']
time_nom_in = dict_nom_in['wp_ib_il_vec']

# Working point Output current data
frequencies_out = dict_nom_out['wp_ob_il_frequency_vec']
currents_out = dict_nom_out['wp_ob_il_amplitude_vec']
phases_out = dict_nom_out['wp_ob_il_phase_rad_vec']
time_nom_out = dict_nom_out['wp_ob_il_vec']

# Find common frequencies
current_pairs_nom, phase_pairs_nom, frequencies = find_common_frequencies(currents_in, phases_in, frequencies_in,
                                                                          currents_out, phases_out, frequencies_out)

# Throw out all harmonics with low impact (f.e. amplitude smaller 1%)
phase_pairs_nom, current_pairs_nom, frequencies = sort_out_small_harmonics(phase_pairs_nom, current_pairs_nom,
                                                                           frequencies, limes=0.05)

# Sort with ascending frequencies
frequencies, phase_pairs_nom, current_pairs_nom = list(zip(*sorted(zip(frequencies,
                                                                       phase_pairs_nom,
                                                                       current_pairs_nom))))

# Change phases from radiant to degree and add 180° phase shift (different conventions)
phase_pairs_nom = [phase_pair*180/np.pi for phase_pair in phase_pairs_nom]
phase_pairs_nom = [[phase_pair[0]+180, phase_pair[1]] for phase_pair in phase_pairs_nom]

# Transform back so simple python lists
phase_pairs_nom = [list(phase_pair) for phase_pair in phase_pairs_nom]
current_pairs_nom = list(current_pairs_nom)
frequencies = list(frequencies)



geo = MagneticComponent(component_type="integrated_transformer")
geo.s = 0.5
geo.core.update(type="EI",
                window_h=parameters["window_h"], window_w=parameters["window_w"], core_w=parameters["core_w"],
                non_linear=False, material=95_100, re_mu_rel=3000)

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
                      radius=geo.core.core_w / 2 + geo.core.window_w - parameters["R_stray"])

geo.update_conductors(n_turns=[[parameters["N"][0, 0], parameters["N"][1, 0]],
                               [parameters["N"][0, 1], parameters["N"][1, 1]]],
                      conductor_type=["litz", "litz"],
                      winding=["interleaved", "interleaved"],
                      scheme=["horizontal", "horizontal"],  # ["square", "square"]  "horizontal"
                      conductor_radii=[parameters["litzes"]["prim"]["conductor_radius"],
                                       parameters["litzes"]["sec"]["conductor_radius"]],
                      litz_para_type=['implicit_ff', 'implicit_ff'],
                      strands_numbers=[parameters["litzes"]["prim"]["N_strands"],
                                       parameters["litzes"]["sec"]["N_strands"]],
                      strand_radii=[parameters["litzes"]["prim"]["strand_radius"],
                                    parameters["litzes"]["sec"]["strand_radius"]],
                      cond_cond_isolation=[2 * parameters["litzes"]["prim"]["isolation"],
                                           2 * parameters["litzes"]["sec"]["isolation"],
                                           parameters["litzes"]["prim"]["isolation"] +
                                           parameters["litzes"]["prim"]["isolation"]],
                      core_cond_isolation=[1e-3, 200e-6])  # 1st argument close to air gaps

# -- Simulation --
I0 = 6  # 6 * 1.25  # 6 is peak of current wave_form
skin_accuracy = 0.5
geo.get_inductances(I0=I0, op_frequency=frequency, skin_mesh_factor=skin_accuracy)

# geo.single_simulation(freq=frequencies[0],
#                       current=current_pairs_nom[0],
#                       phi_deg=phase_pairs_nom[0])

FEM_results = geo.excitation_sweep(frequencies=frequencies,
                                   currents=current_pairs_nom,
                                   phi=phase_pairs_nom,
                                   show_last=True,
                                   return_results=True)
