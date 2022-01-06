from matplotlib import pyplot as plt
from operator import itemgetter
import numpy as np
from femmt import MagneticComponent
from femmt_functions import get_dict_with_unique_keys, get_dicts_with_keys_and_values, find_common_frequencies, \
    sort_out_small_harmonics, store_as_npy_in_directory
from DAB_Input_Data import working_points, non_reluctance_parameters, L_goal, power_nom, power_max
import itertools

# Find optimal frequency
f_opt = None
min_losses = None
results_init_sorted = None
for wp_data in working_points:

    dicts_nom = get_dicts_with_keys_and_values(wp_data, power=power_nom)
    dict_nom_in = get_dict_with_unique_keys(dicts_nom, 'wp_ib_il_phase_rad_vec')
    # Working point switching frequency
    frequency = dict_nom_in["frequency"]

    pathD = f"C:/Users/tillp/OneDrive/Documents/GitHub/FEM_Magnetics_Toolbox/femmt/MA/final/Result_FEM_parameters_{frequency}.npy"

    results = np.load(pathD, allow_pickle=True)

    # print(f"{len(results)=}")
    results = [item for item in results if ("FEM_results" not in item)]

    # Add total losses to each dictionary
    results = [dict(item, **{'total_losses': item["p_hyst_nom"] +
                                             np.sqrt(np.mean(np.array(item["j2H"])**2)) +
                                             np.sqrt(np.mean(np.array(item["j2F"])**2))}) for item in results]
    # Sort with ascending total losses
    results_sorted = sorted(results, key=itemgetter('total_losses'))

    if min_losses is None or min_losses > results_sorted[0]["total_losses"]:
        min_losses = results_sorted[0]["total_losses"]
        f_opt = frequency
        results_init_sorted = results_sorted
    # store_as_npy_in_directory(result_directory, f"results_{frequency}", results)

    # print(f"{len(results)=}")
    # print(f"{results=}")

    # Loss Plots
    total_losses_sorted = [item['total_losses'] for item in results_sorted]
    plt.scatter(np.arange(1, len(total_losses_sorted)+1), total_losses_sorted, label=f"{int(frequency/1000)} kHz")


# Print the optimal loss and corresponding frequency
print(f"{min_losses=} \n"
      f"{f_opt=}")

# show plot
plt.ylabel("Total Losses in W")
plt.xlabel("Sorted Parameter Runs")
plt.legend()
plt.grid()
plt.show()


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

# Maximal Power Data (worst case)
dicts_max = get_dicts_with_keys_and_values(wp_data, power=power_max)
dict_max_in = get_dict_with_unique_keys(dicts_max, 'wp_ib_il_phase_rad_vec')
dict_max_out = get_dict_with_unique_keys(dicts_max, 'wp_ob_il_phase_rad_vec')

time_max_in = dict_max_in['wp_ib_il_vec']
time_max_out = dict_max_out['wp_ob_il_vec']


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create List of Dictionaries for Reluctance Model
midpoint = [-1, 0, 1]  # [25, 30, 35, 40, 45]  # stray_path.midpoint
b_stray_rel_overshoot = [-0.1, 0, 0.1]
N1 = [-1, 0, 1]
N2 = [-1, 0, 1]
Ns1 = [-1, 0, 1]
Ns2 = [-1, 0, 1]
N_flat = list(itertools.product(N1, N2, Ns1, Ns2))
N = [np.reshape(N_flat_single, (2, 2)) for N_flat_single in N_flat]

update_incrementation_categories = ["midpoint", "N", "b_stray_rel_overshoot"]
update_incrementation_values = list(itertools.product(midpoint, N, b_stray_rel_overshoot))

update_incrementation_matrix = []
for objects in update_incrementation_values:
    update_incrementation_matrix.append({key: value for key, value in zip(update_incrementation_categories, objects)})


print(update_incrementation_matrix)
# input("bitte wei")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def new_iteration(init_opt, update_incrementation_matrix, result_history):

    additional_parameters = []
    # Remove Result Data
    keys_to_remove = ["stray_path_width", "R_top_b_peak", "R_bot_b_peak", "R_stray_b_peak", "R_top", "R_bot",
                      "R_stray", "R_stray_real", "p_hyst_nom_1st", "p_hyst_nom", "j2F", "j2H", "p_hyst",
                      "total_losses"]

    init_opt_only_parameters = init_opt
    for key in keys_to_remove:
        del init_opt_only_parameters[key]

    result_history_only_parameters = [dictionary.copy() for dictionary in result_history]
    for part in result_history_only_parameters:
        for key in keys_to_remove:
            del part[key]

    for update_incrementation in update_incrementation_matrix:
        possible_candidate = dict(init_opt_only_parameters)
        valid_candidate = True
        for key, value in update_incrementation.items():
            if key == 'N':
                possible_candidate[key] = np.array(value) + np.array(possible_candidate[key])
                if np.any(possible_candidate[key] < 0):
                    valid_candidate = False

            else:
                possible_candidate[key] += value
                # Check for negative/not allowed values
                if possible_candidate[key] < 0:
                    valid_candidate = False

        if possible_candidate == any(result_history_only_parameters):
            valid_candidate = False

        if valid_candidate:
            additional_parameters.append(possible_candidate)

    return additional_parameters


converged = False
no_iter = 0
while not converged:
    no_iter += 1
    result_directory = f"C:/Users/tillp/OneDrive/Documents/GitHub/FEM_Magnetics_Toolbox/femmt/MA/final/local_{no_iter}"

    # optimum of last iteration
    init_opt = dict(results_init_sorted[0])

    # get new parameter set with local grid
    new_opt_set = new_iteration(init_opt, update_incrementation_matrix, results_init_sorted)

    print(init_opt, "\n\n")
    print(new_opt_set)

    #                                           -- Reluctance Model --
    # -------------------------------------------------------------------------------------------------------------
    store_as_npy_in_directory(result_directory, f"reluctance_parameters_{frequency}", new_opt_set)

    geo = MagneticComponent(component_type="integrated_transformer")

    # valid_reluctance_parameters: List of dictionaries (Core Parameters, Stray Path Parameters, Winding Parameters,
    #                                                    Fluxes, Analytic Core Losses)
    valid_reluctance_parameters = geo.reluctance_model.air_gap_design(L_goal=L_goal[f"{frequency}"],
                                                                      parameters_init=new_opt_set,
                                                                      max_current=[time_max_in, time_max_out],
                                                                      nom_current=[time_nom_in, time_nom_out],
                                                                      nom_current_1st=[current_pairs_nom[0][0],
                                                                                       current_pairs_nom[0][1]],
                                                                      nom_phase_1st=[phase_pairs_nom[0][0],
                                                                                     phase_pairs_nom[0][1]],
                                                                      f_1st=frequency,
                                                                      b_max=0.3, b_stray=0.25,
                                                                      stray_path_parametrization="max_flux",
                                                                      visualize_waveforms=True)

    print(f"{len(valid_reluctance_parameters)=}")

    store_as_npy_in_directory(result_directory, f"valid_reluctance_parameters_{frequency}", valid_reluctance_parameters)

    # Sort out unacceptable hysteresis losses:
    epsilon = 1.2
    hyst = [res["p_hyst_nom"] for res in valid_reluctance_parameters]
    hyst_minimum = min(hyst)
    reluctance_parameters_hyst_good = [item for item in valid_reluctance_parameters if
                                       item["p_hyst_nom"] < hyst_minimum * epsilon]

    print(f"{len(reluctance_parameters_hyst_good)=}")
    store_as_npy_in_directory(result_directory, f"reluctance_parameters_hyst_good{frequency}",
                              reluctance_parameters_hyst_good)

    # input("Press Enter to continue...")

    #                                         -- FEM Simulation --
    # --------------------------------------------------------------------------------------------------------------
    # Bring together valid reluctance parameters and non reluctance parameters
    FEM_parameters = []
    for valid_parameters in reluctance_parameters_hyst_good:
        for non_reluctance_parameter in non_reluctance_parameters:
            FEM_parameters.append(dict(valid_parameters, **non_reluctance_parameter))

    store_as_npy_in_directory(result_directory, f"Init_FEM_parameters_{frequency}", FEM_parameters)

    # Remove:
    geo.s = 0.5

    for n_par, parameters in enumerate(FEM_parameters):
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
        # geo.get_inductances(I0=I0, op_frequency=frequency, skin_mesh_factor=skin_accuracy)

        # geo.single_simulation(freq=frequencies[0],
        #                       current=current_pairs_nom[0],
        #                       phi_deg=phase_pairs_nom[0])

        FEM_results = geo.excitation_sweep(frequencies=frequencies,
                                           currents=current_pairs_nom,
                                           phi=phase_pairs_nom,
                                           show_last=False,
                                           return_results=True)

        FEM_parameters[n_par] = dict(FEM_parameters[n_par], **FEM_results)



    # store_as_npy_in_directory(result_directory, f"Result_FEM_parameters_{frequency}", FEM_parameters)

    print(FEM_parameters)


    print(f"{len(results)=}")
    results = [item for item in FEM_parameters if ("FEM_results" not in item)]

    # Add total losses to each dictionary
    results = [dict(item, **{'total_losses': item["p_hyst_nom"] +
                                             np.sqrt(np.mean(np.array(item["j2H"])**2)) +
                                             np.sqrt(np.mean(np.array(item["j2F"])**2))}) for item in results]
    # Sort with ascending total losses
    results = sorted(results, key=itemgetter('total_losses'))

    # store_as_npy_in_directory(result_directory, f"results_{frequency}", results)
    print(f"{len(results)=}")
    print(f"{results=}")

    store_as_npy_in_directory(result_directory, f"Result_FEM_parameters_{frequency}", results)

    # find optimum
    if results:
        new_opt = results[0]['total_losses']
    else:
        print("Empty Results!")

        print(f"Rewrite last results")
        store_as_npy_in_directory(result_directory, f"Optimal_Result", init_opt)
        break

    # Check for converging
    if init_opt == new_opt:
        print("Converged!")

        print(f"Rewrite last results")
        store_as_npy_in_directory(result_directory, f"Optimal_Result", init_opt)
        break

    else:
        res_tmp = results_init_sorted + results
        print(f"{results_init_sorted=}")
        print(f"{results=}")
        print(f"{res_tmp=}")
        results_init_sorted = sorted(res_tmp, key=itemgetter('total_losses'))
