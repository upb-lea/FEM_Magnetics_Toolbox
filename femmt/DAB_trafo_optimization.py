import numpy as np
from femmt import MagneticComponent
from femmt_functions import get_dict_with_unique_keys, get_dicts_with_keys_and_values, find_common_frequencies, \
    sort_out_small_harmonics
from DAB_Input_Data import working_points, reluctance_parameters, non_reluctance_parameters, L_goal, power_nom, power_max

skin_accuracy = 0.5
mur = 3000


#                                               -- Definitions --
# ----------------------------------------------------------------------------------------------------------------------
for wp_data in working_points:

    #                                           -- Get DAB Data--
    # ------------------------------------------------------------------------------------------------------------------

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Maximal Power Data (worst case)
    dicts_max = get_dicts_with_keys_and_values(wp_data, power=power_max)
    dict_max_in = get_dict_with_unique_keys(dicts_max, 'wp_ib_il_phase_rad_vec')
    dict_max_out = get_dict_with_unique_keys(dicts_max, 'wp_ob_il_phase_rad_vec')

    time_max_in = dict_max_in['wp_ib_il_vec']
    time_max_out = dict_max_out['wp_ob_il_vec']

    #                                           -- Reluctance Model --
    # ------------------------------------------------------------------------------------------------------------------
    geo = MagneticComponent(component_type="integrated_transformer")

    # valid_reluctance_parameters: List of dictionaries (Core Parameters, Stray Path Parameters, Winding Parameters,
    #                                                    Fluxes, Analytic Core Losses)
    valid_reluctance_parameters = geo.reluctance_model.air_gap_design(L_goal=L_goal[f"{frequency}"],
                                                                      parameters_init=reluctance_parameters,
                                                                      max_current=[time_max_in, time_max_out],
                                                                      nom_current=[time_nom_in, time_nom_out],
                                                                      nom_current_1st=[current_pairs_nom[0][0],
                                                                                       current_pairs_nom[0][1]],
                                                                      nom_phase_1st=[phase_pairs_nom[0][0],
                                                                                     phase_pairs_nom[0][1]],
                                                                      f_1st=frequency,
                                                                      b_max=0.3, b_stray=0.25,
                                                                      stray_path_parametrization="max_flux",
                                                                      visualize_waveforms=None)

    print(valid_reluctance_parameters)

    # input("Press Enter to continue...")

    #                                         -- FEM Simulation --
    # --------------------------------------------------------------------------------------------------------------
    # Bring together valid reluctance parameters and non reluctance parameters
    FEM_parameters = []
    for valid_parameters in valid_reluctance_parameters:
        for non_reluctance_parameter in non_reluctance_parameters:
            FEM_parameters.append(dict(valid_parameters, **non_reluctance_parameter))

    # Remove:
    geo.s = 0.5
    geo.ki = 0.53
    geo.alpha = 1.50
    geo.beta = 2.38

    for n_par, parameters in enumerate(FEM_parameters):

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
                              radius=geo.core.core_w / 2 + geo.core.window_w - parameters["R_stray"])

        geo.update_conductors(n_turns=[[parameters["N"][0, 0], parameters["N"][1, 0]],
                                       [parameters["N"][0, 1], parameters["N"][1, 1]]],
                              conductor_type=["litz", "litz"],
                              winding=["interleaved", "interleaved"],
                              scheme=["horizontal", "horizontal"],  # ["square", "square"]  "horizontal"
                              conductor_radii=[0.0005, 0.0009],
                              litz_para_type=['implicit_litz_radius', 'implicit_litz_radius'],
                              strands_numbers=[parameters["N_strands_prim"], parameters["N_strands_sec"]],
                              ff=[0.55, 0.55],
                              strand_radii=[parameters["strand_radius"], parameters["strand_radius"]],
                              cond_cond_isolation=[0.0001, 0.0001, 0.0001],
                              core_cond_isolation=[0.0005, 0.0005])

        # -- Simulation --
        I0 = 6  # 6 * 1.25  # 6 is peak of current wave_form
        # geo.get_inductances(I0=I0, op_frequency=frequency, skin_mesh_factor=skin_accuracy)

        # geo.single_simulation(freq=frequencies[0],
        #                       current=current_pairs_nom[0],
        #                       phi_deg=phase_pairs_nom[0])

        FEM_results = geo.excitation_sweep(frequencies=frequencies,
                                           currents=current_pairs_nom,
                                           phi=phase_pairs_nom,
                                           show_last=True,
                                           return_results=True)

        FEM_parameters[n_par] = dict(FEM_parameters[n_par], **FEM_results)

    print(FEM_parameters)
