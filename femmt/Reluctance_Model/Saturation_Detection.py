import warnings
import numpy as np
import itertools
import matplotlib.pyplot as plt
from femmt import MagneticComponent
from femmt_functions import get_dict_with_unique_keys, get_dicts_with_keys_and_values, find_common_frequencies, \
    sort_out_small_harmonics
from DAB_Input_Data import working_points, reluctance_parameters, non_reluctance_parameters, L, N, power_nom, power_max
from scipy.interpolate import interp1d

mur = 3000
skin_accuracy = 0.5


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
    time_funct_in = dict_nom_in['wp_ib_il_vec']

    # Working point Output current data
    frequencies_out = dict_nom_out['wp_ob_il_frequency_vec']
    currents_out = dict_nom_out['wp_ob_il_amplitude_vec']
    phases_out = dict_nom_out['wp_ob_il_phase_rad_vec']
    time_funct_out = dict_nom_out['wp_ob_il_vec']

    # Find common frequencies
    current_pairs_nom, phase_pairs_nom, frequencies = find_common_frequencies(currents_in, phases_in, frequencies_in,
                                                                              currents_out, phases_out, frequencies_out)

    # Throw out all harmonics with low impact (f.e. amplitude smaller 1%)
    phase_pairs_nom, current_pairs_nom, frequencies = sort_out_small_harmonics(phase_pairs_nom, current_pairs_nom,
                                                                               frequencies, limes=0.01)

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

    time_funct_in = dict_max_in['wp_ib_il_vec']
    time_funct_out = dict_max_out['wp_ob_il_vec']

    # time_funct_in = dict_nom_in['wp_ib_il_vec']
    # time_funct_out = dict_nom_out['wp_ob_il_vec']


    I0 = 6  # 6 * 1.25  # 6 is peak of current wave_form

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    I1 = np.array(time_funct_in[1])  # current_pairs_nom[0]
    I2 = -np.array(time_funct_out[1])  # current_pairs_nom[1]
    # f_I1 = interp1d(time_funct_in[0], time_funct_in[1])  # current_pairs_nom[0]
    # f_I2 = interp1d(time_funct_in[0], time_funct_out[1])  # current_pairs_nom[1]

    # t = np.linspace(min(time_funct_in[0]), max(time_funct_in[0]), 100)
    t = time_funct_in[0]

    #I1 = f_I1(t)
    #I2 = -f_I2(t)

    print(f"{I1=}\n"
          f"{I2=}\n")

    Phi_top = []
    Phi_bot = []
    Phi_stray = []

    for i in range(0, len(I1)):
        I = [I1[i], I2[i]]
        [Phi_top_s, Phi_bot_s] = np.matmul(np.matmul(np.linalg.inv(np.transpose(N[0])), L), np.transpose(I))
        # Phi_stray_s = np.abs(Phi_top_s + Phi_bot_s)
        Phi_stray_s = Phi_bot_s - Phi_top_s

        Phi_top.append(Phi_top_s)
        Phi_bot.append(Phi_bot_s)
        Phi_stray.append(Phi_stray_s)

    print(f"{max(Phi_top)=}\n"
          f"{max(Phi_bot)=}\n"
          f"{max(Phi_stray)=}\n")

    # n = I2 / I1
    [Phi_top_single, Phi_bot_single] = np.matmul(np.matmul(np.linalg.inv(np.transpose(N[0])), L), np.transpose([-I0*1.25, 3.2*I0*1.25]))
    # Phi_stray_single = np.abs(Phi_top_s - Phi_bot_s)
    Phi_stray_single = Phi_bot_single - Phi_top_single

    print(f"{Phi_top_single=}\n"
          f"{Phi_bot_single=}\n"
          f"{Phi_stray_single=}\n")

    figure, axis = plt.subplots(3)


    axis[0].plot(t, Phi_top, label=f"Phi_top")
    axis[0].plot(t, Phi_bot, label=f"Phi_bot")
    axis[0].plot(t, Phi_stray, label=f"Phi_stray")
    axis[1].plot(time_funct_in[0], time_funct_in[1], label=f"Iin")
    axis[1].plot(time_funct_in[0], - np.array(time_funct_out[1]), label=f"Iout")
    # axis[1].plot(time_funct_in[0], np.array(time_funct_out[1])/3.2 - np.array(time_funct_in[1]), label=f"Im")
    axis[2].plot(t, I1, label=f"Iin")
    axis[2].plot(t, I2, label=f"Iout")
    # axis[2].plot(t, I1/3.2 + I2, label=f"Im")

    axis[0].set_ylabel("Magnetic flux in Wb")
    axis[0].set_xlabel("time in s")
    axis[1].set_ylabel("Current in A")
    axis[1].set_xlabel("time in s")
    axis[2].set_ylabel("Current in A")
    axis[2].set_xlabel("time in s")

    axis[0].legend()
    axis[1].legend()
    axis[2].legend()

    axis[0].grid()
    axis[1].grid()
    axis[2].grid()

    plt.show()


    #                                           -- Reluctance Model --
    # ------------------------------------------------------------------------------------------------------------------
    geo = MagneticComponent(component_type="integrated_transformer")

    valid_reluctance_parameters = geo.reluctance_model.air_gap_design(L_goal=L,
                                                                      parameters_init=reluctance_parameters,
                                                                      current=[-1.25 * I0, 1.25 * I0 * 3.2],
                                                                      b_max=0.3, b_stray=0.25,
                                                                      stray_path_parametrization="max_flux")

    input("Press Enter to continue...")
    #                                         -- FEM Simulation --
    # --------------------------------------------------------------------------------------------------------------
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

    geo.s = 1
    geo.ki = 0.53
    geo.alpha = 1.50
    geo.beta = 2.38

    for n_par, parameters in enumerate(FEM_parameters):


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
                              radius=geo.core.core_w / 2 + geo.core.window_w - parameters["R_stray"])

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
        #geo.get_inductances(I0=I0, op_frequency=frequency, skin_mesh_factor=skin_accuracy)


        geo.excitation_sweep(frequencies=frequencies,
                             currents=current_pairs_nom,
                             phi=phase_pairs_nom,
                             show_last=True)

