import numpy as np

from femmt import MagneticComponent
from femmt_functions import get_dict_with_unique_keys, get_dicts_with_keys_and_values, find_common_frequencies, \
    sort_out_small_harmonics, store_as_npy_in_directory
from DAB_Input_Data import working_points, reluctance_parameters, non_reluctance_parameters, L_goal, power_nom, power_max

result_directory = "C:/Users/tillp/OneDrive/Documents/GitHub/FEM_Magnetics_Toolbox/femmt/MA/final"


for wp_data in working_points:
    # Working point switching frequency
    frequency = dict_nom_in["frequency"]

    geo = MagneticComponent(component_type="integrated_transformer")



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
                              conductor_radii=[parameters["litzes"]["prim"]["conductor_radius"], parameters["litzes"]["sec"]["conductor_radius"]],
                              litz_para_type=['implicit_ff', 'implicit_ff'],
                              strands_numbers=[parameters["litzes"]["prim"]["N_strands"], parameters["litzes"]["sec"]["N_strands"]],
                              strand_radii=[parameters["litzes"]["prim"]["strand_radius"], parameters["litzes"]["sec"]["strand_radius"]],
                              cond_cond_isolation=[2*parameters["litzes"]["prim"]["isolation"],
                                                   2*parameters["litzes"]["sec"]["isolation"],
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

    store_as_npy_in_directory(result_directory, f"Result_FEM_parameters_{frequency}", FEM_parameters)

    print(FEM_parameters)
