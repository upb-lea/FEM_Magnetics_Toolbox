"""Stacked transformer optimization."""
# python libraries
import datetime
import logging
import os
import shutil
import json
import pickle

# own libraries
from femmt.optimization.sto_dtos import StoSingleInputConfig, StoTargetAndFixedParameters, FemInput, FemOutput, ReluctanceModelInput, ReluctanceModelOutput
import femmt.functions as ff
import femmt.functions_reluctance as fr
import materialdatabase as mdb
import femmt.optimization.ito_functions as itof
import femmt.optimization.functions_optimization as fo
import femmt as fmt
from materialdatabase.meta.data_enums import Material

# 3rd party libraries
import numpy as np
import optuna
from scipy import optimize
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import magnethub as mh
import pandas as pd
import tqdm

logger = logging.getLogger(__name__)

class StackedTransformerOptimization:
    """Perform optimizations for the stacked transformer."""

    class ReluctanceModel:
        """Perform optimizations using the reluctance model of the stacked transformer."""

        @staticmethod
        def calculate_fix_parameters(config: StoSingleInputConfig) -> StoTargetAndFixedParameters:
            """Calculate fix parameters what can be derived from the input configuration.

            return values are:

                i_rms_1
                i_rms_2
                time_extracted_vec
                current_extracted_1_vec
                current_extracted_2_vec
                material_mu_r_abs_list
                fundamental_frequency
                target_inductance_matrix
                fem_working_directory
                fem_simulation_results_directory
                reluctance_model_results_directory
                fem_thermal_simulation_results_directory

            :param config: configuration file
            :type config: ItoSingleInputConfig
            :return: calculated target and fix parameters
            :rtype: ItoTargetAndFixedParameters
            """
            # currents
            time_extracted, current_extracted_1_vec = fr.time_vec_current_vec_from_time_current_vec(
                config.time_current_1_vec)
            time_extracted, current_extracted_2_vec = fr.time_vec_current_vec_from_time_current_vec(
                config.time_current_2_vec)
            fundamental_frequency = 1 / time_extracted[-1]

            i_rms_1 = fr.i_rms(config.time_current_1_vec)
            i_rms_2 = fr.i_rms(config.time_current_2_vec)

            i_peak_1, i_peak_2 = fr.max_value_from_value_vec(current_extracted_1_vec, current_extracted_2_vec)
            phi_deg_1, phi_deg_2 = fr.phases_deg_from_time_current(time_extracted, current_extracted_1_vec,
                                                                   current_extracted_2_vec)

            # phi_deg_2 = phi_deg_2 - 180

            # target inductances
            target_inductance_matrix = fr.calculate_inductance_matrix_from_ls_lh_n(config.l_s12_target,
                                                                                   config.l_h_target,
                                                                                   config.n_target)

            (fft_frequencies_1, fft_amplitudes_1, fft_phases_1) = ff.fft(
                period_vector_t_i=config.time_current_1_vec, sample_factor=1000, plot='no', mode='time', filter_type='factor', filter_value_factor=0.03)

            (fft_frequencies_2, fft_amplitudes_2, fft_phases_2) = ff.fft(
                period_vector_t_i=config.time_current_2_vec, sample_factor=1000, plot='no', mode='time', filter_type='factor', filter_value_factor=0.03)

            # material properties
            material_db = mdb.Data()

            material_mu_r_abs_list = []
            magnet_model_list = []
            for material_name in config.material_list:
                small_signal_mu_real_over_f_at_T = material_db.get_datasheet_curve(material_name, mdb.DatasheetCurveType.small_signal_mu_real_over_f_at_T)
                small_signal_mu_imag_over_f_at_T = material_db.get_datasheet_curve(material_name, mdb.DatasheetCurveType.small_signal_mu_imag_over_f_at_T)

                mu_real_at_f = np.interp(10_000, small_signal_mu_real_over_f_at_T["f"], small_signal_mu_real_over_f_at_T["mu_real"])
                mu_imag_at_f = np.interp(10_000, small_signal_mu_imag_over_f_at_T["f"], small_signal_mu_imag_over_f_at_T["mu_imag"])

                material_mu_r_abs = np.abs(np.array([mu_real_at_f + 1j * mu_imag_at_f]))
                material_mu_r_abs_list.append(material_mu_r_abs)

                # instantiate material-specific model
                mdl: mh.loss.LossModel = mh.loss.LossModel(material=material_name, team="paderborn")
                magnet_model_list.append(mdl)

            # set up working directories
            working_directories = itof.set_up_folder_structure(config.stacked_transformer_optimization_directory)

            # finalize data to dto
            target_and_fix_parameters = StoTargetAndFixedParameters(
                i_rms_1=i_rms_1,
                i_rms_2=i_rms_2,
                i_peak_1=i_peak_1,
                i_peak_2=i_peak_2,
                i_phase_deg_1=phi_deg_1,
                i_phase_deg_2=phi_deg_2,
                time_extracted_vec=time_extracted,
                magnet_hub_model_list=magnet_model_list,
                current_extracted_1_vec=current_extracted_1_vec,
                current_extracted_2_vec=current_extracted_2_vec,
                material_name_list=config.material_list,
                material_complex_mu_r_list=material_mu_r_abs_list,
                fundamental_frequency=fundamental_frequency,
                target_inductance_matrix=target_inductance_matrix,
                working_directories=working_directories,
                # winding 1
                fft_frequency_list_1=fft_frequencies_1,
                fft_amplitude_list_1=fft_amplitudes_1,
                fft_phases_list_1=fft_phases_1,

                # winding 2
                fft_frequency_list_2=fft_frequencies_2,
                fft_amplitude_list_2=fft_amplitudes_2,
                fft_phases_list_2=fft_phases_2
            )

            return target_and_fix_parameters

        @staticmethod
        def objective(trial: optuna.Trial, config: StoSingleInputConfig, target_and_fixed_parameters: StoTargetAndFixedParameters):
            """
            Objective function to optimize. Uses reluctance model calculation.

            Once core_name_list is not None, the objective function uses fixed core sizes. Cores are picked from the core_database().
            Otherwise, core_inner_diameter_min_max_list, window_w_min_max_list and window_h_bot_min_max_list are used.

            :param trial: Optuna trial
            :type trial: optuna.Trial
            :param config: Stacked transformer optimization configuration file
            :type config: StoSingleInputConfig
            :param target_and_fixed_parameters: Target and fixed parameters
            :type target_and_fixed_parameters: StoTargetAndFixedParameters
            """
            # fixed cores
            if config.core_name_list is not None:
                # using fixed core sizes from the database with flexible height.
                core_name = trial.suggest_categorical("core_name", config.core_name_list)
                core = ff.core_database()[core_name]
                core_inner_diameter = core["core_inner_diameter"]
                window_w = core["window_w"]
                window_h_bot = trial.suggest_float("window_h_bot", 0.3 * core["window_h"], core["window_h"])
                trial.set_user_attr('core_inner_diameter', core_inner_diameter)
                trial.set_user_attr('window_w', window_w)
                trial.set_user_attr('window_h_bot', window_h_bot)

            else:
                # using arbitrary core sizes
                core_inner_diameter = trial.suggest_float("core_inner_diameter", config.core_inner_diameter_min_max_list[0],
                                                          config.core_inner_diameter_min_max_list[1])
                window_w = trial.suggest_float("window_w", config.window_w_min_max_list[0], config.window_w_min_max_list[1])
                window_h_bot = trial.suggest_float('window_h_bot', config.window_h_bot_min_max_list[0], config.window_h_bot_min_max_list[1])

            n_p_top = trial.suggest_int('n_p_top', config.n_p_top_min_max_list[0], config.n_p_top_min_max_list[1])
            n_p_bot = trial.suggest_int('n_p_bot', config.n_p_bot_min_max_list[0], config.n_p_bot_min_max_list[1])
            n_s_bot_min = int(np.round(n_p_bot / config.n_target, decimals=0)) - 2
            if n_s_bot_min < 1:
                n_s_bot_min = 1
            n_s_bot_max = int(np.round(n_p_bot / config.n_target, decimals=0)) + 2
            n_s_bot = trial.suggest_int('n_s_bot', n_s_bot_min, n_s_bot_max)

            primary_litz_name = trial.suggest_categorical('primary_litz_name', config.primary_litz_wire_list)
            primary_litz_dict = ff.litz_database()[primary_litz_name]
            primary_litz_diameter = 2 * primary_litz_dict['conductor_radii']

            secondary_litz_name = trial.suggest_categorical('secondary_litz_name', config.secondary_litz_wire_list)
            secondary_litz_dict = ff.litz_database()[secondary_litz_name]
            secondary_litz_diameter = 2 * secondary_litz_dict['conductor_radii']

            material_name = trial.suggest_categorical('material_name', config.material_list)
            for count, material_mu_r_abs_value in enumerate(target_and_fixed_parameters.material_complex_mu_r_list):
                if target_and_fixed_parameters.material_name_list[count] == material_name:
                    material_mu_r_abs: mdb.MaterialCurve = material_mu_r_abs_value
                    magnet_material_model = target_and_fixed_parameters.magnet_hub_model_list[count]

            # calculate total 2D-axi symmetric volume of the core:
            # formula: number_turns_per_row = (available_width + primary_to_primary) / (wire_diameter + primary_to_primary)
            available_width_top = window_w - config.insulations.iso_window_top_core_left - config.insulations.iso_window_top_core_right
            possible_number_turns_per_column_top_window = int(
                (available_width_top + config.insulations.iso_primary_to_primary) / (primary_litz_diameter + config.insulations.iso_primary_to_primary))
            if possible_number_turns_per_column_top_window < 1:
                return float('nan'), float('nan')
            number_of_rows_needed = np.ceil(n_p_top / possible_number_turns_per_column_top_window)
            needed_height_top_wo_insulation = (number_of_rows_needed * primary_litz_diameter + \
                                               (number_of_rows_needed - 1) * config.insulations.iso_primary_to_primary)
            window_h_top = needed_height_top_wo_insulation + config.insulations.iso_window_top_core_top + config.insulations.iso_window_top_core_bot

            # detailed calculation for the winding window
            # check the area for the primary winding
            available_height_bot = window_h_bot - config.insulations.iso_window_bot_core_top - config.insulations.iso_window_bot_core_bot
            possible_number_prim_turns_per_column_bot_window = int(
                (available_height_bot + config.insulations.iso_primary_to_primary) / (primary_litz_diameter + config.insulations.iso_primary_to_primary))
            if possible_number_prim_turns_per_column_bot_window < 1:
                return float('nan'), float('nan')
            number_of_primary_columns_needed = np.ceil(n_p_bot / possible_number_prim_turns_per_column_bot_window)
            needed_primary_width_bot_wo_insulation = (number_of_primary_columns_needed * primary_litz_diameter + (number_of_primary_columns_needed - 1) * \
                                                      config.insulations.iso_primary_to_primary)
            area_primary_bot = needed_primary_width_bot_wo_insulation * window_h_bot

            # check the area for the secondary winding
            possible_number_sec_turns_per_column_bot_window = int(
                (available_height_bot + config.insulations.iso_secondary_to_secondary) / \
                (secondary_litz_diameter + config.insulations.iso_secondary_to_secondary))
            if possible_number_sec_turns_per_column_bot_window < 1:
                return float('nan'), float('nan')
            number_of_secondary_columns_needed = np.ceil(n_s_bot / possible_number_sec_turns_per_column_bot_window)
            needed_primary_width_bot_wo_insulation = (number_of_secondary_columns_needed * secondary_litz_diameter + \
                                                      (number_of_secondary_columns_needed - 1) * config.insulations.iso_secondary_to_secondary)
            area_secondary_bot = needed_primary_width_bot_wo_insulation * window_h_bot
            area_insulation_prim_sec_bot = 2 * config.insulations.iso_primary_to_secondary * window_h_bot

            total_area_windings_bot = area_primary_bot + area_secondary_bot + area_insulation_prim_sec_bot

            window_bot_available_height = window_h_bot - config.insulations.iso_window_bot_core_top - config.insulations.iso_window_bot_core_bot
            window_bot_available_width = window_w - config.insulations.iso_window_bot_core_left - config.insulations.iso_window_bot_core_right
            window_bot_available_area = window_bot_available_height * window_bot_available_width

            # as the window_h_top is adapted to the number of n_p_top, the top windings always fit into the top window.
            if total_area_windings_bot > window_bot_available_area:
                logger.debug("Winding window too small for too many turns.")
                return float('nan'), float('nan')

            reluctance_model_intput = ReluctanceModelInput(
                target_inductance_matrix=target_and_fixed_parameters.target_inductance_matrix,
                core_inner_diameter=core_inner_diameter,
                window_w=window_w,
                window_h_bot=window_h_bot,
                window_h_top=window_h_top,
                turns_1_top=n_p_top,
                turns_1_bot=n_p_bot,
                turns_2_bot=n_s_bot,
                litz_wire_name_1=primary_litz_name,
                litz_wire_diameter_1=primary_litz_diameter,
                litz_wire_name_2=secondary_litz_name,
                litz_wire_diameter_2=secondary_litz_diameter,

                insulations=config.insulations,
                material_name=material_name,
                material_mu_r_abs=material_mu_r_abs,
                magnet_material_model=magnet_material_model,

                temperature=config.temperature,
                time_extracted_vec=target_and_fixed_parameters.time_extracted_vec,
                current_extracted_vec_1=target_and_fixed_parameters.current_extracted_1_vec,
                current_extracted_vec_2=target_and_fixed_parameters.current_extracted_2_vec,
                fundamental_frequency=target_and_fixed_parameters.fundamental_frequency,
                i_rms_1=target_and_fixed_parameters.i_rms_1,
                i_rms_2=target_and_fixed_parameters.i_rms_2,
                primary_litz_dict=primary_litz_dict,
                secondary_litz_dict=secondary_litz_dict,

                # winding 1
                fft_frequency_list_1=target_and_fixed_parameters.fft_frequency_list_1,
                fft_amplitude_list_1=target_and_fixed_parameters.fft_amplitude_list_1,
                fft_phases_list_1=target_and_fixed_parameters.fft_phases_list_1,
                # winding 2
                fft_frequency_list_2=target_and_fixed_parameters.fft_frequency_list_2,
                fft_amplitude_list_2=target_and_fixed_parameters.fft_amplitude_list_2,
                fft_phases_list_2=target_and_fixed_parameters.fft_phases_list_2,
            )
            try:
                reluctance_output = StackedTransformerOptimization.ReluctanceModel.single_reluctance_model_simulation(reluctance_model_intput)
            except ValueError:
                logger.debug("No fitting air gap length")
                return float('nan'), float('nan')

            # set additional attributes
            trial.set_user_attr('p_hyst', reluctance_output.p_hyst)
            trial.set_user_attr('p_hyst_top', reluctance_output.p_hyst_top)
            trial.set_user_attr('p_hyst_bot', reluctance_output.p_hyst_bot)
            trial.set_user_attr('p_hyst_middle', reluctance_output.p_hyst_middle)
            trial.set_user_attr('b_max_top', reluctance_output.b_max_top)
            trial.set_user_attr('b_max_bot', reluctance_output.b_max_bot)
            trial.set_user_attr('b_max_middle', reluctance_output.b_max_middle)
            trial.set_user_attr('window_h_top', window_h_top)
            trial.set_user_attr('winding_losses', reluctance_output.winding_1_loss + reluctance_output.winding_2_loss)
            trial.set_user_attr('l_top_air_gap', reluctance_output.l_top_air_gap)
            trial.set_user_attr('l_bot_air_gap', reluctance_output.l_bot_air_gap)

            return reluctance_output.volume, reluctance_output.p_loss

        @staticmethod
        def single_reluctance_model_simulation(reluctance_input: ReluctanceModelInput) -> ReluctanceModelOutput:
            """
            Perform a single reluctance model simulation.

            :param reluctance_input: Reluctance model input data.
            :type reluctance_input: ReluctanceModelInput
            :return: Reluctance model output data
            :rtype: ReluctanceModelOutput
            """
            core_total_height = reluctance_input.window_h_top + reluctance_input.window_h_bot + reluctance_input.core_inner_diameter * 3 / 4
            r_outer = fr.calculate_r_outer(reluctance_input.core_inner_diameter, reluctance_input.window_w)
            volume = ff.calculate_cylinder_volume(cylinder_diameter=2 * r_outer, cylinder_height=core_total_height)

            # calculate the reluctance and flux matrix
            winding_matrix = np.array([[reluctance_input.turns_1_top, 0], [reluctance_input.turns_1_bot, reluctance_input.turns_2_bot]])

            reluctance_matrix = fr.calculate_reluctance_matrix(winding_matrix, reluctance_input.target_inductance_matrix)

            current_matrix = np.array([reluctance_input.current_extracted_vec_1, reluctance_input.current_extracted_vec_2])

            flux_matrix = fr.calculate_flux_matrix(reluctance_matrix, winding_matrix, current_matrix)

            flux_top = flux_matrix[0]
            flux_bot = flux_matrix[1]
            flux_middle = flux_bot - flux_top

            core_cross_section = (reluctance_input.core_inner_diameter / 2) ** 2 * np.pi

            flux_density_top = flux_top / core_cross_section
            flux_density_bot = flux_bot / core_cross_section
            flux_density_middle = flux_middle / core_cross_section

            # calculate the core reluctance
            core_inner_cylinder_top = fr.r_core_round(reluctance_input.core_inner_diameter, reluctance_input.window_h_top,
                                                      reluctance_input.material_mu_r_abs)
            core_inner_cylinder_bot = fr.r_core_round(reluctance_input.core_inner_diameter, reluctance_input.window_h_bot,
                                                      reluctance_input.material_mu_r_abs)
            core_top_bot_radiant = fr.r_core_top_bot_radiant(reluctance_input.core_inner_diameter, reluctance_input.window_w,
                                                             reluctance_input.material_mu_r_abs, reluctance_input.core_inner_diameter / 4)

            r_core_top = 2 * core_inner_cylinder_top + core_top_bot_radiant
            r_core_bot = 2 * core_inner_cylinder_bot + core_top_bot_radiant
            r_core_middle = core_top_bot_radiant

            r_core_top_target = reluctance_matrix[0][0] + reluctance_matrix[0][1]  # - r_core_middle
            r_core_bot_target = reluctance_matrix[1][1] + reluctance_matrix[0][1]  # - r_core_middle

            r_air_gap_top_target = r_core_top_target - r_core_top
            r_air_gap_bot_target = r_core_bot_target - r_core_bot

            # calculate air gaps to reach the target parameters
            minimum_air_gap_length = 0.01e-3
            maximum_air_gap_length = 2e-3

            l_top_air_gap = optimize.brentq(
                fr.r_air_gap_round_inf_sct, minimum_air_gap_length, maximum_air_gap_length,
                args=(reluctance_input.core_inner_diameter, reluctance_input.window_h_top, r_air_gap_top_target), full_output=True)[0]

            l_bot_air_gap = optimize.brentq(
                fr.r_air_gap_round_round_sct, minimum_air_gap_length, maximum_air_gap_length,
                args=(reluctance_input.core_inner_diameter, reluctance_input.window_h_bot / 2, reluctance_input.window_h_bot / 2, r_air_gap_bot_target),
                full_output=True)[0]

            # calculate hysteresis losses from mag-net-hub
            interp_points = np.arange(0, 1024) * reluctance_input.time_extracted_vec[-1] / 1024
            flux_density_top_interp = np.interp(interp_points, reluctance_input.time_extracted_vec, flux_density_top)
            flux_density_bot_interp = np.interp(interp_points, reluctance_input.time_extracted_vec, flux_density_bot)
            flux_density_middle_interp = np.interp(interp_points, reluctance_input.time_extracted_vec, flux_density_middle)

            # get power loss in W/m³ and estimated H wave in A/m
            p_density_top, _ = reluctance_input.magnet_material_model(flux_density_top_interp,
                                                                      reluctance_input.fundamental_frequency, reluctance_input.temperature)
            p_density_bot, _ = reluctance_input.magnet_material_model(flux_density_bot_interp,
                                                                      reluctance_input.fundamental_frequency, reluctance_input.temperature)
            p_density_middle, _ = reluctance_input.magnet_material_model(flux_density_middle_interp,
                                                                         reluctance_input.fundamental_frequency, reluctance_input.temperature)

            volume_core_top = (2 * ff.calculate_cylinder_volume(reluctance_input.core_inner_diameter, reluctance_input.window_h_top) - \
                               ff.calculate_cylinder_volume(reluctance_input.core_inner_diameter, l_top_air_gap) + \
                               ff.calculate_cylinder_volume(2 * r_outer, reluctance_input.core_inner_diameter / 4))
            volume_core_bot = (2 * ff.calculate_cylinder_volume(reluctance_input.core_inner_diameter, reluctance_input.window_h_bot) - \
                               ff.calculate_cylinder_volume(reluctance_input.core_inner_diameter, l_bot_air_gap) + \
                               ff.calculate_cylinder_volume(2 * r_outer, reluctance_input.core_inner_diameter / 4))
            volume_core_middle = ff.calculate_cylinder_volume(2 * r_outer, reluctance_input.core_inner_diameter / 4)

            p_top = p_density_top * volume_core_top
            p_bot = p_density_bot * volume_core_bot
            p_middle = p_density_middle * volume_core_middle

            p_hyst = p_top + p_bot + p_middle

            # calculate winding losses
            primary_effective_conductive_cross_section = (
                reluctance_input.primary_litz_dict["strands_numbers"] * reluctance_input.primary_litz_dict["strand_radii"] ** 2 * np.pi)
            primary_effective_conductive_radius = np.sqrt(primary_effective_conductive_cross_section / np.pi)
            primary_resistance_top = fr.resistance_solid_wire(
                reluctance_input.core_inner_diameter, reluctance_input.window_w, reluctance_input.turns_1_top,
                primary_effective_conductive_radius, material='Copper')

            number_bot_prim_turns_per_column = (
                int((reluctance_input.window_h_bot - reluctance_input.insulations.iso_window_bot_core_top - \
                     reluctance_input.insulations.iso_window_bot_core_bot + reluctance_input.insulations.iso_primary_to_primary) / \
                    (2 * reluctance_input.primary_litz_dict["conductor_radii"] + reluctance_input.insulations.iso_primary_to_primary)))
            if number_bot_prim_turns_per_column > reluctance_input.turns_1_bot:
                # single row window only
                primary_resistance_bot_inner = fr.resistance_solid_wire(
                    reluctance_input.core_inner_diameter, reluctance_input.window_w, reluctance_input.turns_1_bot,
                    primary_effective_conductive_radius, material='Copper')
                primary_resistance_bot_outer = 0
            else:
                # multiple row window
                primary_resistance_bot_inner = fr.resistance_solid_wire(
                    reluctance_input.core_inner_diameter, reluctance_input.window_w, number_bot_prim_turns_per_column,
                    primary_effective_conductive_radius, material='Copper')

                primary_resistance_bot_outer = fr.resistance_solid_wire(
                    reluctance_input.core_inner_diameter, reluctance_input.window_w, reluctance_input.turns_1_bot - number_bot_prim_turns_per_column,
                    primary_effective_conductive_radius, material='Copper')

            secondary_effective_conductive_cross_section = (
                reluctance_input.secondary_litz_dict["strands_numbers"] * reluctance_input.secondary_litz_dict["strand_radii"] ** 2 * np.pi)
            secondary_effective_conductive_radius = np.sqrt(secondary_effective_conductive_cross_section / np.pi)
            secondary_resistance = fr.resistance_solid_wire(
                reluctance_input.core_inner_diameter, reluctance_input.window_w, reluctance_input.turns_2_bot, secondary_effective_conductive_radius,
                material='Copper')

            winding_area_1_top = (
                2 * reluctance_input.primary_litz_dict["conductor_radii"] * \
                (reluctance_input.turns_1_top * 2 * reluctance_input.primary_litz_dict["conductor_radii"] + \
                    (reluctance_input.turns_1_top - 1) * reluctance_input.insulations.iso_primary_to_primary))

            p_winding_1_top = 0
            p_winding_1_bot = 0
            for count, fft_frequency in enumerate(reluctance_input.fft_frequency_list_1):
                proximity_factor_1_top = fr.calc_proximity_factor_air_gap(
                    litz_wire_name=reluctance_input.litz_wire_name_1, number_turns=reluctance_input.turns_1_top,
                    r_1=reluctance_input.insulations.iso_window_top_core_left,
                    frequency=fft_frequency, winding_area=winding_area_1_top,
                    litz_wire_material_name='Copper', temperature=reluctance_input.temperature)

                p_winding_1_top += proximity_factor_1_top * primary_resistance_top * reluctance_input.fft_amplitude_list_1[count] ** 2

                if number_bot_prim_turns_per_column > reluctance_input.turns_1_bot:
                    winding_area_1_bot = 2 * reluctance_input.primary_litz_dict["conductor_radii"] * \
                        (reluctance_input.turns_1_bot * 2 * reluctance_input.primary_litz_dict["conductor_radii"] + \
                            (reluctance_input.turns_1_bot - 1) * reluctance_input.insulations.iso_primary_to_primary)

                    proximity_factor_1_bot_inner = fr.calc_proximity_factor_air_gap(
                        litz_wire_name=reluctance_input.litz_wire_name_1, number_turns=reluctance_input.turns_1_bot,
                        r_1=reluctance_input.insulations.iso_window_bot_core_left,
                        frequency=fft_frequency, winding_area=winding_area_1_bot,
                        litz_wire_material_name='Copper', temperature=reluctance_input.temperature)

                    proximity_factor_1_bot_outer = 0
                else:
                    winding_area_1_bot = (2 * reluctance_input.primary_litz_dict["conductor_radii"] * (
                        number_bot_prim_turns_per_column * 2 * reluctance_input.primary_litz_dict["conductor_radii"] + \
                        (number_bot_prim_turns_per_column - 1) * reluctance_input.insulations.iso_primary_to_primary))

                    proximity_factor_1_bot_inner = fr.calc_proximity_factor_air_gap(
                        litz_wire_name=reluctance_input.litz_wire_name_1, number_turns=number_bot_prim_turns_per_column,
                        r_1=reluctance_input.insulations.iso_window_bot_core_left,
                        frequency=fft_frequency, winding_area=winding_area_1_bot,
                        litz_wire_material_name='Copper', temperature=reluctance_input.temperature)

                    proximity_factor_1_bot_outer = fr.calc_proximity_factor(
                        litz_wire_name=reluctance_input.litz_wire_name_1, number_turns=reluctance_input.turns_1_bot - number_bot_prim_turns_per_column,
                        window_h=reluctance_input.window_h_bot,
                        iso_core_top=reluctance_input.insulations.iso_window_bot_core_top, iso_core_bot=reluctance_input.insulations.iso_window_bot_core_bot,
                        frequency=fft_frequency, litz_wire_material_name='Copper', temperature=reluctance_input.temperature)

                p_winding_1_bot_inner = proximity_factor_1_bot_inner * primary_resistance_bot_inner * reluctance_input.fft_amplitude_list_1[count] ** 2
                p_winding_1_bot_outer = proximity_factor_1_bot_outer * primary_resistance_bot_outer * reluctance_input.fft_amplitude_list_1[count] ** 2

                p_winding_1_bot += p_winding_1_bot_inner + p_winding_1_bot_outer

            p_winding_2 = 0
            for count, fft_frequency in enumerate(reluctance_input.fft_frequency_list_2):
                proximity_factor_assumption_2 = fr.calc_proximity_factor(
                    litz_wire_name=reluctance_input.litz_wire_name_2, number_turns=reluctance_input.turns_2_bot, window_h=reluctance_input.window_h_bot,
                    iso_core_top=reluctance_input.insulations.iso_window_bot_core_top, iso_core_bot=reluctance_input.insulations.iso_window_bot_core_bot,
                    frequency=fft_frequency, litz_wire_material_name='Copper', temperature=reluctance_input.temperature)

                p_winding_2 += proximity_factor_assumption_2 * secondary_resistance * reluctance_input.fft_amplitude_list_2[count] ** 2

            p_loss_total = p_hyst + p_winding_1_top + p_winding_1_bot + p_winding_2

            area_to_heat_sink = r_outer ** 2 * np.pi

            reluctance_output = ReluctanceModelOutput(
                # set additional attributes
                p_hyst=p_hyst,
                p_hyst_top=p_top,
                p_hyst_bot=p_bot,
                p_hyst_middle=p_middle,
                b_max_top=np.max(flux_density_top_interp),
                b_max_bot=np.max(flux_density_bot_interp),
                b_max_middle=np.max(flux_density_middle_interp),
                winding_1_loss=p_winding_1_top + p_winding_1_bot,
                winding_2_loss=p_winding_2,
                l_top_air_gap=l_top_air_gap,
                l_bot_air_gap=l_bot_air_gap,
                volume=volume,
                area_to_heat_sink=area_to_heat_sink,
                p_loss=p_loss_total,
            )

            return reluctance_output

        @staticmethod
        def start_proceed_study(config: StoSingleInputConfig, number_trials: int | None = None,
                                target_number_trials: int | None = None, storage: str = 'sqlite',
                                sampler=optuna.samplers.NSGAIIISampler(),
                                ) -> None:
            """
            Proceed a study which is stored as sqlite database.

            :param config: Simulation configuration
            :type config: ItoSingleInputConfig
            :param number_trials: Number of trials adding to the existing study
            :type number_trials: int
            :param storage: storage database, e.g. 'sqlite' or 'mysql'
            :type storage: str
            :param target_number_trials: Number of target trials for the existing study
            :type target_number_trials: int
            :param sampler: optuna.samplers.NSGAIISampler() or optuna.samplers.NSGAIIISampler(). Note about the brackets () !!
            :type sampler: optuna.sampler-object
            """
            if os.path.exists(f"{config.stacked_transformer_optimization_directory}/{config.stacked_transformer_study_name}.sqlite3"):
                logger.info("Existing study found. Proceeding.")

            target_and_fixed_parameters = StackedTransformerOptimization.ReluctanceModel.calculate_fix_parameters(config)

            # introduce study in storage, e.g. sqlite or mysql
            if storage == 'sqlite':
                # Note: for sqlite operation, there needs to be three slashes '///' even before the path '/home/...'
                # Means, in total there are four slashes including the path itself '////home/.../database.sqlite3'
                storage = f"sqlite:///{config.stacked_transformer_optimization_directory}/{config.stacked_transformer_study_name}.sqlite3"
            elif storage == 'mysql':
                storage = "mysql://monty@localhost/mydb",

            # set logging verbosity: https://optuna.readthedocs.io/en/stable/reference/generated/optuna.logging.set_verbosity.html#optuna.logging.set_verbosity
            # .INFO: all messages (default)
            # .WARNING: fails and warnings
            # .ERROR: only errors
            optuna.logging.set_verbosity(optuna.logging.ERROR)

            func = lambda trial: StackedTransformerOptimization.ReluctanceModel.objective(trial, config, target_and_fixed_parameters)

            study_in_storage = optuna.create_study(study_name=config.stacked_transformer_study_name,
                                                   storage=storage,
                                                   directions=['minimize', 'minimize'],
                                                   load_if_exists=True, sampler=sampler)

            if target_number_trials is not None:
                # simulation for a given number of target trials
                if len(study_in_storage.trials) < target_number_trials:
                    study_in_memory = optuna.create_study(directions=['minimize', 'minimize'],
                                                          study_name=config.stacked_transformer_study_name, sampler=sampler)
                    logger.info(f"Sampler is {study_in_memory.sampler.__class__.__name__}")
                    study_in_memory.add_trials(study_in_storage.trials)
                    number_trials = target_number_trials - len(study_in_memory.trials)
                    study_in_memory.optimize(func, n_trials=number_trials, show_progress_bar=True)
                    study_in_storage.add_trials(study_in_memory.trials[-number_trials:])
                    logger.info(f"Finished {number_trials} trials.")
                    logger.info(f"current time: {datetime.datetime.now()}")
                else:
                    logger.info(f"Study has already {len(study_in_storage.trials)} trials, and target is {target_number_trials} trials.")

            else:
                # normal simulation with number_trials
                study_in_memory = optuna.create_study(directions=['minimize', 'minimize'], study_name=config.stacked_transformer_study_name, sampler=sampler)
                logger.info(f"Sampler is {study_in_memory.sampler.__class__.__name__}")
                study_in_memory.add_trials(study_in_storage.trials)
                study_in_memory.optimize(func, n_trials=number_trials, show_progress_bar=True)

                study_in_storage.add_trials(study_in_memory.trials[-number_trials:])
                logger.info(f"Finished {number_trials} trials.")
                logger.info(f"current time: {datetime.datetime.now()}")
            StackedTransformerOptimization.ReluctanceModel.save_config(config)

        @staticmethod
        def show_study_results(config: StoSingleInputConfig) -> None:
            """Show the results of a study.

            A local .html file is generated under config.working_directory to store the interactive plotly plots on disk.

            :param config: Integrated transformer configuration file
            :type config: ItoSingleInputConfig
            """
            study = optuna.load_study(study_name=config.stacked_transformer_study_name,
                                      storage=f"sqlite:///{config.stacked_transformer_optimization_directory}/{config.stacked_transformer_study_name}.sqlite3")

            fig = optuna.visualization.plot_pareto_front(study, targets=lambda t: (t.values[0], t.values[1]), target_names=["volume in m³", "loss in W"])
            fig.update_layout(title=f"{config.stacked_transformer_study_name} <br><sup>{config.stacked_transformer_optimization_directory}</sup>")
            fig.write_html(f"{config.stacked_transformer_optimization_directory}/{config.stacked_transformer_study_name}"
                           f"_{datetime.datetime.now().isoformat(timespec='minutes')}.html")
            fig.show()

        @staticmethod
        def study_to_df(config: StoSingleInputConfig) -> pd.DataFrame:
            """
            Create a Pandas dataframe from a study.

            :param config: configuration
            :type config: InductorOptimizationDTO
            :return: Study results as Pandas Dataframe
            :rtype: pd.DataFrame
            """
            database_url = f'sqlite:///{config.stacked_transformer_optimization_directory}/{config.stacked_transformer_study_name}.sqlite3'
            if os.path.isfile(database_url.replace('sqlite:///', '')):
                logger.info("Existing study found.")
            else:
                raise ValueError(f"Can not find database: {database_url}")
            loaded_study = optuna.load_study(study_name=config.stacked_transformer_study_name, storage=database_url)
            df = loaded_study.trials_dataframe()
            df.to_csv(f'{config.stacked_transformer_optimization_directory}/{config.stacked_transformer_study_name}.csv')
            logger.info(f"Exported study as .csv file: {config.stacked_transformer_optimization_directory}/{config.stacked_transformer_study_name}.csv")
            return df

        @staticmethod
        def load_csv_to_df(csv_filepath: str) -> pd.DataFrame:
            """
            Load a csv file (previously stored from a Pandas dataframe) back to a Pandas dataframe.

            :param csv_filepath: File path of .csv file
            :type csv_filepath: str
            :return: loaded results from the given .csv file
            :rtype: pandas.DataFrame
            """
            df = pd.read_csv(csv_filepath, header=0, index_col=0)
            # reading a pandas dataframe seems to change a global variable in the c subsystem
            # after reading csv values, there are issues running onelab/gmsh, as gmsh writes ',' instead '.' to its own files
            # reading the file again with setting back the delimiter to ';', is a workaround for the mentioned problem.
            pd.read_csv(csv_filepath, header=0, index_col=0, delimiter=';')
            return df

        @staticmethod
        def df_plot_pareto_front(*dataframes: list[pd.DataFrame], label_list: list[str], color_list: list[str] = None,
                                 interactive: bool = True) -> None:
            """
            Plot an interactive Pareto diagram (losses vs. volume) to select the transformers to re-simulate.

            :param dataframes: Dataframe, generated from an optuna study (exported by optuna)
            :type dataframes: pd.Dataframe
            :param label_list: list of labels for the legend. Same order as df.
            :type label_list: list[str]
            :param color_list: list of colors for the points and legend. Same order as df.
            :type color_list: list[str]
            :param interactive: True to show trial numbers if one moves the mouse over it
            :type interactive: bool
            """
            if color_list is None:
                color_list = ['red', 'blue', 'green', 'gray']
            for count, df in enumerate(dataframes):
                # color_list was before list(ff.colors_femmt_default.keys())
                df['color_r'], df['color_g'], df['color_b'] = ff.colors_femmt_default[color_list[count]]

            df_all = pd.concat(dataframes, axis=0)
            color_array = np.transpose(np.array([df_all['color_r'].to_numpy(), df_all['color_g'].to_numpy(), df_all['color_b'].to_numpy()])) / 255

            names = df_all["number"].to_numpy()
            fig, ax = plt.subplots()
            legend_list = []
            for count, label_text in enumerate(label_list):
                legend_list.append(mpatches.Patch(color=np.array(ff.colors_femmt_default[color_list[count]]) / 255, label=label_text))
            plt.legend(handles=legend_list)
            sc = plt.scatter(df_all["values_0"], df_all["values_1"], s=10, c=color_array)

            if interactive:
                annot = ax.annotate("", xy=(0, 0), xytext=(20, 20), textcoords="offset points",
                                    bbox=dict(boxstyle="round", fc="w"),
                                    arrowprops=dict(arrowstyle="->"))
                annot.set_visible(False)

                def update_annot(ind):
                    pos = sc.get_offsets()[ind["ind"][0]]
                    annot.xy = pos
                    text = f"{[names[n] for n in ind['ind']]}"
                    annot.set_text(text)
                    annot.get_bbox_patch().set_alpha(0.4)

                def hover(event):
                    vis = annot.get_visible()
                    if event.inaxes == ax:
                        cont, ind = sc.contains(event)
                        if cont:
                            update_annot(ind)
                            annot.set_visible(True)
                            fig.canvas.draw_idle()
                        else:
                            if vis:
                                annot.set_visible(False)
                                fig.canvas.draw_idle()

                fig.canvas.mpl_connect("motion_notify_event", hover)
            plt.xlabel('Volume in m³')
            plt.ylabel('Losses in W')
            plt.grid()
            plt.show()

        @staticmethod
        def filter_loss_list_df(df: pd.DataFrame, factor_min_dc_losses: float = 1.2, factor_max_dc_losses: float = 10) -> pd.DataFrame:
            """
            Remove designs with too high losses compared to the minimum losses.

            :param df: list of valid DTOs
            :type df: list[ItoSingleResultFile]
            :param factor_min_dc_losses: filter factor for the minimum dc losses
            :type factor_min_dc_losses: float
            :param factor_max_dc_losses: dc_max_loss = factor_max_dc_losses * min_available_dc_losses_in_pareto_front
            :type factor_max_dc_losses: float
            :returns: list with removed objects (too small air gaps)
            :rtype: list[ItoSingleResultFile]
            """
            # figure out pareto front
            # pareto_volume_list, pareto_core_hyst_list, pareto_dto_list = self.pareto_front(volume_list, core_hyst_loss_list, valid_design_list)

            pareto_df: pd.DataFrame = fo.pareto_front_from_df(df)

            vector_to_sort = np.array([pareto_df["values_0"], pareto_df["values_1"]])

            # sorting 2d array by 1st row
            # https://stackoverflow.com/questions/49374253/sort-a-numpy-2d-array-by-1st-row-maintaining-columns
            sorted_vector = vector_to_sort[:, vector_to_sort[0].argsort()]
            x_pareto_vec = sorted_vector[0]
            y_pareto_vec = sorted_vector[1]

            total_losses_list = df["values_1"][~np.isnan(df['values_1'])].to_numpy()

            min_total_dc_losses = total_losses_list[np.argmin(total_losses_list)]
            loss_offset = factor_min_dc_losses * min_total_dc_losses

            ref_loss_max = np.interp(df["values_0"], x_pareto_vec, y_pareto_vec) + loss_offset
            # clip losses to a maximum of the minimum losses
            ref_loss_max = np.clip(ref_loss_max, a_min=-1, a_max=factor_max_dc_losses * min_total_dc_losses)

            pareto_df_offset = df[df['values_1'] < ref_loss_max]

            return pareto_df_offset

        @staticmethod
        def df_trial_numbers(df: pd.DataFrame, trial_number_list: list[int]) -> pd.DataFrame:
            """
            Generate a new dataframe from a given one, just with the trial numbers from the trial_number_list.

            :param df: input dataframe
            :type df: pandas.DataFrame
            :param trial_number_list: list of trials, e.g. [1530, 1870, 3402]
            :type trial_number_list: list[int]
            :return: dataframe with trial numbers from trial_number_list
            :rtype: pandas.DataFrame
            """
            df_trial_numbers = df.loc[df["number"].isin(trial_number_list)]

            return df_trial_numbers

        @staticmethod
        def save_config(config: StoSingleInputConfig) -> None:
            """
            Save the configuration file as pickle file on the disk.

            :param config: configuration
            :type config: InductorOptimizationDTO
            """
            # convert config path to an absolute filepath
            config.inductor_optimization_directory = os.path.abspath(config.stacked_transformer_optimization_directory)
            os.makedirs(config.stacked_transformer_optimization_directory, exist_ok=True)
            with open(f"{config.stacked_transformer_optimization_directory}/{config.stacked_transformer_study_name}.pkl", 'wb') as output:
                pickle.dump(config, output, pickle.HIGHEST_PROTOCOL)

        @staticmethod
        def load_config(config_pickle_filepath: str) -> StoSingleInputConfig:
            """
            Load pickle configuration file from disk.

            :param config_pickle_filepath: filepath to the pickle configuration file
            :type config_pickle_filepath: str
            :return: Configuration file as InductorOptimizationDTO object
            :rtype: InductorOptimizationDTO
            """
            with open(config_pickle_filepath, 'rb') as pickle_file_data:
                return pickle.load(pickle_file_data)

        @staticmethod
        def full_simulation(df_geometry: pd.DataFrame, current_waveform: list, stacked_transformer_config_filepath: str):
            """
            Reluctance model (hysteresis losses) and FEM simulation (winding losses and eddy current losses) for geometries from df_geometry.

            :param df_geometry: Pandas dataframe with geometries
            :type df_geometry: pd.DataFrame
            :param current_waveform: Current waveform to simulate
            :type current_waveform: list
            :param stacked_transformer_config_filepath: Filepath of the inductor optimization configuration file
            :type stacked_transformer_config_filepath: str
            """
            for index, _ in df_geometry.iterrows():

                local_config = StackedTransformerOptimization.ReluctanceModel.load_config(stacked_transformer_config_filepath)

                if local_config.core_name_list is not None:
                    # using fixed core sizes from the database with flexible height.
                    core_name = df_geometry['params_core_name'][index]
                    core = ff.core_database()[core_name]
                    core_inner_diameter = core["core_inner_diameter"]
                    window_w = core["window_w"]
                else:
                    core_inner_diameter = df_geometry['params_core_inner_diameter'][index]
                    window_w = df_geometry['params_window_w'][index]

                # overwrite the old time-current vector with the new one
                local_config.time_current_vec = current_waveform
                target_and_fix_parameters = StackedTransformerOptimization.ReluctanceModel.calculate_fix_parameters(local_config)

                litz_wire_primary_dict = ff.litz_database()[df_geometry['params_primary_litz_name'][index]]
                litz_wire_diameter_primary = 2 * litz_wire_primary_dict["conductor_radii"]
                litz_wire_secondary_dict = ff.litz_database()[df_geometry['params_secondary_litz_name'][index]]
                litz_wire_diameter_secondary = 2 * litz_wire_secondary_dict["conductor_radii"]

                # material properties
                material_db = mdb.Data()
                material_name = Material(df_geometry['params_material_name'][index])
                small_signal_mu_real_over_f_at_T = material_db.get_datasheet_curve(material_name, mdb.DatasheetCurveType.small_signal_mu_real_over_f_at_T)
                small_signal_mu_imag_over_f_at_T = material_db.get_datasheet_curve(material_name, mdb.DatasheetCurveType.small_signal_mu_imag_over_f_at_T)

                mu_real_at_f = np.interp(10_000, small_signal_mu_real_over_f_at_T["f"], small_signal_mu_real_over_f_at_T["mu_real"])
                mu_imag_at_f = np.interp(10_000, small_signal_mu_imag_over_f_at_T["f"], small_signal_mu_imag_over_f_at_T["mu_imag"])

                material_mu_r_abs = np.abs(np.array([mu_real_at_f + 1j * mu_imag_at_f]))

                # instantiate material-specific model
                magnet_material_model: mh.loss.LossModel = mh.loss.LossModel(material=material_name, team="paderborn")

                reluctance_model_input = ReluctanceModelInput(
                    target_inductance_matrix=fr.calculate_inductance_matrix_from_ls_lh_n(local_config.l_s12_target, local_config.l_h_target,
                                                                                         local_config.n_target),
                    core_inner_diameter=core_inner_diameter,
                    window_w=window_w,
                    window_h_bot=df_geometry["params_window_h_bot"][index].item(),
                    window_h_top=df_geometry["user_attrs_window_h_top"][index].item(),
                    turns_1_top=int(df_geometry["params_n_p_top"][index].item()),
                    turns_1_bot=int(df_geometry["params_n_p_bot"][index].item()),
                    turns_2_bot=int(df_geometry["params_n_s_bot"][index].item()),
                    litz_wire_name_1=df_geometry['params_primary_litz_name'][index],
                    litz_wire_diameter_1=litz_wire_diameter_primary,
                    litz_wire_name_2=df_geometry['params_secondary_litz_name'][index],
                    litz_wire_diameter_2=litz_wire_diameter_secondary,

                    insulations=local_config.insulations,
                    material_name=material_name,
                    material_mu_r_abs=material_mu_r_abs,
                    magnet_material_model=magnet_material_model,

                    temperature=local_config.temperature,
                    time_extracted_vec=target_and_fix_parameters.time_extracted_vec,
                    current_extracted_vec_1=target_and_fix_parameters.current_extracted_1_vec,
                    current_extracted_vec_2=target_and_fix_parameters.current_extracted_2_vec,
                    fundamental_frequency=target_and_fix_parameters.fundamental_frequency,

                    i_rms_1=target_and_fix_parameters.i_rms_1,
                    i_rms_2=target_and_fix_parameters.i_rms_2,

                    primary_litz_dict=litz_wire_primary_dict,
                    secondary_litz_dict=litz_wire_secondary_dict,

                    # winding 1
                    fft_frequency_list_1=target_and_fix_parameters.fft_frequency_list_1,
                    fft_amplitude_list_1=target_and_fix_parameters.fft_amplitude_list_1,
                    fft_phases_list_1=target_and_fix_parameters.fft_phases_list_1,

                    # winding 2
                    fft_frequency_list_2=target_and_fix_parameters.fft_frequency_list_2,
                    fft_amplitude_list_2=target_and_fix_parameters.fft_amplitude_list_2,
                    fft_phases_list_2=target_and_fix_parameters.fft_phases_list_2
                )

                reluctance_output: ReluctanceModelOutput = StackedTransformerOptimization.ReluctanceModel.single_reluctance_model_simulation(
                    reluctance_model_input)

                p_total = reluctance_output.p_loss

                return reluctance_output.volume, p_total, reluctance_output.area_to_heat_sink

    class FemSimulation:
        """Contains methods to perform FEM simulations or process their results."""

        @staticmethod
        def fem_simulations_from_reluctance_df(reluctance_df: pd.DataFrame, config: StoSingleInputConfig, show_visual_outputs: bool = False,
                                               process_number: int = 1):
            """
            Perform FEM simulations from a given Pandas dataframe. The dataframe is from the reluctance model results.

            :param reluctance_df: Pandas dataframe containing results from the reluctance model
            :type reluctance_df: pandas.DataFrame
            :param config: Configuration for the optimization of the transformer
            :type config: StoSingleInputConfig
            :param show_visual_outputs: True to show visual outputs like the geometry
            :type show_visual_outputs: bool
            :param process_number: Process number for parallel simulations on multiple CPU cores
            :type process_number: int
            """
            target_and_fix_parameters = StackedTransformerOptimization.ReluctanceModel.calculate_fix_parameters(config)

            working_directory = target_and_fix_parameters.working_directories.fem_working_directory
            if not os.path.exists(working_directory):
                os.mkdir(working_directory)
            working_directory = os.path.join(
                target_and_fix_parameters.working_directories.fem_working_directory, f"process_{process_number}")

            # pd.read_csv(current_waveforms_csv_file, header=0, index_col=0, delimiter=';')
            for index, _ in tqdm.tqdm(reluctance_df.iterrows(), total=reluctance_df.shape[0]):

                destination_json_file = os.path.join(
                    target_and_fix_parameters.working_directories.fem_simulation_results_directory,
                    f'case_{index}.json')

                if os.path.exists(destination_json_file):
                    logger.info(f'case_{index}.json already exists. Skip simulation.')
                else:

                    try:
                        if reluctance_df["params_core_name"] is not None:
                            core_inner_diameter = reluctance_df["user_attrs_core_inner_diameter"][index].item()
                            window_w = reluctance_df["user_attrs_window_w"][index].item()
                            window_h_bot = reluctance_df["user_attrs_window_h_bot"][index].item()
                        else:
                            core_inner_diameter = reluctance_df["params_core_inner_diameter"][index].item()
                            window_w = reluctance_df["params_window_w"][index].item()
                            window_h_bot = reluctance_df["params_window_h_bot"][index].item()

                        fem_input = FemInput(
                            simulation_name='xx',
                            working_directory=working_directory,
                            core_inner_diameter=core_inner_diameter,
                            window_w=window_w,
                            window_h_bot=window_h_bot,
                            window_h_top=reluctance_df['user_attrs_window_h_top'][index].item(),
                            material_name=reluctance_df["params_material_name"][index],
                            temperature=config.temperature,
                            material_data_sources=config.material_data_sources,
                            air_gap_length_top=reluctance_df["user_attrs_l_top_air_gap"][index].item(),
                            air_gap_length_bot=reluctance_df["user_attrs_l_bot_air_gap"][index].item(),
                            insulations=config.insulations,
                            primary_litz_wire_name=reluctance_df['params_primary_litz_name'][index],
                            secondary_litz_wire_name=reluctance_df['params_secondary_litz_name'][index],

                            turns_primary_top=reluctance_df['params_n_p_top'][index].item(),
                            turns_primary_bot=reluctance_df['params_n_p_bot'][index].item(),
                            turns_secondary_bot=int(reluctance_df['params_n_s_bot'][index].item()),

                            fundamental_frequency=target_and_fix_parameters.fundamental_frequency,

                            time_current_1_vec=config.time_current_1_vec,
                            time_current_2_vec=config.time_current_2_vec,
                        )

                        fem_output = StackedTransformerOptimization.FemSimulation.single_fem_simulation(fem_input, show_visual_outputs=show_visual_outputs)

                        reluctance_df.loc[index, 'n'] = fem_output.n_conc
                        reluctance_df.loc[index, 'l_s_conc'] = fem_output.l_s_conc
                        reluctance_df.loc[index, 'l_h_conc'] = fem_output.l_h_conc
                        reluctance_df.loc[index, 'p_loss_winding_1'] = fem_output.p_loss_winding_1
                        reluctance_df.loc[index, 'p_loss_winding_2'] = fem_output.p_loss_winding_2
                        reluctance_df.loc[index, 'eddy_core'] = fem_output.eddy_core
                        reluctance_df.loc[index, 'core'] = fem_output.core

                        # copy result files to result-file folder
                        source_json_file = os.path.join(
                            target_and_fix_parameters.working_directories.fem_working_directory, f'process_{process_number}',
                            "results", "log_electro_magnetic.json")

                        shutil.copy(source_json_file, destination_json_file)

                    except Exception as e:
                        print(f"The following exception was raised: {e}")
                        reluctance_df.loc[index, 'n'] = None
                        reluctance_df.loc[index, 'l_s_conc'] = None
                        reluctance_df.loc[index, 'l_h_conc'] = None
                        reluctance_df.loc[index, 'p_loss_winding_1'] = None
                        reluctance_df.loc[index, 'p_loss_winding_2'] = None
                        reluctance_df.loc[index, 'eddy_core'] = None
                        reluctance_df.loc[index, 'core'] = None
            return reluctance_df

        @staticmethod
        def single_fem_simulation(fem_input: FemInput, show_visual_outputs: bool = False, wire_loss_only: bool = False):
            """
            Perform a single FEM simulation.

            :param fem_input: FEM input DTO
            :type fem_input: FemInput
            :param show_visual_outputs: True to show visual outputs
            :type show_visual_outputs: bool
            :param wire_loss_only: True: only wire loss and inductance calculation, False: Full study (wire loss, hyst. losses, inductance calculation)
            :type wire_loss_only: bool
            :return: FEM output DTO
            :rtype: FemOutput
            """
            time_current_vectors = np.array([fem_input.time_current_1_vec, fem_input.time_current_2_vec])

            # 1. chose simulation type
            geo = fmt.MagneticComponent(component_type=fmt.ComponentType.IntegratedTransformer,
                                        working_directory=fem_input.working_directory, onelab_verbosity=fmt.Verbosity.Silent)

            # 2. set core parameters
            core_dimensions = fmt.dtos.StackedCoreDimensions(core_inner_diameter=fem_input.core_inner_diameter,
                                                             window_w=fem_input.window_w,
                                                             window_h_top=fem_input.window_h_top,
                                                             window_h_bot=fem_input.window_h_bot)

            core_material = fmt.ImportedComplexCoreMaterial(material=fem_input.material_name,
                                                            temperature=fem_input.temperature,
                                                            permeability_datasource=fem_input.material_data_sources.permeability_datasource,
                                                            permittivity_datasource=fem_input.material_data_sources.permittivity_datasource)

            core = fmt.Core(material=core_material,
                            core_type=fmt.CoreType.Stacked,
                            core_dimensions=core_dimensions,
                            detailed_core_model=False)

            geo.set_core(core)

            # 3. set air gap parameters
            air_gaps = fmt.AirGaps(fmt.AirGapMethod.Stacked, core)
            air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, fem_input.air_gap_length_top,
                                 stacked_position=fmt.StackedPosition.Top)
            air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, fem_input.air_gap_length_bot,
                                 stacked_position=fmt.StackedPosition.Bot)
            geo.set_air_gaps(air_gaps)

            # 4. set insulations
            insulation = fmt.Insulation(flag_insulation=False)
            insulation.add_top_section_core_insulations(fem_input.insulations.iso_window_top_core_top, fem_input.insulations.iso_window_top_core_bot,
                                                        fem_input.insulations.iso_window_top_core_left, fem_input.insulations.iso_window_top_core_right)
            insulation.add_bottom_section_core_insulations(fem_input.insulations.iso_window_bot_core_top, fem_input.insulations.iso_window_bot_core_bot,
                                                           fem_input.insulations.iso_window_bot_core_left, fem_input.insulations.iso_window_bot_core_right)
            insulation.add_winding_insulations([[fem_input.insulations.iso_primary_to_primary, fem_input.insulations.iso_secondary_to_secondary],
                                                [fem_input.insulations.iso_secondary_to_secondary, fem_input.insulations.iso_primary_to_primary]])
            geo.set_insulation(insulation)

            winding_window_top, winding_window_bot = fmt.create_stacked_winding_windows(core, insulation)

            vww_top = winding_window_top.split_window(fmt.WindingWindowSplit.NoSplit)
            vww_bot = winding_window_bot.split_window(fmt.WindingWindowSplit.NoSplit)

            # 6. set conductor parameters
            winding1 = fmt.Conductor(0, fmt.ConductorMaterial.Copper)
            primary_litz_wire = fmt.litz_database()[fem_input.primary_litz_wire_name]
            winding1.set_litz_round_conductor(primary_litz_wire['conductor_radii'], primary_litz_wire['strands_numbers'],
                                              primary_litz_wire['strand_radii'], None, fmt.ConductorArrangement.Square)

            winding2 = fmt.Conductor(1, fmt.ConductorMaterial.Copper)
            secondary_litz_wire = fmt.litz_database()[fem_input.secondary_litz_wire_name]
            winding2.set_litz_round_conductor(secondary_litz_wire['conductor_radii'], secondary_litz_wire['strands_numbers'],
                                              secondary_litz_wire['strand_radii'], None, fmt.ConductorArrangement.Square)

            primary_coil_turns = fem_input.turns_primary_top

            # 7. add conductor to vww and add winding window to MagneticComponent
            vww_top.set_interleaved_winding(winding1, primary_coil_turns, winding2, 0, fmt.InterleavedWindingScheme.HorizontalAlternating)
            vww_bot.set_interleaved_winding(winding1, fem_input.turns_primary_bot, winding2,
                                            fem_input.turns_secondary_bot, fmt.InterleavedWindingScheme.HorizontalAlternating)

            geo.set_winding_windows([winding_window_top, winding_window_bot])

            geo.create_model(freq=fem_input.fundamental_frequency, pre_visualize_geometry=show_visual_outputs, save_png=False)

            if wire_loss_only:
                hyst_frequency, _, _ = ff.hysteresis_current_excitation(time_current_vectors)
                inductance_dict = geo.get_inductances(I0=1, skin_mesh_factor=1, op_frequency=hyst_frequency, silent=True)

                study_excitation = geo.stacked_core_study_excitation(time_current_vectors, plot_waveforms=False,
                                                                     fft_filter_value_factor=0.1,
                                                                     transfer_ratio_n=inductance_dict["n_conc"])

                geo.excitation_sweep(study_excitation["linear_losses"]["frequencies"],
                                     study_excitation["linear_losses"]["current_amplitudes"],
                                     study_excitation["linear_losses"]["current_phases_deg"],
                                     inductance_dict=inductance_dict)

                result_dict = geo.read_log()

                fem_output = FemOutput(
                    n_conc=result_dict['inductances']['n_conc'],
                    l_s_conc=result_dict['inductances']['l_s_conc'],
                    l_h_conc=result_dict['inductances']['l_h_conc'],
                    p_loss_winding_1=result_dict['total_losses']['winding1']['total'],
                    p_loss_winding_2=result_dict['total_losses']['winding2']['total'],
                    eddy_core=result_dict['total_losses']['eddy_core'],
                    core=0,
                    volume=result_dict["misc"]["core_2daxi_total_volume"],
                )

            else:
                geo.stacked_core_study(number_primary_coil_turns=primary_coil_turns, time_current_vectors=time_current_vectors,
                                       plot_waveforms=show_visual_outputs, fft_filter_value_factor=0.1)

                result_dict = geo.read_log()

                fem_output = FemOutput(
                    n_conc=result_dict['inductances']['n_conc'],
                    l_s_conc=result_dict['inductances']['l_s_conc'],
                    l_h_conc=result_dict['inductances']['l_h_conc'],
                    p_loss_winding_1=result_dict['total_losses']['winding1']['total'],
                    p_loss_winding_2=result_dict['total_losses']['winding2']['total'],
                    eddy_core=result_dict['total_losses']['eddy_core'],
                    core=result_dict['total_losses']['core'],
                    volume=result_dict["misc"]["core_2daxi_total_volume"],
                )

            return fem_output

        @staticmethod
        def fem_logs_to_df(reluctance_df: pd.DataFrame, fem_results_folder_path: str) -> pd.DataFrame:
            """
            Search the given fem_results_folder_path for .json log files and parse them into a pandas dataframe.

            :param reluctance_df: Pandas Dataframe of reluctance model simulation results
            :type reluctance_df: pd.DataFrame
            :param fem_results_folder_path: Folder with FEM simulation results, simulation results must be named like 'case_123.json'
            :type fem_results_folder_path: str
            :return: Pandas dataframe containing reluctance model results and results from FEM simulations
            :rtype: pd.DataFrame
            """
            files_in_folder = [f for f in os.listdir(fem_results_folder_path) if os.path.isfile(os.path.join(fem_results_folder_path, f))]
            logger.info(files_in_folder)
            logger.info(len(files_in_folder))

            reluctance_df_copy = reluctance_df.copy()

            # add new columns to the dataframe, init values with None
            reluctance_df_copy['fem_inductance'] = None
            reluctance_df_copy['fem_p_loss_winding'] = None
            reluctance_df_copy['fem_eddy_core'] = None
            reluctance_df_copy['fem_core'] = None

            for file in files_in_folder:
                with open(os.path.join(fem_results_folder_path, file), 'r') as log_file:
                    scanned_log_dict = json.load(log_file)

                index = int(file.replace("case_", "").replace(".json", ""))
                # only write FEM simulation log results to df in case of reluctance number simulation is already present there
                if index in reluctance_df_copy["number"]:
                    reluctance_df_copy.at[index, 'fem_inductance'] = scanned_log_dict['single_sweeps'][0]['winding1']['flux_over_current'][0]
                    reluctance_df_copy.at[index, 'fem_p_loss_winding'] = (
                        scanned_log_dict['total_losses']['winding1']['total'] + scanned_log_dict['total_losses']['winding2']['total'])
                    reluctance_df_copy.at[index, 'fem_eddy_core'] = scanned_log_dict['total_losses']['eddy_core']
                    reluctance_df_copy.at[index, 'fem_core'] = scanned_log_dict['total_losses']['core']

            # final loss calculation
            reluctance_df_copy["combined_losses"] = (reluctance_df_copy["fem_eddy_core"] + reluctance_df_copy["fem_p_loss_winding"] + \
                                                     reluctance_df_copy["user_attrs_p_hyst"])

            return reluctance_df_copy

        @staticmethod
        def full_simulation(df_geometry: pd.DataFrame, current_waveform: list, stacked_transformer_config_filepath: str, process_number: int = 1,
                            show_visual_outputs: bool = False, print_derivations: bool = False):
            """
            Reluctance model (hysteresis losses) and FEM simulation (winding losses and eddy current losses) for geometries from df_geometry.

            :param df_geometry: Pandas dataframe with geometries
            :type df_geometry: pd.DataFrame
            :param current_waveform: Current waveform to simulate
            :type current_waveform: list
            :param stacked_transformer_config_filepath: Filepath of the inductor optimization configuration file
            :type stacked_transformer_config_filepath: str
            :param process_number: process number to run the simulation on
            :type process_number: int
            :param show_visual_outputs: True to show visual outputs (geometries)
            :type show_visual_outputs: bool
            :param print_derivations: True to print derivation from FEM simulation to reluctance model
            :type print_derivations: bool
            """
            index_number = df_geometry.first_valid_index()

            if df_geometry.shape[0] > 1:
                raise ValueError("df_geometry must have only one entry.")

            local_config = StackedTransformerOptimization.ReluctanceModel.load_config(stacked_transformer_config_filepath)

            if local_config.core_name_list is not None:
                # using fixed core sizes from the database with flexible height.
                core_name = df_geometry['params_core_name'][index_number]
                core = ff.core_database()[core_name]
                core_inner_diameter = core["core_inner_diameter"]
                window_w = core["window_w"]
            else:
                core_inner_diameter = df_geometry['params_core_inner_diameter'][index_number]
                window_w = df_geometry['params_window_w'][index_number]

            # overwrite the old time-current vector with the new one
            local_config.time_current_vec = current_waveform
            target_and_fix_parameters = StackedTransformerOptimization.ReluctanceModel.calculate_fix_parameters(local_config)

            # material properties
            material_db = mdb.Data()

            material_name = Material(df_geometry['params_material_name'][index_number])

            small_signal_mu_real_over_f_at_T = material_db.get_datasheet_curve(material_name, mdb.DatasheetCurveType.small_signal_mu_real_over_f_at_T)
            small_signal_mu_imag_over_f_at_T = material_db.get_datasheet_curve(material_name, mdb.DatasheetCurveType.small_signal_mu_imag_over_f_at_T)

            mu_real_at_f = np.interp(10_000, small_signal_mu_real_over_f_at_T["f"], small_signal_mu_real_over_f_at_T["mu_real"])
            mu_imag_at_f = np.interp(10_000, small_signal_mu_imag_over_f_at_T["f"], small_signal_mu_imag_over_f_at_T["mu_imag"])

            material_mu_r_abs = np.abs(np.array([mu_real_at_f + 1j * mu_imag_at_f]))

            working_directory = target_and_fix_parameters.working_directories.fem_working_directory
            if not os.path.exists(working_directory):
                os.mkdir(working_directory)
            working_directory = os.path.join(
                target_and_fix_parameters.working_directories.fem_working_directory, f"process_{process_number}")

            fem_input = FemInput(
                # general parameters
                working_directory=working_directory,
                simulation_name="xx",

                # material and geometry parameters
                material_name=material_name,
                primary_litz_wire_name=df_geometry['params_primary_litz_name'][index_number],
                secondary_litz_wire_name=df_geometry['params_secondary_litz_name'][index_number],
                core_inner_diameter=core_inner_diameter,
                window_w=window_w,
                window_h_top=df_geometry["user_attrs_window_h_top"][index_number].item(),
                window_h_bot=df_geometry["params_window_h_bot"][index_number].item(),
                air_gap_length_top=df_geometry["user_attrs_l_top_air_gap"][index_number].item(),
                air_gap_length_bot=df_geometry["user_attrs_l_bot_air_gap"][index_number].item(),
                turns_primary_top=int(df_geometry["params_n_p_top"][index_number].item()),
                turns_primary_bot=int(df_geometry["params_n_p_bot"][index_number].item()),
                turns_secondary_bot=int(df_geometry["params_n_s_bot"][index_number].item()),
                insulations=local_config.insulations,

                # data sources
                material_data_sources=local_config.material_data_sources,

                # operating point conditions
                temperature=local_config.temperature,
                fundamental_frequency=target_and_fix_parameters.fundamental_frequency,
                time_current_1_vec=[list(target_and_fix_parameters.time_extracted_vec), list(target_and_fix_parameters.current_extracted_1_vec)],
                time_current_2_vec=[list(target_and_fix_parameters.time_extracted_vec), list(target_and_fix_parameters.current_extracted_2_vec)]
            )

            # pd.read_csv("~/Downloads/Pandas_trial.csv", header=0, index_col=0, delimiter=';')
            fem_output = StackedTransformerOptimization.FemSimulation.single_fem_simulation(fem_input, show_visual_outputs)

            litz_wire_primary_dict = ff.litz_database()[df_geometry['params_primary_litz_name'][index_number]]
            litz_wire_diameter_primary = 2 * litz_wire_primary_dict["conductor_radii"]
            litz_wire_secondary_dict = ff.litz_database()[df_geometry['params_secondary_litz_name'][index_number]]
            litz_wire_diameter_secondary = 2 * litz_wire_secondary_dict["conductor_radii"]

            # instantiate material-specific model
            magnet_material_model: mh.loss.LossModel = mh.loss.LossModel(material=material_name, team="paderborn")

            reluctance_model_input = ReluctanceModelInput(
                target_inductance_matrix=fr.calculate_inductance_matrix_from_ls_lh_n(local_config.l_s12_target, local_config.l_h_target,
                                                                                     local_config.n_target),
                core_inner_diameter=core_inner_diameter,
                window_w=window_w,
                window_h_bot=df_geometry["params_window_h_bot"][index_number].item(),
                window_h_top=df_geometry["user_attrs_window_h_top"][index_number].item(),
                turns_1_top=int(df_geometry["params_n_p_top"][index_number].item()),
                turns_1_bot=int(df_geometry["params_n_p_bot"][index_number].item()),
                turns_2_bot=int(df_geometry["params_n_s_bot"][index_number].item()),
                litz_wire_name_1=df_geometry['params_primary_litz_name'][index_number],
                litz_wire_diameter_1=litz_wire_diameter_primary,
                litz_wire_name_2=df_geometry['params_secondary_litz_name'][index_number],
                litz_wire_diameter_2=litz_wire_diameter_secondary,

                insulations=local_config.insulations,
                material_name=material_name,
                material_mu_r_abs=material_mu_r_abs,
                magnet_material_model=magnet_material_model,

                temperature=local_config.temperature,
                time_extracted_vec=target_and_fix_parameters.time_extracted_vec,
                current_extracted_vec_1=target_and_fix_parameters.current_extracted_1_vec,
                current_extracted_vec_2=target_and_fix_parameters.current_extracted_2_vec,
                fundamental_frequency=target_and_fix_parameters.fundamental_frequency,

                i_rms_1=target_and_fix_parameters.i_rms_1,
                i_rms_2=target_and_fix_parameters.i_rms_2,

                primary_litz_dict=litz_wire_primary_dict,
                secondary_litz_dict=litz_wire_secondary_dict,

                # winding 1
                fft_frequency_list_1=target_and_fix_parameters.fft_frequency_list_1,
                fft_amplitude_list_1=target_and_fix_parameters.fft_amplitude_list_1,
                fft_phases_list_1=target_and_fix_parameters.fft_phases_list_1,

                # winding 2
                fft_frequency_list_2=target_and_fix_parameters.fft_frequency_list_2,
                fft_amplitude_list_2=target_and_fix_parameters.fft_amplitude_list_2,
                fft_phases_list_2=target_and_fix_parameters.fft_phases_list_2
            )

            reluctance_output: ReluctanceModelOutput = StackedTransformerOptimization.ReluctanceModel.single_reluctance_model_simulation(
                reluctance_model_input)

            p_total = reluctance_output.p_hyst + fem_output.eddy_core + fem_output.p_loss_winding_1 + fem_output.p_loss_winding_2

            if print_derivations:
                logger.info(f"Inductance l_h reluctance: {local_config.l_h_target}")
                logger.info(f"Inductance l_h FEM: {fem_output.l_h_conc}")
                logger.info(f"Inductance l_h derivation: {(fem_output.l_h_conc - local_config.l_h_target) / local_config.l_h_target * 100} %")
                logger.info(f"Inductance l_s reluctance: {local_config.l_s12_target}")
                logger.info(f"Inductance l_s FEM: {fem_output.l_s_conc}")
                logger.info(f"Inductance l_s derivation: {(fem_output.l_s_conc - local_config.l_s12_target) / local_config.l_s12_target * 100} %")
                logger.info(f"Volume reluctance: {reluctance_output.volume}")
                logger.info(f"Volume FEM: {fem_output.volume}")
                logger.info(f"Volume derivation: {(reluctance_output.volume - fem_output.volume) / reluctance_output.volume * 100} %")

                logger.info(f"P_winding_both reluctance: {reluctance_output.winding_1_loss + reluctance_output.winding_2_loss}")
                logger.info(f"P_winding_both FEM: {fem_output.p_loss_winding_1 + fem_output.p_loss_winding_2}")
                winding_derivation = ((fem_output.p_loss_winding_1 + fem_output.p_loss_winding_2 - \
                                      (reluctance_output.winding_1_loss + reluctance_output.winding_2_loss)) / \
                                      (fem_output.p_loss_winding_1 + fem_output.p_loss_winding_2) * 100)
                logger.info(f"P_winding_both derivation: "
                            f"{winding_derivation}")
                logger.info(f"P_hyst reluctance: {reluctance_output.p_hyst}")
                logger.info(f"P_hyst FEM: {fem_output.core}")
                logger.info(f"P_hyst derivation: {(reluctance_output.p_hyst - fem_output.core) / reluctance_output.p_hyst * 100}")

            return reluctance_output.volume, p_total, reluctance_output.area_to_heat_sink
