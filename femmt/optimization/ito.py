"""Integrated transformer optimization."""
# Python libraries
import os
import json
from typing import List, Dict, Optional
import datetime
import pickle
import logging

# 3rd party library import
import materialdatabase as mdb
from scipy import optimize
import optuna
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import magnethub as mh

# femmt import
import femmt.functions as ff
import femmt.functions_reluctance as fr
import femmt.optimization.functions_optimization as fo
from femmt.optimization.ito_dtos import *
import femmt.optimization.ito_functions as itof
import femmt

class MyJSONEncoder(json.JSONEncoder):
    """
    Class to transform dicts with numpy arrays to json.

    This class is used as cls=MyJSONEncoder by json.dump

    See Also
    --------
    https://python-forum.io/thread-35245.html
    """

    def default(self, o: Dict):
        """
        Transform the dictionary to a .json file.

        :param o: Dictionary to transform
        :type o: Dict
        """
        try:
            return o.tolist()  # works with any object that has .tolist() method
        except AttributeError:
            pass
        # Let the base class default method raise the TypeError
        return json.JSONEncoder.default(self, o)


def result_file_dict_to_dto(result_file_dict: Dict) -> ItoSingleResultFile:
    """
    Translate the result file dictionary to a data transfer object (DTO).

    :param result_file_dict: dictionary to translate to the DTO structure
    :type result_file_dict: Dict
    :return: DTO
    :rtype: ItoSingleResultFile
    """
    result_file_dto = ItoSingleResultFile(
        case=result_file_dict["case"],
        air_gap_top=result_file_dict["air_gap_top"],
        air_gap_bot=result_file_dict["air_gap_bot"],
        air_gap_middle=result_file_dict["air_gap_middle"],
        n_p_top=int(result_file_dict["n_p_top"]),
        n_p_bot=int(result_file_dict["n_p_bot"]),
        n_s_top=int(result_file_dict["n_s_top"]),
        n_s_bot=int(result_file_dict["n_s_bot"]),
        window_h_top=result_file_dict["window_h_top"],
        window_h_bot=result_file_dict["window_h_bot"],
        window_w=result_file_dict["window_w"],
        core_material=result_file_dict["core_material"],
        core_inner_diameter=result_file_dict["core_inner_diameter"],
        flux_top_max=result_file_dict["flux_top_max"],
        flux_bot_max=result_file_dict["flux_bot_max"],
        flux_stray_max=result_file_dict["flux_stray_max"],
        flux_density_top_max=result_file_dict["flux_density_top_max"],
        flux_density_bot_max=result_file_dict["flux_density_bot_max"],
        flux_density_stray_max=result_file_dict["flux_density_stray_max"],
        p_hyst=result_file_dict["p_hyst"],
        core_2daxi_total_volume=result_file_dict["core_2daxi_total_volume"],
        litz_wire_name_1=result_file_dict["primary_litz_wire"],
        litz_wire_name_2=result_file_dict["secondary_litz_wire"],
        litz_wire_loss_1=result_file_dict["primary_litz_wire_loss"],
        litz_wire_loss_2=result_file_dict["secondary_litz_wire_loss"],
        total_loss=result_file_dict["total_loss"]
    )
    return result_file_dto


class IntegratedTransformerOptimization:
    """Perform different optimization methods for the integrated transformer."""

    @staticmethod
    def plot(valid_design_list: List[ItoSingleResultFile]) -> None:
        """
        Plot the pareto diagram out of the reluctance model calculation.

        :param valid_design_list:
        :type valid_design_list: List[ItoSingleResultFile]
        """
        volume_list = []
        core_hyst_loss_list = []
        annotation_list = []

        for result in valid_design_list:
            volume_list.append(result.core_2daxi_total_volume)
            core_hyst_loss_list.append(result.total_loss)
            annotation_list.append(result.case)

        fo.plot_2d(volume_list, core_hyst_loss_list, "Volume in m³", "Losses in W", "Pareto Diagram",
                   plot_color="red", annotations=annotation_list)

    # @staticmethod
    # def plot_filtered_pareto_result_list(filter_volume_list, filter_core_hyst_loss_list):
    #     fo.plot_2d(filter_volume_list, filter_core_hyst_loss_list, "volume in m³", "core hysteresis losses in W", "Pareto Diagram", plot_color="red")

    # Very slow for many data points.  Fastest for many costs, most readable
    @staticmethod
    def is_pareto_efficient_dumb(costs: np.array) -> np.array:
        """
        Find the pareto-efficient points.

        :param costs: An (n_points, n_costs) array
        :type costs: np.array
        :return: A (n_points, ) boolean array, indicating whether each point is Pareto efficient
        """
        is_efficient = np.ones(costs.shape[0], dtype=bool)
        for i, c in enumerate(costs):
            is_efficient[i] = np.all(np.any(costs[:i] > c, axis=1)) and np.all(np.any(costs[i + 1:] > c, axis=1))
        return is_efficient

    @staticmethod
    def calculate_fix_parameters(config: ItoSingleInputConfig) -> ItoTargetAndFixedParameters:
        """Calculate fix parameters what can be derived from the input configuration.

        return values are:

            i_rms_1
            i_rms_2
            time_extracted_vec
            current_extracted_1_vec
            current_extracted_2_vec
            material_dto_curve_list
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
        time_extracted, current_extracted_1_vec = fr.time_vec_current_vec_from_time_current_vec(config.time_current_1_vec)
        time_extracted, current_extracted_2_vec = fr.time_vec_current_vec_from_time_current_vec(config.time_current_2_vec)
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
        material_db = mdb.MaterialDatabase(is_silent=True)

        material_data_list = []
        magnet_model_list = []
        for material_name in config.material_list:
            material_dto: mdb.MaterialCurve = material_db.material_data_interpolation_to_dto(material_name, fundamental_frequency, config.temperature)
            material_data_list.append(material_dto)
            # instantiate material-specific model
            mdl: mh.loss.LossModel = mh.loss.LossModel(material=material_name, team="paderborn")
            magnet_model_list.append(mdl)

        # set up working directories
        working_directories = itof.set_up_folder_structure(config.integrated_transformer_optimization_directory)

        # finalize data to dto
        target_and_fix_parameters = ItoTargetAndFixedParameters(
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
            material_dto_curve_list=material_data_list,
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

    class ReluctanceModel:
        """Create and calculate the reluctance model for the integrated transformer."""

        @staticmethod
        def single_reluctance_model_simulation(reluctance_input: ItoReluctanceModelInput) -> ItoReluctanceModelOutput:
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
            winding_matrix = np.array([[reluctance_input.turns_1_top, reluctance_input.turns_2_top],
                                       [reluctance_input.turns_1_bot, reluctance_input.turns_2_bot]])

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
                                                      reluctance_input.material_dto.material_mu_r_abs)
            core_inner_cylinder_bot = fr.r_core_round(reluctance_input.core_inner_diameter, reluctance_input.window_h_bot,
                                                      reluctance_input.material_dto.material_mu_r_abs)
            core_top_bot_radiant = fr.r_core_top_bot_radiant(reluctance_input.core_inner_diameter, reluctance_input.window_w,
                                                             reluctance_input.material_dto.material_mu_r_abs, reluctance_input.core_inner_diameter / 4)

            r_core_top = 2 * core_inner_cylinder_top + core_top_bot_radiant
            r_core_bot = 2 * core_inner_cylinder_bot + core_top_bot_radiant
            r_core_middle = core_top_bot_radiant

            r_top_target = reluctance_matrix[0][0] + reluctance_matrix[0][1]
            r_bot_target = reluctance_matrix[1][1] + reluctance_matrix[0][1]
            r_middle_target = - reluctance_matrix[0][1]

            r_air_gap_top_target = r_top_target - r_core_top
            r_air_gap_bot_target = r_bot_target - r_core_bot
            r_air_gap_middle_target = r_middle_target - r_core_middle

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

            # Note: ideal calculation (360 degree)
            # needs to be translated when it comes to the real setup.
            l_middle_air_gap = optimize.brentq(
                fr.r_air_gap_tablet_cylinder_sct, minimum_air_gap_length, maximum_air_gap_length,
                args=(reluctance_input.core_inner_diameter, reluctance_input.core_inner_diameter/4, reluctance_input.window_w, r_air_gap_middle_target),
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
                    reluctance_input.litz_dict_1["strands_numbers"] * reluctance_input.litz_dict_1["strand_radii"] ** 2 * np.pi)
            primary_effective_conductive_radius = np.sqrt(primary_effective_conductive_cross_section / np.pi)
            primary_resistance_top = fr.resistance_solid_wire(
                reluctance_input.core_inner_diameter, reluctance_input.window_w, reluctance_input.turns_1_top,
                primary_effective_conductive_radius, material='Copper')

            number_bot_prim_turns_per_column = (
                int((reluctance_input.window_h_bot - reluctance_input.insulations.iso_window_bot_core_top - \
                     reluctance_input.insulations.iso_window_bot_core_bot + reluctance_input.insulations.iso_primary_to_primary) / \
                    (2 * reluctance_input.litz_dict_1["conductor_radii"] + reluctance_input.insulations.iso_primary_to_primary)))
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
                    reluctance_input.litz_dict_2["strands_numbers"] * reluctance_input.litz_dict_2["strand_radii"] ** 2 * np.pi)
            secondary_effective_conductive_radius = np.sqrt(secondary_effective_conductive_cross_section / np.pi)
            secondary_resistance = fr.resistance_solid_wire(
                reluctance_input.core_inner_diameter, reluctance_input.window_w, reluctance_input.turns_2_bot, secondary_effective_conductive_radius,
                material='Copper')

            winding_area_1_top = (
                    2 * reluctance_input.litz_dict_1["conductor_radii"] * \
                    (reluctance_input.turns_1_top * 2 * reluctance_input.litz_dict_1["conductor_radii"] + \
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
                    winding_area_1_bot = 2 * reluctance_input.litz_dict_1["conductor_radii"] * \
                                         (reluctance_input.turns_1_bot * 2 * reluctance_input.litz_dict_1["conductor_radii"] + \
                                          (reluctance_input.turns_1_bot - 1) * reluctance_input.insulations.iso_primary_to_primary)

                    proximity_factor_1_bot_inner = fr.calc_proximity_factor_air_gap(
                        litz_wire_name=reluctance_input.litz_wire_name_1, number_turns=reluctance_input.turns_1_bot,
                        r_1=reluctance_input.insulations.iso_window_bot_core_left,
                        frequency=fft_frequency, winding_area=winding_area_1_bot,
                        litz_wire_material_name='Copper', temperature=reluctance_input.temperature)

                    proximity_factor_1_bot_outer = 0
                else:
                    winding_area_1_bot = (2 * reluctance_input.litz_dict_1["conductor_radii"] * (
                            number_bot_prim_turns_per_column * 2 * reluctance_input.litz_dict_1["conductor_radii"] + \
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

            reluctance_output = ItoReluctanceModelOutput(
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
                l_middle_air_gap=l_middle_air_gap,
                volume=volume,
                area_to_heat_sink=area_to_heat_sink,
                p_loss=p_loss_total,
            )

            return reluctance_output

        @staticmethod
        def objective(trial: optuna.Trial, config: ItoSingleInputConfig, target_and_fixed_parameters: ItoTargetAndFixedParameters):
            """
            Objective funktion to optimize. Uses reluctance model calculation.

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
                window_h_half_max = (core["window_h"] - core_inner_diameter) / 2

                trial.set_user_attr('core_inner_diameter', core_inner_diameter)
                trial.set_user_attr('window_w', window_w)

            else:
                # using arbitrary core sizes
                core_inner_diameter = trial.suggest_float("core_inner_diameter", config.core_inner_diameter_min_max_list[0],
                                                          config.core_inner_diameter_min_max_list[1])
                window_w = trial.suggest_float("window_w", config.window_w_min_max_list[0], config.window_w_min_max_list[1])
                window_h_bot = trial.suggest_float('window_h_bot', config.window_h_bot_min_max_list[0], config.window_h_bot_min_max_list[1])

            # suggest turn counts
            n_1_top = trial.suggest_int('n_p_top', config.n_1_top_min_max_list[0], config.n_1_top_min_max_list[1])
            n_1_bot = trial.suggest_int('n_p_bot', config.n_1_bot_min_max_list[0], config.n_1_bot_min_max_list[1])
            n_2_top = trial.suggest_int('n_s_top', config.n_2_top_min_max_list[0], config.n_2_top_min_max_list[1])
            n_2_bot = trial.suggest_int('n_s_bot', config.n_2_bot_min_max_list[0], config.n_2_bot_min_max_list[1])

            # suggest litz wire 1
            litz_name_1 = trial.suggest_categorical('litz_name_1', config.litz_wire_list_1)
            litz_dict_1 = ff.litz_database()[litz_name_1]
            litz_diameter_1 = 2 * litz_dict_1['conductor_radii']

            # suggest litz wire 2
            litz_name_2 = trial.suggest_categorical('litz_name_2', config.litz_wire_list_2)
            litz_dict_2 = ff.litz_database()[litz_name_2]
            litz_diameter_2 = 2 * litz_dict_2['conductor_radii']

            # suggest core material
            material_name = trial.suggest_categorical('material_name', config.material_list)
            for count, material_dto in enumerate(target_and_fixed_parameters.material_dto_curve_list):
                if material_dto.material_name == material_name:
                    material_dto: mdb.MaterialCurve = material_dto
                    magnet_material_model = target_and_fixed_parameters.magnet_hub_model_list[count]

            # calculation of window_h_top: assumption: winding ist stacked
            available_width_top = window_w - config.insulations.iso_window_top_core_left - config.insulations.iso_window_top_core_right
            available_width_bot = window_w - config.insulations.iso_window_bot_core_left - config.insulations.iso_window_bot_core_right

            def calc_winding_height(available_width: float, iso_winding_to_winding: float, litz_diameter: float,
                                    number_turns: float):


                number_turns_per_row = int((available_width + iso_winding_to_winding) / (litz_diameter + iso_winding_to_winding))
                if number_turns_per_row < 1:
                    return float('nan'), float('nan')
                number_of_rows = np.ceil(number_turns / number_turns_per_row)
                winding_height = number_of_rows * litz_diameter + (number_of_rows -1) * iso_winding_to_winding
                return winding_height

            # height of top core needed
            window_h_1_top = calc_winding_height(available_width_top, config.insulations.iso_primary_to_primary, litz_diameter_1, n_1_top)
            window_h_2_top = calc_winding_height(available_width_top, config.insulations.iso_secondary_to_secondary, litz_diameter_2, n_2_top)
            window_h_top = (window_h_1_top + window_h_2_top + config.insulations.iso_primary_to_secondary +
                            config.insulations.iso_window_top_core_top + config.insulations.iso_window_top_core_bot)
            if window_h_top > window_h_half_max:
                return float('nan'), float('nan')

            # height of bot core needed
            window_h_1_bot = calc_winding_height(available_width_bot, config.insulations.iso_primary_to_primary, litz_diameter_1, n_1_bot)
            window_h_2_bot = calc_winding_height(available_width_bot, config.insulations.iso_secondary_to_secondary, litz_diameter_2, n_2_bot)
            window_h_bot = (window_h_1_bot + window_h_2_bot + config.insulations.iso_primary_to_secondary +
                            config.insulations.iso_window_bot_core_top + config.insulations.iso_window_bot_core_bot)
            if window_h_bot > window_h_half_max:
                return float('nan'), float('nan')

            reluctance_model_intput = ItoReluctanceModelInput(
                target_inductance_matrix=target_and_fixed_parameters.target_inductance_matrix,
                core_inner_diameter=core_inner_diameter,
                window_w=window_w,
                window_h_bot=window_h_bot,
                window_h_top=window_h_top,
                turns_1_top=n_1_top,
                turns_1_bot=n_1_bot,
                turns_2_top=n_2_top,
                turns_2_bot=n_2_bot,
                litz_wire_name_1=litz_name_1,
                litz_wire_diameter_1=litz_diameter_1,
                litz_wire_name_2=litz_name_2,
                litz_wire_diameter_2=litz_diameter_2,

                insulations=config.insulations,
                material_dto=material_dto,
                magnet_material_model=magnet_material_model,

                temperature=config.temperature,
                time_extracted_vec=target_and_fixed_parameters.time_extracted_vec,
                current_extracted_vec_1=target_and_fixed_parameters.current_extracted_1_vec,
                current_extracted_vec_2=target_and_fixed_parameters.current_extracted_2_vec,
                fundamental_frequency=target_and_fixed_parameters.fundamental_frequency,
                i_rms_1=target_and_fixed_parameters.i_rms_1,
                i_rms_2=target_and_fixed_parameters.i_rms_2,
                litz_dict_1=litz_dict_1,
                litz_dict_2=litz_dict_2,

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
                reluctance_output = IntegratedTransformerOptimization.ReluctanceModel.single_reluctance_model_simulation(reluctance_model_intput)
            except ValueError:
                logging.warning("No fitting air gap length")
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

        #########################################################
        # set dynamic wire count parameters as optimization parameters
        #########################################################
        # set the winding search space dynamic
        # https://optuna.readthedocs.io/en/stable/faq.html#what-happens-when-i-dynamically-alter-a-search-space

        # if np.linalg.det(t2_reluctance_matrix) != 0 and np.linalg.det(
        #         np.transpose(t2_winding_matrix)) != 0 and np.linalg.det(target_inductance_matrix) != 0:
        #     # calculate the flux
        #     flux_top_vec, flux_bot_vec, flux_stray_vec = fr.flux_vec_from_current_vec(
        #         target_and_fixed_parameters.current_extracted_1_vec,
        #         target_and_fixed_parameters.current_extracted_2_vec,
        #         t2_winding_matrix,
        #         target_inductance_matrix)

        @staticmethod
        def start_proceed_study(config: ItoSingleInputConfig, number_trials: Optional[int] = None,
                                target_number_trials: Optional[int] = None, storage: str = 'sqlite',
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
            if os.path.exists(f"{config.integrated_transformer_optimization_directory}/{config.integrated_transformer_study_name}.sqlite3"):
                print("Existing study found. Proceeding.")

            target_and_fixed_parameters = IntegratedTransformerOptimization.calculate_fix_parameters(config)

            # introduce study in storage, e.g. sqlite or mysql
            if storage == 'sqlite':
                # Note: for sqlite operation, there needs to be three slashes '///' even before the path '/home/...'
                # Means, in total there are four slashes including the path itself '////home/.../database.sqlite3'
                storage = f"sqlite:///{config.integrated_transformer_optimization_directory}/{config.integrated_transformer_study_name}.sqlite3"
            elif storage == 'mysql':
                storage = "mysql://monty@localhost/mydb",

            # set logging verbosity: https://optuna.readthedocs.io/en/stable/reference/generated/optuna.logging.set_verbosity.html#optuna.logging.set_verbosity
            # .INFO: all messages (default)
            # .WARNING: fails and warnings
            # .ERROR: only errors
            optuna.logging.set_verbosity(optuna.logging.ERROR)

            func = lambda trial: IntegratedTransformerOptimization.ReluctanceModel.objective(trial, config, target_and_fixed_parameters)

            study_in_storage = optuna.create_study(study_name=config.integrated_transformer_study_name,
                                                   storage=storage,
                                                   directions=['minimize', 'minimize'],
                                                   load_if_exists=True, sampler=sampler)

            if target_number_trials is not None:
                # simulation for a given number of target trials
                if len(study_in_storage.trials) < target_number_trials:
                    study_in_memory = optuna.create_study(directions=['minimize', 'minimize'],
                                                          study_name=config.integrated_transformer_study_name, sampler=sampler)
                    print(f"Sampler is {study_in_memory.sampler.__class__.__name__}")
                    study_in_memory.add_trials(study_in_storage.trials)
                    number_trials = target_number_trials - len(study_in_memory.trials)
                    study_in_memory.optimize(func, n_trials=number_trials, show_progress_bar=True)
                    study_in_storage.add_trials(study_in_memory.trials[-number_trials:])
                    print(f"Finished {number_trials} trials.")
                    print(f"current time: {datetime.datetime.now()}")
                else:
                    print(f"Study has already {len(study_in_storage.trials)} trials, and target is {target_number_trials} trials.")

            else:
                # normal simulation with number_trials
                study_in_memory = optuna.create_study(directions=['minimize', 'minimize'], study_name=config.integrated_transformer_study_name, sampler=sampler)
                print(f"Sampler is {study_in_memory.sampler.__class__.__name__}")
                study_in_memory.add_trials(study_in_storage.trials)
                study_in_memory.optimize(func, n_trials=number_trials, show_progress_bar=True)

                study_in_storage.add_trials(study_in_memory.trials[-number_trials:])
                print(f"Finished {number_trials} trials.")
                print(f"current time: {datetime.datetime.now()}")
            IntegratedTransformerOptimization.ReluctanceModel.save_config(config)

        @staticmethod
        def study_to_df(config: ItoSingleInputConfig) -> pd.DataFrame:
            """
            Create a Pandas dataframe from a study.

            :param config: configuration
            :type config: InductorOptimizationDTO
            :return: Study results as Pandas Dataframe
            :rtype: pd.DataFrame
            """
            database_url = f'sqlite:///{config.integrated_transformer_optimization_directory}/{config.integrated_transformer_study_name}.sqlite3'
            if os.path.isfile(database_url.replace('sqlite:///', '')):
                print("Existing study found.")
            else:
                raise ValueError(f"Can not find database: {database_url}")
            loaded_study = optuna.load_study(study_name=config.integrated_transformer_study_name, storage=database_url)
            df = loaded_study.trials_dataframe()
            df.to_csv(f'{config.integrated_transformer_optimization_directory}/{config.integrated_transformer_study_name}.csv')
            logging.info(f"Exported study as .csv file: {config.integrated_transformer_optimization_directory}/{config.integrated_transformer_study_name}.csv")
            return df

        #############################
        # filter
        #############################

        @staticmethod
        def filter_loss_list_df(df: pd.DataFrame, factor_min_dc_losses: float = 1.2, factor_max_dc_losses: float = 10) -> pd.DataFrame:
            """
            Remove designs with too high losses compared to the minimum losses.

            :param df: list of valid DTOs
            :type df: List[ItoSingleResultFile]
            :param factor_min_dc_losses: filter factor for the minimum dc losses
            :type factor_min_dc_losses: float
            :param factor_max_dc_losses: dc_max_loss = factor_max_dc_losses * min_available_dc_losses_in_pareto_front
            :type factor_max_dc_losses: float
            :returns: list with removed objects (too small air gaps)
            :rtype: List[ItoSingleResultFile]
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
        def filter_max_air_gap_length(dto_list_to_filter: List[ItoSingleResultFile], max_air_gap_length: float = 1e-6) -> List[ItoSingleResultFile]:
            """
            Remove designs with a too large air gap.

            :param dto_list_to_filter: list of DTOs to filter for too small air gaps
            :type dto_list_to_filter: List[ItoSingleResultFile]
            :param max_air_gap_length: minimum air gap length
            :type max_air_gap_length: float
            :returns: list with removed objects (too small air gaps)
            :rtype: List[ItoSingleResultFile]
            """
            filtered_dtos = []
            for dto in dto_list_to_filter:
                if dto.air_gap_middle < max_air_gap_length and dto.air_gap_top < max_air_gap_length and dto.air_gap_bot < max_air_gap_length:
                    filtered_dtos.append(dto)
            return filtered_dtos

        @staticmethod
        def filter_min_air_gap_length(dto_list_to_filter: List[ItoSingleResultFile], min_air_gap_length: float = 1e-6) -> List[ItoSingleResultFile]:
            """
            Remove designs with a too small air gap.

            :param dto_list_to_filter: list of DTOs to filter for too small air gaps
            :type dto_list_to_filter: List[ItoSingleResultFile]
            :param min_air_gap_length: minimum air gap length
            :type min_air_gap_length: float
            :returns: list with removed objects (too small air gaps)
            :rtype: List[ItoSingleResultFile]
            """
            filtered_dtos = []
            for dto in dto_list_to_filter:
                if dto.air_gap_middle > min_air_gap_length and dto.air_gap_top > min_air_gap_length and dto.air_gap_bot > min_air_gap_length:
                    filtered_dtos.append(dto)
            return filtered_dtos

        #############################
        # show results
        #############################

        @staticmethod
        def show_study_results(config: ItoSingleInputConfig) -> None:
            """Show the results of a study.

            A local .html file is generated under config.working_directory to store the interactive plotly plots on disk.

            :param config: Integrated transformer configuration file
            :type config: ItoSingleInputConfig
            """
            study = optuna.load_study(study_name=config.integrated_transformer_study_name,
                                      storage=f"sqlite:///{config.integrated_transformer_optimization_directory}/"
                                              f"{config.integrated_transformer_study_name}.sqlite3")

            fig = optuna.visualization.plot_pareto_front(study, targets=lambda t: (t.values[0], t.values[1]), target_names=["volume in m³", "loss in W"])
            fig.update_layout(title=f"{config.integrated_transformer_study_name} <br><sup>{config.integrated_transformer_optimization_directory}</sup>")
            fig.write_html(f"{config.integrated_transformer_optimization_directory}/{config.integrated_transformer_study_name}"
                           f"_{datetime.datetime.now().isoformat(timespec='minutes')}.html")
            fig.show()

        @staticmethod
        def df_plot_pareto_front(*dataframes: pd.DataFrame, label_list: list[str], color_list: list[str] = None,
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
                color_list = ['red', 'blue', 'green', 'grey']
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

        #############################
        # save and load
        #############################

        @staticmethod
        def save_config(config: ItoSingleInputConfig) -> None:
            """
            Save the configuration file as pickle file on the disk.

            :param config: configuration
            :type config: InductorOptimizationDTO
            """
            # convert config path to an absolute filepath
            config.inductor_optimization_directory = os.path.abspath(config.integrated_transformer_optimization_directory)
            os.makedirs(config.integrated_transformer_optimization_directory, exist_ok=True)
            with open(f"{config.integrated_transformer_optimization_directory}/{config.integrated_transformer_study_name}.pkl", 'wb') as output:
                pickle.dump(config, output, pickle.HIGHEST_PROTOCOL)

        @staticmethod
        def load_config(config_pickle_filepath: str) -> ItoSingleInputConfig:
            """
            Load pickle configuration file from disk.

            :param config_pickle_filepath: filepath to the pickle configuration file
            :type config_pickle_filepath: str
            :return: Configuration file as InductorOptimizationDTO object
            :rtype: InductorOptimizationDTO
            """
            with open(config_pickle_filepath, 'rb') as pickle_file_data:
                return pickle.load(pickle_file_data)

    class FemSimulation:
        """Group functions to perform FEM simulations."""

        @staticmethod
        def simulate(config_dto: ItoSingleInputConfig, simulation_dto_list: List[ItoSingleResultFile], visualize: bool = False):
            """
            Perform the FEM simulation.

            :param config_dto: Configuration DTO (data transfer object)
            :type config_dto: ItoSingleInputConfig
            :param simulation_dto_list: List of DTOs to simulate
            :type simulation_dto_list: List[ItoSingleResultFile]
            :param visualize: True to visualize the results.
            :type visualize: bool
            """
            femmt.integrated_transformer_fem_simulations_from_result_dtos(config_dto, simulation_dto_list, visualize)

        @staticmethod
        def filter_loss_list(fem_simulations_dict_list: List[Dict], factor_min_dc_losses: float = 0.5) -> List[Dict]:
            """
            Remove too high losses from the given dictionary.

            :param fem_simulations_dict_list: List of dictionaries to filter
            :type fem_simulations_dict_list: List[Dict]
            :param factor_min_dc_losses: factor of the minimum dc losses to filter
            :type factor_min_dc_losses: float
            :return: filtered dictionary list
            :rtype: List[Dict]
            """
            # figure out pareto front
            # pareto_volume_list, pareto_core_hyst_list, pareto_dto_list = self.pareto_front(volume_list, core_hyst_loss_list, valid_design_list)

            x_pareto_vec, y_pareto_vec = fo.pareto_front_from_result_dicts(fem_simulations_dict_list)

            vector_to_sort = np.array([x_pareto_vec, y_pareto_vec])

            # sorting 2d array by 1st row
            # https://stackoverflow.com/questions/49374253/sort-a-numpy-2d-array-by-1st-row-maintaining-columns
            sorted_vector = vector_to_sort[:, vector_to_sort[0].argsort()]
            x_pareto_vec = sorted_vector[0]
            y_pareto_vec = sorted_vector[1]

            total_losses_list = []
            filtered_design_dto_list = []

            for result_dict in fem_simulations_dict_list:
                total_loss = result_dict["total_losses"]["all_windings"] + result_dict["total_losses"]["core"]
                total_losses_list.append(total_loss)

            min_total_dc_losses = total_losses_list[np.argmin(total_losses_list)]
            loss_offset = factor_min_dc_losses * min_total_dc_losses

            for result_dict in fem_simulations_dict_list:
                ref_loss = np.interp(result_dict["misc"]["core_2daxi_total_volume"], x_pareto_vec, y_pareto_vec) + loss_offset
                total_loss = result_dict["total_losses"]["all_windings"] + result_dict["total_losses"]["core"]
                if total_loss < ref_loss:
                    filtered_design_dto_list.append(result_dict)

            return filtered_design_dto_list

        @staticmethod
        def plot(result_log_dict_list: List[Dict]) -> None:
            """
            Plot the pareto diagram out of the fem simulation results.

            :param result_log_dict_list: list of result_log dicts
            :type result_log_dict_list: str
            """
            data_dict_list = []

            volume_list = []
            total_loss_list = []
            total_cost_list = []
            annotation_list = []

            for result_log_dict in result_log_dict_list:

                data_dict_list.append(result_log_dict)
                volume_list.append(result_log_dict["misc"]["core_2daxi_total_volume"])
                total_loss_list.append(result_log_dict["total_losses"]["all_windings"] + result_log_dict["total_losses"]["core"])
                total_cost_list.append(result_log_dict["misc"]["total_cost_incl_margin"])
                annotation_list.append(result_log_dict["case"])

            fo.plot_2d(volume_list, total_loss_list, "Volume in m³", "Losses in W", "Pareto Diagram",
                       plot_color="red", annotations=annotation_list)

        #############################
        # save and load
        #############################

        @staticmethod
        def load(filepath: str) -> List[Dict]:
            """
            Load FEM simulations results from a directory.

            This function appends the case number to the dictionary-list.

            :param filepath: filepath to FEM simulations
            :type: str
            :return: List of result-log dictionaries
            :rtype: List[Dict]
            """
            data_dict_list = []

            for file in os.listdir(filepath):
                if file.endswith(".json"):
                    json_file_path = os.path.join(filepath, file)
                    with open(json_file_path, "r") as fd:
                        loaded_data_dict = json.loads(fd.read())

                    case_number = file.replace("case_", "").replace(".json", "")

                    loaded_data_dict["case"] = f"{case_number}"

                    data_dict_list.append(loaded_data_dict)
            return data_dict_list

        @staticmethod
        def load_unfiltered_results(working_directory: str) -> List[Dict]:
            """
            Load the FEM simulations results as a dict into a list.

            :param working_directory: file path
            :type working_directory: str
            :return:
            :rtype: List[Dict]
            """
            filepath = os.path.join(working_directory, "02_fem_simulation_results")
            return femmt.IntegratedTransformerOptimization.FemSimulation.load(filepath)

        @staticmethod
        def load_filtered_results(working_directory: str) -> List[Dict]:
            """
            Load the FEM simulation results as a dict into a list.

            :param working_directory: file path
            :type working_directory: str
            :return:
            :rtype: List[Dict]
            """
            filepath = os.path.join(working_directory, "02_fem_simulation_results_filtered")
            return femmt.IntegratedTransformerOptimization.FemSimulation.load(filepath)

        @staticmethod
        def save_filtered_results(filtered_dict_list: List[Dict], working_directory: str):
            """
            Save filtered FEM simulation results to hard disk.

            :param filtered_dict_list: list of dictionaries to store
            :type filtered_dict_list: List[Dict]
            :param working_directory: working directory of the optimization
            :param working_directory: str
            """
            for _, result_log in enumerate(filtered_dict_list):
                json_filepath = os.path.join(working_directory, "02_fem_simulation_results_filtered", f"case_{result_log['case']}.json")

                with open(json_filepath, "w") as outfile:
                    json.dump(result_log, outfile)

    class ThermalSimulation:
        """Perform thermal simulations."""

        @staticmethod
        def simulation(config_dto: ItoSingleInputConfig, result_log_dict_list: List[Dict], visualize: bool = False):
            """
            Perform a thermal simulation.

            :param config_dto: Configuration DTO (data transfer object)
            :type config_dto: ItoSingleInputConfig
            :param result_log_dict_list: list of result log dictionaries
            :type result_log_dict_list: List[Dict]
            :param visualize: True to visualize the results.
            :type visualize: bool
            """
            all_filtered_reluctance_dtos = femmt.IntegratedTransformerOptimization.ReluctanceModel.load_filtered_results(
                config_dto.working_directory)

            simulation_dto_list = []
            for result_log in result_log_dict_list:

                case = int(result_log["case"])

                for dto in all_filtered_reluctance_dtos:
                    if case == dto.case:
                        simulation_dto_list.append(dto)

            femmt.integrated_transformer_fem_thermal_simulations_from_result_dtos(config_dto, simulation_dto_list, visualize)

        @staticmethod
        def load_unfiltered_simulations(working_directory: str) -> List[Dict]:
            """
            Load all simulations from a given working directory.

            :param working_directory: directory to load
            :type working_directory: str
            """
            filepath = os.path.join(working_directory, "03_fem_thermal_simulation_results")

            return femmt.IntegratedTransformerOptimization.FemSimulation.load(filepath)

        @staticmethod
        def filter_max_temperature(thermal_result_log_list: List[Dict], max_core_temperature: float, max_winding_temperature: float) -> List[Dict]:
            """
            Filter out designs with too high core and winding temperatures.

            :param thermal_result_log_list: list of thermal result logs
            :type thermal_result_log_list: List[Dict]
            :param max_core_temperature: maximum core temperature
            :type max_core_temperature: float
            :param max_winding_temperature: maximum winding temperature
            :type max_winding_temperature: float
            """
            filtered_thermal_result_log_list = []

            for thermal_result_log in thermal_result_log_list:
                if thermal_result_log["core"]["max"] < max_core_temperature and thermal_result_log["windings"]["total"]["max"] < max_winding_temperature:
                    filtered_thermal_result_log_list.append(thermal_result_log)

            return filtered_thermal_result_log_list

        @staticmethod
        def find_common_cases(fem_results_list: List[Dict], thermal_fem_results_list: List[Dict]) -> List[Dict]:
            """
            Find out common cases in FEM simulation list and thermal FEM simulation list.

            :param fem_results_list: List of FEM result-logs
            :type fem_results_list: List[Dict]
            :param thermal_fem_results_list: List of thermal FEM result-logs
            :type thermal_fem_results_list: List[Dict]
            :return: List of FEM result-logs
            :rtype: List[Dict]
            """
            common_case_list = []
            for thermal_fem_result in thermal_fem_results_list:
                for fem_result in fem_results_list:
                    if thermal_fem_result["case"] == fem_result["case"]:
                        common_case_list.append(fem_result)
                        break

            return common_case_list
