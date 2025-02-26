"""Integrated transformer optimization."""
# Python libraries
import os
import json
import dataclasses
from typing import List, Dict, Tuple, Optional
import datetime
import pickle

# 3rd party library import
import materialdatabase as mdb
from scipy import optimize
import optuna

# femmt import
import femmt.functions as ff
import femmt.functions_reluctance as fr
import femmt.optimization.functions_optimization as fo
from femmt.optimization.ito_dtos import *
import femmt.optimization.optuna_femmt_parser as op
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
        primary_litz_wire=result_file_dict["primary_litz_wire"],
        secondary_litz_wire=result_file_dict["secondary_litz_wire"],
        primary_litz_wire_loss=result_file_dict["primary_litz_wire_loss"],
        secondary_litz_wire_loss=result_file_dict["secondary_litz_wire_loss"],
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
        """
        Calculate fix parameters what can be derived from the input configuration.

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

        # target inductances
        target_inductance_matrix = fr.calculate_inductance_matrix_from_ls_lh_n(config.l_s_target, config.l_h_target,
                                                                               config.n_target)

        # material properties
        material_db = mdb.MaterialDatabase(is_silent=True)

        material_data_list = []
        for material_name in config.material_list:
            material_dto = material_db.material_data_interpolation_to_dto(material_name, fundamental_frequency, config.temperature)
            material_data_list.append(material_dto)

        # set up working directories
        working_directories = itof.set_up_folder_structure(config.integrated_transformer_optimization_directory)

        # finalize data to dto
        target_and_fix_parameters = ItoTargetAndFixedParameters(
            i_rms_1=i_rms_1,
            i_rms_2=i_rms_2,
            time_extracted_vec=time_extracted,
            current_extracted_1_vec=current_extracted_1_vec,
            current_extracted_2_vec=current_extracted_2_vec,
            material_dto_curve_list=material_data_list,
            fundamental_frequency=fundamental_frequency,
            target_inductance_matrix=target_inductance_matrix,
            working_directories=working_directories
        )

        return target_and_fix_parameters

    class ReluctanceModel:
        """Create and calculate the reluctance model for the integrated transformer."""

        @staticmethod
        def objective(trial: optuna.Trial, config: ItoSingleInputConfig,
                      target_and_fixed_parameters: ItoTargetAndFixedParameters) -> Tuple:
            """
            Objective function to optimize.

            Using optuna. Some hints:

             * returning failed trails by using return float('nan'), float('nan'),
               see https://optuna.readthedocs.io/en/stable/faq.html#how-are-nans-returned-by-trials-handled
             * speed up the search for NSGA-II algorithm with dynamic alter the search space, see https://optuna.readthedocs.io/en/stable/faq.html#id10


            :param trial: parameter suggesting by optuna
            :type trial: optuna.Trial
            :param config: input configuration file
            :type config: ItoSingleInputConfig
            :param target_and_fixed_parameters: target and fix parameters
            :type target_and_fixed_parameters: ItoTargetAndFixedParameters
            """
            # pass multiple arguments to the objective function used by optuna
            # https://www.kaggle.com/general/261870

            #########################################################
            # set core geometry optimization parameters
            #########################################################
            core_inner_diameter = trial.suggest_float("core_inner_diameter",
                                                      config.core_inner_diameter_min_max_list[0],
                                                      config.core_inner_diameter_min_max_list[1])
            window_w = trial.suggest_float("window_w", config.window_w_min_max_list[0],
                                           config.window_w_min_max_list[1])
            window_h_top = trial.suggest_float("window_h_top", config.window_h_top_min_max_list[0],
                                               config.window_h_top_min_max_list[1])
            window_h_bot = trial.suggest_float("window_h_bot", config.window_h_bot_min_max_list[0],
                                               config.window_h_bot_min_max_list[1])

            material = trial.suggest_categorical("material", config.material_list)
            primary_litz_wire = trial.suggest_categorical("primary_litz_wire", config.primary_litz_wire_list)
            secondary_litz_wire = trial.suggest_categorical("secondary_litz_wire", config.secondary_litz_wire_list)

            # cross-section comparison is according to a square for round wire.
            # this approximation is more realistic
            # insulation
            insulation_distance = 1e-3
            insulation_cross_section_top = 2 * insulation_distance * (window_w + window_h_top)
            insulation_cross_section_bot = 2 * insulation_distance * (window_w + window_h_bot)

            litz_database = ff.litz_database()

            primary_litz = litz_database[primary_litz_wire]
            secondary_litz = litz_database[secondary_litz_wire]

            total_available_window_cross_section_top = window_h_top * window_w - insulation_cross_section_top
            total_available_window_cross_section_bot = window_h_bot * window_w - insulation_cross_section_bot

            #########################################################
            # set dynamic wire count parameters as optimization parameters
            #########################################################
            # set the winding search space dynamic
            # https://optuna.readthedocs.io/en/stable/faq.html#what-happens-when-i-dynamically-alter-a-search-space

            # n_p_top suggestion
            n_p_top_max = total_available_window_cross_section_top / (2 * primary_litz["conductor_radii"]) ** 2
            n_p_top = trial.suggest_int("n_p_top", 0, n_p_top_max)

            # n_s_top_suggestion
            winding_cross_section_n_p_top_max = n_p_top * (2 * primary_litz["conductor_radii"]) ** 2
            n_s_top_max = int((total_available_window_cross_section_top - winding_cross_section_n_p_top_max) / (
                2 * secondary_litz["conductor_radii"]) ** 2)
            n_s_top = trial.suggest_int("n_s_top", 0, n_s_top_max)

            # n_p_bot suggestion
            n_p_bot_max = total_available_window_cross_section_bot / (2 * primary_litz["conductor_radii"]) ** 2
            n_p_bot = trial.suggest_int("n_p_bot", 0, n_p_bot_max)

            # n_s_bot suggestion
            winding_cross_section_n_p_bot_max = n_p_bot * (2 * primary_litz["conductor_radii"]) ** 2
            n_s_bot_max = int((total_available_window_cross_section_bot - winding_cross_section_n_p_bot_max) / (2 * secondary_litz["conductor_radii"]) ** 2)
            n_s_bot = trial.suggest_int("n_s_bot", 0, n_s_bot_max)

            winding_cross_section_top = n_p_top * (2 * primary_litz["conductor_radii"]) ** 2 + n_s_top * (2 * secondary_litz["conductor_radii"]) ** 2
            winding_cross_section_bot = n_p_bot * (2 * primary_litz["conductor_radii"]) ** 2 + n_s_bot * (2 * secondary_litz["conductor_radii"]) ** 2

            thousand_simulations = trial.number / 1000

            if thousand_simulations.is_integer():
                print(f"simulation count: {trial.number}")

            for material_dto in target_and_fixed_parameters.material_dto_curve_list:
                if material_dto.material_name == material:
                    material_data = material_dto

                material_mu_r_initial = material_data.material_mu_r_abs
                flux_density_data_vec = material_data.material_flux_density_vec
                mu_r_imag_data_vec = material_data.material_mu_r_imag_vec

                core_top_bot_height = core_inner_diameter / 4
                core_cross_section = (core_inner_diameter / 2) ** 2 * np.pi

                t2_winding_matrix = [[n_p_top, n_s_top], [n_p_bot, n_s_bot]]

                target_inductance_matrix = fr.calculate_inductance_matrix_from_ls_lh_n(config.l_s_target,
                                                                                       config.l_h_target,
                                                                                       config.n_target)
                t2_reluctance_matrix = fr.calculate_reluctance_matrix(t2_winding_matrix, target_inductance_matrix)

                core_2daxi_total_volume = fr.calculate_core_2daxi_total_volume(core_inner_diameter,
                                                                               (window_h_bot + window_h_top + core_inner_diameter / 4), window_w)

                if np.linalg.det(t2_reluctance_matrix) != 0 and np.linalg.det(
                        np.transpose(t2_winding_matrix)) != 0 and np.linalg.det(target_inductance_matrix) != 0:
                    # calculate the flux
                    flux_top_vec, flux_bot_vec, flux_stray_vec = fr.flux_vec_from_current_vec(
                        target_and_fixed_parameters.current_extracted_1_vec,
                        target_and_fixed_parameters.current_extracted_2_vec,
                        t2_winding_matrix,
                        target_inductance_matrix)

                    # calculate maximum values
                    flux_top_max, flux_bot_max, flux_stray_max = fr.max_value_from_value_vec(flux_top_vec, flux_bot_vec,
                                                                                             flux_stray_vec)

                    flux_density_top_max = flux_top_max / core_cross_section
                    flux_density_bot_max = flux_bot_max / core_cross_section
                    flux_density_middle_max = flux_stray_max / core_cross_section

                    # calculate target values for r_top and r_bot out of reluctance matrix
                    r_core_middle_cylinder_radial = fr.r_core_top_bot_radiant(core_inner_diameter, window_w,
                                                                              material_data.material_mu_r_abs,
                                                                              core_top_bot_height)

                    r_middle_target = -t2_reluctance_matrix[0][1]
                    r_top_target = t2_reluctance_matrix[0][0] - r_middle_target
                    r_bot_target = t2_reluctance_matrix[1][1] - r_middle_target

                    # calculate the core reluctance of top and bottom and middle part
                    r_core_top_cylinder_inner = fr.r_core_round(core_inner_diameter, window_h_top,
                                                                material_data.material_mu_r_abs)
                    r_core_top = 2 * r_core_top_cylinder_inner + r_core_middle_cylinder_radial
                    r_air_gap_top_target = r_top_target - r_core_top

                    r_core_bot_cylinder_inner = fr.r_core_round(core_inner_diameter, window_h_bot,
                                                                material_data.material_mu_r_abs)
                    r_core_bot = 2 * r_core_bot_cylinder_inner + r_core_middle_cylinder_radial
                    r_air_gap_bot_target = r_bot_target - r_core_bot

                    r_air_gap_middle_target = r_middle_target - r_core_middle_cylinder_radial

                    if r_air_gap_top_target > 0 and r_air_gap_bot_target > 0 and r_air_gap_middle_target > 0:

                        # Note: a minimum air gap length of zero is not allowed. This will lead to failure calculation
                        # when trying to solve (using brentq) r_gap_round_round-function. Calculating an air gap
                        # reluctance with length of zero is not realistic.
                        minimum_air_gap_length = 1e-15
                        maximum_air_gap_length = 5e-3
                        minimum_sort_out_air_gap_length = 0
                        try:
                            # solving brentq needs to be in try/except statement,
                            # as it can be that there is no sign changing in the given interval
                            # to search for the zero.
                            # Note: setting full output to true and taking object [0] is only
                            # to avoid linting error!
                            l_top_air_gap = optimize.brentq(fr.r_air_gap_round_inf_sct, minimum_air_gap_length, maximum_air_gap_length,
                                                            args=(core_inner_diameter, window_h_top, r_air_gap_top_target), full_output=True)[0]

                            l_bot_air_gap = optimize.brentq(fr.r_air_gap_round_round_sct, minimum_air_gap_length,
                                                            maximum_air_gap_length, args=(core_inner_diameter, window_h_bot / 2, window_h_bot / 2,
                                                                                          r_air_gap_bot_target), full_output=True)[0]

                            l_middle_air_gap = optimize.brentq(fr.r_air_gap_tablet_cylinder_sct, minimum_air_gap_length,
                                                               maximum_air_gap_length, args=(core_inner_diameter,
                                                                                             core_inner_diameter / 4, window_w, r_air_gap_middle_target),
                                                               full_output=True)[0]

                        except ValueError:
                            # ValueError is raised in case of an air gap with length of zero
                            return float('nan'), float('nan')

                        if l_top_air_gap > core_inner_diameter or l_bot_air_gap > core_inner_diameter:
                            return float('nan'), float('nan')

                        if l_top_air_gap >= minimum_sort_out_air_gap_length and l_bot_air_gap >= minimum_sort_out_air_gap_length \
                                and l_middle_air_gap >= minimum_sort_out_air_gap_length:
                            p_hyst_top = fr.hyst_losses_core_half_mu_r_imag(core_inner_diameter, window_h_top, window_w,
                                                                            material_data.material_mu_r_abs,
                                                                            flux_top_max,
                                                                            target_and_fixed_parameters.fundamental_frequency,
                                                                            material_data.material_flux_density_vec,
                                                                            material_data.material_mu_r_imag_vec)

                            p_hyst_middle = fr.power_losses_hysteresis_cylinder_radial_direction_mu_r_imag(
                                flux_stray_max,
                                core_inner_diameter / 4,
                                core_inner_diameter / 2,
                                core_inner_diameter / 2 + window_w,
                                target_and_fixed_parameters.fundamental_frequency,
                                material_data.material_mu_r_abs,
                                material_data.material_flux_density_vec,
                                material_data.material_mu_r_imag_vec)

                            p_hyst_bot = fr.hyst_losses_core_half_mu_r_imag(core_inner_diameter, window_h_bot, window_w,
                                                                            material_data.material_mu_r_abs,
                                                                            flux_bot_max,
                                                                            target_and_fixed_parameters.fundamental_frequency,
                                                                            material_data.material_flux_density_vec,
                                                                            material_data.material_mu_r_imag_vec)

                            p_hyst = p_hyst_top + p_hyst_bot + p_hyst_middle

                            primary_effective_conductive_cross_section = primary_litz["strands_numbers"] * primary_litz[
                                "strand_radii"] ** 2 * np.pi
                            primary_effective_conductive_radius = np.sqrt(
                                primary_effective_conductive_cross_section / np.pi)
                            primary_resistance = fr.resistance_solid_wire(core_inner_diameter, window_w,
                                                                          n_p_top + n_p_bot,
                                                                          primary_effective_conductive_radius,
                                                                          material='Copper')
                            primary_dc_loss = primary_resistance * target_and_fixed_parameters.i_rms_1 ** 2

                            secondary_effective_conductive_cross_section = secondary_litz["strands_numbers"] * secondary_litz["strand_radii"] ** 2 * np.pi
                            secondary_effective_conductive_radius = np.sqrt(
                                secondary_effective_conductive_cross_section / np.pi)
                            secondary_resistance = fr.resistance_solid_wire(core_inner_diameter, window_w,
                                                                            n_s_top + n_s_bot,
                                                                            secondary_effective_conductive_radius,
                                                                            material='Copper')
                            secondary_dc_loss = secondary_resistance * target_and_fixed_parameters.i_rms_2 ** 2

                            total_loss = p_hyst + primary_dc_loss + secondary_dc_loss

                            trial.set_user_attr("air_gap_top", l_top_air_gap)
                            trial.set_user_attr("air_gap_bot", l_bot_air_gap)
                            trial.set_user_attr("air_gap_middle", l_middle_air_gap)

                            trial.set_user_attr("flux_top_max", flux_top_max)
                            trial.set_user_attr("flux_bot_max", flux_bot_max)
                            trial.set_user_attr("flux_stray_max", flux_stray_max)
                            trial.set_user_attr("flux_density_top_max", flux_density_top_max)
                            trial.set_user_attr("flux_density_bot_max", flux_density_bot_max)
                            trial.set_user_attr("flux_density_stray_max", flux_density_middle_max)
                            trial.set_user_attr("p_hyst", p_hyst)
                            trial.set_user_attr("primary_litz_wire_loss", primary_dc_loss)
                            trial.set_user_attr("secondary_litz_wire_loss", secondary_dc_loss)

                            print(f"successfully calculated trial {trial.number}")

                            valid_design_dict = ItoSingleResultFile(
                                case=trial.number,
                                air_gap_top=l_top_air_gap,
                                air_gap_bot=l_bot_air_gap,
                                air_gap_middle=l_middle_air_gap,
                                n_p_top=n_p_top,
                                n_p_bot=n_p_bot,
                                n_s_top=n_s_top,
                                n_s_bot=n_s_bot,
                                window_h_top=window_h_top,
                                window_h_bot=window_h_bot,
                                window_w=window_w,
                                core_material=material_data.material_name,
                                core_inner_diameter=core_inner_diameter,
                                primary_litz_wire=primary_litz_wire,
                                secondary_litz_wire=secondary_litz_wire,
                                # results
                                flux_top_max=flux_top_max,
                                flux_bot_max=flux_bot_max,
                                flux_stray_max=flux_stray_max,
                                flux_density_top_max=flux_density_top_max,
                                flux_density_bot_max=flux_density_bot_max,
                                flux_density_stray_max=flux_density_middle_max,
                                p_hyst=p_hyst,
                                core_2daxi_total_volume=core_2daxi_total_volume,
                                primary_litz_wire_loss=primary_dc_loss,
                                secondary_litz_wire_loss=secondary_dc_loss,
                                total_loss=total_loss

                            )

                            return core_2daxi_total_volume, total_loss
                        else:
                            return float('nan'), float('nan')
                    else:
                        return float('nan'), float('nan')
                else:
                    return float('nan'), float('nan')

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

        ##############################
        # load
        ##############################

        @staticmethod
        def load_study_to_dto(study_name: str, config: ItoSingleInputConfig) -> List[ItoSingleResultFile]:
            """
            Load all trials of a study to a DTO-list.

            :param study_name: Name of the study
            :type study_name: str
            :param config: Integrated transformer configuration file
            :type config: ItoSingleInputConfig
            :return: List of all trials
            :rtype: List[ItoSingleResultFile]

            """
            study = optuna.create_study(study_name=study_name,
                                        storage=f"sqlite:///{config.working_directory}/study_{study_name}.sqlite3",
                                        load_if_exists=True)

            dto_list = [op.OptunaFemmtParser.parse(frozen_object) for frozen_object in study.trials \
                        if frozen_object.state == optuna.trial.TrialState.COMPLETE]

            return dto_list

        @staticmethod
        def load_study_best_trials_to_dto(study_name: str, config: ItoSingleInputConfig) -> List[ItoSingleResultFile]:
            """
            Load the best trials (Pareto front) of a study.

            :param study_name: Name of the study
            :type study_name: str
            :param config: Integrated transformer configuration file
            :type config: ItoSingleInputConfig
            :return: List of the best trials.
            :rtype: List[ItoSingleResultFile]

            """
            study = optuna.create_study(study_name=study_name,
                                        storage=f"sqlite:///{config.working_directory}/study_{study_name}.sqlite3",
                                        load_if_exists=True)

            print(study.best_trials[0])

            dto_list = [op.OptunaFemmtParser.parse(frozen_object) for frozen_object in study.best_trials]

            return dto_list

        #############################
        # filter
        #############################

        @staticmethod
        def filter_loss_list(valid_design_list: List[ItoSingleResultFile], factor_min_dc_losses: float = 1.2) -> List[ItoSingleResultFile]:
            """
            Remove designs with too high losses compared to the minimum losses.

            :param valid_design_list: list of valid DTOs
            :type valid_design_list: List[ItoSingleResultFile]
            :param factor_min_dc_losses: filter factor for the minimum dc losses
            :type factor_min_dc_losses: float
            :returns: list with removed objects (too small air gaps)
            :rtype: List[ItoSingleResultFile]
            """
            # figure out pareto front
            # pareto_volume_list, pareto_core_hyst_list, pareto_dto_list = self.pareto_front(volume_list, core_hyst_loss_list, valid_design_list)

            x_pareto_vec, y_pareto_vec = fo.pareto_front_from_dtos(valid_design_list)

            vector_to_sort = np.array([x_pareto_vec, y_pareto_vec])

            # sorting 2d array by 1st row
            # https://stackoverflow.com/questions/49374253/sort-a-numpy-2d-array-by-1st-row-maintaining-columns
            sorted_vector = vector_to_sort[:, vector_to_sort[0].argsort()]
            x_pareto_vec = sorted_vector[0]
            y_pareto_vec = sorted_vector[1]

            total_losses_list = []
            filtered_design_dto_list = []

            for dto in valid_design_list:
                total_losses_list.append(dto.total_loss)

            min_total_dc_losses = total_losses_list[np.argmin(total_losses_list)]
            loss_offset = factor_min_dc_losses * min_total_dc_losses

            for dto in valid_design_list:
                ref_loss = np.interp(dto.core_2daxi_total_volume, x_pareto_vec, y_pareto_vec) + loss_offset
                if dto.total_loss < ref_loss:
                    filtered_design_dto_list.append(dto)

            return filtered_design_dto_list

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

        @staticmethod
        def save_dto_list(result_dto_list: List[ItoSingleResultFile], filepath: str):
            """
            Save the ItoSingleResultFile-List to the file structure.

            :param result_dto_list:
            :type result_dto_list: List[ItoSingleResultFile]
            :param filepath: filepath
            :type filepath: str
            """
            if not os.path.exists(filepath):
                os.mkdir(filepath)

            for _, dto in enumerate(result_dto_list):
                file_name = os.path.join(filepath, f"case_{dto.case}.json")

                result_dict = dataclasses.asdict(dto)
                with open(file_name, "w+", encoding='utf-8') as outfile:
                    json.dump(result_dict, outfile, indent=2, ensure_ascii=False, cls=MyJSONEncoder)

        @staticmethod
        def save_unfiltered_results(config_file: ItoSingleInputConfig, result_file_list: List[ItoSingleResultFile]):
            """
            Save the results of the reluctance model into the file structure.

            :param config_file: integrated transformer configuration file
            :type config_file: ItoSingleInputConfig
            :param result_file_list: list of ItoSingleResultFiles
            :type result_file_list: List[ItoSingleResultFile]
            """
            # generate folder structure
            femmt.set_up_folder_structure(config_file.working_directory)

            # save optimization input parameters
            config_dict = dataclasses.asdict(config_file)

            integrated_transformer_optimization_input_parameters_file = os.path.join(config_file.working_directory, "optimization_input_parameters.json")
            integrated_transformer_reluctance_model_results_directory = os.path.join(config_file.working_directory, "01_reluctance_model_results")

            with open(integrated_transformer_optimization_input_parameters_file, "w+", encoding='utf-8') as outfile:
                json.dump(config_dict, outfile, indent=2, ensure_ascii=False, cls=MyJSONEncoder)

            # save reluctance parameters winning candidates
            femmt.IntegratedTransformerOptimization.ReluctanceModel.save_dto_list(result_file_list, integrated_transformer_reluctance_model_results_directory)

        @staticmethod
        def load_list(filepath: str) -> List[ItoSingleResultFile]:
            """
            Load the list of the reluctance models from the folder structure.

            :param filepath: filepath
            :type filepath: str
            :return: List of ItoSingleResultFiles
            :rtype: List[ItoSingleResultFile]
            """
            valid_design_list = []
            for file in os.listdir(filepath):
                if file.endswith(".json"):
                    json_file_path = os.path.join(filepath, file)
                    with open(json_file_path, "r") as fd:
                        loaded_data_dict = json.loads(fd.read())

                    valid_design_list.append(result_file_dict_to_dto(loaded_data_dict))
            if len(valid_design_list) == 0:
                raise ValueError("Specified file path is empty")

            return valid_design_list

        @staticmethod
        def load_unfiltered_results(working_directory: str) -> List[ItoSingleResultFile]:
            """
            Load the results of the reluctance model and returns the ItoSingleResultFiles as a list.

            :param working_directory: working directory
            :type working_directory: str
            :return: List of ItoSingleResultFiles
            :rtype: List[ItoSingleResultFile]
            """
            integrated_transformer_reluctance_model_results_directory = os.path.join(working_directory, "01_reluctance_model_results")
            print(f"Read results from {integrated_transformer_reluctance_model_results_directory}")
            return femmt.IntegratedTransformerOptimization.ReluctanceModel.load_list(integrated_transformer_reluctance_model_results_directory)

        @staticmethod
        def load_filtered_results(working_directory: str) -> List[ItoSingleResultFile]:
            """
            Load the results of the reluctance model and returns the ItoSingleResultFiles as a list.

            :param working_directory: working directory
            :type working_directory: str
            :return: List of ItoSingleResultFiles
            :rtype: List[ItoSingleResultFile]
            """
            integrated_transformer_reluctance_model_results_directory = os.path.join(working_directory, "01_reluctance_model_results_filtered")
            print(f"Read results from {integrated_transformer_reluctance_model_results_directory}")
            return femmt.IntegratedTransformerOptimization.ReluctanceModel.load_list(integrated_transformer_reluctance_model_results_directory)

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
