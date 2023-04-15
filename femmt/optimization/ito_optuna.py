# python libraries
import os
import itertools
from typing import List, Dict, Tuple

# 3rd party libraries
import pandas as pd
import numpy as np
import optuna
import materialdatabase as mdb
from scipy import optimize

#femmt libraries
from femmt.optimization.integrated_transformer_dtos import *
import femmt.functions_reluctance as fr
import femmt.functions as ff
import femmt.optimization.optuna_femmt_parser as op
import femmt.optimization.ito_functions as itof
import femmt




class ItoOptuna:
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
            material_flux_density_vec, material_mu_r_imag_vec, material_mu_r_real_vec = material_db.permeability_data_to_pro_file(
                config.temperature, fundamental_frequency, material_name, datasource=mdb.MaterialDataSource.ManufacturerDatasheet,
                datatype='permeability_data', plot_interpolation=False)
            material_mu_r_initial = material_db.get_material_property(material_name=material_name, property="initial_permeability")
            saturation_flux_density = material_db.get_saturation_flux_density(material_name=material_name)

            material_data_list.append(MaterialCurve(material_name, material_mu_r_initial, material_flux_density_vec, material_mu_r_imag_vec, material_mu_r_real_vec, saturation_flux_density))

        # working directories
        fem_working_directory = os.path.join(config.working_directory, '00_femmt_simulation')
        reluctance_model_results_directory = os.path.join(config.working_directory, "01_reluctance_model_results")
        fem_simulation_results_directory = os.path.join(config.working_directory, "02_fem_simulation_results")
        fem_thermal_simulation_results_directory = os.path.join(config.working_directory, "03_fem_thermal_simulation_results")


        os.makedirs(fem_working_directory, exist_ok=True)
        os.makedirs(fem_simulation_results_directory, exist_ok=True)
        os.makedirs(reluctance_model_results_directory, exist_ok=True)
        os.makedirs(fem_thermal_simulation_results_directory, exist_ok=True)

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
            fem_working_directory=fem_working_directory,
            fem_simulation_results_directory=fem_simulation_results_directory,
            reluctance_model_results_directory=reluctance_model_results_directory,
            fem_thermal_simulation_results_directory=fem_thermal_simulation_results_directory,
        )

        return target_and_fix_parameters

    class ReluctanceModel:

        @staticmethod
        def objective(trial, config: ItoSingleInputConfig, target_and_fixed_parameters: ItoTargetAndFixedParameters) -> Tuple:
            """
            Objective function to optimize.

            :param config: input configuration file
            :type config: ItoSingleInputConfig
            :param target_and_fixed_parameters: target and fix parameters
            :type target_and_fixed_parameters: ItoTargetAndFixedParameters
            """

            # pass multible arguments to the objective function used by optuna
            # https://www.kaggle.com/general/261870
            # get parameters to tune
            core_inner_diameter = trial.suggest_float("core_inner_diameter", config.core_inner_diameter_min_max_list[0], config.core_inner_diameter_min_max_list[1])
            window_w = trial.suggest_float("window_w", config.window_w_min_max_list[0], config.window_w_min_max_list[1])
            window_h_top = trial.suggest_float("window_h_top", config.window_h_top_min_max_list[0], config.window_h_top_min_max_list[1])
            window_h_bot = trial.suggest_float("window_h_bot", config.window_h_bot_min_max_list[0], config.window_h_bot_min_max_list[1])

            n_p_top = trial.suggest_int("n_p_top", config.n_p_top_min_max_list[0], config.n_p_top_min_max_list[1])
            n_p_bot = trial.suggest_int("n_p_bot", config.n_p_bot_min_max_list[0], config.n_p_bot_min_max_list[1])
            n_s_top = trial.suggest_int("n_s_top", config.n_s_top_min_max_list[0], config.n_s_top_min_max_list[1])
            n_s_bot = trial.suggest_int("n_s_bot", config.n_s_bot_min_max_list[0], config.n_s_bot_min_max_list[1])

            material = trial.suggest_categorical("material", config.material_list)
            primary_litz_wire = trial.suggest_categorical("primary_litz_wire", config.primary_litz_wire_list)
            secondary_litz_wire = trial.suggest_categorical("secondary_litz_wire", config.secondary_litz_wire_list)

            thousand_simulations = trial.number / 1000

            if thousand_simulations.is_integer():
                print(f"simulation count: {trial.number}")

            litz_database = ff.litz_database()


            for material_dto in target_and_fixed_parameters.material_dto_curve_list:
                if material_dto.material_name ==  material:
                    material_data = material_dto

            material_mu_r_initial = material_data.material_mu_r_abs
            flux_density_data_vec = material_data.material_flux_density_vec
            mu_r_imag_data_vec = material_data.material_mu_r_imag_vec

            window_h_total = window_h_top + window_h_bot + core_inner_diameter / 4

            core_top_bot_height = core_inner_diameter / 4
            core_cross_section = (core_inner_diameter / 2) ** 2 * np.pi

            # generate winding matrix
            # note that the t11 winding-matrix will be reshaped later!
            t2_winding_matrix = [[n_p_top, n_s_top], [n_p_bot, n_s_bot]]

            target_inductance_matrix = fr.calculate_inductance_matrix_from_ls_lh_n(config.l_s_target, config.l_h_target,
                                                                                   config.n_target)
            t2_reluctance_matrix = fr.calculate_reluctance_matrix(t2_winding_matrix, target_inductance_matrix)

            if np.linalg.det(t2_reluctance_matrix) != 0 and np.linalg.det(
                    np.transpose(t2_winding_matrix)) != 0 and np.linalg.det(target_inductance_matrix) != 0:
                # calculate the flux
                flux_top_vec, flux_bot_vec, flux_stray_vec = fr.flux_vec_from_current_vec(target_and_fixed_parameters.current_extracted_1_vec,
                                                                                          target_and_fixed_parameters.current_extracted_2_vec,
                                                                                          t2_winding_matrix,
                                                                                          target_inductance_matrix)

                # calculate maximum values
                flux_top_max, flux_bot_max, flux_stray_max = fr.max_value_from_value_vec(flux_top_vec, flux_bot_vec,
                                                                                         flux_stray_vec)

                flux_density_top_max = flux_top_max / core_cross_section
                flux_density_bot_max = flux_bot_max / core_cross_section
                flux_density_middle_max = flux_stray_max / core_cross_section

                if (flux_density_top_max < material_data.saturation_flux_density) and (
                        flux_density_bot_max < material_data.saturation_flux_density) and (
                        flux_density_middle_max < material_data.saturation_flux_density):

                    # calculate target values for r_top and r_bot out of reluctance matrix
                    r_core_middle_cylinder_radial = fr.r_core_top_bot_radiant(core_inner_diameter, window_w, material_data.material_mu_r_abs,
                                                                              core_top_bot_height)

                    r_middle_target = -t2_reluctance_matrix[0][1]
                    r_top_target = t2_reluctance_matrix[0][0] - r_middle_target
                    r_bot_target = t2_reluctance_matrix[1][1] - r_middle_target

                    # calculate the core reluctance of top and bottom and middle part
                    r_core_top_cylinder_inner = fr.r_core_round(core_inner_diameter, window_h_top, material_data.material_mu_r_abs)
                    r_core_top = 2 * r_core_top_cylinder_inner + r_core_middle_cylinder_radial
                    r_air_gap_top_target = r_top_target - r_core_top

                    r_core_bot_cylinder_inner = fr.r_core_round(core_inner_diameter, window_h_bot, material_data.material_mu_r_abs)
                    r_core_bot = 2 * r_core_bot_cylinder_inner + r_core_middle_cylinder_radial
                    r_air_gap_bot_target = r_bot_target - r_core_bot

                    r_air_gap_middle_target = r_middle_target - r_core_middle_cylinder_radial

                    if r_top_target > 0 and r_bot_target > 0 and r_middle_target > 0:

                        minimum_air_gap_length = 0
                        maximum_air_gap_length = 10e-3
                        minimum_sort_out_air_gap_length = 1e-15
                        try:
                            l_top_air_gap = optimize.brentq(fr.r_air_gap_round_inf_sct, minimum_air_gap_length,
                                                            maximum_air_gap_length,
                                                            args=(core_inner_diameter, window_h_top, r_air_gap_top_target))
                            #print(f"{l_top_air_gap = }")
                            l_bot_air_gap = optimize.brentq(fr.r_air_gap_round_round_sct, minimum_air_gap_length,
                                                            maximum_air_gap_length, args=(
                                core_inner_diameter, window_h_bot / 2, window_h_bot / 2, r_air_gap_bot_target))
                            #print(f"{l_bot_air_gap = }")
                            l_middle_air_gap = optimize.brentq(fr.r_air_gap_tablet_cylinder_sct, minimum_air_gap_length,
                                                               maximum_air_gap_length, args=(
                                core_inner_diameter, core_inner_diameter / 4, window_w, r_air_gap_middle_target))
                        except:
                            l_top_air_gap = 0
                            l_bot_air_gap = 0
                            l_middle_air_gap = 0

                        if l_top_air_gap >= minimum_sort_out_air_gap_length and l_bot_air_gap >= minimum_sort_out_air_gap_length and l_middle_air_gap >= minimum_sort_out_air_gap_length:
                            p_hyst_top = fr.hyst_losses_core_half_mu_r_imag(core_inner_diameter, window_h_top, window_w,
                                                                            material_data.material_mu_r_abs, flux_top_max, target_and_fixed_parameters.fundamental_frequency,
                                                                            material_data.material_flux_density_vec, material_data.material_mu_r_imag_vec)

                            p_hyst_middle = fr.power_losses_hysteresis_cylinder_radial_direction_mu_r_imag(flux_stray_max,
                                                                                                           core_inner_diameter / 4,
                                                                                                           core_inner_diameter / 2,
                                                                                                           core_inner_diameter / 2 + window_w,
                                                                                                           target_and_fixed_parameters.fundamental_frequency,
                                                                                                           material_data.material_mu_r_abs,
                                                                                                           material_data.material_flux_density_vec,
                                                                                                           material_data.material_mu_r_imag_vec)

                            p_hyst_bot = fr.hyst_losses_core_half_mu_r_imag(core_inner_diameter, window_h_bot, window_w,
                                                                            material_data.material_mu_r_abs, flux_bot_max,
                                                                            target_and_fixed_parameters.fundamental_frequency,
                                                                            material_data.material_flux_density_vec,
                                                                            material_data.material_mu_r_imag_vec)

                            p_hyst = p_hyst_top + p_hyst_bot + p_hyst_middle

                            core_2daxi_total_volume = fr.calculate_core_2daxi_total_volume(core_inner_diameter, (
                                        window_h_bot + window_h_top + core_inner_diameter / 4), window_w)


                            primary_litz = litz_database[primary_litz_wire]
                            secondary_litz = litz_database[secondary_litz_wire]

                            # cross-section comparison is according to a square for round wire.
                            # this approximation is more realistic
                            # insulation
                            insulation_distance = 1e-3
                            insulation_cross_section_top = 2 * insulation_distance * (window_w + window_h_top)
                            insulation_cross_section_bot = 2 * insulation_distance * (window_w + window_h_bot)

                            total_available_window_cross_section_top = window_h_top * window_w - insulation_cross_section_top
                            winding_cross_section_top = n_p_top * (2 * primary_litz["conductor_radii"]) ** 2 + n_s_top * (
                                        2 * secondary_litz["conductor_radii"]) ** 2

                            total_available_window_cross_section_bot = window_h_bot * window_w - insulation_cross_section_bot
                            winding_cross_section_bot = n_p_bot * (2 * primary_litz["conductor_radii"]) ** 2 + n_s_bot * (
                                        2 * secondary_litz["conductor_radii"]) ** 2



                            if winding_cross_section_top < total_available_window_cross_section_top and winding_cross_section_bot < total_available_window_cross_section_bot:
                                primary_effective_conductive_cross_section = primary_litz["strands_numbers"] * primary_litz[
                                    "strand_radii"] ** 2 * np.pi
                                primary_effective_conductive_radius = np.sqrt(
                                    primary_effective_conductive_cross_section / np.pi)
                                primary_resistance = fr.resistance_solid_wire(core_inner_diameter, window_w,
                                                                              n_p_top + n_p_bot,
                                                                              primary_effective_conductive_radius,
                                                                              material='Copper')
                                primary_dc_loss = primary_resistance * target_and_fixed_parameters.i_rms_1 ** 2

                                secondary_effective_conductive_cross_section = secondary_litz["strands_numbers"] * \
                                                                               secondary_litz["strand_radii"] ** 2 * np.pi
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
                                trial.set_user_attr("flux_density_top_max",flux_density_top_max)
                                trial.set_user_attr("flux_density_bot_max", flux_density_bot_max)
                                trial.set_user_attr("flux_density_stray_max", flux_density_middle_max)
                                trial.set_user_attr("p_hyst", p_hyst)
                                trial.set_user_attr("primary_litz_wire_loss", primary_dc_loss)
                                trial.set_user_attr("secondary_litz_wire_loss", secondary_dc_loss)

                                print(f"succesfully calculated trial {trial.number}")

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
                                return None, None
                        else:
                            return None, None
                    else:
                        return None, None
                else:
                    return None, None
            else:
                return None, None

        @staticmethod
        def start_study(study_name: str, config: ItoSingleInputConfig, number_trials: int, storage: str = None) -> None:
            """
            Start a study to optimize an integrated transformer

            Note: Due to performance reasons, the study is calculated in RAM.
            After finishing the study, the results are copied to sqlite or mysql database by the use of a new study.

            :param study_name: Name of the study
            :type study_name: str
            :param config: simulation configuration
            :type config: ItoSingleInputConfig
            :param number_trials: number of trials
            :type number_trials: int
            :param storage: "sqlite" or "mysql"
            :type storage: str
            """

            target_and_fixed_parameters = femmt.optimization.ItoOptuna.calculate_fix_parameters(config)

            # Wrap the objective inside a lambda and call objective inside it
            func = lambda trial: femmt.optimization.ItoOptuna.objective(trial, config, target_and_fixed_parameters)

            if storage == 'sqlite':
                storage= f"sqlite:///study_{study_name}.sqlite3"
            elif storage == 'mysql':
                storage= "mysql://monty@localhost/mydb",


            # Pass func to Optuna studies
            study_in_memory = optuna.create_study(directions=["minimize", "minimize"],
                                        sampler=optuna.samplers.NSGAIISampler(),
                                        )

            # set logging verbosity: https://optuna.readthedocs.io/en/stable/reference/generated/optuna.logging.set_verbosity.html#optuna.logging.set_verbosity
            # .INFO: all messages (default)
            # .WARNING: fails and warnings
            # .ERROR: only errors
            optuna.logging.set_verbosity(optuna.logging.ERROR)


            print(f"Sampler is {study_in_memory.sampler.__class__.__name__}")
            study_in_memory.optimize(func, n_trials=number_trials, n_jobs=-1)

            # in-memory calculation is shown before saving the data to database
            fig = optuna.visualization.plot_pareto_front(study_in_memory, target_names=["volume", "losses"])
            fig.show()

            study_in_storage = optuna.create_study(directions=["minimize", "minimize"], study_name=study_name,
                                               storage=storage)
            study_in_storage.add_trials(study_in_memory.trials)

        @staticmethod
        def proceed_study(study_name: str, config: ItoSingleInputConfig, number_trials: int) -> None:
            """
            Proceed a study which is stored as sqlite database.

            :param study_name: Name of the study
            :type study_name: str
            :param config: Simulation configuration
            :type ItoSingleInputConfig
            :param number_trials: Number of trials adding to the existing study
            :type number_trials: int
            """
            target_and_fixed_parameters = femmt.optimization.ItoOptuna.calculate_fix_parameters(config)

            # Wrap the objective inside a lambda and call objective inside it
            func = lambda trial: femmt.optimization.ItoOptuna.objective(trial, config, target_and_fixed_parameters)

            study = optuna.create_study(study_name=study_name, storage=f"sqlite:///study_{study_name}.sqlite3", load_if_exists=True)
            study.optimize(func, n_trials=number_trials)

        @staticmethod
        def show_study_results(study_name: str) -> None:
            """
            Show the results of a study.

            :param study_name: Name of the study
            :type study_name: str
            """
            study = optuna.create_study(study_name=study_name, storage=f"sqlite:///study_{study_name}.sqlite3",
                                        load_if_exists=True)

            fig = optuna.visualization.plot_pareto_front(study, target_names=["volume", "losses"])
            fig.show()

        @staticmethod
        def load_study_best_trials(study_name: str) -> List[ItoSingleResultFile]:
            """
            Load the best trials (Pareto front) of a study.

            :param study_name: Name of the study
            :type study_name: str
            :return: List of the best trials.
            :rtype: List[ItoSingleResultFile]
            """
            study = optuna.create_study(study_name=study_name, storage=f"sqlite:///study_{study_name}.sqlite3",
                                        load_if_exists=True)

            print(f"{len(study.trials) = }")
            print(f"{len(study.best_trials) =}")
            dto_list = [op.OptunaFemmtParser.parse(frozen_object) for frozen_object in study.best_trials]

            return dto_list

        @staticmethod
        def load_study_pareto_area(study_name, loss_offset) -> List[ItoSingleResultFile]:
            """
            Load an optuna study and return the pareto front including an offset in losses for further calculations (e.g. FEM simulations)

            :param study_name: name of study to load
            :type study_name: str
            :param loss_offset: loss-offset to be added to the pareto-front losses
            :type loss_offset: float
            """
            study = optuna.create_study(study_name=study_name, storage=f"sqlite:///study_{study_name}.sqlite3", load_if_exists=True)
            full_dto_list = [op.OptunaFemmtParser.parse(frozen_object) for frozen_object in study.trials]

            pareto_dto_list = [op.OptunaFemmtParser.parse(frozen_object) for frozen_object in study.best_trials]

            x_pareto_vec, y_pareto_vec = itof.dto_list_to_vec(pareto_dto_list)

            filtered_design_dto_list = []
            for dto in full_dto_list:
                ref_loss = np.interp(dto.core_2daxi_total_volume, x_pareto_vec, y_pareto_vec) + loss_offset
                if dto.total_loss < ref_loss:
                    filtered_design_dto_list.append(dto)

            return filtered_design_dto_list



    class FemSimulation:
        @staticmethod
        def simulate(config_dto: ItoSingleInputConfig,
                     reluctance_result_list: List[ItoSingleResultFile],
                     visualize = False):
            """
            Perform FEM simulations out of input files.
            """

            itof.integrated_transformer_fem_simulations_from_result_dtos(config_dto,
                                                                reluctance_result_list,
                                                                visualize = visualize)

    class ThermalSimulation:
        pass


