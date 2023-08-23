# python libraries
import os
import shutil
import json

# 3rd party libraries
import optuna

# FEMMT and materialdatabase libraries
from femmt.optimization.sto_dtos import *
import femmt.functions_reluctance as fr
import femmt.functions as ff
import femmt.optimization.ito_functions as itof
import femmt
import materialdatabase as mdb

class StackedTransformerOptimization:

    @staticmethod
    def calculate_fix_parameters(config: StoSingleInputConfig) -> StoTargetAndFixedParameters:
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

        phi_deg_2 = phi_deg_2 - 180

        # target inductances
        target_inductance_matrix = fr.calculate_inductance_matrix_from_ls_lh_n(config.l_s12_target,
                                                                               config.l_h_target,
                                                                               config.n_target)

        # material properties
        material_db = mdb.MaterialDatabase(is_silent=True)

        material_data_list = []
        for material_name in config.material_list:
            material_dto = material_db.material_data_interpolation_to_dto(material_name, fundamental_frequency,
                                                                          config.temperature)
            material_data_list.append(material_dto)

        # set up working directories
        working_directories = itof.set_up_folder_structure(config.working_directory)

        # finalize data to dto
        target_and_fix_parameters = StoTargetAndFixedParameters(
            i_rms_1=i_rms_1,
            i_rms_2=i_rms_2,
            i_peak_1=i_peak_1,
            i_peak_2=i_peak_2,
            i_phase_deg_1=phi_deg_1,
            i_phase_deg_2=phi_deg_2,
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
        pass

    class FemSimulation:
        @staticmethod
        def model_pre_parameterization(trial, config: StoSingleInputConfig):
            parameters = {"core_inner_diameter": trial.suggest_float("core_inner_diameter", config.core_inner_diameter_min_max_list[0], config.core_inner_diameter_min_max_list[1]),
                          "window_w": trial.suggest_float("window_w", config.window_w_min_max_list[0], config.window_w_min_max_list[1]),
                          "air_gap_transformer": trial.suggest_float("air_gap_transformer", 0.1e-3, 5e-3), "inner_coil_insulation": trial.suggest_float("inner_coil_insulation", 0.5e-3, 5e-3),
                          "primary_litz_wire": trial.suggest_categorical("primary_litz_wire", config.primary_litz_wire_list)}

            # suggest core geometry

            parameters["primary_litz_parameters"] = ff.litz_database()[parameters["primary_litz_wire"]]
            parameters["primary_litz_diameter"] = 2 * parameters["primary_litz_parameters"]["conductor_radii"]

            # Will always be calculated from the given parameters
            parameters["available_width"] = parameters["window_w"] - parameters["inner_coil_insulation"] - config.insulations.iso_right_core

            # Suggestion of top window coil
            # Theoretically also 0 coil turns possible (number_rows_coil_winding must then be recalculated to avoid neg. values)
            parameters["primary_coil_turns"] = trial.suggest_int("primary_coil_turns", config.primary_coil_turns_min_max_list[0], config.primary_coil_turns_min_max_list[1])
            # Note: int() is used to round down.
            parameters["number_rows_coil_winding"] = int((parameters["primary_coil_turns"] * (parameters["primary_litz_diameter"] + config.insulations.iso_primary_to_primary) -
                                                          parameters["inner_coil_insulation"]) / parameters["available_width"]) + 1
            parameters["window_h_top"] = config.insulations.iso_top_core + config.insulations.iso_bot_core + \
                                         parameters["number_rows_coil_winding"] * parameters["primary_litz_diameter"] + \
                                         (parameters["number_rows_coil_winding"] - 1) * config.insulations.iso_primary_to_primary

            # Maximum coil air gap depends on the maximum window height top
            parameters["air_gap_coil"] = trial.suggest_float("air_gap_coil", 0.1e-3, parameters["window_h_top"] - 0.1e-3)

            # suggest categorical
            parameters["core_material"] = trial.suggest_categorical("material", config.material_list)
            parameters["foil_thickness"] = trial.suggest_categorical("foil_thickness", config.metal_sheet_thickness_list)

            try:
                if config.max_transformer_total_height is not None:
                    # Maximum transformer height
                    window_h_bot_max = config.max_transformer_total_height - 3 * parameters["core_inner_diameter"] / 4 - parameters["window_h_top"]
                    window_h_bot_min = config.window_h_bot_min_max_list[0]
                    if window_h_bot_min > window_h_bot_max:
                        print(f"{parameters['number_rows_coil_winding'] = }")
                        print(f"{parameters['window_h_top'] = }")
                        raise ValueError(f"{window_h_bot_min = } > {window_h_bot_max = }")

                    parameters["window_h_bot"] = trial.suggest_float("window_h_bot", window_h_bot_min, window_h_bot_max)
                else:
                    parameters["window_h_bot"] = trial.suggest_float("window_h_bot", config.window_h_bot_min_max_list[0], config.window_h_bot_min_max_list[1])
                return parameters

            except Exception as e:
                print(e)
                return None

        @staticmethod
        def fem_simulation(trial, trial_parameters, config: StoSingleInputConfig, target_and_fixed_parameters: StoTargetAndFixedParameters, flag_pre_visualize=False):
            geo = femmt.MagneticComponent(component_type=femmt.ComponentType.IntegratedTransformer,
                                          working_directory=target_and_fixed_parameters.working_directories.fem_working_directory,
                                          silent=True, simulation_name=f"Case_{trial.number}")

            core_dimensions = femmt.dtos.StackedCoreDimensions(core_inner_diameter=trial_parameters["core_inner_diameter"], window_w=trial_parameters["window_w"],
                                                               window_h_top=trial_parameters["window_h_top"], window_h_bot=trial_parameters["window_h_bot"])
            core = femmt.Core(core_type=femmt.CoreType.Stacked, core_dimensions=core_dimensions,
                              # mu_r_abs=3500, phi_mu_deg=12, sigma=1.2,
                              # permeability_datasource=femmt.MaterialDataSource.Custom,
                              # permittivity_datasource=femmt.MaterialDataSource.Custom)
                              material=trial_parameters["core_material"], temperature=config.temperature, frequency=target_and_fixed_parameters.fundamental_frequency,
                              permeability_datasource=femmt.MaterialDataSource.Measurement,
                              permeability_datatype=femmt.MeasurementDataType.ComplexPermeability,
                              permeability_measurement_setup="LEA_LK",
                              permittivity_datasource=femmt.MaterialDataSource.Measurement,
                              permittivity_datatype=femmt.MeasurementDataType.ComplexPermittivity,
                              permittivity_measurement_setup="LEA_LK")

            geo.set_core(core)

            air_gaps = femmt.AirGaps(femmt.AirGapMethod.Stacked, core)
            air_gaps.add_air_gap(femmt.AirGapLegPosition.CenterLeg, trial_parameters["air_gap_coil"], stacked_position=femmt.StackedPosition.Top)
            air_gaps.add_air_gap(femmt.AirGapLegPosition.CenterLeg, trial_parameters["air_gap_transformer"], stacked_position=femmt.StackedPosition.Bot)
            geo.set_air_gaps(air_gaps)

            # set_center_tapped_windings() automatically places the condu
            insulation, coil_window, transformer_window = femmt.functions_topologies.set_center_tapped_windings(
                core=core,

                # primary litz
                primary_additional_bobbin=1e-3,
                primary_turns=14,
                primary_radius=trial_parameters["primary_litz_parameters"]["conductor_radii"],
                primary_number_strands=trial_parameters["primary_litz_parameters"]["strands_numbers"],
                primary_strand_radius=trial_parameters["primary_litz_parameters"]["strand_radii"],

                # secondary foil
                secondary_parallel_turns=2,
                secondary_thickness_foil=trial_parameters["foil_thickness"],

                # insulation
                iso_top_core=config.insulations.iso_top_core, iso_bot_core=config.insulations.iso_bot_core,
                iso_left_core=config.insulations.iso_left_core, iso_right_core=config.insulations.iso_right_core,
                iso_primary_to_primary=config.insulations.iso_primary_to_primary,
                iso_secondary_to_secondary=config.insulations.iso_secondary_to_secondary,
                iso_primary_to_secondary=config.insulations.iso_primary_to_secondary,
                bobbin_coil_top=config.insulations.iso_top_core,
                bobbin_coil_bot=config.insulations.iso_bot_core,
                bobbin_coil_left=trial_parameters["inner_coil_insulation"],
                bobbin_coil_right=config.insulations.iso_right_core,

                # misc
                interleaving_type=femmt.CenterTappedInterleavingType.TypeC,
                primary_coil_turns=trial_parameters["primary_coil_turns"],
                winding_temperature=config.temperature)

            geo.set_insulation(insulation)
            geo.set_winding_windows([coil_window, transformer_window])

            geo.create_model(freq=target_and_fixed_parameters.fundamental_frequency, pre_visualize_geometry=flag_pre_visualize)

            center_tapped_study_excitation = geo.center_tapped_pre_study(plot_waveforms=False,
                                                                         time_current_vectors=[[target_and_fixed_parameters.time_extracted_vec, target_and_fixed_parameters.current_extracted_1_vec],
                                                                                               [target_and_fixed_parameters.time_extracted_vec, target_and_fixed_parameters.current_extracted_2_vec]])

            geo.stacked_core_center_tapped_study(center_tapped_study_excitation=center_tapped_study_excitation)

        @staticmethod
        def load_trial_dict(trial, target_and_fixed_parameters: StoTargetAndFixedParameters):
            # copy result files to result-file folder
            source_json_file = os.path.join(
                target_and_fixed_parameters.working_directories.fem_working_directory, "results",
                "log_electro_magnetic.json")
            destination_json_file = os.path.join(
                target_and_fixed_parameters.working_directories.fem_simulation_results_directory,
                f'case_{trial.number}.json')

            shutil.copy(source_json_file, destination_json_file)

            # read result-log
            with open(source_json_file, "r") as fd:
                loaded_data_dict = json.loads(fd.read())
            return loaded_data_dict

        @staticmethod
        def my_objective4(config, loaded_data_dict):
            total_volume = loaded_data_dict["misc"]["core_2daxi_total_volume"]
            total_loss = loaded_data_dict["total_losses"]["total_losses"]

            # Get inductance values
            difference_l_h = config.l_h_target - loaded_data_dict["inductances"]["L_h"]
            difference_l_s12 = config.l_s12_target - loaded_data_dict["inductances"]["L_s12"]
            return total_volume, total_loss, abs(difference_l_h), abs(difference_l_s12)

        @staticmethod
        def objective(trial, config: StoSingleInputConfig, target_and_fixed_parameters: StoTargetAndFixedParameters, my_objective: str):
            """
            Objective for optuna optimization.

            :param my_objective:
            :param trial: optuna trail objective. Used by optuna
            :param config: simulation configuration file
            :type config: StoSingleInputConfig
            :param target_and_fixed_parameters: contains pre-calculated values
            :type target_and_fixed_parameters: StoTargetAndFixedParameters

            """

            trial_parameters = femmt.optimization.StackedTransformerOptimization.FemSimulation.model_pre_parameterization(trial=trial, config=config)

            if trial_parameters is not None:
                try:
                    femmt.optimization.StackedTransformerOptimization.FemSimulation.fem_simulation(trial, trial_parameters, config, target_and_fixed_parameters)

                    loaded_data_dict = femmt.optimization.StackedTransformerOptimization.FemSimulation.load_trial_dict(trial, target_and_fixed_parameters)

                    if my_objective == "my_objective4":
                        return femmt.optimization.StackedTransformerOptimization.FemSimulation.my_objective4(config, loaded_data_dict)

                except Exception as e:
                    print(e)
                    return float('nan'), float('nan'), float('nan'), float('nan')  # TODO: dependend

            else:
                return float('nan'), float('nan'), float('nan'), float('nan')  # TODO: dependend

        @staticmethod
        def dofs_to_all_parameters(dof_params, config):
            """inout are the necessary parameters (degree of freedom) and all other parameters are returned"""
            all_paramters = dof_params
            all_paramters["core_material"] = all_paramters["material"]  # TODO: ranme to uniform designation
            # calculate dependend (core) parameters
            all_paramters["primary_litz_parameters"] = ff.litz_database()[all_paramters["primary_litz_wire"]]
            primary_litz_diameter = 2 * all_paramters["primary_litz_parameters"]["conductor_radii"]

            # Will always be calculated from the given parameters
            available_width = all_paramters["window_w"] - all_paramters["inner_coil_insulation"] - config.insulations.iso_right_core

            # Note: int() is used to round down.
            number_rows_coil_winding = int((all_paramters["primary_coil_turns"] * (primary_litz_diameter + config.insulations.iso_primary_to_primary) - all_paramters["inner_coil_insulation"]) / available_width) + 1
            all_paramters["window_h_top"] = config.insulations.iso_top_core + config.insulations.iso_bot_core + number_rows_coil_winding * primary_litz_diameter + (
                        number_rows_coil_winding - 1) * config.insulations.iso_primary_to_primary
            return all_paramters


        class NSGAII:

            @staticmethod
            def start_study(study_name: str, config: StoSingleInputConfig, number_trials: int, my_objective: str, my_objective_dimension: int, storage: str = None) -> None:

                if os.path.exists(f"{config.working_directory}/study_{study_name}.sqlite3"):
                    raise Exception(f"study '{study_name}' already availabe. Choose different study name.")

                # calculate the target and fixed parameters
                # and generate the folder structure inside this function
                target_and_fixed_parameters = femmt.optimization.StackedTransformerOptimization.calculate_fix_parameters(config)

                # Wrap the objective inside a lambda and call objective inside it
                func = lambda trial: femmt.StackedTransformerOptimization.FemSimulation.objective(trial, config, target_and_fixed_parameters, my_objective)

                objective_dimension = my_objective_dimension  # TODO: make it generic
                my_sampler = optuna.samplers.NSGAIISampler()
                # my_sampler=optuna.samplers.TPESampler()

                # Pass func to Optuna studies
                study_in_memory = optuna.create_study(directions=["minimize"] * objective_dimension, sampler=my_sampler)

                # set logging verbosity: https://optuna.readthedocs.io/en/stable/reference/generated/optuna.logging.set_verbosity.html#optuna.logging.set_verbosity
                # .INFO: all messages (default)
                # .WARNING: fails and warnings
                # .ERROR: only errors
                # optuna.logging.set_verbosity(optuna.logging.ERROR)

                print(f"Sampler is {study_in_memory.sampler.__class__.__name__}")
                study_in_memory.optimize(func, n_trials=number_trials, gc_after_trial=False, show_progress_bar=True)

                # in-memory calculation is shown before saving the data to database
                # fig = optuna.visualization.plot_pareto_front(study_in_memory, target_names=["volume", "losses", "target_l_h", "target_l_s"])
                # fig.show()

                # introduce study in storage, e.g. sqlite or mysql
                if storage == 'sqlite':
                    # Note: for sqlite operation, there needs to be three slashes '///' even before the path '/home/...'
                    # Means, in total there are four slashes including the path itself '////home/.../database.sqlite3'
                    storage = f"sqlite:///{config.working_directory}/study_{study_name}.sqlite3"
                elif storage == 'mysql':
                    storage = "mysql://monty@localhost/mydb",

                study_in_storage = optuna.create_study(directions=["minimize"] * objective_dimension, study_name=study_name, storage=storage)
                study_in_storage.add_trials(study_in_memory.trials)

            @staticmethod
            def proceed_study(study_name: str, config: StoSingleInputConfig, number_trials: int, my_objective: str, my_objective_dimension: int) -> None:
                """
                Proceed a study which is stored as sqlite database.

                :param study_name: Name of the study
                :type study_name: str
                :param config: Simulation configuration
                :type config: ItoSingleInputConfig
                :param number_trials: Number of trials adding to the existing study
                :type number_trials: int
                """
                target_and_fixed_parameters = femmt.optimization.StackedTransformerOptimization.calculate_fix_parameters(config)

                # Wrap the objective inside a lambda and call objective inside it
                func = lambda trial: femmt.optimization.StackedTransformerOptimization.FemSimulation.objective(trial, config, target_and_fixed_parameters, my_objective)

                objective_dimension = my_objective_dimension  # TODO: make it generic
                my_sampler = optuna.samplers.NSGAIISampler()
                # my_sampler=optuna.samplers.TPESampler()

                study_in_storage = optuna.create_study(study_name=study_name, storage=f"sqlite:///{config.working_directory}/study_{study_name}.sqlite3", load_if_exists=True)

                # Pass func to Optuna studies
                study_in_memory = optuna.create_study(directions=["minimize"] * objective_dimension, sampler=my_sampler)

                study_in_memory.add_trials(study_in_storage.trials)
                study_in_memory.optimize(func, n_trials=number_trials, show_progress_bar=True)

                study_in_storage.add_trials(study_in_memory.trials[(-number_trials - 1):-1])

            @staticmethod
            def show_study_results(study_name: str, config: StoSingleInputConfig,
                                   percent_error_difference_l_h: float = 20,
                                   percent_error_difference_l_s12: float = 20) -> None:
                """
                Show the results of a study.

                :param study_name: Name of the study
                :type study_name: str
                :param config: Integrated transformer configuration file
                :type config: ItoSingleInputConfig
                :param percent_error_difference_l_h: relative error allowed in l_h
                :type percent_error_difference_l_s12: float
                :param percent_error_difference_l_s12: relative error allowed in L_s12
                :type percent_error_difference_l_s12: float
                """
                study = optuna.create_study(study_name=study_name,
                                            storage=f"sqlite:///{config.working_directory}/study_{study_name}.sqlite3",
                                            load_if_exists=True)

                # Order: total_volume, total_loss, difference_l_h, difference_l_s

                l_h_absolute_error =  percent_error_difference_l_h / 100 * config.l_h_target
                print(f"{config.l_h_target = }")
                print(f"{l_h_absolute_error = }")

                l_s_absolute_error = percent_error_difference_l_s12 / 100 * config.l_s12_target
                print(f"{config.l_s12_target = }")
                print(f"{l_s_absolute_error = }")

                # print(study.trials[1414])

                fig = optuna.visualization.plot_pareto_front(study, targets=lambda t: (t.values[0] if -l_h_absolute_error < t.values[2] < l_h_absolute_error else None, t.values[1] if -l_s_absolute_error < t.values[3] < l_s_absolute_error else None), target_names=["volume", "loss"])
                fig.show()

            @staticmethod
            def re_simulate_single_result(study_name: str, config: StoSingleInputConfig, number):
                target_and_fixed_parameters = femmt.optimization.StackedTransformerOptimization.calculate_fix_parameters(config)

                study = optuna.create_study(study_name=study_name,
                                            storage=f"sqlite:///{config.working_directory}/study_{study_name}.sqlite3",
                                            load_if_exists=True)

                # Get the hyperparameters from the trial
                trial = study.trials[number]
                trial_hyperparameters = trial.params

                # Calculate all needed parameters from the hyperparameters
                all_trial_parameters = femmt.optimization.StackedTransformerOptimization.FemSimulation.dofs_to_all_parameters(trial_hyperparameters, config)
                femmt.optimization.StackedTransformerOptimization.FemSimulation.fem_simulation(trial, all_trial_parameters, config, target_and_fixed_parameters, flag_pre_visualize=True)


    class ThermalSimulation:
        pass

