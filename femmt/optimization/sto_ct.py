"""Stacked transformer optimization (STO)."""

# python libraries
import os
import shutil
import json
import datetime
import gc
import logging

# 3rd party libraries
import optuna
import pandas as pd
from matplotlib import pyplot as plt

# FEMMT and materialdatabase libraries
from femmt.optimization.sto_ct_dtos import *
import femmt.functions_reluctance as fr
import femmt.functions as ff
import femmt.optimization.ito_functions as itof
import femmt
import materialdatabase as mdb

logger = logging.getLogger(__name__)

class StackedTransformerCenterTappedOptimization:
    """Optimize a stacked transformer."""

    @staticmethod
    def calculate_fix_parameters(config: StoCtSingleInputConfig) -> StoCtTargetAndFixedParameters:
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
        target_and_fix_parameters = StoCtTargetAndFixedParameters(
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

    @staticmethod
    def objective(trial: optuna.Trial, config: StoCtSingleInputConfig,
                  target_and_fixed_parameters: StoCtTargetAndFixedParameters,
                  number_objectives: int, show_geometries: bool = False, process_number: int = 1):
        """Objective for optuna optimization.

        :param trial: optuna trail objective. Used by optuna
        :type trial: optuna.Trial
        :param config: simulation configuration file
        :type config: StoCtSingleInputConfig
        :param target_and_fixed_parameters: contains pre-calculated values
        :type target_and_fixed_parameters: StoCtTargetAndFixedParameters
        :param number_objectives: number of objectives to give different target output parameters
        :type number_objectives: int
        :param show_geometries: True to display the geometries
        :type show_geometries: bool
        :param process_number: process number. Important for parallel computing, each simulation works in a separate folder
        :type process_number: int
        :return: returns volume, loss and absolute error of inductances in %
            if number_objectives == 3:
                total_volume, total_loss, 100 * (abs(difference_l_h / config.l_h_target) + abs(difference_l_s12 / config.l_s12_target))
            elif number_objectives == 4:
                return total_volume, total_loss, 100 * abs(difference_l_h), 100 * abs(difference_l_s12)
        """
        # suggest core geometry
        core_inner_diameter = trial.suggest_float("core_inner_diameter", config.core_inner_diameter_min_max_list[0],
                                                  config.core_inner_diameter_min_max_list[1])
        window_w = trial.suggest_float("window_w", config.window_w_min_max_list[0], config.window_w_min_max_list[1])
        air_gap_transformer = trial.suggest_float("air_gap_transformer", 0.1e-3, 5e-3)

        # suggest secondary / tertiary inner winding radius
        iso_left_core = trial.suggest_float("iso_left_core", config.insulations.iso_left_core_min,
                                            config.insulations.iso_primary_inner_bobbin)

        primary_additional_bobbin = config.insulations.iso_primary_inner_bobbin - iso_left_core

        primary_litz_wire = trial.suggest_categorical("primary_litz_wire", config.primary_litz_wire_list)

        primary_litz_parameters = ff.litz_database()[primary_litz_wire]
        primary_litz_diameter = 2 * primary_litz_parameters["conductor_radii"]

        # Will always be calculated from the given parameters
        available_width = window_w - iso_left_core - config.insulations.iso_right_core

        # Suggestion of top window coil
        # Theoretically also 0 coil turns possible (number_rows_coil_winding must then be recalculated to avoid neg. values)
        primary_coil_turns = trial.suggest_int("primary_coil_turns", config.primary_coil_turns_min_max_list[0],
                                               config.primary_coil_turns_min_max_list[1])
        # Note: int() is used to round down.
        number_rows_coil_winding = int((primary_coil_turns * (
            primary_litz_diameter + config.insulations.iso_primary_to_primary) - iso_left_core) / available_width) + 1
        window_h_top = config.insulations.iso_top_core + config.insulations.iso_bot_core + number_rows_coil_winding * primary_litz_diameter + (
            number_rows_coil_winding - 1) * config.insulations.iso_primary_to_primary

        # Maximum coil air gap depends on the maximum window height top
        air_gap_coil = trial.suggest_float("air_gap_coil", 0.1e-3, window_h_top - 0.1e-3)

        # suggest categorical
        core_material = trial.suggest_categorical("material", config.material_list)
        foil_thickness = trial.suggest_categorical("foil_thickness", config.metal_sheet_thickness_list)
        interleaving_scheme = trial.suggest_categorical("interleaving_scheme", config.interleaving_scheme_list)
        interleaving_type = trial.suggest_categorical("interleaving_type", config.interleaving_type_list)

        try:
            if config.max_transformer_total_height is not None:
                # Maximum transformer height
                window_h_bot_max = config.max_transformer_total_height - 3 * core_inner_diameter / 4 - window_h_top
                window_h_bot_min = config.window_h_bot_min_max_list[0]
                if window_h_bot_min > window_h_bot_max:
                    print(f"{number_rows_coil_winding=}")
                    print(f"{window_h_top=}")
                    raise ValueError(f"{window_h_bot_min=} > {window_h_bot_max=}")

                window_h_bot = trial.suggest_float("window_h_bot", window_h_bot_min, window_h_bot_max)

            else:
                window_h_bot = trial.suggest_float("window_h_bot", config.window_h_bot_min_max_list[0],
                                                   config.window_h_bot_min_max_list[1])

            if show_geometries:
                verbosity = femmt.Verbosity.ToConsole
            else:
                verbosity = femmt.Verbosity.Silent

            # calculate core volume for stacked core type. Check if this volume is smaller than the given limit
            core_height = window_h_bot + window_h_top + core_inner_diameter * 3 / 4
            r_outer = fr.calculate_r_outer(core_inner_diameter, window_w)
            core_volume = np.pi * r_outer ** 2 * core_height
            if core_volume > config.max_core_volume:
                raise ValueError(f"Core volume of {core_volume} > {config.max_core_volume}.")

            working_directory_single_process = os.path.join(
                target_and_fixed_parameters.working_directories.fem_working_directory, f"process_{process_number}")

            geo = femmt.MagneticComponent(component_type=femmt.ComponentType.IntegratedTransformer,
                                          working_directory=working_directory_single_process,
                                          onelab_verbosity=verbosity, simulation_name=f"Case_{trial.number}")

            geo.update_mesh_accuracies(config.mesh_accuracy, config.mesh_accuracy, config.mesh_accuracy,
                                       config.mesh_accuracy)

            electro_magnetic_directory_single_process = os.path.join(working_directory_single_process,
                                                                     "electro_magnetic")
            strands_coefficients_folder_single_process = os.path.join(electro_magnetic_directory_single_process,
                                                                      "Strands_Coefficients")

            # Update directories for each model
            geo.file_data.update_paths(working_directory_single_process, electro_magnetic_directory_single_process,
                                       strands_coefficients_folder_single_process)

            core_dimensions = femmt.dtos.StackedCoreDimensions(core_inner_diameter=core_inner_diameter,
                                                               window_w=window_w,
                                                               window_h_top=window_h_top, window_h_bot=window_h_bot)
            core = femmt.Core(core_type=femmt.CoreType.Stacked, core_dimensions=core_dimensions,
                              material=core_material, temperature=config.temperature,
                              frequency=target_and_fixed_parameters.fundamental_frequency,
                              permeability_datasource=config.permeability_datasource,
                              permeability_datatype=config.permeability_datatype,
                              permeability_measurement_setup=config.permeability_measurement_setup,
                              permittivity_datasource=config.permittivity_datasource,
                              permittivity_datatype=config.permittivity_datatype,
                              permittivity_measurement_setup=config.permittivity_measurement_setup)

            geo.set_core(core)

            air_gaps = femmt.AirGaps(femmt.AirGapMethod.Stacked, core)
            air_gaps.add_air_gap(femmt.AirGapLegPosition.CenterLeg, air_gap_coil,
                                 stacked_position=femmt.StackedPosition.Top)
            air_gaps.add_air_gap(femmt.AirGapLegPosition.CenterLeg, air_gap_transformer,
                                 stacked_position=femmt.StackedPosition.Bot)
            geo.set_air_gaps(air_gaps)

            # set_center_tapped_windings() automatically places the condu
            insulation, coil_window, transformer_window = femmt.functions_topologies.set_center_tapped_windings(
                core=core,

                # primary litz
                primary_additional_bobbin=primary_additional_bobbin,
                primary_turns=config.n_target,
                primary_radius=primary_litz_parameters["conductor_radii"],
                primary_number_strands=primary_litz_parameters["strands_numbers"],
                primary_strand_radius=primary_litz_parameters["strand_radii"],

                # secondary foil
                secondary_parallel_turns=2,
                secondary_thickness_foil=foil_thickness,

                # insulation
                iso_top_core=config.insulations.iso_top_core,
                iso_bot_core=config.insulations.iso_bot_core,
                iso_left_core=iso_left_core,
                iso_right_core=config.insulations.iso_right_core,
                iso_primary_to_primary=config.insulations.iso_primary_to_primary,
                iso_secondary_to_secondary=config.insulations.iso_secondary_to_secondary,
                iso_primary_to_secondary=config.insulations.iso_primary_to_secondary,
                bobbin_coil_top=config.insulations.iso_top_core,
                bobbin_coil_bot=config.insulations.iso_bot_core,
                bobbin_coil_left=config.insulations.iso_primary_inner_bobbin,
                bobbin_coil_right=config.insulations.iso_right_core,
                center_foil_additional_bobbin=0e-3,
                interleaving_scheme=interleaving_scheme,

                # misc
                interleaving_type=interleaving_type,
                primary_coil_turns=primary_coil_turns,
                winding_temperature=config.temperature)

            geo.set_insulation(insulation)
            geo.set_winding_windows([coil_window, transformer_window])

            geo.create_model(freq=target_and_fixed_parameters.fundamental_frequency,
                             pre_visualize_geometry=show_geometries)

            center_tapped_study_excitation = geo.stacked_core_center_tapped_pre_study(
                time_current_vectors=[[target_and_fixed_parameters.time_extracted_vec,
                                       target_and_fixed_parameters.current_extracted_1_vec],
                                      [target_and_fixed_parameters.time_extracted_vec,
                                       target_and_fixed_parameters.current_extracted_2_vec]],
                fft_filter_value_factor=config.fft_filter_value_factor)

            geo.stacked_core_center_tapped_study(center_tapped_study_excitation,
                                                 number_primary_coil_turns=primary_coil_turns)

            # copy result files to result-file folder
            source_json_file = os.path.join(
                target_and_fixed_parameters.working_directories.fem_working_directory, f'process_{process_number}',
                "results", "log_electro_magnetic.json")
            destination_json_file = os.path.join(
                target_and_fixed_parameters.working_directories.fem_simulation_results_directory,
                f'case_{trial.number}.json')

            shutil.copy(source_json_file, destination_json_file)

            # read result-log
            with open(source_json_file, "r") as fd:
                loaded_data_dict = json.loads(fd.read())

            total_volume = loaded_data_dict["misc"]["core_2daxi_total_volume"]
            total_loss = loaded_data_dict["total_losses"]["total_losses"]
            # total_cost = loaded_data_dict["misc"]["total_cost_incl_margin"]

            # Get inductance values
            difference_l_h = config.l_h_target - geo.L_h
            difference_l_s12 = config.l_s12_target - geo.L_s12

            trial.set_user_attr("l_h", geo.L_h)
            trial.set_user_attr("l_s12", geo.L_s12)

            # TODO: Normalize on goal values here or the whole generation on min and max? for each feature inbetween 0 and 1
            # norm_total_loss, norm_difference_l_h, norm_difference_l_s12 = total_loss/10, #
            #     abs(difference_l_h/config.l_h_target), abs(difference_l_s12/config.l_s12_target)
            # return norm_total_loss, norm_difference_l_h, norm_difference_l_s12
            if number_objectives == 3:
                return total_volume, total_loss, 100 * (abs(difference_l_h / config.l_h_target) + abs(
                    difference_l_s12 / config.l_s12_target))
            elif number_objectives == 4:
                return total_volume, total_loss, 100 * abs(difference_l_h), 100 * abs(difference_l_s12)

        except Exception as e:
            print(e)
            if number_objectives == 3:
                return float('nan'), float('nan'), float('nan')
            elif number_objectives == 4:
                return float('nan'), float('nan'), float('nan'), float('nan')

    @staticmethod
    def start_proceed_study(study_name: str, config: StoCtSingleInputConfig, number_trials: int,
                            number_objectives: int = None,
                            storage: str = 'sqlite',
                            sampler=optuna.samplers.NSGAIISampler(),
                            show_geometries: bool = False,
                            ) -> None:
        """Proceed a study which is stored as sqlite database.

        :param study_name: Name of the study
        :type study_name: str
        :param config: Simulation configuration
        :type config: ItoSingleInputConfig
        :param number_trials: Number of trials adding to the existing study
        :type number_trials: int
        :param number_objectives: number of objectives, e.g. 3 or 4
        :type number_objectives: int
        :param storage: storage database, e.g. 'sqlite' or 'mysql'
        :type storage: str
        :param sampler: optuna.samplers.NSGAIISampler() or optuna.samplers.NSGAIIISampler(). Note about the brackets () !!
        :type sampler: optuna.sampler-object
        :param show_geometries: True to show the geometry of each suggestion (with valid geometry data)
        :type show_geometries: bool
        """

        def objective_directions(configure_number_objectives: int):
            """Check the number of objectives and returns the minimizing targets.

            :param configure_number_objectives: number of objectives
            :type configure_number_objectives: int
            :returns: objective targets and optimization function
            """
            if configure_number_objectives == 3:
                # Wrap the objective inside a lambda and call objective inside it
                return ["minimize", "minimize", "minimize"]
            if configure_number_objectives == 4:
                return ["minimize", "minimize", "minimize", 'minimize']
            else:
                raise ValueError("Invalid objective number.")

        if os.path.exists(f"{config.working_directory}/study_{study_name}.sqlite3"):
            logger.info("Existing study found. Proceeding.")

        target_and_fixed_parameters = femmt.optimization.StackedTransformerCenterTappedOptimization.calculate_fix_parameters(config)

        # introduce study in storage, e.g. sqlite or mysql
        if storage == 'sqlite':
            # Note: for sqlite operation, there needs to be three slashes '///' even before the path '/home/...'
            # Means, in total there are four slashes including the path itself '////home/.../database.sqlite3'
            storage = f"sqlite:///{config.working_directory}/study_{study_name}.sqlite3"
        elif storage == 'mysql':
            storage = "mysql://monty@localhost/mydb",

        # set logging verbosity: https://optuna.readthedocs.io/en/stable/reference/generated/optuna.logging.set_verbosity.html#optuna.logging.set_verbosity
        # .INFO: all messages (default)
        # .WARNING: fails and warnings
        # .ERROR: only errors
        # optuna.logging.set_verbosity(optuna.logging.ERROR)

        directions = objective_directions(number_objectives)

        func = lambda \
            trial: femmt.optimization.StackedTransformerCenterTappedOptimization.objective(trial, config, target_and_fixed_parameters,
                                                                                           number_objectives, show_geometries)

        study_in_storage = optuna.create_study(study_name=study_name,
                                               storage=storage,
                                               directions=directions,
                                               load_if_exists=True, sampler=sampler)

        study_in_memory = optuna.create_study(directions=directions, study_name=study_name, sampler=sampler)
        logger.info(f"Sampler is {study_in_memory.sampler.__class__.__name__}")
        study_in_memory.add_trials(study_in_storage.trials)
        study_in_memory.optimize(func, n_trials=number_trials, show_progress_bar=True,
                                 callbacks=[femmt.StackedTransformerCenterTappedOptimization.run_garbage_collector])

        study_in_storage.add_trials(study_in_memory.trials[-number_trials:])
        logger.info(f"Finished {number_trials} trials.")
        logger.info(f"current time: {datetime.datetime.now()}")

    @staticmethod
    def proceed_multi_core_study(study_name: str, config: StoCtSingleInputConfig, number_trials: int,
                                 number_objectives: int = None,
                                 storage: str = "mysql://monty@localhost/mydb",
                                 sampler=optuna.samplers.NSGAIISampler(),
                                 show_geometries: bool = False,
                                 process_number: int = 1,
                                 ) -> None:
        """Proceed a study which can be paralleled. It is highly recommended to use a mysql-database (or mariadb).

        :param study_name: Name of the study
        :type study_name: str
        :param config: Simulation configuration
        :type config: ItoSingleInputConfig
        :param number_trials: Number of trials adding to the existing study
        :type number_trials: int
        :param number_objectives: number of objectives, e.g. 3 or 4
        :type number_objectives: int
        :param storage: storage database, e.g. 'sqlite' or mysql-storage, e.g. "mysql://monty@localhost/mydb"
        :type storage: str
        :param sampler: optuna.samplers.NSGAIISampler() or optuna.samplers.NSGAIIISampler(). Note about the brackets () !!
        :type sampler: optuna.sampler-object
        :param show_geometries: True to show the geometry of each suggestion (with valid geometry data)
        :type show_geometries: bool
        :type process_number: number of the process, mandatory to split this up for several processes, because they use the same simulation result folder!
        :param process_number: int
        """

        def objective_directions(number_objectives: int):
            """Check if the number of objectives is correct and returns the minimizing targets.

            :param number_objectives: number of objectives
            :type number_objectives: int
            :returns: objective targets and optimization function
            """
            if number_objectives == 3:
                # Wrap the objective inside a lambda and call objective inside it
                return ["minimize", "minimize", "minimize"]
            if number_objectives == 4:
                return ["minimize", "minimize", "minimize", 'minimize']
            else:
                raise ValueError("Invalid objective number.")

        # introduce study in storage, e.g. sqlite or mysql
        if storage == 'sqlite':
            # Note: for sqlite operation, there needs to be three slashes '///' even before the path '/home/...'
            # Means, in total there are four slashes including the path itself '////home/.../database.sqlite3'
            storage = f"sqlite:///{config.working_directory}/study_{study_name}.sqlite3"
        elif storage == 'mysql':
            storage = "mysql://monty@localhost/mydb",

        target_and_fixed_parameters = femmt.optimization.StackedTransformerCenterTappedOptimization.calculate_fix_parameters(config)

        # set logging verbosity: https://optuna.readthedocs.io/en/stable/reference/generated/optuna.logging.set_verbosity.html#optuna.logging.set_verbosity
        # .INFO: all messages (default)
        # .WARNING: fails and warnings
        # .ERROR: only errors
        # optuna.logging.set_verbosity(optuna.logging.ERROR)

        directions = objective_directions(number_objectives)

        func = lambda trial: femmt.optimization.StackedTransformerCenterTappedOptimization.objective(trial, config, target_and_fixed_parameters,
                                                                                                     number_objectives, show_geometries, process_number)

        study_in_database = optuna.create_study(study_name=study_name,
                                                storage=storage,
                                                directions=directions,
                                                load_if_exists=True, sampler=sampler)

        study_in_database.optimize(func, n_trials=number_trials, show_progress_bar=True,
                                   callbacks=[femmt.StackedTransformerCenterTappedOptimization.run_garbage_collector])

    @staticmethod
    def run_garbage_collector(study: optuna.Study, _):
        """Run the garbage collector to prevent high memory consumption.

        :param study: optuna study
        :type study: optuna.Study

        as seen so far, a typical study needs 50.000 to 500.000 trials to have good results when performing
        the optimization. In the range above 200.000 trials, the RAM has a high occupancy rate.
        in case of running 10 or more cores, the garbage collector will run each 10.000 trial
        There could be an improvement to run the garbage collector on every process in the future.
        Now, it is random on which of the processes the garbage collector runs.
        Learn about the garbage collector https://docs.python.org/3/library/gc.html#gc.collect
        Every process runs its own garbage collector. So there is no difference between multiple processes
        https://stackoverflow.com/questions/23272943/how-does-garbage-collection-work-with-multiple-running-processes-threads
        """
        if len(study.trials) % 10000 == 0:
            print("Run garbage collector")
            gc.collect()

    @staticmethod
    def show_study_results(study_name: str, config: StoCtSingleInputConfig,
                           percent_error_difference_l_h: float = 20,
                           percent_error_difference_l_s12: float = 20) -> None:
        """Show the results of a study.

        A local .html file is generated under config.working_directory to store the interactive plotly plots on disk.

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
        l_h_absolute_error = percent_error_difference_l_h / 100 * config.l_h_target
        print(f"{config.l_h_target=}")
        print(f"{l_h_absolute_error=}")

        l_s_absolute_error = percent_error_difference_l_s12 / 100 * config.l_s12_target
        print(f"{config.l_s12_target=}")
        print(f"{l_s_absolute_error=}")

        fig = optuna.visualization.plot_pareto_front(study, targets=lambda t: (
            t.values[0] if -l_h_absolute_error < t.values[2] < l_h_absolute_error else None,
            t.values[1] if -l_s_absolute_error < t.values[3] < l_s_absolute_error else None), target_names=["volume in m³", "loss in W"])
        fig.update_layout(
            title=f"{study_name}: Filtering l_h_absolute_error = {l_h_absolute_error * 100} %, l_s_absolute_error = {l_s_absolute_error * 100} %.")
        fig.write_html(
            f"{config.working_directory}/{study_name}_error_lh_{l_h_absolute_error}_error_ls_{l_s_absolute_error}"
            f"_{datetime.datetime.now().isoformat(timespec='minutes')}.html")
        fig.show()

    @staticmethod
    def show_study_results3(study_name: str, config: StoCtSingleInputConfig,
                            error_difference_inductance_sum_percent, storage: str = 'sqlite') -> None:
        """Show the results of a study.

        A local .html file is generated under config.working_directory to store the interactive plotly plots on disk.

        :param study_name: Name of the study
        :type study_name: str
        :param config: Integrated transformer configuration file
        :type config: ItoSingleInputConfig
        :param error_difference_inductance_sum_percent: |err(L_s12) + err(L_h)| in %
        :type error_difference_inductance_sum_percent: float
        :param storage: storage, e.g. 'sqlite' or path to postgresql-database
        :type storage: str

        """
        if storage == 'sqlite':
            storage = f"sqlite:///{config.working_directory}/study_{study_name}.sqlite3"

        time_start = datetime.datetime.now()
        print(f"Start loading study {study_name} from database")
        study = optuna.load_study(study_name=study_name,
                                  storage=storage)

        print(f"Loaded study {study_name} contains {len(study.trials)} trials.")
        time_stop = datetime.datetime.now()
        print(f"Finished loading study {study_name} from database in time: {time_stop - time_start}")

        print(f"{error_difference_inductance_sum_percent=}")

        time_start = datetime.datetime.now()
        print("start generating Pareto front....")
        fig = optuna.visualization.plot_pareto_front(study, targets=lambda t: (
            t.values[0] if error_difference_inductance_sum_percent > t.values[2] else None,
            t.values[1] if error_difference_inductance_sum_percent > t.values[2] else None), target_names=["volume in m³", "loss in W"])
        fig.update_layout(title=f"{study_name}: Filtering {error_difference_inductance_sum_percent} % of |err(Ls_12)| + |err(L_h)|")
        time_stop = datetime.datetime.now()
        print(f"Finished generating Pareto front in time: {time_stop - time_start}")

        fig.write_html(
            f"{config.working_directory}/{study_name}_error_diff_{error_difference_inductance_sum_percent}%_"
            f"{datetime.datetime.now().isoformat(timespec='minutes')}.html")
        fig.show()

    @staticmethod
    def re_simulate_single_result(study_name: str, config: StoCtSingleInputConfig, number_trial: int,
                                  fft_filter_value_factor: float = 0.01, mesh_accuracy: float = 0.5,
                                  storage: str = "sqlite"):
        """Re-Simulates a FEM simulation from an optuna study.

        Performs a single simulation study (inductance, core loss, winding loss) and shows the geometry of
        number_trial design inside the study 'study_name'. Loads from an optuna study and is very slow.

        Note: This function does not use the fft_filter_value_factor and mesh_accuracy from the config-file.
        The values are given separate. In case of re-simulation, you may want to have more accurate results.

        :param study_name: name of the study
        :type study_name: str
        :param config: stacked transformer configuration file
        :type config: StoCtSingleInputConfig
        :param number_trial: number of trial to simulate
        :type number_trial: int
        :param fft_filter_value_factor: Factor to filter frequencies from the fft. E.g. 0.01 [default] removes all amplitudes below 1 % of
            the maximum amplitude from the result-frequency list
        :type fft_filter_value_factor: float
        :param mesh_accuracy: a mesh_accuracy of 0.5 is recommended. Do not change this parameter, except performing thousands of simulations,
            e.g. a Pareto optimization. In this case, the value can be set e.g. to 0.8
        :type mesh_accuracy: float
        :param storage: storage of the study
        :type storage: str
        """
        target_and_fixed_parameters = femmt.optimization.StackedTransformerCenterTappedOptimization.calculate_fix_parameters(config)

        if storage == "sqlite":
            storage = f"sqlite:///{config.working_directory}/study_{study_name}.sqlite3"

        time_start = datetime.datetime.now()
        print(f"Start loading study {study_name} from database")
        loaded_study = optuna.create_study(study_name=study_name,
                                           storage=storage,
                                           load_if_exists=True)

        print(f"Loaded study {study_name} contains {len(loaded_study.trials)} trials.")
        time_stop = datetime.datetime.now()
        print(f"Finished loading study {study_name} from database in time: {time_stop - time_start}")

        loaded_trial = loaded_study.trials[number_trial]
        loaded_trial_params = loaded_trial.params

        # suggest core geometry
        core_inner_diameter = loaded_trial_params["core_inner_diameter"]
        window_w = loaded_trial_params["window_w"]
        air_gap_transformer = loaded_trial_params["air_gap_transformer"]
        # inner_coil_insulation = trial_params["inner_coil_insulation"]
        iso_left_core = loaded_trial_params["iso_left_core"]

        primary_litz_wire = loaded_trial_params["primary_litz_wire"]

        primary_litz_parameters = ff.litz_database()[primary_litz_wire]
        primary_litz_diameter = 2 * primary_litz_parameters["conductor_radii"]

        # Will always be calculated from the given parameters
        available_width = window_w - iso_left_core - config.insulations.iso_right_core

        # Re-calculation of top window coil
        # Theoretically also 0 coil turns possible (number_rows_coil_winding must then be recalculated to avoid neg. values)
        primary_coil_turns = loaded_trial_params["primary_coil_turns"]
        # Note: int() is used to round down.
        number_rows_coil_winding = int((primary_coil_turns * (
            primary_litz_diameter + config.insulations.iso_primary_to_primary) - config.insulations.iso_primary_inner_bobbin) / available_width) + 1
        window_h_top = config.insulations.iso_top_core + config.insulations.iso_bot_core + number_rows_coil_winding * primary_litz_diameter + (
            number_rows_coil_winding - 1) * config.insulations.iso_primary_to_primary

        primary_additional_bobbin = config.insulations.iso_primary_inner_bobbin - iso_left_core

        # Maximum coil air gap depends on the maximum window height top
        air_gap_coil = loaded_trial_params["air_gap_coil"]

        # suggest categorical
        core_material = Material(loaded_trial_params["material"])
        foil_thickness = loaded_trial_params["foil_thickness"]

        if config.max_transformer_total_height is not None:
            # Maximum transformer height
            window_h_bot_max = config.max_transformer_total_height - 3 * core_inner_diameter / 4 - window_h_top
            window_h_bot_min = config.window_h_bot_min_max_list[0]
            if window_h_bot_min > window_h_bot_max:
                print(f"{number_rows_coil_winding=}")
                print(f"{window_h_top=}")
                raise ValueError(f"{window_h_bot_min=} > {window_h_bot_max=}")

            window_h_bot = loaded_trial_params["window_h_bot"]

        else:
            window_h_bot = loaded_trial_params["window_h_bot"]

        geo = femmt.MagneticComponent(component_type=femmt.ComponentType.IntegratedTransformer,
                                      working_directory=target_and_fixed_parameters.working_directories.fem_working_directory,
                                      onelab_verbosity=femmt.Verbosity.Silent,
                                      simulation_name=f"Single_Case_{loaded_trial._trial_id - 1}")
        # Note: The _trial_id starts counting from 1, while the normal cases count from zero. So a correction needs to be made

        geo.update_mesh_accuracies(mesh_accuracy, mesh_accuracy, mesh_accuracy, mesh_accuracy)

        core_dimensions = femmt.dtos.StackedCoreDimensions(core_inner_diameter=core_inner_diameter, window_w=window_w,
                                                           window_h_top=window_h_top, window_h_bot=window_h_bot)

        core = femmt.Core(core_type=femmt.CoreType.Stacked, core_dimensions=core_dimensions,
                          material=core_material, temperature=config.temperature,
                          frequency=target_and_fixed_parameters.fundamental_frequency,
                          permeability_datasource=config.permeability_datasource,
                          permeability_datatype=config.permeability_datatype,
                          permeability_measurement_setup=config.permeability_measurement_setup,
                          permittivity_datasource=config.permittivity_datasource,
                          permittivity_datatype=config.permittivity_datatype,
                          permittivity_measurement_setup=config.permittivity_measurement_setup)

        geo.set_core(core)

        air_gaps = femmt.AirGaps(femmt.AirGapMethod.Stacked, core)
        air_gaps.add_air_gap(femmt.AirGapLegPosition.CenterLeg, air_gap_coil,
                             stacked_position=femmt.StackedPosition.Top)
        air_gaps.add_air_gap(femmt.AirGapLegPosition.CenterLeg, air_gap_transformer,
                             stacked_position=femmt.StackedPosition.Bot)
        geo.set_air_gaps(air_gaps)

        # set_center_tapped_windings() automatically places the condu
        insulation, coil_window, transformer_window = femmt.functions_topologies.set_center_tapped_windings(
            core=core,

            # primary litz
            primary_additional_bobbin=primary_additional_bobbin,
            primary_turns=config.n_target,
            primary_radius=primary_litz_parameters["conductor_radii"],
            primary_number_strands=primary_litz_parameters["strands_numbers"],
            primary_strand_radius=primary_litz_parameters["strand_radii"],

            # secondary foil
            secondary_parallel_turns=2,
            secondary_thickness_foil=foil_thickness,

            # insulation
            iso_top_core=config.insulations.iso_top_core, iso_bot_core=config.insulations.iso_bot_core,
            iso_left_core=iso_left_core, iso_right_core=config.insulations.iso_right_core,
            iso_primary_to_primary=config.insulations.iso_primary_to_primary,
            iso_secondary_to_secondary=config.insulations.iso_secondary_to_secondary,
            iso_primary_to_secondary=config.insulations.iso_primary_to_secondary,
            bobbin_coil_top=config.insulations.iso_top_core,
            bobbin_coil_bot=config.insulations.iso_bot_core,
            bobbin_coil_left=config.insulations.iso_primary_inner_bobbin,
            bobbin_coil_right=config.insulations.iso_right_core,
            center_foil_additional_bobbin=0e-3,
            interleaving_scheme=InterleavingSchemesFoilLitz(loaded_trial_params['interleaving_scheme']),

            # misc
            interleaving_type=CenterTappedInterleavingType(loaded_trial_params['interleaving_type']),
            primary_coil_turns=primary_coil_turns,
            winding_temperature=config.temperature)

        geo.set_insulation(insulation)
        geo.set_winding_windows([coil_window, transformer_window])

        geo.create_model(freq=target_and_fixed_parameters.fundamental_frequency, pre_visualize_geometry=True)

        center_tapped_study_excitation = geo.stacked_core_center_tapped_pre_study(
            time_current_vectors=[[target_and_fixed_parameters.time_extracted_vec,
                                   target_and_fixed_parameters.current_extracted_1_vec],
                                  [target_and_fixed_parameters.time_extracted_vec,
                                   target_and_fixed_parameters.current_extracted_2_vec]],
            fft_filter_value_factor=fft_filter_value_factor)

        geo.stacked_core_center_tapped_study(center_tapped_study_excitation,
                                             number_primary_coil_turns=primary_coil_turns)

    @staticmethod
    def re_simulate_from_df(df: pd.DataFrame, config: StoCtSingleInputConfig, number_trial: int,
                            fft_filter_value_factor: float = 0.01, mesh_accuracy: float = 0.5,
                            show_simulation_results: bool = False):
        """Perform a single FEM simulation from a given pandas dataframe created by an optuna study.

        Performs a single simulation study (inductance, core loss, winding loss) and shows the geometry of
        number_trial design inside the study 'study_name'. Loads from a pandas dataframe and is very fast.

        Note: This function does not use the fft_filter_value_factor and mesh_accuracy from the config-file.
        The values are given separate. In case of re-simulation, you may want to have more accurate results.

        :param df: pandas dataframe with the loaded study
        :type df: pandas dataframe
        :param config: stacked transformer configuration file
        :type config: StoCtSingleInputConfig
        :param number_trial: number of trial to simulate
        :type number_trial: int
        :param fft_filter_value_factor: Factor to filter frequencies from the fft. E.g. 0.01 [default] removes all amplitudes below 1 % of
            the maximum amplitude from the result-frequency list
        :type fft_filter_value_factor: float
        :param mesh_accuracy: a mesh_accuracy of 0.5 is recommended. Do not change this parameter, except performing thousands of simulations,
            e.g. a Pareto optimization. In this case, the value can be set e.g. to 0.8
        :type mesh_accuracy: float
        :param show_simulation_results: visualize the simulation results of the FEM simulation
        :type show_simulation_results: bool
        """
        target_and_fixed_parameters = femmt.optimization.StackedTransformerCenterTappedOptimization.calculate_fix_parameters(config)

        loaded_trial_params = df.iloc[number_trial]

        # suggest core geometry
        core_inner_diameter = loaded_trial_params["params_core_inner_diameter"]
        window_w = loaded_trial_params["params_window_w"]
        air_gap_transformer = loaded_trial_params["params_air_gap_transformer"]
        # inner_coil_insulation = trial_params["inner_coil_insulation"]
        iso_left_core = loaded_trial_params["params_iso_left_core"]

        primary_litz_wire = loaded_trial_params["params_primary_litz_wire"]

        primary_litz_parameters = ff.litz_database()[primary_litz_wire]
        primary_litz_diameter = 2 * primary_litz_parameters["conductor_radii"]

        # Will always be calculated from the given parameters
        available_width = window_w - iso_left_core - config.insulations.iso_right_core

        # Re-calculation of top window coil
        # Theoretically also 0 coil turns possible (number_rows_coil_winding must then be recalculated to avoid neg. values)
        primary_coil_turns = int(loaded_trial_params["params_primary_coil_turns"])
        # Note: int() is used to round down.
        number_rows_coil_winding = int((primary_coil_turns * (
            primary_litz_diameter + config.insulations.iso_primary_to_primary) - config.insulations.iso_primary_inner_bobbin) / available_width) + 1
        window_h_top = config.insulations.iso_top_core + config.insulations.iso_bot_core + number_rows_coil_winding * primary_litz_diameter + (
            number_rows_coil_winding - 1) * config.insulations.iso_primary_to_primary

        primary_additional_bobbin = config.insulations.iso_primary_inner_bobbin - iso_left_core

        # Maximum coil air gap depends on the maximum window height top
        air_gap_coil = loaded_trial_params["params_air_gap_coil"]

        # suggest categorical
        core_material = Material(loaded_trial_params["params_material"])
        foil_thickness = loaded_trial_params["params_foil_thickness"]

        if config.max_transformer_total_height is not None:
            # Maximum transformer height
            window_h_bot_max = config.max_transformer_total_height - 3 * core_inner_diameter / 4 - window_h_top
            window_h_bot_min = config.window_h_bot_min_max_list[0]
            if window_h_bot_min > window_h_bot_max:
                print(f"{number_rows_coil_winding=}")
                print(f"{window_h_top=}")
                raise ValueError(f"{window_h_bot_min=} > {window_h_bot_max=}")

            window_h_bot = loaded_trial_params["params_window_h_bot"]

        else:
            window_h_bot = loaded_trial_params["params_window_h_bot"]

        geo = femmt.MagneticComponent(component_type=femmt.ComponentType.IntegratedTransformer,
                                      working_directory=target_and_fixed_parameters.working_directories.fem_working_directory,
                                      onelab_verbosity=femmt.Verbosity.Silent,
                                      simulation_name=f"Single_Case_{loaded_trial_params['number']}")

        geo.update_mesh_accuracies(mesh_accuracy_air_gaps=mesh_accuracy, mesh_accuracy_core=mesh_accuracy,
                                   mesh_accuracy_window=mesh_accuracy, mesh_accuracy_conductor=mesh_accuracy)

        geo.update_mesh_accuracies(mesh_accuracy, mesh_accuracy, mesh_accuracy, mesh_accuracy)

        core_dimensions = femmt.dtos.StackedCoreDimensions(core_inner_diameter=core_inner_diameter, window_w=window_w,
                                                           window_h_top=window_h_top, window_h_bot=window_h_bot)

        core = femmt.Core(core_type=femmt.CoreType.Stacked, core_dimensions=core_dimensions,
                          material=core_material, temperature=config.temperature,
                          frequency=target_and_fixed_parameters.fundamental_frequency,
                          permeability_datasource=config.permeability_datasource,
                          permeability_datatype=config.permeability_datatype,
                          permeability_measurement_setup=config.permeability_measurement_setup,
                          permittivity_datasource=config.permittivity_datasource,
                          permittivity_datatype=config.permittivity_datatype,
                          permittivity_measurement_setup=config.permittivity_measurement_setup)

        geo.set_core(core)

        air_gaps = femmt.AirGaps(femmt.AirGapMethod.Stacked, core)
        air_gaps.add_air_gap(femmt.AirGapLegPosition.CenterLeg, air_gap_coil,
                             stacked_position=femmt.StackedPosition.Top)
        air_gaps.add_air_gap(femmt.AirGapLegPosition.CenterLeg, air_gap_transformer,
                             stacked_position=femmt.StackedPosition.Bot)
        geo.set_air_gaps(air_gaps)

        # set_center_tapped_windings() automatically places the condu
        insulation, coil_window, transformer_window = femmt.functions_topologies.set_center_tapped_windings(
            core=core,

            # primary litz
            primary_additional_bobbin=primary_additional_bobbin,
            primary_turns=config.n_target,
            primary_radius=primary_litz_parameters["conductor_radii"],
            primary_number_strands=primary_litz_parameters["strands_numbers"],
            primary_strand_radius=primary_litz_parameters["strand_radii"],

            # secondary foil
            secondary_parallel_turns=2,
            secondary_thickness_foil=foil_thickness,

            # insulation
            iso_top_core=config.insulations.iso_top_core, iso_bot_core=config.insulations.iso_bot_core,
            iso_left_core=iso_left_core, iso_right_core=config.insulations.iso_right_core,
            iso_primary_to_primary=config.insulations.iso_primary_to_primary,
            iso_secondary_to_secondary=config.insulations.iso_secondary_to_secondary,
            iso_primary_to_secondary=config.insulations.iso_primary_to_secondary,
            bobbin_coil_top=config.insulations.iso_top_core,
            bobbin_coil_bot=config.insulations.iso_bot_core,
            bobbin_coil_left=config.insulations.iso_primary_inner_bobbin,
            bobbin_coil_right=config.insulations.iso_right_core,
            center_foil_additional_bobbin=0e-3,
            interleaving_scheme=InterleavingSchemesFoilLitz(loaded_trial_params['params_interleaving_scheme']),

            # misc
            interleaving_type=CenterTappedInterleavingType(loaded_trial_params['params_interleaving_type']),
            primary_coil_turns=primary_coil_turns,
            winding_temperature=config.temperature)

        geo.set_insulation(insulation)
        geo.set_winding_windows([coil_window, transformer_window])

        geo.create_model(freq=target_and_fixed_parameters.fundamental_frequency,
                         pre_visualize_geometry=show_simulation_results)

        center_tapped_study_excitation = geo.stacked_core_center_tapped_pre_study(
            time_current_vectors=[[target_and_fixed_parameters.time_extracted_vec,
                                   target_and_fixed_parameters.current_extracted_1_vec],
                                  [target_and_fixed_parameters.time_extracted_vec,
                                   target_and_fixed_parameters.current_extracted_2_vec]],
            fft_filter_value_factor=fft_filter_value_factor)

        geo.stacked_core_center_tapped_study(center_tapped_study_excitation,
                                             number_primary_coil_turns=primary_coil_turns)

        return geo

    @staticmethod
    def save_png_from_df(df: pd.DataFrame, config: StoCtSingleInputConfig, number_trial: int):
        """
        Create the geometry of specified trials_numbers.

        Note: This function does not use the fft_filter_value_factor and mesh_accuracy from the config-file.
        The values are given separate. In case of re-simulation, you may want to have more accurate results.

        :param df: pandas dataframe with the loaded study
        :type df: pandas dataframe
        :param config: stacked transformer configuration file
        :type config: StoCtSingleInputConfig
        :param number_trial: number of trial to simulate
        :type number_trial: int
        """
        target_and_fixed_parameters = femmt.optimization.StackedTransformerCenterTappedOptimization.calculate_fix_parameters(config)

        loaded_trial_params = df.iloc[number_trial]

        # suggest core geometry
        core_inner_diameter = loaded_trial_params["params_core_inner_diameter"]
        window_w = loaded_trial_params["params_window_w"]
        air_gap_transformer = loaded_trial_params["params_air_gap_transformer"]
        # inner_coil_insulation = trial_params["inner_coil_insulation"]
        iso_left_core = loaded_trial_params["params_iso_left_core"]

        primary_litz_wire = loaded_trial_params["params_primary_litz_wire"]

        primary_litz_parameters = ff.litz_database()[primary_litz_wire]
        primary_litz_diameter = 2 * primary_litz_parameters["conductor_radii"]

        # Will always be calculated from the given parameters
        available_width = window_w - iso_left_core - config.insulations.iso_right_core

        # Re-calculation of top window coil
        # Theoretically also 0 coil turns possible (number_rows_coil_winding must then be recalculated to avoid neg. values)
        primary_coil_turns = int(loaded_trial_params["params_primary_coil_turns"])
        # Note: int() is used to round down.
        number_rows_coil_winding = int((primary_coil_turns * (
            primary_litz_diameter + config.insulations.iso_primary_to_primary) - config.insulations.iso_primary_inner_bobbin) / available_width) + 1
        window_h_top = config.insulations.iso_top_core + config.insulations.iso_bot_core + number_rows_coil_winding * primary_litz_diameter + (
            number_rows_coil_winding - 1) * config.insulations.iso_primary_to_primary

        primary_additional_bobbin = config.insulations.iso_primary_inner_bobbin - iso_left_core

        # Maximum coil air gap depends on the maximum window height top
        air_gap_coil = loaded_trial_params["params_air_gap_coil"]

        # suggest categorical
        core_material = Material(loaded_trial_params["params_material"])
        foil_thickness = loaded_trial_params["params_foil_thickness"]

        if config.max_transformer_total_height is not None:
            # Maximum transformer height
            window_h_bot_max = config.max_transformer_total_height - 3 * core_inner_diameter / 4 - window_h_top
            window_h_bot_min = config.window_h_bot_min_max_list[0]
            if window_h_bot_min > window_h_bot_max:
                print(f"{number_rows_coil_winding=}")
                print(f"{window_h_top=}")
                raise ValueError(f"{window_h_bot_min=} > {window_h_bot_max=}")

            window_h_bot = loaded_trial_params["params_window_h_bot"]

        else:
            window_h_bot = loaded_trial_params["params_window_h_bot"]

        geo = femmt.MagneticComponent(component_type=femmt.ComponentType.IntegratedTransformer,
                                      working_directory=target_and_fixed_parameters.working_directories.fem_working_directory,
                                      onelab_verbosity=femmt.Verbosity.Silent,
                                      simulation_name=f"Single_Case_{loaded_trial_params['number']}")

        core_dimensions = femmt.dtos.StackedCoreDimensions(core_inner_diameter=core_inner_diameter, window_w=window_w,
                                                           window_h_top=window_h_top, window_h_bot=window_h_bot)

        core = femmt.Core(core_type=femmt.CoreType.Stacked, core_dimensions=core_dimensions,
                          material=core_material, temperature=config.temperature,
                          frequency=target_and_fixed_parameters.fundamental_frequency,
                          permeability_datasource=config.permeability_datasource,
                          permeability_datatype=config.permeability_datatype,
                          permeability_measurement_setup=config.permeability_measurement_setup,
                          permittivity_datasource=config.permittivity_datasource,
                          permittivity_datatype=config.permittivity_datatype,
                          permittivity_measurement_setup=config.permittivity_measurement_setup)

        geo.set_core(core)

        air_gaps = femmt.AirGaps(femmt.AirGapMethod.Stacked, core)
        air_gaps.add_air_gap(femmt.AirGapLegPosition.CenterLeg, air_gap_coil,
                             stacked_position=femmt.StackedPosition.Top)
        air_gaps.add_air_gap(femmt.AirGapLegPosition.CenterLeg, air_gap_transformer,
                             stacked_position=femmt.StackedPosition.Bot)
        geo.set_air_gaps(air_gaps)

        # set_center_tapped_windings() automatically places the condu
        insulation, coil_window, transformer_window = femmt.functions_topologies.set_center_tapped_windings(
            core=core,

            # primary litz
            primary_additional_bobbin=primary_additional_bobbin,
            primary_turns=config.n_target,
            primary_radius=primary_litz_parameters["conductor_radii"],
            primary_number_strands=primary_litz_parameters["strands_numbers"],
            primary_strand_radius=primary_litz_parameters["strand_radii"],

            # secondary foil
            secondary_parallel_turns=2,
            secondary_thickness_foil=foil_thickness,

            # insulation
            iso_top_core=config.insulations.iso_top_core, iso_bot_core=config.insulations.iso_bot_core,
            iso_left_core=iso_left_core, iso_right_core=config.insulations.iso_right_core,
            iso_primary_to_primary=config.insulations.iso_primary_to_primary,
            iso_secondary_to_secondary=config.insulations.iso_secondary_to_secondary,
            iso_primary_to_secondary=config.insulations.iso_primary_to_secondary,
            bobbin_coil_top=config.insulations.iso_top_core,
            bobbin_coil_bot=config.insulations.iso_bot_core,
            bobbin_coil_left=config.insulations.iso_primary_inner_bobbin,
            bobbin_coil_right=config.insulations.iso_right_core,
            center_foil_additional_bobbin=0e-3,
            interleaving_scheme=InterleavingSchemesFoilLitz(loaded_trial_params['params_interleaving_scheme']),

            # misc
            interleaving_type=CenterTappedInterleavingType(loaded_trial_params['params_interleaving_type']),
            primary_coil_turns=primary_coil_turns,
            winding_temperature=config.temperature)

        geo.set_insulation(insulation)
        geo.set_winding_windows([coil_window, transformer_window])

        geo.create_model(freq=target_and_fixed_parameters.fundamental_frequency, save_png=True)
        shutil.copy(geo.mesh.hybrid_color_png_file, os.path.join(geo.file_data.working_directory, f'{config.working_directory}/drawings/{number_trial}.png'))

    @staticmethod
    def thermal_simulation_from_geo(geo, thermal_config: ThermalConfig, flag_insulation: bool = False,
                                    show_visual_outputs: bool = True):
        """Perform the thermal simulation for the transformer.

        :param geo: geometry object
        :type geo: femmt.Component
        :param thermal_config: thermal configuration file
        :type thermal_config: ThermalConfig
        :param flag_insulation: True to simulate insulations. As insulations are not fully supported yet,
             the insulation is set to false
        :type flag_insulation: bool
        :param show_visual_outputs: True to show visual simulation result
        :type show_visual_outputs: bool

        """
        geo.thermal_simulation(thermal_config.thermal_conductivity_dict, thermal_config.boundary_temperatures,
                               thermal_config.boundary_flags, thermal_config.case_gap_top,
                               thermal_config.case_gap_right, thermal_config.case_gap_bot,
                               show_visual_outputs, color_scheme=femmt.colors_ba_jonas,
                               colors_geometry=femmt.colors_geometry_ba_jonas, flag_insulation=flag_insulation)

    @staticmethod
    def df_plot_pareto_front(df: pd.DataFrame, sum_inductance_error_percent: float):
        """Plot an interactive Pareto diagram (losses vs. volume) to select the transformers to re-simulate.

        :param df: Dataframe, generated from an optuna study (exported by optuna)
        :type df: pd.Dataframe
        :param sum_inductance_error_percent: maximum allowed error of |error(L_s12)| + |error(L_h)| in %
        :type sum_inductance_error_percent: float
        """
        print(df.head())
        df_pareto = df[df["values_2"] < sum_inductance_error_percent]

        names = df_pareto["number"].to_numpy()
        fig, ax = plt.subplots()
        sc = plt.scatter(df_pareto["values_0"], df_pareto["values_1"], s=10)

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
        plt.title(f"|error(L_s12)| +  |error(L_h)| < {sum_inductance_error_percent} %")
        plt.grid()
        plt.show()

    @staticmethod
    def create_full_report(df: pd.DataFrame, trials_numbers: list[int], config: StoCtSingleInputConfig,
                           thermal_config: ThermalConfig,
                           current_waveforms_operating_points: list[CurrentWorkingPoint],
                           fft_filter_value_factor: float = 0.01, mesh_accuracy: float = 0.5):
        """Create for several geometries and several working points a report.

        Simulates magnetoquasistatic and thermal for all given geometries and current working points.
        Summarizes the losses and temperatures for every working point and geometry.

        :param df: Dataframe, generated from an optuna study (exported by optuna)
        :type df: pd.Dataframe
        :param trials_numbers: List of trial numbers to re-simulate
        :type trials_numbers: list[int]
        :param config: stacked transformer optimization configuration file
        :type config: StoCtSingleInputConfig
        :param thermal_config: thermal configuration file
        :type thermal_config: ThermalConfig
        :param current_waveforms_operating_points: Trial numbers in a list to re-simulate
        :type current_waveforms_operating_points: list[int]
        :param fft_filter_value_factor: Factor to filter frequencies from the fft. E.g. 0.01 [default] removes all
            amplitudes below 1 % of the maximum amplitude from the result-frequency list
        :type fft_filter_value_factor: float
        :param mesh_accuracy: a mesh_accuracy of 0.5 is recommended. Do not change this parameter, except performing thousands of simulations,
            e.g. a Pareto optimization. In this case, the value can be set e.g. to 0.8
        :type mesh_accuracy: float
        """
        report_df = pd.DataFrame()

        for trial_number in trials_numbers:
            simulation_waveforms_dict = {"trial": trial_number}

            for count_current_waveform, current_waveform in enumerate(current_waveforms_operating_points):
                # update currents in the config file
                config.time_current_1_vec = current_waveforms_operating_points[
                    count_current_waveform].time_current_1_vec
                config.time_current_2_vec = current_waveforms_operating_points[
                    count_current_waveform].time_current_2_vec

                # perform the electromagnetic simulation
                geo_sim = femmt.StackedTransformerCenterTappedOptimization.re_simulate_from_df(df, config,
                                                                                               number_trial=trial_number,
                                                                                               show_simulation_results=False,
                                                                                               fft_filter_value_factor=fft_filter_value_factor,
                                                                                               mesh_accuracy=mesh_accuracy)
                # perform the thermal simulation
                femmt.StackedTransformerCenterTappedOptimization.thermal_simulation_from_geo(geo_sim, thermal_config,
                                                                                             show_visual_outputs=False)

                electromagnetoquasistatic_result_dict = geo_sim.read_log()
                thermal_result_dict = geo_sim.read_thermal_log()

                simulation_waveform_dict = {
                    f"total_loss_{current_waveform.name}": electromagnetoquasistatic_result_dict["total_losses"][
                        "total_losses"],
                    f"max_temp_core_{current_waveform.name}": thermal_result_dict["core_parts"]["total"]["max"],
                    f"max_temp_winding_{current_waveform.name}": thermal_result_dict["windings"]["total"]["max"],
                    "volume": electromagnetoquasistatic_result_dict["misc"]["core_2daxi_total_volume"]}

                simulation_waveforms_dict.update(simulation_waveform_dict)

            simulation_df = pd.DataFrame(data=simulation_waveforms_dict, index=[0])

            report_df = pd.concat([report_df, simulation_df], axis=0, ignore_index=True)

        report_df.to_csv(f'{config.working_directory}/summary.csv')
        print(f"Report exported to {config.working_directory}/summary.csv")

    @staticmethod
    def create_pngs(df: pd.DataFrame, trials_numbers: list[int], config: StoCtSingleInputConfig):
        """
        Create the geometry of specified trials_numbers and saves a screenshot in png format for each trial.

        Storage format "trial_number.png" as follows:
        f'{config.working_directory}/drawings/{number_trial}.png')

        :param df: Dataframe, generated from an optuna study (exported by optuna)
        :type df: pd.Dataframe
        :param trials_numbers: List of trial numbers to re-simulate
        :type trials_numbers: list[int]
        :param config: stacked transformer optimization configuration file
        :type config: StoCtSingleInputConfig
        """
        if not os.path.exists(f'{config.working_directory}/drawings'):
            os.mkdir(f'{config.working_directory}/drawings')
        for trial_number in trials_numbers:
            femmt.StackedTransformerCenterTappedOptimization.save_png_from_df(df, config, number_trial=trial_number, show_simulation_results=False)

    @staticmethod
    def study_to_df(study_name: str, database_url: str) -> pd.DataFrame:
        """
        Create a Pandas dataframe from a study.

        :param study_name: name of study
        :type study_name: str
        :param database_url: url of database
        :type database_url: str
        :return: Study results as Pandas Dataframe
        :rtype: pd.DataFrame
        """
        loaded_study = optuna.create_study(study_name=study_name, storage=database_url, load_if_exists=True)
        df = loaded_study.trials_dataframe()
        df.to_csv(f'{study_name}.csv')
        return df
