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
        class NSGAII:
            @staticmethod
            def objective(trial, config: StoSingleInputConfig, target_and_fixed_parameters: StoTargetAndFixedParameters):
                """
                Objective for optuna optimization.

                :param trial: optuna trail objective. Used by optuna
                :param config: simulation configuration file
                :type config: StoSingleInputConfig
                :param target_and_fixed_parameters: contains pre-calculated values
                :type target_and_fixed_parameters: StoTargetAndFixedParameters

                """
                # suggest core geometry
                core_inner_diameter = trial.suggest_float("core_inner_diameter", config.core_inner_diameter_min_max_list[0], config.core_inner_diameter_min_max_list[1])
                window_w = trial.suggest_float("window_w", config.window_w_min_max_list[0], config.window_w_min_max_list[1])
                window_h_top = trial.suggest_float("window_h_top", config.window_h_top_min_max_list[0], config.window_h_top_min_max_list[1])
                window_h_bot = trial.suggest_float("window_h_bot", config.window_h_bot_min_max_list[0], config.window_h_bot_min_max_list[1])
                air_gap_coil = trial.suggest_float("air_gap_coil", 0.1e-3, window_h_top-0.1e-3)
                air_gap_transformer = trial.suggest_float("air_gap_transformer", 0.1e-3, 5e-3)
                foil_thickness = trial.suggest_float("foil_thickness", config.metal_sheet_thickness[0], config.metal_sheet_thickness[1])

                # suggest categorical
                core_material = trial.suggest_categorical("material", config.material_list)
                primary_litz_wire = trial.suggest_categorical("primary_litz_wire", config.primary_litz_wire_list)

                primary_litz_parameters = ff.litz_database()[primary_litz_wire]

                try:
                    geo = femmt.MagneticComponent(component_type=femmt.ComponentType.IntegratedTransformer,
                                                working_directory=target_and_fixed_parameters.working_directories.fem_working_directory,
                                                  silent=False, simulation_name=f"Case_{trial.number}")

                    core_dimensions = femmt.dtos.StackedCoreDimensions(core_inner_diameter=core_inner_diameter, window_w=window_w,
                                                                     window_h_top=window_h_top, window_h_bot=window_h_bot)
                    core = femmt.Core(core_type=femmt.CoreType.Stacked, core_dimensions=core_dimensions,
                                      mu_r_abs=3500, phi_mu_deg=12, sigma=1.2,
                                      permeability_datasource=femmt.MaterialDataSource.Custom,
                                      permittivity_datasource=femmt.MaterialDataSource.Custom)
                                      #material=core_material, temperature=config.temperature, frequency=target_and_fixed_parameters.fundamental_frequency,
                                      #permeability_datasource=femmt.MaterialDataSource.ManufacturerDatasheet,
                                      #permeability_datatype=femmt.MeasurementDataType.ComplexPermeability,
                                      #permittivity_datasource=femmt.MaterialDataSource.ManufacturerDatasheet,
                                      #permittivity_datatype=femmt.MeasurementDataType.ComplexPermittivity)
                    geo.set_core(core)

                    air_gaps = femmt.AirGaps(femmt.AirGapMethod.Stacked, core)
                    air_gaps.add_air_gap(femmt.AirGapLegPosition.CenterLeg, air_gap_coil, stacked_position=femmt.StackedPosition.Top)
                    air_gaps.add_air_gap(femmt.AirGapLegPosition.CenterLeg, air_gap_transformer, stacked_position=femmt.StackedPosition.Bot)
                    geo.set_air_gaps(air_gaps)

                    # set_center_tapped_windings() automatically places the condu
                    insulation, coil_window, transformer_window = femmt.functions_topologies.set_center_tapped_windings(
                        core=core,

                        # primary litz
                        primary_additional_bobbin=1e-3,
                        primary_turns=14,
                        primary_radius=primary_litz_parameters["conductor_radii"],
                        primary_number_strands=primary_litz_parameters["strands_numbers"],
                        primary_strand_radius=primary_litz_parameters["strand_radii"],

                        # secondary foil
                        secondary_parallel_turns=2,
                        secondary_thickness_foil=foil_thickness,

                        # insulation
                        iso_top_core=config.insulations.iso_top_core, iso_bot_core=config.insulations.iso_bot_core,
                        iso_left_core=config.insulations.iso_left_core, iso_right_core=config.insulations.iso_right_core,
                        iso_primary_to_primary=config.insulations.iso_primary_to_primary,
                        iso_secondary_to_secondary=config.insulations.iso_secondary_to_secondary,
                        iso_primary_to_secondary=config.insulations.iso_primary_to_secondary,
                        bobbin_coil_top=0.5e-3,
                        bobbin_coil_bot=core_dimensions.window_h_top / 2 - 0.85e-3,
                        bobbin_coil_left=0.5e-3,
                        bobbin_coil_right=0.2e-3,

                        # misc
                        interleaving_type=femmt.CenterTappedInterleavingType.TypeC,
                        primary_coil_turns=3)

                    geo.set_insulation(insulation)
                    geo.set_winding_windows([coil_window, transformer_window])

                    geo.create_model(freq=target_and_fixed_parameters.fundamental_frequency, pre_visualize_geometry=True)

                    geo.single_simulation(freq=target_and_fixed_parameters.fundamental_frequency, current=[20, 120, 120], phi_deg=[0, 180, 180],
                                          show_fem_simulation_results=False)

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

                    total_volume = loaded_data_dict["misc"]["core_2daxi_total_volume"]
                    total_loss = loaded_data_dict["total_losses"]["total_losses"]
                    total_cost = loaded_data_dict["misc"]["total_cost_incl_margin"]

                    # Getting inductance values should be after the general loss simulation,
                    # otherwise loss simulations will have very small losses (wrong!)
                    geo.get_inductances(I0=1, op_frequency=target_and_fixed_parameters.fundamental_frequency)
                    difference_l_h = config.l_h_target - geo.L_h
                    difference_l_s12 = config.l_s12_target - geo.L_s12

                    return total_volume, total_loss, difference_l_h, difference_l_s12

                except Exception as e:
                    print(e)
                    return float('nan'), float('nan'), float('nan'), float('nan')

            @staticmethod
            def start_study(study_name: str, config: StoSingleInputConfig, number_trials: int, storage: str = None) -> None:

                # calculate the target and fixed parameters
                # and generate the folder structure inside this function
                target_and_fixed_parameters = femmt.optimization.StackedTransformerOptimization.calculate_fix_parameters(config)

                # Wrap the objective inside a lambda and call objective inside it
                func = lambda trial: femmt.StackedTransformerOptimization.FemSimulation.NSGAII.objective(trial, config,
                                                                                   target_and_fixed_parameters)

                # Pass func to Optuna studies
                study_in_memory = optuna.create_study(directions=["minimize", "minimize", "minimize", "minimize"],
                                                      # sampler=optuna.samplers.TPESampler(),
                                                      sampler=optuna.samplers.NSGAIISampler(),
                                                      )

                # set logging verbosity: https://optuna.readthedocs.io/en/stable/reference/generated/optuna.logging.set_verbosity.html#optuna.logging.set_verbosity
                # .INFO: all messages (default)
                # .WARNING: fails and warnings
                # .ERROR: only errors
                #optuna.logging.set_verbosity(optuna.logging.ERROR)

                print(f"Sampler is {study_in_memory.sampler.__class__.__name__}")
                study_in_memory.optimize(func, n_trials=number_trials, gc_after_trial=False, show_progress_bar=True)

                # in-memory calculation is shown before saving the data to database
                #fig = optuna.visualization.plot_pareto_front(study_in_memory, target_names=["volume", "losses", "target_l_h", "target_l_s"])
                #fig.show()

                # introduce study in storage, e.g. sqlite or mysql
                if storage == 'sqlite':
                    # Note: for sqlite operation, there needs to be three slashes '///' even before the path '/home/...'
                    # Means, in total there are four slashes including the path itself '////home/.../database.sqlite3'
                    storage = f"sqlite:///{config.working_directory}/study_{study_name}.sqlite3"
                elif storage == 'mysql':
                    storage = "mysql://monty@localhost/mydb",

                study_in_storage = optuna.create_study(directions=["minimize", "minimize", "minimize", 'minimize'], study_name=study_name,
                                                       storage=storage)
                study_in_storage.add_trials(study_in_memory.trials)

            @staticmethod
            def show_study_results(study_name: str, config: StoSingleInputConfig) -> None:
                """
                Show the results of a study.

                :param study_name: Name of the study
                :type study_name: str
                :param config: Integrated transformer configuration file
                :type config: ItoSingleInputConfig
                """
                study = optuna.create_study(study_name=study_name,
                                            storage=f"sqlite:///{config.working_directory}/study_{study_name}.sqlite3",
                                            load_if_exists=True)

                # Order: total_volume, total_loss, difference_l_h, difference_l_s

                percent_error_difference_l_h = 20
                l_h_absolute_error =  percent_error_difference_l_h / 100 * config.l_h_target
                print(f"{config.l_h_target = }")
                print(f"{l_h_absolute_error = }")

                percent_error_difference_l_s = 20
                l_s_absolute_error = percent_error_difference_l_s / 100 * config.l_s12_target
                print(f"{config.l_s12_target = }")
                print(f"{l_s_absolute_error = }")



                fig = optuna.visualization.plot_pareto_front(study, targets=lambda t: (t.values[0] if -l_h_absolute_error < t.values[2] < l_h_absolute_error else None, t.values[1] if -l_s_absolute_error < t.values[3] < l_s_absolute_error else None), target_names=["volume", "loss"])
                fig.show()


    class ThermalSimulation:
        pass

