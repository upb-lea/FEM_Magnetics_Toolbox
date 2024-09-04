"""Inductor optimization."""
# python libraries
import os
import datetime

# 3rd party libraries
import numpy as np
import optuna
from scipy import optimize
import magnethub as mh

# onw libraries
from femmt.optimization.io_dtos import InductorOptimizationDTO, InductorOptimizationTargetAndFixedParameters
import femmt.functions as ff
import femmt.functions_reluctance as fr
import femmt.optimization.ito_functions as itof
import materialdatabase as mdb

class InductorOptimization:
    """Reluctance model and FEM simulation for the inductor optimization."""

    class ReluctanceModel:
        """Reluctance model methods for the optimization."""

        @staticmethod
        def calculate_fix_parameters(config: InductorOptimizationDTO) -> InductorOptimizationTargetAndFixedParameters:
            """Calculate fix parameters what can be derived from the input configuration.

            return values are:

                i_rms
                i_rms_2
                time_extracted_vec
                current_extracted_vec
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
            time_extracted, current_extracted_vec = fr.time_vec_current_vec_from_time_current_vec(
                config.time_current_vec)
            fundamental_frequency = 1 / time_extracted[-1]

            i_rms = fr.i_rms(config.time_current_vec)

            i_peak = fr.max_value_from_value_vec(current_extracted_vec)

            # material properties
            material_db = mdb.MaterialDatabase(is_silent=True)

            material_data_list = []
            for material_name in config.material_name_list:
                material_dto: mdb.MaterialCurve = material_db.material_data_interpolation_to_dto(material_name, fundamental_frequency, config.temperature)
                material_data_list.append(material_dto)

            # set up working directories
            working_directories = itof.set_up_folder_structure(config.working_directory)

            # finalize data to dto
            target_and_fix_parameters = InductorOptimizationTargetAndFixedParameters(
                i_rms=i_rms,
                i_peak=i_peak,
                time_extracted_vec=time_extracted,
                current_extracted_vec=current_extracted_vec,
                material_dto_curve_list=material_data_list,
                fundamental_frequency=fundamental_frequency,
                working_directories=working_directories
            )

            return target_and_fix_parameters

        @staticmethod
        def objective(trial: optuna.Trial, config: InductorOptimizationDTO, target_and_fixed_parameters: InductorOptimizationTargetAndFixedParameters):
            """
            Optuna objective function.

            :param trial: Trial for the optimization
            :type trial: optuna.Trial
            :param config: inductor optimization configuration file
            :type config: InductorOptimizationDTO
            :param target_and_fixed_parameters: target and fixed pre-calculated values for the inductor optimization
            :type target_and_fixed_parameters: InductorOptimizationTargetAndFixedParameters
            """
            trial.suggest_categorical("material_name", config.material_name_list)

            # Give some core names in core_name_list to use core geometries from the core database.
            if config.core_name_list is not None:
                core_name = trial.suggest_categorical("core_name", config.core_name_list)
                core = ff.core_database()[core_name]
                core_inner_diameter = core["core_inner_diameter"]
                window_w = core["window_w"]
                window_h = trial.suggest_float("window_h", core["window_h"] * 0.3, core["window_h"])

            else:
                core_inner_diameter = trial.suggest_float("core_inner_diameter", config.core_inner_diameter_list[0], config.core_inner_diameter_list[1])
                window_w = trial.suggest_float("window_w", config.window_w_list[0], config.window_w_list[1])
                window_h = trial.suggest_float("window_h", config.window_h_list[0], config.window_h_list[1])

            litz_wire_name = trial.suggest_categorical("litz_wire_name", config.litz_wire_list)
            litz_wire = ff.litz_database()[litz_wire_name]
            litz_wire_diameter = 2 * litz_wire["conductor_radii"]

            # calculate total 2D-axi symmetric volume of the core:
            # formula: number_turns_per_row = (available_width + primary_to_primary) / (wire_diameter + primary_to_primary)
            available_width = window_w - config.insulations.core_left - config.insulations.core_right
            possible_number_turns_per_row = int(
                (available_width + config.insulations.primary_to_primary) / (litz_wire_diameter + config.insulations.primary_to_primary))
            available_height = window_h - config.insulations.core_top - config.insulations.core_bot
            possible_number_turns_per_column = int(
                (available_height + config.insulations.primary_to_primary) / (litz_wire_diameter + config.insulations.primary_to_primary))
            max_turns = possible_number_turns_per_row * possible_number_turns_per_column
            if max_turns < 1:
                print("Max. number of turns per window < 1")
                return float('nan'), float('nan')

            turns = trial.suggest_int('turns', 1, max_turns)

            material_name = trial.suggest_categorical('material_name', config.material_name_list)
            for material_dto in target_and_fixed_parameters.material_dto_curve_list:
                if material_dto.material_name == material_name:
                    material_dto: mdb.MaterialCurve = material_dto

            target_total_reluctance = turns ** 2 / config.target_inductance

            r_core_inner = fr.r_core_round(core_inner_diameter, window_h, material_dto.material_mu_r_abs)
            r_core_top_bot = fr.r_core_top_bot_radiant(core_inner_diameter, window_w, material_dto.material_mu_r_abs, core_inner_diameter/4)
            r_core = 2 * r_core_inner + 2 * r_core_top_bot

            r_air_gap_target = target_total_reluctance - r_core

            flux = turns * target_and_fixed_parameters.current_extracted_vec / target_total_reluctance
            core_cross_section = (core_inner_diameter / 2) ** 2 * np.pi
            flux_density = flux / core_cross_section
            if flux_density.max() > 0.7 * material_dto.saturation_flux_density:
                print(f"Flux density too high (70 % of b_sat): {flux_density} T > 0.7 * {material_dto.saturation_flux_density} T")
                return float('nan'), float('nan')

            # calculate air gaps to reach the target parameters
            minimum_air_gap_length = 0.01e-3
            maximum_air_gap_length = 4e-3

            try:
                l_air_gap = optimize.brentq(
                    fr.r_air_gap_round_round_sct, minimum_air_gap_length, maximum_air_gap_length,
                    args=(core_inner_diameter, window_h / 2, window_h / 2, r_air_gap_target), full_output=True)[0]

            except ValueError:
                print("bot air gap: No fitting air gap length")
                return float('nan'), float('nan')

            # p_loss calculation
            # instantiate material-specific model
            mdl = mh.loss.LossModel(material=material_name, team="paderborn")

            # get power loss in W/m³ and estimated H wave in A/m
            p_density, _ = mdl(flux_density, target_and_fixed_parameters.fundamental_frequency, config.temperature)

            # volume calculation
            r_outer = fr.calculate_r_outer(core_inner_diameter, window_w)
            volume = ff.calculate_cylinder_volume(cylinder_diameter=2 * r_outer, cylinder_height=window_h + core_inner_diameter / 2)

            volume_winding_window = ((core_inner_diameter / 2 + window_w) ** 2 * np.pi - (core_inner_diameter / 2) ** 2 * np.pi) * window_h
            volume_core = volume - volume_winding_window
            p_core = volume_core * p_density

            # calculate winding losses using a proximity factor
            proximity_factor_assumption = 1.3
            effective_conductive_cross_section = litz_wire["strands_numbers"] * litz_wire["strand_radii"] ** 2 * np.pi
            effective_conductive_radius = np.sqrt(effective_conductive_cross_section / np.pi)
            winding_resistance = fr.resistance_solid_wire(
                core_inner_diameter, window_w, turns, effective_conductive_radius, material='Copper')
            winding_dc_loss = winding_resistance * target_and_fixed_parameters.i_rms ** 2

            p_winding = proximity_factor_assumption * winding_dc_loss

            p_loss = p_winding + p_core

            trial.set_user_attr('p_winding', p_winding)
            trial.set_user_attr('p_hyst', p_core)
            trial.set_user_attr('l_air_gap', l_air_gap)

            return p_loss, volume

        @staticmethod
        def start_proceed_study(study_name: str, config: InductorOptimizationDTO, number_trials: int,
                                storage: str = 'sqlite',
                                sampler=optuna.samplers.NSGAIIISampler(),
                                ) -> None:
            """
            Proceed a study which is stored as sqlite database.

            :param study_name: Name of the study
            :type study_name: str
            :param config: Simulation configuration
            :type config: ItoSingleInputConfig
            :param number_trials: Number of trials adding to the existing study
            :type number_trials: int
            :param storage: storage database, e.g. 'sqlite' or 'mysql'
            :type storage: str
            :param sampler: optuna.samplers.NSGAIISampler() or optuna.samplers.NSGAIIISampler(). Note about the brackets () !!
            :type sampler: optuna.sampler-object
            """
            if os.path.exists(f"{config.working_directory}/{study_name}.sqlite3"):
                print("Existing study found. Proceeding.")

            target_and_fixed_parameters = InductorOptimization.ReluctanceModel.calculate_fix_parameters(config)

            # introduce study in storage, e.g. sqlite or mysql
            if storage == 'sqlite':
                # Note: for sqlite operation, there needs to be three slashes '///' even before the path '/home/...'
                # Means, in total there are four slashes including the path itself '////home/.../database.sqlite3'
                storage = f"sqlite:///{config.working_directory}/{study_name}.sqlite3"
            elif storage == 'mysql':
                storage = "mysql://monty@localhost/mydb",

            # set logging verbosity: https://optuna.readthedocs.io/en/stable/reference/generated/optuna.logging.set_verbosity.html#optuna.logging.set_verbosity
            # .INFO: all messages (default)
            # .WARNING: fails and warnings
            # .ERROR: only errors
            # optuna.logging.set_verbosity(optuna.logging.ERROR)

            func = lambda trial: InductorOptimization.ReluctanceModel.objective(trial, config, target_and_fixed_parameters)

            study_in_storage = optuna.create_study(study_name=study_name,
                                                   storage=storage,
                                                   directions=['minimize', 'minimize'],
                                                   load_if_exists=True, sampler=sampler)

            study_in_memory = optuna.create_study(directions=['minimize', 'minimize'], study_name=study_name, sampler=sampler)
            print(f"Sampler is {study_in_memory.sampler.__class__.__name__}")
            study_in_memory.add_trials(study_in_storage.trials)
            study_in_memory.optimize(func, n_trials=number_trials, show_progress_bar=True)

            study_in_storage.add_trials(study_in_memory.trials[-number_trials:])
            print(f"Finished {number_trials} trials.")
            print(f"current time: {datetime.datetime.now()}")

        @staticmethod
        def show_study_results(study_name: str, config: InductorOptimizationDTO) -> None:
            """Show the results of a study.

            A local .html file is generated under config.working_directory to store the interactive plotly plots on disk.

            :param study_name: Name of the study
            :type study_name: str
            :param config: Integrated transformer configuration file
            :type config: ItoSingleInputConfig
            """
            study = optuna.load_study(study_name=study_name, storage=f"sqlite:///{config.working_directory}/{study_name}.sqlite3")

            fig = optuna.visualization.plot_pareto_front(study, targets=lambda t: (t.values[0], t.values[1]), target_names=["volume in m³", "loss in W"])
            fig.update_layout(title=f"{study_name}")
            fig.write_html(f"{config.working_directory}/{study_name}"
                           f"_{datetime.datetime.now().isoformat(timespec='minutes')}.html")
            fig.show()
