"""Inductor optimization."""
# python libraries
import os
import datetime
import pickle
import logging
import shutil
import json

# 3rd party libraries
import numpy as np
import optuna
from scipy import optimize
import magnethub as mh
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import tqdm

# onw libraries
from femmt.optimization.io_dtos import (InductorOptimizationDTO, InductorOptimizationTargetAndFixedParameters, FemInput,
                                        FemOutput, ReluctanceModelInput, ReluctanceModelOutput)
import femmt.functions as ff
import femmt.functions_reluctance as fr
import femmt.optimization.ito_functions as itof
import femmt.optimization.functions_optimization as fo
import femmt as fmt
import materialdatabase as mdb

logger = logging.getLogger(__name__)

class InductorOptimization:
    """Reluctance model and FEM simulation for the inductor optimization."""

    @staticmethod
    def filter_df(df: pd.DataFrame, x: str = "values_0", y: str = "values_1", factor_min_dc_losses: float = 1.2,
                  factor_max_dc_losses: float = 10) -> pd.DataFrame:
        """
        Remove designs with too high losses compared to the minimum losses.

        :param df: pandas dataframe with study results
        :type df: pd.DataFrame
        :param x: x-value name for Pareto plot filtering
        :type x: str
        :param y: y-value name for Pareto plot filtering
        :type y: str
        :param factor_min_dc_losses: filter factor for the minimum dc losses
        :type factor_min_dc_losses: float
        :param factor_max_dc_losses: dc_max_loss = factor_max_dc_losses * min_available_dc_losses_in_pareto_front
        :type factor_max_dc_losses: float
        :returns: pandas dataframe with Pareto front near points
        :rtype: pd.DataFrame
        """
        # figure out pareto front
        # pareto_volume_list, pareto_core_hyst_list, pareto_dto_list = self.pareto_front(volume_list, core_hyst_loss_list, valid_design_list)

        pareto_df: pd.DataFrame = fo.pareto_front_from_df(df)

        vector_to_sort = np.array([pareto_df[x], pareto_df[y]])

        # sorting 2d array by 1st row
        # https://stackoverflow.com/questions/49374253/sort-a-numpy-2d-array-by-1st-row-maintaining-columns
        sorted_vector = vector_to_sort[:, vector_to_sort[0].argsort()]
        x_pareto_vec = sorted_vector[0]
        y_pareto_vec = sorted_vector[1]

        total_losses_list = df[y][~pd.isnull(df[y])].to_numpy()

        min_total_dc_losses = total_losses_list[np.argmin(total_losses_list)]
        loss_offset = factor_min_dc_losses * min_total_dc_losses

        ref_loss_max = np.interp(df[x], x_pareto_vec, y_pareto_vec) + loss_offset
        # clip losses to a maximum of the minimum losses
        ref_loss_max = np.clip(ref_loss_max, a_min=-1, a_max=factor_max_dc_losses * min_total_dc_losses)

        pareto_df_offset = df[df[y] < ref_loss_max]

        return pareto_df_offset

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
            time_extracted, current_extracted_vec = fr.time_vec_current_vec_from_time_current_vec(
                config.time_current_vec)
            fundamental_frequency = 1 / time_extracted[-1]

            i_rms = fr.i_rms(config.time_current_vec)

            i_peak, = fr.max_value_from_value_vec(current_extracted_vec)

            (fft_frequencies, fft_amplitudes, fft_phases) = ff.fft(
                period_vector_t_i=config.time_current_vec, sample_factor=1000, plot='no', mode='time', filter_type='factor', filter_value_factor=0.03)

            # material properties
            material_db = mdb.Data()

            material_mu_r_abs_list = []
            magnet_model_list = []
            for material_name in config.material_name_list:
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
            working_directories = itof.set_up_folder_structure(config.inductor_optimization_directory)

            # finalize data to dto
            target_and_fix_parameters = InductorOptimizationTargetAndFixedParameters(
                i_rms=i_rms,
                i_peak=i_peak,
                time_extracted_vec=time_extracted,
                current_extracted_vec=current_extracted_vec,
                material_name_list=config.material_name_list,
                material_mu_r_abs_list=material_mu_r_abs_list,
                magnet_hub_model_list=magnet_model_list,
                fundamental_frequency=fundamental_frequency,
                working_directories=working_directories,
                fft_frequency_list=fft_frequencies,
                fft_amplitude_list=fft_amplitudes,
                fft_phases_list=fft_phases
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
            # Give some core names in core_name_list to use core geometries from the core database.
            if config.core_name_list is not None:
                # using fixed core sizes from the database with flexible height.
                core_name = trial.suggest_categorical("core_name", config.core_name_list)
                core = ff.core_database()[core_name]
                core_inner_diameter = core["core_inner_diameter"]
                window_w = core["window_w"]
                window_h = trial.suggest_float("window_h", 0.3 * core["window_h"], core["window_h"])

            else:
                # using arbitrary core sizes
                core_inner_diameter = trial.suggest_float("core_inner_diameter",
                                                          config.core_inner_diameter_min_max_list[0], config.core_inner_diameter_min_max_list[1])
                window_w = trial.suggest_float("window_w", config.window_w_min_max_list[0], config.window_w_min_max_list[1])
                window_h = trial.suggest_float("window_h", config.window_h_min_max_list[0], config.window_h_min_max_list[1])

            litz_wire_name = trial.suggest_categorical("litz_wire_name", config.litz_wire_name_list)
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
                logger.debug("Max. number of turns per window < 1")
                return float('nan'), float('nan')

            turns = trial.suggest_int('turns', 1, max_turns)

            material_name = trial.suggest_categorical('material_name', config.material_name_list)
            for count, material_mu_r_abs_value in enumerate(target_and_fixed_parameters.material_mu_r_abs_list):
                if target_and_fixed_parameters.material_name_list[count] == material_name:
                    material_mu_r_abs = material_mu_r_abs_value
                    magnet_material_model = target_and_fixed_parameters.magnet_hub_model_list[count]

            reluctance_model_input = ReluctanceModelInput(
                core_inner_diameter=core_inner_diameter,
                window_w=window_w,
                window_h=window_h,
                turns=turns,
                target_inductance=config.target_inductance,
                litz_wire_name=litz_wire_name,
                litz_wire_diameter=litz_wire_diameter,

                insulations=config.insulations,
                material_mu_r_abs=material_mu_r_abs,
                magnet_material_model=magnet_material_model,

                temperature=config.temperature,
                current_extracted_vec=target_and_fixed_parameters.current_extracted_vec,
                fundamental_frequency=target_and_fixed_parameters.fundamental_frequency,
                fft_frequency_list=target_and_fixed_parameters.fft_frequency_list,
                fft_amplitude_list=target_and_fixed_parameters.fft_amplitude_list,
                fft_phases_list=target_and_fixed_parameters.fft_phases_list
            )
            try:
                reluctance_output: ReluctanceModelOutput = InductorOptimization.ReluctanceModel.single_reluctance_model_simulation(reluctance_model_input)
            except ValueError as e:
                logger.debug("bot air gap: No fitting air gap length")
                return float('nan'), float('nan')

            trial.set_user_attr('p_winding', reluctance_output.p_winding)
            trial.set_user_attr('p_hyst', reluctance_output.p_hyst)
            trial.set_user_attr('l_air_gap', reluctance_output.l_air_gap)
            trial.set_user_attr('core_inner_diameter', core_inner_diameter)
            trial.set_user_attr('window_h', window_h)
            trial.set_user_attr('window_w', window_w)
            trial.set_user_attr('flux_density_peak', reluctance_output.flux_density_peak)

            return reluctance_output.volume, reluctance_output.p_loss_total

        @staticmethod
        def single_reluctance_model_simulation(reluctance_input: ReluctanceModelInput) -> ReluctanceModelOutput:
            """Perform a single reluctance model simulation. E.g. during optimization or during a re-simulation of a single operating point.

            :param reluctance_input: Input parameters for the reluctance model simulation.
            :type reluctance_input: ReluctanceModelInput
            :return: Output parameters of the reluctance model simulation
            :rtype: ReluctanceModelOutput
            """
            target_total_reluctance = reluctance_input.turns ** 2 / reluctance_input.target_inductance

            r_core_inner = fr.r_core_round(reluctance_input.core_inner_diameter, reluctance_input.window_h,
                                           reluctance_input.material_mu_r_abs)
            r_core_top_bot = fr.r_core_top_bot_radiant(reluctance_input.core_inner_diameter, reluctance_input.window_w,
                                                       reluctance_input.material_mu_r_abs, reluctance_input.core_inner_diameter / 4)
            r_core = 2 * r_core_inner + 2 * r_core_top_bot

            r_air_gap_target = target_total_reluctance - r_core

            flux = reluctance_input.turns * reluctance_input.current_extracted_vec / target_total_reluctance
            core_cross_section = (reluctance_input.core_inner_diameter / 2) ** 2 * np.pi
            flux_density = flux / core_cross_section

            # Do not cross out saturation, as the genetic algorithm is missing bad results to improve its suggestions
            # if flux_density.max() > 0.7 * material_mu_r_abs.saturation_flux_density:
            #     print(f"Flux density too high (70 % of b_sat): {flux_density} T > 0.7 * {material_mu_r_abs.saturation_flux_density} T")
            #     return float('nan'), float('nan')

            # calculate air gaps to reach the target parameters
            minimum_air_gap_length = 0.01e-3
            maximum_air_gap_length = 4e-3

            l_air_gap = optimize.brentq(
                fr.r_air_gap_round_round_sct, minimum_air_gap_length, maximum_air_gap_length,
                args=(reluctance_input.core_inner_diameter, reluctance_input.window_h / 2, reluctance_input.window_h / 2, r_air_gap_target),
                full_output=True)[0]

            # p_loss calculation
            # get power loss in W/m³ and estimated H wave in A/m
            p_density, _ = reluctance_input.magnet_material_model(flux_density, reluctance_input.fundamental_frequency, reluctance_input.temperature)

            # volume calculation
            r_outer = fr.calculate_r_outer(reluctance_input.core_inner_diameter, reluctance_input.window_w)
            volume = ff.calculate_cylinder_volume(cylinder_diameter=2 * r_outer,
                                                  cylinder_height=reluctance_input.window_h + reluctance_input.core_inner_diameter / 2)

            volume_winding_window = ((reluctance_input.core_inner_diameter / 2 + reluctance_input.window_w) ** 2 * np.pi - \
                                     (reluctance_input.core_inner_diameter / 2) ** 2 * np.pi) * reluctance_input.window_h
            volume_core = volume - volume_winding_window
            area_to_heat_sink = r_outer ** 2 * np.pi
            p_core = volume_core * p_density

            # winding loss calculation
            winding_dc_resistance = fr.resistance_litz_wire(
                core_inner_diameter=reluctance_input.core_inner_diameter, window_w=reluctance_input.window_w, window_h=reluctance_input.window_h,
                turns_count=reluctance_input.turns, iso_core_top=reluctance_input.insulations.core_top,
                iso_core_bot=reluctance_input.insulations.core_bot, iso_core_left=reluctance_input.insulations.core_left,
                iso_core_right=reluctance_input.insulations.core_right,
                iso_primary_to_primary=reluctance_input.insulations.primary_to_primary, litz_wire_name=reluctance_input.litz_wire_name,
                material="Copper", scheme="vertical_first", temperature=reluctance_input.temperature)

            winding_area = reluctance_input.turns * reluctance_input.litz_wire_diameter ** 2

            p_winding = 0
            for count, fft_frequency in enumerate(reluctance_input.fft_frequency_list):
                # proximity_factor_assumption = fmt.calc_proximity_factor_air_gap(
                #     litz_wire_name=litz_wire_name, number_turns=turns, r_1=config.insulations.core_left,
                #     frequency=fft_frequency, winding_area=winding_area,
                #     litz_wire_material_name='Copper', temperature=config.temperature)

                proximity_factor_assumption = fmt.calc_proximity_factor(
                    litz_wire_name=reluctance_input.litz_wire_name, number_turns=reluctance_input.turns, window_h=reluctance_input.window_h,
                    iso_core_top=reluctance_input.insulations.core_top, iso_core_bot=reluctance_input.insulations.core_bot,
                    frequency=fft_frequency, litz_wire_material_name='Copper', temperature=reluctance_input.temperature)

                p_winding += proximity_factor_assumption * winding_dc_resistance * reluctance_input.fft_amplitude_list[count] ** 2

            p_loss = p_winding + p_core

            reluctance_model_output = ReluctanceModelOutput(
                p_loss_total=p_loss,
                volume=volume,
                area_to_heat_sink=area_to_heat_sink,
                p_winding=p_winding,
                p_hyst=p_core,
                l_air_gap=l_air_gap,
                flux_density_peak=np.max(flux_density),
            )

            return reluctance_model_output

        @staticmethod
        def start_proceed_study(config: InductorOptimizationDTO, number_trials: int | None = None,
                                target_number_trials: int | None = None, storage: str = 'sqlite',
                                sampler=optuna.samplers.NSGAIIISampler(),
                                ) -> None:
            """
            Proceed a study which is stored as sqlite database.

            :param config: Simulation configuration
            :type config: ItoSingleInputConfig
            :param number_trials: Number of trials adding to the existing study
            :type number_trials: int
            :param target_number_trials: Number of target trials for the existing study
            :type target_number_trials: int
            :param storage: storage database, e.g. 'sqlite' or 'mysql'
            :type storage: str
            :param sampler: optuna.samplers.NSGAIISampler() or optuna.samplers.NSGAIIISampler(). Note about the brackets () !!
            :type sampler: optuna.sampler-object
            """
            if os.path.exists(f"{config.inductor_optimization_directory}/{config.inductor_study_name}.sqlite3"):
                logger.info("Existing study found. Proceeding.")

            target_and_fixed_parameters = InductorOptimization.ReluctanceModel.calculate_fix_parameters(config)

            # introduce study in storage, e.g. sqlite or mysql
            if storage == 'sqlite':
                # Note: for sqlite operation, there needs to be three slashes '///' even before the path '/home/...'
                # Means, in total there are four slashes including the path itself '////home/.../database.sqlite3'
                storage = f"sqlite:///{config.inductor_optimization_directory}/{config.inductor_study_name}.sqlite3"
            elif storage == 'mysql':
                storage = "mysql://monty@localhost/mydb",

            # set logging verbosity: https://optuna.readthedocs.io/en/stable/reference/generated/optuna.logging.set_verbosity.html#optuna.logging.set_verbosity
            # .INFO: all messages (default)
            # .WARNING: fails and warnings
            # .ERROR: only errors
            optuna.logging.set_verbosity(optuna.logging.ERROR)

            # check for differences with the old configuration file
            # config_on_disk_filepath = f"{config.inductor_optimization_directory}/{config.inductor_study_name}.pkl"
            # if os.path.exists(config_on_disk_filepath):
            #     config_on_disk = InductorOptimization.ReluctanceModel.load_config(config_on_disk_filepath)
            #     difference = deepdiff.DeepDiff(config, config_on_disk, ignore_order=True, significant_digits=10)
            #     if difference:
            #         print("Configuration file has changed from previous simulation. Do you want to proceed?")
            #         print(f"Difference: {difference}")
            #         read_text = input("'1' or Enter: proceed, 'any key': abort\nYour choice: ")
            #         if read_text == str(1) or read_text == "":
            #             print("proceed...")
            #         else:
            #             print("abort...")
            #             return None

            func = lambda trial: InductorOptimization.ReluctanceModel.objective(trial, config, target_and_fixed_parameters)

            study_in_storage = optuna.create_study(study_name=config.inductor_study_name,
                                                   storage=storage,
                                                   directions=['minimize', 'minimize'],
                                                   load_if_exists=True, sampler=sampler)
            if target_number_trials is not None:
                # simulation for a given number of target trials
                if len(study_in_storage.trials) < target_number_trials:
                    study_in_memory = optuna.create_study(directions=['minimize', 'minimize'], study_name=config.inductor_study_name, sampler=sampler)
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
                study_in_memory = optuna.create_study(directions=['minimize', 'minimize'], study_name=config.inductor_study_name, sampler=sampler)
                logger.info(f"Sampler is {study_in_memory.sampler.__class__.__name__}")
                study_in_memory.add_trials(study_in_storage.trials)
                study_in_memory.optimize(func, n_trials=number_trials, show_progress_bar=True)

                study_in_storage.add_trials(study_in_memory.trials[-number_trials:])
                logger.info(f"Finished {number_trials} trials.")
                logger.info(f"current time: {datetime.datetime.now()}")
            InductorOptimization.ReluctanceModel.save_config(config)
            InductorOptimization.ReluctanceModel.show_study_results(config, show_results=False)

        @staticmethod
        def show_study_results(config: InductorOptimizationDTO, show_results: bool = False) -> None:
            """Show the results of a study.

            A local .html file is generated under config.working_directory to store the interactive plotly plots on disk.

            :param config: Integrated transformer configuration file
            :type config: ItoSingleInputConfig
            :param show_results: True to directly open the browser to view the results.
            :type show_results: bool
            """
            study = optuna.load_study(study_name=config.inductor_study_name,
                                      storage=f"sqlite:///{config.inductor_optimization_directory}/{config.inductor_study_name}.sqlite3")

            fig = optuna.visualization.plot_pareto_front(study, targets=lambda t: (t.values[0], t.values[1]), target_names=["volume in m³", "loss in W"])
            fig.update_layout(title=f"{config.inductor_study_name} <br><sup>{config.inductor_optimization_directory}</sup>")
            fig.write_html(f"{config.inductor_optimization_directory}/{config.inductor_study_name}"
                           f"_{datetime.datetime.now().isoformat(timespec='minutes')}.html")
            if show_results:
                fig.show()

        @staticmethod
        def save_config(config: InductorOptimizationDTO) -> None:
            """
            Save the configuration file as pickle file on the disk.

            :param config: configuration
            :type config: InductorOptimizationDTO
            """
            # convert config path to an absolute filepath
            config.inductor_optimization_directory = os.path.abspath(config.inductor_optimization_directory)
            os.makedirs(config.inductor_optimization_directory, exist_ok=True)
            with open(f"{config.inductor_optimization_directory}/{config.inductor_study_name}.pkl", 'wb') as output:
                pickle.dump(config, output, pickle.HIGHEST_PROTOCOL)

        @staticmethod
        def load_config(config_pickle_filepath: str) -> InductorOptimizationDTO:
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
        def study_to_df(config: InductorOptimizationDTO) -> pd.DataFrame:
            """
            Create a Pandas dataframe from a study.

            :param config: configuration
            :type config: InductorOptimizationDTO
            :return: Study results as Pandas Dataframe
            :rtype: pd.DataFrame
            """
            database_url = f'sqlite:///{config.inductor_optimization_directory}/{config.inductor_study_name}.sqlite3'
            if os.path.isfile(database_url.replace('sqlite:///', '')):
                logger.info("Existing study found.")
            else:
                raise ValueError(f"Can not find database: {database_url}")
            loaded_study = optuna.load_study(study_name=config.inductor_study_name, storage=database_url)
            df = loaded_study.trials_dataframe()
            df.to_csv(f'{config.inductor_optimization_directory}/{config.inductor_study_name}.csv')
            logger.info(f"Exported study as .csv file: {config.inductor_optimization_directory}/{config.inductor_study_name}.csv")
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

            :param df: pandas dataframe with study results
            :type df: pd.DataFrame
            :param factor_min_dc_losses: filter factor for the minimum dc losses
            :type factor_min_dc_losses: float
            :param factor_max_dc_losses: dc_max_loss = factor_max_dc_losses * min_available_dc_losses_in_pareto_front
            :type factor_max_dc_losses: float
            :returns: pandas dataframe with Pareto front near points
            :rtype: pd.DataFrame
            """
            filtered_df = InductorOptimization.filter_df(df, x="values_0", y="values_1", factor_min_dc_losses=factor_min_dc_losses,
                                                         factor_max_dc_losses=factor_max_dc_losses)
            return filtered_df

        @staticmethod
        def df_from_trial_numbers(df: pd.DataFrame, trial_number_list: list[int]) -> pd.DataFrame:
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
        def full_simulation(df_geometry: pd.DataFrame, current_waveform: list, inductor_config_filepath: str) -> tuple:
            """
            Reluctance model (hysteresis losses, winding losses) for geometries from df_geometry.

            :param df_geometry: Pandas dataframe with geometries
            :type df_geometry: pd.DataFrame
            :param current_waveform: Current waveform to simulate
            :type current_waveform: list
            :param inductor_config_filepath: Filepath of the inductor optimization configuration file
            :type inductor_config_filepath: str
            :return: volume, loss
            :rtype: tuple
            """
            for index, _ in df_geometry.iterrows():

                local_config = InductorOptimization.ReluctanceModel.load_config(inductor_config_filepath)

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
                target_and_fix_parameters = InductorOptimization.ReluctanceModel.calculate_fix_parameters(local_config)

                litz_wire = ff.litz_database()[df_geometry['params_litz_wire_name'][index]]
                litz_wire_diameter = 2 * litz_wire["conductor_radii"]

                # material properties
                material_db = mdb.Data()

                material_name = mdb.meta.Material(df_geometry['params_material_name'][index])

                small_signal_mu_real_over_f_at_T = material_db.get_datasheet_curve(material_name, mdb.DatasheetCurveType.small_signal_mu_real_over_f_at_T)
                small_signal_mu_imag_over_f_at_T = material_db.get_datasheet_curve(material_name, mdb.DatasheetCurveType.small_signal_mu_imag_over_f_at_T)

                mu_real_at_f = np.interp(10_000, small_signal_mu_real_over_f_at_T["f"], small_signal_mu_real_over_f_at_T["mu_real"])
                mu_imag_at_f = np.interp(10_000, small_signal_mu_imag_over_f_at_T["f"], small_signal_mu_imag_over_f_at_T["mu_imag"])

                material_mu_r_abs = np.abs(np.array([mu_real_at_f + 1j * mu_imag_at_f]))

                # instantiate material-specific model
                magnet_material_model: mh.loss.LossModel = mh.loss.LossModel(material=material_name, team="paderborn")

                reluctance_model_input = ReluctanceModelInput(
                    core_inner_diameter=core_inner_diameter,
                    window_w=window_w,
                    window_h=df_geometry["params_window_h"][index],
                    turns=int(df_geometry['params_turns'][index].item()),
                    target_inductance=local_config.target_inductance,
                    litz_wire_name=df_geometry['params_litz_wire_name'][index],
                    litz_wire_diameter=litz_wire_diameter,

                    insulations=local_config.insulations,
                    material_mu_r_abs=material_mu_r_abs,
                    magnet_material_model=magnet_material_model,

                    temperature=local_config.temperature,
                    current_extracted_vec=target_and_fix_parameters.current_extracted_vec,
                    fundamental_frequency=target_and_fix_parameters.fundamental_frequency,
                    fft_frequency_list=target_and_fix_parameters.fft_frequency_list,
                    fft_amplitude_list=target_and_fix_parameters.fft_amplitude_list,
                    fft_phases_list=target_and_fix_parameters.fft_phases_list
                )

                reluctance_output: ReluctanceModelOutput = InductorOptimization.ReluctanceModel.single_reluctance_model_simulation(reluctance_model_input)

                p_total = reluctance_output.p_loss_total
                return reluctance_output.volume, p_total, reluctance_output.area_to_heat_sink

    class FemSimulation:
        """Contains methods to perform FEM simulations or process their results."""

        @staticmethod
        def fem_simulations_from_reluctance_df(reluctance_df: pd.DataFrame, config: InductorOptimizationDTO,
                                               show_visual_outputs: bool = False, process_number: int = 1):
            """
            Perform FEM simulations from a given Pandas dataframe. The dataframe is from the reluctance model results.

            :param reluctance_df: Pandas dataframe containing results from the reluctance model
            :type reluctance_df: pandas.DataFrame
            :param config: Configuration for the optimization of the transformer
            :type config: InductorOptimizationDTO
            :param show_visual_outputs: True to show visual outputs like the geometry
            :type show_visual_outputs: bool
            :param process_number: Process number for parallel simulations on multiple cpu cores
            :type process_number: int
            """
            target_and_fix_parameters = InductorOptimization.ReluctanceModel.calculate_fix_parameters(config)

            fem_working_directory = target_and_fix_parameters.working_directories.fem_working_directory
            if not os.path.exists(fem_working_directory):
                os.mkdir(fem_working_directory)
            fem_working_directory = os.path.join(
                target_and_fix_parameters.working_directories.fem_working_directory, f"process_{process_number}")

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
                            window_h = reluctance_df["user_attrs_window_h"][index].item()
                        else:
                            core_inner_diameter = reluctance_df["params_core_inner_diameter"][index].item()
                            window_w = reluctance_df["params_window_w"][index].item()
                            window_h = reluctance_df["params_window_h"][index].item()
                        fem_input = FemInput(
                            simulation_name=f"case_{reluctance_df['number'][index].item()}",
                            working_directory=fem_working_directory,
                            core_inner_diameter=core_inner_diameter,
                            window_w=window_w,
                            window_h=window_h,
                            material_name=reluctance_df["params_material_name"][index],
                            temperature=config.temperature,
                            material_data_sources=config.material_data_sources,
                            insulations=config.insulations,
                            fundamental_frequency=target_and_fix_parameters.fundamental_frequency,
                            air_gap_length=reluctance_df["user_attrs_l_air_gap"][index].item(),
                            litz_wire_name=reluctance_df['params_litz_wire_name'][index],
                            turns=int(reluctance_df["params_turns"][index].item()),
                            fft_frequency_list=target_and_fix_parameters.fft_frequency_list,
                            fft_amplitude_list=target_and_fix_parameters.fft_amplitude_list,
                            fft_phases_list=target_and_fix_parameters.fft_phases_list,
                        )

                        # fem simulation here
                        fem_output = InductorOptimization.FemSimulation.single_fem_simulation(fem_input, False)

                        reluctance_df.loc[index, 'fem_inductance'] = fem_output.fem_inductance
                        reluctance_df.loc[index, 'fem_p_loss_winding'] = fem_output.fem_p_loss_winding
                        reluctance_df.loc[index, 'fem_eddy_core'] = fem_output.fem_eddy_core
                        reluctance_df.loc[index, 'fem_core'] = fem_output.fem_core_total

                        # copy result files to result-file folder
                        source_json_file = os.path.join(
                            target_and_fix_parameters.working_directories.fem_working_directory, f'process_{process_number}',
                            "results", "log_electro_magnetic.json")

                        shutil.copy(source_json_file, destination_json_file)

                    except Exception as e:
                        print(e)
                        reluctance_df.loc[index, 'fem_inductance'] = None
                        reluctance_df.loc[index, 'fem_p_loss_winding'] = None
                        reluctance_df.loc[index, 'fem_eddy_core'] = None
                        reluctance_df.loc[index, 'fem_core'] = None
            return reluctance_df

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

            # add new columns to the dataframe, init values with None
            reluctance_df['fem_inductance'] = None
            reluctance_df['fem_p_loss_winding'] = None
            reluctance_df['fem_eddy_core'] = None
            reluctance_df['fem_core'] = None

            for file in files_in_folder:
                with open(os.path.join(fem_results_folder_path, file), 'r') as log_file:
                    scanned_log_dict = json.load(log_file)

                index = int(file.replace("case_", "").replace(".json", ""))

                # only write FEM simulation log results to df in case of reluctance number simulation is already present there
                if index in reluctance_df["number"]:
                    reluctance_df.at[index, 'fem_inductance'] = scanned_log_dict['single_sweeps'][0]['winding1']['flux_over_current'][0]
                    reluctance_df.at[index, 'fem_p_loss_winding'] = scanned_log_dict['total_losses']['winding1']['total']
                    reluctance_df.at[index, 'fem_eddy_core'] = scanned_log_dict['total_losses']['eddy_core']
                    reluctance_df.at[index, 'fem_core'] = scanned_log_dict['total_losses']['core']

            # final loss calculation
            reluctance_df["combined_losses"] = reluctance_df["fem_eddy_core"] + reluctance_df["fem_p_loss_winding"] + reluctance_df["user_attrs_p_hyst"]

            return reluctance_df

        @staticmethod
        def fem_vs_reluctance_pareto(df: pd.DataFrame) -> None:
            """
            Plot the FEM simulation results vs. the reluctance model results.

            :param df: Pandas Dataframe containing reluctance model results and FEM simulation results
            :type df: pd.DataFrame
            """
            fig, ax = plt.subplots()
            legend_list = []
            plt.legend(handles=legend_list)
            plt.scatter(df["values_0"], df["values_1"], s=10, label='Reluctance Model')  # c=color_array
            df["fem_total_loss"] = df["fem_core"] + df["fem_p_loss_winding"]
            plt.scatter(df["values_0"], df["fem_total_loss"], s=10, label='FEM simulation')  # c=color_array
            plt.scatter(df["values_0"], df["combined_losses"], s=10, label="combined_losses")
            plt.legend()
            plt.xlabel('Volume in m³')
            plt.ylabel('Losses in W')
            plt.grid()
            plt.show()

        @staticmethod
        def fem_vs_reluctance(df: pd.DataFrame, config: InductorOptimizationDTO):
            """
            Add statistics to the dataframe to compare the FEM simulation with the reluctance model results.

            New columns:
             * "winding_loss_fem_vs_reluctance" = df["fem_p_loss_winding"] / df["user_attrs_p_winding"]
             * "core_loss_fem_vs_reluctance" = df["fem_core"] / df["user_attrs_p_hyst"]
             * "inductance_fem_vs_reluctance" = df["fem_inductance"] / config.target_inductance

            :param df: Pandas Dataframe
            :type df: pd.DataFrame
            :param config: Configuration file
            :type config: InductorOptimizationDTO
            :return: Pandas dataframe including the statistics information
            :rtype: pd.DataFrame
            """
            df["winding_loss_fem_vs_reluctance"] = df["fem_p_loss_winding"] / df["user_attrs_p_winding"]
            df["core_loss_fem_vs_reluctance"] = df["fem_core"] / df["user_attrs_p_hyst"]
            df["inductance_fem_vs_reluctance"] = df["fem_inductance"] / config.target_inductance

            return df

        @staticmethod
        def single_fem_simulation(fem_input: FemInput, show_visual_outputs: bool = False) -> FemOutput:
            """
            Perform a single FEM simulation.

            For parallel simulations, use different working directories.

            :param fem_input: FEM input DTO
            :type fem_input: FemInput
            :param show_visual_outputs: True to show visual outputs
            :type show_visual_outputs: bool
            :return: FEM output DTO
            :rtype: FemOutput
            """
            # 1. chose simulation type
            geo = fmt.MagneticComponent(simulation_type=fmt.SimulationType.FreqDomain, component_type=fmt.ComponentType.Inductor,
                                        working_directory=fem_input.working_directory, onelab_verbosity=fmt.Verbosity.Silent,
                                        simulation_name=fem_input.simulation_name)

            # 2. set core parameters
            core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=fem_input.core_inner_diameter,
                                                            window_w=fem_input.window_w,
                                                            window_h=fem_input.window_h, core_h=fem_input.window_h + fem_input.core_inner_diameter / 2)

            core_material = fmt.ImportedComplexCoreMaterial(material=fem_input.material_name,
                                                            temperature=fem_input.temperature,
                                                            permeability_datasource=fmt.DataSource.TDK_MDT,
                                                            permittivity_datasource=fmt.DataSource.LEA_MTB,
                                                            mdb_verbosity=fmt.Verbosity.Silent)

            core = fmt.Core(material=core_material,
                            core_type=fmt.CoreType.Single,
                            core_dimensions=core_dimensions,
                            detailed_core_model=False)

            geo.set_core(core)

            # 3. set air gap parameters
            air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
            air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, fem_input.air_gap_length, 50)
            geo.set_air_gaps(air_gaps)

            # 4. set insulations
            insulation = fmt.Insulation(flag_insulation=True)
            insulation.add_core_insulations(fem_input.insulations.core_top, fem_input.insulations.core_bot,
                                            fem_input.insulations.core_left, fem_input.insulations.core_right)
            insulation.add_winding_insulations([[fem_input.insulations.primary_to_primary]])
            geo.set_insulation(insulation)

            # 5. create winding window and virtual winding windows (vww)
            winding_window = fmt.WindingWindow(core, insulation)
            vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

            # 6. create conductor and set parameters: use solid wires
            winding = fmt.Conductor(0, fmt.ConductorMaterial.Copper, temperature=fem_input.temperature)

            primary_litz_wire = fmt.litz_database()[fem_input.litz_wire_name]
            winding.set_litz_round_conductor(primary_litz_wire['conductor_radii'], primary_litz_wire['strands_numbers'],
                                             primary_litz_wire['strand_radii'], None, fmt.ConductorArrangement.Square)

            # 7. add conductor to vww and add winding window to MagneticComponent

            vww.set_winding(winding, fem_input.turns, None, fmt.Align.CenterOnVerticalAxis,
                            placing_strategy=fmt.ConductorDistribution.VerticalUpward_HorizontalRightward,
                            zigzag=True)
            geo.set_winding_windows([winding_window])

            # 8. create the model
            geo.create_model(freq=fem_input.fundamental_frequency, pre_visualize_geometry=show_visual_outputs, save_png=False)

            current_amplitudes = [[current] for current in fem_input.fft_amplitude_list]
            phases = [[phase] for phase in fem_input.fft_phases_list]

            geo.excitation_sweep(frequency_list=fem_input.fft_frequency_list, current_list_list=current_amplitudes,
                                 phi_deg_list_list=phases, show_last_fem_simulation=show_visual_outputs)

            result_dict = geo.read_log()

            fem_output = FemOutput(
                fem_inductance=result_dict['single_sweeps'][0]['winding1']['flux_over_current'][0],
                fem_p_loss_winding=result_dict['total_losses']['winding1']['total'],
                fem_eddy_core=result_dict['total_losses']['eddy_core'],
                fem_core_total=result_dict['total_losses']['core'],
                volume=result_dict["misc"]["core_2daxi_total_volume"]
            )
            return fem_output

        @staticmethod
        def filter_combined_loss_list_df(df: pd.DataFrame, factor_min_dc_losses: float = 1.2, factor_max_dc_losses: float = 10) -> pd.DataFrame:
            """
            Remove designs with too high losses compared to the minimum losses.

            :param df: pandas dataframe with study results
            :type df: pd.DataFrame
            :param factor_min_dc_losses: filter factor for the minimum dc losses
            :type factor_min_dc_losses: float
            :param factor_max_dc_losses: dc_max_loss = factor_max_dc_losses * min_available_dc_losses_in_pareto_front
            :type factor_max_dc_losses: float
            :returns: pandas dataframe with Pareto front near points
            :rtype: pd.DataFrame
            """
            df = df.drop(df[df["combined_losses"].isna()].index)
            filtered_df = InductorOptimization.filter_df(df, x="values_0", y="combined_losses", factor_min_dc_losses=factor_min_dc_losses,
                                                         factor_max_dc_losses=factor_max_dc_losses)

            return filtered_df

        @staticmethod
        def full_simulation(df_geometry: pd.DataFrame, current_waveform: list, inductor_config_filepath: str, process_number: int = 1,
                            print_derivations: bool = False) -> tuple:
            """
            Reluctance model (hysteresis losses) and FEM simulation (winding losses and eddy current losses) for geometries from df_geometry.

            :param df_geometry: Pandas dataframe with only one single geometries
            :type df_geometry: pd.DataFrame
            :param current_waveform: Current waveform to simulate
            :type current_waveform: list
            :param inductor_config_filepath: Filepath of the inductor optimization configuration file
            :type inductor_config_filepath: str
            :param process_number: process number to run the simulation on
            :type process_number: int
            :param print_derivations: True to print derivation from FEM simulation to reluctance model
            :type print_derivations: bool
            :return: volume, loss
            :rtype: tuple
            """
            index_number = df_geometry.first_valid_index()

            if df_geometry.shape[0] > 1:
                raise ValueError("df_geometry must have only one entry.")

            local_config = InductorOptimization.ReluctanceModel.load_config(inductor_config_filepath)

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
            target_and_fix_parameters = InductorOptimization.ReluctanceModel.calculate_fix_parameters(local_config)

            working_directory = target_and_fix_parameters.working_directories.fem_working_directory
            if not os.path.exists(working_directory):
                os.mkdir(working_directory)
            working_directory = os.path.join(
                target_and_fix_parameters.working_directories.fem_working_directory, f"process_{process_number}")

            # material properties
            material_db = mdb.Data()
            material_name = mdb.meta.Material(df_geometry['params_material_name'][index_number])
            small_signal_mu_real_over_f_at_T = material_db.get_datasheet_curve(material_name, mdb.DatasheetCurveType.small_signal_mu_real_over_f_at_T)
            small_signal_mu_imag_over_f_at_T = material_db.get_datasheet_curve(material_name, mdb.DatasheetCurveType.small_signal_mu_imag_over_f_at_T)

            mu_real_at_f = np.interp(10_000, small_signal_mu_real_over_f_at_T["f"], small_signal_mu_real_over_f_at_T["mu_real"])
            mu_imag_at_f = np.interp(10_000, small_signal_mu_imag_over_f_at_T["f"], small_signal_mu_imag_over_f_at_T["mu_imag"])

            material_mu_r_abs = np.abs(np.array([mu_real_at_f + 1j * mu_imag_at_f]))


            fem_input = FemInput(
                # general parameters
                working_directory=working_directory,
                simulation_name='xx',

                # material and geometry parameters
                material_name=material_name,
                litz_wire_name=df_geometry['params_litz_wire_name'][index_number],
                core_inner_diameter=core_inner_diameter,
                window_w=window_w,
                window_h=df_geometry["params_window_h"][index_number],
                air_gap_length=df_geometry['user_attrs_l_air_gap'][index_number],
                turns=int(df_geometry['params_turns'][index_number].item()),
                insulations=local_config.insulations,

                # data sources
                material_data_sources=local_config.material_data_sources,

                # operating point conditions
                temperature=local_config.temperature,
                fundamental_frequency=target_and_fix_parameters.fundamental_frequency,
                fft_frequency_list=target_and_fix_parameters.fft_frequency_list,
                fft_amplitude_list=target_and_fix_parameters.fft_amplitude_list,
                fft_phases_list=target_and_fix_parameters.fft_phases_list
            )

            fem_output = InductorOptimization.FemSimulation.single_fem_simulation(fem_input, False)

            litz_wire = ff.litz_database()[df_geometry['params_litz_wire_name'][index_number]]
            litz_wire_diameter = 2 * litz_wire["conductor_radii"]

            # material properties
            material_db = mdb.Data()

            # instantiate material-specific model
            magnet_material_model: mh.loss.LossModel = mh.loss.LossModel(material=material_name, team="paderborn")

            reluctance_model_input = ReluctanceModelInput(
                core_inner_diameter=core_inner_diameter,
                window_w=window_w,
                window_h=df_geometry["params_window_h"][index_number],
                turns=int(df_geometry['params_turns'][index_number].item()),
                target_inductance=local_config.target_inductance,
                litz_wire_name=df_geometry['params_litz_wire_name'][index_number],
                litz_wire_diameter=litz_wire_diameter,

                insulations=local_config.insulations,
                material_mu_r_abs=material_mu_r_abs,
                magnet_material_model=magnet_material_model,

                temperature=local_config.temperature,
                current_extracted_vec=target_and_fix_parameters.current_extracted_vec,
                fundamental_frequency=target_and_fix_parameters.fundamental_frequency,
                fft_frequency_list=target_and_fix_parameters.fft_frequency_list,
                fft_amplitude_list=target_and_fix_parameters.fft_amplitude_list,
                fft_phases_list=target_and_fix_parameters.fft_phases_list
            )

            reluctance_output: ReluctanceModelOutput = InductorOptimization.ReluctanceModel.single_reluctance_model_simulation(reluctance_model_input)

            p_total = reluctance_output.p_hyst + fem_output.fem_eddy_core + fem_output.fem_p_loss_winding

            if print_derivations:
                logger.info(f"Inductance reluctance: {local_config.target_inductance}")
                logger.info(f"Inductance FEM: {fem_output.fem_inductance}")
                logger.info(f"Inductance derivation: "
                            f"{(fem_output.fem_inductance - local_config.target_inductance) / local_config.target_inductance * 100} %")
                logger.info(f"Volume reluctance: {reluctance_output.volume}")
                logger.info(f"Volume FEM: {fem_output.volume}")
                logger.info(f"Volume derivation: {(reluctance_output.volume - fem_output.volume) / reluctance_output.volume * 100} %")
                logger.info(f"P_winding reluctance: {reluctance_output.p_winding}")
                logger.info(f"P_winding FEM: {fem_output.fem_p_loss_winding}")
                logger.info(f"P_winding derivation: {(fem_output.fem_p_loss_winding - reluctance_output.p_winding) / fem_output.fem_p_loss_winding * 100}")
                logger.info(f"P_hyst reluctance: {reluctance_output.p_hyst}")
                logger.info(f"P_hyst FEM: {fem_output.fem_core_total}")
                logger.info(f"P_hyst derivation: {(reluctance_output.p_hyst - fem_output.fem_core_total) / reluctance_output.p_hyst * 100}")

            return reluctance_output.volume, p_total, reluctance_output.area_to_heat_sink
