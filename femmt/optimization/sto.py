"""Stacked transformer optimization."""
# python libraries
import datetime
import logging
import os
from typing import List

# own libraries
from femmt.optimization.sto_dtos import StoSingleInputConfig, StoTargetAndFixedParameters
import femmt.functions as ff
import femmt.functions_reluctance as fr
import materialdatabase as mdb
import femmt.optimization.ito_functions as itof
import femmt.optimization.functions_optimization as fo
import femmt as fmt

# 3rd party libraries
import numpy as np
import optuna
from scipy import optimize
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import magnethub as mh
import pandas as pd

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

            # phi_deg_2 = phi_deg_2 - 180

            # target inductances
            target_inductance_matrix = fr.calculate_inductance_matrix_from_ls_lh_n(config.l_s12_target,
                                                                                   config.l_h_target,
                                                                                   config.n_target)

            # material properties
            material_db = mdb.MaterialDatabase(is_silent=True)

            material_data_list = []
            for material_name in config.material_list:
                material_dto: mdb.MaterialCurve = material_db.material_data_interpolation_to_dto(material_name, fundamental_frequency, config.temperature)
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

        @staticmethod
        def objective(trial: optuna.Trial, config: StoSingleInputConfig, target_and_fixed_parameters: StoTargetAndFixedParameters):
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
                window_h_bot = trial.suggest_float("window_h_bot", 0.3 * core["window_h"], core["window_h"])

            else:
                # using arbitrary core sizes
                core_inner_diameter = trial.suggest_float("core_inner_diameter", config.core_inner_diameter_min_max_list[0], config.core_inner_diameter_min_max_list[1])
                window_w = trial.suggest_float("window_w", config.window_w_min_max_list[0], config.window_w_min_max_list[1])
                window_h_bot = trial.suggest_float('window_h_bot', config.window_h_bot_min_max_list[0], config.window_h_bot_min_max_list[1])

            n_p_top = trial.suggest_int('n_p_top', config.n_p_top_min_max_list[0], config.n_p_top_min_max_list[1])
            n_p_bot = trial.suggest_int('n_p_bot', config.n_p_bot_min_max_list[0], config.n_p_bot_min_max_list[1])
            n_s_bot_min = int(np.round(n_p_bot / config.n_target, decimals=0)) - 2
            if n_s_bot_min < 1:
                n_s_bot_min = 1
            n_s_bot_max = int(np.round(n_p_bot / config.n_target, decimals=0)) + 2
            n_s_bot = trial.suggest_int('n_s_bot', n_s_bot_min, n_s_bot_max)

            primary_litz_wire = trial.suggest_categorical('primary_litz_wire', config.primary_litz_wire_list)
            primary_litz = ff.litz_database()[primary_litz_wire]
            primary_litz_wire_diameter = 2 * primary_litz['conductor_radii']

            secondary_litz_wire = trial.suggest_categorical('secondary_litz_wire', config.secondary_litz_wire_list)
            secondary_litz = ff.litz_database()[secondary_litz_wire]
            secondary_litz_wire_diameter = 2 * secondary_litz['conductor_radii']

            material_name = trial.suggest_categorical('material_name', config.material_list)
            for material_dto in target_and_fixed_parameters.material_dto_curve_list:
                if material_dto.material_name == material_name:
                    material_dto: mdb.MaterialCurve = material_dto

            # calculate total 2D-axi symmetric volume of the core:
            # formula: number_turns_per_row = (available_width + primary_to_primary) / (wire_diameter + primary_to_primary)
            available_width_top = window_w - config.insulations.iso_window_top_core_left - config.insulations.iso_window_top_core_right
            possible_number_turns_per_row_top_window = int(
                (available_width_top + config.insulations.iso_primary_to_primary) / (primary_litz_wire_diameter + config.insulations.iso_primary_to_primary))
            number_of_rows_needed = np.ceil(n_p_top / possible_number_turns_per_row_top_window)
            needed_height_top_wo_insulation = (number_of_rows_needed * primary_litz_wire_diameter + \
                                               (number_of_rows_needed - 1) * config.insulations.iso_primary_to_primary)
            window_h_top = needed_height_top_wo_insulation + config.insulations.iso_window_top_core_top + config.insulations.iso_window_top_core_bot

            core_total_height = window_h_top + window_h_bot + core_inner_diameter * 3 / 4
            r_outer = fr.calculate_r_outer(core_inner_diameter, window_w)
            volume = ff.calculate_cylinder_volume(cylinder_diameter=2 * r_outer, cylinder_height=core_total_height)

            window_bot_available_height = window_h_bot - config.insulations.iso_window_bot_core_top - config.insulations.iso_window_bot_core_bot
            window_bot_available_width = window_w - config.insulations.iso_window_bot_core_left - config.insulations.iso_window_bot_core_right
            window_bot_available_area = window_bot_available_height * window_bot_available_width

            # turn area for a single turn is approximated as a rectangle
            window_bot_turns_area = n_p_bot * primary_litz_wire_diameter ** 2 + n_s_bot * secondary_litz_wire_diameter ** 2

            # as the window_h_top is adapted to the number of n_p_top, the top windings always fit into the top window.
            if window_bot_turns_area > window_bot_available_area:
                print("Winding window too small for too many turns.")
                return float('nan'), float('nan')

            # calculate the reluctance and flux matrix
            winding_matrix = np.array([[n_p_top, 0], [n_p_bot, n_s_bot]])

            reluctance_matrix = fr.calculate_reluctance_matrix(winding_matrix, target_and_fixed_parameters.target_inductance_matrix)

            current_matrix = np.array([target_and_fixed_parameters.current_extracted_1_vec, target_and_fixed_parameters.current_extracted_2_vec])

            flux_matrix = fr.calculate_flux_matrix(reluctance_matrix, winding_matrix, current_matrix)

            flux_top = flux_matrix[0]
            flux_bot = flux_matrix[1]
            flux_middle = flux_bot - flux_top

            core_cross_section = (core_inner_diameter / 2) ** 2 * np.pi

            flux_density_top = flux_top / core_cross_section
            flux_density_bot = flux_bot / core_cross_section
            flux_density_middle = flux_middle / core_cross_section

            # calculate the core reluctance
            core_inner_cylinder_top = fr.r_core_round(core_inner_diameter, window_h_top, material_dto.material_mu_r_abs)
            core_inner_cylinder_bot = fr.r_core_round(core_inner_diameter, window_h_bot, material_dto.material_mu_r_abs)
            core_top_bot_radiant = fr.r_core_top_bot_radiant(core_inner_diameter, window_w, material_dto.material_mu_r_abs, core_inner_diameter / 4)

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

            try:
                l_top_air_gap = optimize.brentq(
                    fr.r_air_gap_round_inf_sct, minimum_air_gap_length, maximum_air_gap_length,
                    args=(core_inner_diameter, window_h_top, r_air_gap_top_target), full_output=True)[0]
            except ValueError:
                print("top air gap: No fitting air gap length")
                return float('nan'), float('nan')
            try:
                l_bot_air_gap = optimize.brentq(
                    fr.r_air_gap_round_round_sct, minimum_air_gap_length, maximum_air_gap_length,
                    args=(core_inner_diameter, window_h_bot / 2, window_h_bot / 2, r_air_gap_bot_target), full_output=True)[0]

            except ValueError:
                print("bot air gap: No fitting air gap length")
                return float('nan'), float('nan')

            # calculate hysteresis losses from mag-net-hub
            interp_points = np.arange(0, 1024) * target_and_fixed_parameters.time_extracted_vec[-1] / 1024
            flux_density_top_interp = np.interp(interp_points, target_and_fixed_parameters.time_extracted_vec, flux_density_top)
            flux_density_bot_interp = np.interp(interp_points, target_and_fixed_parameters.time_extracted_vec, flux_density_bot)
            flux_density_middle_interp = np.interp(interp_points, target_and_fixed_parameters.time_extracted_vec, flux_density_middle)

            # instantiate material-specific model
            mdl = mh.loss.LossModel(material=material_name, team="paderborn")

            # get power loss in W/m³ and estimated H wave in A/m
            p_density_top, _ = mdl(flux_density_top_interp, target_and_fixed_parameters.fundamental_frequency, config.temperature)
            p_density_bot, _ = mdl(flux_density_bot_interp, target_and_fixed_parameters.fundamental_frequency, config.temperature)
            p_density_middle, _ = mdl(flux_density_middle_interp, target_and_fixed_parameters.fundamental_frequency, config.temperature)

            volume_core_top = (2 * ff.calculate_cylinder_volume(core_inner_diameter, window_h_top) - \
                               ff.calculate_cylinder_volume(core_inner_diameter, l_top_air_gap) + \
                               ff.calculate_cylinder_volume(2 * r_outer, core_inner_diameter / 4))
            volume_core_bot = (2 * ff.calculate_cylinder_volume(core_inner_diameter, window_h_bot) - \
                               ff.calculate_cylinder_volume(core_inner_diameter, l_bot_air_gap) + \
                               ff.calculate_cylinder_volume(2 * r_outer, core_inner_diameter / 4))
            volume_core_middle = ff.calculate_cylinder_volume(2 * r_outer, core_inner_diameter / 4)

            p_top = p_density_top * volume_core_top
            p_bot = p_density_bot * volume_core_bot
            p_middle = p_density_middle * volume_core_middle

            p_hyst = p_top + p_bot + p_middle

            # calculate winding losses using a proximity factor
            proximity_factor_assumption = 1.3
            primary_effective_conductive_cross_section = primary_litz["strands_numbers"] * primary_litz["strand_radii"] ** 2 * np.pi
            primary_effective_conductive_radius = np.sqrt(primary_effective_conductive_cross_section / np.pi)
            primary_resistance = fr.resistance_solid_wire(
                core_inner_diameter, window_w, n_p_top + n_p_bot, primary_effective_conductive_radius, material='Copper')
            primary_dc_loss = primary_resistance * target_and_fixed_parameters.i_rms_1 ** 2

            secondary_effective_conductive_cross_section = secondary_litz["strands_numbers"] * secondary_litz["strand_radii"] ** 2 * np.pi
            secondary_effective_conductive_radius = np.sqrt(secondary_effective_conductive_cross_section / np.pi)
            secondary_resistance = fr.resistance_solid_wire(
                core_inner_diameter, window_w, n_s_bot, secondary_effective_conductive_radius, material='Copper')
            secondary_dc_loss = secondary_resistance * target_and_fixed_parameters.i_rms_2 ** 2

            winding_losses = proximity_factor_assumption * (primary_dc_loss + secondary_dc_loss)

            p_loss = p_hyst + winding_losses

            # set additional attributes
            trial.set_user_attr('n_s_bot', n_s_bot)
            trial.set_user_attr('p_hyst', p_hyst)
            trial.set_user_attr('p_hyst_top', p_top)
            trial.set_user_attr('p_hyst_bot', p_bot)
            trial.set_user_attr('p_hyst_middle', p_middle)
            trial.set_user_attr('window_h_top', window_h_top)
            trial.set_user_attr('winding_losses', winding_losses)
            trial.set_user_attr('l_top_air_gap', l_top_air_gap)
            trial.set_user_attr('l_bot_air_gap', l_bot_air_gap)

            return volume, p_loss

        @staticmethod
        def start_proceed_study(study_name: str, config: StoSingleInputConfig, number_trials: int,
                                storage: str = 'sqlite',
                                sampler=optuna.samplers.NSGAIIISampler(),
                                show_geometries: bool = False,
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
            :param show_geometries: True to show the geometry of each suggestion (with valid geometry data)
            :type show_geometries: bool
            """
            if os.path.exists(f"{config.working_directory}/{study_name}.sqlite3"):
                print("Existing study found. Proceeding.")

            target_and_fixed_parameters = StackedTransformerOptimization.ReluctanceModel.calculate_fix_parameters(config)

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

            func = lambda trial: StackedTransformerOptimization.ReluctanceModel.objective(trial, config, target_and_fixed_parameters)

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
        def show_study_results(study_name: str, config: StoSingleInputConfig) -> None:
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

        @staticmethod
        def study_to_df(study_name: str, working_directory: str) -> pd.DataFrame:
            """
            Create a Pandas dataframe from a study.

            :param study_name: name of study
            :type study_name: str
            :param working_directory: folder of the located study
            :type working_directory: str
            :return: Study results as Pandas Dataframe
            :rtype: pd.DataFrame
            """
            database_url = f'sqlite:///{working_directory}/{study_name}.sqlite3'
            if os.path.isfile(database_url.replace('sqlite:///', '')):
                print("Existing study found.")
            else:
                raise ValueError(f"Can not find database: {database_url}")
            loaded_study = optuna.load_study(study_name=study_name, storage=database_url)
            df = loaded_study.trials_dataframe()
            df.to_csv(f'{working_directory}/{study_name}.csv')
            logging.info(f"Exported study as .csv file: {working_directory}/{study_name}.csv")
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
        def df_trial_numbers(df: pd.DataFrame, trial_number_list: List[int]) -> pd.DataFrame:
            """
            Generate a new dataframe from a given one, just with the trial numbers from the trial_number_list.

            :param df: input dataframe
            :type df: pandas.DataFrame
            :param trial_number_list: list of trials, e.g. [1530, 1870, 3402]
            :type trial_number_list: List[int]
            :return: dataframe with trial numbers from trial_number_list
            :rtype: pandas.DataFrame
            """
            df_trial_numbers = df.loc[df["number"].isin(trial_number_list)]

            return df_trial_numbers

    class FemSimulation:
        """Contains methods to perform FEM simulations or process their results."""

        @staticmethod
        def fem_simulations_from_reluctance_df(reluctance_df: pd.DataFrame, config: StoSingleInputConfig, show_visual_outputs: bool = False,
                                               process_number: int = 1):
            """
            Perform FEM simulations from a given Pandas dataframe. The dataframe is from the reluctance model results.

            :param reluctance_df: Pandas dataframe containing relults from the relutance model
            :type reluctance_df: pandas.DataFrame
            :param config: Configuration for the optimization of the transformer
            :type config: StoSingleInputConfig
            :param show_visual_outputs: Ture to show visual outputs like the geometry
            :type show_visual_outputs: bool
            :param process_number: Process number for parallel simulations on multiple cpu cores
            :type process_number: int
            """
            target_and_fix_parameters = StackedTransformerOptimization.ReluctanceModel.calculate_fix_parameters(config)

            working_directory = target_and_fix_parameters.working_directories.fem_working_directory
            if not os.path.exists(working_directory):
                os.mkdir(working_directory)
            working_directory = os.path.join(
                target_and_fix_parameters.working_directories.fem_working_directory, f"process_{process_number}")

            time_current_vectors = np.array([config.time_current_1_vec, config.time_current_2_vec])

            # pd.read_csv(current_waveforms_csv_file, header=0, index_col=0, delimiter=';')
            df = pd.DataFrame()

            for index, _ in reluctance_df.iterrows():

                try:
                    # 1. chose simulation type
                    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.IntegratedTransformer,
                                                working_directory=working_directory, verbosity=fmt.Verbosity.ToConsole)

                    # 2. set core parameters
                    core_dimensions = fmt.dtos.StackedCoreDimensions(core_inner_diameter=reluctance_df["params_core_inner_diameter"][index],
                                                                     window_w=reluctance_df["params_window_w"][index],
                                                                     window_h_top=reluctance_df['user_attrs_window_h_top'][index],
                                                                     window_h_bot=reluctance_df["params_window_h_bot"][index])
                    core = fmt.Core(core_type=fmt.CoreType.Stacked, core_dimensions=core_dimensions,
                                    material=reluctance_df["params_material_name"][index], temperature=config.temperature,
                                    frequency=target_and_fix_parameters.fundamental_frequency,
                                    permeability_datasource=config.permeability_datasource,
                                    permeability_datatype=config.permeability_datatype,
                                    permeability_measurement_setup=config.permeability_measurement_setup,
                                    permittivity_datasource=config.permittivity_datasource,
                                    permittivity_datatype=config.permittivity_datatype,
                                    permittivity_measurement_setup=config.permittivity_measurement_setup, mdb_verbosity=fmt.Verbosity.Silent
                                    )
                    geo.set_core(core)

                    # 3. set air gap parameters
                    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Stacked, core)
                    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, reluctance_df["user_attrs_l_top_air_gap"][index],
                                         stacked_position=fmt.StackedPosition.Top)
                    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, reluctance_df["user_attrs_l_bot_air_gap"][index],
                                         stacked_position=fmt.StackedPosition.Bot)
                    geo.set_air_gaps(air_gaps)

                    # 4. set insulations
                    insulation = fmt.Insulation(flag_insulation=False)
                    insulation.add_core_insulations(config.insulations.iso_window_bot_core_top, config.insulations.iso_window_bot_core_bot,
                                                    config.insulations.iso_window_bot_core_left, config.insulations.iso_window_bot_core_right)
                    insulation.add_winding_insulations([[config.insulations.iso_primary_to_primary, config.insulations.iso_secondary_to_secondary],
                                                        [config.insulations.iso_secondary_to_secondary, config.insulations.iso_primary_to_primary]])
                    geo.set_insulation(insulation)

                    winding_window_top, winding_window_bot = fmt.create_stacked_winding_windows(core, insulation)

                    vww_top = winding_window_top.split_window(fmt.WindingWindowSplit.NoSplit)
                    vww_bot = winding_window_bot.split_window(fmt.WindingWindowSplit.NoSplit)

                    # 6. set conductor parameters
                    winding1 = fmt.Conductor(0, fmt.Conductivity.Copper)
                    primary_litz_wire = fmt.litz_database()[reluctance_df['params_primary_litz_wire'][index]]
                    winding1.set_litz_round_conductor(primary_litz_wire['conductor_radii'], primary_litz_wire['strands_numbers'],
                                                      primary_litz_wire['strand_radii'], None, fmt.ConductorArrangement.Square)

                    winding2 = fmt.Conductor(1, fmt.Conductivity.Copper)
                    secondary_litz_wire = fmt.litz_database()[reluctance_df['params_primary_litz_wire'][index]]
                    winding2.set_litz_round_conductor(secondary_litz_wire['conductor_radii'], secondary_litz_wire['strands_numbers'],
                                                      secondary_litz_wire['strand_radii'], None, fmt.ConductorArrangement.Square)

                    primary_coil_turns = reluctance_df['params_n_p_top'][index]
                    print(f"{primary_coil_turns=}")
                    print(f"{reluctance_df['params_n_p_bot'][index]=}")
                    print(f"{reluctance_df['user_attrs_n_s_bot'][index]=}")

                    # 7. add conductor to vww and add winding window to MagneticComponent
                    vww_top.set_interleaved_winding(winding1, primary_coil_turns, winding2, 0, fmt.InterleavedWindingScheme.HorizontalAlternating)
                    vww_bot.set_interleaved_winding(winding1, reluctance_df['params_n_p_bot'][index], winding2,
                                                    int(reluctance_df['user_attrs_n_s_bot'][index]), fmt.InterleavedWindingScheme.HorizontalAlternating)

                    geo.set_winding_windows([winding_window_top, winding_window_bot])

                    geo.create_model(freq=target_and_fix_parameters.fundamental_frequency, pre_visualize_geometry=show_visual_outputs, save_png=False)
                    geo.stacked_core_study(number_primary_coil_turns=primary_coil_turns, time_current_vectors=time_current_vectors,
                                           plot_waveforms=show_visual_outputs, fft_filter_value_factor=0.05)

                    result_dict = geo.read_log()
                    print(f"{result_dict=}")
                    df_single_simulation = pd.DataFrame(result_dict)

                    df = pd.concat([df, df_single_simulation], axis=0)
                except Exception as e:
                    print(e)

            return df
