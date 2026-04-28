"""Example how to perform an optimization using optuna and FEMMT."""
# python libraries
import os
import datetime

# 3rd party libraries
import numpy as np
import pandas as pd

# femmt libraries
import femmt as fmt
from materialdatabase.meta.data_enums import Material, DataSource

i_offset = 6
ac_factor = 1

# exact a single current waveform to optimize the inductor
time_vec: np.ndarray = np.array([0.0, 3.265248131976911e-07, 2.5e-06, 2.8265248131976912e-06, 5e-06])
ac_current_vec: np.ndarray = np.array([-0.9996115022426437, 4.975792579275104, 0.9996115022426446, -4.975792579275103, -0.9996115022426437])*ac_factor
current_vec: np.ndarray = i_offset + ac_current_vec

i_1: list[np.ndarray] = [time_vec, current_vec]

# Initialize the material data source
material_data_sources = fmt.MaterialDataSources(
    permeability_datasource=DataSource("MagNet"),
    # permeability_datasource=DataSource("LEA_MTB"),
    permittivity_datasource=DataSource("LEA_MTB"),
)

fmt_inductor_optimization_dto = fmt.InductorOptimizationDTO(
    inductor_study_name="Example_Inductor",
    core_name_list=["PQ 50/50", "PQ 50/40", "PQ 40/40", "PQ 40/30", "PQ 35/35", "PQ 32/30", "PQ 32/20", "PQ 26/25", "PQ 26/20", "PQ 20/20", "PQ 20/16"],
    material_name_list=[Material("N49")],
    core_inner_diameter_min_max_list=[],
    window_h_min_max_list=[],
    window_w_min_max_list=[],
    litz_wire_name_list=["1.5x105x0.1", "1.4x200x0.071", "1.1x60x0.1"],
    insulations=fmt.InductorInsulationDTO(primary_to_primary=0.2e-3,
                                          core_bot=1e-3,
                                          core_top=1e-3,
                                          core_right=1e-3,
                                          core_left=1e-3),
    target_inductance=5e-6,
    temperature=100,
    time_current_vec=i_1,
    inductor_optimization_directory=os.path.join(os.path.dirname(__file__), "example_results/ito_optuna"),
    material_data_sources=material_data_sources
)

# Set task according your need
task_start_proceed_study = True
task_filter_reluctance_model = True
task_reluctance_model_simulation_from_filtered_reluctance_model_results = True
task_fem_simulation_from_filtered_reluctance_model_results = True
task_plot_study_results = True

# Variable declaration
config_filtered_filepath = os.path.join(fmt_inductor_optimization_dto.inductor_optimization_directory,
                                        f"{fmt_inductor_optimization_dto.inductor_study_name}_filtered.csv")
config_filepath = os.path.join(fmt_inductor_optimization_dto.inductor_optimization_directory, f"{fmt_inductor_optimization_dto.inductor_study_name}.pkl")

if __name__ == '__main__':
    time_start = datetime.datetime.now()

    if task_start_proceed_study:
        fmt.optimization.InductorOptimization.ReluctanceModel.start_proceed_study(fmt_inductor_optimization_dto, 100, storage='sqlite')

    if task_filter_reluctance_model:
        # load trials from reluctance model to a pandas dataframe (df)
        reluctance_result_df = fmt.InductorOptimization.ReluctanceModel.study_to_df(fmt_inductor_optimization_dto)

        # filter for Pareto front
        pareto_reluctance_dto_df = fmt.InductorOptimization.ReluctanceModel.filter_loss_list_df(
            reluctance_result_df, factor_min_dc_losses=0.5)

        fmt.InductorOptimization.ReluctanceModel.df_plot_pareto_front(
            reluctance_result_df, pareto_reluctance_dto_df, label_list=["all", "filtered"], interactive=False)

        inductor_id_list_pareto = pareto_reluctance_dto_df["number"].to_numpy()

        pareto_reluctance_dto_df.to_csv(config_filtered_filepath)

    # load filtered reluctance data
    pareto_reluctance_dto_df = pd.read_csv(config_filtered_filepath)

    if task_reluctance_model_simulation_from_filtered_reluctance_model_results:

        # Filter out 5 simulations
        n = len(pareto_reluctance_dto_df)
        indices: list = []
        for i in range(5):
            indices.append(int(round(i * (n - 1) / 4)))

        # Single index
        indices = [int(round(n/2))]
        # Get filtered subset
        pareto_reluctance_dto_filtered_df = pareto_reluctance_dto_df.iloc[indices]

        df_geometry_re_simulation_number = pareto_reluctance_dto_filtered_df

        volume, combined_losses, area_to_heat_sink, winding_loss, core_loss = fmt.InductorOptimization.ReluctanceModel.full_simulation(
            df_geometry_re_simulation_number, current_waveform=i_1,
            inductor_config_filepath=config_filepath)

    if task_fem_simulation_from_filtered_reluctance_model_results:
        # load filtered reluctance models
        # df_geometry_re_simulation_number = pareto_reluctance_dto_df[pareto_reluctance_dto_df["number"]]

        # Filter out 5 simulations
        n = len(pareto_reluctance_dto_df)
        indices: list = []
        for i in range(5):
            indices.append(int(round(i * (n - 1) / 4)))

        # Single index
        indices = [int(round(n/2))]
        # Get filtered subset
        pareto_reluctance_dto_filtered_df = pareto_reluctance_dto_df.iloc[indices]

        df_geometry_re_simulation_number = pareto_reluctance_dto_filtered_df

        config_filepath = os.path.join(fmt_inductor_optimization_dto.inductor_optimization_directory,
                                       f"{fmt_inductor_optimization_dto.inductor_study_name}.pkl")

        volume, combined_losses, area_to_heat_sink, winding_loss, core_loss = fmt.InductorOptimization.FemSimulation.full_simulation(
            df_geometry_re_simulation_number, current_waveform=i_1,
            inductor_config_filepath=config_filepath, process_number=1, print_derivations=False)

    if task_plot_study_results:
        fmt.InductorOptimization.ReluctanceModel.show_study_results(fmt_inductor_optimization_dto)

    time_stop = datetime.datetime.now()

    time_difference = time_stop - time_start
    print(f"{time_difference=}")
