"""Example how to perform an optimization using optuna and FEMMT."""
# python libraries
import os

# 3rd party libraries
import numpy as np
import datetime

# femmt libraries
import femmt as fmt

core_database = fmt.core_database()
pq3230 = core_database["PQ 32/30"]
pq4040 = core_database["PQ 40/40"]
pq5050 = core_database["PQ 50/50"]
pq6560 = core_database["PQ 65/60"]


i_1 = [[0.0, 3.265248131976911e-07, 2.5e-06, 2.8265248131976912e-06, 5e-06],
       [-0.9996115022426437, 4.975792579275104, 0.9996115022426446, -4.975792579275103, -0.9996115022426437]]
i_2 = [[0.0, 3.265248131976911e-07, 2.5e-06, 2.8265248131976912e-06, 5e-06],
       [-0.9196195846583147, -19.598444313231134, 0.9196195846583122, 19.59844431323113, -0.9196195846583147]]

ito_insulations = fmt.ItoInsulation(
    # insulation for top core window
    iso_window_top_core_top=1e-3,
    iso_window_top_core_bot=1e-3,
    iso_window_top_core_left=1e-3,
    iso_window_top_core_right=1e-3,
    # insulation for bottom core window
    iso_window_bot_core_top=1e-3,
    iso_window_bot_core_bot=1e-3,
    iso_window_bot_core_left=1e-3,
    iso_window_bot_core_right=1e-3,
    # winding-to-winding insulation
    iso_primary_to_primary=10e-6,
    iso_secondary_to_secondary=10e-6,
    iso_primary_to_secondary=10e-6,
)

material_data_sources = fmt.IntegratedTransformerMaterialDataSources(
    permeability_datasource=fmt.MaterialDataSource.Measurement,
    permeability_datatype=fmt.MeasurementDataType.ComplexPermeability,
    permeability_measurement_setup=fmt.MeasurementSetup.MagNet,
    permittivity_datasource=fmt.MaterialDataSource.ManufacturerDatasheet,
    permittivity_datatype=fmt.MeasurementDataType.ComplexPermittivity,
    permittivity_measurement_setup=fmt.MeasurementSetup.LEA_MTB_small_signal
)

dab_transformer_config = fmt.ItoSingleInputConfig(
    integrated_transformer_study_name="2025-02-26",
    integrated_transformer_optimization_directory=os.path.join(os.path.dirname(__file__), "example_results", "optuna_integrated_transformer_optimization"),

    # target parameters
    l_s12_target=85e-6,
    l_h_target=600e-6,
    n_target=2.9,

    # fix parameters
    time_current_1_vec=np.array(i_1),
    time_current_2_vec=np.array(i_2),
    insulations=ito_insulations,

    # optimization parameters
    material_list=["3C95"],
    core_name_list=["PQ 40/40"],
    core_inner_diameter_min_max_list=None,
    window_w_min_max_list=None,
    window_h_top_min_max_list=None,
    window_h_bot_min_max_list=[1 / 6 * pq3230["window_h"], 1 / 6 * pq5050["window_h"]],

    n_1_top_min_max_list=[1, 30],
    n_1_bot_min_max_list=[1, 30],
    n_2_top_min_max_list=[1, 30],
    n_2_bot_min_max_list=[1, 30],
    factor_max_flux_density=1,
    primary_litz_wire_list=["1.4x200x0.071"],
    secondary_litz_wire_list=["1.4x200x0.071"],
    temperature=100,

    material_data_sources=material_data_sources,

)


# task = 'start_study'
# task = 'filter_reluctance_model'
# task = 'fem_simulation_from_filtered_reluctance_model_results'
task = 'plot_study_results'

if __name__ == '__main__':
    time_start = datetime.datetime.now()

    if task == 'start_study':
        fmt.IntegratedTransformerOptimization.ReluctanceModel.start_proceed_study(dab_transformer_config, 10000, storage='sqlite')

    elif task == 'filter_reluctance_model':
        # load trials from reluctance model to a pandas dataframe (df)
        reluctance_result_df = fmt.IntegratedTransformerOptimization.ReluctanceModel.study_to_df(dab_transformer_config)

        # filter air gaps
        # filtered_air_gaps_dto_list = fmt.IntegratedTransformerOptimization.ReluctanceModel.filter_min_air_gap_length(reluctance_result_df)

        # filter for Pareto front
        pareto_reluctance_dto_list = fmt.IntegratedTransformerOptimization.ReluctanceModel.filter_loss_list_df(
            reluctance_result_df, factor_min_dc_losses=0.5)

        fmt.IntegratedTransformerOptimization.ReluctanceModel.df_plot_pareto_front(
            reluctance_result_df, pareto_reluctance_dto_list, label_list=["all", "filtered"], interactive=False)

    elif task == 'fem_simulation_from_filtered_reluctance_model_results':
        # load filtered reluctance models
        pareto_reluctance_dto_list = fmt.IntegratedTransformerOptimization.ReluctanceModel.load_filtered_results(dab_transformer_config.working_directory)
        print(f"{len(pareto_reluctance_dto_list)=}")

        # start FEM simulation
        fmt.IntegratedTransformerOptimization.FemSimulation.simulate(config_dto=dab_transformer_config,
                                                                     simulation_dto_list=pareto_reluctance_dto_list)

    elif task == 'plot_study_results':
        fmt.IntegratedTransformerOptimization.ReluctanceModel.show_study_results(dab_transformer_config)

    time_stop = datetime.datetime.now()

    time_difference = time_stop - time_start
    print(f"{time_difference=}")
