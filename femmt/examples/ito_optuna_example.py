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

dab_transformer_config = fmt.ItoSingleInputConfig(
    l_s_target=85e-6,
    l_h_target=600e-6,
    n_target=2.9,
    time_current_1_vec=np.array(i_1),
    time_current_2_vec=np.array(i_2),
    material_list=["N95"],
    core_inner_diameter_min_max_list=[pq3230["core_inner_diameter"], pq5050["core_inner_diameter"]],
    window_w_min_max_list=[pq3230["window_w"], pq5050["window_w"]],
    window_h_top_min_max_list=[5 / 6 * pq3230["window_h"], 5 / 6 * pq5050["window_h"]],
    window_h_bot_min_max_list=[1 / 6 * pq3230["window_h"], 1 / 6 * pq5050["window_h"]],
    factor_max_flux_density=1,
    primary_litz_wire_list=["1.4x200x0.071"],
    secondary_litz_wire_list=["1.4x200x0.071"],
    temperature=100,
    working_directory=os.path.join(os.path.dirname(__file__), "example_results", "optuna_integrated_transformer_optimization")
)


task = 'start_study'
# task = 'filter_reluctance_model'
# task = 'fem_simulation_from_filtered_reluctance_model_results'
# task = 'plot_study_results'

study_name = "workflow_2023-04-15"

if __name__ == '__main__':
    time_start = datetime.datetime.now()

    if task == 'start_study':
        fmt.IntegratedTransformerOptimization.ReluctanceModel.NSGAII.start_study(study_name, dab_transformer_config, 1000, storage='sqlite')

    elif task == 'filter_reluctance_model':
        # load trials from reluctance model
        reluctance_result_list = fmt.IntegratedTransformerOptimization.ReluctanceModel.NSGAII.load_study_to_dto(study_name, dab_transformer_config)
        print(f"{len(reluctance_result_list) = }")

        # filter air gaps
        filtered_air_gaps_dto_list = fmt.IntegratedTransformerOptimization.ReluctanceModel.filter_min_air_gap_length(reluctance_result_list)
        print(f"{len(filtered_air_gaps_dto_list) = }")

        # filter for Pareto front
        pareto_reluctance_dto_list = fmt.IntegratedTransformerOptimization.ReluctanceModel.filter_loss_list(
            filtered_air_gaps_dto_list, factor_min_dc_losses=0.5)
        print(f"{len(pareto_reluctance_dto_list) = }")

        fmt.IntegratedTransformerOptimization.plot(reluctance_result_list)
        fmt.IntegratedTransformerOptimization.plot(pareto_reluctance_dto_list)

        # save results
        fmt.IntegratedTransformerOptimization.ReluctanceModel.save_dto_list(pareto_reluctance_dto_list, os.path.join(dab_transformer_config.working_directory,
                                                                                                                     '01_reluctance_model_results_filtered'))

    elif task == 'fem_simulation_from_filtered_reluctance_model_results':
        # load filtered reluctance models
        pareto_reluctance_dto_list = fmt.IntegratedTransformerOptimization.ReluctanceModel.load_filtered_results(dab_transformer_config.working_directory)
        print(f"{len(pareto_reluctance_dto_list) = }")

        # start FEM simulation
        fmt.IntegratedTransformerOptimization.FemSimulation.simulate(config_dto=dab_transformer_config,
                                                                     simulation_dto_list=pareto_reluctance_dto_list)

    elif task == 'plot_study_results':
        fmt.IntegratedTransformerOptimization.ReluctanceModel.NSGAII.show_study_results(study_name, dab_transformer_config)

    time_stop = datetime.datetime.now()

    time_difference = time_stop - time_start
    print(f"{time_difference = }")
