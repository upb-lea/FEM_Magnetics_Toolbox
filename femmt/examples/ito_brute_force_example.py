import femmt as fmt
import numpy as np
import os

core_database = fmt.core_database()
pq3230 = core_database["PQ 32/30"]
pq4040 = core_database["PQ 40/40"]
pq5050 = core_database["PQ 50/50"]
pq6560 = core_database["PQ 65/60"]


i_1 = [[0.0, 3.265248131976911e-07, 2.5e-06, 2.8265248131976912e-06, 5e-06], [-0.9996115022426437, 4.975792579275104, 0.9996115022426446, -4.975792579275103, -0.9996115022426437]]
i_2 = [[0.0, 3.265248131976911e-07, 2.5e-06, 2.8265248131976912e-06, 5e-06], [-0.9196195846583147, -19.598444313231134, 0.9196195846583122, 19.59844431323113, -0.9196195846583147]]

dab_transformer_config = fmt.ItoSingleInputConfig(
    l_s_target = 85e-6,
    l_h_target= 600e-6,
    n_target= 2.9,
    time_current_1_vec = np.array(i_1),
    time_current_2_vec = np.array(i_2),
    material_list = ["N95"],
    core_inner_diameter_min_max_list= [pq3230["core_inner_diameter"], pq6560["core_inner_diameter"], 4],
    window_w_min_max_list= [pq3230["window_w"], pq6560["window_w"], 3],
    window_h_top_min_max_list= [5 / 6 * pq3230["window_h"], 5 / 6 * pq6560["window_h"], 2],
    window_h_bot_min_max_list= [1 / 6 * pq3230["window_h"], 1 / 6 * pq6560["window_h"], 2],
    factor_max_flux_density = 1,
    n_p_top_min_max_list = [15,30],
    n_p_bot_min_max_list = [0,8],
    n_s_top_min_max_list = [0,7],
    n_s_bot_min_max_list = [0,7],
    primary_litz_wire_list= ["1.4x200x0.071"],
    secondary_litz_wire_list= ["1.4x200x0.071"],
    temperature=100,
    working_directory= os.path.join(os.path.dirname(__file__), "example_results", "integrated_transformer_optimization")
)

#task = 'simulation_reluctance'
#task = 'load_reluctance_and_filter'
#task = 'load_reluctance_filter_and_simulate_fem'
#task = 'plot_fem_simulations_results'
#task = "single_fem_simulation_from_reluctance_result"
task = 'load_fem_simulation_results_and_perform_thermal_simulations'


if task == 'simulation_reluctance':
    # reluctance model calculation
    #valid_reluctance_model_designs = fmt.integrated_transformer_optimization(dab_transformer_config)

    valid_reluctance_model_designs = fmt.IntegratedTransformerOptimization.ReluctanceModel.brute_force_calculation(dab_transformer_config)

    fmt.IntegratedTransformerOptimization.ReluctanceModel.save_unfiltered_results(config_file=dab_transformer_config, result_file_list=valid_reluctance_model_designs)

elif task == 'load_reluctance_and_filter':
    valid_reluctance_model_designs = fmt.IntegratedTransformerOptimization.ReluctanceModel.load_unfiltered_results(dab_transformer_config.working_directory)
    fmt.IntegratedTransformerOptimization.plot_reluctance_model_result_list(valid_reluctance_model_designs)
    print(f"{len(valid_reluctance_model_designs) = }")

    filtered_air_gaps_dto_list = fmt.IntegratedTransformerOptimization.ReluctanceModel.filter_min_air_gap_length(valid_reluctance_model_designs)
    #opt.plot_reluctance_model_result_list(filtered_air_gaps_dto_list)
    print(f"{len(filtered_air_gaps_dto_list) = }")

    pareto_reluctance_dto_list = fmt.IntegratedTransformerOptimization.ReluctanceModel.filter_loss_list(filtered_air_gaps_dto_list, factor_min_dc_losses=0.5)
    print(f"{len(pareto_reluctance_dto_list) = }")
    fmt.IntegratedTransformerOptimization.plot_reluctance_model_result_list(pareto_reluctance_dto_list)

    fmt.IntegratedTransformerOptimization.ReluctanceModel.save_dto_list(pareto_reluctance_dto_list, os.path.join(dab_transformer_config.working_directory, '01_reluctance_model_results_filtered'))

elif task == 'load_reluctance_filter_and_simulate_fem':
    pareto_reluctance_dto_list = fmt.IntegratedTransformerOptimization.ReluctanceModel.load_filtered_results(dab_transformer_config.working_directory)
    print(f"{len(pareto_reluctance_dto_list) = }")

    fmt.IntegratedTransformerOptimization.FemSimulation.simulate(config_dto=dab_transformer_config,
                                                                simulation_dto_list=pareto_reluctance_dto_list)

elif task == 'plot_fem_simulations_results':
    fmt.IntegratedTransformerOptimization.FemSimulation.plot_unfiltered_results(dab_transformer_config.working_directory)

elif task == 'load_fem_simulation_results_and_perform_thermal_simulations':

    valid_reluctance_model_designs = fmt.IntegratedTransformerOptimization.ReluctanceModel.load_filtered_results(dab_transformer_config.working_directory)

    fmt.integrated_transformer_fem_thermal_simulations_from_result_dtos(config_dto=dab_transformer_config, simulation_dto_list = valid_reluctance_model_designs, visualize = False)

elif task == 'plot_thermal_fem_simulation_results':
    pass