import femmt as fmt
import numpy as np
import os

core_database = fmt.core_database()
pq2016 = core_database["PQ 20/16"]
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
    core_inner_diameter_min_max_list= [pq2016["core_inner_diameter"], pq4040["core_inner_diameter"], 2],
    window_w_min_max_list= [pq2016["window_w"], pq4040["window_w"], 3],
    window_h_top_min_max_list= [5 / 6 * pq2016["window_h"], 5 / 6 * pq3230["window_h"], 2],
    window_h_bot_min_max_list= [1 / 6 * pq2016["window_h"], 1 / 6 * pq3230["window_h"], 2],
    factor_max_flux_density = 1,
    primary_litz_wire_list= ["1.4x200x0.071"],
    secondary_litz_wire_list= ["1.4x200x0.071"],
    temperature=100,
    working_directory= os.path.join(os.path.dirname(__file__), "example_results", "integrated_transformer_optimization")
)

task = 'simulation_reluctance'
#task = 'load_reluctance_and_filter'
#task = 'load_reluctance_filter_and_simulate_fem'
#task = 'plot_and_filter_fem_simulations_results'
#task = "single_fem_simulation_from_reluctance_result"
#task = 'load_fem_simulation_results_and_perform_thermal_simulations'
#task = 'load_and_filter_thermal_simulations'
#task = 'plot_thermal_fem_simulation_results'


if task == 'simulation_reluctance':

    # reluctance model brute force calculation
    valid_reluctance_model_designs = fmt.IntegratedTransformerOptimization.ReluctanceModel.BruteForce.brute_force_calculation(dab_transformer_config)

    # save all calculated and valid reluctance model calculations
    fmt.IntegratedTransformerOptimization.ReluctanceModel.save_unfiltered_results(config_file=dab_transformer_config, result_file_list=valid_reluctance_model_designs)

elif task == 'load_reluctance_and_filter':
    # load all calculated and valid reluctance model calculations
    valid_reluctance_model_designs = fmt.IntegratedTransformerOptimization.ReluctanceModel.load_unfiltered_results(dab_transformer_config.working_directory)
    print(f"{len(valid_reluctance_model_designs) = }")

    # filter air gaps
    filtered_air_gaps_dto_list = fmt.IntegratedTransformerOptimization.ReluctanceModel.filter_min_air_gap_length(valid_reluctance_model_designs)
    print(f"{len(filtered_air_gaps_dto_list) = }")

    # filter for Pareto front
    pareto_reluctance_dto_list = fmt.IntegratedTransformerOptimization.ReluctanceModel.filter_loss_list(filtered_air_gaps_dto_list, factor_min_dc_losses=0.5)
    print(f"{len(pareto_reluctance_dto_list) = }")

    # save results
    fmt.IntegratedTransformerOptimization.ReluctanceModel.save_dto_list(pareto_reluctance_dto_list, os.path.join(dab_transformer_config.working_directory, '01_reluctance_model_results_filtered'))

    # plot unfiltered and filtered Pareto planes
    fmt.IntegratedTransformerOptimization.plot(valid_reluctance_model_designs)
    fmt.IntegratedTransformerOptimization.plot(pareto_reluctance_dto_list)


elif task == 'load_reluctance_filter_and_simulate_fem':

    # load filtered reluctance models
    pareto_reluctance_dto_list = fmt.IntegratedTransformerOptimization.ReluctanceModel.load_filtered_results(dab_transformer_config.working_directory)
    print(f"{len(pareto_reluctance_dto_list) = }")

    # start FEM simulation
    fmt.IntegratedTransformerOptimization.FemSimulation.simulate(config_dto=dab_transformer_config,
                                                                simulation_dto_list=pareto_reluctance_dto_list)

elif task == 'plot_and_filter_fem_simulations_results':

    # load unfiltered results
    unfiltered_fem_results = fmt.IntegratedTransformerOptimization.FemSimulation.load_unfiltered_results(dab_transformer_config.working_directory)
    print(f"{len(unfiltered_fem_results) = }")

    # plot unfiltered results
    fmt.IntegratedTransformerOptimization.FemSimulation.plot(unfiltered_fem_results)

    # filter results
    filtered_fem_results = fmt.IntegratedTransformerOptimization.FemSimulation.filter_loss_list(unfiltered_fem_results, factor_min_dc_losses=0.1)
    print(f"{len(filtered_fem_results) = }")

    # save filtered results
    # fmt.IntegratedTransformerOptimization.FemSimulation.save_filtered_results(filtered_fem_results, dab_transformer_config.working_directory)

    # plot filtered results
    fmt.IntegratedTransformerOptimization.FemSimulation.plot(filtered_fem_results)

elif task == 'load_fem_simulation_results_and_perform_thermal_simulations':

    # load filtered FEM simulation results
    filtered_fem_simulation_results = fmt.IntegratedTransformerOptimization.FemSimulation.load_filtered_results(dab_transformer_config.working_directory)

    # start thermal FEM simulation
    fmt.IntegratedTransformerOptimization.ThermalSimulation.simulation(config_dto=dab_transformer_config, result_log_dict_list=filtered_fem_simulation_results)

elif task == 'plot_thermal_fem_simulation_results':
    # load thermal simulation results
    unfiltered_thermal_fem_results = fmt.IntegratedTransformerOptimization.ThermalSimulation.load_unfiltered_simulations(dab_transformer_config.working_directory)
    print(f"{len(unfiltered_thermal_fem_results) = }")

    # filter thermal FEM simulations
    filtered_thermal_fem_results = fmt.IntegratedTransformerOptimization.ThermalSimulation.filter_max_temperature(unfiltered_thermal_fem_results, 125, 125)
    print(f"{len(filtered_thermal_fem_results) = }")

    # load filtered FEM simulation results
    unfiltered_fem_simulation_results = fmt.IntegratedTransformerOptimization.FemSimulation.load_unfiltered_results(
        dab_transformer_config.working_directory)

    # find common cases of loss vs. volume and temperature vs. volume
    # du to the plot will contain the loss vs. volume area, so non-working cases of temperature vs. volume need
    # to be removed in loss vs. volume plane.
    valid_thermal_simulations = fmt.IntegratedTransformerOptimization.ThermalSimulation.find_common_cases(unfiltered_fem_simulation_results, filtered_thermal_fem_results)
    print(f"{len(valid_thermal_simulations) = }")

    # plot all thermal simulation results, and plot
    fmt.IntegratedTransformerOptimization.FemSimulation.plot(unfiltered_thermal_fem_results)
    fmt.IntegratedTransformerOptimization.FemSimulation.plot(valid_thermal_simulations)