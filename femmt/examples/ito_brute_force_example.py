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
    core_inner_diameter_min_max_list= [pq3230["core_inner_diameter"], pq6560["core_inner_diameter"], 3],
    window_w_min_max_list= [pq3230["window_w"], pq6560["window_w"], 3],
    window_h_top_min_max_list= [5 / 6 * pq3230["window_h"], 5 / 6 * pq6560["window_h"], 3],
    window_h_bot_min_max_list= [1 / 6 * pq3230["window_h"], 1 / 6 * pq6560["window_h"], 3],
    factor_max_flux_density = 1,
    n_p_top_min_max_list = [15,30],
    n_p_bot_min_max_list = [0,8],
    n_s_top_min_max_list = [0,7],
    n_s_bot_min_max_list = [0,7],
    primary_litz_wire_list= ["1.4x200x0.071", "1.1x60x0.1", "1.35x200x0.071", "1.8x512x0.05"],
    secondary_litz_wire_list= ["1.4x200x0.071", "1.35x200x0.071"],
    temperature=100,
    working_directory=''
)

#task = 'simulation_reluctance'
task = 'load_reluctance_and_filter'
#task = 'load_reluctance_filter_and_simulate_fem'
#task = 'load_fem_simulations'
#task = "single_fem_simulation_from_reluctance_result"
#task = 'fem_without_filter'
#task = 'visualize_case_geometry'
#task = 'load_reluctance_simulations'

example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
if not os.path.exists(example_results_folder):
    os.mkdir(example_results_folder)

# Working directory can be set arbitrarily
working_directory = os.path.join(example_results_folder, "integrated_transformer_optimization")
if not os.path.exists(working_directory):
    os.mkdir(working_directory)



opt = fmt.optimization.IntegratedTransformerOptimization(working_directory=working_directory)

if task == 'simulation_reluctance':
    # reluctance model calculation
    valid_reluctance_model_designs = opt.integrated_transformer_optimization(dab_transformer_config)

    opt.save_reluctance_model_result_list(config_file=dab_transformer_config, result_file_list=valid_reluctance_model_designs)

elif task == 'fem_without_filter':
    valid_reluctance_model_designs = opt.load_reluctance_model_result_list()
    opt.fem_simulation(config_dto=dab_transformer_config, simulation_dto_list=valid_reluctance_model_designs)

elif task == 'load_reluctance_and_filter':
    valid_reluctance_model_designs = opt.load_reluctance_model_result_list()
    opt.plot_reluctance_model_result_list(valid_reluctance_model_designs)
    print(f"{len(valid_reluctance_model_designs) = }")

    filtered_air_gaps_dto_list = opt.filter_min_air_gap_length(valid_reluctance_model_designs)
    #opt.plot_reluctance_model_result_list(filtered_air_gaps_dto_list)
    print(f"{len(filtered_air_gaps_dto_list) = }")

    pareto_reluctance_dto_list = opt.filter_reluctance_model_list(filtered_air_gaps_dto_list, factor_min_dc_losses=0.5)
    print(f"{len(pareto_reluctance_dto_list) = }")
    opt.plot_reluctance_model_result_list(pareto_reluctance_dto_list)

    opt.save_dto_list(pareto_reluctance_dto_list, os.path.join(opt.optimization_working_directory, 'filtered_results'))
    #opt.plot_filtered_pareto_result_list(filter_volume_list, filter_core_hyst_loss_list)

elif task == 'load_reluctance_filter_and_simulate_fem':
    valid_reluctance_model_designs = opt.load_reluctance_model_result_list()
    #opt.plot_reluctance_model_result_list(valid_reluctance_model_designs)
    print(f"{len(valid_reluctance_model_designs) = }")

    filtered_air_gaps_dto_list = opt.filter_min_air_gap_length(valid_reluctance_model_designs)
    #opt.plot_reluctance_model_result_list(filtered_air_gaps_dto_list)
    print(f"{len(filtered_air_gaps_dto_list) = }")

    pareto_reluctance_dto_list = opt.filter_reluctance_model_list(filtered_air_gaps_dto_list, factor_min_dc_losses=0.5)
    print(f"{len(pareto_reluctance_dto_list) = }")
    #opt.plot_reluctance_model_result_list(pareto_dto_list)

    opt.fem_simulation(config_dto=dab_transformer_config, simulation_dto_list=pareto_reluctance_dto_list)

    #opt.plot_filtered_pareto_result_list(filter_volume_list, filter_core_hyst_loss_list)

elif task == "single_fem_simulation_from_reluctance_result":
    pareto_reluctance_dto_list = opt.load_reluctance_model_list(os.path.join(opt.optimization_working_directory, 'reluctance_model_results'))
    print(len(pareto_reluctance_dto_list))
    for count, reluctance_dto in (enumerate(pareto_reluctance_dto_list)):
        if reluctance_dto.case == 0:
            index = count

    single_dto = pareto_reluctance_dto_list[index]
    print(single_dto)
    opt.fem_simulation(config_dto=dab_transformer_config, simulation_dto_list=[single_dto], visualize=True)

elif task == 'load_fem_simulations':
    opt.plot_fem_simulation_result_list(os.path.join(opt.optimization_working_directory, 'fem_simulation_results'))

elif task == 'load_reluctance_simulations':
    valid_reluctance_model_designs = opt.load_reluctance_model_result_list()
    opt.plot_reluctance_model_result_list(valid_reluctance_model_designs)

elif task == 'visualize_case_geometry':

    file = os.path.join(opt.integrated_transformer_fem_simulations_results_directory, "case_0.json")

    geo = fmt.MagneticComponent.decode_settings_from_log(file, opt.femmt_working_directory)

    geo.create_model(freq=100000, visualize_before=True, save_png=False)

    #geo.single_simulation(freq=100000, current=[4.5], show_results=True)