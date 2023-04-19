# python libraries
import os

# 3rd party libraries
import pandas as pd
import numpy as np
import datetime

#femmt libraries
import femmt as fmt



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
    core_inner_diameter_min_max_list= [pq3230["core_inner_diameter"], pq5050["core_inner_diameter"]],
    window_w_min_max_list= [pq3230["window_w"], pq5050["window_w"]],
    window_h_top_min_max_list= [5 / 6 * pq3230["window_h"], 5 / 6 * pq5050["window_h"]],
    window_h_bot_min_max_list= [1 / 6 * pq3230["window_h"], 1 / 6 * pq5050["window_h"]],
    factor_max_flux_density = 1,
    n_p_top_min_max_list = [1, 100],
    n_p_bot_min_max_list = [0,30],
    n_s_top_min_max_list = [0,70],
    n_s_bot_min_max_list = [0,80],
    primary_litz_wire_list= ["1.4x200x0.071"],
    secondary_litz_wire_list= ["1.4x200x0.071"],
    temperature=100,
    working_directory=os.path.join(os.path.dirname(__file__), "example_results", "optuna_integrated_transformer_optimization")
)


task = 'start_study'
#task = 'fem_simulation_best_trials'
#task = 'fem_simulation_best_trials_offset'
#task = 'plot_study_results'

study_name = "workflow_2023-04-15"

if __name__ == '__main__':
    time_start = datetime.datetime.now()

    object1 = fmt.ItoOptuna()

    if task == 'start_study':
        object1.ReluctanceModel.start_study(study_name, dab_transformer_config, 20000, storage='sqlite')

    elif task == 'fem_simulation_best_trials':
        reluctance_result_list = object1.ReluctanceModel.load_study_best_trials(study_name)

        pareto_reluctance_result_list = object1.ReluctanceModel.load_study_pareto_area(study_name, 5)


        fmt.ItoOptuna.FemSimulation.simulate(dab_transformer_config,
                                             reluctance_result_list,
                                             visualize = False)

    elif task == 'fem_simulation_best_trials_offset':
        reluctance_result_list = object1.ReluctanceModel.load_study_pareto_area(study_name, 10)

        pareto_reluctance_result_list = object1.ReluctanceModel.load_study_pareto_area(study_name, 1)


        fmt.ItoOptuna.FemSimulation.simulate(dab_transformer_config,
                                             reluctance_result_list,
                                             visualize = True)

    elif task == 'plot_study_results':
        fmt.ItoOptuna.ReluctanceModel.show_study_results(study_name)


    time_stop = datetime.datetime.now()

    time_difference = time_stop - time_start
    print(f"{time_difference = }")