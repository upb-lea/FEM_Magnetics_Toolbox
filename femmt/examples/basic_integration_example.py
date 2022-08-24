from femmt.femmt_enumerations import AirGapLegPosition, AirGapMethod
import femmt as fmt
import numpy as np
import os

core_db1 = fmt.core_database()["PQ 20/16"]
core_db2 = fmt.core_database()["PQ 20/20"]
core_db3 = fmt.core_database()["PQ 26/20"]
core_db4 = fmt.core_database()["PQ 26/25"]

core_width_list = [core_db1["core_w"], core_db2["core_w"], core_db3["core_w"], core_db4["core_w"]]
core_window_w_list = [core_db1["window_w"], core_db2["window_w"], core_db3["window_w"], core_db4["window_w"]]
core_window_h_list = [core_db1["window_h"], core_db2["window_h"], core_db3["window_h"], core_db4["window_h"]]

no_of_turns = np.linspace(5, 10, 1)
n_air_gaps = np.linspace(1, 5, 1)
air_gap_length = np.linspace(0.0001, 0.0005, 0.0001)
# air_gap_position = np.linspace(0.001, 0.007, 5)

goal_inductance = 0.0002

mc1 = MagneticCircuit(core_width_list, core_window_w_list, core_window_h_list, no_of_turns, n_air_gaps, air_gap_length)
mc1.core_reluctance()
mc1.air_gap_reluctance()


def aut_action_run_simulation(sim_value):
    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor)
    core = fmt.Core(core_w=sim_value[0], window_h=sim_value[1], window_w=sim_value[2],
                    mu_rel=3100, phi_mu_deg=12,
                    sigma=0.6)
    geo.set_core(core)

    # 3. set air gap parameters
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Manually, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.001, sim_value[5])
    geo.set_air_gaps(air_gaps)

    # 4. set conductor parameters: use solid wires
    winding = fmt.Winding(sim_value[3], 0, fmt.Conductivity.Copper,
                          fmt.WindingType.Primary,
                          fmt.WindingScheme.Square)
    winding.set_litz_conductor(None, 600, 35.5e-6, 0.6)
    # winding.set_solid_conductor(0.0015)
    geo.set_windings([winding])

    # 5. set isolations
    isolation = fmt.Isolation()
    isolation.add_core_isolations(0.001, 0.001, 0.002, 0.001)
    isolation.add_winding_isolations(0.0001)
    geo.set_isolation(isolation)

    # 5. create the model
    geo.create_model(freq=100000, visualize_before=False, save_png=False)

    # 6.a. start single simulation
    #geo.single_simulation(freq=100000, current=[3], show_results=False)

    # 6.b. start sweep simulation
    fs = [0, 10000, 30000, 60000, 100000, 150000]
    amplitude_list = [[10], [2], [1], [0.5], [0.2], [0.1]]
    phase_list = [[0], [10], [20], [30], [40], [50]]
    geo.excitation_sweep(frequency_list=fs, current_list_list=amplitude_list, phi_deg_list_list=phase_list)

"""
FEM_data_matrix = [core_w, window_h, window_w, no_of_turns, n_air_gaps, air_gap_h, inductance]
FEM_data_matrix = [[0.0149, 0.0295, 0.01105, 8, 1, 0.0001, 0.0002],
                   [0.0149, 0.0295, 0.01105, 8, 1, 0.0002, 0.0002],
                   [0.0149, 0.0295, 0.01105, 8, 1, 0.0003, 0.0002],
                   [0.0149, 0.0295, 0.01105, 8, 1, 0.0004, 0.0002]]"""

num_r = len(FEM_data_matrix)
m = []
for i in range(num_r):
    m.append([float(x) for x in FEM_data_matrix[i]])
    for i in m:
        aut_action_run_simulation(i)


