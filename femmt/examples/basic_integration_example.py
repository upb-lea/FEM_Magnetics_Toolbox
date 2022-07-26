from femmt.femmt_enumerations import AirGapLegPosition, AirGapMethod
import femmt as fmt
import numpy as np
import os

geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor)

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
#air_gap_position = np.linspace(0.001, 0.007, 5)

goal_inductance = 0.0002

mc1 = MagneticCircuit(core_width_list, core_window_w_list, core_window_h_list, no_of_turns, n_air_gaps, air_gap_length)
mc1.core_reluctance()
mc1.air_gap_reluctance()







