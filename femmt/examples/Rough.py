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
