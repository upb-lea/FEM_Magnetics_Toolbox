from LogParser import * 
import os
import matplotlib.pyplot as plt
import numpy as np

directory = r"C:\Uni\Bachelorarbeit\sweeps\center core mesh"

axis_dirs = [os.path.join(directory, f"air_gap_axis_{i}") for i in range(5)]
rect_dirs = [os.path.join(directory, f"air_gap_rect_{i}") for i in range(5)]

axis_logs = FEMMTLogParser.get_log_files_from_working_directories(axis_dirs)
rect_logs = FEMMTLogParser.get_log_files_from_working_directories(rect_dirs)

parser_axis = FEMMTLogParser(axis_logs)
parser_rect = FEMMTLogParser(rect_logs)

inductance_axis = []
inductance_rect = []

for file_data in parser_axis.data.values():
    inductance_axis.append(file_data.sweeps[0].windings[0].self_inductance)

for file_data in parser_rect.data.values():
    inductance_rect.append(file_data.sweeps[0].windings[0].self_inductance)

air_gap_heights = [0.0020, 0.0010, 0.0005, 0.000025, 0.0000125]

plt.plot(air_gap_heights, inductance_axis, "bo", label="axis_points")
plt.plot(air_gap_heights, inductance_rect, "ro", label="rect")

plt.xlabel("air gap heights")
plt.ylabel("self inductance")

plt.legend()
plt.show()

