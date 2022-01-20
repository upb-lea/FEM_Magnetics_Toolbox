import json
import numpy as np

def calculate_heat_flux_round_wire(power, radius):
    """
    :param power: losses in Watts
    :param radius: wire thickness in m
    """

    wire_thickness = radius

    return power/(np.pi * wire_thickness**2)

def calculate_heat_flux_core(power, dim_core, dim_window, dim_air_gaps):
    """
    :param power: losses in Watts
    :dimCore: dimensions of the core in the format [width, height]
    :dimWindow: dimensions of the winding window in the format [width, height]
    :dimAirGaps: dimensions of the air gaos in the format [[air_gap_1_width, air_gap_1_height], [air_gap_2_width, air_gap_2_height], ..]
    """

    return power/(dim_core[0]*dim_core[1] - dim_window[0]*dim_window[1]*sum([x[0]*x[1] for x in dim_air_gaps]))

def read_results_log(results_log_file_path):
    losses = {}
    with open(results_log_file_path, "r") as fd:
        content = json.loads(fd.read())
        losses = content["Losses"]

    return losses