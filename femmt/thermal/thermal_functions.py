import json
import numpy as np
import os

def calculate_heat_flux_round_wire(power, radius):
    """
    :param power: losses in Watts
    :param radius: wire thickness in m
    """

    wire_thickness = radius

    # Power density for volumes W/m^3
    return power/(2*np.pi**2 * wire_thickness**2)

    # Power density for surfaces W/m^2
    #return power/(np.pi * wire_thickness**2)

def read_results_log(results_log_file_path):
    losses = {}
    if not os.path.exists(results_log_file_path):
        raise Exception(f"Losses file not found {results_log_file_path}.")
    with open(results_log_file_path, "r") as fd:
        content = json.loads(fd.read())
        losses = content["Losses"]

    return losses