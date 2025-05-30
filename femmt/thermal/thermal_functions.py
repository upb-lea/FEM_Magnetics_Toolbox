"""Functions for the thermal simulations."""
# Python standard libraries
import json
import os

# 3rd party libraries
import numpy as np


def calculate_heat_flux_round_wire(power, wire_radius, wire_distance):
    """
    Calculate teat flux for round wire.

    :param power: losses in Watts
    :type power: float
    :param wire_radius: wire thickness in m
    :type wire_radius: float
    :param wire_distance: length of the wire
    :type wire_distance: float
    """
    # Power density for volumes W/m^3
    volume = 2 * np.pi**2 * wire_radius**2 * wire_distance

    return power/volume


def read_results_log(results_log_file_path: str) -> dict:
    """
    Read loss results from logfile to get the losses to perform a thermal simulation.

    :param results_log_file_path: file path to result-log
    :type results_log_file_path: str
    """
    if not os.path.exists(results_log_file_path):
        raise Exception(f"Losses file not found {results_log_file_path}.")
    with open(results_log_file_path, "r") as fd:
        content = json.loads(fd.read())
        losses = content["total_losses"]

    return losses
