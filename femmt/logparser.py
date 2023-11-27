import os
import json
import matplotlib.pyplot as plt
from typing import List, Dict
from enum import Enum
from dataclasses import dataclass


class SweepTypes(Enum):
    SingleSweep = "single_sweeps"


@dataclass
class WindingData:
    flux: complex
    turns: int
    self_inductance: complex
    voltage: complex
    current: complex
    active_power: float
    reactive_power: float
    apparent_power: float


@dataclass
class SweepData:
    frequency: float
    core_eddy_losses: float
    core_hyst_losses: float
    winding_losses: float
    windings: List[WindingData]


@dataclass
class FileData:
    file_path: str
    sweeps: List[SweepData]
    total_winding_losses: float
    total_core_eddy_losses: float
    total_core_hyst_losses: float
    total_core_losses: float
    core_2daxi_total_volume: float
    total_cost: float


class FEMMTLogParser:
    """Class to parse the electromagnetic_results_log file created by FEMMT.
    Creates a class structure from the file in order to easy access the data and create plots.
    """
    # Contains the complete data
    data: Dict[str, FileData]

    def __init__(self, file_paths_dict: Dict):
        """Creates the data dict out of the given file_paths.

        :param file_paths_dict: List of paths to every log file that should be added to the data.
        :type file_paths_dict: List[str]
        """
        self.data = {}

        for name, file_path in file_paths_dict.items():
            if not os.path.isfile(file_path):
                raise Exception(f"File {file_path} does not exist.")
        
            self.data[name] = self.parse_file(file_path, SweepTypes.SingleSweep)

    def plot_frequency_sweep_losses(self, data_names: List[str], loss_parameter: str, plot_label: str = "") -> None:
        """
        Example function for a possible sweep plot. Sweeps over the frequency of different simulations from
        one or multiple files.

        :param data_names: Name of the data (keys of data dict). If the list is empty every key will be taken.
        :param loss_parameter: Name of the variable from SweepData as str which will be set on the y-axis.
        :param plot_label: Title of the plot.
        """
        if len(data_names) == 0:
            data_names = self.data.keys()

        for data_name in data_names:
            freq_data = []
            for sweep in self.data[data_name].sweeps:
                freq_data.append([sweep.frequency, getattr(sweep, loss_parameter)])

            freq_data.sort(key=lambda x_vector: x_vector[0])

            x = [x[0] for x in freq_data]
            y = [y[1] for y in freq_data]

            plt.plot(x, y, "o", label=f"{data_name}")
        
        plt.legend()
        plt.title(plot_label)
        plt.xlabel("frequency")
        plt.ylabel(loss_parameter)
        plt.show()

    def plot_frequency_sweep_winding_params(self, data_names: str, winding_number: int, winding_parameter: str,
                                            plot_label: str = "") -> None:
        """
        Example function for a possible sweep plot. Sweeps over the frequency of different simulations from
        one or multiple files.

        :param data_names: Name of the data (keys of data dict). If the list is empty every key will be taken.
        :type data_names: str
        :param winding_number: Number of winding which shall be compared.
        :type winding_number: int
        :param plot_label: Title of the plot.
        :type plot_label: str
        :param winding_parameter:
        :type winding_parameter: str
        """
        if len(data_names) == 0:
            data_names = self.data.keys()

        for data_name in data_names:
            freq_data = []
            for sweep in self.data[data_name].sweeps:
                if len(sweep.windings) < winding_number:
                    raise Exception(f"Winding number {winding_number} is too high for the data")
                
                data_value = getattr(sweep.windings[0], winding_parameter)
                if type(data_value) == complex:
                    data_value = abs(data_value)

                freq_data.append([sweep.frequency, data_value])

            freq_data.sort(key=lambda x_vector: x_vector[0])

            x = [x[0] for x in freq_data]
            y = [y[1] for y in freq_data]
            
            plt.plot(x, y, "o", label=f"{data_name}")
        
        plt.legend()
        plt.title(plot_label)
        plt.xlabel("frequency")
        plt.ylabel(winding_parameter)
        plt.show()

    @staticmethod
    def get_log_files_from_working_directories(working_directories: List[str]) -> Dict:
        """
        Returns a dict containing the log files for each given working directory together with the name of the
        directory as key. For every working directory the local path to the log file
        is working_directory/results/log_electro_magnetic.json

        :param working_directories: Working directories.
        :return: Dictionary with the name of the directory as key and the log.json as value.
        """
        log_files = {}
        for working_directory in working_directories:
            name = os.path.basename(working_directory)
            log_files[name] = os.path.join(working_directory, "results", "log_electro_magnetic.json")

        return log_files

    @staticmethod
    def parse_complex(data):
        """Returns complex number if data is list with two elements.

        :param data: List of 2 values or single value
        :return: Either float or complex value
        """
        if type(data) is list:
            if len(data) == 2:
                return complex(data[0], data[1])
            else:
                raise Exception(f"In order to parse complex values the list needs a size of 2: {data}")
        return float(data)

    @staticmethod
    def parse_file(file_path: str, sweep_type: SweepTypes) -> FileData:
        """Internal function used to parse the JSON-File to a Class structure.

        :param file_path: Full path to file
        :type file_path: str
        :param sweep_type: Sweep type to parse from (most parent element in JSON-File)
        :type sweep_type: SweepTypes
        :raises Exception: _description_
        :return: Data stored in a class
        :rtype: FileData
        """
        if not sweep_type == SweepTypes.SingleSweep:
            raise Exception("Currently only single sweep is supported")

        with open(file_path, "r") as fd:
            full_data = json.loads(fd.read())

        sweep = full_data[sweep_type.value]

        sweeps_data = []
        for item in sweep:
            windings = []

            index = 1
            while f"winding{index}" in item:
                current_winding = item[f"winding{index}"]
                winding_data = {
                    "flux": FEMMTLogParser.parse_complex(current_winding["flux"]),
                    "turns": current_winding["number_turns"],
                    "self_inductance": FEMMTLogParser.parse_complex(current_winding["self_inductance"]),
                    "voltage": FEMMTLogParser.parse_complex(current_winding["V"]),
                    "current": FEMMTLogParser.parse_complex(current_winding["I"]),
                    "active_power": current_winding["P"],
                    "reactive_power": current_winding["Q"],
                    "apparent_power": current_winding["S"],
                }
                winding_data_class = WindingData(**winding_data)

                windings.append(winding_data_class)
                index += 1

            sweep_data = {
                "frequency": item["f"],
                "core_eddy_losses": item["core_eddy_losses"],
                "core_hyst_losses": item["core_hyst_losses"],
                "winding_losses": item["all_winding_losses"],
                "windings": windings,
            }

            sweep_data_class = SweepData(**sweep_data)

            sweeps_data.append(sweep_data_class)

        total = {
            "file_path": file_path,
            "sweeps": sweeps_data,
            "total_winding_losses": full_data["total_losses"]["all_windings"],
            "total_core_eddy_losses": full_data["total_losses"]["eddy_core"],
            "total_core_hyst_losses": full_data["total_losses"]["hyst_core_fundamental_freq"],
            "total_core_losses": full_data["total_losses"]["core"],
            "core_2daxi_total_volume": full_data["misc"]["core_2daxi_total_volume"],
            "total_cost": full_data["misc"]["total_cost_incl_margin"],
        }

        return FileData(**total)
