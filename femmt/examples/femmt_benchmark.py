"""Speed up of FEMMT by using parallel processing. File to generate benchmarks of different speed-up techniques.

It contains multiple benchmarking functions in order to analyse the runtime and the accuracy of the simulation results.
"""

# Python standard libraries
from typing import List
from dataclasses import dataclass
import os
import json
import numpy as np
import matplotlib.pyplot as plt

# Local libraries
import femmt as fmt

# ------ Classes ------

@dataclass
class MeshAccuracies:
    """Mesh accuracies for core, window, air gaps and conductors."""

    mesh_accuracy_core: float
    mesh_accuracy_window: float
    mesh_accuracy_conductor: float
    mesh_accuracy_air_gaps: float

class SingleBenchmark:
    """Define a single benchmark."""
    
    # Benchmark times
    high_level_geo_gen_time: float              # Time to create the model in gmsh
    generate_hybrid_mesh_time: float           # Time to generate the hybrid mesh
    generate_electro_magnetic_mesh_time: float  # Time to generate the electro magnetic mesh
    prepare_simulation_time: float              # Time to prepare the simulation (create GetDP files and set simulation options)
    real_simulation_time: float                 # Time to run the simulation in GetDP
    logging_time: float                         # Time for FEMMT to write the results in a log

    create_model_time: float                    # high_level_geo_gen_time + generate_hybrid_mesh_time
    simulation_time: float                      # generate_electro_magnetic_mesh_time + prepare_simulation_time + real_simulation_time + logging_time

    execution_time_: float                      # all measurements summed up

    flux_over_current: List[float]
    total_losses: float
    total_winding_losses: float

    def __init__(self, working_directory: str, model: fmt.MagneticComponent):
        self.working_directory = working_directory
        self.model = model

    def benchmark_simulation(self, inductor_frequency: int):
        """
        Start the simulation for a single benchmark.

        :param inductor_frequency: frequency in Hz
        :type inductor_frequency: int
        """
        # Simulate
        self.high_level_geo_gen_time, self.generate_hybrid_mesh_time = self.model.create_model(freq=inductor_frequency, pre_visualize_geometry=False,
                                                                                               save_png=False, benchmark=True)

        self.generate_electro_magnetic_mesh_time, self.prepare_simulation_time, self.real_simulation_time, self.logging_time = self.model.single_simulation(
            freq=inductor_frequency, current=[4.5], plot_interpolation=False, show_fem_simulation_results=False, benchmark=True)
       
        self.create_model_time = self.high_level_geo_gen_time + self.generate_hybrid_mesh_time
        self.simulation_time = self.generate_electro_magnetic_mesh_time + self.prepare_simulation_time + self.real_simulation_time + self.logging_time
        self.execution_time = self.create_model_time + self.simulation_time

        # Get simulation results
        results_log = os.path.join(self.working_directory, "results", "log_electro_magnetic.json")

        with open(results_log, "r") as fd:
            content = json.load(fd)

            self.flux_over_current = content["single_sweeps"][0]["winding1"]["flux_over_current"]
            self.total_losses = content["total_losses"]["total_losses"]
            self.winding_losses = content["total_losses"]["all_windings"]

    def to_dict(self):
        """Export the single benchmark results to a dict."""
        return {
            "create_model_time": self.create_model_time,
            "simulation_time": self.simulation_time,
            "high_level_geo_gen_time": self.simulation_time,
            "generate_hybrid_mesh_time": self.generate_hybrid_mesh_time,
            "generate_electro_magnetic_mesh_time": self.generate_electro_magnetic_mesh_time,
            "prepare_simulation_time": self.prepare_simulation_time,
            "real_simulation_time": self.real_simulation_time,
            "logging_time": self.logging_time,
            "execution_time": self.execution_time,
            "flux_over_current": self.flux_over_current,
            "total_losses": self.total_losses,
            "total_winding_losses": self.total_winding_losses
        }

class Benchmark:
    """Benchmarks different mesh accuracies for a general comparison."""

    benchmarks: SingleBenchmark

    def __init__(self, benchmarks):
        self.benchmarks = benchmarks

    def to_json(self, file_path):
        """Export the benchmark results to a .json file."""
        output = {
            "averages": {
                "create_model_time": np.mean([x.create_model_time for x in self.benchmarks]),
                "simulation_time": np.mean([x.simulation_time for x in self.benchmarks]),
                "high_level_geo_gen_time": np.mean([x.high_level_geo_gen_time for x in self.benchmarks]),
                "generate_hybrid_mesh_time": np.mean([x.generate_hybrid_mesh_time for x in self.benchmarks]),
                "generate_electro_magnetic_mesh_time": np.mean([x.generate_electro_magnetic_mesh_time for x in self.benchmarks]),
                "prepare_simulation_time": np.mean([x.prepare_simulation_time for x in self.benchmarks]),
                "real_simulation_time": np.mean([x.real_simulation_time for x in self.benchmarks]),
                "logging_time": np.mean([x.logging_time for x in self.benchmarks]),
                "execution_time": np.mean([x.execution_time for x in self.benchmarks]),
            },
            "variances": {
                "create_model_time": np.var([x.create_model_time for x in self.benchmarks]),
                "simulation_time": np.var([x.simulation_time for x in self.benchmarks]),
                "high_level_geo_gen_time": np.var([x.high_level_geo_gen_time for x in self.benchmarks]),
                "generate_hybrid_mesh_time": np.var([x.generate_hybrid_mesh_time for x in self.benchmarks]),
                "generate_electro_magnetic_mesh_time": np.var([x.generate_electro_magnetic_mesh_time for x in self.benchmarks]),
                "prepare_simulation_time": np.var([x.prepare_simulation_time for x in self.benchmarks]),
                "real_simulation_time": np.var([x.real_simulation_time for x in self.benchmarks]),
                "logging_time": np.var([x.logging_time for x in self.benchmarks]),
                "execution_time": np.var([x.execution_time for x in self.benchmarks]),
            },
            "benchmarks": [x.to_dict() for x in self.benchmarks] 
        }

        with open(file_path, "w") as fd:
            json.dump(output, fd, indent=2)

# ------ Generic Functions ------


def create_model(working_directory, mesh_accuracies=None, aspect_ratio=10, wwr_enabled=True, number_of_conductors: int = 9):
    """Create the model for benchmark."""
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)
    inductor_frequency = 270000
    
    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=working_directory,
                                verbosity=fmt.Verbosity.Silent, wwr_enabled=wwr_enabled)

    core_db = fmt.core_database()["PQ 40/40"]
    core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                    window_w=core_db["window_w"],
                                                    window_h=core_db["window_h"],
                                                    core_h=core_db["core_h"])
    core = fmt.Core(core_type=fmt.CoreType.Single,
                    core_dimensions=core_dimensions,
                    material=fmt.Material.N49, temperature=45, frequency=inductor_frequency,
                    permeability_datasource=fmt.MaterialDataSource.Measurement,
                    permeability_datatype=fmt.MeasurementDataType.ComplexPermeability,
                    permeability_measurement_setup=fmt.MeasurementSetup.LEA_LK,
                    permittivity_datasource=fmt.MaterialDataSource.Measurement,
                    permittivity_datatype=fmt.MeasurementDataType.ComplexPermittivity,
                    permittivity_measurement_setup=fmt.MeasurementSetup.LEA_LK, mdb_verbosity=fmt.Verbosity.Silent)
    
    if mesh_accuracies:
        if isinstance(mesh_accuracies, MeshAccuracies):
            geo.update_mesh_accuracies(mesh_accuracies.mesh_accuracy_core, mesh_accuracies.mesh_accuracy_window,
                                       mesh_accuracies.mesh_accuracy_conductor, mesh_accuracies.mesh_accuracy_air_gaps)
        elif isinstance(mesh_accuracies, float):
            geo.update_mesh_accuracies(mesh_accuracies, mesh_accuracies, mesh_accuracies, mesh_accuracies)
    geo.set_core(core)
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
    geo.set_air_gaps(air_gaps)
    insulation = fmt.Insulation(aspect_ratio, aspect_ratio > 0)
    insulation.add_core_insulations(0.001, 0.001, 0.004, 0.001)
    insulation.add_winding_insulations([[0.0005]])
    geo.set_insulation(insulation)
    winding_window = fmt.WindingWindow(core, insulation)
    vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)
    winding = fmt.Conductor(0, fmt.Conductivity.Copper)
    winding.set_solid_round_conductor(conductor_radius=0.0013, conductor_arrangement=fmt.ConductorArrangement.Square)
    winding.parallel = False
    vww.set_winding(winding, number_of_conductors, None)
    geo.set_winding_windows([winding_window])

    return geo


def create_rectangular_conductor_model(working_directory, mesh_accuracies, thickness, center_factor, left_bound_delta=None, wwr_enabled=True, aspect_ratio=10):
    """Create a model with rectangular condutors."""
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)
    inductor_frequency = 270000
    
    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=working_directory,
                                verbosity=fmt.Verbosity.Silent, wwr_enabled=wwr_enabled)

    core_db = fmt.core_database()["PQ 40/40"]
    core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                    window_w=core_db["window_w"],
                                                    window_h=core_db["window_h"],
                                                    core_h=core_db["core_h"])
    # core = fmt.Core(core_type=fmt.CoreType.Single,
    #                core_dimensions=core_dimensions,
    #                material=fmt.Material.N49, temperature=45, frequency=inductor_frequency,
    #                permeability_datasource=fmt.MaterialDataSource.Measurement,
    #                permeability_datatype=fmt.MeasurementDataType.ComplexPermeability,
    #                permeability_measurement_setup=fmt.MeasurementSetup.LEA_LK,
    #                permittivity_datasource=fmt.MaterialDataSource.Measurement,
    #                permittivity_datatype=fmt.MeasurementDataType.ComplexPermittivity,
    #                permittivity_measurement_setup=fmt.MeasurementSetup.LEA_LK, mdb_verbosity=fmt.Verbosity.Silent)
    
    # Used for FEMM comparison
    core = fmt.Core(core_type=fmt.CoreType.Single,
                    core_dimensions=core_dimensions,
                    material=fmt.Material.N49, temperature=45, frequency=inductor_frequency,
                    # permeability_datasource="manufacturer_datasheet",
                    permeability_datasource=fmt.MaterialDataSource.Custom,
                    # permeability_datatype=fmt.MeasurementDataType.ComplexPermeability,
                    # permeability_measurement_setup="LEA_LK",
                    mu_r_abs=1500, phi_mu_deg=0,
                    permittivity_datasource=fmt.MaterialDataSource.Custom,
                    # permittivity_datatype=fmt.MeasurementDataType.ComplexPermittivity,
                    # permittivity_measurement_setup="LEA_LK",
                    sigma=2.38732212414336e-001,
                    mdb_verbosity=fmt.Verbosity.Silent,)

    if mesh_accuracies:
        if isinstance(mesh_accuracies, MeshAccuracies):
            geo.update_mesh_accuracies(mesh_accuracies.mesh_accuracy_core, mesh_accuracies.mesh_accuracy_window,
                                       mesh_accuracies.mesh_accuracy_conductor, mesh_accuracies.mesh_accuracy_air_gaps)
        elif isinstance(mesh_accuracies, float):
            geo.update_mesh_accuracies(mesh_accuracies, mesh_accuracies, mesh_accuracies, mesh_accuracies)

    geo.set_core(core)
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
    geo.set_air_gaps(air_gaps)
    insulation = fmt.Insulation(aspect_ratio, not left_bound_delta)
    if left_bound_delta:
        insulation.add_core_insulations(0.001, 0.001, 0, 0)
    else:
        insulation.add_core_insulations(0.001, 0.001, 0.004, 0.001)
    insulation.add_winding_insulations([[0.0005]])
    geo.set_insulation(insulation)
    winding_window = fmt.WindingWindow(core, insulation)
    vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)
    winding = fmt.Conductor(0, fmt.Conductivity.Copper)
    winding.set_rectangular_conductor(thickness=thickness)
    winding.parallel = False
    vww.set_winding(winding, 14, fmt.WindingScheme.FoilHorizontal, fmt.WrapParaType.FixedThickness)
    
    if left_bound_delta:
        conductor_length = 0.006
        vww.left_bound += left_bound_delta
        vww.right_bound = vww.left_bound + conductor_length
        winding_window.max_left_bound += left_bound_delta
        winding_window.max_right_bound = winding_window.max_left_bound + conductor_length

    geo.set_winding_windows([winding_window])
    geo.mesh_data.center_factor = center_factor

    return geo

def plot_mesh_over_precision(benchmarks, x_list, x_label, title):
    """Plot the mesh over precision."""
    total_losses = [bm.total_losses for bm in benchmarks]
    self_inductance = [np.sqrt(bm.flux_over_current[0]**2+bm.flux_over_current[1]**2) for bm in benchmarks]
    execution_times = [bm.execution_time for bm in benchmarks]

    print(f"{x_label}: {x_list}")
    print(f"Total Losses: {total_losses}")
    print(f"Self inductance: {self_inductance}")
    print(f"Execution times: {execution_times}")

    figure, axis = plt.subplots(3, sharex=True)

    figure.suptitle(title, fontsize=20)
    axis[0].plot(x_list, total_losses, "bo", markersize=7)
    axis[0].set_ylabel("Total losses", fontsize=15)
    axis[0].set_xticks(x_list, fontsize=15)
    axis[0].grid()

    axis[1].plot(x_list, self_inductance, "bo", markersize=7)
    axis[1].set_ylabel("|Self inductance|", fontsize=15)
    axis[1].set_xticks(x_list, ddontsize=15)
    axis[1].grid()
    
    axis[2].plot(x_list, execution_times, "bo", markersize=7)
    axis[2].set_ylabel("Execution time", fontsize=15)
    axis[2].set_xlabel(x_label, fontsize=15)
    axis[2].set_xticks(x_list, fontsize=15)
    axis[2].grid()
    plt.show()

# ------ Benchmarks ------

def benchmark_rectangular_conductor_offset(working_directory):
    """Benchmark the conductor offset."""
    left_bound_deltas = np.linspace(0.001, 0.005, num=5)
    default_mesh_accuracy = 0.5
    thickness = 0.0015
    mesh_accuracies = [1, 0.8, 0.6, 0.4, 0.2]

    self_inductance = []
    winding_losses = []
    execution_times = []

    for mesh_accuracy in mesh_accuracies:
        current_mesh_accuracy_benchmarks = []
        for left_bound_delta in left_bound_deltas:
            wd = os.path.join(working_directory, f"{mesh_accuracy}")
            current_mesh_accuracy = MeshAccuracies(default_mesh_accuracy, default_mesh_accuracy, mesh_accuracy, default_mesh_accuracy)
            geo = create_rectangular_conductor_model(wd, current_mesh_accuracy, thickness, 1, left_bound_delta)
            benchmark = SingleBenchmark(wd, geo)
            benchmark.benchmark_simulation(270000)
            current_mesh_accuracy_benchmarks.append(benchmark)

        self_inductance.append([np.sqrt(bm.flux_over_current[0]**2+bm.flux_over_current[1]**2) for bm in current_mesh_accuracy_benchmarks])
        winding_losses.append([bm.winding_losses for bm in current_mesh_accuracy_benchmarks])
        execution_times.append([bm.execution_time for bm in current_mesh_accuracy_benchmarks])

    figure, axis = plt.subplots(3, sharex=True)
    for index, mesh_accuracy in enumerate(mesh_accuracies):
        current_self_inductance = [x/self_inductance[-1][i] for i, x in enumerate(self_inductance[index])]
        current_winding_losses = [x/winding_losses[-1][i] for i, x in enumerate(winding_losses[index])]
        current_execution_times = [x/execution_times[-1][i] for i, x in enumerate(execution_times[index])]
        axis[0].plot(left_bound_deltas, current_self_inductance, "o")
        axis[1].plot(left_bound_deltas, current_winding_losses, "o")
        axis[2].plot(left_bound_deltas, current_execution_times, "o", label=f"Mesh accuracy: {mesh_accuracy}")

    axis[0].set_ylabel("|Self indutance|")
    axis[0].set_xticks(left_bound_deltas)
    axis[0].grid()
    
    axis[1].set_ylabel("Total winding losses")
    axis[1].set_xticks(left_bound_deltas)
    axis[1].grid()

    axis[2].set_ylabel("Execution time")
    axis[2].set_xlabel("Left bound deltas")
    axis[2].set_xticks(left_bound_deltas)
    axis[2].grid()

    axis[2].legend()
    plt.show()


def benchmark_rectangular_conductor(working_directory):
    """Benchmark mesh accuracies inside a rectangular condutor."""
    # mesh_accuracies = np.arange(0.1, 1, step=0.1)
    # mesh_accuracies_conductor = np.linspace(0.001, 0.00005, num=10)
    default_mesh_accuracy = 0.4
    thickness = 0.0015
    # num_vertical_mesh_points = np.arange(1, 15, step=1.5)
    # mesh_accuracies_conductor = [thickness/x for x in num_vertical_mesh_points]
    mesh_accuracies_conductor = [0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
    winding_losses = []
    execution_times = []

    skin_depth = 0

    mesh_sizes = [thickness / 4 * x for x in mesh_accuracies_conductor]

    center_factors = [1, 2, 3, 4]
    
    femm_simulation_done = False
    femm_simulation_path = None

    for center_factor in center_factors:
        current_center_factor_benchmarks = []
        for mesh_accuracy_conductor in mesh_accuracies_conductor:
            wd = os.path.join(working_directory, f"{mesh_accuracy_conductor}")
            mesh_accuracy = MeshAccuracies(default_mesh_accuracy, default_mesh_accuracy, mesh_accuracy_conductor, default_mesh_accuracy)
            geo = create_rectangular_conductor_model(wd, mesh_accuracy, thickness, center_factor)
            benchmark = SingleBenchmark(wd, geo)
            benchmark.benchmark_simulation(270000)
            current_center_factor_benchmarks.append(benchmark)
            skin_depth = geo.mesh_data.delta
            
            if not femm_simulation_done:
                geo.femm_reference(freq=270000, current=[4.5], sign=[1], non_visualize=1, mesh_size=0.0, mesh_size_conductor=mesh_sizes[-1])
                femm_simulation_done = True
                femm_simulation_path = os.path.join(wd, "femm", "result_log_femm.json")

        winding_losses.append([bm.winding_losses for bm in current_center_factor_benchmarks])
        execution_times.append([bm.execution_time for bm in current_center_factor_benchmarks])

    figure, axis = plt.subplots(2, sharex=True)
    print(skin_depth)
    x_values = [mesh_size/skin_depth for mesh_size in mesh_sizes]

    figure.suptitle(f"Different meshes sizes for rectangular conductor; FEMM at mesh_size/skin_depth={mesh_sizes[-1]/skin_depth}")
    for index, center_factor in enumerate(center_factors):
        current_winding_losses = winding_losses[index]
        current_execution_times = execution_times[index]
        print(current_winding_losses)
        axis[0].plot(x_values, current_winding_losses, "o")
        axis[1].plot(x_values, current_execution_times, "o", label=f"Center factor: {center_factor}")

    femm_log_content = None
    with open(femm_simulation_path, "r") as fd:
        femm_log_content = json.load(fd)

    # Add FEMM horizontal line
    axis[0].axhline(y=femm_log_content["Winding1 Losses"], linestyle="--", color="black", label="FEMM Comparison")

    axis[0].set_ylabel("Total Winding Losses")
    axis[0].set_xticks(x_values)
    axis[0].grid()
    
    axis[1].set_ylabel("Execution time")
    axis[1].set_xlabel("Mesh size conductor/Skin depth")
    axis[1].set_xticks(x_values)
    axis[1].legend()
    axis[1].grid()
    plt.show()


def benchmark_different_mesh_accuracies(working_directory):
    """Simulate for different mesh accuracies and benchmark the results."""
    # comparison_data = {
    #    # Data from simulation with 0.05 mesh accuracy (everywhere) and no wwr and no insulations
    #    "total_losses": 7.681125170876381,
    #    "self_inductance": 3.536549787506707e-05,
    #    "execution_time": 239.89776849746704
    # }

    mesh_accuracy_values = np.linspace(0.3, 1, num=8)
    # mesh_accuracy_values = [1]
    default_accuracy = 0.3

    insulation_deltas = 0.00003
    wwr_enabled = False

    core_benchmarks = []
    window_benchmarks = []
    winding_benchmarks = []
    air_gaps_benchmarks = []
    comparison_benchmarks = []

    frequency = 270000
    
    # Core
    for accuracy_value in mesh_accuracy_values:
        wd = os.path.join(working_directory, "core")
        geo = create_model(wd, MeshAccuracies(accuracy_value, default_accuracy, default_accuracy, default_accuracy), 5, wwr_enabled)
        benchmark = SingleBenchmark(wd, geo)
        benchmark.benchmark_simulation(frequency)
        core_benchmarks.append(benchmark)
    
    # Window
    for accuracy_value in mesh_accuracy_values:
        wd = os.path.join(working_directory, "window")
        geo = create_model(wd, MeshAccuracies(default_accuracy, accuracy_value, default_accuracy, default_accuracy), 5, wwr_enabled)
        benchmark = SingleBenchmark(wd, geo)
        benchmark.benchmark_simulation(frequency)
        window_benchmarks.append(benchmark)
    
    # Winding
    for accuracy_value in mesh_accuracy_values:
        wd = os.path.join(working_directory, "winding")
        geo = create_model(wd, MeshAccuracies(default_accuracy, default_accuracy, accuracy_value, default_accuracy), 5, wwr_enabled)
        benchmark = SingleBenchmark(wd, geo)
        benchmark.benchmark_simulation(frequency)
        winding_benchmarks.append(benchmark)

    # Air gaps
    for accuracy_value in mesh_accuracy_values:
        wd = os.path.join(working_directory, "air_gaps")
        geo = create_model(wd, MeshAccuracies(default_accuracy, default_accuracy, default_accuracy, accuracy_value), 5, wwr_enabled)
        benchmark = SingleBenchmark(wd, geo)
        benchmark.benchmark_simulation(frequency)
        air_gaps_benchmarks.append(benchmark)
    
    # Comparison
    for accuracy_value in mesh_accuracy_values:
        wd = os.path.join(working_directory, "comparison")
        geo = create_model(wd, MeshAccuracies(accuracy_value, accuracy_value, accuracy_value, accuracy_value), 5, wwr_enabled)
        benchmark = SingleBenchmark(wd, geo)
        benchmark.benchmark_simulation(frequency)
        comparison_benchmarks.append(benchmark)
    
    font_size_default = 15
    markersize = 7
    font_size_title = 20

    figure, axis = plt.subplots(3, sharex=True)

    figure.suptitle(f"Benchmark: Different mesh accuracies based on position (specific region is coarser); default mesh size = {default_accuracy}",
                    fontsize=font_size_title)
    axis[0].plot(mesh_accuracy_values, [bm.total_losses for bm in core_benchmarks], "ro", label="Fine winding window mesh", markersize=markersize)
    axis[0].plot(mesh_accuracy_values, [bm.total_losses for bm in window_benchmarks], "yo", label="Coarse winding window mesh", markersize=markersize)
    axis[0].plot(mesh_accuracy_values, [bm.total_losses for bm in winding_benchmarks], "go", label="Winding", markersize=markersize)
    axis[0].plot(mesh_accuracy_values, [bm.total_losses for bm in air_gaps_benchmarks], "bo", label="Air Gaps", markersize=markersize)
    axis[0].plot(mesh_accuracy_values, [bm.total_losses for bm in comparison_benchmarks], "co", label="Comparison", markersize=markersize)
    # axis[0].axhline(y=comparison_data["total_losses"], linestyle="--", label="Ideal value")
    axis[0].set_ylabel("Total losses", fontsize=font_size_default)
    axis[0].set_xticks(mesh_accuracy_values, fontsize=font_size_default)
    axis[0].grid()

    axis[1].plot(mesh_accuracy_values, [np.sqrt(bm.flux_over_current[0]**2+bm.flux_over_current[1]**2) for bm in core_benchmarks], "ro",
                 label="Fine winding window mesh", markersize=markersize)
    axis[1].plot(mesh_accuracy_values, [np.sqrt(bm.flux_over_current[0]**2+bm.flux_over_current[1]**2) for bm in window_benchmarks], "yo",
                 label="Coarse winding window mesh", markersize=markersize)
    axis[1].plot(mesh_accuracy_values, [np.sqrt(bm.flux_over_current[0]**2+bm.flux_over_current[1]**2) for bm in winding_benchmarks], "go",
                 label="Winding", markersize=markersize)
    axis[1].plot(mesh_accuracy_values, [np.sqrt(bm.flux_over_current[0]**2+bm.flux_over_current[1]**2) for bm in air_gaps_benchmarks], "bo",
                 label="Air Gaps", markersize=markersize)
    axis[1].plot(mesh_accuracy_values, [np.sqrt(bm.flux_over_current[0]**2+bm.flux_over_current[1]**2) for bm in comparison_benchmarks], "co",
                 label="Comparison", markersize=markersize)
    # axis[1].axhline(y=comparison_data["self_inductance"], linestyle="--", label="Ideal value")
    axis[1].set_ylabel("|Self inductance|", fontsize=font_size_default)
    axis[1].set_xticks(mesh_accuracy_values, fontsize=font_size_default)
    axis[1].grid()
    
    axis[2].plot(mesh_accuracy_values, [bm.execution_time for bm in core_benchmarks], "ro", label="Fine winding window mesh", markersize=markersize)
    axis[2].plot(mesh_accuracy_values, [bm.execution_time for bm in window_benchmarks], "yo", label="Coarse winding window mesh", markersize=markersize)
    axis[2].plot(mesh_accuracy_values, [bm.execution_time for bm in winding_benchmarks], "go", label="Winding", markersize=markersize)
    axis[2].plot(mesh_accuracy_values, [bm.execution_time for bm in air_gaps_benchmarks], "bo", label="Air Gaps", markersize=markersize)
    axis[2].plot(mesh_accuracy_values, [bm.execution_time for bm in comparison_benchmarks], "co", label="Comparison", markersize=markersize)
    # axis[2].axhline(y=comparison_data["execution_time"], linestyle="--", label="Ideal value")
    axis[2].set_ylabel("Execution time", fontsize=font_size_default)
    axis[2].set_xlabel("Mesh accuracy for set region", fontsize=font_size_default)
    axis[2].set_xticks(mesh_accuracy_values, fontsize=font_size_default)
    axis[2].grid()

    plt.legend(fontsize=font_size_default)
    plt.show()

def benchmark_aspect_ratios(folder: str):
    """Create a benchmark for the different insulation_deltas. Opens a plot in the end."""
    aspect_ratios = np.arange(0, 10, step=1)
    frequency = 270000
    benchmarks = []

    for aspect_ratio in aspect_ratios:
        working_directory = os.path.join(folder, f"{aspect_ratio}")
        model = create_model(working_directory, 0.4, aspect_ratio)
        benchmark = SingleBenchmark(working_directory, model)
        benchmark.benchmark_simulation(frequency)
        benchmarks.append(benchmark)

    plot_mesh_over_precision(benchmarks, aspect_ratios, "aspect ratio", "Benchmark: Different aspect ratios")

def benchmark_mesh_accuracy(folder: str):
    """Create a benchmark for the different mesh_accuracies. Opens a plot at the end."""
    mesh_accuracies = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1]
    # mesh_accuracies = [1]
    frequency = 270000 
    font_size_default = 15
    markersize = 7
    font_size_title = 20

    # With insulation and wwr
    current_folder = os.path.join(folder, "insulation_and_wwr")

    if not os.path.exists(current_folder):
        os.mkdir(current_folder)

    benchmarks_insulation_and_wwr = []

    for mesh_accuracy in mesh_accuracies:
        working_directory = os.path.join(current_folder, f"{mesh_accuracy}")
        model = create_model(working_directory, mesh_accuracy, 5)
        benchmark = SingleBenchmark(working_directory, model)
        benchmark.benchmark_simulation(frequency)
        benchmarks_insulation_and_wwr.append(benchmark)

    total_losses_insulation_and_wwr = [bm.total_losses for bm in benchmarks_insulation_and_wwr]
    self_inductance_insulation_and_wwr = [np.sqrt(bm.flux_over_current[0]**2+bm.flux_over_current[1]**2) for bm in benchmarks_insulation_and_wwr]
    execution_time_insulation_and_wwr = [bm.execution_time for bm in benchmarks_insulation_and_wwr]

    print("Insulation and WWR <Total losses>:", total_losses_insulation_and_wwr)
    print("Insulation and WWR <Self inductance>:", self_inductance_insulation_and_wwr)
    print("Insulation and WWR <Execution time>:", execution_time_insulation_and_wwr)

    # Without insulation but with wwr
    current_folder = os.path.join(folder, "wwr")

    if not os.path.exists(current_folder):
        os.mkdir(current_folder)
        
    benchmarks_wwr = []

    for mesh_accuracy in mesh_accuracies:
        working_directory = os.path.join(current_folder, f"{mesh_accuracy}")
        model = create_model(working_directory, mesh_accuracy, 5, True)
        benchmark = SingleBenchmark(working_directory, model)
        benchmark.benchmark_simulation(frequency)
        benchmarks_wwr.append(benchmark)

    total_losses_wwr = [bm.total_losses for bm in benchmarks_wwr]
    self_inductance_wwr = [np.sqrt(bm.flux_over_current[0]**2+bm.flux_over_current[1]**2) for bm in benchmarks_wwr]
    execution_time_wwr = [bm.execution_time for bm in benchmarks_wwr]

    print("WWR <Total losses>:", total_losses_wwr)
    print("WWR <Self inductance>:", self_inductance_wwr)
    print("WWR <Execution time>:", execution_time_wwr)

    # With insulation but without wwr
    current_folder = os.path.join(folder, "insulation")

    if not os.path.exists(current_folder):
        os.mkdir(current_folder)
        
    benchmarks_insulation = []

    for mesh_accuracy in mesh_accuracies:
        working_directory = os.path.join(current_folder, f"{mesh_accuracy}")
        model = create_model(working_directory, mesh_accuracy, 5, False)
        benchmark = SingleBenchmark(working_directory, model)
        benchmark.benchmark_simulation(frequency)
        benchmarks_insulation.append(benchmark)

    total_losses_insulation = [bm.total_losses for bm in benchmarks_insulation]
    self_inductance_insulation = [np.sqrt(bm.flux_over_current[0]**2+bm.flux_over_current[1]**2) for bm in benchmarks_insulation]
    execution_time_insulation = [bm.execution_time for bm in benchmarks_insulation]

    print("Insulation <Total losses>:", total_losses_insulation)
    print("Insulation <Self inductance>:", self_inductance_insulation)
    print("Insulation <Execution time>:", execution_time_insulation)

    # Without both
    current_folder = os.path.join(folder, "without_both")

    if not os.path.exists(current_folder):
        os.mkdir(current_folder)
        
    benchmarks = []

    for mesh_accuracy in mesh_accuracies:
        working_directory = os.path.join(current_folder, f"{mesh_accuracy}")
        model = create_model(working_directory, mesh_accuracy, 5, False)
        benchmark = SingleBenchmark(working_directory, model)
        benchmark.benchmark_simulation(frequency)
        benchmarks.append(benchmark)

    total_losses = [bm.total_losses for bm in benchmarks]
    self_inductance = [np.sqrt(bm.flux_over_current[0]**2+bm.flux_over_current[1]**2) for bm in benchmarks]
    execution_time = [bm.execution_time for bm in benchmarks]

    print("Without both <Total losses>:", total_losses)
    print("Without both <Self inductance>:", self_inductance)
    print("Without both <Execution time>:", execution_time)

    figure, axis = plt.subplots(3, sharex=True)

    figure.suptitle("Benchmark: Different mesh precisions", fontsize=font_size_title)
    axis[0].plot(mesh_accuracies, total_losses_insulation_and_wwr, "bo", label="Insulation+WWR", markersize=markersize)
    axis[0].plot(mesh_accuracies, total_losses_insulation, "ro", label="Insulation", markersize=markersize)
    axis[0].plot(mesh_accuracies, total_losses_wwr, "yo", label="WWR", markersize=markersize)
    axis[0].plot(mesh_accuracies, total_losses, "go", label="No Insulation, No WWR", markersize=markersize)
    axis[0].set_ylabel("Total losses", fontsize=font_size_default)
    axis[0].set_xticks(mesh_accuracies, fontsize=font_size_default)
    axis[0].grid()

    axis[1].plot(mesh_accuracies, self_inductance_insulation_and_wwr, "bo", label="Insulation+WWR", markersize=markersize)
    axis[1].plot(mesh_accuracies, self_inductance_insulation, "ro", label="Insulation", markersize=markersize)
    axis[1].plot(mesh_accuracies, self_inductance_wwr, "yo", label="WWR", markersize=markersize)
    axis[1].plot(mesh_accuracies, self_inductance, "go", label="No Insulation, No WWR", markersize=markersize)
    axis[1].set_ylabel("|Self inductance|", fontsize=font_size_default)
    axis[1].set_xticks(mesh_accuracies, fontsize=font_size_default)
    axis[1].grid()
    
    axis[2].plot(mesh_accuracies, execution_time_insulation_and_wwr, "bo", label="Insulation+WWR", markersize=markersize)
    axis[2].plot(mesh_accuracies, execution_time_insulation, "ro", label="Insulation", markersize=markersize)
    axis[2].plot(mesh_accuracies, execution_time_wwr, "yo", label="WWR", markersize=markersize)
    axis[2].plot(mesh_accuracies, execution_time, "go", label="No Insulation, No WWR", markersize=markersize)
    axis[2].set_ylabel("Execution time", fontsize=font_size_default)
    axis[2].set_xlabel("Mesh accuracy", fontsize=font_size_default)
    axis[2].set_xticks(mesh_accuracies, fontsize=font_size_default)
    axis[2].set_yscale("log", base=2)
    axis[2].grid()
    plt.legend()
    plt.show()

def single_benchmark(benchmark_results_folder, frequency):
    """Run many simulations of the same model and writes the execution times as well as a multiple other measured times in a log file."""
    benchmark_results_file = os.path.join(benchmark_results_folder, "benchmark_test_without_insulation_bug.json")
    benchmarks = []

    for _ in range(4):
        working_directory = os.path.join(benchmark_results_folder, "without_bug")
        model = create_model(working_directory)
        benchmark = SingleBenchmark(working_directory, model)
        benchmark.benchmark_simulation(frequency)
        benchmarks.append(benchmark)
    benchmark = Benchmark(benchmarks)
    benchmark.to_json(benchmark_results_file)


def benchmark_winding_window_rasterization(folder):
    """Introduce different winding counts to benachmark the winding window rasterization."""
    numbers_of_conductors = np.arange(2, 10, step=2)
    # numbers_of_conductors = [1]
    benchmarks_no_wwr = []
    benchmarks_wwr = []
    frequency = 270000

    # No WWR
    for conductors in numbers_of_conductors:
        working_directory = os.path.join(benchmark_results_folder, f"{conductors}")
        model = create_model(working_directory, 0.5, 5, False)
        benchmark = SingleBenchmark(working_directory, model)
        benchmark.benchmark_simulation(frequency)
        benchmarks_no_wwr.append(benchmark)
        
    # WWR
    for conductors in numbers_of_conductors:
        working_directory = os.path.join(benchmark_results_folder, f"{conductors}")
        model = create_model(working_directory, 0.5, 5, True)
        benchmark = SingleBenchmark(working_directory, model)
        benchmark.benchmark_simulation(frequency)
        benchmarks_wwr.append(benchmark)

    font_size_default = 15
    markersize = 7
    font_size_title = 20

    figure, axis = plt.subplots(3, sharex=True)

    figure.suptitle("Comparing coarser winding window with default winding window", fontsize=font_size_title)
    axis[0].plot(numbers_of_conductors, [bm.total_losses for bm in benchmarks_no_wwr], "ro", label="No WWR", markersize=markersize)
    axis[0].plot(numbers_of_conductors, [bm.total_losses for bm in benchmarks_wwr], "yo", label="WWR", markersize=markersize)
    # axis[0].axhline(y=comparison_data["total_losses"], linestyle="--", label="Ideal value")
    axis[0].set_ylabel("Total losses", fontsize=font_size_default)
    axis[0].set_xticks(numbers_of_conductors, fontsize=font_size_default)
    axis[0].grid()

    axis[1].plot(numbers_of_conductors, [np.sqrt(bm.flux_over_current[0]**2+bm.flux_over_current[1]**2) for bm in benchmarks_no_wwr], "ro", label="No WWR",
                 markersize=markersize)
    axis[1].plot(numbers_of_conductors, [np.sqrt(bm.flux_over_current[0]**2+bm.flux_over_current[1]**2) for bm in benchmarks_wwr], "yo", label="WWR",
                 markersize=markersize)
    # axis[1].axhline(y=comparison_data["self_inductance"], linestyle="--", label="Ideal value")
    axis[1].set_ylabel("|Self inductance|", fontsize=font_size_default)
    axis[1].set_xticks(numbers_of_conductors, fontsize=font_size_default)
    axis[1].grid()
    
    axis[2].plot(numbers_of_conductors, [bm.execution_time for bm in benchmarks_no_wwr], "ro", label="No WWR", markersize=markersize)
    axis[2].plot(numbers_of_conductors, [bm.execution_time for bm in benchmarks_wwr], "yo", label="WWR", markersize=markersize)
    # axis[2].axhline(y=comparison_data["execution_time"], linestyle="--", label="Ideal value")
    axis[2].set_ylabel("Execution time", fontsize=font_size_default)
    axis[2].set_xlabel("Number of conductors", fontsize=font_size_default)
    axis[2].set_xticks(numbers_of_conductors, fontsize=font_size_default)
    axis[2].grid()

    plt.legend(fontsize=font_size_default)
    plt.show()

def general_comparison(folder: str):
    """
    Comparison with and without meshing techniques.

    :param folder: file path to folder
    :type folder: str
    """
    mesh_accuracies = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    # mesh_accuracies = [0.4, 0.5, 1]
    frequency = 270000

    font_size_default = 15
    markersize = 7
    font_size_title = 20

    all_on = []
    for mesh_accuracy in mesh_accuracies:
        working_directory = os.path.join(folder, f"{mesh_accuracy}_all_on")
        model = create_rectangular_conductor_model(working_directory, mesh_accuracy, 0.0015, 4, None, True, 1)
        benchmark = SingleBenchmark(working_directory, model)
        benchmark.benchmark_simulation(frequency)
        all_on.append(benchmark)
        
    all_off = []
    for mesh_accuracy in mesh_accuracies:
        working_directory = os.path.join(folder, f"{mesh_accuracy}_all_off")
        model = create_rectangular_conductor_model(working_directory, mesh_accuracy, 0.0015, 1, None, False, 15)
        benchmark = SingleBenchmark(working_directory, model)
        benchmark.benchmark_simulation(frequency)
        all_off.append(benchmark)

    total_losses_all_on = [bm.total_losses for bm in all_on]
    self_inductance_all_on = [np.sqrt(bm.flux_over_current[0]**2+bm.flux_over_current[1]**2) for bm in all_on]
    execution_time_all_on = [bm.execution_time for bm in all_on]

    total_losses_all_off = [bm.total_losses for bm in all_off]
    self_inductance_all_off = [np.sqrt(bm.flux_over_current[0]**2+bm.flux_over_current[1]**2) for bm in all_off]
    execution_time_all_off = [bm.execution_time for bm in all_off]

    figure, axis = plt.subplots(3, sharex=True)

    figure.suptitle("Comparison with and without meshing techniques", fontsize=font_size_title)
    axis[0].plot(mesh_accuracies, total_losses_all_on, "bo", label="With meshing techniques", markersize=markersize)
    axis[0].plot(mesh_accuracies, total_losses_all_off, "ro", label="Without meshing techniques", markersize=markersize)
    axis[0].set_ylabel("Total losses", fontsize=font_size_default)
    axis[0].set_xticks(mesh_accuracies, fontsize=font_size_default)
    axis[0].grid()

    axis[1].plot(mesh_accuracies, self_inductance_all_on, "bo", label="With meshing techniques", markersize=markersize)
    axis[1].plot(mesh_accuracies, self_inductance_all_off, "ro", label="Without meshing techniques", markersize=markersize)
    axis[1].set_ylabel("|Self inductance|", fontsize=font_size_default)
    axis[1].set_xticks(mesh_accuracies, fontsize=font_size_default)
    axis[1].grid()
    
    axis[2].plot(mesh_accuracies, execution_time_all_on, "bo", label="With meshing techniques", markersize=markersize)
    axis[2].plot(mesh_accuracies, execution_time_all_off, "ro", label="Without meshing techniques", markersize=markersize)
    axis[2].set_ylabel("Execution time", fontsize=font_size_default)
    axis[2].set_xlabel("Mesh accuracy", fontsize=font_size_default)
    axis[2].set_xticks(mesh_accuracies, fontsize=font_size_default)
    axis[2].set_yscale("log")
    axis[2].grid()
    plt.legend()
    plt.show()


if __name__ == "__main__":
    example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
    benchmark_results_folder = os.path.join(example_results_folder, "benchmarks")
    benchmark_single_results_folder = os.path.join(benchmark_results_folder, "single_benchmarks")
    benchmark_accuracy_results_folder = os.path.join(benchmark_results_folder, "benchmark_accuracy")
    benchmark_aspect_ratio_folder = os.path.join(benchmark_results_folder, "benchmark_aspect_ratio")
    benchmark_different_accuracy_folder = os.path.join(benchmark_results_folder, "benchmark_different_accuracy")
    benchmark_rect_conductor_folder = os.path.join(benchmark_results_folder, "benchmark_rect_conductor")
    benchmark_rect_conductor_offset_folder = os.path.join(benchmark_results_folder, "benchmark_rect_conductor_offset")
    benchmark_winding_window_rasterization_folder = os.path.join(benchmark_results_folder, "benchmark_winding_window_rasterization")
    general_comparison_folder = os.path.join(benchmark_results_folder, "general_comparison")

    if not os.path.exists(benchmark_results_folder):
        os.mkdir(benchmark_results_folder)
    if not os.path.exists(benchmark_accuracy_results_folder):
        os.mkdir(benchmark_accuracy_results_folder)
    if not os.path.exists(benchmark_aspect_ratio_folder):
        os.mkdir(benchmark_aspect_ratio_folder)
    if not os.path.exists(benchmark_different_accuracy_folder):
        os.mkdir(benchmark_different_accuracy_folder)
    if not os.path.exists(benchmark_rect_conductor_folder):
        os.mkdir(benchmark_rect_conductor_folder)
    if not os.path.exists(benchmark_rect_conductor_offset_folder):
        os.mkdir(benchmark_rect_conductor_offset_folder)
    if not os.path.exists(benchmark_winding_window_rasterization_folder):
        os.mkdir(benchmark_winding_window_rasterization_folder)
    if not os.path.exists(general_comparison_folder):
        os.mkdir(general_comparison_folder)

    # single_benchmark(benchmark_single_results_folder, 270000)
    # benchmark_mesh_accuracy(benchmark_accuracy_results_folder)
    # benchmark_aspect_ratios(benchmark_aspect_ratio_folder)
    # benchmark_different_mesh_accuracies(benchmark_different_accuracy_folder)
    # benchmark_rectangular_conductor(benchmark_rect_conductor_folder)
    # benchmark_rectangular_conductor_offset(benchmark_rect_conductor_offset_folder)
    # benchmark_winding_window_rasterization(benchmark_winding_window_rasterization_folder)
    general_comparison(general_comparison_folder)
