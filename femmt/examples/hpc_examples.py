# Python standard libraries
from typing import Dict
from itertools import product
import os
import time
import matplotlib.pyplot as plt
import statistics


# Local libraries
import femmt as fmt
import materialdatabase as mdb

"""This file contains examples for the use of the hpc function. Internally the multprocessing package is used.

IMPORTANT:
The __name__ == __main__ is necessary for the hpc function to work properly. So if a custom file is created from which the run_hpc function is
called please make sure to also include a __name__ == __main__ otherwise errors will be thrown.
For more information have a look here: https://docs.python.org/2/library/multiprocessing.html#windows
"""

# ---- Utility functions ----
def create_parallel_example_transformer() -> fmt.MagneticComponent:
    """Creates an example model which is used for the parallel execution example. This does implement a simple transformer.
    """ 
    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Transformer, working_directory=working_directory, verbosity=fmt.Verbosity.ToFile)
    core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=0.015, window_w=0.012, window_h=0.0295, core_h=0.015/2)
    core = fmt.Core(core_dimensions=core_dimensions, non_linear=False, sigma=1, re_mu_rel=3200, phi_mu_deg=10,
                    permeability_datasource=fmt.MaterialDataSource.Custom, permittivity_datasource=fmt.MaterialDataSource.Custom,
                    mdb_verbosity=fmt.Verbosity.Silent)
    geo.set_core(core)
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
    geo.set_air_gaps(air_gaps)
    insulation = fmt.Insulation()
    insulation.add_core_insulations(0.001, 0.001, 0.002, 0.001)
    insulation.add_winding_insulations([[0.0002, 0.001],
                                        [0.001, 0.0002]])
    geo.set_insulation(insulation)
    winding_window = fmt.WindingWindow(core, insulation)
    vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)
    winding1 = fmt.Conductor(0, fmt.Conductivity.Copper)
    winding1.set_solid_round_conductor(0.0011, None)
    winding2 = fmt.Conductor(1, fmt.Conductivity.Copper)
    winding2.set_solid_round_conductor(0.0011, None)
    vww.set_interleaved_winding(winding1, 21, winding2, 7, fmt.InterleavedWindingScheme.HorizontalAlternating)
    geo.set_winding_windows([winding_window])

    return geo

def create_parallel_example_inductor(inductor_frequency: int, air_gap_height: float = 0.0005, air_gap_position: int = 50) -> fmt.MagneticComponent:
    """Creates an example model which is used for the parallel execution example. This does implement a simple inductor with given inductor_frequency.

    :param inductor_frequency: Frequency for the inductor.
    :type inductor_frequency: int
    """ 
    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=None,  # Can be set to None since it will be overwritten anyways
                                clean_previous_results=False, verbosity=fmt.Verbosity.ToFile)
    core_db = fmt.core_database()["PQ 40/40"]
    core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                    window_w=core_db["window_w"],
                                                    window_h=core_db["window_h"],
                                                    core_h=core_db["core_h"])
    core = fmt.Core(core_type=fmt.CoreType.Single,
                    core_dimensions=core_dimensions,
                    material=mdb.Material.N49, temperature=45, frequency=inductor_frequency,
                    permeability_datasource=fmt.MaterialDataSource.Measurement,
                    permeability_datatype=fmt.MeasurementDataType.ComplexPermeability,
                    permeability_measurement_setup=mdb.MeasurementSetup.LEA_LK,
                    permittivity_datasource=fmt.MaterialDataSource.Measurement,
                    permittivity_datatype=fmt.MeasurementDataType.ComplexPermittivity,
                    permittivity_measurement_setup=mdb.MeasurementSetup.LEA_LK, mdb_verbosity=fmt.Verbosity.Silent)
    geo.set_core(core)
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, air_gap_height, air_gap_position)
    geo.set_air_gaps(air_gaps)
    insulation = fmt.Insulation()
    insulation.add_core_insulations(0.001, 0.001, 0.004, 0.001)
    insulation.add_winding_insulations([[0.0005]])
    geo.set_insulation(insulation)
    winding_window = fmt.WindingWindow(core, insulation)
    vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)
    winding = fmt.Conductor(0, fmt.Conductivity.Copper)
    winding.set_solid_round_conductor(conductor_radius=0.0013, conductor_arrangement=fmt.ConductorArrangement.Square)
    # winding.set_litz_round_conductor(conductor_radius=0.0013, number_strands=150, strand_radius=100e-6,fill_factor=None,
    # conductor_arrangement=fmt.ConductorArrangement.Square)
    winding.parallel = False
    vww.set_winding(winding, 9, None)
    geo.set_winding_windows([winding_window])

    return geo

def custom_hpc(parameters: Dict):
    """Very simple example for a custom hpc_function which can be given to the hpc.run() function.

    :param parameters: Dictionary containing the model and the given simulation_parameters.
    :type parameters: Dict
    """
    model = parameters["model"]
    simulation_parameters = parameters["simulation_parameters"]

    if "current" not in simulation_parameters:
        print("'current' argument is missing. Simulation will be skipped.")
        return
    if "phi_deg" not in simulation_parameters:
        print("'phi_deg' argument is missing. Simulation will be skipped.")
        return

    current = simulation_parameters["current"]
    phi_deg = simulation_parameters["phi_deg"]

    model.create_model(freq=250000, pre_visualize_geometry=False)
    model.single_simulation(freq=250000, current=current, phi_deg=phi_deg, show_fem_simulation_results=False)

def parallel_simulation_study(averaging_count):
    example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
    parallel_results_folder = os.path.join(example_results_folder, "parallel")
    study_results_folder = os.path.join(parallel_results_folder, "study")

    if not os.path.exists(parallel_results_folder):
        os.mkdir(parallel_results_folder)

    if not os.path.exists(study_results_folder):
        os.mkdir(study_results_folder)

    process_counts = [6, 7, 8]
    frequencies = [100000, 150000, 200000]
    air_gap_heights = [0.0002, 0.0005, 0.0007]
    air_gap_positions = [20, 40, 60, 80]

    models = []
    simulation_parameters = []

    runtimes = []

    for frequency, air_gap_height, air_gap_position in product(frequencies, air_gap_heights, air_gap_positions):
        models.append(create_parallel_example_inductor(frequency, air_gap_height, air_gap_position))
        simulation_parameters.append({
            "freq": frequency,
            "current": [1]})

    for process_count in process_counts:
        working_directory = os.path.join(study_results_folder, f"{process_count}")
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        simulation_times = []
        for count in range(averaging_count):
            start_time = time.time()
            fmt.run_hpc(process_count, models, simulation_parameters, working_directory)
            simulation_times.append(time.time() - start_time)

        runtimes.append(statistics.fmean(simulation_times))

    print("Process counts:", process_counts)
    print("Runtimes:", runtimes)
    
    plt.plot(process_counts, runtimes, "bo")
    plt.title(f"Parallel study ({len(models)} different models)")
    plt.xlabel("Number of processes")
    if averaging_count > 1:
        plt.ylabel(f"Runtime (mean of {averaging_count} simulations)")
    else:
        plt.ylabel(f"Runtime")
    plt.show()


if __name__ == "__main__":
    # ---- Choosing the execution ----
    execution_type = "default_example"
    # execution_type = "custom_hpc"
    # execution_type = "parallel_study"

    if execution_type == "default_example":
        example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
        if not os.path.exists(example_results_folder):
            os.mkdir(example_results_folder)
        parallel_results_folder = os.path.join(example_results_folder, "parallel")
        if not os.path.exists(parallel_results_folder):
            os.mkdir(parallel_results_folder)
        working_directory = os.path.join(parallel_results_folder, "default")
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        number_of_models = 10
        number_of_processes = 5

        inductor_frequency = 270000

        geos = []
        simulation_parameters = []
        working_directories = []
        for i in range(number_of_models):
            geos.append(create_parallel_example_inductor(inductor_frequency))
            simulation_parameters.append({
                "freq": inductor_frequency,
                "current": [4.5*(1+i/10)]
            })

        start_time = time.time()
        fmt.run_hpc(number_of_processes, geos, simulation_parameters, working_directory)
        execution_time = time.time() - start_time

        print(f"Execution time: {execution_time}")

    elif execution_type == "custom_hpc":
        example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
        if not os.path.exists(example_results_folder):
            os.mkdir(example_results_folder)
        parallel_results_folder = os.path.join(example_results_folder, "parallel")
        if not os.path.exists(parallel_results_folder):
            os.mkdir(parallel_results_folder)
        working_directory = os.path.join(parallel_results_folder, "default")
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        number_of_models = 8
        number_of_processes = 4

        inductor_frequency = 270000

        geos = []
        simulation_parameters = []
        working_directories = []
        for i in range(number_of_models):
            geos.append(create_parallel_example_inductor(inductor_frequency))
            simulation_parameters.append({
                "current": [4, 12],
                "phi_deg": [0, 180]
            })

        start_time = time.time()
        fmt.run_hpc(number_of_processes, geos, simulation_parameters, working_directory, custom_hpc)
        execution_time = time.time() - start_time

        print(f"Execution time: {execution_time}")

    elif execution_type == "parallel_study":
        example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
        if not os.path.exists(example_results_folder):
            os.mkdir(example_results_folder)
        parallel_results_folder = os.path.join(example_results_folder, "parallel")
        if not os.path.exists(parallel_results_folder):
            os.mkdir(parallel_results_folder)
        working_directory = os.path.join(parallel_results_folder, "default")
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        parallel_simulation_study(3)
    else:
        raise Exception(f"Execution type {execution_type} not found.")
