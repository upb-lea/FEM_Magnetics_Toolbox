"""Examples for the parallel simulation."""
# Python standard libraries
from typing import Dict
from itertools import product
import os
import time
import matplotlib.pyplot as plt
import statistics

# File to generate Data
import Testdata_Generator

# Local libraries
import femmt as fmt
import materialdatabase as mdb

"""This file contains examples for the use of the hpc function. Internally the multprocessing package is used."""

def create_parallel_example_transformer() -> fmt.MagneticComponent:
    """Creates an example model which is used for the parallel execution example. This does implement a simple transformer.
    """
    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Transformer, working_directory=working_directory, verbosity=fmt.Verbosity.ToFile)
    core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=Testdata_Generator.corediameter(), window_w=Testdata_Generator.windowwidth(), window_h=Testdata_Generator.windowheight(), core_h= Testdata_Generator.coreheight())
    core = fmt.Core(core_dimensions=core_dimensions, non_linear=False, correct_outer_leg=True,
                    material=mdb.Material.N95, temperature=Testdata_Generator.temperatures(), frequency=Testdata_Generator.frequencies(),
                    permeability_datasource=fmt.MaterialDataSource.Measurement,
                    permeability_datatype=fmt.MeasurementDataType.ComplexPermeability,
                    permeability_measurement_setup=mdb.MeasurementSetup.LEA_LK,
                    permittivity_datasource=fmt.MaterialDataSource.Measurement,
                    permittivity_datatype=fmt.MeasurementDataType.ComplexPermittivity,
                    permittivity_measurement_setup=mdb.MeasurementSetup.LEA_LK, mdb_verbosity=fmt.Verbosity.Silent)

    geo.set_core(core)
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, Testdata_Generator.airgapheight(),
                         Testdata_Generator.airgapposition(1))

    geo.set_air_gaps(air_gaps)
    insulation = fmt.Insulation()
    insulation.add_core_insulations(Testdata_Generator.coreinsulation())
    insulation.add_winding_insulations(Testdata_Generator.innerwindinginsulation())
    geo.set_insulation(insulation)
    winding_window = fmt.WindingWindow(core, insulation)
    vww = winding_window.split_window(Testdata_Generator.windingwindowsplit())
    winding1 = fmt.Conductor(0, fmt.Conductivity.Copper)
    winding1.set_solid_round_conductor(conductor_radius=Testdata_Generator.conductorradius(),
                                      conductor_arrangement=None)
    winding2 = fmt.Conductor(1, fmt.Conductivity.Copper)
    winding2.set_solid_round_conductor(conductor_radius=Testdata_Generator.conductorradius(),
                                      conductor_arrangement=None )
    vww.set_interleaved_winding(winding1, Testdata_Generator.windingturns(1), winding2, Testdata_Generator.windingturns(2))
    geo.set_winding_windows([winding_window])

    return geo


def create_parallel_example_inductor(inductor_frequency: int, temperature: int, turns: int, number_of_air_gaps: int, air_gap_position: list) -> fmt.MagneticComponent:
    # Creates an example model which is used for the parallel execution example. This does implement a simple inductor with given inductor_frequency.

    #:param inductor_frequency: Frequency for the inductor.
    #:type inductor_frequency: int

    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=None,
                                # Can be set to None since it will be overwritten anyways
                                clean_previous_results=False, verbosity=fmt.Verbosity.Silent)
    core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=Testdata_Generator.corediameter(),
                                                    window_w= Testdata_Generator.windowwidth(),
                                                    window_h= Testdata_Generator.windowheight(),
                                                    core_h= Testdata_Generator.coreheight())
    core = fmt.Core(fmt.CoreType.Single, core_dimensions=core_dimensions,
                    material=mdb.Material.custom_material, temperature= temperature, non_linear=False, correct_outer_leg=True,
                    mu_r_abs=3100, phi_mu_deg=12, sigma=1.2,
                    permeability_datasource=fmt.MaterialDataSource.Custom,
                    permittivity_datasource=fmt.MaterialDataSource.Custom,
                    frequency=Testdata_Generator.frequencies(),
                    permeability_datatype=fmt.MeasurementDataType.ComplexPermeability,
                    permeability_measurement_setup=mdb.MeasurementSetup.LEA_LK,
                    permittivity_datatype=fmt.MeasurementDataType.ComplexPermittivity,
                    permittivity_measurement_setup=mdb.MeasurementSetup.LEA_LK, mdb_verbosity=fmt.Verbosity.Silent)

    geo.set_core(core)
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    # sets first air gap
    if number_of_air_gaps in [1, 2, 3]:
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, Testdata_Generator.airgapheight(), air_gap_position[0])
    # sets second air gap
    if number_of_air_gaps in [2, 3]:
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, Testdata_Generator.airgapheight(), air_gap_position[1])
    # sets third air gap
    if number_of_air_gaps in [3]:
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, Testdata_Generator.airgapheight(), air_gap_position[2])
    geo.set_air_gaps(air_gaps)
    insulation = fmt.Insulation()
    insulation.add_core_insulations(Testdata_Generator.coreinsulation()[0],Testdata_Generator.coreinsulation()[1],
                                    Testdata_Generator.coreinsulation()[2],Testdata_Generator.coreinsulation()[3])
    insulation.add_winding_insulations([[Testdata_Generator.innerwindinginsulation()]])
    geo.set_insulation(insulation)
    winding_window = fmt.WindingWindow(core, insulation)
    vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)
    winding = fmt.Conductor(0, fmt.Conductivity.Copper)
    winding.set_solid_round_conductor(conductor_radius=Testdata_Generator.conductorradius(), conductor_arrangement=Testdata_Generator.conductorarrangement())
    winding.parallel = False
    vww.set_winding(winding, turns, None)
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

    model.create_model(freq=Testdata_Generator.frequencies(), pre_visualize_geometry=False)
    model.single_simulation(freq=Testdata_Generator.frequencies(), current=current, phi_deg=phi_deg, show_fem_simulation_results=False)

def parallel_simulation_study(averaging_count):
    """Perform several parallel simulations."""
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
                "current": [1]
    })

    for process_count in process_counts:
        working_directory = os.path.join(study_results_folder, f"{process_count}")
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        simulation_times = []
        for _ in range(averaging_count):
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
        plt.ylabel("Runtime")
    plt.show()


if __name__ == "__main__":

    # ---- Instances for Simulation loop ----
    Model_instance = Testdata_Generator.Model_Values
    Results_instance = Testdata_Generator.ChangeFolder

    # ---- Preparing for Simulation loop ----
    winding_turns = 1
    print("1")
    Model_instance.next(Testdata_Generator.Model_Values)
    print("3")
    # ---- Choosing the execution ----
    execution_type = "default_example"
    #execution_type = "custom_hpc"
    #execution_type = "parallel_study"


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

        inductor_frequencies = Testdata_Generator.frequencies()
        temperatures = Testdata_Generator.temperatures()
        number_of_air_gaps, air_gap_position = Testdata_Generator.airgapnumber_and_position()

        # This while loop checks if the windings fit into the winding window and if the magnetic table flow is too high. Every Iteration adds one winding turn
        while (Testdata_Generator.check_windings_fit(winding_turns, Testdata_Generator.conductorradius(), Testdata_Generator.coreinsulation()[0],Testdata_Generator.coreinsulation()[1],
                                                     Testdata_Generator.coreinsulation()[2],Testdata_Generator.coreinsulation()[3], Testdata_Generator.windowheight(), Testdata_Generator.windowwidth())
               and Testdata_Generator.is_magnetic_flux_density_below_limit(winding_turns, Testdata_Generator.coreheight(), Testdata_Generator.current())):

            print("Winding:", winding_turns)



            combinations = list(product(inductor_frequencies, temperatures))

            number_of_models = len(combinations)
            number_of_processes = 5

            geos = []
            simulation_parameters = []
            working_directories = []

            for inductor_frequency, temperature in combinations:
                    geos.append(create_parallel_example_inductor(inductor_frequency, temperature, winding_turns, number_of_air_gaps, air_gap_position))
                    simulation_parameters.append({
                        "freq": inductor_frequency,
                        "current": [Testdata_Generator.current()]
                    })

            start_time = time.time()

            fmt.run_hpc(number_of_processes, geos, simulation_parameters, working_directory)
            execution_time = time.time() - start_time

            Results_instance.copy_result_to_new_folder(Testdata_Generator.ChangeFolder, number_of_models)

            print(f"Execution time: {execution_time}")

            winding_turns += 1


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

        inductor_frequency = Testdata_Generator.frequencies()

        geos = []
        simulation_parameters = []
        working_directories = []
        for i in range(number_of_models):
            geos.append(create_parallel_example_inductor(inductor_frequency, 40,10))
            simulation_parameters.append({
                "current": Testdata_Generator.current(),
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
