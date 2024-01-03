"""Advanced example to demonstrate an air gap sweep for an inductor."""
import matplotlib.pyplot as plt
import femmt as fmt
import os
from typing import Optional

if not os.path.exists(os.path.join(os.path.dirname(__file__), "sweep_examples")):
    os.mkdir(os.path.join(os.path.dirname(__file__), "sweep_examples"))


def basic_example_sweep(onelab_folder: Optional[str] = None, show_visual_outputs: bool = True, is_test: bool = False):
    """Advanced example to demonstrate an air gap sweep for an inductor."""
    # In this sweep an inductor with variable air gap height is simulated
    air_gap_heights = [0.0020, 0.0010, 0.0005, 0.000025, 0.0000125]

    working_directories = []

    example_results_folder = os.path.join(os.path.dirname(__file__), "example_results", "sweep_examples")
    if not os.path.exists(example_results_folder):
        os.mkdir(example_results_folder)

    for i, height in enumerate(air_gap_heights):
        # In order to save the results for every simulation the working directory is changed
        # Working directory can be set arbitrarily
        directory = os.path.join(example_results_folder, f"air_gap_{i}")
        if not directory:
            os.mkdir(directory)

        working_directories.append(directory)

        geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=directory, verbosity=fmt.Verbosity.Silent, is_gui=is_test)

        # This line is for automated pytest running on GitHub only. Please ignore this line!
        if onelab_folder is not None:
            geo.file_data.onelab_folder_path = onelab_folder

        core_db = fmt.core_database()["PQ 40/40"]
        core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"], window_w=core_db["window_w"],
                                                        window_h=core_db["window_h"], core_h=core_db["core_h"])
        core = fmt.Core(core_dimensions=core_dimensions,
                        material=fmt.Material.N95, temperature=25, frequency=100000,
                        permeability_datasource=fmt.MaterialDataSource.Measurement,
                        permeability_datatype=fmt.MeasurementDataType.ComplexPermeability,
                        permeability_measurement_setup=fmt.MeasurementSetup.LEA_LK,
                        permittivity_datasource=fmt.MaterialDataSource.Measurement,
                        permittivity_datatype=fmt.MeasurementDataType.ComplexPermittivity,
                        permittivity_measurement_setup=fmt.MeasurementSetup.LEA_LK
                        )
        geo.set_core(core)

        air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, height, 50)
        geo.set_air_gaps(air_gaps)

        insulation = fmt.Insulation()
        insulation.add_core_insulations(0.001, 0.001, 0.004, 0.001)
        insulation.add_winding_insulations([[0.0005]])
        geo.set_insulation(insulation)

        winding_window = fmt.WindingWindow(core, insulation)
        complete = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

        conductor = fmt.Conductor(0, fmt.Conductivity.Copper)
        conductor.set_litz_round_conductor(conductor_radius=0.0013, number_strands=150, strand_radius=100e-6, 
                                           fill_factor=None, conductor_arrangement=fmt.ConductorArrangement.Square)
        complete.set_winding(conductor, 9, None)
        geo.set_winding_windows([winding_window])

        geo.create_model(freq=100000, pre_visualize_geometry=False, save_png=False)
        geo.single_simulation(freq=100000, current=[4.5], show_fem_simulation_results=False)

    # After the simulations the sweep can be analyzed
    # This could be done using the FEMMTLogParser:
    logs = fmt.FEMMTLogParser.get_log_files_from_working_directories(working_directories)
    log_parser = fmt.FEMMTLogParser(logs)

    # In this case the self inductivity of winding1 will be analyzed
    inductivities = []
    for _, data in log_parser.data.items():
        inductivities.append(data.sweeps[0].windings[0].flux_over_current)

    if not is_test:
        plt.plot(air_gap_heights, inductivities, "ro")
        plt.title("Air gap height sweep")
        plt.xlabel("air gap height")
        plt.ylabel("self inductance")
        plt.grid()
        plt.show()


if __name__ == '__main__':
    basic_example_sweep()
