"""Advanced example to show an excitation sweep of multiple frequencies for an inductor."""
import numpy as np

import femmt as fmt
import materialdatabase as mdb
import os


def advanced_example_inductor_sweep(onelab_folder: str = None, show_visual_outputs: bool = True, is_test: bool = False):
    """
    Advanced example to show an excitation sweep of multiple frequencies for an inductor.

    :param onelab_folder: onelab folder path
    :type onelab_folder: str
    :param show_visual_outputs: True to show visual outputs (simulation results)
    :type show_visual_outputs: bool
    :param is_test: True for pytest usage. Defaults to False.
    :type is_test: bool
    """
    example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
    if not os.path.exists(example_results_folder):
        os.mkdir(example_results_folder)

    # Working directory can be set arbitrarily
    working_directory = os.path.join(example_results_folder, os.path.splitext(os.path.basename(__file__))[0])
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    # 1. chose simulation type
    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=working_directory, is_gui=is_test)

    # This line is for automated pytest running on GitHub only. Please ignore this line!
    if onelab_folder is not None:
        geo.file_data.onelab_folder_path = onelab_folder

    inductor_frequency = 5e5

    # 2. set core parameters
    core_db = fmt.core_database()["PQ 50/50"]
    core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                    window_w=core_db["window_w"],
                                                    window_h=core_db["window_h"],
                                                    core_h=core_db["core_h"])

    core = fmt.Core(core_type=fmt.CoreType.Single,
                    core_dimensions=core_dimensions,
                    material=mdb.Material.N95, temperature=30,
                    permeability_datasource=fmt.MaterialDataSource.Measurement,
                    permeability_datatype=fmt.MeasurementDataType.ComplexPermeability,
                    permeability_measurement_setup=mdb.MeasurementSetup.LEA_LK,
                    permittivity_datasource=fmt.MaterialDataSource.Measurement,
                    permittivity_datatype=fmt.MeasurementDataType.ComplexPermittivity,
                    permittivity_measurement_setup=mdb.MeasurementSetup.LEA_LK)

    geo.set_core(core)

    # 3. set air gap parameters
    # air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    # air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
    # air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0002, 90)
    # geo.set_air_gaps(air_gaps)

    # 4. set insulations
    insulation = fmt.Insulation(flag_insulation=False)
    insulation.add_core_insulations(0.001, 0.014, 0.006, 0.001)
    insulation.add_winding_insulations([[0.0005]])
    geo.set_insulation(insulation)

    # 5. create winding window and virtual winding windows (vww)
    winding_window = fmt.WindingWindow(core, insulation)
    vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

    # 6. create conductor and set parameters: use solid wires
    winding = fmt.Conductor(0, fmt.ConductorMaterial.Copper, temperature=30)
    winding.set_solid_round_conductor(conductor_radius=0.0013, conductor_arrangement=fmt.ConductorArrangement.Square)
    # winding.parallel = False  # set True to make the windings parallel, currently only for solid conductors
    # winding.set_litz_round_conductor(conductor_radius=0.0013, number_strands=150, strand_radius=100e-6,
    #                                  fill_factor=None, conductor_arrangement=fmt.ConductorArrangement.Square)

    # 7. add conductor to vww and add winding window to MagneticComponent
    vww.set_winding(winding, 1, None, fmt.Align.ToEdges, fmt.ConductorDistribution.VerticalUpward_HorizontalRightward, zigzag=False)
    geo.set_winding_windows([winding_window])

    # 8. create the model
    geo.create_model(freq=inductor_frequency, pre_visualize_geometry=show_visual_outputs, save_png=False)

    # 6.a. start simulation
    # geo.single_simulation(freq=inductor_frequency, current=[0.01],
    #                       plot_interpolation=True, show_fem_simulation_results=show_visual_outputs)

    # geo.femm_reference(freq=inductor_frequency, current=[4.5], sign=[1], non_visualize=0)

    # # 6.b. Excitation Sweep Example
    # # Perform a sweep using more than one frequency
    fs = list(np.linspace(100e3, 500e3, 5))
    amplitude_list = [[0.01] for _ in range(5)]
    # currents = np.linspace(0.5, 1/30, 10)
    # amplitude_list = [[x] for x in currents]
    phase_list = [[0] for _ in range(5)]
    geo.excitation_sweep(frequency_list=fs, current_list_list=amplitude_list, phi_deg_list_list=phase_list)


if __name__ == "__main__":
    advanced_example_inductor_sweep(show_visual_outputs=True)
