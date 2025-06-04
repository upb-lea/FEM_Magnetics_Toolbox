"""
Basic example to show how to simulate an inductor with electrostatic analysis.

After starting the program, the geometry dimensions are displayed. Verify this geometry, close the window to continue the simulation.
This example simulates the electrostatic properties of an inductor, such as capacitance and electric fields.

Once the geometry is verified, an electrostatic simulation will be run with voltages applied to the turns of the winding.
The simulation results will include electrostatic potential distributions, electric field maps, and capacitance data.
These results can help you understand how the electric field is distributed within the component and where high-field regions may occur.

The simulation results will be visualized. In these visual outputs, you will be able to see the distribution of electrostatic potential in different turns of
the winding and the influence of the core and other materials in the geometry.
"""
import femmt as fmt
import os


def basic_example_inductor_electrostatic(onelab_folder: str = None, show_visual_outputs: bool = True, is_test: bool = False):
    """
    Run the example code for the inductor.

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
    geo = fmt.MagneticComponent(simulation_type=fmt.SimulationType.ElectroStatic, component_type=fmt.ComponentType.Inductor,
                                working_directory=working_directory, is_gui=is_test)

    # This line is for automated pytest running on GitHub only. Please ignore this line!
    if onelab_folder is not None:
        geo.file_data.onelab_folder_path = onelab_folder

    inductor_frequency = 2700000

    # 2. set core parameters
    core_db = fmt.core_database()["PQ 40/40"]
    core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                    window_w=core_db["window_w"],
                                                    window_h=core_db["window_h"],
                                                    core_h=core_db["core_h"])

    core = fmt.Core(core_type=fmt.CoreType.Single,
                    core_dimensions=core_dimensions,
                    detailed_core_model=False,
                    material=fmt.Material.N49, temperature=45, frequency=inductor_frequency,
                    # permeability_datasource="manufacturer_datasheet",
                    permeability_datasource=fmt.MaterialDataSource.Measurement,
                    permeability_datatype=fmt.MeasurementDataType.ComplexPermeability,
                    permeability_measurement_setup=fmt.MeasurementSetup.LEA_LK,
                    permittivity_datasource=fmt.MaterialDataSource.Measurement,
                    permittivity_datatype=fmt.MeasurementDataType.ComplexPermittivity,
                    permittivity_measurement_setup=fmt.MeasurementSetup.LEA_LK, mdb_verbosity=fmt.Verbosity.Silent)

    geo.set_core(core)
    # 3. set air gap parameters
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.001, 50)
    # air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0002, 90)
    geo.set_air_gaps(air_gaps)

    # 4. set insulations
    bobbin_db = fmt.bobbin_database()["PQ 40/40"]
    bobbin_dimensions = fmt.dtos.BobbinDimensions(bobbin_inner_diameter=bobbin_db["bobbin_inner_diameter"],
                                                  bobbin_window_w=bobbin_db["bobbin_window_w"],
                                                  bobbin_window_h=bobbin_db["bobbin_window_h"],
                                                  bobbin_h=bobbin_db["bobbin_h"])
    insulation = fmt.Insulation(flag_insulation=True, bobbin_dimensions=bobbin_dimensions)
    insulation.add_core_insulations(1.5e-3, 1.5e-3, 1e-3, 1e-3)
    insulation.add_turn_insulation([0.5e-4], add_turn_insulations=True)
    insulation.add_winding_insulations([[0.6e-3, 0.1e-3]], per_layer_of_turns=True)
    insulation.add_kapton_insulation(add_kapton_material=True, thickness=0.6e-3)
    geo.set_insulation(insulation)

    # 5. create winding window and virtual winding windows (vww)
    winding_window = fmt.WindingWindow(core, insulation)
    vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

    # 6. create conductor and set parameters: use solid wires
    winding = fmt.Conductor(0, fmt.Conductivity.Copper, winding_material_temperature=45)
    winding.set_solid_round_conductor(conductor_radius=1.1506e-3, conductor_arrangement=fmt.ConductorArrangement.Square)
    winding.parallel = False  # set True to make the windings parallel, currently only for solid conductors
    # 7. add conductor to vww and add winding window to MagneticComponent
    vww.set_winding(winding, 14, None, fmt.Align.ToEdges, placing_strategy=fmt.ConductorDistribution.VerticalUpward_HorizontalRightward,
                    zigzag=True)
    geo.set_winding_windows([winding_window])
    # 8. create the model
    geo.create_model(freq=inductor_frequency, pre_visualize_geometry=show_visual_outputs, save_png=False, skin_mesh_factor=0.5)
    # 8. run electrostatic simulation
    num_turns_w1 = 14
    # Create a linear voltage distribution along winding 1 from V_A to V_B ( the first turn to the last turn)
    V_A = 1
    V_B = 0
    voltages_winding_1 = [
        V_A - (V_A - V_B) * i / (num_turns_w1 - 1)
        for i in range(num_turns_w1)
    ]
    # geo.electrostatic_simulation(voltage=[voltages_winding_1], ground_outer_boundary=False, core_voltage=0,
    #                              show_fem_simulation_results=show_visual_outputs, save_to_excel=False)
    # Run simulation in FEMM
    geo.femm_reference_electrostatic(voltages=[voltages_winding_1], ground_core=True, ground_outer_boundary=True,
                                     non_visualize=0, save_to_excel=False, compare_excel_files_to_femmt=False, mesh_size_conductor=0.0)


if __name__ == "__main__":
    basic_example_inductor_electrostatic(show_visual_outputs=True)
