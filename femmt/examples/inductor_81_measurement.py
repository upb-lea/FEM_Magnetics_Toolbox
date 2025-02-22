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


def basic_example_inductor_measurement(onelab_folder: str = None, show_visual_outputs: bool = True, is_test: bool = False):
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
    geo = fmt.MagneticComponent(simulation_type=fmt.SimulationType.ElectroStatic, component_type=fmt.ComponentType.Inductor, working_directory=working_directory,
                                verbosity=fmt.Verbosity.ToConsole, is_gui=is_test)

    # This line is for automated pytest running on GitHub only. Please ignore this line!
    if onelab_folder is not None:
        geo.file_data.onelab_folder_path = onelab_folder

    inductor_frequency = 100000

    # 2. set core parameters
    # core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=0.02, window_w=0.01, window_h=0.03,
    #                                                 core_h=0.02)
    # core = fmt.Core(core_dimensions=core_dimensions, mu_r_abs=3100, phi_mu_deg=12, sigma=1.2,
    #                 permeability_datasource=fmt.MaterialDataSource.Custom,
    #                 permittivity_datasource=fmt.MaterialDataSource.Custom,
    #                 detailed_core_model=False)
    # geo.set_core(core)
    # 2. set core parameters
    core_db = fmt.core_database()["PQ 40/40"]
    core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                    window_w=core_db["window_w"],
                                                    window_h=core_db["window_h"],
                                                    core_h=core_db["core_h"])
    bobbin_db = fmt.bobbin_database()["PQ 40/40"]
    bobbin_dimensions = fmt.dtos.BobbinDimensions(bobbin_inner_diameter=bobbin_db["bobbin_inner_diameter"],
                                                  bobbin_window_w=bobbin_db["bobbin_window_w"],
                                                  bobbin_window_h=bobbin_db["bobbin_window_h"],
                                                  bobbin_h=bobbin_db["bobbin_h"])

    core = fmt.Core(core_type=fmt.CoreType.Single,
                    core_dimensions=core_dimensions,
                    bobbin_dimensions=bobbin_dimensions,
                    detailed_core_model=False,
                    material=fmt.Material.N95, temperature=45, frequency=inductor_frequency,
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
    insulation = fmt.Insulation(flag_insulation=True)
    insulation.add_core_insulations(1.55e-3, 1.55e-3, 0.9e-3, 1.5e-4)
    insulation.add_winding_insulations([[0.025e-3]])
    insulation.add_conductor_air_conductor_insulation([[1.0837e-4, 1.54e-4, 1.9285e-4, 2.2777e-4],
                                                       [3.06e-4, 2.653e-4]])
    insulation.add_kapton_insulation(add_kapton=True, thickness=0.07e-3)
    geo.set_insulation(insulation)

    # 5. create winding window and virtual winding windows (vww)
    winding_window = fmt.WindingWindow(core, insulation)
    vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

    # 6. create conductor and set parameters: use solid wires
    winding = fmt.Conductor(0, fmt.Conductivity.Copper, winding_material_temperature=45)
    # winding.set_solid_round_conductor(conductor_radius=1.1506e-3, conductor_arrangement=fmt.ConductorArrangement.Square)
    winding.set_solid_round_conductor(conductor_radius=0.225e-3, conductor_arrangement=fmt.ConductorArrangement.Square)
    winding.parallel = False  # set True to make the windings parallel, currently only for solid conductors
    # winding.set_litz_round_conductor(conductor_radius=0.0013, number_strands=150, strand_radius=100e-6,
    # fill_factor=None, conductor_arrangement=fmt.ConductorArrangement.Square)
    # 7. add conductor to vww and add winding window to MagneticComponent
    vww.set_winding(winding, 81, None, fmt.Align.ToEdges, placing_strategy=fmt.ConductorDistribution.VerticalUpward_HorizontalRightward,
                    zigzag=True)
    geo.set_winding_windows([winding_window])
    num_turns_w1 = 81
    # Create a linear voltage distribution along winding 1 from V_A to V_B
    V_A = 1
    V_B = 0
    voltages_winding_1 = [
        V_A - (V_A - V_B) * i / (num_turns_w1 - 1)
        for i in range(num_turns_w1)
    ]
    # voltages_winding_1 = [0] * 80 + [1] * 1
    # 8. create the model
    geo.create_model(freq=inductor_frequency, pre_visualize_geometry=show_visual_outputs, save_png=False, skin_mesh_factor=0.5)
    # 8. run electrostatic simulation
    geo.electrostatic_simulation(voltage=[voltages_winding_1], core_voltage=None, ground_outer_boundary=False,
                                 show_fem_simulation_results=show_visual_outputs, save_to_excel=False)
    geo.get_total_charges()
    # Call the electrostatic FEMM simulation function
    # voltages = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140]
    geo.femm_reference_electrostatic(voltages=[voltages_winding_1], ground_core=True, ground_outer_boundary=False,
                                     non_visualize=0, save_to_excel=False, compare_excel_files_to_femmt=False, mesh_size_conductor=0.0)


if __name__ == "__main__":
    basic_example_inductor_measurement(show_visual_outputs=True)
