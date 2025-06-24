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

def basic_example_transformer_electrostatic(onelab_folder: str = None, show_visual_outputs: bool = True, is_test: bool = False):
    """
    Run the example code for the transformer.

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

    # Example for a transformer with multiple virtual winding windows.
    working_directory = os.path.join(example_results_folder, os.path.splitext(os.path.basename(__file__))[0])
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    # 1. chose simulation type
    geo = fmt.MagneticComponent(simulation_type=fmt.SimulationType.ElectroStatic, component_type=fmt.ComponentType.Transformer,
                                working_directory=working_directory, is_gui=is_test)

    # This line is for automated pytest running on GitHub only. Please ignore this line!
    if onelab_folder is not None:
        geo.file_data.onelab_folder_path = onelab_folder

    # 2. set core parameters
    # core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=0.02, window_w=0.01, window_h=0.03,
    #                                                 core_h=0.02)
    # core = fmt.Core(core_dimensions=core_dimensions, mu_r_abs=3100, phi_mu_deg=12, sigma=1.2,
    #                 permeability_datasource=fmt.MaterialDataSource.Custom,
    #                 permittivity_datasource=fmt.MaterialDataSource.Custom,
    #                 detailed_core_model=False)
    # geo.set_core(core)
    core_db = fmt.core_database()["PQ 40/40"]
    core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                   window_w=core_db["window_w"],
                                                   window_h=core_db["window_h"],
                                                   core_h=core_db["core_h"])
    core = fmt.Core(core_type=fmt.CoreType.Single,
                    core_dimensions=core_dimensions,
                    detailed_core_model=False,
                    material=fmt.Material.N49, temperature=45, frequency=0,
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
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
    geo.set_air_gaps(air_gaps)

    # 4. set insulation

    bobbin_db = fmt.bobbin_database()["PQ 40/40"]
    bobbin_dimensions = fmt.dtos.BobbinDimensions(bobbin_inner_diameter=bobbin_db["bobbin_inner_diameter"],
                                                  bobbin_window_w=bobbin_db["bobbin_window_w"],
                                                  bobbin_window_h=bobbin_db["bobbin_window_h"],
                                                  bobbin_h=bobbin_db["bobbin_h"])
    insulation = fmt.Insulation(flag_insulation=True, bobbin_dimensions=bobbin_dimensions)
    bobbin_material = fmt.insulation_materials_database()["core_insulation"]["bobbins"]["Thermoset"]["Phenolic"]
    insulation.add_core_insulations(1.55e-3, 1.55e-3, 0.9e-3, 1.5e-4,
                                    dielectric_constant=bobbin_material["dielectric_constant"])
    insulation.add_winding_insulations([[0.0002, 0.095e-3],
                                        [0.095e-3, 0.0002]], per_layer_of_turns=True)
    turn_insulation_material = fmt.insulation_materials_database()["wire_insulation"]["plastic_insulation"]["Plenum Polyvinyl Chloride (Plenum PVC)"]
    insulation.add_turn_insulation([0.25e-5, 0.25e-5],
                                   dielectric_constant=[turn_insulation_material["dielectric_constant"], turn_insulation_material["dielectric_constant"]],
                                   add_turn_insulations=False)
    layer_insulation = fmt.insulation_materials_database()["film_insulation"]["Kapton"]
    insulation.add_insulation_between_layers(add_insulation_material=True, thickness=0.5e-3, dielectric_constant=layer_insulation["dielectric_constant"])
    geo.set_insulation(insulation)

    # 5. create winding window and virtual winding windows (vww)
    winding_window = fmt.WindingWindow(core, insulation)
    # bot, top = winding_window.split_window(fmt.WindingWindowSplit.HorizontalSplit, split_distance=0.001)
    # 109-49
    cells = winding_window.NHorizontalAndVerticalSplit(horizontal_split_factors=[0.29],
                                                       vertical_split_factors=None)
    # cells = winding_window.NHorizontalAndVerticalSplit(horizontal_split_factors=[0.],
    #                                                    vertical_split_factors=None)

    # 6. create conductors and set parameters
    winding1 = fmt.Conductor(0, fmt.Conductivity.Copper)
    # winding1.set_solid_round_conductor(1.1506e-3, fmt.ConductorArrangement.Square)
    #winding1.set_solid_round_conductor(0.35e-3, fmt.ConductorArrangement.Square)
    winding1.set_solid_round_conductor(1.1506e-3, fmt.ConductorArrangement.Square)

    # winding1 = fmt.Conductor(0, fmt.Conductivity.Copper)
    # winding1.set_litz_round_conductor(0.0011, 50, 0.00011, None, fmt.ConductorArrangement.Square)

    # winding2 = fmt.Conductor(1, fmt.Conductivity.Copper)
    # winding2.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

    winding2 = fmt.Conductor(1, fmt.Conductivity.Copper)
    # winding2.set_solid_round_conductor(1.1506e-3, fmt.ConductorArrangement.Square)
    winding2.set_solid_round_conductor(1.1506e-3, fmt.ConductorArrangement.Square)
    winding2.parallel = False
    # winding2.set_litz_round_conductor(0.0011, 50, 0.00011, None, fmt.ConductorArrangement.Square)

    # 7. add conductor to vww and add winding window to MagneticComponent
    # top.set_winding(winding2, 15, None, fmt.Align.ToEdges, fmt.ConductorDistribution.VerticalUpward_HorizontalRightward, zigzag=True)
    cells[1].set_winding(winding2, 10, None, fmt.Align.ToEdges, fmt.ConductorDistribution.VerticalUpward_HorizontalRightward, zigzag=True)
    # bot.set_winding(winding2, 29, None, fmt.Align.ToEdges, fmt.ConductorDistribution.VerticalUpward_HorizontalRightward, zigzag=True)
    cells[0].set_winding(winding1, 10, None, fmt.Align.ToEdges, fmt.ConductorDistribution.VerticalUpward_HorizontalRightward, zigzag=True)
    # top.set_winding(winding1, 109, None, fmt.Align.ToEdges, fmt.ConductorDistribution.VerticalUpward_HorizontalRightward, zigzag=True)
    geo.set_winding_windows([winding_window])



    # 8. start simulation with given frequency, currents and phases
    geo.create_model(freq=0, pre_visualize_geometry=show_visual_outputs)
    # Number of turns in each winding
    num_turns_w1 = 10
    num_turns_w2 = 10

    """---------------------------------The definition of the 10 simulations needed for extracting the capacitance of the transformer------------------- """
    # Note that the user can run these 10 simulations by simply calling the get_capacitance_of_transformer() function
    # and capacitance extraction is shown in the terminal

    # --------------------------------------------------------------------------------------
    # Simulation 1 (V_A, V_B, V_C, V_D = 1, 0, 0 , 0) --- (V_1, V_2, V_3, V_4 = 1, 0, 0, 0)
    V_A = 1
    V_B = 0
    voltages_winding_1 = [
        V_A - (V_A - V_B) * i / (num_turns_w1 - 1)
        for i in range(num_turns_w1)
    ]
    voltages_winding_2 = [1] * num_turns_w2

    # --------------------------------------------------------------------------------------
    # # Simulation 2 (V_A, V_B, V_C, V_D = 0, 0, 1 , 0) --- (V_1, V_2, V_3, V_4 = 0, 1, 0, 0)
    # voltages_winding_1 = [0] * num_turns_w1
    # V_C = 1.0
    # V_D = 0.0
    # voltages_winding_2 = [
    #     V_C - (V_C - V_D) * i / (num_turns_w2 - 1)
    #     for i in range(num_turns_w2)
    # ]

    # ---------------------------------------------------------------------------------------------
    # # Simulation 3 (V_A, V_B, V_C, V_D = 0, 0, 1 , 1) --- (V_1, V_2, V_3, V_4 = 0, 0, 1, 0)
    # voltages_winding_1 = [0] * num_turns_w1
    # voltages_winding_2 = [1] * num_turns_w2

    # -----------------------------------------------------------------------------------------------
    # # Simulation 4 (V_1, V_2, V_3, V_4 = 1, 1, 1 , 1) --- (V_1, V_2, V_3, V_4 = 0, 0, 0, 1)
    # voltages_winding_1 = [1] * num_turns_w1
    # voltages_winding_2 = [1] * num_turns_w2

    # -----------------------------------------------------------------------------------------------
    # # Simulation 5 (V_A, V_B, V_C, V_D = 1, 0, 1 , 0) --- (V_1, V_2, V_3, V_4 = 1, 1, 0, 0)
    # V_A = 1.0
    # V_B = 0.0
    # V_C = 1.0
    # V_D = 0.0
    # voltages_winding_1 = [
    #     V_A - (V_A - V_B) * i / (num_turns_w1 - 1)
    #     for i in range(num_turns_w1)
    # ]
    # voltages_winding_2 = [
    #     V_C - (V_C - V_D) * i / (num_turns_w2 - 1)
    #     for i in range(num_turns_w2)
    # ]

    # ------------------------------------------------------------------------------------------------
    # # Simulation 6 (V_A, V_B, V_C, V_D = 1, 0, 1 , 1) --- (V_1, V_2, V_3, V_4 = 1, 0, 1, 0)
    # V_A = 1.0
    # V_B = 0.0
    # voltages_winding_1 = [
    #     V_A - (V_A - V_B) * i / (num_turns_w1 - 1)
    #     for i in range(num_turns_w1)
    # ]
    # voltages_winding_2 = [1] * num_turns_w2

    # --------------------------------------------------------------------------------------------------
    # # Simulation 7 (V_A, V_B, V_C, V_D = 2, 1, 1 , 1) --- (V_1, V_2, V_3, V_4 = 1, 0, 0, 1)
    # V_A = 2.0
    # V_B = 1.0
    # voltages_winding_1 = [
    #     V_A - (V_A - V_B) * i / (num_turns_w1 - 1)
    #     for i in range(num_turns_w1)
    # ]
    # voltages_winding_2 = [1] * num_turns_w2

    # ------------------------------------------------------------------------------------------------------
    # # Simulation 8 (V_A, V_B, V_C, V_D = 0, 0, 2 , 1) --- (V_1, V_2, V_3, V_4 = 0, 1, 1, 1)
    # voltages_winding_1 = [0] * num_turns_w1
    # V_C = 2.0
    # V_D = 1.0
    # voltages_winding_2 = [
    #     V_C - (V_C - V_D) * i / (num_turns_w2 - 1)
    #     for i in range(num_turns_w2)
    # ]

    # ---------------------------------------------------------------
    # # Simulation 9 (V_A, V_B, V_C, V_D = 1, 1, 2 , 1) --- (V_1, V_2, V_3, V_4 = 0, 1, 0, 1)
    # voltages_winding_1 = [1] * num_turns_w1
    # V_C = 2.0
    # V_D = 1.0
    # voltages_winding_2 = [
    #     V_C - (V_C - V_D) * i / (num_turns_w2 - 1)
    #     for i in range(num_turns_w2)
    # ]

    # ---------------------------------------------------------------
    # # Simulation 10 (V_A, V_B, V_C, V_D = 1, 1, 2 , 2) --- (V_1, V_2, V_3, V_4 = 0, 1, 0, 1)
    # # Create a fixed voltage from C to D
    # voltages_winding_1 = [1] * num_turns_w1
    # voltages_winding_2 = [2] * num_turns_w2

    geo.electrostatic_simulation( voltage=[voltages_winding_1, voltages_winding_2], core_voltage=0, ground_outer_boundary=False,
                                 show_fem_simulation_results=show_visual_outputs, save_to_excel=False)
    # geo.get_total_charges()
    # geo.femm_reference_electrostatic(voltages=[voltages_winding_1, voltages_winding_2], ground_core=True, ground_outer_boundary=True, non_visualize=0, save_to_excel=False,
    #                                  compare_excel_files_to_femmt=False)

if __name__ == "__main__":
    basic_example_transformer_electrostatic(show_visual_outputs=True)
