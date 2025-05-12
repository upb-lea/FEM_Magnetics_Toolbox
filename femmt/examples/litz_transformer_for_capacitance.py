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

def basic_example_transformer_litz(onelab_folder: str = None, show_visual_outputs: bool = True, is_test: bool = False):
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
    geo = fmt.MagneticComponent(simulation_type=fmt.SimulationType.FreqDomain, component_type=fmt.ComponentType.Transformer,
                                working_directory=working_directory, verbosity=fmt.Verbosity.ToConsole, is_gui=is_test)

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
    bobbin_db = fmt.bobbin_database()["PQ 40/40"]
    bobbin_dimensions = fmt.dtos.BobbinDimensions(bobbin_inner_diameter=bobbin_db["bobbin_inner_diameter"],
                                                  bobbin_window_w=bobbin_db["bobbin_window_w"],
                                                  bobbin_window_h=bobbin_db["bobbin_window_h"],
                                                  bobbin_h=bobbin_db["bobbin_h"])
    core = fmt.Core(core_type=fmt.CoreType.Single,
                    core_dimensions=core_dimensions,
                    bobbin_dimensions=bobbin_dimensions,
                    detailed_core_model=False,
                    material=fmt.Material.N95, temperature=45, frequency=0,
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
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 1e-3, 50)
    geo.set_air_gaps(air_gaps)

    # 4. set insulation
    insulation = fmt.Insulation(flag_insulation=True)
    # insulation.add_core_insulations(2.2e-3, 2.2e-3, 1.25e-3, 1.25e-3)
    # # core_insulation (measured)
    # insulation.add_core_insulations(1.7e-3, 1.7e-3, 1.25e-3, 1.25e-3)
    # core_insulation (from datasheet)
    insulation.add_core_insulations(1.55e-3, 1.55e-3, 0.9e-3, 1.5e-4)
    # # # # 109-49 transformer
    # insulation.add_winding_insulations([[0.025e-3, 0.095e-3],
    #                                     [0.095e-3, 0.025e-3]])
    # insulation.add_winding_insulations([[0.025e-3, 0.095e-3],
    #                                     [0.095e-3, 0.025e-3]])
    # insulation.add_winding_insulations([[1.4946e-4, 0.095e-3],
    #                                     [0.095e-3, 1.4946e-4]])
    insulation.add_winding_insulations([[0.29e-3, 0.095e-3],
                                        [0.095e-3, 0.29e-3]])
    insulation.add_conductor_air_conductor_insulation([[2.5e-4, 3e-4, 2.281e-4, 3.064e-4],
                                                       [2.5e-4, 3e-4]])
    # othman
    insulation.add_kapton_insulation(add_kapton=True, thickness=0.025e-3)
    # till
    #insulation.add_kapton_insulation(add_kapton=True, thickness=0.9e-3)
    # 59-56 transformer
    # insulation.add_winding_insulations([[1.33125e-4, 0.095e-3],
    #                                     [0.095e-3, 1.88667e-4]])

    geo.set_insulation(insulation)

    # 5. create winding window and virtual winding windows (vww)
    winding_window = fmt.WindingWindow(core, insulation)
    # bot, top = winding_window.split_window(fmt.WindingWindowSplit.HorizontalSplit, split_distance=0.001)
    # 109-49
    #othman
    cells = winding_window.NHorizontalAndVerticalSplit(horizontal_split_factors=[0.43],
                                                       vertical_split_factors=None)
    # till
    # cells = winding_window.NHorizontalAndVerticalSplit(horizontal_split_factors=[0.29],
    #                                                    vertical_split_factors=None)
    # cells = winding_window.NHorizontalAndVerticalSplit(horizontal_split_factors=[0.],
    #                                                    vertical_split_factors=None)

    # 6. create conductors and set parameters
    winding1 = fmt.Conductor(0, fmt.Conductivity.Copper)
    # winding1.set_solid_round_conductor(1.1506e-3, fmt.ConductorArrangement.Square)
    # winding1.set_solid_round_conductor(0.0008728074876212745, fmt.ConductorArrangement.Square)
    #winding1.set_litz_round_conductor(None, 405, 35.5e-6, 0.67, fmt.ConductorArrangement.Square)
    winding1.set_solid_round_conductor(0.71e-3, fmt.ConductorArrangement.Square)

    # winding1 = fmt.Conductor(0, fmt.Conductivity.Copper)
    # winding1.set_litz_round_conductor(0.0011, 50, 0.00011, None, fmt.ConductorArrangement.Square)

    # winding2 = fmt.Conductor(1, fmt.Conductivity.Copper)
    # winding2.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

    winding2 = fmt.Conductor(1, fmt.Conductivity.Copper)
    # winding2.set_solid_round_conductor(1.1506e-3, fmt.ConductorArrangement.Square)
    # winding2.set_solid_round_conductor(0.0012696221085454768, fmt.ConductorArrangement.Square)
    winding2.set_solid_round_conductor(0.71e-3, fmt.ConductorArrangement.Square)
    winding2.parallel = False
    #winding2.set_litz_round_conductor(None, 1200, 30e-6, 0.67, fmt.ConductorArrangement.Square)

    # 7. add conductor to vww and add winding window to MagneticComponent
    # top.set_winding(winding2, 15, None, fmt.Align.ToEdges, fmt.ConductorDistribution.VerticalUpward_HorizontalRightward, zigzag=True)
    cells[1].set_winding(winding2, 21, None, fmt.Align.ToEdges, fmt.ConductorDistribution.VerticalUpward_HorizontalRightward, zigzag=True)
    # bot.set_winding(winding2, 29, None, fmt.Align.ToEdges, fmt.ConductorDistribution.VerticalUpward_HorizontalRightward, zigzag=True)
    cells[0].set_winding(winding1, 21, None, fmt.Align.ToEdges, fmt.ConductorDistribution.VerticalUpward_HorizontalRightward, zigzag=True)
    # top.set_winding(winding1, 109, None, fmt.Align.ToEdges, fmt.ConductorDistribution.VerticalUpward_HorizontalRightward, zigzag=True)
    geo.set_winding_windows([winding_window])



    # 8. start simulation with given frequency, currents and phases
    # [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    # [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150]
    # [30, 35, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]
    # [20, 40, 60, 80, 100]
    geo.create_model(freq=0, pre_visualize_geometry=True)

    # geo.get_total_charges()
    # geo.femm_reference_electrostatic(voltages=[voltages_winding_1, voltages_winding_2], ground_core=True, ground_outer_boundary=True, non_visualize=0, save_to_excel=False,
    #                                  compare_excel_files_to_femmt=False)
    # geo.get_inductances(I0=2, op_frequency=20000, skin_mesh_factor=0.5)
    # geo.single_simulation(freq=200000, current=[10, 14], phi_deg=[0, 180],
    #                       plot_interpolation=False, show_fem_simulation_results=show_visual_outputs)
    # geo.get_transformer_capacitance(flag_cd=True)
    meas = (
        6.05915e-11, 3.581451055984E-11, 6.35685E-11,
        2.66E-11, 4.79E-11, 5.97E-11,
        4.33E-11, None, None,
        None
    )

    geo.get_transformer_capacitance(c_meas_short= None, c_meas_open=4.268050496E-11
, measured_values = meas, show_plot_comparison=True)

if __name__ == "__main__":
    basic_example_transformer_litz(show_visual_outputs=True)
