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



def basic_example_transformer_electrostatic_ansys(onelab_folder: str = None, show_visual_outputs: bool = True, is_test: bool = False):
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
                                working_directory=working_directory, verbosity=fmt.Verbosity.ToConsole, is_gui=is_test)

    # This line is for automated pytest running on GitHub only. Please ignore this line!
    if onelab_folder is not None:
        geo.file_data.onelab_folder_path = onelab_folder

    # 2. set core parameters
    core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=0.02, window_w=0.01, window_h=0.03,
                                                    core_h=0.02)
    core = fmt.Core(core_dimensions=core_dimensions, mu_r_abs=3100, phi_mu_deg=12, sigma=1.2,
                    permeability_datasource=fmt.MaterialDataSource.Custom,
                    permittivity_datasource=fmt.MaterialDataSource.Custom,
                    detailed_core_model=False)
    geo.set_core(core)
    # 3. set air gap parameters
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.001, 50)
    geo.set_air_gaps(air_gaps)

    # 4. set insulation
    insulation = fmt.Insulation(flag_insulation=True)
    insulation.add_core_insulations(0.001, 0.001, 0.001, 0.001)
    # This is actually now the insulation of the winding itself (insulation around the conductor)
    insulation.add_winding_insulations([[0.0002, 0.001],
                                        [0.001, 0.0002]])
    # An Air-Insulation in every layer for every winding
    insulation.add_air_conductor_insulations([[0.0002, 0.001, 0.002],
                                              [0.001, 0.0002, 0.001]])
    insulation.add_kapton_insulation(add_kapton=True, thickness=0.05e-3)
    geo.set_insulation(insulation)

    # 5. create winding window and virtual winding windows (vww)
    winding_window = fmt.WindingWindow(core, insulation)
    bot, top = winding_window.split_window(fmt.WindingWindowSplit.HorizontalSplit, split_distance=0)

    # 6. create conductors and set parameters
    winding1 = fmt.Conductor(0, fmt.Conductivity.Copper)
    # winding1.set_solid_round_conductor(1.1506e-3, fmt.ConductorArrangement.Square)
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
    bot.set_winding(winding2, 7, None, fmt.Align.ToEdges, fmt.ConductorDistribution.VerticalDownward_HorizontalLeftward, zigzag=False)
    # bot.set_winding(winding1, 29, None, fmt.Align.ToEdges, fmt.ConductorDistribution.VerticalUpward_HorizontalRightward, zigzag=True)
    top.set_winding(winding1, 8, None, fmt.Align.ToEdges, fmt.ConductorDistribution.HorizontalRightward_VerticalUpward, zigzag=False)
    geo.set_winding_windows([winding_window])

    # 8. start simulation with given frequency, currents and phases
    geo.create_model(freq=0, pre_visualize_geometry=show_visual_outputs)
    geo.electrostatic_simulation( voltage=[[1, 2, 3, 4 , 5, 6, 7, 8], [3, 6, 9, 12, 15, 14, 12]], core_voltage=0, ground_outer_boundary=False,
                                 show_fem_simulation_results=show_visual_outputs, save_to_excel=False)
    geo.femm_reference_electrostatic(voltages=[[1, 2, 3, 4 , 5, 6, 7, 8], [3, 6, 9, 12, 15, 14, 12]], ground_core=True, ground_outer_boundary=False, non_visualize=0, save_to_excel=False,
                                     compare_excel_files_to_femmt=False)

if __name__ == "__main__":
    basic_example_transformer_electrostatic_ansys(show_visual_outputs=True)
