"""
Basic example to show how to simulate a stacked transformer.

After starting the program, the geometry dimensions are displayed. Verify this geometry, close the window, to continue the simulation.
After a short time, B-Field and winding losses simulation results are shown. Winding losses are shown as a colormap.
In the core, the magnitude B-Field in Tesla is shown. With the gmsh window, one can move the picture in the 3D way (not recommended).
If you close this window, the thermal simulation will be continued, if programmed. If true, the thermal heat distribution will be displayed.
To continue with the next simulation (or end the program), you need to close this window. All results are written to the result
folder .../femmt/examples/example_results/simulation_file_name/results/log_electro_magnetic.json. and .../results_thermal.json.
"""
import femmt as fmt
import os


def basic_example_transformer_stacked_electrostatic(onelab_folder: str = None, show_visual_outputs: bool = True,
                                      is_test: bool = False):
    """
    Run the example code for the stacked transformer.

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

    working_directory = os.path.join(example_results_folder, os.path.splitext(os.path.basename(__file__))[0])
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    # 1. chose simulation type
    geo = fmt.MagneticComponent(simulation_type=fmt.SimulationType.ElectroStatic, component_type=fmt.ComponentType.IntegratedTransformer,
                                working_directory=working_directory, verbosity=fmt.Verbosity.ToConsole, is_gui=is_test)

    # This line is for automated pytest running on GitHub only. Please ignore this line!
    if onelab_folder is not None:
        geo.file_data.onelab_folder_path = onelab_folder

    # 2. set core parameters
    core_dimensions = fmt.dtos.StackedCoreDimensions(core_inner_diameter=0.02, window_w=0.02, window_h_top=0.01,
                                                     window_h_bot=0.03)
    core = fmt.Core(core_type=fmt.CoreType.Stacked, core_dimensions=core_dimensions, mu_r_abs=3100, phi_mu_deg=12,
                    sigma=0.6,
                    permeability_datasource=fmt.MaterialDataSource.Custom,
                    permittivity_datasource=fmt.MaterialDataSource.Custom)
    geo.set_core(core)

    # 3. set air gap parameters
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Stacked, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.002, stacked_position=fmt.StackedPosition.Top)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.001, stacked_position=fmt.StackedPosition.Bot)
    geo.set_air_gaps(air_gaps)

    # 4. set insulations
    insulation = fmt.Insulation(flag_insulation=False)
    insulation.add_core_insulations(0.001, 0.001, 0.001, 0.001)  # [bot, top, left, right]
    insulation.add_winding_insulations([[0.0002, 0.001],
                                        [0.001, 0.0002]])
    geo.set_insulation(insulation)

    winding_window_top, winding_window_bot = fmt.create_stacked_winding_windows(core, insulation)

    vww_top = winding_window_top.split_window(fmt.WindingWindowSplit.NoSplit)
    vww_bot = winding_window_bot.split_window(fmt.WindingWindowSplit.NoSplit)

    # 6. set conductor parameters
    winding1 = fmt.Conductor(0, fmt.Conductivity.Copper)
    winding1.set_solid_round_conductor(1e-3, fmt.ConductorArrangement.Square)
    # winding1.set_litz_round_conductor(None, 100, 70e-6, 0.5, fmt.ConductorArrangement.Square)

    winding2 = fmt.Conductor(1, fmt.Conductivity.Copper)
    winding2.set_solid_round_conductor(1e-3, fmt.ConductorArrangement.Square)
    # winding2.set_litz_round_conductor(None, 120, 70e-6, 0.5, fmt.ConductorArrangement.Square)

    # 7. add conductor to vww and add winding window to MagneticComponent
    vww_top.set_interleaved_winding(winding1, 16, winding2, 0, fmt.InterleavedWindingScheme.HorizontalAlternating)
    vww_bot.set_interleaved_winding(winding1, 9, winding2, 20, fmt.InterleavedWindingScheme.HorizontalAlternating)

    geo.set_winding_windows([winding_window_top, winding_window_bot])

    # 8. start simulation with given frequency, currents and phases
    geo.create_model(freq=0, pre_visualize_geometry=show_visual_outputs)
    geo.electrostatic_simulation(voltage=[[10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250],
                                          [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200]],
                                 ground_core=False, ground_outer_boundary=True, show_fem_simulation_results=show_visual_outputs, save_to_excel=True)

if __name__ == "__main__":
    basic_example_transformer_stacked_electrostatic(show_visual_outputs=True)
