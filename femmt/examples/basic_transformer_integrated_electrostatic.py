"""
Basic example to show how to simulate an integrated transformer.

After starting the program, the geometry dimensions are displayed. Verify this geometry, close the window, to continue the simulation.
After a short time, B-Field and winding losses simulation results are shown. Winding losses are shown as a colormap.
In the core, the magnitude B-Field in Tesla is shown. With the gmsh window, one can move the picture in the 3D way (not recommended).
If you close this window, the thermal simulation will be continued, if programmed. If true, the thermal heat distribution will be displayed.
To continue with the next simulation (or end the program), you need to close this window. All results are written to the result
folder .../femmt/examples/example_results/simulation_file_name/results/log_electro_magnetic.json. and .../results_thermal.json.
"""
import os
import femmt as fmt


def basic_example_transformer_integrated_electrostatic(onelab_folder: str = None, show_visual_outputs: bool = True,
                                         is_test: bool = False):
    """
    Run the example code for the integrated transformer.

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
    core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=0.02, window_w=0.011, window_h=0.03,
                                                    core_h=0.08)
    core = fmt.Core(core_dimensions=core_dimensions, mu_r_abs=3100, phi_mu_deg=12, sigma=0.6,
                    permeability_datasource=fmt.MaterialDataSource.Custom,
                    permittivity_datasource=fmt.MaterialDataSource.Custom)
    geo.set_core(core)

    # 2.1 set stray path parameters
    stray_path = fmt.StrayPath(start_index=0, length=geo.core.core_inner_diameter / 2 + geo.core.window_w - 0.001)
    geo.set_stray_path(stray_path)
    print(stray_path)

    # 3. set air gap parameters
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.003, 30)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.003, 60)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.001, 80)
    geo.set_air_gaps(air_gaps)

    # 4. set insulations
    insulation = fmt.Insulation(flag_insulation=False)
    insulation.add_core_insulations(0.001, 0.001, 0.002, 0.001)
    insulation.add_winding_insulations([[0.0002, 0.001],
                                        [0.001, 0.0002]])
    geo.set_insulation(insulation)

    # 5. create winding window and virtual winding windows (vww)
    # For an integrated transformer it is not necessary to set horizontal and vertical split factors
    # since this is determined by the stray_path
    winding_window = fmt.WindingWindow(core, insulation, stray_path, air_gaps)
    top, bot = winding_window.split_window(fmt.WindingWindowSplit.HorizontalSplit)

    # 6. set conductor parameters
    winding1 = fmt.Conductor(0, fmt.Conductivity.Copper)
    # winding1.set_litz_round_conductor(None, 100, 70e-6, 0.5, fmt.ConductorArrangement.Square)
    winding1.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

    winding2 = fmt.Conductor(1, fmt.Conductivity.Copper)
    # winding2.set_litz_round_conductor(None, 100, 70e-6, 0.5, fmt.ConductorArrangement.Square)
    winding2.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

    # 7. add conductor to vww and add winding window to MagneticComponent
    top.set_interleaved_winding(winding1, 3, winding2, 4, fmt.InterleavedWindingScheme.HorizontalAlternating)
    bot.set_interleaved_winding(winding1, 1, winding2, 2, fmt.InterleavedWindingScheme.HorizontalAlternating)
    geo.set_winding_windows([winding_window])

    # 8. start simulation with given frequency, currents and phases
    geo.create_model(freq=250000, pre_visualize_geometry=show_visual_outputs)
    geo.electrostatic_simulation(voltage=[[10, 20, 30, 40], [10, 20, 30, 40, 50, 60]],
                                 ground_core=False, ground_outer_boundary=False, show_fem_simulation_results=show_visual_outputs, save_to_excel=True)

if __name__ == "__main__":
    basic_example_transformer_integrated_electrostatic(show_visual_outputs=True)
