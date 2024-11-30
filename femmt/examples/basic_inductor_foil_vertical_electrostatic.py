"""
Basic example to show how to simulate an inductor with vertical foil winding.

After starting the program, the geometry dimensions are displayed. Verify this geometry, close the window, to continue the simulation.
After a short time, B-Field and winding losses simulation results are shown. Winding losses are shown as a colormap.
In the core, the magnitude B-Field in Tesla is shown. With the gmsh window, one can move the picture in the 3D way (not recommended).
If you close this window, the thermal simulation will be continued, if programmed. If true, the thermal heat distribution will be displayed.
To continue with the next simulation (or end the program), you need to close this window. All results are written to the result
folder .../femmt/examples/example_results/simulation_file_name/results/log_electro_magnetic.json. and .../results_thermal.json.
"""
import femmt as fmt
import os


def basic_example_inductor_foil_vertical_electrostatic(onelab_folder: str = None, show_visual_outputs: bool = True,
                                         is_test: bool = False):
    """
    Run the example code for the inductor with vertical foil winding.

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

    # Choose wrap para type
    wrap_para_type = fmt.WrapParaType.FixedThickness

    # Set is_gui = True so FEMMt won't ask for the onelab path if no config is found.
    geo = fmt.MagneticComponent(simulation_type=fmt.SimulationType.ElectroStatic, component_type=fmt.ComponentType.Inductor, working_directory=working_directory,
                                verbosity=fmt.Verbosity.ToConsole, is_gui=is_test)

    # This line is for automated pytest running on GitHub only. Please ignore this line!
    if onelab_folder is not None:
        geo.file_data.onelab_folder_path = onelab_folder

    core_db = fmt.core_database()["PQ 40/40"]
    core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                    window_w=core_db["window_w"],
                                                    window_h=core_db["window_h"],
                                                    core_h=core_db["core_h"])

    core = fmt.Core(core_type=fmt.CoreType.Single, core_dimensions=core_dimensions,
                    mu_r_abs=3100, phi_mu_deg=12,
                    sigma=0.6, permeability_datasource=fmt.MaterialDataSource.Custom,
                    permittivity_datasource=fmt.MaterialDataSource.Custom)
    geo.set_core(core)

    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Center, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005)
    geo.set_air_gaps(air_gaps)

    insulation = fmt.Insulation(flag_insulation=False)
    insulation.add_core_insulations(0.001, 0.001, 0.001, 0.001)
    insulation.add_winding_insulations([[0.0005]])
    geo.set_insulation(insulation)

    winding_window = fmt.WindingWindow(core, insulation)
    vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

    winding = fmt.Conductor(0, fmt.Conductivity.Copper, winding_material_temperature=25)
    winding.set_rectangular_conductor(thickness=1e-3)

    vww.set_winding(winding, 3, fmt.WindingScheme.FoilHorizontal, fmt.Align.ToEdges, wrap_para_type=wrap_para_type,
                    foil_horizontal_placing_strategy=fmt.FoilHorizontalDistribution.VerticalUpward)
    geo.set_winding_windows([winding_window])

    geo.create_model(freq=100000, pre_visualize_geometry=show_visual_outputs, save_png=False)

    # 8. run electrostatic simulation
    geo.electrostatic_simulation(voltage=[[0, 0, 5]], ground_core=False, ground_outer_boundary=False,
                                 show_fem_simulation_results=show_visual_outputs, save_to_excel=False)
    # Call the electrostatic FEMM simulation function
    # voltages = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140]
    # geo.femm_reference_electrostatic(voltages=[[5, 0]], ground_core=True, ground_outer_boundary=False,
    #                                  non_visualize=0, save_to_excel=True, compare_excel_files_to_femmt=True, mesh_size_conductor=0.0)


if __name__ == "__main__":
    basic_example_inductor_foil_vertical_electrostatic(show_visual_outputs=True)
