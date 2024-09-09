"""
Example how to use the split winding method. Run this file, to see the different winding orders.

Just shows different split methods. No thermal simulation.

After starting the program, the geometry dimensions are displayed. Verify this geometry, close the window, to continue the simulation.
After a short time, B-Field and winding losses simulation results are shown. Winding losses are shown as a colormap.
In the core, the magnitude B-Field in Tesla is shown. With the gmsh window, one can move the picture in the 3D way (not recommended).
If you close this window, the thermal simulation will be continued, if programmed. If true, the thermal heat distribution will be displayed.
To continue with the next simulation (or end the program), you need to close this window. All results are written to the result
folder .../femmt/examples/example_results/simulation_file_name/results/log_electro_magnetic.json. and .../results_thermal.json.
"""
import femmt as fmt
import os

def run_transformer_vvw_split_examples(num_windings: int, onelab_folder: str = None, show_visual_outputs: bool = True, is_test: bool = False):
    """
    Run the example code for the transformer.

    :param num_windings: Number of windings to simulate
    :type num_windings: int
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
    working_directory_group_folder = os.path.join(example_results_folder, os.path.splitext(os.path.basename(__file__))[0])
    if not os.path.exists(working_directory_group_folder):
        os.mkdir(working_directory_group_folder)

    def setup_simulation(working_directory, horizontal_split_factors, vertical_split_factors):
        geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Transformer, working_directory=working_directory,
                                    verbosity=fmt.Verbosity.Silent, is_gui=is_test)

        # This line is for automated pytest running on GitHub only. Please ignore this line!
        if onelab_folder is not None:
            geo.file_data.onelab_folder_path = onelab_folder

        core_dimensions = fmt.dtos.SingleCoreDimensions(window_h=16.1e-3, window_w=(22.5 - 12) / 2 * 1e-3,
                                                        core_inner_diameter=12e-3, core_h=22e-3)
        core = fmt.Core(core_dimensions=core_dimensions, material=fmt.Material.N95, temperature=60, frequency=100000,
                        permeability_datasource=fmt.MaterialDataSource.Measurement,
                        permeability_datatype=fmt.MeasurementDataType.ComplexPermeability,
                        permeability_measurement_setup=fmt.MeasurementSetup.LEA_LK,
                        permittivity_datasource=fmt.MaterialDataSource.Measurement,
                        permittivity_datatype=fmt.MeasurementDataType.ComplexPermittivity,
                        permittivity_measurement_setup=fmt.MeasurementSetup.LEA_LK)
        geo.set_core(core)

        air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.00016, 50)
        geo.set_air_gaps(air_gaps)

        insulation = fmt.Insulation()
        insulation.add_core_insulations(0.0008, 0.0008, 0.0001, 0.00001)
        iso_self = 0.0001
        iso_against = 0.0002
        insulation.add_winding_insulations(
            [[iso_self, iso_against, iso_against, iso_against, iso_against, iso_against, iso_against],
             [iso_against, iso_self, iso_against, iso_against, iso_against, iso_against, iso_against],
             [iso_against, iso_against, iso_self, iso_against, iso_against, iso_against, iso_against],
             [iso_against, iso_against, iso_against, iso_self, iso_against, iso_against, iso_against],
             [iso_against, iso_against, iso_against, iso_against, iso_self, iso_against, iso_against],
             [iso_against, iso_against, iso_against, iso_against, iso_against, iso_self, iso_against],
             [iso_against, iso_against, iso_against, iso_against, iso_against, iso_against, iso_self]])
        geo.set_insulation(insulation)

        winding_window = fmt.WindingWindow(core, insulation)

        cells = winding_window.flexible_split(
            horizontal_split_factors=horizontal_split_factors,
            vertical_split_factors=vertical_split_factors
        )

        windings = []
        for i in range(num_windings):
            winding = fmt.Conductor(i, fmt.Conductivity.Copper)
            winding.set_litz_round_conductor(0.85e-3 / 2, 40, 0.1e-3 / 2, None, fmt.ConductorArrangement.Square)
            windings.append(winding)

        for i in range(num_windings):
            cells[i].set_winding(windings[i], 7 - i, fmt.WindingType.Single, fmt.Align.ToEdges,
                                 fmt.ConductorDistribution.VerticalUpward_HorizontalRightward)

        geo.set_winding_windows([winding_window])

        geo.create_model(freq=100000, pre_visualize_geometry=show_visual_outputs)

        return geo

    if num_windings == 2:
        # Run with vertical split
        working_directory = os.path.join(working_directory_group_folder, "2-windings-vertical-split-only")
        setup_simulation(working_directory, horizontal_split_factors=[], vertical_split_factors=[[0.5]])

        # Run with horizontal split
        working_directory = os.path.join(working_directory_group_folder, "2-windings-horizontal-split-only")
        setup_simulation(working_directory, horizontal_split_factors=[0.5], vertical_split_factors=[])

    elif num_windings == 3:
        working_directory = os.path.join(working_directory_group_folder, "3-windings-vertical-split-only")
        setup_simulation(working_directory, horizontal_split_factors=[], vertical_split_factors=[[0.33, 0.66]])

        working_directory = os.path.join(working_directory_group_folder, "3-windings-horizontal-split-only")
        setup_simulation(working_directory, horizontal_split_factors=[0.33, 0.66], vertical_split_factors=[])

    elif num_windings == 5:
        working_directory = os.path.join(working_directory_group_folder, "5-windings")
        setup_simulation(working_directory, horizontal_split_factors=[0.48, 0.75], vertical_split_factors=[[0.5], [0.5], None])

    elif num_windings == 6:
        working_directory = os.path.join(working_directory_group_folder, "6-windings")
        setup_simulation(working_directory, horizontal_split_factors=[0.48, 0.75], vertical_split_factors=[[0.5], [0.5], [0.5]])

    else:
        raise ValueError("Unsupported number of windings")


if __name__ == "__main__":
    # Run simulations for different numbers of windings
    for num_windings in [2, 3, 5, 6]:
        print(f"Running simulation for {num_windings} windings")
        run_transformer_vvw_split_examples(num_windings, show_visual_outputs=True)
