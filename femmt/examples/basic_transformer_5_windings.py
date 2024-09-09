"""
Basic example to show how to simulate a 5-winding transformer.

After starting the program, the geometry dimensions are displayed. Verify this geometry, close the window, to continue the simulation.
After a short time, B-Field and winding losses simulation results are shown. Winding losses are shown as a colormap.
In the core, the magnitude B-Field in Tesla is shown. With the gmsh window, one can move the picture in the 3D way (not recommended).
If you close this window, the thermal simulation will be continued, if programmed. If true, the thermal heat distribution will be displayed.
To continue with the next simulation (or end the program), you need to close this window. All results are written to the result
folder .../femmt/examples/example_results/simulation_file_name/results/log_electro_magnetic.json. and .../results_thermal.json.
"""
import femmt as fmt
import os


def basic_example_transformer_5_windings(onelab_folder: str = None, show_visual_outputs: bool = True,
                                         is_test: bool = False):
    """
    Run the example code for the 5-winding transformer.

    :param onelab_folder: onelab folder path
    :type onelab_folder: str
    :param show_visual_outputs: True to show visual outputs (simulation results)
    :type show_visual_outputs: bool
    :param is_test: True for pytest usage. Defaults to False.
    :type is_test: bool
    """
    def example_thermal_simulation(show_thermal_visual_outputs: bool = True, flag_insulation: bool = True):
        # Thermal simulation:
        # The losses calculated by the magnetics simulation can be used to calculate the heat distribution of the
        # given magnetic component. In order to use the thermal simulation, thermal conductivities for each material
        # can be entered as well as a boundary temperature which will be applied on the boundary of the
        # simulation (dirichlet boundary condition).

        # The case parameter sets the thermal conductivity for a case which will be set around the core.
        # This could model some case in which the transformer is placed in together with a set potting material.
        thermal_conductivity_dict = {
            "air": 0.0263,
            "case": {  # epoxy resign
                "top": 1.54,
                "top_right": 1.54,
                "right": 1.54,
                "bot_right": 1.54,
                "bot": 1.54
            },
            "core": 5,  # ferrite
            "winding": 400,  # copper
            "air_gaps": 180,  # aluminium nitride
            "insulation": 0.42 if flag_insulation else None  # polyethylene
        }

        # Here the case size can be determined
        case_gap_top = 0.002
        case_gap_right = 0.0025
        case_gap_bot = 0.002

        # Here the boundary temperatures can be set, currently it is set to 20°C (around 293°K).
        # This does not change the results of the simulation (at least when every boundary is set equally)
        # but will set the temperature offset.
        boundary_temperatures = {
            "value_boundary_top": 20,
            "value_boundary_top_right": 20,
            "value_boundary_right_top": 20,
            "value_boundary_right": 20,
            "value_boundary_right_bottom": 20,
            "value_boundary_bottom_right": 20,
            "value_boundary_bottom": 20
        }

        # Here the boundary sides can be turned on (1) or off (0)
        # By turning off the flag a neumann boundary will be applied at this point with heat flux = 0
        boundary_flags = {
            "flag_boundary_top": 0,
            "flag_boundary_top_right": 0,
            "flag_boundary_right_top": 1,
            "flag_boundary_right": 1,
            "flag_boundary_right_bottom": 1,
            "flag_boundary_bottom_right": 1,
            "flag_boundary_bottom": 1
        }

        # In order for the thermal simulation to work an electro_magnetic simulation has to run before.
        # The em-simulation will create a file containing the losses.
        # When the losses file is already created and contains the losses for the current model, it is enough to
        # run geo.create_model in order for the thermal simulation to work (geo.single_simulation is not needed).
        # Obviously when the model is modified and the losses can be out of date and therefore the
        # geo.single_simulation needs to run again.
        geo.thermal_simulation(thermal_conductivity_dict, boundary_temperatures, boundary_flags, case_gap_top,
                               case_gap_right, case_gap_bot, show_thermal_visual_outputs,
                               color_scheme=fmt.colors_ba_jonas, colors_geometry=fmt.colors_geometry_ba_jonas,
                               flag_insulation=flag_insulation)

    example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
    if not os.path.exists(example_results_folder):
        os.mkdir(example_results_folder)

    # Example for a transformer with multiple virtual winding windows.
    working_directory = os.path.join(example_results_folder, os.path.splitext(os.path.basename(__file__))[0])
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    # 1. chose simulation type
    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Transformer, working_directory=working_directory,
                                verbosity=fmt.Verbosity.ToConsole, is_gui=is_test)

    # This line is for automated pytest running on GitHub only. Please ignore this line!
    if onelab_folder is not None:
        geo.file_data.onelab_folder_path = onelab_folder

    # 2. set core parameters
    core_dimensions = fmt.dtos.SingleCoreDimensions(window_h=16.1e-3, window_w=(22.5 - 12) / 2 * 1e-3,
                                                    core_inner_diameter=12e-3, core_h=22e-3)
    core = fmt.Core(core_dimensions=core_dimensions, material=fmt.Material.N95, temperature=60, frequency=100000,
                    # permeability_datasource="manufacturer_datasheet",
                    permeability_datasource=fmt.MaterialDataSource.Measurement,
                    permeability_datatype=fmt.MeasurementDataType.ComplexPermeability,
                    permeability_measurement_setup=fmt.MeasurementSetup.LEA_LK,
                    permittivity_datasource=fmt.MaterialDataSource.Measurement,
                    permittivity_datatype=fmt.MeasurementDataType.ComplexPermittivity,
                    permittivity_measurement_setup=fmt.MeasurementSetup.LEA_LK)
    geo.set_core(core)

    # 3. set air gap parameters
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.00016, 50)
    geo.set_air_gaps(air_gaps)

    # 4. set insulation
    insulation = fmt.Insulation()
    insulation.add_core_insulations(0.0008, 0.0008, 0.001, 0.0001)
    iso_self = 0.0001
    iso_against = 0.0002
    insulation.add_winding_insulations(
        [[iso_self, iso_against, iso_against, iso_against, iso_against],
         [iso_against, iso_self, iso_against, iso_against, iso_against],
         [iso_against, iso_against, iso_self, iso_against, iso_against],
         [iso_against, iso_against, iso_self, iso_against, iso_against],
         [iso_against, iso_against, iso_against, iso_against, iso_self]])
    geo.set_insulation(insulation)

    # 5. create winding window and virtual winding windows (vww)
    winding_window = fmt.WindingWindow(core, insulation)
    cells = winding_window.NHorizontalAndVerticalSplit(horizontal_split_factors=[0.48, 0.75],
                                                       vertical_split_factors=[None, [0.5, 0.85], None])

    # 6. create windings and assign conductors
    winding1 = fmt.Conductor(0, fmt.Conductivity.Copper)
    winding1.set_litz_round_conductor(0.85e-3 / 2, 40, 0.1e-3 / 2, None, fmt.ConductorArrangement.Square)

    winding2 = fmt.Conductor(1, fmt.Conductivity.Copper)
    winding2.set_litz_round_conductor(1.0e-3 / 2, 60, 0.1e-3 / 2, None, fmt.ConductorArrangement.Square)

    winding3 = fmt.Conductor(2, fmt.Conductivity.Copper)
    winding3.set_litz_round_conductor(0.75e-3 / 2, 40, 0.1e-3 / 2, None, fmt.ConductorArrangement.Square)

    winding4 = fmt.Conductor(3, fmt.Conductivity.Copper)
    winding4.set_litz_round_conductor(0.95e-3 / 2, 40, 0.1e-3 / 2, None, fmt.ConductorArrangement.Square)

    winding5 = fmt.Conductor(4, fmt.Conductivity.Copper)
    winding5.set_litz_round_conductor(0.75e-3 / 2, 40, 0.1e-3 / 2, None, fmt.ConductorArrangement.Square)

    # 7. assign windings to virtual winding windows (cells)
    cells[0].set_winding(winding1, 22, fmt.WindingType.Single, fmt.Align.ToEdges, fmt.ConductorDistribution.VerticalUpward_HorizontalRightward)
    cells[1].set_winding(winding2, 6, fmt.WindingType.Single, fmt.Align.ToEdges, fmt.ConductorDistribution.VerticalUpward_HorizontalRightward)
    cells[2].set_winding(winding3, 6, fmt.WindingType.Single, fmt.Align.ToEdges, fmt.ConductorDistribution.VerticalUpward_HorizontalRightward)
    cells[3].set_winding(winding4, 1, fmt.WindingType.Single, fmt.Align.ToEdges, fmt.ConductorDistribution.VerticalUpward_HorizontalRightward)
    cells[4].set_winding(winding5, 2, fmt.WindingType.Single, fmt.Align.ToEdges, fmt.ConductorDistribution.VerticalUpward_HorizontalRightward)
    geo.set_winding_windows([winding_window])

    # 8. perform an FEM simulation
    geo.create_model(freq=100000, pre_visualize_geometry=show_visual_outputs)
    geo.single_simulation(freq=100000, current=[1.625, 6.0, 4.9, 0.5, 1],
                          phi_deg=[0, 180, 180, 90, 89], show_fem_simulation_results=show_visual_outputs)

    # 9. prepare and start thermal simulation
    example_thermal_simulation(show_thermal_visual_outputs=show_visual_outputs, flag_insulation=True)


if __name__ == "__main__":
    basic_example_transformer_5_windings(show_visual_outputs=True)
