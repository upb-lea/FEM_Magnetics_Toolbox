"""Basic example to show how to simulate an inductor with vertical foil winding."""
import femmt as fmt
import os


def basic_example_inductor_foil_vertical(onelab_folder: str = None, show_visual_outputs: bool = True,
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
    working_directory = os.path.join(example_results_folder, "basic_inductor_foil_vertical")
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    # Choose wrap para type
    wrap_para_type = fmt.WrapParaType.FixedThickness

    # Set is_gui = True so FEMMt won't ask for the onelab path if no config is found.
    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=working_directory,
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

    insulation = fmt.Insulation()
    insulation.add_core_insulations(0.001, 0.001, 0.002, 0.001)
    insulation.add_winding_insulations([[0.0005]])
    geo.set_insulation(insulation)

    winding_window = fmt.WindingWindow(core, insulation)
    vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

    winding = fmt.Conductor(0, fmt.Conductivity.Copper, winding_material_temperature=25)
    winding.set_rectangular_conductor(thickness=1e-3)

    vww.set_winding(winding, 5, fmt.WindingScheme.FoilVertical, fmt.Align.ToEdges, wrap_para_type=wrap_para_type,
                    foil_vertical_placing_strategy=fmt.FoilVerticalDistribution.HorizontalRightward)
    geo.set_winding_windows([winding_window])

    geo.create_model(freq=100000, pre_visualize_geometry=show_visual_outputs, save_png=False)

    geo.single_simulation(freq=100000, current=[3], show_fem_simulation_results=show_visual_outputs)

    example_thermal_simulation(show_visual_outputs, flag_insulation=True)


if __name__ == "__main__":
    basic_example_inductor_foil_vertical(show_visual_outputs=True)
