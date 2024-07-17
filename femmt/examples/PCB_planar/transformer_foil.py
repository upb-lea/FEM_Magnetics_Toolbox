import femmt as fmt
import os

from datetime import datetime
Date = datetime.now().strftime("%Y%m%d-%H%M%S")
def basic_example_transformer_foil(onelab_folder: str = None, show_visual_outputs: bool = True, is_test: bool = False):
    """Run the example code for the transformer."""

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
    working_directory = os.path.join(example_results_folder, "transformer" + Date)
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    # Choose wrap para type
    wrap_para_type = fmt.WrapParaType.FixedThickness

    # 1. chose simulation type
    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Transformer, working_directory=working_directory,
                                verbosity=fmt.Verbosity.ToConsole, is_gui=is_test)

    # This line is for automated pytest running on GitHub only. Please ignore this line!
    if onelab_folder is not None:
        geo.file_data.onelab_folder_path = onelab_folder

    # 2. set core parameters
    core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=5e-3,
                                                    window_w=6e-3,
                                                    window_h=2e-3,
                                                    core_h=8e-3)


    core = fmt.Core(core_type=fmt.CoreType.Single, core_dimensions=core_dimensions,
                    mu_r_abs=3100, phi_mu_deg=12,
                    sigma=0.6, permeability_datasource=fmt.MaterialDataSource.Custom,
                    permittivity_datasource=fmt.MaterialDataSource.Custom)
    geo.set_core(core)

    # 3. set air gap parameters
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
    geo.set_air_gaps(air_gaps)

    # 4. set insulation
    prepeg_thickness = 0.00023
    wind_insulation = prepeg_thickness / 2
    insulation = fmt.Insulation(flag_insulation=True)
    insulation.add_core_insulations(0.0005, 0.0005, 0.0005, 0.0005)
    insulation.add_winding_insulations([[wind_insulation, wind_insulation],
                                        [wind_insulation, wind_insulation]])
    geo.set_insulation(insulation)
    # 4 virtuelle Wicklungs fenster
    # 5. create winding window and virtual winding windows (vww)
    winding_window = fmt.WindingWindow(core, insulation)
    #vww_bot, vww_top = winding_window.split_window(fmt.WindingWindowSplit.HorizontalSplit, split_distance=0.001)
    vww_top_left, vww_top_right, vww_bot_left, vww_bot_right = winding_window.split_window(fmt.WindingWindowSplit.HorizontalAndVerticalSplit, split_distance=0.0001)
    # 6. create conductors and set parameters
    winding1 = fmt.Conductor(0, fmt.Conductivity.Copper)
    winding1.set_rectangular_conductor(thickness=35e-6)

    winding2 = fmt.Conductor(1, fmt.Conductivity.Copper)
    winding2.set_rectangular_conductor(thickness=35e-6)
    winding2.parallel = False

    # 7. add conductor to vww and add winding window to MagneticComponent

    #PlacingStrategy sorgt für Fehler: "Zu viele Leiter",
    vww_top_left.set_winding(winding1, 1, fmt.WindingScheme.FoilHorizontal, wrap_para_type=wrap_para_type, foil_horizontal_placing_strategy=fmt.FoilHorizontalDistribution.VerticalDownward)
    vww_top_right.set_winding(winding1, 1, fmt.WindingScheme.FoilHorizontal, wrap_para_type=wrap_para_type, foil_horizontal_placing_strategy=fmt.FoilHorizontalDistribution.VerticalDownward)
    vww_bot_left.set_winding(winding2, 1, fmt.WindingScheme.FoilHorizontal, wrap_para_type=wrap_para_type, foil_horizontal_placing_strategy=fmt.FoilHorizontalDistribution.VerticalDownward)
    vww_bot_right.set_winding(winding2, 1, fmt.WindingScheme.FoilHorizontal, wrap_para_type=wrap_para_type, foil_horizontal_placing_strategy=fmt.FoilHorizontalDistribution.VerticalDownward)

    # Oben Wicklung 1, unten Wicklung 2
    #vww_top_left.set_winding(winding1, 3, fmt.WindingScheme.FoilHorizontal,wrap_para_type=wrap_para_type,)
    #vww_top_right.set_winding(winding1, 3, fmt.WindingScheme.FoilHorizontal,wrap_para_type=wrap_para_type)
    #vww_bot_left.set_winding(winding2, 3, fmt.WindingScheme.FoilHorizontal, wrap_para_type=wrap_para_type)
    #vww_bot_right.set_winding(winding2, 3, fmt.WindingScheme.FoilHorizontal, wrap_para_type=wrap_para_type)

    #links Wicklung 1, rechts wicklung 2
    #vww_top_left.set_winding(winding1, 3, fmt.WindingScheme.FoilHorizontal,wrap_para_type=wrap_para_type,)
    #vww_top_right.set_winding(winding2, 3, fmt.WindingScheme.FoilHorizontal,wrap_para_type=wrap_para_type)
    #vww_bot_left.set_winding(winding1, 3, fmt.WindingScheme.FoilHorizontal, wrap_para_type=wrap_para_type)
    #vww_bot_right.set_winding(winding2, 3, fmt.WindingScheme.FoilHorizontal, wrap_para_type=wrap_para_type)


    #Wicklungen Diagonal
    #vww_top_left.set_winding(winding1, 3, fmt.WindingScheme.FoilHorizontal,wrap_para_type=wrap_para_type,)
    #vww_top_right.set_winding(winding2, 3, fmt.WindingScheme.FoilHorizontal,wrap_para_type=wrap_para_type)
    #vww_bot_left.set_winding(winding2, 3, fmt.WindingScheme.FoilHorizontal, wrap_para_type=wrap_para_type)
    #vww_bot_right.set_winding(winding1, 3, fmt.WindingScheme.FoilHorizontal, wrap_para_type=wrap_para_type)
    geo.set_winding_windows([winding_window])

    # NxM Virtuelle Wicklungsfenster
    #cells = winding_window.NHorizontalAndVerticalSplit(horizontal_split_factors=[0.33, 0.67],vertical_split_factors=[[1/6, 2/6, 3/6, 4/6, 5/6],[1/6, 2/6, 3/6, 4/6, 5/6], [1/6, 2/6, 3/6, 4/6, 5/6]])


    # 8. start simulation with given frequency, currents and phases
    geo.create_model(freq=1000000, pre_visualize_geometry=show_visual_outputs)
    geo.single_simulation(freq=1000000, current=[7, 6], phi_deg=[0, 180],
                          show_fem_simulation_results=show_visual_outputs)

    example_thermal_simulation(show_visual_outputs, flag_insulation=True)


if __name__ == "__main__":
    basic_example_transformer_foil(show_visual_outputs=True)