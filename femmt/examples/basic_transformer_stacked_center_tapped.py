"""
Basic example to show how to build up the geometry for a stacked center tapped transformer.

After starting the program, the geometry dimensions are displayed. Verify this geometry, close the window, to continue the simulation.
After a short time, B-Field and winding losses simulation results are shown. Winding losses are shown as a colormap.
In the core, the magnitude B-Field in Tesla is shown. With the gmsh window, one can move the picture in the 3D way (not recommended).
If you close this window, the thermal simulation will be continued, if programmed. If true, the thermal heat distribution will be displayed.
To continue with the next simulation (or end the program), you need to close this window. All results are written to the result
folder .../femmt/examples/example_results/simulation_file_name/results/log_electro_magnetic.json. and .../results_thermal.json.
"""
import femmt as fmt
import os
import logging

# configure logging to show femmt terminal output
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

def basic_example_transformer_stacked_center_tapped(onelab_folder: str = None, show_visual_outputs: bool = True,
                                                    is_test: bool = False):
    """
    Show how to build up the geometry for a stacked center tapped transformer.

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

    working_directory = os.path.join(example_results_folder, os.path.splitext(os.path.basename(__file__))[0])
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.IntegratedTransformer,
                                working_directory=working_directory, is_gui=is_test)

    # This line is for automated pytest running on GitHub only. Please ignore this line!
    if onelab_folder is not None:
        geo.file_data.onelab_folder_path = onelab_folder

    core_dimensions = fmt.dtos.StackedCoreDimensions(core_inner_diameter=0.02, window_w=0.015, window_h_top=0.005,
                                                     window_h_bot=0.017)
    core = fmt.Core(core_type=fmt.CoreType.Stacked, core_dimensions=core_dimensions, mu_r_abs=3100, phi_mu_deg=12,
                    sigma=1.2,
                    permeability_datasource=fmt.MaterialDataSource.Custom,
                    permittivity_datasource=fmt.MaterialDataSource.Custom)
    geo.set_core(core)

    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Stacked, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.002, stacked_position=fmt.StackedPosition.Top)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.001, stacked_position=fmt.StackedPosition.Bot)
    geo.set_air_gaps(air_gaps)

    # set_center_tapped_windings() automatically places the condu
    insulation, coil_window, transformer_window = fmt.functions_topologies.set_center_tapped_windings(
        core=core,
        primary_turns=14,
        primary_radius=1.1e-3,
        primary_number_strands=50,
        primary_strand_radius=0.00011,
        secondary_parallel_turns=2,
        secondary_thickness_foil=1e-3,
        iso_top_core=0.001,
        iso_bot_core=0.001,
        iso_left_core=0.002,
        iso_right_core=0.001,
        iso_primary_to_primary=2e-4,
        iso_secondary_to_secondary=2e-4,
        iso_primary_to_secondary=4e-4,
        interleaving_type=fmt.CenterTappedInterleavingType.TypeC,
        interleaving_scheme=fmt.InterleavingSchemesFoilLitz.ter_3_4_sec_ter_4_3_sec,
        primary_coil_turns=3,
        primary_additional_bobbin=1e-3,
        winding_temperature=100,
        bobbin_coil_left=3e-3,
        center_foil_additional_bobbin=0e-3,
        wrap_para_type=fmt.WrapParaType.FixedThickness,
        foil_horizontal_placing_strategy=fmt.FoilHorizontalDistribution.VerticalUpward)

    insulation.add_turn_insulation([0.25e-5, 0.25e-5], add_turn_insulations=False)
    insulation.add_insulation_between_layers(add_kapton_material=False, thickness=0.0005)

    geo.set_insulation(insulation)
    geo.set_winding_windows([coil_window, transformer_window])

    geo.create_model(freq=200000, pre_visualize_geometry=show_visual_outputs)

    geo.single_simulation(freq=200000, current=[20, 120, 120], phi_deg=[0, 180, 180],
                          show_fem_simulation_results=show_visual_outputs)

    geo.get_inductances(I0=1, op_frequency=200000)
    # 7. prepare and start thermal simulation
    example_thermal_simulation(show_visual_outputs, flag_insulation=False)


if __name__ == "__main__":
    basic_example_transformer_stacked_center_tapped(show_visual_outputs=True)
