import femmt as fmt
import os

def basic_example_transformer_center_tapped(onelab_folder: str = None, show_visual_outputs: bool = True, is_test: bool = False):
    def example_thermal_simulation(show_visual_outputs: bool = True):
        # Thermal simulation:
        # The losses calculated by the magnetics simulation can be used to calculate the heat distribution of the given magnetic component
        # In order to use the thermal simulation, thermal conductivities for each material can be entered as well as a boundary temperature
        # which will be applied on the boundary of the simulation (dirichlet boundary condition).

        # The case parameter sets the thermal conductivity for a case which will be set around the core.
        # This could model some case in which the transformer is placed in together with a set potting material.
        thermal_conductivity_dict = {
            "air": 0.0263,
            "case": { # epoxy resign
                "top": 1.54,
                "top_right": 1.54,
                "right": 1.54,
                "bot_right": 1.54,
                "bot": 1.54
            },
            "core": 5, # ferrite
            "winding": 400, # copper
            "air_gaps": 180, # aluminiumnitride
            "insulation": 0.42 # polyethylen
        }

        # Here the case size can be determined
        case_gap_top = 0.002
        case_gap_right = 0.0025
        case_gap_bot = 0.002

        # Here the boundary temperatures can be set, currently it is set to 20°C (around 293°K).
        # This does not change the results of the simulation (at least when every boundary is set equally) but will set the temperature offset.
        boundary_temperatures = {
            "value_boundary_top": 20,
            "value_boundary_top_right": 20,
            "value_boundary_right_top": 20,
            "value_boundary_right": 20,
            "value_boundary_right_bottom": 20,
            "value_boundary_bottom_right": 20,
            "value_boundary_bottom": 20
        }

        # In order to compare the femmt thermal simulation with a femm heat flow simulation the same boundary temperature should be applied.
        # Currently only one temperature can be applied which will be set on every boundary site.
        femm_boundary_temperature = 20

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
        # When the losses file is already created and contains the losses for the current model, it is enough to run geo.create_model in
        # order for the thermal simulation to work (geo.single_simulation is not needed).
        # Obviously when the model is modified and the losses can be out of date and therefore the geo.single_simulation needs to run again.
        geo.thermal_simulation(thermal_conductivity_dict, boundary_temperatures, boundary_flags, case_gap_top,
                               case_gap_right, case_gap_bot, show_visual_outputs, color_scheme=fmt.colors_ba_jonas, colors_geometry=fmt.colors_geometry_ba_jonas)

    example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
    if not os.path.exists(example_results_folder):
        os.mkdir(example_results_folder)



    working_directory = os.path.join(example_results_folder, "center-tapped-transformer")
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Transformer, working_directory=working_directory,
                                verbosity=fmt.Verbosity.ToConsole, is_gui=is_test)

    # This line is for automated pytest running on github only. Please ignore this line!
    if onelab_folder is not None: geo.file_data.onelab_folder_path = onelab_folder

    core_dimensions = fmt.dtos.SingleCoreDimensions(window_h=0.025, window_w=0.02, core_inner_diameter=0.015)
    core = fmt.Core(core_dimensions=core_dimensions, mu_r_abs=3100, phi_mu_deg=12, sigma=1.2,
                    permeability_datasource=fmt.MaterialDataSource.Custom, permittivity_datasource=fmt.MaterialDataSource.Custom)
    geo.set_core(core)

    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
    geo.set_air_gaps(air_gaps)

    # set_center_tapped_windings() automatically places the conductors
    insulation, winding_window = fmt.functions_topologies.set_center_tapped_windings(core=core,
                                                                                     primary_turns=12, primary_radius=1.1e-3, primary_number_strands=50, primary_strand_radius=0.00011,
                                                                                     secondary_parallel_turns=3, secondary_thickness_foil=1e-3,
                                                                                     iso_top_core=0.001, iso_bot_core=0.001, iso_left_core=0.002, iso_right_core=0.001,
                                                                                     iso_primary_to_primary=1e-4, iso_secondary_to_secondary=2e-4, iso_primary_to_secondary=5e-4,
                                                                                     interleaving_type=fmt.CenterTappedInterleavingType.TypeA,
                                                                                     primary_additional_bobbin=100, winding_temperature=100)


    geo.set_insulation(insulation)
    geo.set_winding_windows([winding_window])

    geo.create_model(freq=200000, pre_visualize_geometry=show_visual_outputs)
    geo.single_simulation(freq=200000, current=[20, 120, 120], phi_deg=[0, 180, 180], show_fem_simulation_results=show_visual_outputs)

if __name__ == "__main__":
    basic_example_transformer_center_tapped(show_visual_outputs=True)