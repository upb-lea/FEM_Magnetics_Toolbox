import femmt as fmt
import os



def basic_example_transformer_n_winding(onelab_folder: str = None, show_visual_outputs: bool = True, is_test: bool = False):
    def example_thermal_simulation(show_visual_outputs: bool = True):


        # Thermal simulation:
        # The losses calculated by the magnetics simulation can be used to calculate the heat distribution of the given magnetic component
        # In order to use the thermal simulation, thermal conductivities for each material can be entered as well as a boundary temperature
        # which will be applied on the boundary of the simulation (dirichlet boundary condition).

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
            "air_gaps": 180,  # aluminiumnitride
            "insulation": 0.42  # polyethylen
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
                               case_gap_right, case_gap_bot, show_visual_outputs, color_scheme=fmt.colors_ba_jonas,
                               colors_geometry=fmt.colors_geometry_ba_jonas)


    example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
    if not os.path.exists(example_results_folder):
        os.mkdir(example_results_folder)

    # Example for a transformer with multiple virtual winding windows.
    working_directory = os.path.join(example_results_folder, "n-winding-transformer")
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    # 1. chose simulation type
    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Transformer, working_directory=working_directory,
                                verbosity=fmt.Verbosity.Silent)

    # 2. set core parameters
    core_dimensions = fmt.dtos.SingleCoreDimensions(window_h=0.12, window_w=0.09, core_inner_diameter=0.050, core_h= 0.2)
    core = fmt.Core(core_dimensions=core_dimensions, mu_r_abs=3100, phi_mu_deg=12, sigma=1.2,
                    permeability_datasource=fmt.MaterialDataSource.Custom,
                    permittivity_datasource=fmt.MaterialDataSource.Custom)
    geo.set_core(core)

    # 3. set air gap parameters
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
    geo.set_air_gaps(air_gaps)

    # 4. set insulation
    insulation = fmt.Insulation()
    insulation.add_core_insulations(0.001, 0.001, 0.002, 0.001)
    # insulation.add_winding_insulations([0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002], 0.0005)
    insulation.add_winding_insulations(
        [[0.0002, 0.0004, 0.0004, 0.0002, 0.0002, 0.0004, 0.0002, 0.0002, 0.0002, 0.0004, 0.0002, 0.0004],
         [0.0004, 0.0002, 0.0004, 0.0002, 0.0002, 0.0004, 0.0002, 0.0002, 0.0002, 0.0004, 0.0002, 0.0004],
         [0.0004, 0.0002, 0.0004, 0.0002, 0.0002, 0.0004, 0.0002, 0.0002, 0.0002, 0.0004, 0.0002, 0.0004],
         [0.0002, 0.0004, 0.0004, 0.0002, 0.0002, 0.0004, 0.0002, 0.0002, 0.0002, 0.0004, 0.0002, 0.0004],
         [0.0004, 0.0002, 0.0004, 0.0002, 0.0002, 0.0004, 0.0002, 0.0002, 0.0002, 0.0004, 0.0002, 0.0004],
         [0.0004, 0.0002, 0.0004, 0.0002, 0.0002, 0.0004, 0.0002, 0.0002, 0.0002, 0.0004, 0.0002, 0.0004],
         [0.0002, 0.0004, 0.0004, 0.0002, 0.0002, 0.0004, 0.0002, 0.0002, 0.0002, 0.0004, 0.0002, 0.0004],
         [0.0004, 0.0002, 0.0004, 0.0002, 0.0002, 0.0004, 0.0002, 0.0002, 0.0002, 0.0004, 0.0002, 0.0004],
         [0.0004, 0.0002, 0.0004, 0.0002, 0.0002, 0.0004, 0.0002, 0.0002, 0.0002, 0.0004, 0.0002, 0.0004],
         [0.0002, 0.0004, 0.0004, 0.0002, 0.0002, 0.0004, 0.0002, 0.0002, 0.0002, 0.0004, 0.0002, 0.0004],
         [0.0004, 0.0002, 0.0004, 0.0002, 0.0002, 0.0004, 0.0002, 0.0002, 0.0002, 0.0004, 0.0002, 0.0004],
         [0.0004, 0.0002, 0.0004, 0.0002, 0.0002, 0.0004, 0.0002, 0.0002, 0.0002, 0.0004, 0.0002, 0.0004]])
    geo.set_insulation(insulation)

    # 5. create winding window and virtual winding windows (vww)
    winding_window = fmt.WindingWindow(core, insulation)
    cells = winding_window.NCellsSplit(0, [1 / 6, 2 / 6, 3 / 6, 4 / 6, 5 / 6])
    # 6. create conductors and set parameters

    winding1 = fmt.Conductor(0, fmt.Conductivity.Copper)
    winding1.set_litz_round_conductor(0.0011, 50, 0.00011, None, fmt.ConductorArrangement.Square)

    # winding1 = fmt.Conductor(0, fmt.Conductivity.Copper)
    # winding1.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

    # winding2 = fmt.Conductor(1, fmt.Conductivity.Copper)
    # winding2.set_litz_round_conductor(0.0011, 50, 0.00011, None, fmt.ConductorArrangement.Square)

    winding2 = fmt.Conductor(1, fmt.Conductivity.Copper)
    winding2.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

    # winding3 = fmt.Conductor(2, fmt.Conductivity.Copper)
    # winding3.set_litz_round_conductor(0.0011, 50, 0.00011, None, fmt.ConductorArrangement.Square)

    winding3 = fmt.Conductor(2, fmt.Conductivity.Copper)
    winding3.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

    winding4 = fmt.Conductor(3, fmt.Conductivity.Copper)
    winding4.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)
    # winding4 = fmt.Conductor(3, fmt.Conductivity.Copper)
    # winding4.set_litz_round_conductor(0.0011, 50, 0.000011, None, fmt.ConductorArrangement.Square)

    winding5 = fmt.Conductor(4, fmt.Conductivity.Copper)
    winding5.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

    winding6 = fmt.Conductor(5, fmt.Conductivity.Copper)
    winding6.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

    winding7 = fmt.Conductor(6, fmt.Conductivity.Copper)
    winding7.set_litz_round_conductor(0.0011, 5, 0.00011, None, fmt.ConductorArrangement.Square)

    winding7 = fmt.Conductor(6, fmt.Conductivity.Copper)
    winding7.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

    winding8 = fmt.Conductor(7, fmt.Conductivity.Copper)
    winding8.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

    # winding8 = fmt.Conductor(7, fmt.Conductivity.Copper)
    # winding8.set_litz_round_conductor(0.0011, 50, 0.000011, None, fmt.ConductorArrangement.Square)

    winding9 = fmt.Conductor(8, fmt.Conductivity.Copper)
    winding9.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

    winding10 = fmt.Conductor(9, fmt.Conductivity.Copper)
    winding10.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

    winding11 = fmt.Conductor(10, fmt.Conductivity.Copper)
    winding11.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

    winding12 = fmt.Conductor(11, fmt.Conductivity.Copper)
    winding12.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

    # 5. create winding window and virtual winding windows (vww)
    # winding_window = fmt.WindingWindow(core, insulation)
    # cells = winding_window.NCellsSplit(0, [1/6, 2/6, 3/6, 4/6, 5/6])
    # 7. add conductor to vww and add winding window to MagneticComponent
    cells[0].set_winding(winding1, 9, fmt.WindingType.Single)
    cells[1].set_winding(winding2, 10, fmt.WindingType.Single)
    cells[2].set_winding(winding3, 9, fmt.WindingType.Single)
    cells[3].set_winding(winding4, 7, fmt.WindingType.Single)
    cells[4].set_winding(winding5, 12, fmt.WindingType.Single)
    cells[5].set_winding(winding6, 9, fmt.WindingType.Single)
    cells[6].set_winding(winding7, 9, fmt.WindingType.Single)
    cells[7].set_winding(winding8, 12, fmt.WindingType.Single)
    cells[8].set_winding(winding9, 13, fmt.WindingType.Single)
    cells[9].set_winding(winding10, 15, fmt.WindingType.Single)
    cells[10].set_winding(winding11, 11, fmt.WindingType.Single)
    cells[11].set_winding(winding12, 15, fmt.WindingType.Single)
    geo.set_winding_windows([winding_window])
    # 8. start simulation with given frequency, currents and phases
    geo.create_model(freq=250000, pre_visualize_geometry=show_visual_outputs)
    geo.single_simulation(freq=250000, current=[4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
                          phi_deg=[0, 180, 0, 180, 0, 180, 0, 180, 0, 180, 0, 180], show_fem_simulation_results=show_visual_outputs)
    # 7. prepare and start thermal simulation
    example_thermal_simulation(show_visual_outputs=show_visual_outputs)

    # read inductances
    # geo.get_inductances(I0=4, op_frequency=250000, skin_mesh_factor=0.5, visualize=False)
    # Reference simulation using FEMM
    # geo.femm_reference(freq=250000, current=[4, 4, 4, 4, 4, 4, 4, 4, 4, 4], sign=[1, -1, -1, -1, -1, -1, -1, -1, -1, -1], non_visualize=0)

if __name__ == "__main__":
    basic_example_transformer_n_winding(show_visual_outputs=True)