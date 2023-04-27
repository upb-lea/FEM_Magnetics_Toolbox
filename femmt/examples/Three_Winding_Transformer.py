import femmt as fmt
import os

def example_thermal_simulation():
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
                           case_gap_right, case_gap_bot, True, color_scheme=fmt.colors_ba_jonas, colors_geometry=fmt.colors_geometry_ba_jonas)

example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
if not os.path.exists(example_results_folder):
    os.mkdir(example_results_folder)

#component = "inductor"
#component = "transformer-interleaved"
component = "transformer"
#component = "integrated_transformer"
# component = "load_from_file"

# Create Object
if component == "inductor":
    # Working directory can be set arbitrarily
    working_directory = os.path.join(example_results_folder, "inductor")
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    # 1. chose simulation type
    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=working_directory, silent=True)

    # 2. set core parameters
    core_db = fmt.core_database()["PQ 40/40"]
    inductor_frequency = 270000
    core = fmt.Core(core_inner_diameter=core_db["core_inner_diameter"], window_w=core_db["window_w"], window_h=core_db["window_h"],
                    material="N49", temperature=45, frequency=inductor_frequency,
                    # permeability_datasource="manufacturer_datasheet",
                    permeability_datasource=fmt.MaterialDataSource.Measurement,
                    permeability_datatype=fmt.MeasurementDataType.ComplexPermeability,
                    permeability_measurement_setup="LEA_LK",
                    permittivity_datasource=fmt.MaterialDataSource.Measurement,
                    permittivity_datatype=fmt.MeasurementDataType.ComplexPermittivity,
                    permittivity_measurement_setup="LEA_LK")
                    # mu_rel=3000, phi_mu_deg=10,
                    # sigma=0.5)
    geo.set_core(core)

    # 3. set air gap parameters
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
    #air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 90)
    geo.set_air_gaps(air_gaps)

    # 4. set insulations
    insulation = fmt.Insulation()
    insulation.add_core_insulations(0.001, 0.001, 0.004, 0.001)
    insulation.add_winding_insulations([0.0005], 0.0001)
    geo.set_insulation(insulation)

    # 5. create winding window and virtual winding windows (vww)
    winding_window = fmt.WindingWindow(core, insulation)
    vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

    # 6. create conductor and set parameters: use solid wires
    winding = fmt.Conductor(0, fmt.Conductivity.Copper)
    winding.set_solid_round_conductor(conductor_radius=0.0013, conductor_arrangement=fmt.ConductorArrangement.Square)
    winding.parallel = False  # set True to make the windings parallel, currently only for solid conductors
    #winding.set_litz_round_conductor(conductor_radius=0.0013, number_strands=150, strand_radius=100e-6,
    # fill_factor=None, conductor_arrangement=fmt.ConductorArrangement.Square)

    # 7. add conductor to vww and add winding window to MagneticComponent
    vww.set_winding(winding, 9, None)
    geo.set_winding_window(winding_window)

    # 8. create the model
    geo.create_model(freq=inductor_frequency, pre_visualize_geometry=True, save_png=False)

    # 6.a. start simulation
    geo.single_simulation(freq=inductor_frequency, current=[4.5],
                          plot_interpolation=False, show_fem_simulation_results=True)

    # geo.femm_reference(freq=inductor_frequency, current=[4.5], sign=[1], non_visualize=0)

    # 6.b. Excitation Sweep Example
    # Perform a sweep using more than one frequency
    # fs = [0, 10000, 30000, 60000, 100000, 150000]
    # amplitude_list = [[10], [2], [1], [0.5], [0.2], [0.1]]
    # phase_list = [[0], [10], [20], [30], [40], [50]]
    # geo.excitation_sweep(frequency_list=fs, current_list_list=amplitude_list, phi_deg_list_list=phase_list)
    # 9. start simulation

    # 7. prepare and start thermal simulation
    # example_thermal_simulation()

if component == "transformer-interleaved":
    working_directory = os.path.join(example_results_folder, "transformer-interleaved")
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    # 1. chose simulation type
    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Transformer, working_directory=working_directory)

    # 2. set core parameters
    core = fmt.Core(window_h=0.0295, window_w=0.012, core_inner_diameter=0.015,
                    non_linear=False, sigma=1, re_mu_rel=3200, phi_mu_deg=10, permeability_datasource = fmt.MaterialDataSource.Custom, permittivity_datasource = fmt.MaterialDataSource.Custom)

    geo.set_core(core)

    # 3. set air gap parameters
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
    geo.set_air_gaps(air_gaps)

    # 4. set insulations
    insulation = fmt.Insulation()
    insulation.add_core_insulations(0.001, 0.001, 0.002, 0.001)
    insulation.add_winding_insulations([0.0002, 0.0002], 0.0001)
    geo.set_insulation(insulation)

    # 5. create winding window and virtual winding windows (vww)
    winding_window = fmt.WindingWindow(core, insulation)
    vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

    # 6. create conductors and set parameters
    winding1 = fmt.Conductor(0, fmt.Conductivity.Copper)
    winding1.set_solid_round_conductor(0.0011, None)

    winding2 = fmt.Conductor(1, fmt.Conductivity.Copper)
    winding2.set_solid_round_conductor(0.0011, None)

    # 7. add conductor to vww and add winding window to MagneticComponent
    vww.set_interleaved_winding(winding1, 21, winding2, 7, fmt.InterleavedWindingScheme.HorizontalAlternating, 0.0005)
    geo.set_winding_window(winding_window)

    # 8. start simulation with given frequency, currents and phases
    geo.create_model(freq=250000, pre_visualize_geometry=True)
    geo.single_simulation(freq=250000, current=[4, 12], phi_deg=[0, 180], show_fem_simulation_results=True)

    # other simulation options:
    # ------------------------
    # read inductances
    # geo.get_inductances(I0=8, op_frequency=250000, skin_mesh_factor=0.5)

    # 9. start thermal simulation
    #example_thermal_simulation()
    
if component == "transformer":
    # Example for a transformer with multiple virtual winding windows.
    working_directory = os.path.join(example_results_folder, "transformer")
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    # 1. chose simulation type
    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Transformer, working_directory=working_directory, silent=True)

    # 2. set core parameters
    core_dimensions = fmt.dtos.SingleCoreDimensions(window_h=0.06, window_w=0.03, core_inner_diameter=0.015)
    core = fmt.Core(core_dimensions=core_dimensions, mu_r_abs=3100, phi_mu_deg=12, sigma=1.2,
                    permeability_datasource=fmt.MaterialDataSource.Custom, permittivity_datasource=fmt.MaterialDataSource.Custom)
    geo.set_core(core)

    # 3. set air gap parameters
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
    geo.set_air_gaps(air_gaps)

    # 4. set insulation
    insulation = fmt.Insulation()
    insulation.add_core_insulations(0.001, 0.001, 0.002, 0.001)
    insulation.add_winding_insulations([0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002], 0.0005)
    geo.set_insulation(insulation)

    # 5. create winding window and virtual winding windows (vww)
    winding_window = fmt.WindingWindow(core, insulation)
    top_left, top_right, mid_top_left, mid_top_right, mid_left, mid_right, mid_bot_left, mid_bot_right,  bot_left, bot_right = winding_window.split_window(fmt.WindingWindowSplit.TenCells_Split, horizontal_split_factor=0.2, horizontal_split_factor_2=0.4,
                                                                    horizontal_split_factor_3=0.6, horizontal_split_factor_4=0.8 )
    #mid_left = winding_window.combine_vww(mid_bot_left, mid_left)
    #top_left = winding_window.combine_vww(top_left, bot_left)
    #top_right = winding_window.combine_vww(bot_right, top_right)
    #top_right = winding_window.combine_vww(top_right, bot_right)
    # 6. create conductors and set parameters

    #winding1 = fmt.Conductor(0, fmt.Conductivity.Copper)
    #winding1.set_litz_round_conductor(0.0011, 50, 0.00011, None, fmt.ConductorArrangement.Square)

    winding1 = fmt.Conductor(0, fmt.Conductivity.Copper)
    winding1.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

    winding2 = fmt.Conductor(1, fmt.Conductivity.Copper)
    winding2.set_litz_round_conductor(0.0011, 50, 0.00011, None, fmt.ConductorArrangement.Square)

    #winding2 = fmt.Conductor(1, fmt.Conductivity.Copper)
    #winding2.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

    winding3 = fmt.Conductor(2, fmt.Conductivity.Copper)
    winding3.set_litz_round_conductor(0.0011, 50, 0.00011, None, fmt.ConductorArrangement.Square)

    #winding3 = fmt.Conductor(2, fmt.Conductivity.Copper)
    #winding3.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

    winding4 = fmt.Conductor(3, fmt.Conductivity.Copper)
    winding4.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Hexagonal)
    #winding4 = fmt.Conductor(3, fmt.Conductivity.Copper)
    #winding4.set_litz_round_conductor(0.0011, 50, 0.000011, None, fmt.ConductorArrangement.Square)

    winding5 = fmt.Conductor(4, fmt.Conductivity.Copper)
    winding5.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

    winding6 = fmt.Conductor(5, fmt.Conductivity.Copper)
    winding6.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

    #winding7 = fmt.Conductor(6, fmt.Conductivity.Copper)
    #winding7.set_litz_round_conductor(0.0011, 5, 0.00011, None, fmt.ConductorArrangement.Square)

    winding7 = fmt.Conductor(6, fmt.Conductivity.Copper)
    winding7.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

    winding8 = fmt.Conductor(7, fmt.Conductivity.Copper)
    winding8.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

    #winding8 = fmt.Conductor(7, fmt.Conductivity.Copper)
    #winding8.set_litz_round_conductor(0.0011, 50, 0.000011, None, fmt.ConductorArrangement.Square)

    winding9 = fmt.Conductor(8, fmt.Conductivity.Copper)
    winding9.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

    winding10 = fmt.Conductor(9, fmt.Conductivity.Copper)
    winding10.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)



    # 7. add conductor to vww and add winding window to MagneticComponent
    top_left.set_winding(winding1, 15, fmt.WindingType.Single)
    top_right.set_winding(winding2, 10, fmt.WindingType.Single)
    mid_top_left.set_winding(winding3, 10, fmt.WindingType.Single)
    mid_top_right.set_winding(winding4, 10, fmt.WindingType.Single)
    mid_left.set_winding(winding5, 10, fmt.WindingType.Single)
    mid_right.set_winding(winding6, 10, fmt.WindingType.Single)
    mid_bot_left.set_winding(winding7, 10, fmt.WindingType.Single)
    mid_bot_right.set_winding(winding8, 10, fmt.WindingType.Single)
    bot_left.set_winding(winding9, 8, fmt.WindingType.Single)
    bot_right.set_winding(winding10, 12, fmt.WindingType.Single)
    geo.set_winding_windows([winding_window])
    # 8. start simulation with given frequency, currents and phases
    geo.create_model(freq=250000, pre_visualize_geometry=True)
    geo.single_simulation(freq=250000, current=[4, 4, 4, 4, 4, 4, 4, 4, 4, 4], phi_deg=[0, 180, 180, 180, 180, 180, 180, 180, 180, 180])

    # read inductances
    #geo.get_inductances(I0=4, op_frequency=250000, skin_mesh_factor=0.5, visualize=False)
    # Reference simulation using FEMM
    geo.femm_reference(freq=250000, current=[4, 4, 4, 4, 4, 4, 4, 4, 4, 4], sign=[1, -1, -1, -1, -1, -1, -1, -1, -1, -1], non_visualize=0)
if component == "integrated_transformer":
    working_directory = os.path.join(example_results_folder, "integrated-transformer")
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    # 1. chose simulation type
    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.IntegratedTransformer, working_directory=working_directory)

    # 2. set core parameters
    core = fmt.Core(window_h=0.03, window_w=0.011, core_inner_diameter=0.02,
                    mu_r_abs=3100, phi_mu_deg=12,
                    sigma=0.6, permeability_datasource = fmt.MaterialDataSource.Custom, permittivity_datasource = fmt.MaterialDataSource.Custom)
    geo.set_core(core)

    # 2.1 set stray path parameters
    stray_path = fmt.StrayPath(start_index=0, length=geo.core.core_inner_diameter / 2 + geo.core.window_w - 0.001)
    geo.set_stray_path(stray_path)

    # 3. set air gap parameters
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.001, 30)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.001, 40)
    geo.set_air_gaps(air_gaps)

    # 4. set insulations
    insulation = fmt.Insulation()
    insulation.add_core_insulations(0.001, 0.001, 0.002, 0.001)
    insulation.add_winding_insulations([0.0002, 0.0002], 0.0001)
    geo.set_insulation(insulation)

    # 5. create winding window and virtual winding windows (vww)
    # For an integrated transformer it is not necessary to set horizontal and vertical split factors
    # since this is determined by the stray_path
    winding_window = fmt.WindingWindow(core, insulation, stray_path, air_gaps)
    top, bot = winding_window.split_window(fmt.WindingWindowSplit.HorizontalSplit)

    # 6. set conductor parameters
    winding1 = fmt.Conductor(0, fmt.Conductivity.Copper)
    winding1.set_litz_round_conductor(None, 100, 70e-6, 0.5, fmt.ConductorArrangement.Square)

    winding2 = fmt.Conductor(1, fmt.Conductivity.Copper)
    winding2.set_litz_round_conductor(None, 100, 70e-6, 0.5, fmt.ConductorArrangement.Square)

    # 7. add conductor to vww and add winding window to MagneticComponent
    top.set_interleaved_winding(winding1, 3, winding2, 6, fmt.InterleavedWindingScheme.HorizontalAlternating, 0.0005)
    bot.set_interleaved_winding(winding1, 1, winding2, 2, fmt.InterleavedWindingScheme.HorizontalAlternating, 0.0005)
    geo.set_winding_window(winding_window)

    # 8. start simulation with given frequency, currents and phases
    geo.create_model(freq=250000, pre_visualize_geometry=True)
    geo.single_simulation(freq=250000, current=[8.0, 4.0], phi_deg=[0, 180])

    # other simulation options:
    # -------------------------
    # geo.get_inductances(I0=10, op_frequency=100000, skin_mesh_factor=0.5)

if component == "load_from_file":
    working_directory = os.path.join(example_results_folder, "from-file")
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    file = os.path.join(os.path.dirname(__file__), "example_log.json")

    geo = fmt.MagneticComponent.decode_settings_from_log(file, working_directory)

    geo.create_model(freq=100000, pre_visualize_geometry=False, save_png=False)

    geo.single_simulation(freq=100000, current=[4.5], show_fem_simulation_results=True)
