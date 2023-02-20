import femmt as fmt
import os

# simulation = "lab_model"
# simulation = "pq4040_ansys_comparison"
simulation = "pq4040axisymmetric"

working_directory = os.path.join(os.path.dirname(__file__), "example_results")

if not os.path.exists(working_directory):
    os.mkdir(working_directory)

if simulation == "lab_model":
    # This simulation is used for the lab model

    cwd = os.path.join(working_directory, "lab_model")
    if not os.path.exists(cwd):
        os.mkdir(cwd)

    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=cwd)

    # 2. set core parameters
    core_db = fmt.core_database()["PQ 40/40"]
    
    core = fmt.Core(core_inner_diameter=core_db["core_inner_diameter"], window_w=core_db["window_w"], window_h=core_db["window_h"],
                    mu_r_abs=3100, phi_mu_deg=12, sigma=0.6)
    geo.set_core(core)

    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
    geo.set_air_gaps(air_gaps)

    insulation = fmt.Insulation()
    insulation.add_core_insulations(0.001, 0.001, 0.002, 0.001)
    insulation.add_winding_insulations([0.0001])
    geo.set_insulation(insulation)

    winding_window = fmt.WindingWindow(core, insulation)
    vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

    winding = fmt.Conductor(0, fmt.Conductivity.Copper)
    winding.set_litz_round_conductor(conductor_radius=0.0015, number_strands=150, strand_radius=100e-6, 
                                        fill_factor=None, conductor_arrangement=fmt.ConductorArrangement.Square)

    vww.set_winding(winding, 8, None)
    geo.set_winding_window(winding_window)

    geo.create_model(freq=100000, visualize_before=True, save_png=False)
    geo.single_simulation(freq=100000, current=[3], show_results=False)

    thermal_conductivity_dict = {
            "air": 0.122, # potting epoxy resign
            "case": {
                "top": 0.122,
                "top_right": 0.122,
                "right": 0.122,
                "bot_right": 0.122,
                "bot": 0.122
            },
            "core": 5, # ferrite
            "winding": 0.517, # copper
            "air_gaps": 1.57,
            "insulation": 1.57
    }

    case_gap_top = 0.0004
    case_gap_right = 0.01
    case_gap_bot = 0.03

    boundary_temperatures = {
        "value_boundary_top": 20,
        "value_boundary_top_right": 20,
        "value_boundary_right_top": 20,
        "value_boundary_right": 20,
        "value_boundary_right_bottom": 20,
        "value_boundary_bottom_right": 20,
        "value_boundary_bottom": 20
    }
    
    femm_boundary_temperature = 20

    boundary_flags = {
        "flag_boundary_top": 1,
        "flag_boundary_top_right": 1,
        "flag_boundary_right_top": 1,
        "flag_boundary_right": 1,
        "flag_boundary_right_bottom": 1,
        "flag_boundary_bottom_right": 1,
        "flag_boundary_bottom": 1
    }

    #color_scheme = fmt.colors_ba_jonas
    #colors_geometry = fmt.colors_geometry_ba_jonas
    color_scheme = fmt.colors_ba_jonas
    colors_geometry = fmt.colors_geometry_draw_only_lines

    geo.thermal_simulation(thermal_conductivity_dict, boundary_temperatures, boundary_flags, case_gap_top, case_gap_right, 
        case_gap_bot, show_results=True, visualize_before=True, color_scheme=color_scheme, colors_geometry=colors_geometry)

if simulation == "pq4040_ansys_comparison":
    # This simulation is used for the ansys simulation comparison

    cwd = os.path.join(working_directory, "pq4040_ansys_comparison")
    if not os.path.exists(cwd):
        os.mkdir(cwd)

    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=cwd)

    core = fmt.Core(core_h=0.04, core_inner_diameter=0.0149, window_h=0.0278, window_w=0.01105, mu_r_abs=3100, phi_mu_deg=12, sigma=0.6)
    geo.set_core(core)

    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Center, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005)
    geo.set_air_gaps(air_gaps)

    insulation = fmt.Insulation()
    insulation.add_core_insulations(0.001, 0.001, 0.002, 0.001)
    insulation.add_winding_insulations([0.0001])
    geo.set_insulation(insulation)

    winding_window = fmt.WindingWindow(core, insulation)
    vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

    winding = fmt.Conductor(0, fmt.Conductivity.Copper)
    winding.set_solid_round_conductor(conductor_radius=0.0015, conductor_arrangement=fmt.ConductorArrangement.Square)

    vww.set_winding(winding, 8, None)
    geo.set_winding_window(winding_window)

    geo.create_model(freq=100000, visualize_before=True, save_png=False)
    geo.single_simulation(freq=100000, current=[3], show_results=False)

    thermal_conductivity_dict = {
            "air": 1.57, # potting epoxy resign
            "case": {
                "top": 1.57,
                "top_right": 1.57,
                "right": 1.57,
                "bot_right": 1.57,
                "bot": 1.57
            },
            "core": 5, # ferrite
            "winding": 400, # copper
            "air_gaps": 1.57,
            "insulation": 1.57
    }

    case_gap_top = 0.0004
    case_gap_right = 0.000655
    case_gap_bot = 0.0004

    boundary_temperatures = {
        "value_boundary_top": 20,
        "value_boundary_top_right": 20,
        "value_boundary_right_top": 20,
        "value_boundary_right": 20,
        "value_boundary_right_bottom": 20,
        "value_boundary_bottom_right": 20,
        "value_boundary_bottom": 20
    }
    
    femm_boundary_temperature = 20

    boundary_flags = {
        "flag_boundary_top": 1,
        "flag_boundary_top_right": 1,
        "flag_boundary_right_top": 1,
        "flag_boundary_right": 1,
        "flag_boundary_right_bottom": 1,
        "flag_boundary_bottom_right": 1,
        "flag_boundary_bottom": 1
    }

    #color_scheme = fmt.colors_ba_jonas
    #colors_geometry = fmt.colors_geometry_ba_jonas
    color_scheme = fmt.colors_ba_jonas
    colors_geometry = fmt.colors_geometry_draw_only_lines


    geo.thermal_simulation(thermal_conductivity_dict, boundary_temperatures, boundary_flags, case_gap_top, case_gap_right, 
        case_gap_bot, show_results=True, visualize_before=True, color_scheme=color_scheme, colors_geometry=colors_geometry)
    #geo.femm_thermal_validation(thermal_conductivity_dict, femm_boundary_temperature, case_gap_top, case_gap_right, case_gap_bot)

if simulation == "pq4040axisymmetric":
    cwd = os.path.join(working_directory, "pq4040axisymmetric")
    if not os.path.exists(cwd):
        os.mkdir(cwd)

    core_db = fmt.core_database()["PQ 40/40"]

    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=cwd)
    
    core = fmt.Core(core_inner_diameter=core_db["core_inner_diameter"], window_w=core_db["window_w"], window_h=core_db["window_h"],
                    mu_r_abs=3100, phi_mu_deg=12,
                    sigma=0.)
    geo.set_core(core)

    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Center, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005)
    geo.set_air_gaps(air_gaps)

    insulation = fmt.Insulation()
    insulation.add_core_insulations(0.001, 0.001, 0.002, 0.001)
    insulation.add_winding_insulations([0.0001])
    geo.set_insulation(insulation)

    winding_window = fmt.WindingWindow(core, insulation)
    vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

    winding = fmt.Conductor(0, fmt.Conductivity.Copper)
    winding.set_solid_round_conductor(conductor_radius=0.0015, conductor_arrangement=fmt.ConductorArrangement.Square)

    vww.set_winding(winding, 8, None)
    geo.set_winding_window(winding_window)

    geo.create_model(freq=100000, visualize_before=True, save_png=False)
    geo.single_simulation(freq=100000, current=[3], show_results=True)

    thermal_conductivity_dict = {
            "air": 1.57, # potting epoxy resign
            "case": {
                "top": 1.57,
                "top_right": 1.57,
                "right": 1.57,
                "bot_right": 1.57,
                "bot": 1.57
            },
            "core": 5, # ferrite
            "winding": 400, # copper
            "air_gaps": 1.57,
            "insulation": 1.57
    }

    case_gap_top = 0.002
    case_gap_right = 0.0025
    case_gap_bot = 0.002

    boundary_temperatures = {
        "value_boundary_top": 20,
        "value_boundary_top_right": 20,
        "value_boundary_right_top": 20,
        "value_boundary_right": 20,
        "value_boundary_right_bottom": 20,
        "value_boundary_bottom_right": 20,
        "value_boundary_bottom": 20
    }

    boundary_flags = {
        "flag_boundary_top": 1,
        "flag_boundary_top_right": 1,
        "flag_boundary_right_top": 1,
        "flag_boundary_right": 1,
        "flag_boundary_right_bottom": 1,
        "flag_boundary_bottom_right": 1,
        "flag_boundary_bottom": 1
    }

    #color_scheme = fmt.colors_ba_jonas
    #colors_geometry = fmt.colors_geometry_ba_jonas
    color_scheme = fmt.colors_ba_jonas
    colors_geometry = fmt.colors_geometry_draw_only_lines

    geo.thermal_simulation(thermal_conductivity_dict, boundary_temperatures, boundary_flags, case_gap_top, case_gap_right, 
        case_gap_bot, show_results=True, visualize_before=True, color_scheme=color_scheme, colors_geometry=colors_geometry)
