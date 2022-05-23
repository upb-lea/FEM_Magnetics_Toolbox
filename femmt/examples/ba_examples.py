import femmt as fmt

def pq4040():
    # This simulation is used for the ansys simulation comparison

    geo = fmt.MagneticComponent(component_type="inductor")

    geo.core.update(core_h=0.04, core_w=0.0149, window_h=0.0278, window_w=0.01105, mu_rel=3100, phi_mu_deg=12, sigma=0.6)
    geo.air_gaps.update(method="center", n_air_gaps=1, air_gap_h=[0.0005], position_tag=[0])

    geo.update_conductors(n_turns=[[8]], conductor_type=["solid"], conductor_radii=[0.0015],
                        winding=["primary"], scheme=["square"],
                        core_cond_isolation=[0.001, 0.001, 0.001, 0.001], cond_cond_isolation=[0.0001],
                        conductivity_sigma=["copper"])

    geo.create_model(freq=100000, visualize_before=True, save_png=False)
    #geo.single_simulation(freq=100000, current=[3], show_results=False)

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
            "isolation": 1.57
    }

    case_gap_top = 0.001
    case_gap_right = 0.001
    case_gap_bot = 0.001

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

def basic_example():
    geo = fmt.MagneticComponent(component_type="transformer")

    # Update geometry
    geo.core.update(window_h=0.0295, window_w=0.012, core_w=0.015)

    # Add air gaps
    geo.air_gaps.update(method="percent", n_air_gaps=1, air_gap_h=[0.0005],
                    air_gap_position=[50], position_tag=[0])

    # Add conductors
    geo.update_conductors(n_turns=[[21], [7]], conductor_type=["solid", "solid"],
                        litz_para_type=['implicit_litz_radius', 'implicit_litz_radius'],
                        ff=[None, None], strands_numbers=[None, None], strand_radii=[None, None],
                        conductor_radii=[0.0011, 0.0011],
                        winding=["interleaved"], scheme=["horizontal"],
                        core_cond_isolation=[0.001, 0.001, 0.002, 0.001], cond_cond_isolation=[0.0002, 0.0002, 0.0005],
                        conductivity_sigma=["copper", "copper"])

    # Create model
    geo.create_model(freq=250000, visualize_before=True)

    #geo.single_simulation(freq=250000, current=[4, 12], phi_deg=[0, 180], show_results=True)

    thermal_conductivity_dict = {
        "air": 0.0263,
        "case": {
            "top": 1.54,
            "top_right": 1.54,
            "right": 1.54,
            "bot_right": 1.54,
            "bot": 1.54
        },
        "core": 5,
        "winding": 400,
        "air_gaps": 180,
        "isolation": 0.42
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

    geo.thermal_simulation(thermal_conductivity_dict, boundary_temperatures, boundary_flags, case_gap_top, case_gap_right, case_gap_bot, show_results=True, visualize_before=True)


if __name__ == "__main__":
    pq4040()
    #basic_example()