import femmt as fmt

def example_thermal_simulation(geo):
    # Thermal simulation:
    # The losses calculated by the magnetics simulation can be used to calculate the heat distribution of the given magnetic component
    # In order to use the thermal simulation, thermal conductivities for each material can be entered as well as a boundary temperature
    # which will be applied on the boundary of the simulation (dirichlet boundary condition).
    
    # The case parameter sets the thermal conductivity for a case which will be set around the core.
    # This could model some case in which the transformer is placed in together with a set potting material.
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

    # Here the case size can be determined
    case_gap_top = 0.001
    case_gap_right = 0.001
    case_gap_bot = 0.001

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
    boundary_flags = {
        "flag_boundary_top": 1,
        "flag_boundary_top_right": 1,
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
    geo.thermal_simulation(thermal_conductivity_dict, boundary_temperatures, boundary_flags, case_gap_top, case_gap_right, case_gap_bot, True)
    #geo.femm_thermal_validation(thermal_conductivity_dict, femm_boundary_temperature, case_gap_top, case_gap_right, case_gap_bot)

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
    example_thermal_simulation(geo)

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

    geo.thermal_simulation(thermal_conductivity_dict, boundary_temperatures, boundary_flags, case_gap_top, case_gap_right, case_gap_bot, True)


if __name__ == "__main__":
    pq4040()
    #basic_example()