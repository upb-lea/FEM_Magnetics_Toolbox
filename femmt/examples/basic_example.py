# minimal example for github clone
from femmt import *
import numpy as np

# component = "inductor"
component = "transformer"
# component = "integrated_transformer"


# Create Object
if component == "inductor":
    geo = MagneticComponent(component_type="inductor")

    # Update Geometry
    geo.core.update(window_h=0.03, window_w=0.011)

    # geo.air_gaps.update(method="percent", n_air_gaps=4, air_gap_h=[0.0005, 0.0005, 0.0005, 0.0005],
    #                     position_tag=[0, 0, 0, 0], air_gap_position=[20, 40, 60, 80])
    geo.air_gaps.update(method="center", n_air_gaps=1, air_gap_h=[0.0005], position_tag=[0])

    geo.update_conductors(n_turns=[[14]], conductor_type=["solid"], conductor_radii=[0.0015],
                          winding=["primary"], scheme=["square"],
                          core_cond_isolation=[0.0005], cond_cond_isolation=[0.0001])

    # geo.single_simulation(freq=100000, current=[1])
    geo.femm_reference(freq=100000, current=[1], sigma_cu=58, sign=[1], non_visualize=0)

if component == "transformer":
    simulate_before_thermal = True

    geo = MagneticComponent(component_type="transformer")
    geo.visualize_before = False

    # Update Geometry
    geo.core.update(window_h=0.0295, window_w=0.012, core_w=0.015)

    # geo.air_gaps.update(n_air_gaps=0)
    geo.air_gaps.update(method="percent", n_air_gaps=1, air_gap_h=[0.0005],
                        air_gap_position=[50], position_tag=[0])

    geo.update_conductors(n_turns=[[36], [11]], conductor_type=["solid", "litz"],
                        litz_para_type=['implicit_litz_radius', 'implicit_litz_radius'],
                        ff=[None, 0.6], strands_numbers=[None, 600], strand_radii=[70e-6, 35.5e-6],
                        conductor_radii=[0.0011, None],
                        winding=["interleaved"], scheme=["horizontal"],
                        core_cond_isolation=[0.0005, 0.0005], cond_cond_isolation=[0.0002, 0.0002, 0.0005])

    if simulate_before_thermal:
        # Perform a single simulation
        geo.single_simulation(freq=250000, current=[4.14723021, 14.58960019], phi_deg=[- 1.66257715/np.pi*180, 170])
        # geo.single_simulation(freq=250000, current=[4.18368713, 4.28975166], phi_deg=[-1.09710805/np.pi*180,
        #                                                                               - 1.47917789/np.pi*180 + 180])

        # geo.get_inductances(I0=8, op_frequency=250000, skin_mesh_factor=0.5)
        # geo.femm_reference(freq=100000, current=[1, 2], sigma_cu=58, sign=[1, -1], non_visualize=0)

    thermal_conductivity_dict = {
            "air": 0.0263,
            "case": 0.3,
            "core": 5,
            "winding": 400 
    }

    geo.thermal_simulation(thermal_conductivity_dict)
    geo.femm_thermal_validation(thermal_conductivity_dict)

if component == "integrated_transformer":
    geo = MagneticComponent(component_type="integrated_transformer")

    geo.update_conductors(n_turns=[[1, 3], [2, 6]], conductor_type=["litz", "litz"],
                          litz_para_type=['implicit_litz_radius', 'implicit_litz_radius'],
                          ff=[0.5, 0.5], strands_numbers=[100, 100], strand_radii=[70e-6, 70e-6],
                          winding=["interleaved", "interleaved"], scheme=["horizontal", "horizontal"],
                          core_cond_isolation=[0.0005, 0.0005], cond_cond_isolation=[0.0002, 0.0002, 0.0005])

    geo.core.update(window_h=0.03, window_w=0.011)

    geo.air_gaps.update(method="percent",
                        n_air_gaps=2,
                        position_tag=[0, 0],
                        air_gap_h=[0.001, 0.001],
                        air_gap_position=[30, 40])

    geo.stray_path.update(start_index=0,
                          radius=geo.core.core_w / 2 + geo.core.window_w - 0.001)

    # Perform a single simulation
    geo.single_simulation(freq=250000, current=[8.0, 4.0], phi_deg=[0, 180])
    # geo.get_inductances(I0=10, op_frequency=100000, skin_mesh_factor=0.5)

