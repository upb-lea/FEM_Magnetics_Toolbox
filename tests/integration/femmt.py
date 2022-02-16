import pytest
import numpy as np
from femmt import MagneticComponent

@pytest.mark.integtest
def test_femmt():

    # Build model and create test files
    geo = MagneticComponent(component_type="transformer")
    geo.visualize_before = False
    geo.core.update(window_h=0.0295, window_w=0.012, core_w=0.015)
    geo.air_gaps.update(method="percent", n_air_gaps=1, air_gap_h=[0.0005],
                        air_gap_position=[50], position_tag=[0])
    geo.update_conductors(n_turns=[[36], [11]], conductor_type=["solid", "litz"],
                        litz_para_type=['implicit_litz_radius', 'implicit_litz_radius'],
                        ff=[None, 0.6], strands_numbers=[None, 600], strand_radii=[70e-6, 35.5e-6],
                        conductor_radii=[0.0011, None],
                        winding=["interleaved"], scheme=["horizontal"],
                        core_cond_isolation=[0.0005, 0.0005], cond_cond_isolation=[0.0002, 0.0002, 0.0005])

    geo.create_model(freq=250000)
    geo.single_simulation(freq=250000, current=[4.14723021, 14.58960019], phi_deg=[- 1.66257715/np.pi*180, 170])

    thermal_conductivity_dict = {
            "air": 0.0263,
            "case": 0.3,
            "core": 5,
            "winding": 400,
            "air_gaps": 180
    }

    case_gap_top = 0.0015
    case_gap_right = 0.0025
    case_gap_bot = 0.002
    boundary_temperatures = {
        "value_boundary_top": 293,
        "value_boundary_top_right": 293,
        "value_boundary_right_top": 293,
        "value_boundary_right": 293,
        "value_boundary_right_bottom": 293,
        "value_boundary_bottom_right": 293,
        "value_boundary_bottom": 293
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

    geo.thermal_simulation(thermal_conductivity_dict, boundary_temperatures, boundary_flags, case_gap_top, case_gap_right, case_gap_bot)