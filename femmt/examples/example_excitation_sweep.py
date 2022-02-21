# minimal example for github clone
from femmt import *
import numpy as np

# component = "inductor"
component = "transformer"

# Create Object
if component == "inductor":
    geo = MagneticComponent(component_type="inductor", working_directory="")

    frequencies = [100000, 200000]
    current_amplitudes = [[10], [4]]
    phases = [[0], [179]]

    # Update Geometry
    geo.core.update(window_h=0.03, window_w=0.011)

    # geo.air_gaps.update(method="percent", n_air_gaps=4, air_gap_h=[0.0005, 0.0005, 0.0005, 0.0005],
    #                     position_tag=[0, 0, 0, 0], air_gap_position=[20, 40, 60, 80])
    geo.air_gaps.update(method="center", n_air_gaps=1, air_gap_h=[0.0005], position_tag=[0])

    geo.update_conductors(n_turns=[[17]], conductor_type=["litz"], conductor_radii=[0.0012],
                          litz_para_type=['implicit_ff'], strands_numbers=[500], strand_radii=[35e-6],
                          winding=["primary"], scheme=["square_full_width"],
                          core_cond_isolation=[0.0005, 0.0005, 0.0005, 0.0005], cond_cond_isolation=[0.0002])

    # Perform a frequency sweep simulation
    geo.visualize_before = False
    geo.excitation_sweep(frequency_list=frequencies, current_list_list=current_amplitudes, phi_deg_list_list=phases, show_last=True)


if component == "transformer":
    geo = MagneticComponent(component_type="transformer", working_directory="")

    frequencies = [100000, 200000]
    current_amplitudes = [[4, 14.5], [2, 6]]
    phases = [[0, 176], [0, 163]]

    # Update Geometry
    geo.core.update(window_h=0.0295, window_w=0.012, core_w=0.015)

    # geo.air_gaps.update(n_air_gaps=0)
    geo.air_gaps.update(method="percent", n_air_gaps=1, air_gap_h=[0.0005],
                        air_gap_position=[50], position_tag=[0])

    geo.update_conductors(n_turns=[[36], [11]], conductor_type=["litz", "litz"],
                          litz_para_type=['implicit_ff', 'implicit_ff'],
                          strands_numbers=[400, 600], strand_radii=[35e-6, 35e-6],
                          conductor_radii=[0.0010, 0.0012],
                          winding=["interleaved"], scheme=["horizontal"],
                          core_cond_isolation=[0.0005, 0.0005, 0.0005, 0.0005], cond_cond_isolation=[0.0002, 0.0002, 0.0005])

    # Perform a frequency sweep simulation
    geo.visualize_before = True
    geo.excitation_sweep(frequency_list=frequencies, current_list_list=current_amplitudes, phi_deg_list_list=phases, show_last=True)
