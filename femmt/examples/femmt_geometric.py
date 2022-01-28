from femmt import MagneticComponent
import numpy as np

# -- Definitions --
frequency = 100000
flag_reference = 1
inductance_calculation = 0
skin_accuracy = 1 # 0.5
N1 = 6
N2 = 6
I0 = 20

# -- Component Selection --
geo = MagneticComponent(component_type="transformer")
# geo = MagneticComponent(component_type="inductor")

# -- Core --
geo.update_core(core_type="EI", window_h=0.03, window_w=0.01, core_w=0.02)

# -- Air Gaps --
"""
geo.update_air_gaps(method="percent",
                    n_air_gaps=4,
                    position_tag=[0, 0, 0, 0],
                    air_gap_h=[0.0005, 0.0005, 0.0005, 0.0005],
                    air_gap_position=[20, 40, 60, 80],
                    stray_path=[[2, 3], geo.core_w/2+geo.window_w-0.001 ])
"""
geo.update_air_gaps(method="center", n_air_gaps=1, air_gap_h=[0.002])


# -- Simulation --
if geo.component_type == "transformer":
    geo.update_conductors(n_turns=[N1, N2],
                          conductor_type=["litz", "litz"],
                          conductor_radix=[0.0015, 0.0015],
                          scheme=["square", "square"],
                          layer_numbers=[10, 10],
                          strand_radix=[0.05e-3, 0.05e-3],
                          thickness=[0.000975, 0.0015],
                          wrap_para=["fixed_thickness", "interpolate"],
                          cond_cond_isolation=[0.0002, 0.0002, 0.0008],
                          core_cond_isolation=[0.001, 0.001],
                          )

    # Start Simulation
    if inductance_calculation and not flag_reference:
        geo.get_inductances(I0=I0, op_frequency=frequency, skin_mesh_factor=skin_accuracy)
    else:
        if flag_reference:
            geo.femm_reference(freq=frequency, current=[0, I0], sign=[1, -1], sigma=58, non_visualize=0)
        else:
            # geo.adaptive_single_simulation(freq=1000000, current=[10], max_iter=10, local=1)
            geo.single_simulation(freq=frequency, current=[I0, N1/N2*I0], phi=[0, 0], skin_mesh_factor=skin_accuracy)
            # geo.get_Core_Loss(Ipeak=[I0, N1/N2*I0], ki=0.53, alpha=1.50, beta=2.38, t_rise=0.3/frequency,
            #                  t_fall=0.3/frequency, f_switch=frequency, skin_mesh_factor=skin_accuracy)
            geo.get_Core_Loss(Ipeak=[I0, N1/N2*I0], ki=0.23, alpha=1.80, beta=2.88, t_rise=0.3/frequency,
                              t_fall=0.3/frequency, f_switch=frequency, skin_mesh_factor=skin_accuracy)

if geo.component_type == "inductor":
    geo.update_conductors(n_turns=[1],
                          conductor_type=["full"],
                          conductor_radix=[0.0015],
                          scheme=["hexa"],
                          wrap_para=["fixed_thickness"],
                          thickness=[0.000975],
                          cond_cond_isolation=[0.0002],
                          core_cond_isolation=[0.0004],
                          )
    """
    geo.update_conductors(n_turns=[6],
                          conductor_type=["litz"],
                          conductor_radix=[0.0015],
                          layer_numbers=[11],
                          strand_radix=[0.05e-3])
    """
    # Start Simulation
    if flag_reference:
        geo.femm_reference(freq=frequency, current=[1], sigma=58, non_visualize=0)
    else:
        geo.single_simulation(freq=frequency, current=[1], skin_mesh_factor=accuracy)
