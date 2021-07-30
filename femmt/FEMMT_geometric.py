from FEMMT import MagneticComponent

# -- Definitions --
frequency = 100000
flag_reference = 0
inductance_calculation = 1
accuracy = 0.5
N1 = 8
N2 = 8
I0 = 1

# -- Component Selection --
geo = MagneticComponent(component_type="transformer")
# geo = MagneticComponent(component_type="inductor")

# -- Core --
geo.update_core(core_type="EI", window_h=0.03)

# -- Air Gaps --
# geo.update_air_gaps(method="percent", n_air_gaps=0, position_tag=[0, 0],
# air_gap_h=[0.001, 0.001], air_gap_position=[20, 80])
geo.update_air_gaps(method="center", n_air_gaps=1, air_gap_h=[0.002])

# -- Simulation --
if geo.component_type == "transformer":
    geo.update_conductors(n_turns=[N1, N2],
                          conductor_type=["solid", "solid"],
                          conductor_radix=[0.0015, 0.0015],
                          layer_numbers=[10, 10],
                          strand_radix=[0.05e-3, 0.05e-3])
    """
    geo.update_conductors(n_turns=[10, 10],      # 0.0011879; 0.0015
                          conductor_type=["solid", "solid"],
                          conductor_radix=[0.0015, 0.0015],
                          layer_numbers=[10, 10],
                          strand_radix=[0.05e-3, 0.05e-3])
    """
    # Start Simulation
    if inductance_calculation and not flag_reference:
        geo.get_inductances(I0=I0, op_frequency=frequency, mesh_accuracy=accuracy)
    else:
        if flag_reference:
            geo.femm_reference(freq=frequency, current=[0, I0], sign=[1, -1], sigma=58, non_visualize=0)
        else:
            # geo.adaptive_single_simulation(freq=1000000, current=[10], max_iter=10, local=1)
            geo.single_simulation(freq=frequency, current=[0, I0], phi=[0, 0], skin_mesh_factor=accuracy)


if geo.component_type == "inductor":
    geo.update_conductors(n_turns=[20],
                          conductor_type=["solid"],
                          conductor_radix=[0.0015])
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
