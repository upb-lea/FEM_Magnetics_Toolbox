from femmt import MagneticComponent

# -- Control Parameters --
inductance_calculation = 0
accuracy = 0.5

# -- Definitions --
frequency = 100000
N1 = 6
N2 = 14
I0 = 10

# -- Component Selection --
geo = MagneticComponent(component_type="transformer")

# -- Core --
geo.update_core(core_type="EI", window_h=0.03, window_w=0.01, core_w=0.02)

# -- Air Gaps --
geo.update_air_gaps(method="center", n_air_gaps=1, air_gap_h=[0.002])

# -- Conductors --
geo.update_conductors(n_turns=[N1, N2],
                      conductor_type=["solid", "solid"],
                      conductor_radix=[0.0016, 0.0015],
                      scheme=["square", "hexa"],
                      cond_cond_isolation=[0.0002, 0.0002, 0.0008],
                      core_cond_isolation=[0.0004, 0.0004]
                      )

# -- Simulation --
if inductance_calculation:
    geo.get_inductances(I0=I0, op_frequency=frequency, mesh_accuracy=accuracy)
else:
    # Attention: ideal coupling assumed (only for demonstration)
    # geo.femm_reference(freq=frequency, current=[I0, N1/N2*I0], sign=[1, -1], sigma=58, non_visualize=0)
    geo.single_simulation(freq=frequency, current=[I0, N1/N2*I0], phi=[0, 180], skin_mesh_factor=accuracy)
