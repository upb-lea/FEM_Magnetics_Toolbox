from FEMMT import MagneticComponent

"""
# Transformer
geo = MagneticComponent(component_type="transformer")

geo.update_conductors(n_turns=[6, 8],
                      conductor_type=["solid", "solid"],
                      conductor_radix=[0.0015, 0.001])

geo.update_conductors(n_turns=[10, 10],
                      conductor_type=["solid", "litz"],
                      conductor_radix=[0.0015, 0.0015],
                      layer_numbers=[10, 10],
                      strand_radix=[0.05e-3, 0.05e-3])
# 0.0011879
# 0.0015
"""

# Inductor
geo = MagneticComponent(component_type="inductor")

geo.update_core(core_type="EI",
                window_h=0.03)

geo.update_air_gaps(n_air_gaps=0)

geo.update_conductors(n_turns=[6],
                      conductor_type=["solid"],
                      conductor_radix=[0.0015])

"""
geo.update_conductors(n_turns=[6],
                      conductor_type=["litz"],
                      conductor_radix=[0.0015],
                      layer_numbers=[10],
                      strand_radix=[0.05e-3])
"""

geo.single_simulation(freq=100000, current=[10])
