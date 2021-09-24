# minimal example for github clone
from FEMMT import MagneticComponent

# Create Object
geo = MagneticComponent(component_type="inductor")
#geo = MagneticComponent(component_type="transformer")

# Update Geometry
geo.update_core(core_type="EI", window_h=0.03)

geo.update_air_gaps(method="center", n_air_gaps=1, air_gap_h=[0.001])

geo.update_conductors(n_turns=[[14]], conductor_type=["solid"], conductor_radii=[0.0015],
                      winding=["primary"], scheme=["square"],
                      core_cond_isolation=[0.0005], cond_cond_isolation=[0.0001])

#geo.update_conductors(n_turns=[[6, 0], [0, 6]], conductor_type=["solid", "solid"],
#                      conductor_radii=[0.0015, 0.0015], winding=["interleaved", "interleaved"],
#                      scheme=["horizontal", "horizontal"],
#                      cond_cond_isolation=[0.0001, 0.0001, 0.0003], core_cond_isolation=[0.0005])

# Perform a single simulation
geo.single_simulation(freq=1000000, current=[10])
#geo.single_simulation(freq=1000000, current=[10, 10])