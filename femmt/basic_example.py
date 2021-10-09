# minimal example for github clone
from femmt import MagneticComponent

# Create Object
#geo = MagneticComponent(component_type="inductor")
geo = MagneticComponent(component_type="transformer")

# Update Geometry
geo.core.update(type="EI", window_h=0.03)

# geo.air_gaps.update(n_air_gaps=0)
geo.air_gaps.update(method="center", n_air_gaps=1, air_gap_h=[0.002], air_gap_position=[0])

#geo.update_conductors(n_turns=[[14]], conductor_type=["solid"], conductor_radii=[0.0015],
#                      winding=["primary"], scheme=["square"],
#                      core_cond_isolation=[0.0005], cond_cond_isolation=[0.0001])

geo.update_conductors(n_turns=[[3, 0], [0, 5]], conductor_type=["litz", "litz"],
                      litz_para_type=['implicit_litz_radius', 'implicit_litz_radius'],
                      ff=[0.5, 0.5], strands_numbers=[400, 200], strand_radii=[70e-6, 70e-6],
                      winding=["interleaved", "interleaved"], scheme=["horizontal", "horizontal"],
                      core_cond_isolation=[0.0005, 0.0005], cond_cond_isolation=[0.0001, 0.0001, 0.0001])


#geo.update_conductors(n_turns=[[6, 0], [0, 6]], conductor_type=["solid", "solid"],
#                      conductor_radii=[0.0015, 0.0015], winding=["interleaved", "interleaved"],
#                      scheme=["horizontal", "horizontal"],
#                      cond_cond_isolation=[0.0001, 0.0001, 0.0003], core_cond_isolation=[0.0005])

# Perform a single simulation
#geo.single_simulation(freq=1000000, current=[10])
#geo.single_simulation(freq=1000000, current=[10, 10])
#geo.get_inductances(I0=10, op_frequency=100000, skin_mesh_factor=0.5)

geo.femm_reference(freq=100000, current=[1, 2], sigma_cu=58, sign=[1, -1], non_visualize=0)
