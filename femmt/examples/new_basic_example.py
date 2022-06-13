import femmt
import femmt as fmt
import numpy as np
from femmt_windings import *

geo = fmt.MagneticComponent(component_type="transformer")

geo.core.update(window_h=0.0295, window_w=0.012, core_w=0.015,
                mu_rel=3100, phi_mu_deg=12,
                sigma=0.6)

geo.air_gaps.update(method="percent", n_air_gaps=1, air_gap_h=[0.0005],
                    air_gap_position=[50], position_tag=[0])

windings = Windings()
windings.add_core_isolations(0.001, 0.001, 0.002, 0.001)
windings.add_winding_isolations(0.0002, 0.0002, 0.0005)

winding1 = Winding(10, 0, Conductivities.Copper, WindingTypes.Primary, WindingSchemes.Square)
winding1.set_solid_conductor(0.0011)

winding2 = Winding(0, 10, Conductivities.Copper, WindingTypes.Secondary, WindingSchemes.Square)
winding2.set_litz_conductor(None, 600, 35.5e-6, 0.6)

windings.add_windings([winding1, winding2])


print(windings.get_all_parameters())
geo.update_conductors(**windings.get_all_parameters())

"""
geo.update_conductors(n_turns=[[10, 0], [0, 10]], conductor_type=["solid", "litz"],
                    litz_para_type=['implicit_litz_radius', 'implicit_litz_radius'],
                    ff=[None, 0.6], strands_numbers=[None, 600], strand_radii=[70e-6, 35.5e-6],
                    conductor_radii=[0.0011, None],
                    winding=["primary", "secondary"], scheme=["square", "square"],
                    core_cond_isolation=[0.001, 0.001, 0.002, 0.001], cond_cond_isolation=[0.0002, 0.0002, 0.0005],
                    conductivity_sigma=["copper", "copper"])
"""

geo.create_model(freq=250000, visualize_before=True)
geo.single_simulation(freq=250000, current=[4.14723021, 14.58960019], phi_deg=[- 1.66257715/np.pi*180, 170])