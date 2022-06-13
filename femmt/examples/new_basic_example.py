import femmt
import femmt as fmt
import numpy as np
from femmt_windings import *
from femmt_air_gaps import *

geo = fmt.MagneticComponent(component_type="transformer")

geo.core.update(window_h=0.0295, window_w=0.012, core_w=0.015,
                mu_rel=3100, phi_mu_deg=12,
                sigma=0.6)

"""
geo.air_gaps.update(method="percent", n_air_gaps=1, air_gap_h=[0.0005],
                    air_gap_position=[50], position_tag=[0])
"""

# Add air gaps
air_gaps = AirGaps(AirGapMethods.Percent)
air_gaps.add_air_gap(AirGapLegPositions.CenterLeg, 50, 0.0005)
geo.air_gaps.update(**air_gaps.get_all_parameters())

# Add windings
windings = Windings()
windings.add_core_isolations(0.001, 0.001, 0.002, 0.001)
windings.add_winding_isolations(0.0002, 0.0002, 0.0005)

winding1 = Winding(10, 0, Conductivities.Copper, WindingTypes.Primary, WindingSchemes.Square)
winding1.set_solid_conductor(0.0011)

winding2 = Winding(0, 10, Conductivities.Copper, WindingTypes.Secondary, WindingSchemes.Square)
winding2.set_litz_conductor(None, 600, 35.5e-6, 0.6)

windings.add_windings([winding1, winding2])
geo.update_conductors(**windings.get_all_parameters())

# Create model
geo.create_model(freq=250000, visualize_before=True)
geo.single_simulation(freq=250000, current=[4.14723021, 14.58960019], phi_deg=[- 1.66257715/np.pi*180, 170])