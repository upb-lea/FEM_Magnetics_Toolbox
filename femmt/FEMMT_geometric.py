from FEMMT import MagneticComponent
import numpy as np

# Create Object
geo = MagneticComponent()

"""
# Perform a single simulation
geo.single_simulation(freq=100000,current=1)
#geo.pre_simulation()
"""

# ===============
# Generate Mesh
geo.mesh()
# Iterate on frequency
frequencies = np.linspace(0, 250000, 6)
currents = [10, 2, 1, 0.5, 0.2, 0.1]
geo.excitation_sweep(currents=currents, frequencies=frequencies)
