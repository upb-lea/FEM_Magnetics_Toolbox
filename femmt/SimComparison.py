import pandas as pd
from FEMMT import MagneticComponent
import numpy as np

# Create Object
geo = MagneticComponent()

# Generate Mesh
geo.mesh()
# Iterate on frequency
frequencies = np.linspace(0, 250000, 6)
currents = [2, 2, 2, 2, 2, 2]
geo.excitation_sweep(currents=currents, frequencies=frequencies)

# Reference simulation with FEMM
losses = []
for f in frequencies:
    geo.femm_reference(freq=f, current=2, sigma=58)
    losses.append(geo.tot_loss_femm.real)


"""
Pv_sim = np.array([0.12,  2.16, 7.83, 16.66, 28.82, 43.93])
freq_sim = np.array([0, 50, 100, 150, 200, 250])




fig, ax = plt.subplots()
ax.axis
ax.plot(freq, Pv, label='litz (6 layer, 0.1 mm), 3A')
ax.set_xlabel('$f$ in Hz')
ax.set_ylabel('$P_v$ in W')
ax.legend()
plt.show()
"""