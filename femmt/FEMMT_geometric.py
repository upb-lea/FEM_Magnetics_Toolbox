from FEMMT import MagneticComponent
import numpy as np

# Folders must exits - currently not created automatically
tag='/log/litz_eva/layers1'


layers = [4, 5, 6, 9, 10, 11, 12]  # [9]
radix =  [0.05e-3]  # [0.04e-3, 0.05e-3, 0.06e-3]  #
fillfactors = [0.5]  # [0.4, 0.5, 0.6, 0.7]


# Create object instance
geo = MagneticComponent(conductor_type="litz")

# Perform sweeps
for layer in layers:
    for radius in radix:
        for ff in fillfactors:
            geo.litz_loss_comparison(FF=ff, n_layers=layer, strand_radius=radius, sim_choice='both', nametag=tag)

# -- Load Results and plot'em --
geo.load_litz_loss_logs(tag=tag)