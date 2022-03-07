# from BH_extractor import *
from mu_extractor import *


# Relative Permeability
directory = f"C:/Users/tillp/sciebo/Exchange_FEMMT/05_Materials/data/Loss_Data_PHD_Keuck/mu_r_Plot"
temperatures = [30, 60, 80, 100]
frequencies = [300000]
# load_and_extract(parameter="mu_r", read_directory=directory, write_directory="materials/N95/mu_r/", temperatures=temperatures, frequencies=frequencies, material="N95", do_plot=True)

# Permeability loss angle
directory = f"C:/Users/tillp/sciebo/Exchange_FEMMT/05_Materials/data/Loss_Data_PHD_Keuck/mu_phi_Plot"
temperatures = [30, 60, 80, 100]
frequencies = [100000, 200000, 300000, 400000]
# load_and_extract(parameter="mu_phi", read_directory=directory, write_directory="materials/N95/mu_phi/", temperatures=temperatures, frequencies=frequencies, material="N95", do_plot=True)


# create_arithmetic_form
temperatures = [30, 60, 80, 100]
frequencies = [100000, 200000, 300000, 400000]
create_arithmetic_form(temperatures=temperatures, frequencies=frequencies, material="N95", mu_r_frequency=300000)
