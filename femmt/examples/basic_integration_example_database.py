import materialdatabase as mdb
import matplotlib.pyplot as plt
database = mdb.MaterialDatabase()

mdb.compare_core_loss_flux_density_data(material_list=["N95", "N87"], temperature=25)
mdb.compare_core_loss_temperature(material_list=["N95", "N87"], flux=200e-3)