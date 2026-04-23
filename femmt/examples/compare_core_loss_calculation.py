"""Example how to perform an optimization using optuna and FEMMT."""
# python libraries

# 3rd party libraries
import numpy as np
import pandas as pd

# femmt libraries
import femmt as fmt
from femmt.optimization import InductorOptimization as io
from femmt import ImportedComplexCoreMaterial
import magnethub as mh
from materialdatabase.meta.data_enums import Material, DataSource
from materialdatabase import Data
i_offset = 0

# Exact a single current waveform to optimize the inductor
time_vec: np.ndarray = np.array([0.0, 3.265248131976911e-07, 2.5e-06, 2.8265248131976912e-06, 5e-06])
current_vec: np.ndarray = i_offset + np.array([-0.9996115022426437, 4.975792579275104, 0.9996115022426446, -4.975792579275103, -0.9996115022426437])

i_1: list[np.ndarray] = [time_vec, current_vec]

# Initialize the material data source
material_data_sources = fmt.MaterialDataSources(
    permeability_datasource=DataSource("MagNet"),
    permittivity_datasource=DataSource("LEA_MTB"),
)

material_db = Data()
material = Material("N49")

if __name__ == '__main__':

    # target_total_reluctance 10/20=0.5 -> 0.5=Bsat
    target_total_reluctance = 10

    # Prepare equidistant flux vector for the magnet loss simulation
    time_interp = np.linspace(time_vec[0], time_vec[-1], 1024)
    flux_density_vec = np.interp(time_interp, time_vec, current_vec) / target_total_reluctance
    fundamental_frequency = 1/time_vec[-1]

    # Debug
    flux_avg_density = 0

    # Variable declaration
    temperature = 80
    # instantiate material-specific model
    magnet_material_model: mh.loss.LossModel = mh.loss.LossModel(material=material.name, team="paderborn")
    # Initial magnetization curve
    initial_mag_curve = material_db.get_initial_magnetization_curve(material, 500, temperature)
    imported_complex_material = ImportedComplexCoreMaterial(material, temperature,
                                                            material_data_sources.permeability_datasource,
                                                            material_data_sources.permittivity_datasource)

    result_list = []

    for count in range(100):
        factor = count/100
        # Calculate the flux without DC-offset
        flux_density = flux_density_vec * factor
        flux_peak_density = (flux_density.max() - flux_density.min())/2

        # p_loss calculation
        # get power loss in W/m³ and estimated H wave in A/m
        p_density_cmp, _ = magnet_material_model(flux_density, fundamental_frequency, temperature)
        p_density = io.ReluctanceModel.get_power_density(
            flux_avg_density, flux_peak_density, imported_complex_material,
            initial_mag_curve, float(fundamental_frequency), temperature)

        p_density_factor = p_density_cmp / p_density * 170 / 21
        act_k = p_density / p_density_cmp

        result_list.append((flux_peak_density, p_density, p_density_cmp, act_k))

    result: pd.DataFrame = pd.DataFrame(result_list, columns=["flux_peak_density", "p_density", "p_density_magnet", "act_k"])
    result.to_csv("/home/andreas/Workspace/Projekt/Test/Testdata/p_loss_vgl.csv")
