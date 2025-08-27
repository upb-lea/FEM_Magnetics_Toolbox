"""Common data classes uses by multiple optimization classes and methods."""

# python libraries
import dataclasses

# own libraries
from materialdatabase.meta.data_enums import DataSource

@dataclasses.dataclass
class MaterialDataSources:
    """Data sources for the FEM simulation."""

    permeability_datasource: DataSource
    permittivity_datasource: DataSource

@dataclasses.dataclass
class WorkingDirectories:
    """Working directories for an integrated transformer optimization."""

    fem_working_directory: str
    reluctance_model_results_directory: str
    fem_simulation_results_directory: str
    fem_simulation_filtered_results_directory: str
    fem_thermal_simulation_results_directory: str
    fem_thermal_filtered_simulation_results_directory: str
