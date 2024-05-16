"""Contains information about the file structure."""
# Python standard libraries
import os

import numpy as np
from typing import List

# Local libraries
from femmt.enumerations import ConductorType
from femmt.model import Conductor
from typing import Optional


class FileData:
    """Contains paths to every folder and file needed in femmt."""

    def __init__(self, working_directory: str, electro_magnetic_folder_path: str = None, strands_coefficients_folder_path: str = None):
        if working_directory is not None:
            self.update_paths(working_directory, electro_magnetic_folder_path, strands_coefficients_folder_path)
        self.onelab_folder_path: Optional[str] = None

    @staticmethod
    def create_folders(*args) -> None:
        """Create folders for every given folder path (if it does not exist)."""
        for folder in list(args):
            if not os.path.exists(folder):
                os.mkdir(folder)

    def clear_previous_simulation_results(self):
        """
        Clear all simulation results from previous simulations.

        Therefore, the result-folder structure as well as some temporary files
        (Parameter.pro, core_materials_temp.pro) are cleaned up.
        """
        self.clean_folder_structure(self.results_folder_path)
        if os.path.exists(os.path.join(self.electro_magnetic_folder_path, "core_materials_temp.pro")):
            os.remove(os.path.join(self.electro_magnetic_folder_path, "core_materials_temp.pro"))
        if os.path.exists(os.path.join(self.electro_magnetic_folder_path, "Parameter.pro")):
            os.remove(os.path.join(self.electro_magnetic_folder_path, "Parameter.pro"))

    @staticmethod
    def clean_folder_structure(folder_path: str):
        """Clean all files from a folder structure. The folder structure remains intact."""
        try:
            for root, _, files in os.walk(folder_path):
                for file in files:
                    file_path = os.path.join(root, file)
                    os.remove(file_path)
                    # print(f"remove {file_path}")
            # print("All simulation results from previous simulations have been deleted successfully.")
        except OSError:
            print("Error occurred while deleting files and subdirectories.")

    def update_paths(self, working_directory: str, electro_magnetic_folder_path: str = None, strands_coefficients_folder_path: str = None) -> None:
        """Set the local path based on the given working directory.

        :param working_directory: working directory folder path
        :type working_directory: str
        :param electro_magnetic_folder_path: folder path to electro magnetic simulation folders
        :type electro_magnetic_folder_path: str
        :param strands_coefficients_folder_path: folder path to strand coefficients
        :type strands_coefficients_folder_path: str
        """
        # Setup folder paths 
        self.working_directory = working_directory
        self.femmt_folder_path = os.path.dirname(__file__)
        self.mesh_folder_path = os.path.join(self.working_directory, "mesh")
        if electro_magnetic_folder_path:
            self.electro_magnetic_folder_path = electro_magnetic_folder_path
        else:
            self.electro_magnetic_folder_path = os.path.join(self.femmt_folder_path, "electro_magnetic")
        self.results_folder_path = os.path.join(self.working_directory, "results")
        self.e_m_values_folder_path = os.path.join(self.results_folder_path, "values")
        self.e_m_fields_folder_path = os.path.join(self.results_folder_path, "fields")
        self.e_m_circuit_folder_path = os.path.join(self.results_folder_path, "circuit")
        if strands_coefficients_folder_path:
            self.e_m_strands_coefficients_folder_path = strands_coefficients_folder_path
        else:
            self.e_m_strands_coefficients_folder_path = os.path.join(self.electro_magnetic_folder_path, "Strands_Coefficients")
        self.femm_folder_path = os.path.join(self.working_directory, "femm")
        self.reluctance_model_folder_path = os.path.join(self.working_directory, "reluctance_model")
        self.thermal_results_folder_path = os.path.join(self.results_folder_path, "thermal")

        # Setup file paths
        self.e_m_results_log_path = os.path.join(self.results_folder_path, "log_electro_magnetic.json")
        self.coordinates_description_log_path = os.path.join(self.results_folder_path, "log_coordinates_description.json")
        self.material_log_path = os.path.join(self.results_folder_path, "log_material.json")
        self.femm_results_log_path = os.path.join(self.femm_folder_path, "result_log_femm.json")
        self.config_path = os.path.join(self.femmt_folder_path, "config.json")
        self.e_m_mesh_file = os.path.join(self.mesh_folder_path, "electro_magnetic.msh")
        self.model_geo_file = os.path.join(self.mesh_folder_path, "model.geo_unrolled")
        self.hybrid_mesh_file = os.path.join(self.mesh_folder_path, "hybrid.msh")
        self.hybrid_color_mesh_file = os.path.join(self.mesh_folder_path, "hybrid_color.msh")
        self.hybrid_color_visualize_file = os.path.join(self.mesh_folder_path, "hybrid_color.png")
        self.thermal_mesh_file = os.path.join(self.mesh_folder_path, "thermal.msh")
        self.results_em_simulation = os.path.join(self.mesh_folder_path, "results.png")
        self.gmsh_log = os.path.join(self.results_folder_path, "log_gmsh.txt")
        self.getdp_log = os.path.join(self.results_folder_path, "log_getdp.txt")
        self.femmt_log = os.path.join(self.results_folder_path, "log_femmt.txt")

        # Create necessary folders
        self.create_folders(self.femmt_folder_path, self.mesh_folder_path, self.electro_magnetic_folder_path, 
                            self.results_folder_path, self.e_m_values_folder_path, self.e_m_fields_folder_path,
                            self.e_m_circuit_folder_path, self.e_m_strands_coefficients_folder_path)


class MeshData:
    """Contains data which is needed for the mesh generation. Is updated by high_level_geo_gen."""

    padding: float           # > 1
    skin_mesh_factor: float
    c_core: float
    c_window: float
    c_conductor = List[float]
    c_center_conductor = List[float]
    c_air_gaps: float
    
    center_factor: float

    mu0: float
    core_w: float
    window_w: float
    windings: List["Conductor"]  # This is written as string because it is a forward import

    def __init__(self, mesh_accuracy_core: float,
                 mesh_accuracy_window: float,
                 mesh_accuracy_conductor: float,
                 mesh_accuracy_air_gaps: float,
                 padding: float,
                 mu0: float):
        self.mesh_accuracy_core = mesh_accuracy_core
        self.mesh_accuracy_window = mesh_accuracy_window
        self.mesh_accuracy_conductor = mesh_accuracy_conductor
        self.mesh_accuracy_air_gaps = mesh_accuracy_air_gaps
        self.padding = padding
        self.mu0 = mu0

        # TODO This value should be set from user/outside?
        # The value 4 has good impact on the runtime and does not affect the simulation results too much.
        self.center_factor = 4

    def update_spatial_data(self, core_w: float, window_w: float, windings: List["Conductor"]):
        """Update geometry data of the core of the magnetic component."""
        self.core_w = core_w
        self.window_w = window_w
        self.windings = windings

        # Empty lists
        self.c_conductor = [None] * len(windings)
        self.c_center_conductor = [None] * len(windings)

        self.c_core = core_w / 10 * self.mesh_accuracy_core
        self.c_window = window_w / 30 * self.mesh_accuracy_window
        self.c_air_gaps = window_w / 20 * self.mesh_accuracy_air_gaps

    def update_data(self, frequency: float, skin_mesh_factor: float) -> None:
        """Update the mesh data according to the given frequency and skin_mesh_factor.

        :param frequency: Frequency of the model (updates skin depth which affects the mesh)
        :type frequency: float
        :param skin_mesh_factor: Factor for skin mesh
        :type skin_mesh_factor: float
        """
        self.skin_mesh_factor = skin_mesh_factor

        # Update Skin Depth (needed for meshing)
        if frequency is not None:
            if frequency == 0:
                self.delta = 1e9
            else:
                self.delta = np.sqrt(2 / (2 * frequency * np.pi * self.windings[0].cond_sigma * self.mu0))
            for i in range(len(self.windings)):
                if self.windings[i].conductor_type == ConductorType.RoundSolid:
                    self.c_conductor[i] = min([self.delta * self.skin_mesh_factor, self.windings[i].conductor_radius / 4 * self.mesh_accuracy_conductor])
                    # * self.mesh.skin_mesh_factor])
                    self.c_center_conductor[i] = self.windings[i].conductor_radius / 4 * self.mesh_accuracy_conductor  # * self.mesh.skin_mesh_factor
                elif self.windings[i].conductor_type == ConductorType.RoundLitz:
                    self.c_conductor[i] = self.windings[i].conductor_radius / 4 * self.mesh_accuracy_conductor
                    self.c_center_conductor[i] = self.windings[i].conductor_radius / 4 * self.mesh_accuracy_conductor
                else:
                    self.c_conductor[i] = self.windings[i].thickness / 4 * self.mesh_accuracy_conductor  # TODO: dynamic implementation
                    self.c_center_conductor[i] = self.center_factor * self.windings[i].thickness / 4 * self.mesh_accuracy_conductor
