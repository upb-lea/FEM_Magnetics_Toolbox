# Python standard libraries
import os
import shutil

import numpy as np
from typing import List

import femmt
# Local libraries
from femmt.enumerations import ConductorType
from femmt.model import Conductor

class FileData:
    """Contains paths to every folder and file needed in femmt.
    """
    def __init__(self, working_directory: str, electro_magnetic_folder_path: str = None, strands_coefficients_folder_path: str = None):
        if working_directory is not None:
            self.update_paths(working_directory, electro_magnetic_folder_path, strands_coefficients_folder_path)

    @staticmethod
    def create_folders(*args) -> None:
        """
        Creates folder for every given folder path (if it does not exist).
        """
        for folder in list(args):
            if not os.path.exists(folder):
                os.mkdir(folder)

    def clear_previous_simulation_results(self):
        self.clean_folder_structure(self.results_folder_path)

    @staticmethod
    def clean_folder_structure(folder_path: str):
        """
        Cleans all files from a folder structure. The folder structure remains intact.
        """
        try:
            for root, dirs, files in os.walk(folder_path):
                for file in files:
                    file_path = os.path.join(root, file)
                    os.remove(file_path)
                    # print(f"remove {file_path}")
            # print("All simulation results from previous simulations have been deleted successfully.")
        except OSError:
            print("Error occurred while deleting files and subdirectories.")

    def update_paths(self, working_directory: str, electro_magnetic_folder_path: str = None, strands_coefficients_folder_path: str = None) -> None:
        """Sets the local path based on the given working directory

        :param working_directory: working directory folder path
        :type working_directory: str
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
    """Contains data which is needed for the mesh generation.
    Is updated by high_level_geo_gen.
    """

    global_accuracy: float  # Parameter for mesh-accuracy
    padding: float           # > 1
    skin_mesh_factor: float
    c_core : float
    c_window: float
    c_conductor = List[float]
    c_center_conductor = List[float]

    mu0: float
    core_w: float
    window_w: float
    windings: List["Conductor"] # This is written as string because it is a forward import

    def __init__(self, global_accuracy: float, padding: float, mu0: float, core_w: float, window_w: float, windings: List["Conductor"]):
        self.global_accuracy = global_accuracy
        self.padding = padding
        self.mu0 = mu0
        self.core_w = core_w
        self.window_w = window_w
        self.windings = windings

        # Empty lists
        self.c_conductor = [None] * len(windings)
        self.c_center_conductor = [None] * len(windings)

    def update_data(self, frequency: float, skin_mesh_factor: float) -> None:
        """Updates the mesh data according to the given frequency and skin_mesh_factor.

        :param frequency: Frequency of the model (updates skin depth which affects the mesh)
        :type frequency: float
        :param skin_mesh_factor: Factor for skin mesh
        :type skin_mesh_factor: float
        """

        # Mesh-Parameters must be updated depending on geometry size
        self.c_core = self.core_w / 10. * self.global_accuracy
        self.c_window = self.window_w / 30 * self.global_accuracy
        self.skin_mesh_factor = skin_mesh_factor

        # Update Skin Depth (needed for meshing)
        if frequency is not None:
            if frequency == 0:
                self.delta = 1e9
            else:
                self.delta = np.sqrt(2 / (2 * frequency * np.pi * self.windings[0].cond_sigma * self.mu0))
            for i in range(len(self.windings)):
                if self.windings[i].conductor_type == ConductorType.RoundSolid:
                    self.c_conductor[i] = min([self.delta * self.skin_mesh_factor, self.windings[i].conductor_radius / 4 * self.global_accuracy]) #* self.mesh.skin_mesh_factor])
                    self.c_center_conductor[i] = self.windings[i].conductor_radius / 4 * self.global_accuracy  # * self.mesh.skin_mesh_factor
                elif self.windings[i].conductor_type == ConductorType.RoundLitz:
                    self.c_conductor[i] = self.windings[i].conductor_radius / 4 * self.global_accuracy
                    self.c_center_conductor[i] = self.windings[i].conductor_radius / 4 * self.global_accuracy
                else:
                    self.c_conductor[i] = 0.0001  # TODO: dynamic implementation

