# Python libraries
import inspect
import os
import json
import itertools
import shutil
from dataclasses import dataclass
from typing import List

# 3rd party library import
import numpy as np
import materialdatabase as mdb
from scipy import optimize

# femmt import
import femmt.functions as ff
import femmt.functions_reluctance as fr
import femmt.optimization.functions_optimization as fof
from femmt.optimization.integrated_transformer_dtos import *



def config_file_dto_to_dict(config_file_dto: InputConfig):

    config_file_dict = {
        "l_s_target": config_file_dto.l_s_target,
        "l_h_target": config_file_dto.l_h_target,
        "n_target": config_file_dto.n_target,
        "time_current_1_vec": config_file_dto.time_current_1_vec.tolist(),
        "time_current_2_vec": config_file_dto.time_current_2_vec.tolist(),
        "material_list": config_file_dto.material_list,
        "core_inner_diameter_min_max_count_list": config_file_dto.core_inner_diameter_min_max_count_list,
        "window_w_min_max_count_list": config_file_dto.window_w_min_max_count_list,
        "window_h_top_min_max_count_list": config_file_dto.window_h_top_min_max_count_list,
        "window_h_bot_min_max_count_list": config_file_dto.window_h_bot_min_max_count_list,
        "factor_max_flux_density": config_file_dto.factor_max_flux_density,
        "n_p_top_min_max_list": config_file_dto.n_p_top_min_max_list,
        "n_p_bot_min_max_list": config_file_dto.n_p_bot_min_max_list,
        "n_s_top_min_max_list": config_file_dto.n_s_top_min_max_list,
        "n_s_bot_min_max_list": config_file_dto.n_s_bot_min_max_list
    }
    return config_file_dict






def result_file_dto_to_dict(reluctance_model_result_dto):
    result_file_dict = {
        "case": reluctance_model_result_dto.case,
        "air_gap_top": reluctance_model_result_dto.air_gap_top,
        "air_gap_bot": reluctance_model_result_dto.air_gap_bot,
        "air_gap_middle": reluctance_model_result_dto.air_gap_middle,
        "n_p_top": reluctance_model_result_dto.n_p_top,
        "n_p_bot": reluctance_model_result_dto.n_p_bot,
        "n_s_top": reluctance_model_result_dto.n_s_top,
        "n_s_bot": reluctance_model_result_dto.n_s_bot,
        "window_h_top": reluctance_model_result_dto.window_h_top,
        "window_h_bot": reluctance_model_result_dto.window_h_bot,
        "window_w": reluctance_model_result_dto.window_w,
        "mu_r_abs": reluctance_model_result_dto.mu_r_abs,
        "core_inner_diameter": reluctance_model_result_dto.core_inner_diameter,
        "flux_top_max": reluctance_model_result_dto.flux_top_max,
        "flux_bot_max": reluctance_model_result_dto.flux_bot_max,
        "flux_stray_max": reluctance_model_result_dto.flux_stray_max,
        "p_hyst": reluctance_model_result_dto.p_hyst,
        "core_2daxi_total_volume": reluctance_model_result_dto.core_2daxi_total_volume,
    }
    return result_file_dict

def result_file_dict_to_dto(result_file_dict):
    result_file_dto = ResultFile(
        case = result_file_dict["case"],
        air_gap_top = result_file_dict["air_gap_top"],
        air_gap_bot= result_file_dict["air_gap_bot"],
        air_gap_middle=result_file_dict["air_gap_middle"],
        n_p_top=result_file_dict["n_p_top"],
        n_p_bot=result_file_dict["n_p_bot"],
        n_s_top=result_file_dict["n_s_top"],
        n_s_bot=result_file_dict["n_s_bot"],
        window_h_top=result_file_dict["window_h_top"],
        window_h_bot=result_file_dict["window_h_bot"],
        window_w=result_file_dict["window_w"],
        mu_r_abs=result_file_dict["mu_r_abs"],
        core_inner_diameter=result_file_dict["core_inner_diameter"],
        flux_top_max=result_file_dict["flux_top_max"],
        flux_bot_max=result_file_dict["flux_bot_max"],
        flux_stray_max=result_file_dict["flux_stray_max"],
        p_hyst = result_file_dict["p_hyst"],
        core_2daxi_total_volume=result_file_dict["core_2daxi_total_volume"]
    )
    return result_file_dto



class IntegratedTransformerOptimization:


    def __init__(self, working_directory: str = None, silent: bool = False):

        # Variable to set silent mode
        ff.set_silent_status(silent)

        self.set_up_folder_structure(working_directory)
        self.material_db = mdb.MaterialDatabase(is_silent=silent)

        ff.femmt_print(f"\n"
              f"Start optimization \n"
              f"--- --- --- ---")

    def set_up_folder_structure(self, working_directory):
        if working_directory is None:
            caller_filename = inspect.stack()[1].filename
            working_directory = os.path.join(os.path.dirname(caller_filename), "integrated_transformer_optimization")

        # generate new and empty working directory
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # set up folders for optimization
        self.optimization_working_directory = working_directory
        self.femmt_working_directory = os.path.join(self.optimization_working_directory, "femmt_simulation")
        self.integrated_transformer_reluctance_model_results_directory = os.path.join(self.optimization_working_directory, "reluctance_model_results")
        self.integrated_transformer_fem_simulations_results_directory = os.path.join(self.optimization_working_directory, "fem_simulation_results")
        self.integrated_transformer_optimization_input_parameters_file = os.path.join(self.optimization_working_directory, "optimization_input_parameters.json")

        self.create_folders(self.optimization_working_directory, self.integrated_transformer_reluctance_model_results_directory, self.integrated_transformer_fem_simulations_results_directory, self.femmt_working_directory)



    @staticmethod
    def create_folders(*args) -> None:
        """Creates folder for every given folder path (if it does not exist).
        """
        for folder in list(args):
            if not os.path.exists(folder):
                os.mkdir(folder)

    @staticmethod
    def calculate_sweep_tensors(input_parameters_dto: InputConfig):

        sweep_tensor = SweepTensor(
            # tensors: windings
            t1_n_p_top = np.arange(input_parameters_dto.n_p_top_min_max_list[0], input_parameters_dto.n_p_top_min_max_list[1] + 1),
            t1_n_p_bot = np.arange(input_parameters_dto.n_p_bot_min_max_list[0], input_parameters_dto.n_p_bot_min_max_list[1] + 1),
            t1_n_s_top = np.arange(input_parameters_dto.n_s_top_min_max_list[0], input_parameters_dto.n_s_top_min_max_list[1] + 1),
            t1_n_s_bot = np.arange(input_parameters_dto.n_s_bot_min_max_list[0], input_parameters_dto.n_s_bot_min_max_list[1] + 1),

            # tensors: outer core geometry and material
            t1_window_h_top = np.linspace(input_parameters_dto.window_h_top_min_max_count_list[0], input_parameters_dto.window_h_top_min_max_count_list[1],
                                          input_parameters_dto.window_h_top_min_max_count_list[2]),
            t1_window_h_bot = np.linspace(input_parameters_dto.window_h_bot_min_max_count_list[0], input_parameters_dto.window_h_bot_min_max_count_list[1],
                                          input_parameters_dto.window_h_bot_min_max_count_list[2]),
            t1_window_w = np.linspace(input_parameters_dto.window_w_min_max_count_list[0], input_parameters_dto.window_w_min_max_count_list[1],
                                      input_parameters_dto.window_w_min_max_count_list[2]),
            t1_core_material = input_parameters_dto.material_list,
            t1_core_inner_diameter = np.linspace(input_parameters_dto.core_inner_diameter_min_max_count_list[0],
                                                 input_parameters_dto.core_inner_diameter_min_max_count_list[1],
                                                 input_parameters_dto.core_inner_diameter_min_max_count_list[2]),


            time_current_1_vec = input_parameters_dto.time_current_1_vec,
            time_current_2_vec = input_parameters_dto.time_current_2_vec,

            l_s_target_value = input_parameters_dto.l_s_target,
            l_h_target_value = input_parameters_dto.l_h_target,
            n_target_value = input_parameters_dto.n_target,
            factor_max_flux_density = input_parameters_dto.factor_max_flux_density
        )
        return sweep_tensor

    def save_reluctance_model_result_list(self, config_file: InputConfig, result_file_list: List[ResultFile]):
        # save optimization input parameters
        config_dict = config_file_dto_to_dict(config_file)

        print(config_dict)

        with open(self.integrated_transformer_optimization_input_parameters_file, "w+", encoding='utf-8') as outfile:
            json.dump(config_dict, outfile, indent=2, ensure_ascii=False)

        # save reluctance parameters winning candidates
        for count, reluctance_model_result_dto in enumerate(result_file_list):
            file_name = os.path.join(self.integrated_transformer_reluctance_model_results_directory, f"valid_reluctance_design_{count}.json")
            reluctance_model_result_dict = result_file_dto_to_dict(reluctance_model_result_dto)
            with open(file_name, "w+", encoding='utf-8') as outfile:
                json.dump(reluctance_model_result_dict, outfile, indent=2, ensure_ascii=False)


    def load_optimization_input_parameters(self, filepath: str = None):
        if filepath is None:
            with open(self.integrated_transformer_optimization_input_parameters_file, "r") as fd:
                self.optimization_input_parameters = json.loads(fd.read())
        else:
            pass

        self.calculate_sweep_tensors(self.optimization_input_parameters)

    @staticmethod
    def t2_calculate_reluctance_matrix(t2_inductance_matrix, t2_winding_matrix, t2_winding_matrix_transpose):
        """
        Calculates the inductance matrix out of reluctance matrix and winding matrix.

        :param t2_inductance_matrix: matrix of transformer inductance
        :param t2_winding_matrix: matrix of transformer windings
        :return: reluctance matrix

        winding matrix e.g.
        N = [ [N_1a, N_2b], [N_1b, N_2b] ]

        inductance matrix e.g.
        L = [ [L_11, M], [M, L_22] ]

        returns reluctance matrix e.g.
        r = [ [], [] ]
        """

        # invert inductance matrix
        t2_inductance_matrix_invert = np.linalg.inv(t2_inductance_matrix)

        # Formular: L = N^T * R^-1 * N
        # Note: Be careful when trying to multiply the matrices in one single step. Some pre-tests failed.
        # The following commented example returns a different result as the code-version. The code-version is
        # verified with a 2x2 example.
        # So this line is not correct!
        # return np.einsum('...ij, ...jj, ...jk -> ...ik', t11_winding_matrix_transpose, t11_reluctance_matrix_invert,
        #                  t11_winding_matrix), t9_valid_design_mask

        # minimal example to understand the operation
        # matrix1 = np.array([[1, 2], [3, 4]])
        # matrix2 = np.array([[5, 6], [7, 8]])
        # matrix3 = np.array([[9, 10], [11, 12]])
        #
        # # reference code
        # normal_multiplication = np.matmul(np.matmul(matrix1, matrix2), matrix3)
        # print(f"{normal_multiplication = }")
        #
        # # This does not macht to the reference code!!!
        # einsum_multiplication = np.einsum('...ij, ...jj, ...jk -> ...ik', matrix1, matrix2, matrix3)
        # print(f"{einsum_multiplication = }")
        #
        # # two einsum multiplications: same result as reference code
        # einsum_multiplication_part_1 = np.einsum('...ij, ...jh -> ...ih', matrix1, matrix2)
        # einsum_multiplication_part_2 = np.einsum('...ij, ...jh -> ...ih', einsum_multiplication_part_1, matrix3)
        # print(f"{einsum_multiplication_part_2 = }")
        einsum_multiplication_part_1 = np.einsum('...ij, ...jh -> ...ih', t2_winding_matrix, t2_inductance_matrix_invert)
        einsum_multiplication_part_2 = np.einsum('...ij, ...jh -> ...ih', einsum_multiplication_part_1, t2_winding_matrix_transpose)

        return einsum_multiplication_part_2


    def load_reluctance_model_result_list(self) -> List[ResultFile]:
        ff.femmt_print(f"Read results from {self.integrated_transformer_reluctance_model_results_directory}")

        valid_design_list = []

        for file in os.listdir(self.integrated_transformer_reluctance_model_results_directory):
            if file.endswith(".json"):
                json_file_path = os.path.join(self.integrated_transformer_reluctance_model_results_directory, file)
                with open(json_file_path, "r") as fd:
                    loaded_data_dict = json.loads(fd.read())

                valid_design_list.append(result_file_dict_to_dto(loaded_data_dict))
        return valid_design_list

    @staticmethod
    def plot_reluctance_model_result_list(valid_design_list: List[ResultFile]):

        volume_list = []
        core_hyst_loss_list = []
        annotation_list = []

        for result in valid_design_list:
            volume_list.append(result.core_2daxi_total_volume)
            core_hyst_loss_list.append(result.p_hyst)
            annotation_list.append(result.case)


        fof.plot_2d(volume_list, core_hyst_loss_list, "volume in m³", "core hysteresis losses in W", "Pareto Diagram",
                    plot_color="red", annotations=annotation_list)


    def plot_pareto_result_list(self):

        fof.plot_2d(self.pareto_volume_list, self.pareto_core_hyst_list, "volume in m³", "core hysteresis losses in W", "Pareto Diagram", plot_color="red")


    @staticmethod
    def plot_filtered_pareto_result_list(filter_volume_list, filter_core_hyst_loss_list):
        fof.plot_2d(filter_volume_list, filter_core_hyst_loss_list, "volume in m³", "core hysteresis losses in W", "Pareto Diagram", plot_color="red")


    # Very slow for many datapo ints.  Fastest for many costs, most readable
    @staticmethod
    def is_pareto_efficient_dumb(costs):
        """
        Find the pareto-efficient points
        :param costs: An (n_points, n_costs) array
        :return: A (n_points, ) boolean array, indicating whether each point is Pareto efficient
        """
        is_efficient = np.ones(costs.shape[0], dtype=bool)
        for i, c in enumerate(costs):
            is_efficient[i] = np.all(np.any(costs[:i] > c, axis=1)) and np.all(np.any(costs[i + 1:] > c, axis=1))
        return is_efficient



    @staticmethod
    # Faster than is_pareto_efficient_simple, but less readable.
    def is_pareto_efficient(costs, return_mask=True):
        """
        Find the pareto-efficient points
        :param costs: An (n_points, n_costs) array
        :param return_mask: True to return a mask
        :return: An array of indices of pareto-efficient points.
            If return_mask is True, this will be an (n_points, ) boolean array
            Otherwise it will be a (n_efficient_points, ) integer array of indices.
        """
        is_efficient = np.arange(costs.shape[0])
        n_points = costs.shape[0]
        next_point_index = 0  # Next index in the is_efficient array to search for
        while next_point_index < len(costs):
            nondominated_point_mask = np.any(costs < costs[next_point_index], axis=1)
            nondominated_point_mask[next_point_index] = True
            is_efficient = is_efficient[nondominated_point_mask]  # Remove dominated points
            costs = costs[nondominated_point_mask]
            next_point_index = np.sum(nondominated_point_mask[:next_point_index]) + 1
        if return_mask:
            is_efficient_mask = np.zeros(n_points, dtype=bool)
            is_efficient_mask[is_efficient] = True
            return is_efficient_mask
        else:
            return is_efficient

    def pareto_front(self, x_vec, y_vec, dto_list):


        print(f"{len(x_vec) = }")
        tuple_vec = []

        for count_y, y in enumerate(y_vec):
            tuple_vec.append((x_vec[count_y], y_vec[count_y]))

        tuple_vec = np.array(tuple_vec)

        print(f"{tuple_vec = }")

        pareto_tuple_mask_vec = self.is_pareto_efficient(tuple_vec)

        x_pareto_vec = []
        y_pareto_vec = []
        pareto_dto_vec = []

        for count_mask, mask in enumerate(pareto_tuple_mask_vec):
            if mask:
                x_pareto_vec.append(x_vec[count_mask])
                y_pareto_vec.append(y_vec[count_mask])
                pareto_dto_vec.append(dto_list[count_mask])

        print(f"{len(x_pareto_vec) = }")

        return np.array(x_pareto_vec), np.array(y_pareto_vec), pareto_dto_vec



    def filter_reluctance_model_list(self, valid_design_list: List[ResultFile], factor_min_hyst_losses: float = 1.2) -> List[ResultFile]:
        # figure out pareto front
        #pareto_volume_list, pareto_core_hyst_list, pareto_dto_list = self.pareto_front(volume_list, core_hyst_loss_list, valid_design_list)

        volume_list = np.array([])
        min_hyst_losses_list = np.array([])
        pareto_dto_list = []
        # figure out pareto front with minimum loss-per-volume approach
        for result_dto in valid_design_list:
            # new volume
            if result_dto.core_2daxi_total_volume not in volume_list:
                volume_list = np.append(volume_list, result_dto.core_2daxi_total_volume)
                min_hyst_losses_list = np.append(min_hyst_losses_list, result_dto.p_hyst)
                pareto_dto_list.append(result_dto)
            # existing volume but lower losses
            elif np.shape(min_hyst_losses_list) == (1,) and np.isclose(result_dto.core_2daxi_total_volume, volume_list):
                if result_dto.p_hyst < min_hyst_losses_list:
                    min_hyst_losses_list = np.array([result_dto.p_hyst])
                    pareto_dto_list = [result_dto]
            else:
                reference_index = int(np.where(np.isclose(volume_list, result_dto.core_2daxi_total_volume))[0])
                if result_dto.p_hyst < min_hyst_losses_list[reference_index]:
                    min_hyst_losses_list[reference_index] = result_dto.p_hyst
                    pareto_dto_list[reference_index] = result_dto

        min_hyst_losses = min_hyst_losses_list[np.argmin(min_hyst_losses_list)]

        print(f"{volume_list = }")
        print(f"{min_hyst_losses_list = }")

        loss_offset = factor_min_hyst_losses * min_hyst_losses

        filtered_design_dto_list = []
        for count, dto in enumerate(valid_design_list):
            reference_index = np.where(np.isclose(dto.core_2daxi_total_volume, volume_list))[0]
            if dto.p_hyst < (min_hyst_losses_list[reference_index] + loss_offset):
                filtered_design_dto_list.append(dto)

        return filtered_design_dto_list

    @staticmethod
    def filter_air_gap(dto_list_to_filter: List[ResultFile], max_air_gap_length = 1e-9) -> List[ResultFile]:
        filtered_dtos = []

        for dto in dto_list_to_filter:
            if dto.air_gap_middle < max_air_gap_length:
                filtered_dtos.append(dto)
        return filtered_dtos


    def integrated_transformer_optimization(self, config_file: InputConfig):
        case_number = 0

        # 0. Empty folder
        if os.path.exists(self.optimization_working_directory):
            shutil.rmtree(self.optimization_working_directory)
            self.set_up_folder_structure(self.optimization_working_directory)


        sweep_dto = self.calculate_sweep_tensors(config_file)



        # 1. Extract fundamental frequency from current vectors
        time_extracted, current_extracted_1_vec = fr.time_vec_current_vec_from_time_current_vec(sweep_dto.time_current_1_vec)
        time_extracted, current_extracted_2_vec = fr.time_vec_current_vec_from_time_current_vec(sweep_dto.time_current_2_vec)
        fundamental_frequency = 1 / time_extracted[-1]
        print(f"{fundamental_frequency = }")

        # generate list of all parameter combinations
        t2_parameter_sweep = np.array(list(itertools.product(sweep_dto.t1_window_w, sweep_dto.t1_window_h_top,
                                                             sweep_dto.t1_window_h_bot, sweep_dto.t1_core_inner_diameter,
                                                             sweep_dto.t1_n_p_top, sweep_dto.t1_n_p_bot, sweep_dto.t1_n_s_top,
                                                             sweep_dto.t1_n_s_bot)))

        # report simulation progress
        number_of_simulations = len(t2_parameter_sweep) * len(sweep_dto.t1_core_material)
        ff.femmt_print(f"Simulation count: {number_of_simulations}")
        simulations_per_percent = int(number_of_simulations / 99)
        ff.femmt_print(f"{simulations_per_percent = }")
        simulation_progress_percent = 0

        valid_design_list = []

        # initialize parameters staying same form simulation
        t2_inductance_matrix = [
            [sweep_dto.l_s_target_value + sweep_dto.l_h_target_value, sweep_dto.l_h_target_value / sweep_dto.n_target_value],
            [sweep_dto.l_h_target_value / sweep_dto.n_target_value, sweep_dto.l_h_target_value / (sweep_dto.n_target_value ** 2)]]


        for count_core_material, material_name in enumerate(sweep_dto.t1_core_material):
            """
            outer core material loop loads material properties from material database
             * mu_r_abs
             * saturation_flux_density and calculates the dimensioning_flux_density from it
             * material vectors for mu_r_real and mu_r_imag depending on flux_density
             
            """

            mu_r_abs = self.material_db.get_material_property(material_name=material_name, property="initial_permeability")

            saturation_flux_density = self.material_db.get_saturation_flux_density(material_name=material_name)
            dimensioning_max_flux_density = saturation_flux_density * sweep_dto.factor_max_flux_density

            print(f"{saturation_flux_density = }")
            print(f"{dimensioning_max_flux_density = }")

            # get material properties, especially mu_r_imag
            temperature = 100
            material_flux_density_vec, material_mu_r_imag_vec, material_mu_r_real_vec = self.material_db.permeability_data_to_pro_file(temperature,
                                                                                                            fundamental_frequency,
                                                                                                            material_name,
                                                                                                            datasource=mdb.MaterialDataSource.ManufacturerDatasheet,
                                                                                                            datatype='permeability_data',
                                                                                                            plot_interpolation=False)
            mu_r_abs_vec = np.sqrt(material_mu_r_real_vec ** 2 + material_mu_r_imag_vec ** 2)



            for count, t1d_core_geometry_material in enumerate(t2_parameter_sweep):
                """
                inner optimization loop
                 * report about simulation progress
                 * generate winding matrix from sweep parameters
                 * calculate reluctance matrix
                
                """
                # report about simulation progress
                if count == simulations_per_percent * simulation_progress_percent:
                    simulation_progress_percent += 1
                    print(f"{simulation_progress_percent} simulation_progress_percent")

                # assign materials and values from sweep-tensor
                window_w = t1d_core_geometry_material[0]
                window_h_top = t1d_core_geometry_material[1]
                window_h_bot = t1d_core_geometry_material[2]
                core_inner_diameter = t1d_core_geometry_material[3]
                n_p_top = t1d_core_geometry_material[4]
                n_p_bot = t1d_core_geometry_material[5]
                n_s_top = t1d_core_geometry_material[6]
                n_s_bot = t1d_core_geometry_material[7]

                core_top_bot_height = core_inner_diameter / 4
                core_cross_section = (core_inner_diameter / 2) ** 2 * np.pi

                # generate winding matrix
                # note that the t11 winding-matrix will be reshaped later!
                t2_winding_matrix = [[n_p_top, n_s_top], [n_p_bot, n_s_bot]]

                # matrix reshaping
                t2_winding_matrix_transpose = np.transpose(t2_winding_matrix, (1, 0))

                t2_reluctance_matrix = self.t2_calculate_reluctance_matrix(t2_inductance_matrix, t2_winding_matrix, t2_winding_matrix_transpose)

                if np.linalg.det(t2_reluctance_matrix) != 0 and np.linalg.det(np.transpose(t2_winding_matrix)) != 0 and np.linalg.det(t2_inductance_matrix) != 0:
                    # calculate the flux
                    flux_top_vec, flux_bot_vec, flux_stray_vec = fr.flux_vec_from_current_vec(current_extracted_1_vec,
                                                                                              current_extracted_2_vec,
                                                                                              t2_winding_matrix,
                                                                                              t2_inductance_matrix)

                    # calculate maximum values
                    flux_top_max, flux_bot_max, flux_stray_max = fr.max_flux_from_flux_vec(flux_top_vec, flux_bot_vec,
                                                                                           flux_stray_vec)

                    flux_density_top_max = flux_top_max / core_cross_section
                    flux_density_bot_max = flux_bot_max / core_cross_section
                    flux_density_middle_max = flux_stray_max / core_cross_section

                    if (flux_density_top_max < dimensioning_max_flux_density) and (flux_density_bot_max < dimensioning_max_flux_density) and (flux_density_middle_max < dimensioning_max_flux_density):

                        # calculate target values for r_top and r_bot out of reluctance matrix
                        r_core_middle_cylinder_radial = fr.r_core_top_bot_radiant(core_inner_diameter, window_w, mu_r_abs, core_top_bot_height)

                        r_middle_target = -t2_reluctance_matrix[0][1]
                        r_top_target = t2_reluctance_matrix[0][0] - r_middle_target
                        r_bot_target = t2_reluctance_matrix[1][1] - r_middle_target

                        # calculate the core reluctance of top and bottom and middle part
                        r_core_top_cylinder_inner = fr.r_core_round(core_inner_diameter, window_h_top, mu_r_abs)
                        r_core_top = 2 * r_core_top_cylinder_inner + r_core_middle_cylinder_radial
                        r_air_gap_top_target = r_top_target - r_core_top

                        r_core_bot_cylinder_inner = fr.r_core_round(core_inner_diameter, window_h_bot, mu_r_abs)
                        r_core_bot = 2 * r_core_bot_cylinder_inner + r_core_middle_cylinder_radial
                        r_air_gap_bot_target = r_bot_target - r_core_bot

                        r_air_gap_middle_target = r_middle_target - r_core_middle_cylinder_radial


                        minimum_air_gap_length = 0
                        maximum_air_gap_length = 1e-3

                        l_top_air_gap = optimize.brentq(fr.r_air_gap_round_inf_sct, minimum_air_gap_length, maximum_air_gap_length, args=(core_inner_diameter, window_h_top, r_air_gap_top_target))
                        l_bot_air_gap = optimize.brentq(fr.r_air_gap_round_round_sct, minimum_air_gap_length, maximum_air_gap_length, args=(core_inner_diameter, window_h_bot / 2, window_h_bot / 2, r_air_gap_bot_target))
                        l_middle_air_gap = optimize.brentq(fr.r_air_gap_tablet_cylinder_sct, minimum_air_gap_length, maximum_air_gap_length, args=(core_inner_diameter, core_inner_diameter / 4, window_w, r_air_gap_middle_target))


                        if l_bot_air_gap > 0 and l_bot_air_gap > 0 and l_middle_air_gap > 0:

                            p_hyst_top = fr.hyst_losses_core_half_mu_r_imag(core_inner_diameter, window_h_top, window_w, mu_r_abs, flux_top_max, fundamental_frequency, material_flux_density_vec, material_mu_r_imag_vec)

                            p_hyst_middle = fr.power_losses_hysteresis_cylinder_radial_direction_mu_r_imag(flux_stray_max, core_inner_diameter/4,
                                                                                        core_inner_diameter/2, core_inner_diameter/2 + window_w, fundamental_frequency,
                                                                                        mu_r_abs, material_flux_density_vec, material_mu_r_imag_vec)

                            p_hyst_bot = fr.hyst_losses_core_half_mu_r_imag(core_inner_diameter, window_h_bot, window_w,
                                                                            mu_r_abs, flux_bot_max,
                                                                            fundamental_frequency,
                                                                            material_flux_density_vec,
                                                                            material_mu_r_imag_vec)

                            p_hyst = p_hyst_top + p_hyst_bot + p_hyst_middle

                            core_2daxi_total_volume = fr.calculate_core_2daxi_total_volume(core_inner_diameter, (window_h_bot + window_h_top + core_inner_diameter / 4), window_w)

                            valid_design_dict = ResultFile(
                                case= case_number,
                                air_gap_top= l_top_air_gap,
                                air_gap_bot= l_bot_air_gap,
                                air_gap_middle = l_middle_air_gap,
                                n_p_top= n_p_top,
                                n_p_bot= n_p_bot,
                                n_s_top= n_s_top,
                                n_s_bot= n_s_bot,
                                window_h_top= window_h_top,
                                window_h_bot= window_h_bot,
                                window_w= window_w,
                                mu_r_abs= mu_r_abs,
                                core_inner_diameter= core_inner_diameter,
                                flux_top_max= flux_top_max,
                                flux_bot_max= flux_bot_max,
                                flux_stray_max= flux_stray_max,
                                p_hyst= p_hyst,
                                core_2daxi_total_volume= core_2daxi_total_volume,
                            )



                            # Add dict to list of valid designs
                            valid_design_list.append(valid_design_dict)
                            case_number += 1

        print(f"Number of valid designs: {len(valid_design_list)}")
        return valid_design_list


