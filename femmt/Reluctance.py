import femmt as fmt
import numpy as np
import matplotlib
import os
from os import listdir
from os.path import isfile, join
import json
matplotlib.rc('xtick', labelsize=24)
matplotlib.rc('ytick', labelsize=24)
from itertools import product

# 3rd library imports
import numpy as np
from matplotlib import pyplot as plt

# femmt imports
import femmt.Functions as ff


def load_design(working_directory: str):
    """
    Load FEM simulation results from given working directory

    param working_directory: Sets the working directory
    :type working_directory: str
    """
    working_directories = []
    labels = []
    # working_directory = os.path.join(working_directory, 'fem_simulation_data')
    print("##########################")
    print(f"{working_directory =}")
    print("##########################")
    file_names = [f for f in listdir(working_directory) if isfile(join(working_directory, f))]

    counter = 0
    for name in file_names:
        temp_var = os.path.join(working_directory, name)
        working_directories.append(temp_var)
        labels.append(name.removesuffix('.json'))
        counter = counter + 1

    zip_iterator = zip(file_names, working_directories)
    logs = dict(zip_iterator)

    # After the simulations the sweep can be analyzed
    # This could be done using the FEMMTLogParser:
    log_parser = fmt.FEMMTLogParser(logs)

    # In this case the self inductivity of winding1 will be analyzed
    inductivities = []
    total_loss = []
    total_volume = []
    total_cost = []
    for name, data in log_parser.data.items():
        inductivities.append(data.sweeps[0].windings[0].self_inductance)
        total_loss.append(data.total_core_losses + data.total_winding_losses)
        total_volume.append(data.core_2daxi_total_volume)
        total_cost.append(data.total_cost)

    real_inductance = []
    for i in range(len(total_loss)):
        real_inductance.append(inductivities[i].real)

    return real_inductance, total_loss, total_volume, total_cost, labels

def plot_limitation():
    length = 15
    width = 100 * length
    height = 101 - length

    r_mx = 1 / (fmt.mu0 * (width / 2 / length + 2 / np.pi * (
                1 + np.log(np.pi * height / 4 / length))))
    print(height)

    # width_c = 100
    # length_c = 0.5
    # height1_c = np.linspace(50, 100, 1000)
    # height2_c =100 - height1_c
    # h_l = height2_c / length
    #
    # r_m1 = 1 / (fmt.mu0 * (width_c / 2 / length_c + 2 / np.pi * (
    #         1 + np.log(np.pi * height1_c / 4 / length_c))))
    #
    # r_m2 = 1 / (fmt.mu0 * (width_c / 2 / length_c + 2 / np.pi * (
    #         1 + np.log(np.pi * height2_c / 4 / length_c))))
    #
    # r_m = r_m1 + r_m2
    # combined = np.vstack((h_l, r_m)).T
    # print(combined)

    width_c = 100
    length_c = 0.5
    height1_c = np.linspace(50, 100, 1000)
    height2_c =100 - height1_c
    h_l = height2_c / length
    # print(h_l)
    r_m1 = 1 / (fmt.mu0 * (width_c / 2 / length_c + 2 / np.pi * (
            1 + np.log(np.pi * height1_c / 4 / length_c))))

    r_m2 = 1 / (fmt.mu0 * (width_c / 2 / length_c + 2 / np.pi * (
            1 + np.log(np.pi * height2_c / 4 / length_c))))

    r_m = r_m1 + r_m2
    ratio = r_mx / r_m
    # print(ratio)
    fig, ax = fmt.plt.subplots()  # Create a figure containing a single axes.
    fmt.plt.title("R_basic vs h/l")
    fmt.plt.xlabel("h/l")
    fmt.plt.ylabel("R_basic")

    ax.plot(h_l, r_m)
    # ax.hlines(y=r_m, xmin=0, xmax=50, linewidth=2, color='g')
    ax.hlines(y=r_mx, xmin=-1, xmax=51, linewidth=2, color='r')
    ax.invert_xaxis()
    ax.grid()
    plt.show()


def plot_r_basis():
    # width = 100
    # length = 1
    # height = np.linspace(100, 0, 1000)
    # h_l = height / length
    #
    # r_m = 1 / (fmt.mu0 * (width / 2 / length + 2 / np.pi * (
    #             1 + np.log(np.pi * height / 4 / length))))
    #
    # combined = np.vstack((h_l, r_m)).T
    # print(combined)
    # fig, ax = fmt.plt.subplots(figsize=(3.54, 3.54), dpi=150)  # Create a figure containing a single axes.
    # # fmt.plt.title("R_basic vs h/l")
    # fmt.plt.xlabel("$\dfrac{h}{l}$", fontsize=24)
    # fmt.plt.ylabel("$R_{\mathrm{basic}}^{\prime}$ / AT/Wb", fontsize=24)
    # ax.plot(h_l, r_m, linewidth=4, label=f'w/l ={width}')
    # ax.invert_xaxis()
    # ax.legend(fontsize=24)
    # ax.grid()
    # fmt.plt.show()

    # width = np.linspace(100, 0, 10000)
    # length = 1
    # height = 100
    # w_l = width / length
    #
    # r_m = 1 / (fmt.mu0 * (width / 2 / length + 2 / np.pi * (
    #         1 + np.log(np.pi * height / 4 / length))))
    #
    # combined = np.vstack((w_l, r_m)).T
    # print(combined)
    # fig, ax = fmt.plt.subplots(figsize=(3.54, 3.54), dpi=150)  # Create a figure containing a single axes.
    # # fmt.plt.title("R_basic vs w/l")
    # fmt.plt.xlabel("$\dfrac{w}{l}$", fontsize=24)
    # fmt.plt.ylabel("$R_{\mathrm{basic}}^{\prime}$ / AT/Wb", fontsize=24)
    # ax.plot(w_l, r_m, linewidth=4, label=f'h/l ={height}')
    # ax.invert_xaxis()
    # ax.legend(fontsize=24)
    # ax.grid()
    # fmt.plt.show()

    width = np.linspace(100, 0, 5)
    length = 1
    height = np.linspace(100, 0, 1000)
    h_l = height / length

    w_l = width / length
    fig, ax = fmt.plt.subplots(figsize=(3.54, 3.54), dpi=150)  # Create a figure containing a single axes.
    # fmt.plt.title("$R_{basic}$ vs $\dfrac{h}{l}$", fontsize=20)
    fmt.plt.xlabel("$\dfrac{h}{l}$", fontsize=24)
    fmt.plt.ylabel("$R_{\mathrm{basic}}^{\prime}$ / AT/Wb", fontsize=24)

    for i, wid in enumerate(width):
        r_m = 1 / (fmt.mu0 * (wid / 2 / length + 2 / np.pi * (
                1 + np.log(np.pi * height / 4 / length))))

        combined = np.vstack((h_l, r_m)).T
        # print(combined)

        ax.plot(h_l, r_m, linewidth=4, label=f'w/l ={w_l[i]}')

    ax.invert_xaxis()
    ax.set_yscale('log')
    ax.legend(fontsize=24)
    ax.grid()
    plt.show()





class MagneticCircuit:
    """This is a class for calculating the reluctance and inductance for 2D axis symmetric components
        and visualising magnetic circuit"""

    def __init__(self, core_inner_diameter: list, window_h: list, window_w: list, no_of_turns: list, n_air_gaps: list,
                 air_gap_h: list, air_gap_position: list, mu_rel: list, mult_air_gap_type: list = None,
                 air_gap_method: str = 'Percent', component_type: str = 'inductor', sim_type: str = 'single'):
        """
        :param core_inner_diameter: Diameter of center leg of the core in meter
        :type core_inner_diameter: list
        :param window_h: Height of the core window [in meter]
        :type window_h: list
        :param window_w: Width of the core window [in meter]
        :type window_w: list
        :param no_of_turns: Number of turns
        :type no_of_turns: list
        :param n_air_gaps: Number of air-gaps in the center leg of the core
        :type n_air_gaps: list
        :param air_gap_h: Air-gap height [in meter]
        :type air_gap_h: list
        :param air_gap_position: Position of the air-gap in the percentage with respect to window_h
        :type air_gap_position: list
        :param mu_rel: Relative permeability of the core [in F/m]
        :type mu_rel: list
        :param mult_air_gap_type: Two types of equally distributed air-gaps (used only for air-gaps more than 1)
            Type 1: Equally distributed air-gaps including corner air-gaps (eg: air-gaps-position = [0, 50, 100] for 3 air-gaps)
            Type 2: Equally distributed air-gaps excluding corner air-gaps (eg: air-gaps-position = [25, 50, 75] for 3 air-gaps)
        :type mult_air_gap_type: list
                :param air_gap_h: Air-gap height [in meter]
        :param air_gap_method: Input method of air gap position ( either in 'Percent', 'Center' or 'Manually')
        :type air_gap_method: str
        :param component_type: Position of the air-gap in the percentage with respect to window_h
        :type component_type: str
        :param sim_type: Relative permeability of the core [in F/m]
        :type sim_type: str
        """

        self.row_num = 0
        self.single_air_gap_len = None
        self.data_matrix_len = None
        self.data_matrix = None

        if not all(isinstance(item, int) for item in no_of_turns):
            raise Exception("no_of_turns list elements should be integer")
        if not all(isinstance(item, int) for item in n_air_gaps):
            raise Exception("n_air_gaps list elements should be integer")
        if not (air_gap_method == 'Center' or air_gap_method == 'Percent' or air_gap_method == 'Manually'):
            raise Exception("string value wrong for air_gap_method argument")
        if not (sim_type == 'single' or sim_type == 'sweep'):
            raise Exception("string value wrong for sim_type argument")
        if not (component_type == 'inductor' or component_type == 'integrated_transformer'):
            raise Exception("string value wrong for component_type argument")
        # if any(item > 0.0005 for item in air_gap_h):
        #     raise Exception("Model accuracy is not good for air_gap_h more than 0.0005")
        if sim_type == 'single':
            if not (len(core_inner_diameter) == 1 and len(window_h) == 1 and len(window_w) == 1 and len(no_of_turns) == 1
                    and len(n_air_gaps) == 1 and len(mu_rel) == 1):
                raise Exception("single sim_type requires single list elements")
            if not (n_air_gaps[0] == len(air_gap_h) and n_air_gaps[0] == len(air_gap_position)):
                raise Exception("No. of elements of air_gap_h and air_gap_position should match n_air_gaps")

        # Sort air_gap_position and air_gap_h based on air_gap_position


        self.core_inner_diameter = np.array(core_inner_diameter)
        self.window_h = np.array(window_h)
        self.window_w = np.array(window_w)
        self.mu_rel = np.array(mu_rel)
        self.no_of_turns = np.array(no_of_turns)
        self.n_air_gaps = np.array(n_air_gaps)

        self.mult_air_gap_type = None
        self.sim_type = sim_type
        self.data_matrix = np.zeros((1, 10))

        if not len(n_air_gaps):
            self.n_air_gaps = np.zeros((1, 1))

        if self.n_air_gaps[0] != 0:
            zipped_lists = zip(air_gap_position, air_gap_h)
            sorted_pairs = sorted(zipped_lists)

            tuples = zip(*sorted_pairs)
            self.air_gap_position, self.air_gap_h = [list(tpl) for tpl in tuples]
            self.air_gap_h = np.array(self.air_gap_h)
            self.percent_position_air_gap = np.array(self.air_gap_position)
        else:
            self.air_gap_h = np.array(air_gap_h)
            self.percent_position_air_gap = np.array(air_gap_position)

        # Creates the data matrix with all the input parameter combinations for sim_type = 'sweep'
        if sim_type == 'sweep':
            self.create_data_matrix(core_inner_diameter, window_h, window_w, no_of_turns, n_air_gaps, air_gap_h,
                                    air_gap_position, mu_rel, mult_air_gap_type)

            self.core_inner_diameter = self.data_matrix[:, 0]
            self.window_h = self.data_matrix[:, 1]
            self.window_w = self.data_matrix[:, 2]
            self.mu_rel = self.data_matrix[:, 3]
            self.no_of_turns = self.data_matrix[:, 4]
            self.n_air_gaps = self.data_matrix[:, 5]
            self.air_gap_h = self.data_matrix[:, 6]
            self.percent_position_air_gap = self.data_matrix[:, 7]
            self.mult_air_gap_type = self.data_matrix[:, 8]

        self.core_h = self.window_h + self.core_inner_diameter / 2
        self.r_outer = None
        self.r_inner = None
        self.core_h_middle = None  # height of upper and lower part of the window in the core
        self.outer_w = None  # Outer leg width
        self.mu_0 = 4 * np.pi * 1e-7

        self.abs_position_air_gap = None
        self.air_gap_method = air_gap_method

        self.cal_inductance = None
        self.section = None
        self.length = None
        self.area = None
        self.reluctance = None

        self.max_percent_position = None
        self.min_percent_position = None
        self.param_pos_dict = None

        if component_type == 'inductor':
            self.core_reluctance()
            if sim_type == "single":
                self.air_gap_reluctance_single()
            elif sim_type == "sweep":
                self.air_gap_reluctance_sweep()

    def create_data_matrix(self, core_inner_diameter: list, window_h: list, window_w: list, no_of_turns: list, n_air_gaps: list,
                           air_gap_h: list, air_gap_position: list, mu_rel: list, mult_air_gap_type: list):

        """ Creates matrix consisting of input design parameters with all their combinations

        :param core_inner_diameter: Diameter of center leg of the core in meter
        :type core_inner_diameter: list
        :param window_h: Height of the core window [in meter]
        :type window_h: list
        :param window_w: Width of the core window [in meter]
        :type window_w: list
        :param no_of_turns: Number of turns
        :type no_of_turns: list
        :param n_air_gaps: Number of air-gaps in the center leg of the core
        :type n_air_gaps: list
        :param air_gap_h: Air-gap height [in meter]
        :type air_gap_h: list
        :param air_gap_position: Position of the air-gap in the percentage with respect to window_h
        :type air_gap_position: list
        :param mu_rel: Relative permeability of the core [in F/m]
        :type mu_rel: list
        :param mult_air_gap_type: Two types of equally distributed air-gaps (used only for air-gaps more than 1)
            Type 1: Equally distributed air-gaps including corner air-gaps (eg: air-gaps-position = [0, 50, 100] for 3 air-gaps)
            Type 2: Equally distributed air-gaps excluding corner air-gaps (eg: air-gaps-position = [25, 50, 75] for 3 air-gaps)
        :type mult_air_gap_type: list
        """

        # example: data_matrix = [core_inner_diameter, window_h, window_w, mu_rel, no_of_turns, n_air_gaps, air_gap_h,
        #                      air_gap_position, mult_air_gap_type, inductance]
        clone_n_air_gaps = n_air_gaps

        if 1 in clone_n_air_gaps:
            self.data_matrix = np.zeros((len(core_inner_diameter) * len(no_of_turns) * len(air_gap_h) * len(mu_rel) * (
                    len(air_gap_position) + (len(n_air_gaps) - 1) * len(mult_air_gap_type)), 10))
        else:
            self.data_matrix = np.zeros((len(core_inner_diameter) * len(no_of_turns) * len(air_gap_h) * len(n_air_gaps) * len(
                mu_rel) * len(mult_air_gap_type), 10))

        if 1 in clone_n_air_gaps:
            clone_n_air_gaps.remove(1)
            for index_1 in range(len(core_inner_diameter)):
                for index_2 in range(len(mu_rel)):
                    for index_3 in range(len(no_of_turns)):
                        for index_4 in range(len(air_gap_h)):
                            for index_5 in range(len(air_gap_position)):
                                self.data_matrix[self.row_num, 0] = core_inner_diameter[index_1]
                                self.data_matrix[self.row_num, 1] = window_h[index_1]
                                self.data_matrix[self.row_num, 2] = window_w[index_1]
                                self.data_matrix[self.row_num, 3] = mu_rel[index_2]
                                self.data_matrix[self.row_num, 4] = no_of_turns[index_3]
                                self.data_matrix[self.row_num, 5] = 1
                                self.data_matrix[self.row_num, 6] = air_gap_h[index_4]
                                self.data_matrix[self.row_num, 7] = air_gap_position[index_5]
                                self.data_matrix[self.row_num, 8] = 0
                                self.row_num = self.row_num + 1

        self.single_air_gap_len = self.row_num

        if len(clone_n_air_gaps) != 0:
            for index_1 in range(len(core_inner_diameter)):
                for index_2 in range(len(mu_rel)):
                    for index_3 in range(len(no_of_turns)):
                        for index_4 in range(len(clone_n_air_gaps)):
                            for index_5 in range(len(air_gap_h)):
                                for index_6 in range(len(mult_air_gap_type)):
                                    self.data_matrix[self.row_num, 0] = core_inner_diameter[index_1]
                                    self.data_matrix[self.row_num, 1] = window_h[index_1]
                                    self.data_matrix[self.row_num, 2] = window_w[index_1]
                                    self.data_matrix[self.row_num, 3] = mu_rel[index_2]
                                    self.data_matrix[self.row_num, 4] = no_of_turns[index_3]
                                    self.data_matrix[self.row_num, 5] = n_air_gaps[index_4]
                                    self.data_matrix[self.row_num, 6] = air_gap_h[index_5]
                                    self.data_matrix[self.row_num, 7] = 0
                                    self.data_matrix[self.row_num, 8] = mult_air_gap_type[index_6]
                                    self.row_num = self.row_num + 1

        self.data_matrix_len = self.row_num

    def core_reluctance(self):
        """Calculates the core reluctance along with length and area of each section of the core geometry"""

        self.core_h_middle = (self.core_h - self.window_h) / 2
        self.r_inner = self.core_inner_diameter / 2 + self.window_w
        self.r_outer = np.sqrt((self.core_inner_diameter / 2) ** 2 + self.r_inner ** 2)
        self.outer_w = self.r_outer - self.r_inner

        self.section = [0, 1, 2, 3,
                        4]  # Section [0]: Center leg; [1]: Central corner; [2]: Winding window; [3]: Outer corners; [4]: Outer leg
        self.length = np.zeros((len(self.data_matrix), len(self.section)))
        self.area = np.zeros((len(self.data_matrix), len(self.section)))
        self.reluctance = np.zeros((len(self.data_matrix), len(self.section) + 1))

        if self.sim_type == 'sweep':
            self.length[:, 0] = self.window_h - (self.n_air_gaps * self.air_gap_h)
        else:
            self.length[:, 0] = self.window_h - sum(self.air_gap_h)

        self.area[:, 0] = np.pi * ((self.core_inner_diameter / 2) ** 2)

        self.length[:, 1] = (np.pi / 8) * (self.core_inner_diameter / 2 + self.core_h_middle)
        # self.area[:, 1] = ((self.core_inner_diameter / 2 + self.core_h_middle) / 2) * 2 * np.pi * (self.core_inner_diameter / 2)
        self.area[:, 1] = np.pi * self.core_inner_diameter / 2 * ((self.core_inner_diameter / 2 + self.core_h_middle) / 2) * (
                2 + ((self.core_inner_diameter / 2 + self.core_h_middle) / 2) / (
            np.sqrt((self.core_inner_diameter / 2) ** 2 + self.core_h_middle ** 2)))

        self.length[:, 2] = self.window_w
        self.area[:, 2] = np.nan
        self.reluctance[:, 2] = ((self.mu_0 * self.mu_rel * 2 * np.pi * self.core_h_middle) ** -1) * np.log(
            (2 * self.r_inner) / self.core_inner_diameter)

        self.length[:, 3] = (np.pi / 8) * (self.outer_w + self.core_h_middle)
        # self.area[:, 3] = ((self.outer_w + self.core_h_middle) / 2) * 2 * np.pi * self.r_inner
        self.area[:, 3] = 2 * np.pi * self.r_inner * self.core_h_middle

        self.length[:, 4] = self.window_h
        self.area[:, 4] = np.pi * (self.r_outer ** 2 - self.r_inner ** 2)

        self.reluctance[:, 0] = self.length[:, 0] / (self.mu_0 * self.mu_rel * self.area[:, 0])
        self.reluctance[:, 1] = self.length[:, 1] / (self.mu_0 * self.mu_rel * self.area[:, 1])
        self.reluctance[:, 3] = self.length[:, 3] / (self.mu_0 * self.mu_rel * self.area[:, 3])
        self.reluctance[:, 4] = self.length[:, 4] / (self.mu_0 * self.mu_rel * self.area[:, 4])
        self.reluctance[:, 1:4] = 2 * self.reluctance[:, 1:4]

    def air_gap_reluctance_sweep(self):
        """Calculates air-gap reluctance and the inductance of the given geometry"""

        # Single air-gap reluctance calculations
        self.max_percent_position = ((self.window_h[0:self.single_air_gap_len] - (
                self.air_gap_h[0:self.single_air_gap_len] / 2)) / self.window_h[0:self.single_air_gap_len]) * 100
        self.min_percent_position = ((self.air_gap_h[0:self.single_air_gap_len] / 2) / self.window_h[
                                                                                       0:self.single_air_gap_len]) * 100
        self.abs_position_air_gap = (self.percent_position_air_gap[0:self.single_air_gap_len] * self.window_h[
                                                                                                0:self.single_air_gap_len]) / 100  # Convert percent position to absolute value position
        h = np.zeros((len(self.abs_position_air_gap), 2))

        h[:, 0] = np.where((self.percent_position_air_gap[0:self.single_air_gap_len] <= self.min_percent_position) | (
                self.percent_position_air_gap[0:self.single_air_gap_len] >= self.max_percent_position),
                           self.window_h[0:self.single_air_gap_len] - self.air_gap_h[0:self.single_air_gap_len],
                           self.abs_position_air_gap - (self.air_gap_h[0:self.single_air_gap_len] / 2))
        h[:, 1] = np.where((self.percent_position_air_gap[0:self.single_air_gap_len] <= self.min_percent_position) | (
                self.percent_position_air_gap[0:self.single_air_gap_len] >= self.max_percent_position),
                           0,
                           self.window_h[0:self.single_air_gap_len] - self.abs_position_air_gap - (
                                   self.air_gap_h[0:self.single_air_gap_len] / 2))

        self.reluctance[0:self.single_air_gap_len, 5] = np.where(h[:, 1] == 0,
                                                                 single_round_inf(
                                                                     self.air_gap_h[0:self.single_air_gap_len],
                                                                     self.core_inner_diameter[0:self.single_air_gap_len], h[:, 0]),
                                                                 single_round_round(
                                                                     self.air_gap_h[0:self.single_air_gap_len],
                                                                     self.core_inner_diameter[0:self.single_air_gap_len], h[:, 0],
                                                                     h[:, 1]))

        # Distributed air-gaps reluctance calculations
        h_multiple = np.where(self.mult_air_gap_type[self.single_air_gap_len:self.data_matrix_len] == 1,
                              (self.window_h[self.single_air_gap_len:self.data_matrix_len] - (
                                      self.n_air_gaps[self.single_air_gap_len:self.data_matrix_len] * self.air_gap_h[
                                                                                                      self.single_air_gap_len:self.data_matrix_len])) / (
                                      (self.n_air_gaps[self.single_air_gap_len:self.data_matrix_len] - 1) * 2),
                              (self.window_h[self.single_air_gap_len:self.data_matrix_len] - (
                                      self.n_air_gaps[self.single_air_gap_len:self.data_matrix_len] * self.air_gap_h[
                                                                                                      self.single_air_gap_len:self.data_matrix_len])) / (
                                      self.n_air_gaps[self.single_air_gap_len:self.data_matrix_len] + 1))

        self.reluctance[self.single_air_gap_len:self.data_matrix_len, 5] = np.where(
            self.mult_air_gap_type[self.single_air_gap_len:self.data_matrix_len] == 1,
            distributed_type_1(self.air_gap_h[self.single_air_gap_len:self.data_matrix_len],
                               self.core_inner_diameter[self.single_air_gap_len:self.data_matrix_len],
                               self.n_air_gaps[self.single_air_gap_len:self.data_matrix_len], h_multiple),
            distributed_type_2(self.air_gap_h[self.single_air_gap_len:self.data_matrix_len],
                               self.core_inner_diameter[self.single_air_gap_len:self.data_matrix_len],
                               self.n_air_gaps[self.single_air_gap_len:self.data_matrix_len], h_multiple))

        # Inductance calculation
        self.data_matrix[:, 9] = (self.no_of_turns ** 2) / np.sum(self.reluctance, axis=1)
        self.data_matrix = np.hstack(
            (self.data_matrix, np.reshape(self.core_h_middle, (self.data_matrix_len, 1))))  # position: 10
        self.data_matrix = np.hstack(
            (self.data_matrix, np.reshape(self.r_inner, (self.data_matrix_len, 1))))  # position: 11
        self.data_matrix = np.hstack(
            (self.data_matrix, np.reshape(self.r_outer, (self.data_matrix_len, 1))))  # position: 12
        self.data_matrix = np.hstack(
            (self.data_matrix, np.reshape(self.area[:, 0], (self.data_matrix_len, 1))))  # position: 13
        self.data_matrix = np.hstack(
            (self.data_matrix, np.reshape(self.area[:, 4], (self.data_matrix_len, 1))))  # position: 14
        self.data_matrix = np.hstack(
            (self.data_matrix, np.reshape(self.core_h, (self.data_matrix_len, 1))))  # position: 15

    def get_parameters_position_dict(self):
        self.param_pos_dict = {"core_inner_diameter": 0, "window_h": 1, "window_w": 2, "mu_rel": 3, "no_of_turns": 4,
                               "n_air_gaps": 5,
                               "air_gap_h": 6, "air_gap_position": 7, "mult_air_gap_type": 8, "inductance": 9,
                               "core_h_middle": 10,
                               "r_inner": 11, "r_outer": 12, "center_leg_area": 13, "outer_leg_area": 14, "core_h": 15}

        return self.param_pos_dict

    def air_gap_reluctance_single(self):
        flag_0 = 0
        flag_1 = 0
        flag_2 = 0
        if self.n_air_gaps[0] != 0:
            if self.air_gap_method == 'Center':
                self.section.append(6)  # round-round type airgap
                temp1 = ff.r_basis([self.air_gap_h[0] / 2], [self.core_inner_diameter], [(self.window_h - self.air_gap_h[0]) / 2])
                temp2 = ff.sigma([self.air_gap_h[0]], [self.core_inner_diameter / 2], 2 * temp1)
                temp3 = ff.r_round_round([self.air_gap_h[0]], temp2, [self.core_inner_diameter / 2])
                temp4 = self.air_gap_h[0] / (self.mu_0 * np.pi * (self.core_inner_diameter / 2) ** 2)  # classical reluctance formula
                self.reluctance[:, 5] = temp3

            elif self.air_gap_method == 'Percent' or self.air_gap_method == 'Manually':
                self.max_percent_position = ((self.window_h - self.air_gap_h[self.n_air_gaps - 1] / 2) / self.window_h) * 100
                self.min_percent_position = ((self.air_gap_h[0] / 2) / self.window_h) * 100
                if self.air_gap_method == 'Percent':
                    self.position = np.array(self.percent_position_air_gap) / 100 * self.window_h  # Convert percent position to absolute value position
                print(f"Max percent: {self.max_percent_position}")
                print(f"Min percent: {self.min_percent_position}")

                if self.percent_position_air_gap[0] <= self.min_percent_position:
                    flag_0 = 1
                    self.section.append(8)
                    if self.n_air_gaps == 1:
                        h = self.window_h - self.air_gap_h[0]
                    else:
                        h = ((self.position[1] - self.air_gap_h[1] / 2) - self.air_gap_h[0]) / 2

                    temp1 = ff.r_basis([self.air_gap_h[0]], [self.core_inner_diameter], [h])
                    temp2 = ff.sigma([self.air_gap_h[0]], [self.core_inner_diameter / 2], temp1)
                    temp3 = ff.r_round_inf([self.air_gap_h[0]], temp2, [self.core_inner_diameter / 2])
                    self.reluctance[:, 5] = self.reluctance[:, 5] + temp3
                    print('air gap is at lower corner')

                if self.percent_position_air_gap[self.n_air_gaps - 1] >= self.max_percent_position:
                    flag_1 = 1
                    self.section.append(8)
                    if self.n_air_gaps == 1:
                        h = self.window_h - self.air_gap_h[self.n_air_gaps - 1]
                    else:
                        h = (self.position[self.n_air_gaps - 1] - self.position[self.n_air_gaps - 2] - self.air_gap_h[
                            self.n_air_gaps - 1] / 2 - self.air_gap_h[self.n_air_gaps - 2] / 2) / 2

                    temp1 = ff.r_basis([self.air_gap_h[self.n_air_gaps - 1]], [self.core_inner_diameter], [h])
                    temp2 = ff.sigma([self.air_gap_h[self.n_air_gaps - 1]], [self.core_inner_diameter / 2], temp1)
                    temp3 = ff.r_round_inf([self.air_gap_h[self.n_air_gaps - 1]], temp2, [self.core_inner_diameter / 2])
                    self.reluctance[:, 5] = self.reluctance[:, 5] + temp3
                    print('air gap is at upper corner')

                for i in range(self.n_air_gaps[0]):
                    if self.min_percent_position < self.percent_position_air_gap[i] < self.max_percent_position:
                        self.section.append(7)
                        if flag_2 == 0:
                            if flag_0 == 0 and flag_1 == 0:  # No corner air-gaps
                                self.position = np.append(self.position,
                                                          self.window_h + (self.window_h - self.position[self.n_air_gaps - 1]))
                                self.position = np.insert(self.position, 0, -self.position[0])
                                self.air_gap_h = np.append(self.air_gap_h, self.air_gap_h[self.n_air_gaps - 1])
                                self.air_gap_h = np.insert(self.air_gap_h, 0, self.air_gap_h[0])
                            elif flag_0 == 1 and flag_1 == 0:  # Only lower air-gap is present
                                self.position = np.append(self.position,
                                                          self.window_h + (self.window_h - self.position[self.n_air_gaps - 1]))
                                self.air_gap_h = np.append(self.air_gap_h, self.air_gap_h[self.n_air_gaps - 1])
                            elif flag_0 == 0 and flag_1 == 1:  # Only Upper air-gap is present
                                self.position = np.insert(self.position, 0, -self.position[0])
                                self.air_gap_h = np.insert(self.air_gap_h, 0, self.air_gap_h[0])
                            flag_2 = 1

                        if flag_0 == 0 and flag_1 == 0:
                            h1 = (self.position[i + 1] - self.position[i] - self.air_gap_h[i + 1] / 2 - self.air_gap_h[
                                i] / 2) / 2
                            h2 = (self.position[i + 2] - self.position[i + 1] - self.air_gap_h[i + 2] / 2 - self.air_gap_h[
                                i + 1] / 2) / 2
                            print('No corner air gap detected')
                        elif flag_0 == 1 and flag_1 == 0:
                            h1 = (self.position[i] - self.position[i - 1] - self.air_gap_h[i] / 2 - self.air_gap_h[
                                i - 1] / 2) / 2
                            h2 = (self.position[i + 1] - self.position[i] - self.air_gap_h[i + 1] / 2 - self.air_gap_h[
                                i] / 2) / 2
                            print('Lower air gap detected')
                        elif flag_0 == 0 and flag_1 == 1:
                            h1 = (self.position[i + 1] - self.position[i] - self.air_gap_h[i + 1] / 2 - self.air_gap_h[
                                i] / 2) / 2
                            h2 = (self.position[i + 2] - self.position[i + 1] - self.air_gap_h[i + 2] / 2 - self.air_gap_h[
                                i + 1] / 2) / 2
                            print('Upper air gap detected')
                        else:
                            h1 = (self.position[i] - self.position[i - 1] - self.air_gap_h[i] / 2 - self.air_gap_h[
                                i - 1] / 2) / 2
                            h2 = (self.position[i + 1] - self.position[i] - self.air_gap_h[i + 1] / 2 - self.air_gap_h[
                                i] / 2) / 2
                            print('Both air gap detected')

                        r_basis_1 = ff.r_basis([self.air_gap_h[i] / 2], [self.core_inner_diameter], [h1])
                        r_basis_2 = ff.r_basis([self.air_gap_h[i] / 2], [self.core_inner_diameter], [h2])
                        temp2 = ff.sigma([self.air_gap_h[i]], [self.core_inner_diameter / 2], r_basis_1 + r_basis_2)
                        temp3 = ff.r_round_round([self.air_gap_h[i]], temp2, [self.core_inner_diameter / 2])
                        self.reluctance[:, 5] = self.reluctance[:, 5] + temp3

    def calculate_inductance(self):
        # self.section, self.orientation = set_orientation(self.section, len(self.section))
        self.cal_inductance = (self.no_of_turns * self.no_of_turns) / np.sum(self.reluctance, axis=1)
        self.data_matrix[:, 9] = self.cal_inductance


def single_round_inf(air_gap_h, core_inner_diameter, height_core_material):
    """Returns reluctance of a single air-gap at the corner

    :param air_gap_h: Air-gap height [in meter]
    :type air_gap_h: list
    :param core_inner_diameter: Diameter of center leg of the core [in meter]
    :type core_inner_diameter: list
    :param height_core_material: Core distance between air-gap and other end of the window-h [in meter]
    :type height_core_material: list
    :return: Reluctance of a single air-gap at the corner
    :rtype: list"""

    r_basis_round_inf = ff.r_basis(air_gap_h, core_inner_diameter, height_core_material)
    sigma_round_inf = ff.sigma(air_gap_h, core_inner_diameter / 2, r_basis_round_inf)
    reluctance_round_inf = ff.r_round_inf(air_gap_h, sigma_round_inf, core_inner_diameter / 2)

    return reluctance_round_inf


def single_round_round(air_gap_total_hight, core_inner_diameter, height_core_material_0, height_core_material_1):
    """Returns reluctance of a single air-gap at position other than corner on the center leg

    :param air_gap_total_hight: Air-gap total height [in meter]
    :type air_gap_total_hight: list
    :param core_inner_diameter: Diameter of center leg of the core [in meter]
    :type core_inner_diameter: list
    :param height_core_material_0: Distance between window_h and air_gap_position for a single air-gap [in meter]
    :type height_core_material_0: list
    :param height_core_material_1: Height of air-gap from the base of the core window [in meter]
    :type height_core_material_1: list
    :return: Reluctance of a single air-gap at position other than corner on the center leg
    :rtype: list"""

    r_basis_1 = ff.r_basis(air_gap_total_hight / 2, core_inner_diameter, height_core_material_0)
    r_basis_2 = ff.r_basis(air_gap_total_hight / 2, core_inner_diameter, height_core_material_1)
    sigma_round_round = ff.sigma(air_gap_total_hight, core_inner_diameter / 2, r_basis_1 + r_basis_2)
    reluctance_round_round = ff.r_round_round(air_gap_total_hight, sigma_round_round, core_inner_diameter / 2)

    return reluctance_round_round


def distributed_type_1(air_gap_hight_single_air_gap, core_inner_diameter, n_air_gaps, h_multiple):
    """Returns distributed air-gap reluctance of Type 1 (Where corner air-gaps are present)

    :param air_gap_hight_single_air_gap: Air-gap height [in meter]
    :type air_gap_hight_single_air_gap: list
    :param core_inner_diameter: Diameter of center leg of the core [in meter]
    :type core_inner_diameter: list
    :param n_air_gaps: Number of air-gaps in the center leg of the core
    :type n_air_gaps: list
    :param h_multiple: Half of core height between two consecutive air-gaps in an equally distributed air-gaps [in meter]
    :type h_multiple: list
    :return: Distributed air-gap reluctance of Type 1 (Where corner air-gaps are present)
    :rtype: list"""

    # ToDo: Raise Error for less than two air gaps

    # first part calculates the two outer air gaps (very top and very bottom)
    r_basis = ff.r_basis(air_gap_hight_single_air_gap, core_inner_diameter, h_multiple)
    sigma_round_inf = ff.sigma(air_gap_hight_single_air_gap, core_inner_diameter / 2, r_basis)
    reluctance_round_inf = ff.r_round_inf(air_gap_hight_single_air_gap, sigma_round_inf, core_inner_diameter / 2)
    reluctance = (2 * reluctance_round_inf)

    # second part calculates the inner air gaps between top and bottom air gaps (if available)
    r_basis_1 = ff.r_basis(air_gap_hight_single_air_gap / 2, core_inner_diameter, h_multiple)
    r_basis_2 = ff.r_basis(air_gap_hight_single_air_gap / 2, core_inner_diameter, h_multiple)
    sigma_round_round = ff.sigma(air_gap_hight_single_air_gap, core_inner_diameter / 2, r_basis_1 + r_basis_2)
    reluctance_round_round = ff.r_round_round(air_gap_hight_single_air_gap, sigma_round_round, core_inner_diameter / 2)
    reluctance = reluctance + ((n_air_gaps - 2) * reluctance_round_round)

    return reluctance


def distributed_type_2(air_gap_hight_single_air_gap, core_inner_diameter, n_air_gaps, h_multiple):
    """Returns distributed air-gap reluctance of Type 2 (Where corner air-gaps are absent)

    :param air_gap_hight_single_air_gap: Air-gap height [in meter]
    :type air_gap_hight_single_air_gap: list
    :param core_inner_diameter: Diameter of center leg of the core [in meter]
    :type core_inner_diameter: list
    :param n_air_gaps: Number of air-gaps in the center leg of the core
    :type n_air_gaps: list
    :param h_multiple: Core height between two consecutive air-gaps in an equally distributed air-gaps [in meter]
    :type h_multiple: list
    :return: Distributed air-gap reluctance of Type 2 (Where corner air-gaps are absent)
    :rtype: list"""

    #ToDo: Raise Error for less than two air gaps

    # First part calculates two outer air gaps (very top and very bottom)
    r_basis_airgap_airgap = ff.r_basis(air_gap_hight_single_air_gap / 2, core_inner_diameter, h_multiple)
    r_basis_airgap_corner = ff.r_basis(air_gap_hight_single_air_gap / 2, core_inner_diameter, h_multiple / 2)
    sigma_round_round = ff.sigma(air_gap_hight_single_air_gap, core_inner_diameter / 2, r_basis_airgap_airgap + r_basis_airgap_corner)
    reluctance_round_round = ff.r_round_round(air_gap_hight_single_air_gap, sigma_round_round, core_inner_diameter / 2)
    reluctance = (2 * reluctance_round_round)

    # second part calculates air gaps between the outer air gaps (if available)
    r_basis_airgap_airgap_1 = ff.r_basis(air_gap_hight_single_air_gap / 2, core_inner_diameter, h_multiple / 2)
    r_basis_airgap_airgap_2 = ff.r_basis(air_gap_hight_single_air_gap / 2, core_inner_diameter, h_multiple / 2)
    sigma_round_round = ff.sigma(air_gap_hight_single_air_gap, core_inner_diameter / 2, r_basis_airgap_airgap_1 + r_basis_airgap_airgap_2)
    reluctance_round_round = ff.r_round_round(air_gap_hight_single_air_gap, sigma_round_round, core_inner_diameter / 2)
    reluctance = reluctance + ((n_air_gaps - 2) * reluctance_round_round)

    return reluctance


if __name__ == '__main__':
    mc1 = MagneticCircuit(core_inner_diameter=[0.0149], window_h=[0.0295], window_w=[0.01105], no_of_turns=[9], n_air_gaps=[3],
                          air_gap_h=[0.0005, 0.0005, 0.0005], air_gap_position=[23.3, 25, 26.7], mu_rel=[3000], mult_air_gap_type=[1, 2],
                          air_gap_method='Percent', component_type='inductor', sim_type='single')  # 0.0149
    #
    # mc1 = MagneticCircuit(core_inner_diameter=[0.0149], window_h=[0.0295], window_w=[0.01105], no_of_turns=[9],
    #                       n_air_gaps=[1],
    #                       air_gap_h=[0.0005], air_gap_position=[26.7], mu_rel=[3000],
    #                       mult_air_gap_type=[1, 2],
    #                       air_gap_method='Percent', component_type='inductor', sim_type='single')  # 0.0149
    # # print(np.sum(mc1.reluctance, axis=1) - mc1.reluctance[0, 5])
    # print(mc1.reluctance[0, 5])
    mc1.calculate_inductance()
    #
    # print(f"Inductance is {mc1.cal_inductance}")
    # plot_r_basis()
    # plot_limitation()

    # # Sweep of air-gap and winding position and compare it with FEM simulation
    # sweep_air_gap_h = np.linspace(0.0001, 0.001, 5)
    # sweep_wndg_pos = np.linspace(0.001, 0.007, 7)
    # # sweep_current = np.linspace(0.1, 10, 5)
    # fem_ind = np.zeros((len(sweep_air_gap_h), len(sweep_wndg_pos)))
    # cal_ind = np.zeros((len(sweep_air_gap_h), len(sweep_wndg_pos)))
    # # fem_ind = np.zeros(len(sweep_wndg_pos))
    # # cal_ind = np.zeros(len(sweep_wndg_pos))
    # # cal_fringe_dist = np.zeros((len(sweep_air_gap_h), len(sweep_wndg_pos)))
    #
    # example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
    # if not os.path.exists(example_results_folder):
    #     os.mkdir(example_results_folder)
    #
    # working_directory = os.path.join(example_results_folder, "inductor")
    # if not os.path.exists(working_directory):
    #     os.mkdir(working_directory)
    #
    # for j in range(len(sweep_air_gap_h)):
    #     for k in range(len(sweep_wndg_pos)):
    #             mc1 = MagneticCircuit(core_inner_diameter=[0.0149], window_h=[0.0295], window_w=[0.01105], no_of_turns=[9], n_air_gaps=[1],
    #                           air_gap_h=[sweep_air_gap_h[j]], air_gap_position=[50], mu_rel=[3000], mult_air_gap_type=[1, 2],
    #                           air_gap_method='Percent', component_type='inductor', sim_type='single')  # 0.0149
    #             mc1.calculate_inductance()
    #             cal_ind[j, k] = mc1.data_matrix[:, 9]  # - (2.38295 * 1e-6 - 0.000326175 * sweep_wndg_pos[k])
    #             # cal_fringe_dist[j, k] = mc1.fringe_dist
    #             # mc1.max_percent = ((mc1.window_h - (sweep_air_gap_h[j] / 2)) / mc1.window_h) * 100
    #             # mc1.min_percent = ((sweep_air_gap_h[j] / 2) / mc1.window_h) * 100
    #
    #             # 1. chose simulation type
    #             geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor,
    #                                         working_directory=working_directory,
    #                                         silent=True)
    #
    #             # 2. set core parameters
    #             core_db = fmt.core_database()["PQ 40/40"]
    #
    #             core = fmt.Core(core_inner_diameter=core_db["core_inner_diameter"], window_w=core_db["window_w"],
    #                             window_h=core_db["window_h"],
    #                             material="N95", temperature=25, frequency=100000, datasource="manufacturer_datasheet")
    #             # mu_rel=3000, phi_mu_deg=10,
    #             # sigma=0.5)
    #             geo.set_core(core)
    #
    #             # 3. set air gap parameters
    #             air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    #             air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, sweep_air_gap_h[j], 50)
    #             # air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 90)
    #             geo.set_air_gaps(air_gaps)
    #
    #             # 4. set insulations
    #             insulation = fmt.Insulation()
    #             insulation.add_core_insulations(0.001, 0.001, sweep_wndg_pos[k], 0.001)
    #             insulation.add_winding_insulations([0.0005], 0.0001)
    #             geo.set_insulation(insulation)
    #
    #             # 5. create winding window and virtual winding windows (vww)
    #             winding_window = fmt.WindingWindow(core, insulation)
    #             vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)
    #
    #             # 6. create conductor and set parameters: use solid wires
    #             winding = fmt.Conductor(0, fmt.Conductivity.Copper)
    #             # winding.set_solid_round_conductor(conductor_radius=0.0013,
    #             #                                   conductor_arrangement=fmt.ConductorArrangement.Square)
    #             winding.set_litz_round_conductor(conductor_radius=0.0013, number_strands=150, strand_radius=100e-6,
    #             fill_factor=None, conductor_arrangement=fmt.ConductorArrangement.Square)
    #
    #             # 7. add conductor to vww and add winding window to MagneticComponent
    #             vww.set_winding(winding, 9, None)
    #             geo.set_winding_window(winding_window)
    #
    #             # 8. create the model
    #             geo.create_model(freq=100000, visualize_before=False, save_png=False)
    #
    #             # 6.a. start simulation
    #             geo.single_simulation(freq=100000, current=[4.5], show_results=False)
    #
    #             # basic_example_func(sweep_air_gap_h[j], [mc1.min_percent + 0.1, 25, 50, 75, mc1.max_percent - 0.1], 9, sweep_wndg_pos[k], 10)
    #             # basic_example_func(sweep_air_gap_h[j], [50], 9, sweep_wndg_pos[k], 10)
    #             fem_ind[j, k] = geo.read_log()["single_sweeps"][0]["winding1"]["self_inductance"][0]
    #             # Load design and plot various plots for analysis
    #             # working_directory_load = os.path.join(working_directory, "results")
    #             # inductance, total_loss, total_volume, total_cost, annotation_list = load_design \
    #             #     (working_directory=working_directory_load)
    #             # print(inductance)
    #             # fem_ind[j, k] = inductance
    #             # fem_ind[j, k] = geo.read_log()["single_sweeps"][0]["winding1"]["Q"]
    #
    #
    # # print(f"Air-gap length: {sweep_air_gap_h}")
    # print(f"Winding position: {sweep_wndg_pos}")
    # print(f"FEM inductance: {fem_ind}")
    # print(f"Calculated inductance: {cal_ind}")
    # # print(f"Calculated fringe dist: {cal_fringe_dist}")
    #
    # # h_by_l = ((0.0295 - sweep_air_gap_h) / 2) / sweep_air_gap_h
    #
    # abs_error = cal_ind - fem_ind
    # error = (abs_error / fem_ind) * 100
    # # avg_abs_error = np.sum(abs_error, axis=0) / 5
    #
    # # print(f"h_by_l: {h_by_l}")
    # print(f"abs_error: {abs_error}")
    # print(f"percent error: {error}")
    # # print(f"average actual error: {avg_abs_error}")
    #
    # # np.savetxt('absolute_error.txt', abs_error)
    #
    # print(type(error))
    # print(type(sweep_wndg_pos))
    #
    #
    # # Plotting tools
    # fig, ax = fmt.plt.subplots(figsize=(3.54, 3.54), dpi=150)  # Create a figure containing a single axes.
    # fmt.plt.title("(9 litz-conductors) (Center air-gap)", fontsize=24)
    # fmt.plt.xlabel("Winding position (in m)", fontsize=24)
    # fmt.plt.ylabel("Percent inductance error (in %)", fontsize=24)
    # # ax.plot(sweep_wndg_pos, error, 'o', c='#%02x%02x%02x' % fmt.colors_femmt_default['red'], linewidth=4)
    # # ax.plot(sweep_wndg_pos, error, linewidth=4)
    # for j in range(len(sweep_air_gap_h)):
    #     # ax.plot(sweep_wndg_pos, fem_ind[j, :], sweep_wndg_pos, cal_ind[j, :], label=str(sweep_air_gap_h[j]))
    #     ax.plot(sweep_wndg_pos, error[j, :], label=f"{str(round(sweep_air_gap_h[j], 5))} m", linewidth=4)
    # ax.legend(loc='best', fontsize=18, title="air-gap length", title_fontsize=18)
    # ax.grid()
    # fmt.plt.show()

