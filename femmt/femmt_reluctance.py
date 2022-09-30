import femmt as fmt
import numpy as np
from itertools import product


def plot_r_basis():
    # width = 1
    # length = 1
    # height = np.linspace(10, 0.1, 1000)
    width = 0.0149
    length = 0.0005
    height = np.linspace(0.005, 0, 1000)
    h_l = height / length

    r_m = 1 / (fmt.mu0 * (width / 2 / length + 2 / np.pi * (
                1 + np.log(np.pi * height / 4 / length))))

    combined = np.vstack((h_l, r_m)).T
    print(combined)
    fig, ax = fmt.plt.subplots()  # Create a figure containing a single axes.
    fmt.plt.title("R_basic vs h/l")
    fmt.plt.xlabel("h/l")
    fmt.plt.ylabel("R_basic")
    ax.plot(h_l, r_m)
    ax.invert_xaxis()
    ax.grid()
    fmt.plt.show()


class MagneticCircuit:
    """This is a class for calculating the reluctance and inductance for 2D axis symmetric components
        and visualising magnetic circuit"""

    def __init__(self, core_w: list, window_h: list, window_w: list, no_of_turns: list, n_air_gaps: list,
                 air_gap_h: list, air_gap_position: list, mu_rel: list, mult_air_gap_type: list = None,
                 air_gap_method='percent', component_type='inductor', sim_type='single'):
        """
        :param core_w: Diameter of center leg of the core in meter
        :type core_w: list
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

        self.row_num = 0
        self.single_air_gap_len = None
        self.data_matrix_len = None
        self.data_matrix = None

        if not all(isinstance(item, int) for item in no_of_turns):
            raise Exception("no_of_turns list elements should be integer")
        if not all(isinstance(item, int) for item in n_air_gaps):
            raise Exception("n_air_gaps list elements should be integer")
        if not (air_gap_method == 'center' or air_gap_method == 'percent' or air_gap_method == 'manual'):
            raise Exception("string value wrong for air_gap_method argument")
        if not (sim_type == 'single' or sim_type == 'sweep'):
            raise Exception("string value wrong for sim_type argument")
        if not (component_type == 'inductor' or component_type == 'integrated_transformer'):
            raise Exception("string value wrong for component_type argument")
        if any(item > 0.0005 for item in air_gap_h):
            raise Exception("Model accuracy is not good for air_gap_h more than 0.0005")
        if sim_type == 'single':
            if not (len(core_w) == 1 and len(window_h) == 1 and len(window_w) == 1 and len(no_of_turns) == 1
                    and len(n_air_gaps) == 1 and len(mu_rel) == 1):
                raise Exception("single sim_type requires single list elements")
            if not (n_air_gaps[0] == len(air_gap_h) and n_air_gaps[0] == len(air_gap_position)):
                raise Exception("No. of elements of air_gap_h and air_gap_position should match n_air_gaps")

        # Sort air_gap_position and air_gap_h based on air_gap_position
        zipped_lists = zip(air_gap_position, air_gap_h)
        sorted_pairs = sorted(zipped_lists)

        tuples = zip(*sorted_pairs)
        self.air_gap_position, self.air_gap_h = [list(tpl) for tpl in tuples]

        self.core_w = np.array(core_w)
        self.window_h = np.array(window_h)
        self.window_w = np.array(window_w)
        self.mu_rel = np.array(mu_rel)
        self.no_of_turns = np.array(no_of_turns)
        self.n_air_gaps = np.array(n_air_gaps)
        self.air_gap_h = np.array(self.air_gap_h)
        self.percent_position_air_gap = np.array(self.air_gap_position)
        self.mult_air_gap_type = None
        self.sim_type = sim_type
        self.data_matrix = np.zeros((1, 10))

        if not len(n_air_gaps):
            self.n_air_gaps = np.zeros((1, 1))

        # Creates the data matrix with all the input parameter combinations for sim_type = 'sweep'
        if sim_type == 'sweep':
            self.create_data_matrix(core_w, window_h, window_w, no_of_turns, n_air_gaps, air_gap_h,
                                    air_gap_position, mu_rel, mult_air_gap_type)

            self.core_w = self.data_matrix[:, 0]        # Diameter of center leg
            self.window_h = self.data_matrix[:, 1]
            self.window_w = self.data_matrix[:, 2]
            self.mu_rel = self.data_matrix[:, 3]        # 3000
            self.no_of_turns = self.data_matrix[:, 4]
            self.n_air_gaps = self.data_matrix[:, 5]
            self.air_gap_h = self.data_matrix[:, 6]
            self.percent_position_air_gap = self.data_matrix[:, 7]
            self.mult_air_gap_type = self.data_matrix[:, 8]

        self.core_h = self.window_h + self.core_w / 2
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

    def create_data_matrix(self, core_w: list, window_h: list, window_w: list, no_of_turns: list, n_air_gaps: list,
                           air_gap_h: list, air_gap_position: list, mu_rel: list, mult_air_gap_type: list):

        """ Creates matrix consisting of input design parameters with all their combinations

        :param core_w: Diameter of center leg of the core in meter
        :type core_w: list
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

        # example: data_matrix = [core_w, window_h, window_w, mu_rel, no_of_turns, n_air_gaps, air_gap_h,
        #                      air_gap_position, mult_air_gap_type, inductance]
        clone_n_air_gaps = n_air_gaps

        if 1 in clone_n_air_gaps:
            self.data_matrix = np.zeros((len(core_w) * len(no_of_turns) * len(air_gap_h) * len(mu_rel) * (
                    len(air_gap_position) + (len(n_air_gaps) - 1) * len(mult_air_gap_type)), 10))
        else:
            self.data_matrix = np.zeros((len(core_w) * len(no_of_turns) * len(air_gap_h) * len(n_air_gaps) * len(
                mu_rel) * len(mult_air_gap_type), 10))

        if 1 in clone_n_air_gaps:
            clone_n_air_gaps.remove(1)
            for index_1 in range(len(core_w)):
                for index_2 in range(len(mu_rel)):
                    for index_3 in range(len(no_of_turns)):
                        for index_4 in range(len(air_gap_h)):
                            for index_5 in range(len(air_gap_position)):
                                self.data_matrix[self.row_num, 0] = core_w[index_1]
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
            for index_1 in range(len(core_w)):
                for index_2 in range(len(mu_rel)):
                    for index_3 in range(len(no_of_turns)):
                        for index_4 in range(len(clone_n_air_gaps)):
                            for index_5 in range(len(air_gap_h)):
                                for index_6 in range(len(mult_air_gap_type)):
                                    self.data_matrix[self.row_num, 0] = core_w[index_1]
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
        self.r_inner = self.core_w / 2 + self.window_w
        self.r_outer = np.sqrt((self.core_w / 2) ** 2 + self.r_inner ** 2)
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

        self.area[:, 0] = np.pi * ((self.core_w / 2) ** 2)

        self.length[:, 1] = (np.pi / 8) * (self.core_w / 2 + self.core_h_middle)
        # self.area[:, 1] = ((self.core_w / 2 + self.core_h_middle) / 2) * 2 * np.pi * (self.core_w / 2)
        self.area[:, 1] = np.pi * self.core_w / 2 * ((self.core_w / 2 + self.core_h_middle) / 2) * (
                2 + ((self.core_w / 2 + self.core_h_middle) / 2) / (
            np.sqrt((self.core_w / 2) ** 2 + self.core_h_middle ** 2)))

        self.length[:, 2] = self.window_w
        self.area[:, 2] = np.nan
        self.reluctance[:, 2] = ((self.mu_0 * self.mu_rel * 2 * np.pi * self.core_h_middle) ** -1) * np.log(
            (2 * self.r_inner) / self.core_w)

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
                                                                     self.core_w[0:self.single_air_gap_len], h[:, 0]),
                                                                 single_round_round(
                                                                     self.air_gap_h[0:self.single_air_gap_len],
                                                                     self.core_w[0:self.single_air_gap_len], h[:, 0],
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
                               self.core_w[self.single_air_gap_len:self.data_matrix_len],
                               self.n_air_gaps[self.single_air_gap_len:self.data_matrix_len], h_multiple),
            distributed_type_2(self.air_gap_h[self.single_air_gap_len:self.data_matrix_len],
                               self.core_w[self.single_air_gap_len:self.data_matrix_len],
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

    def get_param_pos_dict(self):
        self.param_pos_dict = {"core_w": 0, "window_h": 1, "window_w": 2, "mu_rel": 3, "no_of_turns": 4,
                               "n_air_gaps": 5,
                               "air_gap_h": 6, "air_gap_position": 7, "mult_air_gap_type": 8, "inductance": 9,
                               "core_h_middle": 10,
                               "r_inner": 11, "r_outer": 12, "center_leg_area": 13, "outer_leg_area": 14}

        return self.param_pos_dict

    def air_gap_reluctance_single(self):
        flag_0 = 0
        flag_1 = 0
        flag_2 = 0
        if self.n_air_gaps[0] != 0:
            if self.air_gap_method == 'center':
                self.section.append(6)  # round-round type airgap
                temp1 = fmt.r_basis([self.air_gap_h[0] / 2], [self.core_w], [(self.window_h - self.air_gap_h[0]) / 2])
                temp2 = fmt.sigma([self.air_gap_h[0]], [self.core_w / 2], 2 * temp1)
                temp3 = fmt.r_round_round([self.air_gap_h[0]], temp2, [self.core_w / 2])
                temp4 = self.air_gap_h[0] / (self.mu_0 * np.pi * (self.core_w / 2) ** 2)  # classical reluctance formula
                self.reluctance[:, 5] = temp3

            elif self.air_gap_method == 'percent' or self.air_gap_method == 'manual':
                self.max_percent_position = ((self.window_h - self.air_gap_h[self.n_air_gaps - 1] / 2) / self.window_h) * 100
                self.min_percent_position = ((self.air_gap_h[0] / 2) / self.window_h) * 100
                if self.air_gap_method == 'percent':
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

                    temp1 = fmt.r_basis([self.air_gap_h[0]], [self.core_w], [h])
                    temp2 = fmt.sigma([self.air_gap_h[0]], [self.core_w / 2], temp1)
                    temp3 = fmt.r_round_inf([self.air_gap_h[0]], temp2, [self.core_w / 2])
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

                    temp1 = fmt.r_basis([self.air_gap_h[self.n_air_gaps - 1]], [self.core_w], [h])
                    temp2 = fmt.sigma([self.air_gap_h[self.n_air_gaps - 1]], [self.core_w / 2], temp1)
                    temp3 = fmt.femmt_functions.r_round_inf([self.air_gap_h[self.n_air_gaps - 1]], temp2, [self.core_w / 2])
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

                        r_basis_1 = fmt.r_basis([self.air_gap_h[i] / 2], [self.core_w], [h1])
                        r_basis_2 = fmt.r_basis([self.air_gap_h[i] / 2], [self.core_w], [h2])
                        temp2 = fmt.sigma([self.air_gap_h[i]], [self.core_w / 2], r_basis_1 + r_basis_2)
                        temp3 = fmt.femmt_functions.r_round_round([self.air_gap_h[i]], temp2, [self.core_w / 2])
                        self.reluctance[:, 5] = self.reluctance[:, 5] + temp3

    def calculate_inductance(self):
        # self.section, self.orientation = set_orientation(self.section, len(self.section))
        self.cal_inductance = (self.no_of_turns * self.no_of_turns) / np.sum(self.reluctance, axis=1)
        self.data_matrix[:, 9] = self.cal_inductance


def single_round_inf(air_gap_h, core_w, h):
    """Returns reluctance of a single air-gap at the corner

    :param air_gap_h: Air-gap height [in meter]
    :type air_gap_h: list
    :param core_w: Diameter of center leg of the core [in meter]
    :type core_w: list
    :param h: Core distance between air-gap and other end of the window-h [in meter]
    :type h: list
    :return: Reluctance of a single air-gap at the corner
    :rtype: list"""

    temp1 = fmt.r_basis(air_gap_h, core_w, h)
    temp2 = fmt.sigma(air_gap_h, core_w / 2, temp1)
    temp3 = fmt.r_round_inf(air_gap_h, temp2, core_w / 2)

    return temp3


def single_round_round(air_gap_h, core_w, h0, h1):
    """Returns reluctance of a single air-gap at position other than corner on the center leg

    :param air_gap_h: Air-gap height [in meter]
    :type air_gap_h: list
    :param core_w: Diameter of center leg of the core [in meter]
    :type core_w: list
    :param h0: Distance between window_h and air_gap_position for a single air-gap [in meter]
    :type h0: list
    :param h1: Height of air-gap from the base of the core window [in meter]
    :type h1: list
    :return: Reluctance of a single air-gap at position other than corner on the center leg
    :rtype: list"""

    r_basis_1 = fmt.r_basis(air_gap_h / 2, core_w, h0)
    r_basis_2 = fmt.r_basis(air_gap_h / 2, core_w, h1)
    temp2 = fmt.sigma(air_gap_h, core_w / 2, r_basis_1 + r_basis_2)
    temp3 = fmt.r_round_round(air_gap_h, temp2, core_w / 2)

    return temp3


def distributed_type_1(air_gap_h, core_w, n_air_gaps, h_multiple):
    """Returns distributed air-gap reluctance of Type 1 (Where corner air-gaps are present)

    :param air_gap_h: Air-gap height [in meter]
    :type air_gap_h: list
    :param core_w: Diameter of center leg of the core [in meter]
    :type core_w: list
    :param n_air_gaps: Number of air-gaps in the center leg of the core
    :type n_air_gaps: list
    :param h_multiple: Half of core height between two consecutive air-gaps in an equally distributed air-gaps [in meter]
    :type h_multiple: list
    :return: Distributed air-gap reluctance of Type 1 (Where corner air-gaps are present)
    :rtype: list"""

    temp1 = fmt.r_basis(air_gap_h, core_w, h_multiple)
    temp2 = fmt.sigma(air_gap_h, core_w / 2, temp1)
    temp3 = fmt.r_round_inf(air_gap_h, temp2, core_w / 2)
    reluctance = (2 * temp3)

    r_basis_1 = fmt.r_basis(air_gap_h / 2, core_w, h_multiple)
    r_basis_2 = fmt.r_basis(air_gap_h / 2, core_w, h_multiple)
    temp2 = fmt.sigma(air_gap_h, core_w / 2, r_basis_1 + r_basis_2)
    temp3 = fmt.r_round_round(air_gap_h, temp2, core_w / 2)
    reluctance = reluctance + ((n_air_gaps - 2) * temp3)

    return reluctance


def distributed_type_2(air_gap_h, core_w, n_air_gaps, h_multiple):
    """Returns distributed air-gap reluctance of Type 2 (Where corner air-gaps are absent)

    :param air_gap_h: Air-gap height [in meter]
    :type air_gap_h: list
    :param core_w: Diameter of center leg of the core [in meter]
    :type core_w: list
    :param n_air_gaps: Number of air-gaps in the center leg of the core
    :type n_air_gaps: list
    :param h_multiple: Core height between two consecutive air-gaps in an equally distributed air-gaps [in meter]
    :type h_multiple: list
    :return: Distributed air-gap reluctance of Type 2 (Where corner air-gaps are absent)
    :rtype: list"""

    r_basis_1 = fmt.r_basis(air_gap_h / 2, core_w, h_multiple)
    r_basis_2 = fmt.r_basis(air_gap_h / 2, core_w, h_multiple / 2)
    temp2 = fmt.sigma(air_gap_h, core_w / 2, r_basis_1 + r_basis_2)
    temp3 = fmt.r_round_round(air_gap_h, temp2, core_w / 2)
    reluctance = (2 * temp3)

    r_basis_1 = fmt.r_basis(air_gap_h / 2, core_w, h_multiple / 2)
    r_basis_2 = fmt.r_basis(air_gap_h / 2, core_w, h_multiple / 2)
    temp2 = fmt.sigma(air_gap_h, core_w / 2, r_basis_1 + r_basis_2)
    temp3 = fmt.r_round_round(air_gap_h, temp2, core_w / 2)
    reluctance = reluctance + ((n_air_gaps - 2) * temp3)

    return reluctance


if __name__ == '__main__':
    mc1 = MagneticCircuit(core_w=[0.0149], window_h=[0.0295], window_w=[0.01105], no_of_turns=[8], n_air_gaps=[3],
                          air_gap_h=[0.0005, 0.0002, 0.0001], air_gap_position=[0, 50, 100], mu_rel=[3000], mult_air_gap_type=[1, 2],
                          air_gap_method='percent', component_type='inductor', sim_type='single')  # 0.0149

    mc1.calculate_inductance()

    print(f"Inductance is {mc1.cal_inductance}")
    # plot_r_basis()

# def basic_example_func(f_height, f_position, f_n_turns, f_core_cond_iso, f_current):
#     # 2. set core parameters
#     core_db = fmt.core_database()["PQ 40/40"]
#
#     core = fmt.Core(core_w=core_db["core_w"], window_w=core_db["window_w"], window_h=core_db["window_h"],
#                     material="95_100")
#     # mu_rel=3000, phi_mu_deg=10,
#     # sigma=0.5)
#     geo.set_core(core)
#
#     # 3. set air gap parameters
#     air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
#     air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, f_position[0], f_height)
#     # air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, f_position[1], f_height)
#     # air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, f_position[2], f_height)
#     # air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, f_position[3], f_height)
#     # air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, f_position[4], f_height)
#     geo.set_air_gaps(air_gaps)
#
#     # 4. set conductor parameters: use solid wires
#     winding = fmt.Winding(f_n_turns, 0, fmt.Conductivity.Copper, fmt.WindingType.Primary, fmt.WindingScheme.Square)
#     winding.set_solid_conductor(0.0013)
#     # winding.set_litz_conductor(conductor_radius=0.0013, number_strands=150, strand_radius=100e-6, fill_factor=None)
#     geo.set_windings([winding])
#
#     # 5. set isolations
#     isolation = fmt.Isolation()
#     isolation.add_core_isolations(0.001, 0.001, f_core_cond_iso, 0.001)
#     isolation.add_winding_isolations(0.0005)
#     geo.set_isolation(isolation)
#
#     # 5. create the model
#     geo.create_model(freq=100000, visualize_before=False, save_png=False)
#
#     # 6. start simulation
#     geo.single_simulation(freq=100000, current=[f_current], show_results=False)
#
#
# # Sweep of air-gap and winding position and compare it with FEM simulation
# sweep_air_gap_h = np.linspace(0.0001, 0.0005, 3)
# sweep_wndg_pos = np.linspace(0.001, 0.007, 5)
# # sweep_current = np.linspace(0.1, 10, 5)
# fem_ind = np.zeros((len(sweep_air_gap_h), len(sweep_wndg_pos)))
# cal_ind = np.zeros((len(sweep_air_gap_h), len(sweep_wndg_pos)))
# # cal_fringe_dist = np.zeros((len(sweep_air_gap_h), len(sweep_wndg_pos)))
#
# # Working directory can be set arbitrarily
# working_directory = fmt.os.path.join(fmt.os.path.dirname(__file__), '..')
#
# # 1. chose simulation type
# geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=working_directory)
#
# for j in range(len(sweep_air_gap_h)):
#     for k in range(len(sweep_wndg_pos)):
#         mc1 = MagneticCircuit([0.0149], [0.0295], [0.01105], [8], [1], [0.0001, 0.0005])  # 0.0149
#         mc1.core_reluctance()
#         mc1.air_gap_reluctance()
#         cal_ind[j, k] = mc1.L  # - (2.38295 * 1e-6 - 0.000326175 * sweep_wndg_pos[k])
#         # cal_fringe_dist[j, k] = mc1.fringe_dist
#         mc1.max_percent = ((mc1.window_h - (sweep_air_gap_h[j] / 2)) / mc1.window_h) * 100
#         mc1.min_percent = ((sweep_air_gap_h[j] / 2) / mc1.window_h) * 100
#
#         # basic_example_func(sweep_air_gap_h[j], [mc1.min_percent + 0.1, 25, 50, 75, mc1.max_percent - 0.1], 9, sweep_wndg_pos[k], 10)
#         # basic_example_func(sweep_air_gap_h[j], [50], 9, sweep_wndg_pos[k], 10)
#         fem_ind[j, k] = geo.read_log()["single_sweeps"][0]["winding1"]["self_inductivity"][0]
#         # fem_ind[j, k] = geo.read_log()["single_sweeps"][0]["winding1"]["Q"]
#
# print(f"Air-gap length: {sweep_air_gap_h}")
# print(f"Winding position: {sweep_wndg_pos}")
# print(f"FEM inductance: {fem_ind}")
# print(f"Calculated inductance: {cal_ind}")
# # print(f"Calculated fringe dist: {cal_fringe_dist}")
#
# h_by_l = ((0.0295 - sweep_air_gap_h) / 2) / sweep_air_gap_h
# abs_error = cal_ind - fem_ind
# error = (abs_error / fem_ind) * 100
# avg_abs_error = np.sum(abs_error, axis=0) / 5
#
# print(f"h_by_l: {h_by_l}")
# print(f"abs_error: {abs_error}")
# print(f"percent error: {error}")
# print(f"average actual error: {avg_abs_error}")
#
# np.savetxt('absolute_error.txt', abs_error)
#
# # Plotting tools
# fig, ax = fmt.plt.subplots()  # Create a figure containing a single axes.
# fmt.plt.title("Inductance % error vs winding pos (3 conductors) (center airgaps)")
# fmt.plt.xlabel("Winding position (in m)")
# fmt.plt.ylabel("Percent error (in %)")
# # ax.plot(Air_gap_length, FEM_inductance, 'o')
# for j in range(len(sweep_air_gap_h)):
#     # ax.plot(sweep_wndg_pos, fem_ind[j, :], sweep_wndg_pos, cal_ind[j, :], label=str(sweep_air_gap_h[j]))
#     ax.plot(sweep_wndg_pos, error[j, :], label=str(sweep_air_gap_h[j]))
# ax.legend(loc='best')
# ax.grid()
# fmt.plt.show()

# mc1 = MagneticCircuit(0.0398, 0.0149, 0.0295, 0.01105, 8, 3, 'center', 1, [0.0005], [50])
# mc1 = MagneticCircuit(0.0149, 0.0295, 0.01105, 8, 3, 'center', 1, [0.02], [50])
# mc1 = MagneticCircuit(0.0149, 0.0295, 0.01105, 8, 3, 'percent', 2, [0.0005, 0.0005], [20, 50])
# mc1.core_reluctance()
# mc1.air_gap_reluctance()
# mc1.draw_schematic()
# plot_error()

# import schemdraw
# import schemdraw.elements as elm

# def set_orientation(section: list, n: int) -> tuple:    # Function to dynamically define circuit component orientation: 'up', 'right', 'down', 'left'
#     temp1 = ['up', 'right', 'down', 'left']
#     temp2 = []
#     if n % 2 != 0:
#         section.append('l')
#         n += 1
#     if n % 2 == 0:
#         if n % 4 == 0:
#             for i in range(4):
#                 for j in range(n // 4):
#                     temp2.append(temp1[i])
#         else:
#             for i in range(4):
#                 for j in range((n-2) // 4):
#                     temp2.append(temp1[i])
#
#             temp2.insert((n-2) // 4, 'right')
#             temp2.append('left')
#     return section, temp2


# def plot_r_basis():
#     width = 1
#     length = 1
#     height = np.linspace(10, 0.1, 1000)
#     h_l = height / length
#
#     r_m = 1 / (mu0 * (width / 2 / length + 2 / np.pi * (
#                 1 + np.log(np.pi * height / 4 / length))))
#
#     combined = np.vstack((h_l, r_m)).T
#     print(combined)
#     fig, ax = plt.subplots()  # Create a figure containing a single axes.
#     plt.title("R_basic vs h/l")
#     plt.xlabel("h/l")
#     plt.ylabel("R_basic")
#     ax.plot(h_l, r_m)
#     ax.invert_xaxis()
#     ax.grid()
#     plt.show()


# def plot_error():
#     fem_ind = np.array([116.0773e-6, 32.1e-6, 18.71e-6, 6.08575e-6, 4.324969e-6,  3.81392e-6, 2.345389e-6])
#     cal_ind = np.array([116.518e-6, 32.12472e-6, 18.72275e-6, 6.25557e-6, 4.2895602e-6, 3.611239e-6, 0.71733329e-6])
#     air_gap_l = np.array([0.0001, 0.0005, 0.001, 0.005, 0.00833, 0.01, 0.02])
#     h_by_l = ((0.0295 - air_gap_l) / 2) / air_gap_l
#     error = ((fem_ind - cal_ind) / fem_ind) * 100
#     fig, ax = plt.subplots()  # Create a figure containing a single axes.
#     plt.title("inductance vs h/l")
#     plt.xlabel("h/l")
#     plt.ylabel("Error in %")
#     # ax.plot(h_by_l, fem_ind, h_by_l, cal_ind)
#     ax.plot(h_by_l, error)
#     # ax.invert_xaxis()
#     ax.grid()
#     plt.show()


# def draw_schematic(self):

# print('Section:', self.section)
# print('Orientation:', self.orientation)
# print('List of length:', self.length)
# print('List of cross-section area:', self.area)
# print('List of reluctance:', self.reluctance)
# print('Inductance:', self.L)

# schemdraw.theme('monokai')
# d = schemdraw.Drawing()
# for i in range(len(self.section)):
#     if self.section[i] == 9:
#         d += getattr(elm.SourceV().dot().label(str(round(self.reluctance[i], 2)) + ' AT'), self.orientation[i])
#     elif self.section[i] == 'l':
#         d += getattr(elm.Line().dot(), self.orientation[i])
#     elif self.section[i] == 6 or self.section[i] == 7 or self.section[i] == 8:
#         d += getattr(elm.ResistorIEC().dot().label(str(round(self.reluctance[i], 2)) + ' AT/Wb'), self.orientation[i])
#     else:
#         d += getattr(elm.Resistor().dot().label(str(round(self.reluctance[i], 2)) + ' AT/Wb'), self.orientation[i])
#
# d.draw()
# d.save('my_circuit.svg')

# Air_gap_length = [1.00000000e-07, 4.08261224e-04, 8.16422449e-04, 1.22458367e-03,
#  1.63274490e-03, 2.04090612e-03, 2.44906735e-03, 2.85722857e-03,
#  3.26538980e-03, 3.67355102e-03, 4.08171224e-03, 4.48987347e-03,
#  4.89803469e-03, 5.30619592e-03, 5.71435714e-03, 6.12251837e-03,
#  6.53067959e-03, 6.93884082e-03, 7.34700204e-03, 7.75516327e-03,
#  8.16332449e-03, 8.57148571e-03, 8.97964694e-03, 9.38780816e-03,
#  9.79596939e-03, 1.02041306e-02, 1.06122918e-02, 1.10204531e-02,
#  1.14286143e-02, 1.18367755e-02, 1.22449367e-02, 1.26530980e-02,
#  1.30612592e-02, 1.34694204e-02, 1.38775816e-02, 1.42857429e-02,
#  1.46939041e-02, 1.51020653e-02, 1.55102265e-02, 1.59183878e-02,
#  1.63265490e-02, 1.67347102e-02, 1.71428714e-02, 1.75510327e-02,
#  1.79591939e-02, 1.83673551e-02, 1.87755163e-02, 1.91836776e-02,
#  1.95918388e-02, 2.00000000e-02]
#
# FEM_inductance = [[5.25733873e-04],
#  [3.56509342e-05],
#  [1.96273371e-05],
#  [1.38607988e-05],
#  [1.08431919e-05],
#  [8.96908451e-06],
#  [7.68081665e-06],
#  [6.73635898e-06],
#  [6.00979154e-06],
#  [5.42909853e-06],
#  [4.95256678e-06],
#  [4.55325727e-06],
#  [4.21369509e-06],
#  [3.92151073e-06],
#  [3.66711744e-06],
#  [3.44366581e-06],
#  [3.24587203e-06],
#  [3.06932936e-06],
#  [2.91034587e-06],
#  [2.76613127e-06],
#  [2.63475469e-06],
#  [2.51460888e-06],
#  [2.40401003e-06],
#  [2.30159406e-06],
#  [2.20646338e-06],
#  [2.11756703e-06],
#  [2.03425524e-06],
#  [1.95578553e-06],
#  [1.88179771e-06],
#  [1.81227984e-06],
#  [1.74712510e-06],
#  [1.68603210e-06],
#  [1.62871517e-06],
#  [1.57488291e-06],
#  [1.52419216e-06],
#  [1.47629414e-06],
#  [1.43075069e-06],
#  [1.38734476e-06],
#  [1.34601111e-06],
#  [1.30665534e-06],
#  [1.26911897e-06],
#  [1.23328032e-06],
#  [1.19909825e-06],
#  [1.16647940e-06],
#  [1.13532254e-06],
#  [1.10558445e-06],
#  [1.07726437e-06],
#  [1.05039705e-06],
#  [1.02504119e-06],
#  [1.00107838e-06]]

# def air_gap_reluctance(self):
#     flag_0 = 0
#     flag_1 = 0
#     flag_2 = 0
#
#     if self.method == 'center':
#         self.section.append(6)  # round-round
#         temp1 = fmt.r_basis(self.air_gap_h[0] / 2, self.core_w, (self.window_h - self.air_gap_h[0]) / 2)
#         temp2 = fmt.sigma(self.air_gap_h[0], self.core_w / 2, 2 * temp1)
#         temp3 = fmt.r_round_round(self.air_gap_h[0], temp2, self.core_w / 2)
#         temp4 = self.air_gap_h[0] / (self.mu_0 * np.pi * (self.core_w / 2) ** 2)  # classical reluctance formula
#         self.reluctance = np.append(self.reluctance, temp4)
#         self.fringe_area = (np.pi * (self.core_w / 2) ** 2) * (1 / (temp2 * temp2))
#         self.fringe_dist = np.sqrt(self.fringe_area / np.pi) - (self.core_w / 2)
#     # Assuming only equally distributed airgaps
#     elif self.method == 'percent':
#         self.max_percent = ((self.window_h - self.air_gap_h[self.n_air_gaps - 1] / 2) / self.window_h) * 100
#         self.min_percent = ((self.air_gap_h[0] / 2) / self.window_h) * 100
#         self.position = np.array(
#             self.air_gap_position) / 100 * self.window_h  # Convert percent position to absolute value position
#         print(f"Max percent: {self.max_percent}")
#         print(f"Min percent: {self.min_percent}")
#
#         if self.air_gap_position[0] <= self.min_percent:
#             flag_0 = 1
#             self.section.append(8)
#             if self.n_air_gaps == 1:
#                 h = self.window_h - self.air_gap_h[0]
#             else:
#                 h = ((self.position[1] - self.air_gap_h[1] / 2) - self.air_gap_h[0]) / 2
#
#             temp1 = fmt.r_basis(self.air_gap_h[0], self.core_w, h)
#             temp2 = fmt.sigma(self.air_gap_h[0], self.core_w / 2, temp1)
#             temp3 = fmt.r_round_inf(self.air_gap_h[0], temp2, self.core_w / 2)
#             self.reluctance = np.append(self.reluctance, temp3)
#             print('air gap is at lower corner')
#
#         if self.air_gap_position[self.n_air_gaps - 1] >= self.max_percent:
#             flag_1 = 1
#             self.section.append(8)
#             if self.n_air_gaps == 1:
#                 h = self.window_h - self.air_gap_h[self.n_air_gaps - 1]
#             else:
#                 h = (self.position[self.n_air_gaps - 1] - self.position[self.n_air_gaps - 2] - self.air_gap_h[
#                     self.n_air_gaps - 1] / 2 - self.air_gap_h[self.n_air_gaps - 2] / 2) / 2
#
#             temp1 = fmt.r_basis(self.air_gap_h[self.n_air_gaps - 1], self.core_w, h)
#             temp2 = fmt.sigma(self.air_gap_h[self.n_air_gaps - 1], self.core_w / 2, temp1)
#             temp3 = fmt.femmt_functions.r_round_inf(self.air_gap_h[self.n_air_gaps - 1], temp2, self.core_w / 2)
#             self.reluctance = np.append(self.reluctance, temp3)
#             print('air gap is at upper corner')
#
#         for i in range(self.n_air_gaps):
#             if self.min_percent < self.air_gap_position[i] < self.max_percent:
#                 self.section.append(7)
#                 if flag_2 == 0:
#
#                     if flag_0 == 0 and flag_1 == 0:  # No corner air-gaps
#                         self.position = np.append(self.position,
#                                                   self.window_h + (self.window_h - self.position[self.n_air_gaps - 1]))
#                         self.position = np.insert(self.position, 0, -self.position[0])
#                         self.air_gap_h = np.append(self.air_gap_h, self.air_gap_h[self.n_air_gaps - 1])
#                         self.air_gap_h = np.insert(self.air_gap_h, 0, self.air_gap_h[0])
#                     elif flag_0 == 1 and flag_1 == 0:  # Only lower air-gap is present
#                         self.position = np.append(self.position,
#                                                   self.window_h + (self.window_h - self.position[self.n_air_gaps - 1]))
#                         self.air_gap_h = np.append(self.air_gap_h, self.air_gap_h[self.n_air_gaps - 1])
#                     elif flag_0 == 0 and flag_1 == 1:  # Only Upper air-gap is present
#                         self.position = np.insert(self.position, 0, -self.position[0])
#                         self.air_gap_h = np.insert(self.air_gap_h, 0, self.air_gap_h[0])
#                     flag_2 = 1
#
#                 if flag_0 == 0 and flag_1 == 0:
#                     h1 = (self.position[i + 1] - self.position[i] - self.air_gap_h[i + 1] / 2 - self.air_gap_h[
#                         i] / 2) / 2
#                     h2 = (self.position[i + 2] - self.position[i + 1] - self.air_gap_h[i + 2] / 2 - self.air_gap_h[
#                         i + 1] / 2) / 2
#                     print('No corner air gap detected')
#                 elif flag_0 == 1 and flag_1 == 0:
#                     h1 = (self.position[i] - self.position[i - 1] - self.air_gap_h[i] / 2 - self.air_gap_h[
#                         i - 1] / 2) / 2
#                     h2 = (self.position[i + 1] - self.position[i] - self.air_gap_h[i + 1] / 2 - self.air_gap_h[
#                         i] / 2) / 2
#                     print('Lower air gap detected')
#                 elif flag_0 == 0 and flag_1 == 1:
#                     h1 = (self.position[i + 1] - self.position[i] - self.air_gap_h[i + 1] / 2 - self.air_gap_h[
#                         i] / 2) / 2
#                     h2 = (self.position[i + 2] - self.position[i + 1] - self.air_gap_h[i + 2] / 2 - self.air_gap_h[
#                         i + 1] / 2) / 2
#                     print('Upper air gap detected')
#                 else:
#                     h1 = (self.position[i] - self.position[i - 1] - self.air_gap_h[i] / 2 - self.air_gap_h[
#                         i - 1] / 2) / 2
#                     h2 = (self.position[i + 1] - self.position[i] - self.air_gap_h[i + 1] / 2 - self.air_gap_h[
#                         i] / 2) / 2
#                     print('Both air gap detected')
#
#                 r_basis_1 = fmt.r_basis(self.air_gap_h[i] / 2, self.core_w, h1)
#                 r_basis_2 = fmt.r_basis(self.air_gap_h[i] / 2, self.core_w, h2)
#                 temp2 = fmt.sigma(self.air_gap_h[i], self.core_w / 2, r_basis_1 + r_basis_2)
#                 temp3 = fmt.femmt_functions.r_round_round(self.air_gap_h[i], temp2, self.core_w / 2)
#                 self.reluctance = np.append(self.reluctance, temp3)
#
#     # self.section, self.orientation = set_orientation(self.section, len(self.section))
#     self.L = (self.no_of_turns * self.no_of_turns) / sum(self.reluctance)

# self.data_matrix = np.zeros(
#             (len(core_w) * len(no_of_turns) * len(n_air_gaps) * len(air_gap_h) * len(air_gap_position), 8))
#         row_num = 0
#         for index_1 in range(len(core_w)):
#             for index_2 in range(len(no_of_turns)):
#                 for index_3 in range(len(n_air_gaps)):
#                     for index_4 in range(len(air_gap_h)):
#                         for index_5 in range(len(air_gap_position)):
#                             self.data_matrix[row_num, 0] = core_w[index_1]
#                             self.data_matrix[row_num, 1] = window_h[index_1]
#                             self.data_matrix[row_num, 2] = window_w[index_1]
#                             self.data_matrix[row_num, 3] = no_of_turns[index_2]
#                             self.data_matrix[row_num, 4] = n_air_gaps[index_3]
#                             self.data_matrix[row_num, 5] = air_gap_h[index_4]
#                             self.data_matrix[row_num, 6] = air_gap_position[index_5]
#                             row_num = row_num + 1

# def air_gap_reluctance(self):
#      # self.section.append(6)  # acknowledge addition of air-gap
#      self.max_percent_position = ((self.window_h - (self.air_gap_h / 2)) / self.window_h) * 100
#      self.min_percent_position = ((self.air_gap_h / 2) / self.window_h) * 100
#      self.abs_position_air_gap = (self.percent_position_air_gap * self.window_h) / 100  # Convert percent position to absolute value position
#      h = np.zeros((len(self.data_matrix), 2))
#
#      h[:, 0] = np.where((self.percent_position_air_gap <= self.min_percent_position) | (
#              self.percent_position_air_gap >= self.max_percent_position),
#                         self.window_h - self.air_gap_h,
#                         self.abs_position_air_gap - (self.air_gap_h / 2))
#      h[:, 1] = self.window_h - (self.abs_position_air_gap - (self.air_gap_h / 2))
#
#      self.reluctance[:, 5] = np.where(h[:, 1] == 0, fmt.r_round_inf(self.air_gap_h, fmt.sigma(self.air_gap_h,
#                                                                                               self.core_w / 2,
#                                                                                               fmt.r_basis(
#                                                                                                   self.air_gap_h,
#                                                                                                   self.core_w,
#                                                                                                   h[:, 0])),
#                                                                     self.core_w / 2),
#                                       fmt.femmt_functions.r_round_round(self.air_gap_h, fmt.sigma(self.air_gap_h,
#                                                                                                   self.core_w / 2,
#                                                                                                   (fmt.r_basis(
#                                                                                                       self.air_gap_h / 2,
#                                                                                                       self.core_w,
#                                                                                                       h[:, 0])) +
#                                                                                                   (fmt.r_basis(
#                                                                                                       self.air_gap_h / 2,
#                                                                                                       self.core_w,
#                                                                                                       h[:, 1]))),
#                                                                         self.core_w / 2))
#
#      self.data_matrix[:, 7] = (self.no_of_turns ** 2) / np.sum(self.reluctance, axis=1)
#      print(self.data_matrix)
# for i in range(self.single_air_gap_len, self.row_num):
#     if self.mult_air_gap_type[i] == 1:
#         h_multiple = (self.window_h[i] - (self.n_air_gaps[i] * self.air_gap_h[i])) / ((self.n_air_gaps[i] - 1) * 2)
#
#         temp1 = fmt.r_basis([self.air_gap_h[i]], [self.core_w[i]], [h_multiple])
#         temp2 = fmt.sigma([self.air_gap_h[i]], [self.core_w[i] / 2], [temp1])
#         temp3 = fmt.r_round_inf([self.air_gap_h[i]], [temp2], [self.core_w[i] / 2])
#         self.reluctance[i, 5] = self.reluctance[i, 5] + (2 * temp3)
#
#         r_basis_1 = fmt.r_basis([self.air_gap_h[i] / 2], [self.core_w[i]], [h_multiple])
#         r_basis_2 = fmt.r_basis([self.air_gap_h[i] / 2], [self.core_w[i]], [h_multiple])
#         temp2 = fmt.sigma([self.air_gap_h[i]], [self.core_w[i] / 2], [r_basis_1 + r_basis_2])
#         temp3 = fmt.femmt_functions.r_round_round([self.air_gap_h[i]], [temp2], [self.core_w[i] / 2])
#         self.reluctance[i, 5] = self.reluctance[i, 5] + ((self.n_air_gaps[i] - 2) * temp3)
#     else:
#         h_multiple = (self.window_h[i] - (self.n_air_gaps[i] * self.air_gap_h[i])) / (self.n_air_gaps[i] + 1)
#
#         r_basis_1 = fmt.r_basis([self.air_gap_h[i] / 2], [self.core_w[i]], [h_multiple])
#         r_basis_2 = fmt.r_basis([self.air_gap_h[i] / 2], [self.core_w[i]], [h_multiple / 2])
#         temp2 = fmt.sigma([self.air_gap_h[i]], [self.core_w[i] / 2], [r_basis_1 + r_basis_2])
#         temp3 = fmt.femmt_functions.r_round_round([self.air_gap_h[i]], [temp2], [self.core_w[i] / 2])
#         self.reluctance[i, 5] = self.reluctance[i, 5] + (2 * temp3)
#
#         r_basis_1 = fmt.r_basis([self.air_gap_h[i] / 2], [self.core_w[i]], [h_multiple / 2])
#         r_basis_2 = fmt.r_basis([self.air_gap_h[i] / 2], [self.core_w[i]], [h_multiple / 2])
#         temp2 = fmt.sigma([self.air_gap_h[i]], [self.core_w[i] / 2], [r_basis_1 + r_basis_2])
#         temp3 = fmt.femmt_functions.r_round_round([self.air_gap_h[i]], [temp2], [self.core_w[i] / 2])
#         self.reluctance[i, 5] = self.reluctance[i, 5] + ((self.n_air_gaps[i] - 2) * temp3)


# self.reluctance[self.single_air_gap_len:self.data_matrix_len, 5] = np.where(self.mult_air_gap_type[self.single_air_gap_len:self.data_matrix_len] == 1,
#     (2 * fmt.r_round_inf(self.air_gap_h[self.single_air_gap_len:self.data_matrix_len],
#                          fmt.sigma(self.air_gap_h[self.single_air_gap_len:self.data_matrix_len],
#                                    self.core_w[self.single_air_gap_len:self.data_matrix_len] / 2,
#                                    fmt.r_basis(self.air_gap_h[self.single_air_gap_len:self.data_matrix_len],
#                                                self.core_w[self.single_air_gap_len:self.data_matrix_len], h_multiple)),
#                          self.core_w[self.single_air_gap_len:self.data_matrix_len] / 2)) + (
#                 (self.n_air_gaps[self.single_air_gap_len:self.data_matrix_len] - 2) * fmt.r_round_round(
#             self.air_gap_h[self.single_air_gap_len:self.data_matrix_len],
#             fmt.sigma(self.air_gap_h[self.single_air_gap_len:self.data_matrix_len],
#                       self.core_w[self.single_air_gap_len:self.data_matrix_len] / 2,
#                       2 * fmt.r_basis([self.air_gap_h[self.single_air_gap_len:self.data_matrix_len] / 2],
#                                       [self.core_w[self.single_air_gap_len:self.data_matrix_len]], [h_multiple])),
#             self.core_w[self.single_air_gap_len:self.data_matrix_len] / 2)),
#     (2 * fmt.r_round_round(self.air_gap_h[self.single_air_gap_len:self.data_matrix_len],
#                                            fmt.sigma(self.air_gap_h[self.single_air_gap_len:self.data_matrix_len],
#                                                      self.core_w[self.single_air_gap_len:self.data_matrix_len] / 2,
#                                                      fmt.r_basis(
#                                                          self.air_gap_h[self.single_air_gap_len:self.data_matrix_len] / 2,
#                                                          self.core_w[self.single_air_gap_len:self.data_matrix_len],
#                                                          h_multiple) + fmt.r_basis(
#                                                          self.air_gap_h[self.single_air_gap_len:self.data_matrix_len] / 2,
#                                                          self.core_w[self.single_air_gap_len:self.data_matrix_len],
#                                                          h_multiple / 2)),
#                                            self.core_w[self.single_air_gap_len:self.data_matrix_len] / 2)) + (
#                 (self.n_air_gaps[self.single_air_gap_len:self.data_matrix_len] - 2) * fmt.femmt_functions.r_round_round(
#             self.air_gap_h[self.single_air_gap_len:self.data_matrix_len],
#             fmt.sigma(self.air_gap_h[self.single_air_gap_len:self.data_matrix_len],
#                       self.core_w[self.single_air_gap_len:self.data_matrix_len] / 2,
#                       2 * fmt.r_basis(self.air_gap_h[self.single_air_gap_len:self.data_matrix_len] / 2,
#                                       self.core_w[self.single_air_gap_len:self.data_matrix_len], h_multiple / 2)),
#             self.core_w[self.single_air_gap_len:self.data_matrix_len] / 2)))
