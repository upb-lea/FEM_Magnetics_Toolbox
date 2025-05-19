"""Create and calculate reluctance models."""
# python imports

# 3rd library imports
import numpy as np
from matplotlib import pyplot as plt

# femmt imports
import femmt.functions_reluctance as fr
from femmt.constants import mu_0


def plot_limitation():
    """Plot limitation."""
    length = 15
    width = 100 * length
    height = 101 - length

    r_mx = 1 / (mu_0 * (width / 2 / length + 2 / np.pi * (
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
    height2_c = 100 - height1_c
    h_l = height2_c / length
    # print(h_l)
    r_m1 = 1 / (mu_0 * (width_c / 2 / length_c + 2 / np.pi * (1 + np.log(np.pi * height1_c / 4 / length_c))))

    r_m2 = 1 / (mu_0 * (width_c / 2 / length_c + 2 / np.pi * (1 + np.log(np.pi * height2_c / 4 / length_c))))

    r_m = r_m1 + r_m2
    ratio = r_mx / r_m
    # print(ratio)
    fig, ax = plt.subplots()  # Create a figure containing a single axes.
    plt.title("R_basic vs h/l")
    plt.xlabel("h/l")
    plt.ylabel("R_basic")

    ax.plot(h_l, r_m)
    # ax.hlines(y=r_m, xmin=0, xmax=50, linewidth=2, color='g')
    ax.hlines(y=r_mx, xmin=-1, xmax=51, linewidth=2, color='r')
    ax.invert_xaxis()
    ax.grid()
    plt.show()


def plot_r_basic():
    """
    Plot the 2D Reluctance of the basic geometry described in Muelethaler thesis.

    (using Schwarz-Christoffel transformation) at page no. 35.
    It plots the Reluctance formula with respect to its variables (h/l and w/l).
    It is an independent function and has been used to analyze the expression and its limitation.
    """
    # # Uncomment this section of code if Reluctance change with respect to h/l is desired
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
    # fig, ax = fmt.plt.subplots(figsize=(3.54, 3.54), dpi=150)
    # # fmt.plt.title("R_basic vs h/l")
    # fmt.plt.xlabel("$\dfrac{h}{l}$", fontsize=24)
    # fmt.plt.ylabel("$R_{\mathrm{basic}}^{\prime}$ / AT/Wb", fontsize=24)
    # ax.plot(h_l, r_m, linewidth=4, label=f'w/l ={width}')
    # ax.invert_xaxis()
    # ax.legend(fontsize=24)
    # ax.grid()
    # fmt.plt.show()

    # # Uncomment this section of code if Reluctance change with respect to w/l is desired
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

    # Uncomment this section of code if Reluctance change with respect to h/l and w/l is desired
    width = np.linspace(100, 0, 5)
    length = 1
    height = np.linspace(100, 0, 1000)
    h_l = height / length
    w_l = width / length
    fig, ax = plt.subplots(figsize=(3.54, 3.54), dpi=150)  # Create a figure containing a single axes.
    # fmt.plt.title("$R_{basic}$ vs $\dfrac{h}{l}$", fontsize=20)
    plt.xlabel(r"$\dfrac{h}{l}$", fontsize=24)
    plt.ylabel(r"$R_{\mathrm{basic}}^{\prime}$ / AT/Wb", fontsize=24)

    for i, wid in enumerate(width):
        r_m = 1 / (mu_0 * (wid / 2 / length + 2 / np.pi * (1 + np.log(np.pi * height / 4 / length))))

        combined = np.vstack((h_l, r_m)).T
        # print(combined)

        ax.plot(h_l, r_m, linewidth=4, label=f'w/l ={w_l[i]}')

    ax.invert_xaxis()
    ax.set_yscale('log')
    ax.legend(fontsize=24)
    ax.grid()
    plt.show()


def distributed_type_1(air_gap_height_single_air_gap, core_inner_diameter, n_air_gaps, h_multiple):
    """Calculate distributed air-gap reluctance of Type 1 (Where corner air-gaps are present).

    :param air_gap_height_single_air_gap: Air-gap height [in meter]
    :type air_gap_height_single_air_gap: list
    :param core_inner_diameter: Diameter of center leg of the core [in meter]
    :type core_inner_diameter: list
    :param n_air_gaps: Number of air-gaps in the center leg of the core
    :type n_air_gaps: list
    :param h_multiple: Half of core height between two consecutive air-gaps in an equally distributed air-gaps [in meter]
    :type h_multiple: ndarray
    :return: Distributed air-gap reluctance of Type 1 (Where corner air-gaps are present)
    :rtype: list
    """
    # ToDo: Raise Error for less than two air gaps

    # first part calculates the two outer air gaps (very top and very bottom)
    reluctance_round_inf = fr.r_air_gap_round_inf(air_gap_height_single_air_gap, core_inner_diameter, h_multiple)
    reluctance_total = (2 * reluctance_round_inf)

    # second part calculates the inner air gaps between top and bottom air gaps (if available)
    reluctance_round_round = fr.r_air_gap_round_round(air_gap_height_single_air_gap, core_inner_diameter, h_multiple,
                                                      h_multiple)
    reluctance_total = reluctance_total + ((n_air_gaps - 2) * reluctance_round_round)

    return reluctance_total


def distributed_type_2(air_gap_height_single_air_gap, core_inner_diameter, n_air_gaps, h_multiple):
    """Calculate distributed air-gap reluctance of Type 2 (Where corner air-gaps are absent).

    :param air_gap_height_single_air_gap: Air-gap height [in meter]
    :type air_gap_height_single_air_gap: list
    :param core_inner_diameter: Diameter of center leg of the core [in meter]
    :type core_inner_diameter: list
    :param n_air_gaps: Number of air-gaps in the center leg of the core
    :type n_air_gaps: list
    :param h_multiple: Core height between two consecutive air-gaps in an equally distributed air-gaps [in meter]
    :type h_multiple: ndarray
    :return: Distributed air-gap reluctance of Type 2 (Where corner air-gaps are absent)
    :rtype: list
    """
    # ToDo: Raise Error for less than two air gaps

    # First part calculates two outer air gaps (very top and very bottom)
    reluctance_round_round = fr.r_air_gap_round_round(air_gap_height_single_air_gap, core_inner_diameter, h_multiple,
                                                      h_multiple / 2)
    reluctance = (2 * reluctance_round_round)

    # second part calculates air gaps between the outer air gaps (if available)
    reluctance_round_round = fr.r_air_gap_round_round(air_gap_height_single_air_gap, core_inner_diameter, h_multiple / 2,
                                                      h_multiple / 2)
    reluctance = reluctance + ((n_air_gaps - 2) * reluctance_round_round)

    return reluctance


def create_data_matrix(core_inner_diameter: list, window_h: list, window_w: list, no_of_turns: list,
                       n_air_gaps: list,
                       air_gap_h: list, air_gap_position: list, mu_rel: list, mult_air_gap_type: list):
    """Create matrix consisting of input design parameters with all their combinations.

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
        Type 1: Equally distributed air-gaps including corner air-gaps (e.g.: air-gaps-position = [0, 50, 100] for 3 air-gaps)
        Type 2: Equally distributed air-gaps excluding corner air-gaps (e.g.: air-gaps-position = [25, 50, 75] for 3 air-gaps)
    :type mult_air_gap_type: list
    """
    # Structure: data_matrix = [core_inner_diameter, window_h, window_w, mu_rel, no_of_turns, n_air_gaps, air_gap_h,
    #                      air_gap_position, mult_air_gap_type, inductance]
    clone_n_air_gaps = n_air_gaps
    row_num = 0

    if 1 in clone_n_air_gaps:
        data_matrix = np.zeros((len(core_inner_diameter) * len(no_of_turns) * len(air_gap_h) * len(mu_rel) * (
            len(air_gap_position) + (len(n_air_gaps) - 1) * len(mult_air_gap_type)), 10))
    else:
        data_matrix = np.zeros(
            (len(core_inner_diameter) * len(no_of_turns) * len(air_gap_h) * len(n_air_gaps) * len(
                mu_rel) * len(mult_air_gap_type), 10))

    if 1 in clone_n_air_gaps:
        clone_n_air_gaps.remove(1)
        for index_1 in range(len(core_inner_diameter)):
            for index_2 in range(len(mu_rel)):
                for index_3 in range(len(no_of_turns)):
                    for index_4 in range(len(air_gap_h)):
                        for index_5 in range(len(air_gap_position)):
                            data_matrix[row_num, 0] = core_inner_diameter[index_1]
                            data_matrix[row_num, 1] = window_h[index_1]
                            data_matrix[row_num, 2] = window_w[index_1]
                            data_matrix[row_num, 3] = mu_rel[index_2]
                            data_matrix[row_num, 4] = no_of_turns[index_3]
                            data_matrix[row_num, 5] = 1
                            data_matrix[row_num, 6] = air_gap_h[index_4]
                            data_matrix[row_num, 7] = air_gap_position[index_5]
                            data_matrix[row_num, 8] = np.nan
                            row_num = row_num + 1

    single_air_gap_len = row_num

    if len(clone_n_air_gaps) != 0:
        for index_1 in range(len(core_inner_diameter)):
            for index_2 in range(len(mu_rel)):
                for index_3 in range(len(no_of_turns)):
                    for index_4 in range(len(clone_n_air_gaps)):
                        for index_5 in range(len(air_gap_h)):
                            for index_6 in range(len(mult_air_gap_type)):
                                data_matrix[row_num, 0] = core_inner_diameter[index_1]
                                data_matrix[row_num, 1] = window_h[index_1]
                                data_matrix[row_num, 2] = window_w[index_1]
                                data_matrix[row_num, 3] = mu_rel[index_2]
                                data_matrix[row_num, 4] = no_of_turns[index_3]
                                data_matrix[row_num, 5] = n_air_gaps[index_4]
                                data_matrix[row_num, 6] = air_gap_h[index_5]
                                data_matrix[row_num, 7] = np.nan
                                data_matrix[row_num, 8] = mult_air_gap_type[index_6]
                                row_num = row_num + 1

    data_matrix_len = row_num

    return data_matrix, single_air_gap_len, data_matrix_len


class MagneticCircuit:
    """Class object for calculating the reluctance and inductance of 2D-axis symmetric inductor."""

    def __init__(self, core_inner_diameter: list, window_h: list, window_w: list, no_of_turns: list, n_air_gaps: list,
                 air_gap_h: list, air_gap_position: list, mu_r_abs: list, mult_air_gap_type: list = None,
                 air_gap_method: str = 'Percent', component_type: str = 'inductor', sim_type: str = 'single'):
        """
        Initialize the MagneticCircuit class.

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
        :param mu_r_abs: Relative permeability of the core [in F/m]
        :type mu_r_abs: list
        :param mult_air_gap_type: Two types of equally distributed air-gaps (used only for air-gaps more than 1)
            Type 1: Equally distributed air-gaps including corner air-gaps (e.g.: air-gaps-position = [0, 50, 100] for 3 air-gaps)
            Type 2: Equally distributed air-gaps excluding corner air-gaps (e.g.: air-gaps-position = [25, 50, 75] for 3 air-gaps)
        :type mult_air_gap_type: list
                :param air_gap_h: Air-gap height [in meter]
        :param air_gap_method: Input method of air gap position ( either in 'Percent', 'Center' or 'Manually')
        :type air_gap_method: str
        :param component_type: Position of the air-gap in the percentage with respect to window_h
        :type component_type: str
        :param sim_type: Relative permeability of the core [in F/m]
        :type sim_type: str
        """
        # Storing input arguments into object variables
        self.core_inner_diameter = core_inner_diameter
        self.window_h = window_h
        self.window_w = window_w
        self.no_of_turns = no_of_turns
        self.n_air_gaps = n_air_gaps
        self.air_gap_h = air_gap_h
        self.air_gap_position = air_gap_position
        self.mu_r_abs = mu_r_abs
        self.mult_air_gap_type = mult_air_gap_type
        self.air_gap_method = air_gap_method
        self.sim_type = sim_type
        self.component_type = component_type

        # Checking input variables
        self.input_pre_check()

        # Definition of data matrix (which stores all inputs and results) and its length variables
        self.data_matrix = np.zeros((1, 10))
        self.single_air_gap_len = None
        self.data_matrix_len = None

        # Definition of geometry variables which are calculated using the basic 3 available dimensions
        self.core_h = None
        self.r_outer = None
        self.r_inner = None
        self.core_h_middle = None  # height of upper/lower part of the window in the core
        self.outer_leg_w = None

        # Definition of other variables required for various calculations
        self.abs_position_air_gap = None
        self.cal_inductance = None
        self.section = None
        self.length = None
        self.area = None
        self.reluctance = None

        self.max_percent_position = None
        self.min_percent_position = None

        # Dictionary to store the column names of the data_matrix (2D array)
        self.param_pos_dict = {"core_inner_diameter": 0, "window_h": 1, "window_w": 2, "mu_r_abs": 3, "no_of_turns": 4,
                               "n_air_gaps": 5, "air_gap_h": 6, "air_gap_position": 7, "mult_air_gap_type": 8,
                               "inductance": 9}

        if self.component_type == 'inductor':
            if sim_type == 'single':
                # Conversion of list to numpy array
                self.core_inner_diameter = np.array(self.core_inner_diameter)
                self.window_h = np.array(self.window_h)
                self.window_w = np.array(self.window_w)
                self.no_of_turns = np.array(self.no_of_turns)
                self.n_air_gaps = np.array(self.n_air_gaps)
                self.air_gap_h = np.array(self.air_gap_h)
                self.air_gap_position = np.array(self.air_gap_position)
                self.mu_r_abs = np.array(self.mu_r_abs)

                # Call to functions for single inductance calculation
                self.core_reluctance()
                self.air_gap_reluctance_single()
                self.calculate_inductance()
            else:
                # Creates the data matrix with all the input parameter combinations for sim_type = 'sweep'
                self.data_matrix, self.single_air_gap_len, self.data_matrix_len = \
                    create_data_matrix(self.core_inner_diameter, self.window_h, self.window_w, self.no_of_turns,
                                       self.n_air_gaps, self.air_gap_h, self.air_gap_position, self.mu_r_abs,
                                       self.mult_air_gap_type)

                # Mapping variables to data_matrix columns
                self.core_inner_diameter = self.data_matrix[:, 0]
                self.window_h = self.data_matrix[:, 1]
                self.window_w = self.data_matrix[:, 2]
                self.mu_r_abs = self.data_matrix[:, 3]
                self.no_of_turns = self.data_matrix[:, 4]
                self.n_air_gaps = self.data_matrix[:, 5]
                self.air_gap_h = self.data_matrix[:, 6]
                self.air_gap_position = self.data_matrix[:, 7]
                self.mult_air_gap_type = self.data_matrix[:, 8]

                # Call to functions for sweep inductance calculation
                self.core_reluctance()
                self.air_gap_reluctance_sweep()
                self.calculate_inductance()

                # Adding other important variables to data_matrix
                self.data_matrix = self.add_column_to_data_matrix(self.data_matrix, self.core_h_middle, 'core_h_middle')  # 10
                self.data_matrix = self.add_column_to_data_matrix(self.data_matrix, self.r_inner, 'r_inner')  # 11
                self.data_matrix = self.add_column_to_data_matrix(self.data_matrix, self.r_outer, 'r_outer')  # 12
                self.data_matrix = self.add_column_to_data_matrix(self.data_matrix, self.area[:, 0], 'center_leg_area')  # 13
                self.data_matrix = self.add_column_to_data_matrix(self.data_matrix, self.area[:, 4], 'outer_leg_area')  # 14
                self.data_matrix = self.add_column_to_data_matrix(self.data_matrix, self.core_h, 'core_h')  # 15

    def input_pre_check(self):
        """Check the correctness of the inputs provided to class MagneticCircuit."""
        if not (len(self.core_inner_diameter) and len(self.window_h) and len(self.window_w) and len(self.no_of_turns) \
                and len(self.n_air_gaps) and len(self.mu_r_abs)):
            raise Exception("one of the passed arguments are empty list")
        if not all(isinstance(item, int) for item in self.no_of_turns):
            raise Exception("no_of_turns list elements should be integer")
        if not all(isinstance(item, int) for item in self.n_air_gaps):
            raise Exception("n_air_gaps list elements should be integer")
        if not (self.air_gap_method == 'Center' or self.air_gap_method == 'Percent' or self.air_gap_method == 'Manually'):
            raise Exception("string value wrong for air_gap_method argument")
        if not (self.sim_type == 'single' or self.sim_type == 'sweep'):
            raise Exception("string value wrong for sim_type argument")
        if not (self.component_type == 'inductor' or self.component_type == 'integrated_transformer'):
            raise Exception("string value wrong for component_type argument")
        # if any(item > 0.0005 for item in self.air_gap_h):
        #    raise Exception("Model accuracy is not good for air_gap_h more than 0.0005")
        if self.sim_type == 'single':
            if not (len(self.core_inner_diameter) == 1 and len(self.window_h) == 1 and len(self.window_w) == 1 and len(
                    self.no_of_turns) == 1 and len(self.n_air_gaps) == 1 and len(self.mu_r_abs) == 1):
                raise Exception("single sim_type requires single list elements")
            if not (self.n_air_gaps[0] == len(self.air_gap_h) and self.n_air_gaps[0] == len(self.air_gap_position)):
                raise Exception("No. of elements of air_gap_h and air_gap_position should match n_air_gaps")

            # Sort air_gap_position and air_gap_h based on air_gap_position
            if self.n_air_gaps[0] != 0:
                self.air_gap_position, self.air_gap_h = [list(tpl) for tpl in
                                                         zip(*sorted(zip(self.air_gap_position, self.air_gap_h)))]
                self.air_gap_h = np.array(self.air_gap_h)
                self.air_gap_position = np.array(self.air_gap_position)
                print(self.air_gap_position)
                print(self.air_gap_h)
        if self.sim_type == 'sweep':
            if not self.air_gap_method == 'Percent':
                raise Exception("For sim_type = 'sweep', 'air_gap_method' should be only provided in percent")

    def core_reluctance(self):
        """
        Calculate the core reluctance along with length and area of each section of the core geometry.

        Core reluctance are referred from Appendix B of book "E. C. Snelling. Soft Ferrites,
        Properties and Applications. 2nd edition. Butterworths, 1988". This book is referred in Muelethaler thesis
        at page no. 26
        """
        # Geometry definitions
        self.core_h = self.window_h + self.core_inner_diameter / 2
        self.core_h_middle = self.core_inner_diameter / 2 / 2
        self.r_inner = self.core_inner_diameter / 2 + self.window_w
        self.r_outer = np.sqrt((self.core_inner_diameter / 2) ** 2 + self.r_inner ** 2)
        self.outer_leg_w = self.r_outer - self.r_inner

        # Section [0]: Center leg; [1]: Inner corner; [2]: Winding window; [3]: Outer corner; [4]: Outer leg
        self.section = [0, 1, 2, 3, 4]
        self.length = np.zeros((len(self.data_matrix), len(self.section)))
        self.area = np.zeros((len(self.data_matrix), len(self.section)))

        # One more column is added for reluctance matrix to accommodate air-gap reluctance later on
        self.reluctance = np.zeros((len(self.data_matrix), len(self.section) + 1))

        # Section 0: Center leg
        if self.sim_type == 'sweep':
            self.length[:, 0] = self.window_h - (self.n_air_gaps * self.air_gap_h)
        else:
            self.length[:, 0] = self.window_h - sum(self.air_gap_h)
        self.area[:, 0] = np.pi * ((self.core_inner_diameter / 2) ** 2)

        # Section 1: Inner corner
        s_1 = (self.core_inner_diameter / 2) - (self.core_inner_diameter / (2 * np.sqrt(2)))
        self.length[:, 1] = (np.pi / 4) * (s_1 + (self.core_h_middle / 2))
        self.area[:, 1] = (np.pi / 2) * ((2 * (self.core_inner_diameter / 2) * self.core_h_middle) + ((self.core_inner_diameter / 2) ** 2))

        # Section 2: Winding window section
        self.length[:, 2] = self.window_w
        self.area[:, 2] = np.nan

        # Section 3: Outer corner
        s_2 = np.sqrt(((self.r_inner ** 2) + (self.r_outer ** 2)) / 2) - self.r_inner
        self.length[:, 3] = (np.pi / 4) * (s_2 + (self.core_h_middle / 2))
        self.area[:, 3] = (np.pi / 2) * ((2 * self.r_inner * self.core_h_middle) + ((self.r_outer ** 2) - (self.r_inner ** 2)))

        # Section 4: Outer leg
        self.length[:, 4] = self.window_h
        self.area[:, 4] = np.pi * (self.r_outer ** 2 - self.r_inner ** 2)

        # Reluctance calculation
        self.reluctance[:, 0] = self.length[:, 0] / (mu_0 * self.mu_r_abs * self.area[:, 0])
        self.reluctance[:, 1] = self.length[:, 1] / (mu_0 * self.mu_r_abs * self.area[:, 1])
        self.reluctance[:, 2] = ((mu_0 * self.mu_r_abs * 2 * np.pi * self.core_h_middle) ** -1) * np.log((2 * self.r_inner) / self.core_inner_diameter)
        self.reluctance[:, 3] = self.length[:, 3] / (mu_0 * self.mu_r_abs * self.area[:, 3])
        self.reluctance[:, 4] = self.length[:, 4] / (mu_0 * self.mu_r_abs * self.area[:, 4])

        # Doubling the reluctance of section 1, 2, 3, to accommodate horizontal symmetry
        self.reluctance[:, 1:4] = 2 * self.reluctance[:, 1:4]

    def air_gap_reluctance_sweep(self):
        """
        Calculate the air-gap reluctance for a sweep simulation with single/distributed air-gaps.

        Method according to the following paper:
        ["A Novel Approach for 3D Air Gap Reluctance Calculations" - J. Mühlethaler, J.W. Kolar, A. Ecklebe]

        It is called when input variable sim_type == 'sweep'.
        Its calculation for multiple air-gap is based on
        series connection. That is, multiple air-gaps are divided equally and common height (h) is calculated.
        Then their Reluctances are calculated and added such that they are connected in series.
        """
        # Single air-gap reluctance calculations
        self.max_percent_position = ((self.window_h[0:self.single_air_gap_len] - (
            self.air_gap_h[0:self.single_air_gap_len] / 2)) / self.window_h[0:self.single_air_gap_len]) * 100
        self.min_percent_position = ((self.air_gap_h[0:self.single_air_gap_len] / 2) / self.window_h[0:self.single_air_gap_len]) * 100
        # Convert percent position to absolute value position
        self.abs_position_air_gap = (self.air_gap_position[0:self.single_air_gap_len] * self.window_h[0:self.single_air_gap_len]) / 100
        h = np.zeros((len(self.abs_position_air_gap), 2))

        h[:, 0] = np.where((self.air_gap_position[0:self.single_air_gap_len] <= self.min_percent_position) | (
            self.air_gap_position[0:self.single_air_gap_len] >= self.max_percent_position),
            self.window_h[0:self.single_air_gap_len] - self.air_gap_h[0:self.single_air_gap_len],
            self.abs_position_air_gap - (self.air_gap_h[0:self.single_air_gap_len] / 2))
        h[:, 1] = np.where((self.air_gap_position[0:self.single_air_gap_len] <= self.min_percent_position) | (
            self.air_gap_position[0:self.single_air_gap_len] >= self.max_percent_position), 0,
            self.window_h[0:self.single_air_gap_len] - self.abs_position_air_gap - (self.air_gap_h[0:self.single_air_gap_len] / 2))

        self.reluctance[0:self.single_air_gap_len, 5] = np.where(h[:, 1] == 0, fr.r_air_gap_round_inf(
            self.air_gap_h[0:self.single_air_gap_len], self.core_inner_diameter[0:self.single_air_gap_len],
            h[:, 0]), fr.r_air_gap_round_round(self.air_gap_h[0:self.single_air_gap_len], self.core_inner_diameter[0:self.single_air_gap_len],
                                               h[:, 0], h[:, 1]))

        # Distributed air-gaps reluctance calculations
        h_multiple = np.where(self.mult_air_gap_type[self.single_air_gap_len:self.data_matrix_len] == 1,
                              (self.window_h[self.single_air_gap_len:self.data_matrix_len] - (
                               self.n_air_gaps[self.single_air_gap_len:self.data_matrix_len] * \
                               self.air_gap_h[self.single_air_gap_len:self.data_matrix_len])) / (
                                   (self.n_air_gaps[self.single_air_gap_len:self.data_matrix_len] - 1) * 2),
                              (self.window_h[self.single_air_gap_len:self.data_matrix_len] - (
                                  self.n_air_gaps[self.single_air_gap_len:self.data_matrix_len] * \
                                  self.air_gap_h[self.single_air_gap_len:self.data_matrix_len])) / (
                                  self.n_air_gaps[self.single_air_gap_len:self.data_matrix_len] + 1))

        self.reluctance[self.single_air_gap_len:self.data_matrix_len, 5] = np.where(
            self.mult_air_gap_type[self.single_air_gap_len:self.data_matrix_len] == 1,
            distributed_type_1(self.air_gap_h[self.single_air_gap_len:self.data_matrix_len],
                               self.core_inner_diameter[self.single_air_gap_len:self.data_matrix_len],
                               self.n_air_gaps[self.single_air_gap_len:self.data_matrix_len], h_multiple),
            distributed_type_2(self.air_gap_h[self.single_air_gap_len:self.data_matrix_len],
                               self.core_inner_diameter[self.single_air_gap_len:self.data_matrix_len],
                               self.n_air_gaps[self.single_air_gap_len:self.data_matrix_len], h_multiple))

    def air_gap_reluctance_sweep_new(self):
        """
        Calculate the air-gap reluctance for a sweep simulation with single/distributed air-gaps.

        Method according to the following paper:
        ["A Novel Approach for 3D Air Gap Reluctance Calculations" - J. Mühlethaler, J.W. Kolar, A. Ecklebe]

        It is called when input variable sim_type == 'sweep'.
        Its calculation for multiple air-gap is based on superposition.
        That is, multiple air-gap's reluctance is calculated by taking one at a time and then adding them together
        (like in superposition theorem)
        """
        # Single air-gap reluctance calculations
        # Convert percent position to absolute position
        air_gap_position_absolute = (self.air_gap_position[0:self.single_air_gap_len] / 100) * self.window_h[0:self.single_air_gap_len]

        core_height_lower = air_gap_position_absolute - (self.air_gap_h[0:self.single_air_gap_len] / 2)
        core_height_upper = self.window_h[0:self.single_air_gap_len] - air_gap_position_absolute - (self.air_gap_h[0:self.single_air_gap_len] / 2)

        core_height_lower = np.where(core_height_lower <= 0, 0, core_height_lower)
        core_height_upper = np.where(core_height_upper <= 0, 0, core_height_upper)

        self.reluctance[0:self.single_air_gap_len, 5] = np.where(core_height_upper * core_height_lower == 0,
                                                                 fr.r_air_gap_round_inf(
                                                                     self.air_gap_h[0:self.single_air_gap_len],
                                                                     self.core_inner_diameter[0:self.single_air_gap_len],
                                                                     core_height_lower + core_height_upper),
                                                                 fr.r_air_gap_round_round(
                                                                     self.air_gap_h[0:self.single_air_gap_len],
                                                                     self.core_inner_diameter[0:self.single_air_gap_len],
                                                                     core_height_upper,
                                                                     core_height_lower))

        # Distributed air-gaps reluctance calculations
        for i in range(self.single_air_gap_len, self.data_matrix_len):
            if self.mult_air_gap_type[i] == 1:
                air_gap_position_percent = np.arange(0, 101, 100 / (self.n_air_gaps[i] - 1))
            else:
                air_gap_position_percent = np.arange(100 / (self.n_air_gaps[i] + 1), 100,
                                                     100 / (self.n_air_gaps[i] + 1))

            air_gap_position_abs = (air_gap_position_percent / 100) * self.window_h[i]

            h_upper = self.window_h[i] - air_gap_position_abs - (self.air_gap_h[i] / 2)
            h_lower = air_gap_position_abs - (self.air_gap_h[i] / 2)

            h_upper = np.where(h_upper <= 0, 0, h_upper)
            h_lower = np.where(h_lower <= 0, 0, h_lower)

            self.reluctance[i, 5] = np.sum(np.where(h_lower * h_upper == 0,
                                                    fr.r_air_gap_round_inf(self.air_gap_h[i],
                                                                           self.core_inner_diameter[i],
                                                                           h_lower + h_upper),
                                                    fr.r_air_gap_round_round(self.air_gap_h[i],
                                                                             self.core_inner_diameter[i],
                                                                             h_upper, h_lower)))
        # print(f"reluctances:{self.reluctance[:, 5]}")

    def get_parameters_position_dict(self):
        """
        Return dictionary 'param_pos_dict'.

        Used to refer the column number of data_matrix by using the column names.
        """
        return self.param_pos_dict

    def air_gap_reluctance_single(self):
        """
        Calculate the air-gap reluctance for a single simulation with single/multiple air-gaps.

        Method according to the following paper:
        ["A Novel Approach for 3D Air Gap Reluctance Calculations" - J. Mühlethaler, J.W. Kolar, A. Ecklebe]

        It is called when input variable sim_type == 'single'.
        Its calculation for multiple air-gap is based on
        series connection. That is, multiple air-gaps are divided equally and common height (h) is calculated.
        Then their reluctances are calculated and added such that they are connected in series.
        """
        flag_0 = 0
        flag_1 = 0
        flag_2 = 0
        if self.n_air_gaps[0] != 0:
            if self.air_gap_method == 'Center':
                self.section.append(6)  # round-round type air gap

                self.reluctance[:, 5] = fr.r_air_gap_round_round([self.air_gap_h[0]], [self.core_inner_diameter],
                                                                 [(self.window_h - self.air_gap_h[0]) / 2],
                                                                 [(self.window_h - self.air_gap_h[0]) / 2])

            elif self.air_gap_method == 'Percent' or self.air_gap_method == 'Manually':
                self.max_percent_position = ((self.window_h - self.air_gap_h[
                    self.n_air_gaps - 1] / 2) / self.window_h) * 100
                self.min_percent_position = ((self.air_gap_h[0] / 2) / self.window_h) * 100
                if self.air_gap_method == 'Percent':
                    self.position = np.array(
                        self.air_gap_position) / 100 * self.window_h  # Convert percent position to absolute value position
                print(f"Max percent: {self.max_percent_position}")
                print(f"Min percent: {self.min_percent_position}")

                if self.air_gap_position[0] <= self.min_percent_position:
                    flag_0 = 1
                    self.section.append(8)
                    if self.n_air_gaps == 1:
                        h = self.window_h - self.air_gap_h[0]
                    else:
                        h = ((self.position[1] - self.air_gap_h[1] / 2) - self.air_gap_h[0]) / 2

                    self.reluctance[:, 5] = self.reluctance[:, 5] + fr.r_air_gap_round_inf([self.air_gap_h[0]],
                                                                                           [self.core_inner_diameter],
                                                                                           [h])
                    print('air gap is at lower corner')

                if self.air_gap_position[self.n_air_gaps - 1] >= self.max_percent_position:
                    flag_1 = 1
                    self.section.append(8)
                    if self.n_air_gaps == 1:
                        h = self.window_h - self.air_gap_h[self.n_air_gaps - 1]
                    else:
                        h = (self.position[self.n_air_gaps - 1] - self.position[self.n_air_gaps - 2] - self.air_gap_h[
                            self.n_air_gaps - 1] / 2 - self.air_gap_h[self.n_air_gaps - 2] / 2) / 2

                    self.reluctance[:, 5] = self.reluctance[:, 5] + fr.r_air_gap_round_inf(
                        [self.air_gap_h[self.n_air_gaps - 1]], [self.core_inner_diameter], [h])
                    print('air gap is at upper corner')

                for i in range(self.n_air_gaps[0]):
                    if self.min_percent_position < self.air_gap_position[i] < self.max_percent_position:
                        self.section.append(7)
                        if flag_2 == 0:
                            if flag_0 == 0 and flag_1 == 0:  # No corner air-gaps
                                self.position = np.append(self.position,
                                                          self.window_h + (self.window_h - self.position[
                                                              self.n_air_gaps - 1]))
                                self.position = np.insert(self.position, 0, -self.position[0])
                                self.air_gap_h = np.append(self.air_gap_h, self.air_gap_h[self.n_air_gaps - 1])
                                self.air_gap_h = np.insert(self.air_gap_h, 0, self.air_gap_h[0])
                            elif flag_0 == 1 and flag_1 == 0:  # Only lower air-gap is present
                                self.position = np.append(self.position,
                                                          self.window_h + (self.window_h - self.position[
                                                              self.n_air_gaps - 1]))
                                self.air_gap_h = np.append(self.air_gap_h, self.air_gap_h[self.n_air_gaps - 1])
                            elif flag_0 == 0 and flag_1 == 1:  # Only Upper air-gap is present
                                self.position = np.insert(self.position, 0, -self.position[0])
                                self.air_gap_h = np.insert(self.air_gap_h, 0, self.air_gap_h[0])
                            flag_2 = 1

                        if flag_0 == 0 and flag_1 == 0:
                            h1 = (self.position[i + 1] - self.position[i] - self.air_gap_h[i + 1] / 2 - self.air_gap_h[i] / 2) / 2
                            h2 = (self.position[i + 2] - self.position[i + 1] - self.air_gap_h[i + 2] / 2 - self.air_gap_h[i + 1] / 2) / 2
                            print('No corner air gap detected')
                        elif flag_0 == 1 and flag_1 == 0:
                            h1 = (self.position[i] - self.position[i - 1] - self.air_gap_h[i] / 2 - self.air_gap_h[
                                i - 1] / 2) / 2
                            h2 = (self.position[i + 1] - self.position[i] - self.air_gap_h[i + 1] / 2 - self.air_gap_h[
                                i] / 2) / 2
                            print('Lower air gap detected')
                        elif flag_0 == 0 and flag_1 == 1:
                            h1 = (self.position[i + 1] - self.position[i] - self.air_gap_h[i + 1] / 2 - self.air_gap_h[i] / 2) / 2
                            h2 = (self.position[i + 2] - self.position[i + 1] - self.air_gap_h[i + 2] / 2 - self.air_gap_h[i + 1] / 2) / 2
                            print('Upper air gap detected')
                        else:
                            h1 = (self.position[i] - self.position[i - 1] - self.air_gap_h[i] / 2 - self.air_gap_h[
                                i - 1] / 2) / 2
                            h2 = (self.position[i + 1] - self.position[i] - self.air_gap_h[i + 1] / 2 - self.air_gap_h[
                                i] / 2) / 2
                            print('Both air gap detected')

                        self.reluctance[:, 5] = self.reluctance[:, 5] + fr.r_air_gap_round_round([self.air_gap_h[i]], [
                            self.core_inner_diameter], [h1], [h2])

    def air_gap_reluctance_single_new(self):
        """
        Calculate the air-gap reluctance for a single simulation with single/multiple air-gaps.

        Method is according to the following paper:
        ["A Novel Approach for 3D Air Gap Reluctance Calculations" - J. Mühlethaler, J.W. Kolar, A. Ecklebe]

        It is called when input variable sim_type == 'single'.
        Its calculation for multiple air-gap is based on superposition.
        That is, multiple air-gap's reluctance is calculated by taking one at a time and then adding them together
        (like in superposition theorem)
        """
        # Conversion of different air_gap_method to manual method
        if self.air_gap_method == "Percent":
            self.air_gap_position = (np.array(self.air_gap_position) / 100) * self.window_h[0]
        elif self.air_gap_method == "Center":
            self.air_gap_position = [self.window_h[0] / 2]

        core_height_upper = self.window_h - self.air_gap_position - (self.air_gap_h / 2)
        core_height_lower = self.air_gap_position - (self.air_gap_h / 2)

        core_height_upper = np.where(core_height_upper <= 0, 0, core_height_upper)
        core_height_lower = np.where(core_height_lower <= 0, 0, core_height_lower)

        self.reluctance[:, 5] = np.sum(np.where(core_height_upper * core_height_lower == 0,
                                       fr.r_air_gap_round_inf(self.air_gap_h, self.core_inner_diameter, core_height_lower + core_height_upper),
                                       fr.r_air_gap_round_round(self.air_gap_h, self.core_inner_diameter, core_height_upper, core_height_lower)))

    def calculate_inductance(self):
        """Calculate the inductance from Number of turns and total reluctance (L = N^2 / R_m)."""
        self.cal_inductance = (self.no_of_turns * self.no_of_turns) / np.sum(self.reluctance, axis=1)
        self.data_matrix[:, 9] = self.cal_inductance
        print(f"Inductance:{self.cal_inductance}")

    def add_column_to_data_matrix(self, data_matrix, column_value, column_name: str):
        """
        Add column to the given matrix.

        :param data_matrix: Matrix containing the design parameters
        :type data_matrix: ndarray
        :param column_value: Column to be added
        :type column_value: ndarray
        :param column_name: Identifier of the column
        :type column_name: str
        """
        size = len(data_matrix[0])
        data_matrix = np.hstack((data_matrix, np.reshape(column_value, (len(column_value), 1))))
        self.param_pos_dict[column_name] = size

        return data_matrix


if __name__ == '__main__':
    mc1 = MagneticCircuit(core_inner_diameter=[0.0149], window_h=[0.0295], window_w=[0.01105], no_of_turns=[1],
                          n_air_gaps=[1],
                          air_gap_h=[0.0005], air_gap_position=[50], mu_r_abs=[3000],
                          mult_air_gap_type=[1, 2],
                          air_gap_method='Percent', component_type='inductor', sim_type='single')  # 0.0149

    # mc1 = MagneticCircuit(core_inner_diameter=[0.0149], window_h=[0.0295], window_w=[0.01105], no_of_turns=[9],
    #                       n_air_gaps=[1, 3],
    #                       air_gap_h=[0.0005], air_gap_position=[0, 50, 100], mu_r_abs=[3000],
    #                       mult_air_gap_type=[1, 2],
    #                       air_gap_method='Percent', component_type='inductor', sim_type='sweep')  # 0.0149
    # plot_r_basic()
