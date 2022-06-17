# 2D-axis symmetric core reluctance calculations
import femmt as fmt
import numpy as np
import matplotlib.pyplot as plt
import schemdraw
import schemdraw.elements as elm
from femmt.femmt_functions import *


def set_orientation(section: list, n: int) -> tuple:    # Function to dynamically define circuit component orientation: 'up', 'right', 'down', 'left'
    temp1 = ['up', 'right', 'down', 'left']
    temp2 = []
    if n % 2 != 0:
        section.append('l')
        n += 1
    if n % 2 == 0:
        if n % 4 == 0:
            for i in range(4):
                for j in range(n // 4):
                    temp2.append(temp1[i])
        else:
            for i in range(4):
                for j in range((n-2) // 4):
                    temp2.append(temp1[i])

            temp2.insert((n-2) // 4, 'right')
            temp2.append('left')
    return section, temp2


def plot_r_basis():
    width = 1
    length = 1
    height = np.linspace(10, 0.1, 1000)
    h_l = height / length

    r_m = 1 / (mu0 * (width / 2 / length + 2 / np.pi * (
                1 + np.log(np.pi * height / 4 / length))))

    combined = np.vstack((h_l, r_m)).T
    print(combined)
    fig, ax = plt.subplots()  # Create a figure containing a single axes.
    plt.title("R_basic vs h/l")
    plt.xlabel("h/l")
    plt.ylabel("R_basic")
    ax.plot(h_l, r_m)
    ax.invert_xaxis()
    ax.grid()
    plt.show()


class MagneticCircuit:

    """This is a class for calculating the reluctance and inductance and visualising magnetic circuit"""

    def __init__(self, core_h, core_w, window_h, window_w, r_outer,
                 no_of_turns, current, method, n_air_gaps, air_gap_h, air_gap_position):

        self.core_h = core_h
        self.core_w = core_w
        self.window_h = window_h
        self.window_w = window_w
        self.r_outer = r_outer

        self.no_of_turns = no_of_turns
        self.current = current
        self.method = method
        self.n_air_gaps = n_air_gaps
        self.air_gap_h = air_gap_h
        self.air_gap_position = air_gap_position

        self.middle_h = None   # height of upper and lower part of the window in the core
        self.outer_w = None    # Outer leg width
        self.L = None          # Total inductance
        self.section = None
        self.orientation = None
        self.length = None
        self.area = None
        self.reluctance = None
        self.mu_rel = 3100
        self.mu_0 = 4 * np.pi * 1e-7

        self.max_percent = None
        self.min_percent = None
        self.position = None

    # position_tag = 0 for 2D axis symmetric cores
    def core_reluctance(self):
        self.middle_h = (self.core_h - self.window_h)/2
        r_inner = self.core_w/2 + self.window_w
        self.outer_w = self.r_outer - r_inner

        self.section = [0, 1, 2, 3, 4, 5, 4, 3, 2]
        self.orientation = []
        self.length = np.zeros(len(self.section))
        self.area = np.zeros(len(self.section))
        self.reluctance = np.zeros(len(self.section))

        for i in range(len(self.section)):
            if self.section[i] == 0:
                self.reluctance[i] = self.no_of_turns * self.current

            if self.section[i] == 1:
                self.length[i] = self.window_h - sum(self.air_gap_h)
                self.area[i] = np.pi * ((self.core_w/2) ** 2)
                self.reluctance[i] = self.length[i] / (self.mu_0 * self.mu_rel * self.area[i])

            elif self.section[i] == 2:
                self.length[i] = (np.pi / 8) * (self.core_w / 2 + self.middle_h)
                self.area[i] = ((self.core_w / 2 + self.middle_h) / 2) * 2 * np.pi * (self.core_w / 2)
                self.reluctance[i] = self.length[i] / (self.mu_0 * self.mu_rel * self.area[i])

            elif self.section[i] == 3:
                self.length[i] = self.window_w
                self.area[i] = np.nan
                self.reluctance[i] = ((self.mu_0 * self.mu_rel * 2 * np.pi * self.middle_h) ** -1) * np.log((2 * r_inner)/self.core_w)

            elif self.section[i] == 4:
                self.length[i] = (np.pi / 8) * (self.outer_w + self.middle_h)
                self.area[i] = ((self.outer_w + self.middle_h) / 2) * 2 * np.pi * r_inner
                self.reluctance[i] = self.length[i] / (self.mu_0 * self.mu_rel * self.area[i])

            elif self.section[i] == 5:
                self.length[i] = self.window_h
                self.area[i] = np.pi * (self.r_outer ** 2 - r_inner ** 2)
                self.reluctance[i] = self.length[i] / (self.mu_0 * self.mu_rel * self.area[i])

    def air_gap_reluctance(self):
        flag_0 = 0
        flag_1 = 0
        flag_2 = 0
        if self.method == 'center':
            self.section.append(6)  # round-round
            temp1 = r_basis(self.air_gap_h[0] / 2, self.core_w / 2, (self.window_h - self.air_gap_h[0]) / 2)
            temp2 = sigma(self.air_gap_h[0], self.core_w / 2, 2 * temp1)
            temp3 = fmt.femmt_functions.r_round_round(self.air_gap_h[0], temp2, self.core_w / 2)
            self.reluctance = np.append(self.reluctance, temp3)

        elif self.method == 'percent':
            self.max_percent = ((self.window_h - self.air_gap_h[self.n_air_gaps - 1] / 2) / self.window_h) * 100
            self.min_percent = ((self.air_gap_h[0] / 2) / self.window_h) * 100
            self.position = np.array(self.air_gap_position) / 100 * self.window_h  # Convert percent to absolute value

            if self.air_gap_position[0] <= self.min_percent:
                flag_0 = 1
                self.section.append(8)
                if self.n_air_gaps == 1:
                    h = self.window_h - self.air_gap_h[0]
                else:
                    h = (self.position[1] - self.position[0]) / 2

                temp1 = r_basis(self.air_gap_h[0], self.core_w / 2, h)
                temp2 = sigma(self.air_gap_h[0], self.core_w / 2, temp1)
                temp3 = fmt.femmt_functions.r_round_inf(self.air_gap_h[0], temp2, self.core_w / 2)
                self.reluctance = np.append(self.reluctance, temp3)
                print('air gap is at lower corner')

            if self.air_gap_position[self.n_air_gaps - 1] >= self.max_percent:
                flag_1 = 1
                self.section.append(8)
                if self.n_air_gaps == 1:
                    h = self.window_h - self.air_gap_h[0]
                else:
                    h = (self.position[self.n_air_gaps - 1] - self.position[self.n_air_gaps - 2]) / 2

                temp1 = r_basis(self.air_gap_h[0], self.core_w / 2, h)
                temp2 = sigma(self.air_gap_h[0], self.core_w / 2, temp1)
                temp3 = fmt.femmt_functions.r_round_inf(self.air_gap_h[0], temp2, self.core_w / 2)
                self.reluctance = np.append(self.reluctance, temp3)
                print('air gap is at upper corner')

            for i in range(self.n_air_gaps):
                if self.min_percent < self.air_gap_position[i] < self.max_percent:
                    self.section.append(7)
                    if flag_2 == 0:
                        if flag_0 == 0 and flag_1 == 0:
                            self.position = np.append(self.position, self.window_h + (self.window_h - self.position[self.n_air_gaps - 1]))
                            self.position = np.insert(self.position, 0, self.position[0])
                        elif flag_0 == 1 and flag_1 == 0:
                            self.position = np.append(self.position, self.window_h + (self.window_h - self.position[self.n_air_gaps - 1]))
                        elif flag_0 == 0 and flag_1 == 1:
                            self.position = np.insert(self.position, 0, self.position[0])
                        flag_2 = 1

                    if flag_0 == 0 and flag_1 == 0:
                        h1 = (self.position[i + 1] - self.position[i]) / 2
                        h2 = (self.position[i + 2] - self.position[i + 1]) / 2
                        print('No corner air gap detected')
                    elif flag_0 == 1 and flag_1 == 0:
                        h1 = (self.position[i] - self.position[i - 1]) / 2
                        h2 = (self.position[i + 1] - self.position[i]) / 2
                        print('Lower air gap detected')
                    elif flag_0 == 0 and flag_1 == 1:
                        h1 = (self.position[i + 1] - self.position[i]) / 2
                        h2 = (self.position[i + 2] - self.position[i + 1]) / 2
                        print('Upper air gap detected')
                    else:
                        h1 = (self.position[i] - self.position[i - 1]) / 2
                        h2 = (self.position[i + 1] - self.position[i]) / 2
                        print('Both air gap detected')

                    r_basis_1 = r_basis(self.air_gap_h[0] / 2, self.core_w / 2, h1)
                    r_basis_2 = r_basis(self.air_gap_h[0] / 2, self.core_w / 2, h2)
                    temp2 = sigma(self.air_gap_h[0], self.core_w / 2, r_basis_1 + r_basis_2)
                    temp3 = fmt.femmt_functions.r_round_round(self.air_gap_h[0], temp2, self.core_w / 2)
                    self.reluctance = np.append(self.reluctance, temp3)

            self.section, self.orientation = set_orientation(self.section, len(self.section))
            self.L = (self.no_of_turns * self.no_of_turns) / sum(self.reluctance)

            print('Section:', self.section)
            print(self.max_percent)
            print(self.min_percent)
            print('Orientation:', self.orientation)
            print('List of length:', self.length)
            print('List of cross-section area:', self.area)
            print('List of reluctance:', self.reluctance)
            print('Inductance:', self.L)

    def draw_schematic(self):
        d = schemdraw.Drawing()
        for i in range(len(self.section)):
            if self.section[i] == 0:
                d += getattr(elm.SourceV().dot().label(str(round(self.reluctance[i], 2)) + 'AT'), self.orientation[i])
            elif self.section[i] == 'l':
                d += getattr(elm.Line().dot(), self.orientation[i])
            elif self.section[i] == 6 or self.section[i] == 7 or self.section[i] == 8:
                d += getattr(elm.ResistorIEC().dot().label(str(round(self.reluctance[i], 2)) + 'AT/Wb'), self.orientation[i])
            else:
                d += getattr(elm.Resistor().dot().label(str(round(self.reluctance[i], 2)) + 'AT/Wb'), self.orientation[i])

        d.draw()
        d.save('my_circuit.svg')


mc1 = MagneticCircuit(0.0398, 0.0149, 0.0295, 0.01105, 0.01994, 8, 3, 'percent', 2, [0.0005, 0.0005], [1, 70])
mc1.core_reluctance()
mc1.air_gap_reluctance()
mc1.draw_schematic()
# plot_r_basis()

