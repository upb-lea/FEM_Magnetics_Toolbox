# 2D-axis symmetric core reluctance calculations


import femmt as fmt
from femmt.femmt_functions import *
import numpy as np
import schemdraw
import schemdraw.elements as elm


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


def basic_example_func(f_method, f_n, f_h, f_pos, f_n_turns, f_core_cond_iso, f_current):
    # 2. set core parameters
    core = fmt.core_database()["PQ 40/40"]

    geo.core.update(core_w=core["core_w"], window_w=core["window_w"], window_h=core["window_h"],
                    mu_rel=3100, phi_mu_deg=12,
                    sigma=0.6)

    # 3. set air gap parameters
    geo.air_gaps.update(method=f_method, n_air_gaps=f_n, air_gap_h=f_h, air_gap_position=f_pos,
                        position_tag=[0])

    geo.update_conductors(n_turns=[f_n_turns], conductor_type=["solid"], conductor_radii=[0.0015],
                          winding=["primary"], scheme=["square"],
                          core_cond_isolation=[0.001, 0.001, f_core_cond_iso, 0.001], cond_cond_isolation=[0.0001],
                          conductivity_sigma=["copper"])

    # geo.update_conductors(n_turns=[f_n_turns], conductor_type=["litz"], conductor_radii=[0.0015],
    #                       litz_para_type=['implicit_ff'],
    #                       ff=[None], strands_numbers=[100], strand_radii=[70e-6],
    #                       winding=["primary"], scheme=["square"],
    #                       core_cond_isolation=[0.001, 0.001, f_core_cond_iso, 0.001], cond_cond_isolation=[0.0001],
    #                       conductivity_sigma=["copper"])

    geo.create_model(freq=100000, visualize_before=False, do_meshing=True, save_png=False)

    geo.single_simulation(freq=100000, current=f_current, phi_deg=[0], show_results=False)


class MagneticCircuit:

    """This is a class for calculating the reluctance and inductance and visualising magnetic circuit"""

    # position_tag = 0 for 2D axis symmetric cores
    def __init__(self, core_w, window_h, window_w,
                 no_of_turns, current, method, n_air_gaps, air_gap_h, air_gap_position):

        self.core_h = window_h + core_w / 2
        self.core_w = core_w
        self.window_h = window_h
        self.window_w = window_w
        self.r_outer = None
        self.r_inner = None
        self.core_h_middle = None  # height of upper and lower part of the window in the core
        self.outer_w = None  # Outer leg width

        self.no_of_turns = no_of_turns
        self.current = current
        self.method = method
        self.n_air_gaps = n_air_gaps
        self.air_gap_h = air_gap_h
        self.air_gap_position = air_gap_position

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

    def core_reluctance(self):
        self.core_h_middle = (self.core_h - self.window_h)/2
        self.r_inner = self.core_w/2 + self.window_w
        self.r_outer = np.sqrt((self.core_w / 2) ** 2 + self.r_inner ** 2)
        self.outer_w = self.r_outer - self.r_inner

        self.section = [0, 1, 2, 3, 4]
        self.length = np.zeros(len(self.section))
        self.area = np.zeros(len(self.section))
        self.reluctance = np.zeros(len(self.section))

        self.length[0] = self.window_h - sum(self.air_gap_h)
        self.area[0] = np.pi * ((self.core_w/2) ** 2)

        self.length[1] = (np.pi / 8) * (self.core_w / 2 + self.core_h_middle)
        self.area[1] = ((self.core_w / 2 + self.core_h_middle) / 2) * 2 * np.pi * (self.core_w / 2)

        self.length[2] = self.window_w
        self.area[2] = np.nan
        self.reluctance[2] = ((self.mu_0 * self.mu_rel * 2 * np.pi * self.core_h_middle) ** -1) * np.log((2 * self.r_inner) / self.core_w)

        self.length[3] = (np.pi / 8) * (self.outer_w + self.core_h_middle)
        self.area[3] = ((self.outer_w + self.core_h_middle) / 2) * 2 * np.pi * self.r_inner

        self.length[4] = self.window_h
        self.area[4] = np.pi * (self.r_outer ** 2 - self.r_inner ** 2)

        self.reluctance[~np.isnan(self.area)] = self.length[~np.isnan(self.area)] / (self.mu_0 * self.mu_rel * self.area[~np.isnan(self.area)])
        self.reluctance[1:3] = 2 * self.reluctance[1:3]

    def air_gap_reluctance(self):
        flag_0 = 0
        flag_1 = 0
        flag_2 = 0
        if self.method == 'center':
            self.section.append(6)  # round-round
            temp1 = r_basis(self.air_gap_h[0] / 2, self.core_w, (self.window_h - self.air_gap_h[0]) / 2)
            temp2 = sigma(self.air_gap_h[0], self.core_w / 2, 2 * temp1)
            temp3 = fmt.femmt_functions.r_round_round(self.air_gap_h[0], temp2, self.core_w / 2)
            temp4 = self.air_gap_h[0] / (self.mu_0 * np.pi * (self.core_w / 2) ** 2)    # classical reluctance formula
            self.reluctance = np.append(self.reluctance, temp3)
            # print(f"New approach reluctance: {temp3}")
            # print(f"Classical reluctance: {temp4}")

        elif self.method == 'percent':
            self.max_percent = ((self.window_h - self.air_gap_h[self.n_air_gaps - 1] / 2) / self.window_h) * 100
            self.min_percent = ((self.air_gap_h[0] / 2) / self.window_h) * 100
            self.position = np.array(self.air_gap_position) / 100 * self.window_h  # Convert percent to absolute value
            print(f"Max percent: {self.max_percent}")
            print(f"Min percent: {self.min_percent}")

            if self.air_gap_position[0] <= self.min_percent:
                flag_0 = 1
                self.section.append(8)
                if self.n_air_gaps == 1:
                    h = self.window_h - self.air_gap_h[0]
                else:
                    h = (self.position[1] - self.position[0]) / 2

                temp1 = r_basis(self.air_gap_h[0], self.core_w, h)
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

                temp1 = r_basis(self.air_gap_h[0], self.core_w, h)
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
                            self.position = np.insert(self.position, 0, -self.position[0])
                        elif flag_0 == 1 and flag_1 == 0:
                            self.position = np.append(self.position, self.window_h + (self.window_h - self.position[self.n_air_gaps - 1]))
                        elif flag_0 == 0 and flag_1 == 1:
                            self.position = np.insert(self.position, 0, -self.position[0])
                        flag_2 = 1

                    if flag_0 == 0 and flag_1 == 0:
                        h1 = (self.position[i + 1] - self.position[i]) / 2
                        h2 = (self.position[i + 2] - self.position[i + 1]) / 2
                        # h1 = self.position[i + 1] - (self.air_gap_h[i] / 2)
                        # h2 = self.window_h - (self.position[i + 1] + (self.air_gap_h[i] / 2))
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

                    r_basis_1 = r_basis(self.air_gap_h[0] / 2, self.core_w, h1)
                    r_basis_2 = r_basis(self.air_gap_h[0] / 2, self.core_w, h2)
                    temp2 = sigma(self.air_gap_h[0], self.core_w / 2, r_basis_1 + r_basis_2)
                    temp3 = fmt.femmt_functions.r_round_round(self.air_gap_h[0], temp2, self.core_w / 2)
                    self.reluctance = np.append(self.reluctance, temp3)

        # self.section, self.orientation = set_orientation(self.section, len(self.section))
        self.L = (self.no_of_turns * self.no_of_turns) / sum(self.reluctance)

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


# Sweep of air-gap and winding position and compare it with FEM simulation
sweep_air_gap_h = np.linspace(20, 80, 4)
sweep_wndg_pos = np.linspace(0.001, 0.007, 7)
# sweep_current = np.linspace(0.1, 10, 5)
fem_ind = np.zeros((len(sweep_air_gap_h), len(sweep_wndg_pos)))
cal_ind = np.zeros((len(sweep_air_gap_h), len(sweep_wndg_pos)))

# 1. chose simulation type
geo = fmt.MagneticComponent(component_type="inductor")

for j in range(len(sweep_air_gap_h)):
    for i in range(len(sweep_wndg_pos)):
        mc1 = MagneticCircuit(0.0149, 0.0295, 0.01105, 8, 10, 'percent', 1, [0.001], [sweep_air_gap_h[j]])  # 0.0149
        mc1.core_reluctance()
        mc1.air_gap_reluctance()
        cal_ind[j, i] = mc1.L
        # mc1.max_percent = ((mc1.window_h - (sweep_air_gap_h[j] / 2)) / mc1.window_h) * 100
        # mc1.min_percent = ((sweep_air_gap_h[j] / 2) / mc1.window_h) * 100

        basic_example_func("percent", 1, [0.001], [sweep_air_gap_h[j]], [8], sweep_wndg_pos[i], [10])
        fem_ind[j, i] = geo.read_log()["single_sweeps"][0]["winding1"]["self_inductivity"][0]
        # fem_ind[j, i] = geo.read_log()["single_sweeps"][0]["winding1"]["Q"]

print(f"Air-gap length: {sweep_air_gap_h}")
print(f"Winding position: {sweep_wndg_pos}")
print(f"FEM inductance: {fem_ind}")
print(f"Calculated inductance: {cal_ind}")

h_by_l = ((0.0295 - sweep_air_gap_h) / 2) / sweep_air_gap_h
error = ((fem_ind - cal_ind) / fem_ind) * 100

# Plotting tools
fig, ax = plt.subplots()  # Create a figure containing a single axes.
plt.title("inductance vs winding pos (8 conductor only) (variable airgap pos) (1mm)")
plt.xlabel("Winding position (in m)")
plt.ylabel("error (in %)")
for j in range(len(sweep_air_gap_h)):
    ax.plot(sweep_wndg_pos, error[j, :], label=str(sweep_air_gap_h[j]))
ax.legend(loc='best')
ax.grid()
plt.show()

# mc1 = MagneticCircuit(0.0398, 0.0149, 0.0295, 0.01105, 8, 3, 'center', 1, [0.0005], [50])
# mc1 = MagneticCircuit(0.0149, 0.0295, 0.01105, 8, 3, 'center', 1, [0.02], [50])
# mc1 = MagneticCircuit(0.0149, 0.0295, 0.01105, 8, 3, 'percent', 2, [0.0005, 0.0005], [20, 50])
# mc1.core_reluctance()
# mc1.air_gap_reluctance()
# mc1.draw_schematic()
# plot_error()