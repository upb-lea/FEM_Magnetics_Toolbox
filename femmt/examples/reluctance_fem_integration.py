import femmt as fmt
import numpy as np
import re
import os
from os import listdir
from os.path import isfile, join
import shutil
import matplotlib.pyplot as plt
from itertools import product
import logging
from scipy.interpolate import interp1d
from mpl_toolkits.mplot3d import proj3d
import materialdatabase as mdb

material_db = mdb.MaterialDatabase()


def plot_2d(x_value: list, y_value: list, x_label: str, y_label: str, title: str, annotations: list):
    """
        Visualize data in 2d plot with popover next to mouse position.

        param x_value: Data points for x-axis
        :type x_value: list
        :param y_value: Data points for y-axis
        :type y_value: list
        :param x_label: x-axis label
        :type x_label: str
        :param y_label: y-axis label
        :type y_label: str
        :param title: Title of the graph
        :type title: str
        :param annotations: Annotations corresponding to the 3D points
        :type annotations: list
        """
    names = np.array(annotations)
    x_value_str = [str(round(x, 6)) for x in x_value]
    y_value_str = [str(round(y, 6)) for y in y_value]

    c = np.random.randint(1, 5, size=len(y_value))

    norm = plt.Normalize(1, 4)
    cmap = plt.cm.RdYlGn

    fig, ax = plt.subplots()
    fmt.plt.title(title)
    fmt.plt.xlabel(x_label)
    fmt.plt.ylabel(y_label)
    sc = plt.scatter(x_value, y_value, c=c, s=50, cmap=cmap, norm=norm)

    annot = ax.annotate("", xy=(0, 0), xytext=(20, 20), textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)

    def update_annot(ind):
        """Create popover annotations in 2d plot"""

        pos = sc.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        text = "{}\nVolume: {}\nLoss:{}".format(" ".join([names[n] for n in ind["ind"]]),
                                                " ".join([x_value_str[n] for n in ind["ind"]]),
                                                " ".join([y_value_str[n] for n in ind["ind"]]))
        annot.set_text(text)
        annot.get_bbox_patch().set_facecolor(cmap(norm(c[ind["ind"][0]])))
        annot.get_bbox_patch().set_alpha(0.4)

    def hover(event):
        """Event that is triggered when mouse is hovered.
        Shows text annotation over data point closest to mouse."""
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = sc.contains(event)
            if cont:
                update_annot(ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

    fig.canvas.mpl_connect("motion_notify_event", hover)
    ax.grid()
    plt.show()


def plot_3d(x_value: list, y_value: list, z_value: list, x_label: str, y_label: str, z_label: str, x_limit: list,
            y_limit: list, z_limit: list, title: str, annotations: list):
    """
    Visualize data in 3d plot with popover next to mouse position.

    param x_value: Data points for x-axis
    :type x_value: list
    :param y_value: Data points for y-axis
    :type y_value: list
    :param z_value: Data points for z-axis
    :type z_value: list
    :param x_label: x-axis label
    :type x_label: str
    :param y_label: y-axis label
    :type y_label: str
    :param z_label: z-axis label
    :type z_label: str
    :param x_limit: Min and max limit of x-axis
    :type x_limit: list
    :param y_limit: Min and max limit of y-axis
    :type y_limit: list
    :param z_limit: Min and max limit of z-axis
    :type z_limit: list
    :param title: Title of the graph
    :type title: str
    :param annotations: Annotations corresponding to the 3D points
    :type annotations: list
    """
    names = np.array(annotations)
    num = [re.findall(r'\d+', item) for item in names]
    case_num_list = [int(item[0]) for item in num]

    X = np.zeros((len(x_value), 1))
    X = np.hstack((X, np.reshape(x_value, (len(x_value), 1))))
    X = np.hstack((X, np.reshape(y_value, (len(y_value), 1))))
    X = np.hstack((X, np.reshape(z_value, (len(z_value), 1))))
    X = np.hstack((X, np.reshape(case_num_list, (len(case_num_list), 1))))
    X = np.delete(X, 0, 1)
    X = X[X[:, 3].argsort()]

    x_value_str = [str(round(x, 6)) for x in X[:, 0]]
    y_value_str = [str(round(y, 6)) for y in X[:, 1]]
    z_value_str = [str(round(z, 6)) for z in X[:, 2]]

    fig = plt.figure(figsize=(16, 10))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_zlabel(z_label)
    ax.set_xlim(xmin=x_limit[0], xmax=x_limit[1])
    ax.set_ylim(ymin=y_limit[0], ymax=y_limit[1])
    ax.set_zlim(zmin=z_limit[0], zmax=z_limit[1])
    XX = X[:, 0]
    XX[XX > x_limit[1]] = np.nan
    YY = X[:, 1]
    YY[YY > y_limit[1]] = np.nan
    ZZ = X[:, 2]
    ZZ[ZZ > z_limit[1]] = np.nan
    X[:, 0] = XX
    X[:, 1] = YY
    X[:, 2] = ZZ
    print(X)
    c = np.random.randint(1, 5, size=len(X))

    norm = plt.Normalize(1, 4)
    cmap = plt.cm.RdYlGn

    ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=c, s=10, cmap=cmap, norm=norm, depthshade=False, picker=True)

    def distance(point, event):
        """Return distance between mouse position and given data point

        Args:
            point (np.array): np.array of shape (3,), with x,y,z in data coords
            event (MouseEvent): mouse event (which contains mouse position in .x and .xdata)
        Returns:
            distance (np.float64): distance (in screen coords) between mouse pos and data point
        """
        assert point.shape == (3,), "distance: point.shape is wrong: %s, must be (3,)" % point.shape

        # Project 3d data space to 2d data space
        x2, y2, _ = proj3d.proj_transform(point[0], point[1], point[2], plt.gca().get_proj())
        # Convert 2d data space to 2d screen space
        x3, y3 = ax.transData.transform((x2, y2))

        return np.sqrt((x3 - event.x) ** 2 + (y3 - event.y) ** 2)

    def calcClosestDatapoint(X, event):
        """Calculate which data point is closest to the mouse position.

        Args:
            X (np.array) - array of points, of shape (numPoints, 3)
            event (MouseEvent) - mouse event (containing mouse position)
        Returns:
            smallestIndex (int) - the index (into the array of points X) of the element closest to the mouse position
        """
        X_modified = X[np.all(~np.isnan(X), axis=1), 0:3]
        distances = [distance(X_modified[i, 0:3], event) for i in range(X_modified.shape[0])]
        return np.argmin(distances)

    def annotatePlot(X, index):
        """Create popover label in 3d chart

        Args:
            X (np.array) - array of points, of shape (numPoints, 3)
            index (int) - index (into points array X) of item which should be printed
        Returns:
            None
        """
        # If we have previously displayed another label, remove it first
        if hasattr(annotatePlot, 'label'):
            annotatePlot.label.remove()
        # Get data point from array of points X, at position index
        x2, y2, _ = proj3d.proj_transform(X[index, 0], X[index, 1], X[index, 2], ax.get_proj())
        annotatePlot.label = plt.annotate(
            f'case{index}\nVolume:{x_value_str[index]}\nLoss:{y_value_str[index]}\nCost:{z_value_str[index]}',
            xy=(x2, y2), xytext=(-20, 20), textcoords='offset points', ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
            arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
        fig.canvas.draw()

    def onMouseMotion(event):
        """Event that is triggered when mouse is moved. Shows text annotation over data point closest to mouse."""
        closestIndex = calcClosestDatapoint(X, event)
        annotatePlot(X, closestIndex)

    fig.canvas.mpl_connect('motion_notify_event', onMouseMotion)  # on mouse motion
    plt.show()


def load_design(working_directory: str):
    """
    Load FEM simulation results from given working directory

    param working_directory: Sets the working directory
    :type working_directory: str
    """
    working_directories = []
    labels = []
    working_directory = os.path.join(working_directory, 'fem_simulation_data')
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
        total_loss.append(data.sweeps[0].windings[0].active_power)
        total_volume.append(data.core_2daxi_total_volume)
        total_cost.append(data.total_cost)

    real_inductance = []
    for i in range(len(total_loss)):
        real_inductance.append(inductivities[i].real)
    print(f'labels: {labels}')
    print(total_volume)
    print(total_loss)
    print(total_cost)

    return real_inductance, total_loss, total_volume, total_cost, labels


class AutomatedDesign:
    """
    AutomatedDesign class implements brute force optimization for magnetic component design.
    It consists of input parameters sweep, filtration, FEM simulation and plotting of relevant results.
    """

    # copper conductivity (sigma) @ 20 degree celsius (in siemens/meter)
    copper_conductivity = 5.96 * 1e7
    winding_factor_dict = {'Square': 0.785, 'Hexagonal': 0.907}
    winding_scheme_dict = {'Square': fmt.ConductorArrangement.Square, 'Hexagonal': fmt.ConductorArrangement.Hexagonal}

    def __init__(self, working_directory: str, component: str, goal_inductance: float, frequency: float,
                 inductance_percent_tolerance: int, winding_scheme: str, current_max: float, percent_of_b_sat: int,
                 percent_of_total_loss: int, database_core_names: list, database_litz_names: list,
                 solid_conductor_r: list, manual_core_w: list, manual_window_h: list, manual_window_w: list,
                 no_of_turns: list, n_air_gaps: list, air_gap_height: list, air_gap_position: list, material_list: list,
                 mult_air_gap_type: list, top_core_insulation: float, bot_core_insulation: float,
                 left_core_insulation: float, right_core_insulation: float, inner_winding_insulations: list,
                 temperature: float):

        """
                :param working_directory: Sets the working directory
                :type working_directory: str
                :param component: Sets magnetic component: 'inductor', 'integrated transformer'
                :type component: str
                :param goal_inductance: Sets goal inductance for design [in Henry]
                :type goal_inductance: float
                :param frequency: Operating frequency [in Hz]
                :type frequency: float
                :param inductance_percent_tolerance: percent tolerance with respect to goal inductance [in percent]
                :type inductance_percent_tolerance: float
                :param winding_scheme: Winding scheme: 'Square' or 'Hexagonal'
                :type winding_scheme: str
                :param current_max: Max current amplitude with assumption of sinusoidal current waveform [in Ampere]
                :type current_max: float
                :param percent_of_b_sat: percent of saturation of magnetic flux density [in percent]
                :type percent_of_b_sat: int
                :param percent_of_total_loss: percentage of the total loss [in percent]
                :type percent_of_total_loss: int
                :param database_core_names: list of core names from the database
                :type database_core_names: list
                :param database_litz_names: list of litz wire names from the database
                :type database_litz_names: list
                :param solid_conductor_r: Solid conductor radius [in meter]
                :type solid_conductor_r: list
                :param manual_core_w: Diameter of center leg of the core [in meter]
                :type manual_core_w: list
                :param manual_window_h: Height of the core window [in meter]
                :type manual_window_h: list
                :param manual_window_w: Width of the core window [in meter]
                :type manual_window_w: list
                :param no_of_turns: Number of turns
                :type no_of_turns: list
                :param n_air_gaps: Number of air-gaps in the center leg of the core
                :type n_air_gaps: list
                :param air_gap_height: Air-gap height [in meter]
                :type air_gap_height: list
                :param air_gap_position: Position of the air-gap in the percentage with respect to window_h [in percent]
                :type air_gap_position: list
                :param material_list: Relative permeability of the core [in F/m]
                :type material_list: list
                :param mult_air_gap_type: Two types of equally distributed air-gaps (used only for air-gaps more than 1)
                    Type 1: Edge distributed
                    Type 2: Center distributed
                :type mult_air_gap_type: list
                :param top_core_insulation: top_core_insulation [in meter]
                :type top_core_insulation: float
                :param bot_core_insulation: bot_core_insulation [in meter]
                :type bot_core_insulation: float
                :param left_core_insulation: left_core_insulation [in meter]
                :type left_core_insulation: float
                :param right_core_insulation: right_core_insulation [in meter]
                :type right_core_insulation: float
                :param inner_winding_insulations: inner_winding_insulations [in meter]
                :type inner_winding_insulations: float
                :param temperature: core temperature [in degree Celsius]
                :type temperature: float
        """

        self.working_directory = working_directory
        if not os.path.exists(self.working_directory):
            os.mkdir(self.working_directory)
        self.component = component

        self.goal_inductance = goal_inductance
        self.inductance_percent_tolerance = inductance_percent_tolerance
        self.current_max = current_max
        self.percent_of_b_sat = percent_of_b_sat
        self.percent_of_total_loss = percent_of_total_loss
        self.frequency = frequency
        self.temperature = temperature

        # Set core-geometry from core database or/and manual entry
        self.db_core_names = database_core_names
        self.manual_core_w = manual_core_w
        self.manual_window_h = manual_window_h
        self.manual_window_w = manual_window_w

        # Set winding settings (Solid and Litz winding type)
        self.db_litz_names = database_litz_names
        self.solid_conductor_r = solid_conductor_r

        # Set air-gap and core parameters
        self.no_of_turns = no_of_turns
        self.n_air_gaps = n_air_gaps
        self.air_gap_height = air_gap_height
        self.air_gap_position = air_gap_position

        # Set core material
        self.material_list = material_list

        # Multiple air-gap type ('edge_distributed' and/or 'center_distributed')
        self.mult_air_gap_type = mult_air_gap_type

        # Set windng and insulation data
        self.winding_scheme = winding_scheme
        self.winding_factor = self.winding_factor_dict[winding_scheme]
        self.top_core_insulation = top_core_insulation
        self.bot_core_insulation = bot_core_insulation
        self.left_core_insulation = left_core_insulation
        self.right_core_insulation = right_core_insulation
        self.inner_winding_insulations = inner_winding_insulations

        self.core_w_list, self.window_h_list, self.window_w_list, self.litz_conductor_r, self.litz_strand_r, \
        self.litz_strand_n, self.mu_rel, self.mult_air_gap_type_list = self.input_pre_process()

        self.min_conductor_r = min(self.litz_conductor_r + self.solid_conductor_r)
        self.conductor_r_list = self.litz_conductor_r + self.solid_conductor_r

        # Call to Reluctance model (Class MagneticCircuit)
        mc = fmt.MagneticCircuit(core_w=self.core_w_list, window_h=self.window_h_list, window_w=self.window_w_list,
                                 no_of_turns=self.no_of_turns, n_air_gaps=self.n_air_gaps,
                                 air_gap_h=self.air_gap_height, air_gap_position=self.air_gap_position,
                                 mu_rel=self.mu_rel, mult_air_gap_type=self.mult_air_gap_type_list,
                                 air_gap_method='percent', component_type=self.component, sim_type='sweep')
        self.param = mc.get_param_pos_dict()

        # Filtration of the design cases which are not important
        data_matrix_1 = self.filter1_geometry(mc.data_matrix)
        data_matrix_2 = self.filter2_inductance(data_matrix_1)
        data_matrix_3 = self.filter3_flux_saturation(data_matrix_2)
        data_matrix_4 = self.filter4_losses(data_matrix_3)
        self.plot_volume_loss(data_matrix_4)
        self.data_matrix_fem = data_matrix_4
        #     # FEM_data_matrix = data_matrix_4[np.where(data_matrix_4[:, param["normalized_total_loss"]]
        #     <= ((0.01 / data_matrix_4[:, param["normalized_total_volume"]]) + 0.6))]

    def input_pre_process(self):
        """ Pre-process the user input to prepare lists for reluctance model"""
        all_manual_combinations = list(product(self.manual_core_w, self.manual_window_h, self.manual_window_w))
        manual_core_w = [item[0] for item in all_manual_combinations]
        manual_window_h = [item[1] for item in all_manual_combinations]
        manual_window_w = [item[2] for item in all_manual_combinations]

        core_db = fmt.core_database()
        db_core_w = [core_db[core_name]["core_inner_diameter"] for core_name in self.db_core_names]
        db_window_h = [core_db[core_name]["window_h"] for core_name in self.db_core_names]
        db_window_w = [core_db[core_name]["window_w"] for core_name in self.db_core_names]

        core_w_list = db_core_w + manual_core_w
        window_h_list = db_window_h + manual_window_h
        window_w_list = db_window_w + manual_window_w

        litz_db = fmt.litz_database()
        litz_conductor_r = [litz_db[litz_name]["conductor_radii"] for litz_name in self.db_litz_names]
        litz_strand_r = [litz_db[litz_name]["strand_radii"] for litz_name in self.db_litz_names]
        litz_strand_n = [litz_db[litz_name]["strands_numbers"] for litz_name in self.db_litz_names]

        mu_rel = [material_db.get_material_property(material_name=material_name, property="initial_permeability")
                  for material_name in self.material_list]

        mult_air_gap_type_list = []
        for item in self.mult_air_gap_type:
            if item == 'edge_distributed':
                mult_air_gap_type_list.append(1)
            elif item == 'center_distributed':
                mult_air_gap_type_list.append(2)
            else:
                raise Exception('Wrong string input for multiple air-gap type')

        return core_w_list, window_h_list, window_w_list, litz_conductor_r, litz_strand_r, litz_strand_n, mu_rel, mult_air_gap_type_list

    def filter1_geometry(self, data_matrix):
        """
        Filter out design cases which are not physical possible based on no_of_turns and winding area

        param data_matrix: Matrix containing the design parameters
        :type data_matrix: array
        """
        window_area = data_matrix[:, self.param["window_h"]] * data_matrix[:, self.param["window_w"]]
        insulation_area = ((self.left_core_insulation + self.right_core_insulation) *
                           data_matrix[:, self.param["window_h"]]) + \
                          ((self.top_core_insulation + self.bot_core_insulation) *
                           (data_matrix[:, self.param["window_w"]] -
                            (self.left_core_insulation + self.right_core_insulation)))

        data_matrix = data_matrix[
            np.where((data_matrix[:, self.param["no_of_turns"]] * np.pi * self.min_conductor_r ** 2)
                     < (self.winding_factor * (window_area - insulation_area)))]

        return data_matrix

    def filter2_inductance(self, data_matrix):
        """
        Filter out design cases which are in between the given goal inductance tolerance band

        param data_matrix: Matrix containing the design parameters
        :type data_matrix: array
        """
        data_matrix = data_matrix[np.where((data_matrix[:, self.param["inductance"]] >
                                            ((100 - self.inductance_percent_tolerance) / 100) *
                                            self.goal_inductance) &
                                           (data_matrix[:, self.param["inductance"]] <
                                            ((100 + self.inductance_percent_tolerance) / 100) * self.goal_inductance))]
        return data_matrix

    def filter3_flux_saturation(self, data_matrix):
        """
        Filter out design cases based on the maximum magnetic flux allowed in the magnetic core.

        param data_matrix: Matrix containing the design parameters
        :type data_matrix: array
        """
        b_sat_dict = {}
        for material_name in self.material_list:
            b_sat_key = material_db.get_material_property(material_name=material_name, property="initial_permeability")
            b_sat_dict[b_sat_key] = material_db.get_material_property(material_name=material_name,
                                                                      property="max_flux_density")

        # Creating B_saturated array corresponding to the material type
        b_sat = np.zeros((len(data_matrix), 1))
        for index in range(len(data_matrix)):
            b_sat[index] = b_sat_dict[data_matrix[index, self.param["mu_rel"]]]

        # flux_max = L * i_max / N
        total_flux_max = (data_matrix[:, self.param["inductance"]] * self.current_max) / data_matrix[:,
                                                                                         self.param["no_of_turns"]]
        b_max_center = total_flux_max / data_matrix[:, self.param["center_leg_area"]]
        b_max_middle = total_flux_max / (
                np.pi * data_matrix[:, self.param["core_w"]] * data_matrix[:, self.param["core_h_middle"]])
        b_max_outer = total_flux_max / data_matrix[:, self.param["outer_leg_area"]]

        data_matrix = self.add_column_to_data_matrix(data_matrix, total_flux_max, 'total_flux_max')
        data_matrix = self.add_column_to_data_matrix(data_matrix, b_max_center, 'b_max_center')
        data_matrix = self.add_column_to_data_matrix(data_matrix, b_max_middle, 'b_max_middle')
        data_matrix = self.add_column_to_data_matrix(data_matrix, b_max_outer, 'b_max_outer')

        data_matrix_temp = np.zeros((0, len(data_matrix[0])))
        for index in range(len(data_matrix)):
            if (data_matrix[index, self.param["b_max_center"]] < (self.percent_of_b_sat / 100) * b_sat[index]) & \
                    (data_matrix[index, self.param["b_max_outer"]] < (self.percent_of_b_sat / 100) * b_sat[index]):
                data_matrix_temp = np.vstack([data_matrix_temp, data_matrix[index, :]])

        return data_matrix_temp

    def filter4_losses(self, data_matrix):
        """
       Filter out design cases based on the calculated hysteresis and DC loss

       param data_matrix: Matrix containing the design parameters
       :type data_matrix: array
        """
        # Filter out data-matrix according to calculated hysteresis loss + DC winding loss
        # Volume chosen as per "Masterthesis_Till_Piepenbrock" pg-45
        mu_imag_dict = {}
        counter = 0
        for material_name in self.material_list:
            mu_imag_key = material_db.get_material_property(material_name=material_name,
                                                            property="initial_permeability")
            mu_imag_dict[mu_imag_key] = counter
            counter = counter + 1

        material_data_list = [material_db.permeability_data_to_pro_file(self.temperature, self.frequency,
                                                                        material_name, "manufacturer_datasheet")
                              for material_name in self.material_list]

        mu_imag_interpol_func = [interp1d(material_data_list[i][0], material_data_list[i][1], kind="cubic") for i in
                                 range(len(self.material_list))]

        # Creating mu_imag array corresponding to the material type
        mu_imag = np.zeros(len(data_matrix))
        for index in range(len(data_matrix)):
            mu_imag[index] = mu_imag_interpol_func[mu_imag_dict[data_matrix[index, self.param["mu_rel"]]]] \
                (data_matrix[index, self.param["b_max_center"]])

        volume_center = (np.pi * (data_matrix[:, self.param["core_w"]] / 2) ** 2) * \
                        (data_matrix[:, self.param["window_h"]] + data_matrix[:, self.param["core_h_middle"]] -
                         (data_matrix[:, self.param["n_air_gaps"]] * data_matrix[:, self.param["air_gap_h"]]))
        volume_outer = (np.pi * ((data_matrix[:, self.param["r_outer"]] ** 2) -
                                 (data_matrix[:, self.param["r_inner"]] ** 2))) * \
                       (data_matrix[:, self.param["window_h"]] + data_matrix[:, self.param["core_h_middle"]])

        P_hyst_center = 0.5 * (2 * np.pi * self.frequency) * fmt.mu0 * mu_imag * (
                    (data_matrix[:, self.param["b_max_center"]] /
                     (fmt.mu0 * data_matrix[:,
                                self.param["mu_rel"]])) ** 2)

        P_hyst_outer = 0.5 * (2 * np.pi * self.frequency) * mu_imag * fmt.mu0 * (
                    (data_matrix[:, self.param["b_max_outer"]] /
                     (fmt.mu0 * data_matrix[:, self.param["mu_rel"]])) ** 2)

        P_hyst_density_center = P_hyst_center * volume_center
        P_hyst_density_middle = 0.5 * (2 * np.pi * self.frequency) * mu_imag * fmt.mu0 * \
                                ((data_matrix[:, self.param["total_flux_max"]] / (
                                        fmt.mu0 * data_matrix[:, self.param["mu_rel"]])) ** 2) * \
                                (1 / (2 * np.pi * data_matrix[:, self.param["core_h_middle"]])) * \
                                np.log(
                                    (data_matrix[:, self.param["r_inner"]] * 2) / data_matrix[:, self.param["core_w"]])
        P_hyst_density_outer = P_hyst_outer * volume_outer
        total_hyst_loss = P_hyst_density_center + (2 * P_hyst_density_middle) + P_hyst_density_outer

        # Winding loss (only DC loss)
        Resistance = (data_matrix[:, self.param["no_of_turns"]] * 2 * np.pi *
                      (data_matrix[:, self.param["core_w"]] / 2 + self.min_conductor_r)) / \
                     (self.copper_conductivity * (np.pi * (self.min_conductor_r ** 2)))

        DC_loss = ((self.current_max ** 2) / 2) * Resistance

        total_loss = DC_loss + total_hyst_loss
        max_total_loss = max(total_loss)
        normalized_total_loss = total_loss / max_total_loss

        total_volume = np.pi * (data_matrix[:, self.param["r_outer"]] ** 2) * data_matrix[:, self.param["core_h"]]
        max_volume = max(total_volume)
        normalized_total_volume = total_volume / max_volume

        data_matrix = self.add_column_to_data_matrix(data_matrix, total_hyst_loss, 'total_hyst_loss')
        data_matrix = self.add_column_to_data_matrix(data_matrix, DC_loss, 'DC_loss')
        data_matrix = self.add_column_to_data_matrix(data_matrix, total_loss, 'total_loss')
        data_matrix = self.add_column_to_data_matrix(data_matrix, normalized_total_loss, 'normalized_total_loss')
        data_matrix = self.add_column_to_data_matrix(data_matrix, total_volume, 'total_volume')
        data_matrix = self.add_column_to_data_matrix(data_matrix, normalized_total_volume, 'normalized_total_volume')

        data_matrix = data_matrix[data_matrix[:, self.param["total_loss"]].argsort()]
        data_matrix = data_matrix[0:int((self.percent_of_total_loss / 100) * len(data_matrix)), :]

        return data_matrix

    def fem_simulation(self):
        """
        FEM simulation of the design cases and saving the result in the given working directory for later analysis
        """
        example_results_folder = os.path.join(self.working_directory, "example_results")
        if not os.path.exists(example_results_folder):
            os.mkdir(example_results_folder)

        working_directory = os.path.join(example_results_folder, "inductor")
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        data_folder = os.path.join(self.working_directory, 'fem_simulation_data')
        if not os.path.exists(data_folder):
            os.mkdir(data_folder)

        src_path = os.path.join(self.working_directory, "example_results/inductor/results/log_electro_magnetic.json")

        data_files = []
        file_names = []
        counter3 = 0
        for j in range(len(self.conductor_r_list)):
            for i in range(len(self.data_matrix_fem)):
                window_area = self.data_matrix_fem[i, self.param["window_h"]] * self.data_matrix_fem[
                    i, self.param["window_w"]]
                insulation_area = ((self.left_core_insulation + self.right_core_insulation) *
                                   self.data_matrix_fem[i, self.param["window_h"]]) + \
                                  ((self.top_core_insulation + self.bot_core_insulation) *
                                   (self.data_matrix_fem[i, self.param["window_w"]] -
                                    (self.left_core_insulation + self.right_core_insulation)))

                if not ((self.data_matrix_fem[i, self.param["no_of_turns"]] * np.pi * self.conductor_r_list[j] ** 2)
                        < (self.winding_factor * (window_area - insulation_area))):
                    continue

                # MagneticComponent class object
                geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor,
                                            working_directory=working_directory, silent=True)

                core = fmt.Core(core_inner_diameter=self.data_matrix_fem[i, self.param["core_w"]],
                                window_w=self.data_matrix_fem[i, self.param["window_w"]],
                                window_h=self.data_matrix_fem[i, self.param["window_h"]],
                                # material="N95", temperature=25, frequency=freq, datasource="manufacturer_datasheet")
                                # material="95_100")
                                mu_rel=self.data_matrix_fem[i, self.param["mu_rel"]], phi_mu_deg=10,
                                sigma=0.5)
                # TODO: (completed)
                # mu_rel=3000, phi_mu_deg=10,
                # sigma=0.5)
                geo.set_core(core)

                # 3. set air gap parameters
                air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
                if int(self.data_matrix_fem[i, self.param["n_air_gaps"]]) == 1:
                    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg,
                                         self.data_matrix_fem[i, self.param["air_gap_h"]],
                                         self.data_matrix_fem[i, self.param["air_gap_position"]])
                else:
                    if int(self.data_matrix_fem[i, self.param["mult_air_gap_type"]]) == 1:
                        position_list = list(
                            np.linspace(0, 100, int(self.data_matrix_fem[i, self.param["n_air_gaps"]])))
                        for position in position_list:
                            air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg,
                                                 self.data_matrix_fem[i, self.param["air_gap_h"]],
                                                 position)

                    elif int(self.data_matrix_fem[i, self.param["mult_air_gap_type"]]) == 2:
                        position_list = list(
                            np.linspace(0, 100, int(self.data_matrix_fem[i, self.param["n_air_gaps"]]) + 2))
                        position_list.remove(0.0)
                        position_list.remove(100.0)
                        for position in position_list:
                            air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg,
                                                 self.data_matrix_fem[i, self.param["air_gap_h"]],
                                                 position)
                geo.set_air_gaps(air_gaps)

                # 4. set insulations
                insulation = fmt.Insulation()
                insulation.add_core_insulations(self.top_core_insulation, self.bot_core_insulation,
                                                self.left_core_insulation, self.right_core_insulation)
                insulation.add_winding_insulations(self.inner_winding_insulations, 0.0001)
                geo.set_insulation(insulation)

                # 5. create winding window and virtual winding windows (vww)
                winding_window = fmt.WindingWindow(core, insulation)
                vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

                # 6. create conductor and set parameters: use solid wires
                winding = fmt.Conductor(0, fmt.Conductivity.Copper)
                if j < len(self.litz_conductor_r):
                    winding.set_litz_round_conductor(conductor_radius=self.litz_conductor_r[j],
                                                     number_strands=self.litz_strand_n[j],
                                                     strand_radius=self.litz_strand_r[j], fill_factor=None,
                                                     conductor_arrangement=self.winding_scheme_dict[self.winding_scheme])
                else:
                    winding.set_solid_round_conductor(conductor_radius=self.conductor_r_list[j],
                                                      conductor_arrangement=self.winding_scheme_dict[self.winding_scheme])
                # winding.set_litz_round_conductor(conductor_radius=0.0013, number_strands=150, strand_radius=100e-6,
                # fill_factor=None, conductor_arrangement=fmt.ConductorArrangement.Square)

                # 7. add conductor to vww and add winding window to MagneticComponent
                vww.set_winding(winding, int(self.data_matrix_fem[i, self.param["no_of_turns"]]), None)
                geo.set_winding_window(winding_window)

                try:
                    # 5. create the model
                    geo.create_model(freq=self.frequency, visualize_before=False, save_png=False)

                    # 6. start simulation
                    geo.single_simulation(freq=self.frequency, current=[self.current_max], show_results=False)

                    shutil.copy2(src_path, data_folder)
                    old_filename = os.path.join(data_folder, "log_electro_magnetic.json")
                    new_filename = os.path.join(data_folder, f"case{counter3}.json")
                    os.rename(old_filename, new_filename)
                    data_files.append(new_filename)
                    file_names.append(f"case{counter3}")
                    # print(f"{counter3} of {n_cases_FEM * 2}")
                    counter3 = counter3 + 1

                except (Exception,) as e:
                    print("next iteration")
                    logging.exception(e)

    def add_column_to_data_matrix(self, data_matrix, column_value: list, column_name: str):
        """
        Adds column to the given matrix

        param data_matrix: Matrix containing the design parameters
       :type data_matrix: array
       :param column_value: Column to be added
       :type column_value: list
       :param column_name: Identifier of the column
       :type column_name: str
        """
        size = len(data_matrix[0])
        data_matrix = np.hstack((data_matrix, np.reshape(column_value, (len(column_value), 1))))
        self.param[column_name] = size

        return data_matrix

    def plot_volume_loss(self, data_matrix):
        """
        Plots estimated normalised volume vs loss graph from reluctance model results

        param data_matrix: Matrix containing the design parameters
       :type data_matrix: array
        """
        fig, ax = fmt.plt.subplots()  # Create a figure containing a single axes.
        fmt.plt.title("Normalised volume vs Normalised losses")
        fmt.plt.xlabel("Normalised Volume")
        fmt.plt.ylabel("Normalised Losses")
        ax.plot(data_matrix[:, self.param['normalized_total_volume']],
                data_matrix[:, self.param['normalized_total_loss']], 'o')
        ax.grid()
        fmt.plt.show()


if __name__ == '__main__':
    ad = AutomatedDesign(working_directory='D:/Personal_data/MS_Paderborn/Sem4/Project_2/sweep3d',
                         component='inductor',
                         goal_inductance=120 * 1e-6,
                         frequency=100000,
                         inductance_percent_tolerance=10,
                         winding_scheme='Square',
                         current_max=8,
                         percent_of_b_sat=70,
                         percent_of_total_loss=1,
                         database_core_names=[],
                         database_litz_names=['1.5x105x0.1'],
                         solid_conductor_r=[0.0013],
                         manual_core_w=list(np.linspace(0.005, 0.05, 10)),
                         manual_window_h=list(np.linspace(0.01, 0.08, 5)),
                         manual_window_w=list(np.linspace(0.005, 0.04, 10)),
                         no_of_turns=[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
                         n_air_gaps=[1, 2],
                         air_gap_height=list(np.linspace(0.0001, 0.0005, 5)),
                         air_gap_position=list(np.linspace(20, 80, 2)),
                         material_list=['N95'],
                         mult_air_gap_type=['center_distributed'],
                         top_core_insulation=0.0001,
                         bot_core_insulation=0.0001,
                         left_core_insulation=0.0004,
                         right_core_insulation=0.0001,
                         inner_winding_insulations=[0.0005],
                         temperature=25.0)

    ad.fem_simulation()
    real_inductance, total_loss, total_volume, total_cost, labels = load_design \
        (working_directory='D:/Personal_data/MS_Paderborn/Sem4/Project_2/sweep3d')
    # plot_2d(x_value=total_volume, y_value=total_loss, x_label='Volume / m\u00b3', y_label='Loss / W',
    #         title='Volume vs Loss', annotations=labels)
    # plot_2d(x_value=total_volume, y_value=total_cost, x_label='Volume / m\u00b3', y_label='Cost / \u20ac',
    #         title='Volume vs Cost', annotations=labels)
    plot_3d(x_value=total_volume, y_value=total_loss, z_value=total_cost, x_label='Volume / m\u00b3',
            y_label='Loss / W', z_label='Cost / \u20ac', x_limit=[0, 0.0005], y_limit=[0, 20],
            z_limit=[4.5, 5.5], title='Volume vs Loss vs Cost', annotations=labels)

# def automated_design_func():
#     # ########################################   {DESIGN PARAMETERS}   #################################################
#     save_directory_name = "sweep_new_3d_with_cost"  # New directory is created in FEM_Magnetics_Toolbox/femmt/examples/
#     goal_inductance = 120 * 1e-6
#     L_tolerance_percent = 10                # Inductance tolerance of +-10% is applied
#     winding_factor = 0.91
#     i_max = 8                               # Max current amplitude with assumption of sinusoidal current waveform 8
#     percent_of_B_sat = 70                   # Percent of B_sat allowed in the designed core
#     percent_of_total_loss = 30              # Percent of total_loss allowed in FEM simulation 30
#
#     freq = 100 * 1e3                        # Switching frequency in Hz
#     # mu_imag = 100                           # TODO: (completed)
#     Cu_sigma = 5.96 * 1e7                   # copper conductivity (sigma) @ 20 degree celsius
#
#
#     # Set core-geometry from core database or/and manual entry
#     db_core_names = []  # "PQ 40/40", "PQ 40/30"    blank
#
#     manual_core_w = list(np.linspace(0.005, 0.05, 10))# 10
#     manual_window_h = list(np.linspace(0.01, 0.08, 5))# 5
#     manual_window_w = list(np.linspace(0.005, 0.04, 10))# 10
#
#     all_manual_combinations = list(product(manual_core_w, manual_window_h, manual_window_w))
#     manual_core_w = [item[0] for item in all_manual_combinations]
#     manual_window_h = [item[1] for item in all_manual_combinations]
#     manual_window_w = [item[2] for item in all_manual_combinations]
#
#     core_db = fmt.core_database()
#     db_core_w = [core_db[core_name]["core_inner_diameter"] for core_name in db_core_names]
#     db_window_h = [core_db[core_name]["window_h"] for core_name in db_core_names]
#     db_window_w = [core_db[core_name]["window_w"] for core_name in db_core_names]
#
#     core_w_list = db_core_w + manual_core_w
#     window_h_list = db_window_h + manual_window_h
#     window_w_list = db_window_w + manual_window_w
#
#     # Set winding settings (Solid and Litz winding type)
#     solid_conductor_r = [0.0013]  # blank
#     litz_names = ["1.5x105x0.1"]  # "1.5x105x0.1", "1.4x200x0.071"
#
#     litz_db = fmt.litz_database()
#     litz_conductor_r = [litz_db[litz_name]["conductor_radii"] for litz_name in litz_names]
#     litz_strand_r = [litz_db[litz_name]["strand_radii"] for litz_name in litz_names]
#     litz_strand_n = [litz_db[litz_name]["strands_numbers"] for litz_name in litz_names]
#
#     min_conductor_r = min(litz_conductor_r + solid_conductor_r)
#
#     # Set air-gap and core parameters
#     no_of_turns = [2,3,4,5,6,7,8, 9, 10, 11, 12, 13, 14,15,16,17,18,19,20]  # Set No. of turns (N)   2,3,4,5,6,7,8, 9, 10, 11, 12, 13, 14,15,16,17,18,19,20
#     print(f"{no_of_turns = }")
#     n_air_gaps = [1, 2]  # Set No. of air-gaps (n)
#     air_gap_height = list(np.linspace(0.0001, 0.0005, 5))  # Set air-gap length in metre (l)
#     air_gap_position = list(np.linspace(20, 80, 2))  # Set air-gap position in percent w.r.t. core window height
#
#     material_names = ["N95"]  # Set relative permeability in F/m (u) , "N87"
#
#     mu_rel = [material_db.get_material_property(material_name=material_name, property="initial_permeability")
#               for material_name in material_names]
#
#     temp_var = [material_db.permeability_data_to_pro_file(25, freq, material_name, "manufacturer_datasheet")
#                 for material_name in material_names]
#     print(temp_var[0][1])
#     mu_imag_interpol_func = [interp1d(temp_var[i][0], temp_var[i][1], kind="cubic") for i in range(len(material_names))]
#     # print(f[0](0.2))
#     # Set two types of equally distributed air-gaps (used only for air-gaps more than 1):
#     # Type 1: Equally distributed air-gaps including corner air-gaps (eg: air-gaps-position = [0, 50, 100])
#     # Type 2: Equally distributed air-gaps excluding corner air-gaps (eg: air-gaps-position = [25, 50, 75])
#     # 'Type1 = with corner air-gaps; 'Type2' = without corner air-gaps; 'Type0' = single air-gap
#     mult_air_gap_type = [2] #Type1-Edge, Type2: Centre
#     # TODO: check if the issue has been resolved
#
#     # ######################################   {RELUCTANCE_CALCULATION}   ##############################################
#     # Call to Reluctance model (Class MagneticCircuit)
#     mc1 = fmt.MagneticCircuit(core_w=core_w_list, window_h=window_h_list, window_w=window_w_list,
#                               no_of_turns=no_of_turns, n_air_gaps=n_air_gaps, air_gap_h=air_gap_height,
#                               air_gap_position=air_gap_position, mu_rel=mu_rel, mult_air_gap_type=mult_air_gap_type,
#                               air_gap_method='percent', component_type='inductor', sim_type='sweep')
#     param = mc1.get_param_pos_dict()
#     n_cases_0 = len(mc1.data_matrix)
#     print(n_cases_0)
#     # Calculate total reluctance and creates data_matrix to access all corresponding parameters and results
#     # To access all/any data from MagneticCircuit class, use self.data_matrix[:, param["parameter_name"]].
#     # The parameters are arranged as shown below:
#     # Example: If you want to access inductance, type self.data_matrix[:, param["inductance"]]
#
#     # ############################################   {FILTRATION}   ####################################################
#     # 1st Filter: ------------------------------------------------------------------------------------------------------
#     # Filter out cases where physical geometry is not possible
#     data_matrix_1 = mc1.data_matrix[np.where((mc1.data_matrix[:, param["no_of_turns"]] * np.pi * min_conductor_r ** 2)
#                                              < (winding_factor * mc1.data_matrix[:,
#                                                                  param["window_h"]] * mc1.data_matrix[:,
#                                                                                       param[
#                                                                                           "window_w"]]))]
#     n_cases_1 = len(data_matrix_1)
#     print(n_cases_1)
#     # 2nd Filter:-------------------------------------------------------------------------------------------------------
#     # Based on +-10% goal inductance tolerance band
#     data_matrix_2 = data_matrix_1[
#         np.where((data_matrix_1[:, param["inductance"]] > ((100 - L_tolerance_percent) / 100) * goal_inductance) &
#                  (data_matrix_1[:, param["inductance"]] < ((100 + L_tolerance_percent) / 100) * goal_inductance))]
#     n_cases_2 = len(data_matrix_2)
#     print(n_cases_2)
#     # 3rd Filter:-------------------------------------------------------------------------------------------------------
#     # Filter out cases where B_max is greater than B_sat
#     B_sat_dict = {}
#     counter1 = 0
#     for material_name in material_names:
#         B_sat_key = material_db.get_material_property(material_name=material_name, property="initial_permeability")
#         B_sat_dict[B_sat_key] = material_db.get_material_property(material_name=material_name, property="max_flux_density")
#         counter1 = counter1 + 1
#
#     # Creating B_saturated array corresponding to the material type
#     B_sat = np.zeros((len(data_matrix_2), 1))
#     for index in range(len(data_matrix_2)):
#         B_sat[index] = B_sat_dict[data_matrix_2[index, param["mu_rel"]]]
#
#     # flux_max = L * i_max / N
#     total_flux_max = (data_matrix_2[:, param["inductance"]] * i_max) / data_matrix_2[:, param["no_of_turns"]]
#     data_matrix_2 = np.hstack((data_matrix_2, np.reshape(total_flux_max, (len(total_flux_max), 1))))
#     param["total_flux_max"] = 15
#
#     B_max_center = total_flux_max / data_matrix_2[:, param["center_leg_area"]]
#     B_max_middle = total_flux_max / (
#             np.pi * data_matrix_2[:, param["core_w"]] * data_matrix_2[:, param["core_h_middle"]])
#     B_max_outer = total_flux_max / data_matrix_2[:, param["outer_leg_area"]]
#
#     data_matrix_2 = np.hstack((data_matrix_2, np.reshape(B_max_center, (len(B_max_center), 1))))
#     param["B_max_center"] = 16
#     data_matrix_2 = np.hstack((data_matrix_2, np.reshape(B_max_outer, (len(B_max_outer), 1))))
#     param["B_max_outer"] = 17
#
#     data_matrix_3 = np.zeros((0, 18))
#     for index in range(len(data_matrix_2)):
#         if (data_matrix_2[index, param["B_max_center"]] < (percent_of_B_sat / 100) * B_sat[index]) & \
#                 (data_matrix_2[index, param["B_max_outer"]] < (percent_of_B_sat / 100) * B_sat[index]):
#             data_matrix_3 = np.vstack([data_matrix_3, data_matrix_2[index, :]])
#     n_cases_3 = len(data_matrix_3)
#     print(n_cases_3)
#     # 4th Filter:-------------------------------------------------------------------------------------------------------
#     # Filter out data-matrix according to calculated hysteresis loss + DC winding loss
#     # Volume chosen as per "Masterthesis_Till_Piepenbrock" pg-45
#     mu_imag_dict = {}
#     counter_x = 0
#     for material_name in material_names:
#         mu_imag_key = material_db.get_material_property(material_name=material_name, property="initial_permeability")
#         mu_imag_dict[mu_imag_key] = counter_x
#         counter_x = counter_x + 1
#
#     # Creating B_saturated array corresponding to the material type
#     mu_imag = np.zeros(len(data_matrix_3))
#     for index in range(len(data_matrix_3)):
#         mu_imag[index] = mu_imag_interpol_func[mu_imag_dict[data_matrix_3[index, param["mu_rel"]]]](data_matrix_3[index, param["B_max_center"]])
#
#
#     volume_center = (np.pi * (data_matrix_3[:, param["core_w"]] / 2) ** 2) * \
#                     (data_matrix_3[:, param["window_h"]] + data_matrix_3[:, param["core_h_middle"]] -
#                      (data_matrix_3[:, param["n_air_gaps"]] * data_matrix_3[:, param["air_gap_h"]]))
#     volume_outer = (np.pi * ((data_matrix_3[:, param["r_outer"]] ** 2) - (data_matrix_3[:, param["r_inner"]] ** 2))) * \
#                    (data_matrix_3[:, param["window_h"]] + data_matrix_3[:, param["core_h_middle"]])
#
#     P_hyst_center = 0.5 * (2 * np.pi * freq) * fmt.mu0 * mu_imag * ((data_matrix_3[:, param["B_max_center"]] /
#                                                                      (fmt.mu0 * data_matrix_3[:,
#                                                                                 param["mu_rel"]])) ** 2)
#     # P_hyst_center = np.array(mu_imag) * np.array(P_hyst_center)
#     P_hyst_outer = 0.5 * (2 * np.pi * freq) * mu_imag * fmt.mu0 * ((data_matrix_3[:, param["B_max_outer"]] /
#                                                                     (fmt.mu0 * data_matrix_3[:, param["mu_rel"]])) ** 2)
#
#     P_hyst_density_center = P_hyst_center * volume_center
#     P_hyst_density_middle = 0.5 * (2 * np.pi * freq) * mu_imag * fmt.mu0 * \
#                             ((data_matrix_3[:, param["total_flux_max"]] / (
#                                     fmt.mu0 * data_matrix_3[:, param["mu_rel"]])) ** 2) * \
#                             (1 / (2 * np.pi * data_matrix_3[:, param["core_h_middle"]])) * \
#                             np.log((data_matrix_3[:, param["r_inner"]] * 2) / data_matrix_3[:, param["core_w"]])
#     P_hyst_density_outer = P_hyst_outer * volume_outer
#
#     total_hyst_loss = P_hyst_density_center + (2 * P_hyst_density_middle) + P_hyst_density_outer
#     data_matrix_3 = np.hstack((data_matrix_3, np.reshape(total_hyst_loss, (len(total_hyst_loss), 1))))  # position: 18
#     param["P_hyst_density_total"] = 18
#
#     # Winding loss (only DC loss)
#     Resistance = (data_matrix_3[:, param["no_of_turns"]] * 2 * np.pi *
#                   (data_matrix_3[:, param["core_w"]] / 2 + min_conductor_r)) / \
#                  (Cu_sigma * (np.pi * (min_conductor_r ** 2)))
#
#     DC_loss = ((i_max ** 2) / 2) * Resistance
#     data_matrix_3 = np.hstack((data_matrix_3, np.reshape(DC_loss, (len(DC_loss), 1))))  # position: 19
#     param["DC_loss"] = 19
#
#     total_loss = DC_loss + total_hyst_loss
#     data_matrix_3 = np.hstack((data_matrix_3, np.reshape(total_loss, (len(total_loss), 1))))  # position: 20
#     param["total_loss"] = 20
#     data_matrix_3 = data_matrix_3[data_matrix_3[:, param["total_loss"]].argsort()]
#     data_matrix_4 = data_matrix_3[0:int((percent_of_total_loss / 100) * len(data_matrix_3)), :]
#
#     total_loss = data_matrix_4[:, param["total_loss"]]
#     max_total_loss = max(total_loss)
#     normalized_total_loss = total_loss / max_total_loss
#     data_matrix_4 = np.hstack((data_matrix_4, np.reshape(normalized_total_loss, (len(normalized_total_loss), 1))))  # position: 20
#     param["normalized_total_loss"] = 21
#
#     total_volume = np.pi * (data_matrix_4[:, param["r_outer"]] ** 2) * data_matrix_4[:, param["core_h_middle"]]
#     max_volume = max(total_volume)
#     normalized_total_volume = total_volume / max_volume
#     data_matrix_4 = np.hstack((data_matrix_4, np.reshape(total_volume, (len(total_volume), 1))))  # position: 20
#     param["total_volume"] = 22
#     data_matrix_4 = np.hstack((data_matrix_4, np.reshape(normalized_total_volume, (len(normalized_total_volume), 1))))  # position: 20
#     param["normalized_total_volume"] = 23
#
#     # Sort the data_matrix with respect to total losses column----------------------------------------------------------
#     # FEM_data_matrix = data_matrix_4[np.where(data_matrix_4[:, param["normalized_total_loss"]] <= ((0.01 / data_matrix_4[:, param["normalized_total_volume"]]) + 0.6))]
#     FEM_data_matrix = data_matrix_4
#     n_cases_FEM = len(FEM_data_matrix)
#     print(n_cases_FEM)
#
#     fig, ax = fmt.plt.subplots()  # Create a figure containing a single axes.
#     fmt.plt.title("Normalised volume vs Normalised losses")
#     fmt.plt.xlabel("Normalised Volume")
#     fmt.plt.ylabel("Normalised Losses")
#     ax.plot(FEM_data_matrix[:, param['normalized_total_volume']], FEM_data_matrix[:, param['normalized_total_loss']], 'o')
#     ax.grid()
#     fmt.plt.show()
#
#     # ##########################################   {FEM_SIMULATION}   ##################################################
#     qwerty = 1
#
#     example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
#     if not os.path.exists(example_results_folder):
#         os.mkdir(example_results_folder)
#
#     working_directory = os.path.join(example_results_folder, "inductor")
#     if not os.path.exists(working_directory):
#         os.mkdir(working_directory)
#
#     if not os.path.exists(os.path.join(os.path.dirname(__file__), save_directory_name)):
#         os.mkdir(os.path.join(os.path.dirname(__file__), save_directory_name))
#
#     working_directories = []
#     file_names = []
#
#     src_path = os.path.join(os.path.dirname(__file__), "example_results/inductor/results/log_electro_magnetic.json")
#
#     counter3 = 0
#     for j in range(len(solid_conductor_r) + len(litz_conductor_r)):
#         conductor_r_list = litz_conductor_r + solid_conductor_r
#         for i in range(len(FEM_data_matrix)):
#             print(f"value of i:{i}")
#             print(f"value of j:{j}")
#             if not ((FEM_data_matrix[i, param["no_of_turns"]] * np.pi * conductor_r_list[j] ** 2)
#                     < (winding_factor * FEM_data_matrix[i, param["window_h"]] * mc1.data_matrix[i, param["window_w"]])):
#                 continue
#
#             # MagneticComponent class object
#             geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=working_directory, silent=True)
#
#             core = fmt.Core(core_inner_diameter=FEM_data_matrix[i, param["core_w"]], window_w=FEM_data_matrix[i, param["window_w"]],
#                             window_h=FEM_data_matrix[i, param["window_h"]],
#                             # material="N95", temperature=25, frequency=freq, datasource="manufacturer_datasheet")
#                             # material="95_100")
#             mu_rel=3000, phi_mu_deg=10,
#             sigma=0.5)
#             # TODO: (completed)
#             # mu_rel=3000, phi_mu_deg=10,
#             # sigma=0.5)
#             geo.set_core(core)
#
#             # 3. set air gap parameters
#             air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
#             if int(FEM_data_matrix[i, param["n_air_gaps"]]) == 1:
#                 air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, FEM_data_matrix[i, param["air_gap_h"]],
#                                      FEM_data_matrix[i, param["air_gap_position"]])
#             else:
#                 if int(FEM_data_matrix[i, param["mult_air_gap_type"]]) == 1:
#                     position_list = list(np.linspace(0, 100, int(FEM_data_matrix[i, param["n_air_gaps"]])))
#                     for position in position_list:
#                         air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, FEM_data_matrix[i, param["air_gap_h"]],
#                                              position)
#
#                 elif int(FEM_data_matrix[i, param["mult_air_gap_type"]]) == 2:
#                     position_list = list(np.linspace(0, 100, int(FEM_data_matrix[i, param["n_air_gaps"]]) + 2))
#                     position_list.remove(0.0)
#                     position_list.remove(100.0)
#                     for position in position_list:
#                         air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, FEM_data_matrix[i, param["air_gap_h"]],
#                                              position)
#             geo.set_air_gaps(air_gaps)
#
#             # 4. set insulations
#             insulation = fmt.Insulation()
#             insulation.add_core_insulations(0.001, 0.001, 0.004, 0.001)
#             insulation.add_winding_insulations([0.0005], 0.0001)
#             geo.set_insulation(insulation)
#
#             # 5. create winding window and virtual winding windows (vww)
#             winding_window = fmt.WindingWindow(core, insulation)
#             vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)
#
#             # 6. create conductor and set parameters: use solid wires
#             winding = fmt.Conductor(0, fmt.Conductivity.Copper)
#             if j < len(litz_conductor_r):
#                 winding.set_litz_round_conductor(conductor_radius=litz_conductor_r[j], number_strands=litz_strand_n[j],
#                                                  strand_radius=litz_strand_r[j], fill_factor=None,
#                                                  conductor_arrangement=fmt.ConductorArrangement.Square)
#             else:
#                 winding.set_solid_round_conductor(conductor_radius=conductor_r_list[j],
#                                                   conductor_arrangement=fmt.ConductorArrangement.Square)
#             # winding.set_litz_round_conductor(conductor_radius=0.0013, number_strands=150, strand_radius=100e-6,
#             # fill_factor=None, conductor_arrangement=fmt.ConductorArrangement.Square)
#
#             # 7. add conductor to vww and add winding window to MagneticComponent
#             vww.set_winding(winding, int(FEM_data_matrix[i, param["no_of_turns"]]), None)
#             geo.set_winding_window(winding_window)
#
#             try:
#                 # 5. create the model
#                 geo.create_model(freq=freq, visualize_before=False, save_png=False)
#
#                 # 6. start simulation
#                 geo.single_simulation(freq=freq, current=[i_max], show_results=False)
#
#                 shutil.copy2(src_path, os.path.join(os.path.dirname(__file__), save_directory_name))
#                 old_filename = os.path.join(os.path.dirname(__file__), save_directory_name, "log_electro_magnetic.json")
#                 new_filename = os.path.join(os.path.dirname(__file__), save_directory_name, f"case{counter3}.json")
#                 os.rename(old_filename, new_filename)
#                 working_directories.append(new_filename)
#                 file_names.append(f"case{counter3}")
#                 counter3 = counter3 + 1
#                 print(f"{counter3} of {n_cases_FEM * 2}")
#
#             except (Exception,) as e:
#                 print("next iteration")
#                 logging.exception(e)
#
#             # # 5. create the model
#             # geo.create_model(freq=freq, visualize_before=False, save_png=False)
#             #
#             # # 6. start simulation
#             # geo.single_simulation(freq=freq, current=[i_max], show_results=False)
#
#
# def load_design(load_directory_name):
#     working_directories = []
#     labels = []
#     working_directory = os.path.join(os.path.dirname(__file__), load_directory_name)
#     print("##########################")
#     print(f"{working_directory =}")
#     print("##########################")
#     file_names = [f for f in listdir(working_directory) if isfile(join(working_directory, f))]
#     file_names.sort()
#     counter2 = 0
#     for name in file_names:
#         temp_var = os.path.join(os.path.dirname(__file__), load_directory_name, name)
#         working_directories.append(temp_var)
#         # labels.append(f"case{counter2}")
#         labels.append(name)
#         counter2 = counter2 + 1
#
#     zip_iterator = zip(file_names, working_directories)
#     logs = dict(zip_iterator)
#
#     # After the simulations the sweep can be analyzed
#     # This could be done using the FEMMTLogParser:
#     log_parser = fmt.FEMMTLogParser(logs)
#
#     # In this case the self inductivity of winding1 will be analyzed
#     inductivities = []
#     active_power = []
#     total_volume = []
#     total_cost = []
#     for name, data in log_parser.data.items():
#         inductivities.append(data.sweeps[0].windings[0].self_inductance)
#         active_power.append(data.sweeps[0].windings[0].active_power)
#         total_volume.append(data.core_2daxi_total_volume)
#         total_cost.append(data.total_cost)
#
#     real_inductance = []
#     for i in range(len(active_power)):
#         real_inductance.append(inductivities[i].real)
#     print(file_names)
#     print(total_volume)
#     print(active_power)
#     print(total_cost)
#
#     names = np.array(labels)
#     c = np.random.randint(1, 5, size=len(active_power))
#
#     norm = plt.Normalize(1, 4)
#     cmap = plt.cm.RdYlGn
#
#     fig, ax = plt.subplots()
#     fmt.plt.title("Volume vs loss")
#     fmt.plt.xlabel("Volume (in cubic m)")
#     fmt.plt.ylabel("Total loss (in W)")
#     sc = plt.scatter(total_volume, active_power, c=c, s=50, cmap=cmap, norm=norm)
#
#     annot = ax.annotate("", xy=(0, 0), xytext=(20, 20), textcoords="offset points",
#                         bbox=dict(boxstyle="round", fc="w"),
#                         arrowprops=dict(arrowstyle="->"))
#     annot.set_visible(False)
#
#     def update_annot(ind):
#         pos = sc.get_offsets()[ind["ind"][0]]
#         annot.xy = pos
#         # text = "{}, {}".format(" ".join(list(map(str, ind["ind"]))),
#         #                       " ".join([names[n] for n in ind["ind"]]))
#         text = "{}".format(" ".join([names[n] for n in ind["ind"]]))
#         annot.set_text(text)
#         annot.get_bbox_patch().set_facecolor(cmap(norm(c[ind["ind"][0]])))
#         annot.get_bbox_patch().set_alpha(0.4)
#
#     def hover(event):
#         vis = annot.get_visible()
#         if event.inaxes == ax:
#             cont, ind = sc.contains(event)
#             if cont:
#                 update_annot(ind)
#                 annot.set_visible(True)
#                 fig.canvas.draw_idle()
#             else:
#                 if vis:
#                     annot.set_visible(False)
#                     fig.canvas.draw_idle()
#
#     fig.canvas.mpl_connect("motion_notify_event", hover)
#     ax.grid()
#     plt.show()
#
#
# def compare_graph(directory_names):
#     temp = 0
#
#     for directory_name in directory_names:
#         working_directories = []
#         labels = []
#         working_directory = os.path.join(os.path.dirname(__file__), directory_name)
#         print("##########################")
#         print(f"{working_directory =}")
#         print("##########################")
#         file_names = [f for f in listdir(working_directory) if isfile(join(working_directory, f))]
#         file_names.sort()
#         counter2 = 0
#         for name in file_names:
#             temp_var = os.path.join(os.path.dirname(__file__), directory_name, name)
#             working_directories.append(temp_var)
#             # labels.append(f"case{counter2}")
#             labels.append(name)
#             counter2 = counter2 + 1
#
#         zip_iterator = zip(file_names, working_directories)
#         logs = dict(zip_iterator)
#
#         # After the simulations the sweep can be analyzed
#         # This could be done using the FEMMTLogParser:
#         log_parser = fmt.FEMMTLogParser(logs)
#
#         # In this case the self inductivity of winding1 will be analyzed
#         inductivities = []
#         active_power = []
#         total_volume = []
#         for name, data in log_parser.data.items():
#             inductivities.append(data.sweeps[0].windings[0].self_inductance)
#             active_power.append(data.sweeps[0].windings[0].active_power)
#             total_volume.append(data.core_2daxi_total_volume)
#
#         real_inductance = []
#         for i in range(len(active_power)):
#             real_inductance.append(inductivities[i].real)
#
#         print(real_inductance)
#         print(active_power)
#
#         # fig, ax = fmt.plt.subplots()  # Create a figure containing a single axes.
#         # fmt.plt.title("Normalised volume vs Normalised losses")
#         # fmt.plt.xlabel("Normalised Volume")
#         # fmt.plt.ylabel("Normalised Losses")
#         # ax.grid()
#         if temp == 0:
#             plt.plot(total_volume, active_power, 'o')
#             temp =1
#         else:
#             plt.plot(total_volume, active_power, '*')
#     plt.title("volume vs losses")
#     plt.xlabel("Total Volume")
#     plt.ylabel("Total Losses")
#     plt.grid()
#     fmt.plt.show()
#
#
# def plot_3d(load_directory_name):
#     working_directories = []
#     labels = []
#     working_directory = os.path.join(os.path.dirname(__file__), load_directory_name)
#     print("##########################")
#     print(f"{working_directory =}")
#     print("##########################")
#     file_names = [f for f in listdir(working_directory) if isfile(join(working_directory, f))]
#     print(file_names)
#     num = [re.findall(r'\d+', item) for item in file_names]
#     case_num_list = [int(item[0]) for item in num]
#     file_names.sort()
#     print(file_names)
#     counter2 = 0
#     for name in file_names:
#         temp_var = os.path.join(os.path.dirname(__file__), load_directory_name, name)
#         working_directories.append(temp_var)
#         # labels.append(f"case{counter2}")
#         labels.append(name)
#         counter2 = counter2 + 1
#
#     zip_iterator = zip(file_names, working_directories)
#     logs = dict(zip_iterator)
#
#     # After the simulations the sweep can be analyzed
#     # This could be done using the FEMMTLogParser:
#     log_parser = fmt.FEMMTLogParser(logs)
#
#     # In this case the self inductivity of winding1 will be analyzed
#     inductivities = []
#     active_power = []
#     total_volume = []
#     total_cost = []
#     for name, data in log_parser.data.items():
#         inductivities.append(data.sweeps[0].windings[0].self_inductance)
#         active_power.append(data.sweeps[0].windings[0].active_power)
#         total_volume.append(data.core_2daxi_total_volume)
#         total_cost.append(data.total_cost)
#
#     real_inductance = []
#     for i in range(len(active_power)):
#         real_inductance.append(inductivities[i].real)
#     print(file_names)
#     print(total_volume)
#     print(active_power)
#     print(total_cost)
#
#     names = np.array(labels)
#     case_num_list = np.array(case_num_list)
#
#     X = np.zeros((len(total_volume), 1))
#     X = np.hstack((X, np.reshape(total_volume, (len(total_volume), 1))))
#     X = np.hstack((X, np.reshape(active_power, (len(active_power), 1))))
#     X = np.hstack((X, np.reshape(total_cost, (len(total_cost), 1))))
#     X = np.hstack((X, np.reshape(case_num_list, (len(case_num_list), 1))))
#     X = np.delete(X, 0, 1)
#     X = X[X[:, 3].argsort()]
#     print(len(X))
#     return X, names
#
#
# def visualize3DData (K):
#     names = K[1]
#     X = K[0]
#     """Visualize data in 3d plot with popover next to mouse position.
#
#     Args:
#         X (np.array) - array of points, of shape (numPoints, 3)
#     Returns:
#         None
#     """
#     fig = plt.figure(figsize = (16,10))
#     ax = fig.add_subplot(111, projection = '3d')
#     ax.set_title('Volume vs loss vs cost')
#     ax.set_xlabel('volume (in cubic m)')
#     ax.set_ylabel('losses (in W)')
#     ax.set_zlabel('cost (in Euro)')
#     ax.set_xlim(0, 0.0005)
#     ax.set_ylim(0, 20)
#     ax.set_zlim(4.5, 5.5)
#     XX = X[:, 0]
#     XX[XX > 0.0005] = np.nan
#     Y = X[:, 1]
#     Y[Y > 20] = np.nan
#     Z = X[:, 2]
#     Z[Z > 5.5] = np.nan
#     c = np.random.randint(1, 5, size=len(X))
#
#     norm = plt.Normalize(1, 4)
#     cmap = plt.cm.RdYlGn
#
#
#     ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=c, s=10, cmap=cmap, norm=norm, depthshade = False, picker = True)
#     # ax.scatter(XX, Y, Z, c=c, s=10, cmap=cmap, norm=norm, depthshade=False, picker=True)
#
#     def distance(point, event):
#         """Return distance between mouse position and given data point
#
#         Args:
#             point (np.array): np.array of shape (3,), with x,y,z in data coords
#             event (MouseEvent): mouse event (which contains mouse position in .x and .xdata)
#         Returns:
#             distance (np.float64): distance (in screen coords) between mouse pos and data point
#         """
#         assert point.shape == (3,), "distance: point.shape is wrong: %s, must be (3,)" % point.shape
#
#         # Project 3d data space to 2d data space
#         x2, y2, _ = proj3d.proj_transform(point[0], point[1], point[2], plt.gca().get_proj())
#         # Convert 2d data space to 2d screen space
#         x3, y3 = ax.transData.transform((x2, y2))
#
#         return np.sqrt ((x3 - event.x)**2 + (y3 - event.y)**2)
#
#     def calcClosestDatapoint(X, event):
#         """"Calculate which data point is closest to the mouse position.
#
#         Args:
#             X (np.array) - array of points, of shape (numPoints, 3)
#             event (MouseEvent) - mouse event (containing mouse position)
#         Returns:
#             smallestIndex (int) - the index (into the array of points X) of the element closest to the mouse position
#         """
#         distances = [distance (X[i, 0:3], event) for i in range(X.shape[0])]
#         return np.argmin(distances)
#
#     def annotatePlot(X, index):
#         """Create popover label in 3d chart
#
#         Args:
#             X (np.array) - array of points, of shape (numPoints, 3)
#             index (int) - index (into points array X) of item which should be printed
#         Returns:
#             None
#         """
#         # If we have previously displayed another label, remove it first
#         if hasattr(annotatePlot, 'label'):
#             annotatePlot.label.remove()
#         # Get data point from array of points X, at position index
#         x2, y2, _ = proj3d.proj_transform(X[index, 0], X[index, 1], X[index, 2], ax.get_proj())
#         annotatePlot.label = plt.annotate( "case %d" % index,
#             xy = (x2, y2), xytext = (-20, 20), textcoords = 'offset points', ha = 'right', va = 'bottom',
#             bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
#             arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
#         fig.canvas.draw()
#
#     def onMouseMotion(event):
#         """Event that is triggered when mouse is moved. Shows text annotation over data point closest to mouse."""
#         closestIndex = calcClosestDatapoint(X, event)
#         annotatePlot (X, closestIndex)
#
#     fig.canvas.mpl_connect('motion_notify_event', onMouseMotion)  # on mouse motion
#     plt.show()
#
#
# if __name__ == '__main__':
#
#     # automated_design_func()
#
#     # design_name = "sweep_new_3d_with_cost"
#     # load_design(design_name)
#
#     design_name = "sweep_new_3d_with_cost"
#     K = plot_3d(design_name)
#     print(K)
#     visualize3DData(K)
#
#     # directory_names = ["sweep_examples_5", "hyperbolic_filtered_0_01_0_6"]
#     # compare_graph(directory_names)
