import femmt as fmt
import numpy as np
import re
import os
from os import listdir
from os.path import isfile, join
import shutil
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import proj3d
from itertools import product
import logging
from scipy.interpolate import interp1d
import materialdatabase as mdb


material_db = mdb.MaterialDatabase()


def load_from_single_file(working_directory: str, file_name: str):
    """
        Load from a single FEM simulation case for checking the result in detail

        param working_directory: Working directory where all the simulated cases have been saved from the automated design
        :type working_directory: str
        :param file_name: Log file which needs to be simulated (e.g. 'case1344.json')
        :type file_name: str
    """

    file_path = os.path.join(working_directory, "fem_simulation_data", file_name)

    example_results_folder = os.path.join(working_directory, "example_results")
    if not os.path.exists(example_results_folder):
        os.mkdir(example_results_folder)

    working_directory = os.path.join(example_results_folder, "from-file")
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    file_path_dict = {file_name: file_path}
    # After the simulations the sweep can be analyzed
    # This could be done using the FEMMTLogParser:
    log_parser = fmt.FEMMTLogParser(file_path_dict)

    frequency = 0
    current = 0
    for name, data in log_parser.data.items():
        frequency = data.sweeps[0].frequency
        current = data.sweeps[0].windings[0].current.real

    geo = fmt.MagneticComponent.decode_settings_from_log(file_path, working_directory)
    # TODO: Update the decode_settings_from_log as per the new core database

    geo.create_model(freq=frequency, visualize_before=False, save_png=False)

    geo.single_simulation(freq=frequency, current=[current], show_results=True)


def plot_2d(x_value: list, y_value: list, x_label: str, y_label: str, title: str, plot_color: str, annotations: list = None):
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
        :param plot_color: Color of the plot (the colors are based on 'fmt.colors_femmt_default')
        :type annotations: str
    """
    if annotations is None:
        names = [str(x) for x in list(range(len(x_value)))]
    else:
        names = np.array(annotations)

    x_value_str = [str(round(x, 6)) for x in x_value]
    y_value_str = [str(round(y, 6)) for y in y_value]

    fig, ax = plt.subplots()
    fmt.plt.title(title)
    fmt.plt.xlabel(x_label)
    fmt.plt.ylabel(y_label)

    # c = np.random.randint(1, 5, size=len(y_value))
    # norm = plt.Normalize(1, 4)
    # cmap = plt.cm.RdYlGn

    # sc = plt.scatter(x_value, y_value, c=c, s=50, cmap=cmap, norm=norm)

    sc = plt.scatter(x_value, y_value, c='#%02x%02x%02x' % fmt.colors_femmt_default[plot_color])

    annot = ax.annotate("", xy=(0, 0), xytext=(20, 20), textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)

    def update_annot(ind):
        """Create popover annotations in 2d plot"""

        pos = sc.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        text = "{}\n{}: {}\n{}:{}".format(" ".join([names[n] for n in ind["ind"]]), x_label,
                                                " ".join([x_value_str[n] for n in ind["ind"]]), y_label,
                                                " ".join([y_value_str[n] for n in ind["ind"]]))
        annot.set_text(text)
        # annot.get_bbox_patch().set_facecolor(cmap(norm(c[ind["ind"][0]])))
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


def plot_3d(x_value: list, y_value: list, z_value: list, x_label: str, y_label: str, z_label: str,
            title: str, annotations: list, plot_color: str):
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
        :param title: Title of the graph
        :type title: str
        :param annotations: Annotations corresponding to the 3D points
        :type annotations: list
        :param plot_color: Color of the plot (the colors are based on 'fmt.colors_femmt_default')
        :type annotations: str
    """
    names = np.array(annotations)
    num = [re.findall(r'\d+', item) for item in names]
    print(num)
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

    fig = plt.figure()  # figsize=(16, 10)
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_zlabel(z_label)

    # c = np.random.randint(1, 5, size=len(X))
    #
    # norm = plt.Normalize(1, 4)
    # cmap = plt.cm.RdYlGn
    # ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=c, s=10, cmap=cmap, norm=norm, depthshade=False, picker=True)

    ax.scatter(X[:, 0], X[:, 1], X[:, 2], c='#%02x%02x%02x' % fmt.colors_femmt_default[plot_color],
               depthshade=True, picker=True)

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
        total_loss.append(data.total_winding_losses + data.total_core_losses)
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
    It consists of input parameters sweep, filtration, FEM simulation and plotting of the relevant results.
    """

    # copper conductivity (sigma) @ 20 degree celsius (in siemens/meter)
    copper_conductivity = 5.96 * 1e7
    winding_scheme_dict = {'Square': [0.785, fmt.ConductorArrangement.Square],
                           'Hexagonal': [0.907, fmt.ConductorArrangement.Hexagonal]}
    constant_convert_dict = {np.nan: None}

    def __init__(self, working_directory: str, magnetic_component: str, goal_inductance: float, frequency: float,
                 goal_inductance_percent_tolerance: int, winding_scheme: str, peak_current: float, percent_of_b_sat: int,
                 percent_of_total_loss: int, database_core_names: list, database_litz_names: list,
                 solid_conductor_r: list, manual_core_inner_diameter: list, manual_window_h: list, manual_window_w: list,
                 no_of_turns: list, n_air_gaps: list, air_gap_height: list, air_gap_position: list, core_material: list,
                 mult_air_gap_type: list, top_core_insulation: float, bot_core_insulation: float,
                 left_core_insulation: float, right_core_insulation: float, inner_winding_insulation: float,
                 temperature: float, manual_litz_conductor_r: list, manual_litz_strand_r: list,
                 manual_litz_strand_n: list, manual_litz_fill_factor: list):

        """
                :param working_directory: Sets the working directory
                :type working_directory: str
                :param magnetic_component: Sets magnetic component: 'inductor', 'integrated transformer'
                :type magnetic_component: str
                :param goal_inductance: Sets goal inductance for design [in Henry]
                :type goal_inductance: float
                :param frequency: Operating frequency [in Hz]
                :type frequency: float
                :param goal_inductance_percent_tolerance: percent tolerance with respect to goal inductance [in percent]
                :type goal_inductance_percent_tolerance: float
                :param winding_scheme: Winding scheme: 'Square' or 'Hexagonal'
                :type winding_scheme: str
                :param peak_current: Max current amplitude with assumption of sinusoidal current waveform [in Ampere]
                :type peak_current: float
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
                :param manual_core_inner_diameter: Diameter of center leg of the core [in meter]
                :type manual_core_inner_diameter: list
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
                :param core_material: Relative permeability of the core [in F/m]
                :type core_material: list
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
                :param inner_winding_insulation: inner_winding_insulation [in meter]
                :type inner_winding_insulation: float
                :param temperature: core temperature [in degree Celsius]
                :type temperature: float
        """

        self.working_directory = working_directory
        if not os.path.exists(self.working_directory):
            os.mkdir(self.working_directory)
        self.magnetic_component = magnetic_component

        self.goal_inductance = goal_inductance
        self.goal_inductance_percent_tolerance = goal_inductance_percent_tolerance
        self.peak_current = peak_current
        self.percent_of_b_sat = percent_of_b_sat
        self.percent_of_total_loss = percent_of_total_loss
        self.frequency = frequency
        self.temperature = temperature

        # Set core-geometry from core database or/and manual entry
        self.database_core_names = database_core_names
        self.manual_core_inner_diameter = manual_core_inner_diameter
        self.manual_window_h = manual_window_h
        self.manual_window_w = manual_window_w

        # Set winding settings (Solid and Litz winding type)
        self.database_litz_names = database_litz_names
        self.manual_litz_conductor_r = manual_litz_conductor_r
        self.manual_litz_strand_r = manual_litz_strand_r
        self.manual_litz_strand_n = manual_litz_strand_n
        self.manual_litz_fill_factor = manual_litz_fill_factor
        self.solid_conductor_r = solid_conductor_r

        # Set air-gap and core parameters
        self.no_of_turns = no_of_turns
        self.n_air_gaps = n_air_gaps
        self.air_gap_height = air_gap_height
        self.air_gap_position = air_gap_position

        # Set core material
        self.core_material = core_material
        self.core_material_dict = {}

        # Multiple air-gap type ('edge_distributed' and/or 'center_distributed')
        self.mult_air_gap_type = mult_air_gap_type

        # Set windng and insulation data
        self.winding_scheme = winding_scheme
        self.winding_factor = self.winding_scheme_dict[winding_scheme][0]
        self.top_core_insulation = top_core_insulation
        self.bot_core_insulation = bot_core_insulation
        self.left_core_insulation = left_core_insulation
        self.right_core_insulation = right_core_insulation
        self.inner_winding_insulation = inner_winding_insulation

        # Pre-process function to prepare inputs for Reluctance model
        self.core_inner_diameter_list, self.window_h_list, self.window_w_list, self.litz_conductor_r, \
        self.litz_strand_r, self.litz_strand_n, self.litz_fill_factor, self.mu_rel, \
        self.mult_air_gap_type_list = self.input_pre_process()

        # self.min_conductor_r = min(self.litz_conductor_r + self.solid_conductor_r)  # TODO:
        self.conductor_r_list = self.litz_conductor_r + self.solid_conductor_r  # TODO:

        # Call to Reluctance model (Class MagneticCircuit)
        mc = fmt.MagneticCircuit(core_inner_diameter=self.core_inner_diameter_list, window_h=self.window_h_list, window_w=self.window_w_list,
                                 no_of_turns=self.no_of_turns, n_air_gaps=self.n_air_gaps,
                                 air_gap_h=self.air_gap_height, air_gap_position=self.air_gap_position,
                                 mu_rel=self.mu_rel, mult_air_gap_type=self.mult_air_gap_type_list,
                                 air_gap_method='Percent', component_type=self.magnetic_component, sim_type='sweep')
        self.param = mc.get_parameters_position_dict()

        # Filtration of the design cases which are not important
        self.data_matrix_0 = mc.data_matrix

        self.data_matrix_1 = self.filter_goal_inductance(self.data_matrix_0)
        self.data_matrix_2 = self.filter_flux_saturation(self.data_matrix_1)
        self.data_matrix_3 = self.filter_geometry(self.data_matrix_2)
        self.data_matrix_4 = self.filter_losses(self.data_matrix_3)

        #self.plot_volume_loss(self.data_matrix_4)
        #   TODO: remove plot from here
        self.data_matrix_fem = self.data_matrix_4
        #     # FEM_data_matrix = data_matrix_4[np.where(data_matrix_4[:, param["normalized_total_loss"]]
        #     <= ((0.01 / data_matrix_4[:, param["normalized_total_volume"]]) + 0.6))]

    def input_pre_process(self):
        """ Pre-process the user input to prepare lists for reluctance model"""
        all_manual_combinations = list(product(self.manual_core_inner_diameter, self.manual_window_h, self.manual_window_w))
        manual_core_inner_diameter = [item[0] for item in all_manual_combinations]
        manual_window_h = [item[1] for item in all_manual_combinations]
        manual_window_w = [item[2] for item in all_manual_combinations]

        core_db = fmt.core_database()
        db_core_inner_diameter = [core_db[core_name]["core_inner_diameter"] for core_name in self.database_core_names]
        db_window_h = [core_db[core_name]["window_h"] for core_name in self.database_core_names]
        db_window_w = [core_db[core_name]["window_w"] for core_name in self.database_core_names]

        core_inner_diameter_list = db_core_inner_diameter + manual_core_inner_diameter
        window_h_list = db_window_h + manual_window_h
        window_w_list = db_window_w + manual_window_w

        litz_db = fmt.litz_database()
        db_litz_conductor_r = [np.nan if litz_db[litz_name]["conductor_radii"] == ""
                               else litz_db[litz_name]["conductor_radii"] for litz_name in self.database_litz_names]
        db_litz_strand_r = [np.nan if litz_db[litz_name]["strand_radii"] == ""
                            else litz_db[litz_name]["strand_radii"] for litz_name in self.database_litz_names]
        db_litz_strand_n = [np.nan if litz_db[litz_name]["strands_numbers"] == ""
                            else litz_db[litz_name]["strands_numbers"] for litz_name in self.database_litz_names]
        db_litz_fill_factor = [np.nan if litz_db[litz_name]["ff"] == ""
                               else litz_db[litz_name]["ff"] for litz_name in self.database_litz_names]

        litz_conductor_r = db_litz_conductor_r + self.manual_litz_conductor_r
        litz_strand_r = db_litz_strand_r + self.manual_litz_strand_r
        litz_strand_n = db_litz_strand_n + self.manual_litz_strand_n
        litz_fill_factor = db_litz_fill_factor + self.manual_litz_fill_factor

        mu_rel = [material_db.get_material_property(material_name=material_name, property="initial_permeability")
                  for material_name in self.core_material]

        mult_air_gap_type_list = []
        for item in self.mult_air_gap_type:
            if item == 'edge_distributed':
                mult_air_gap_type_list.append(1)
            elif item == 'center_distributed':
                mult_air_gap_type_list.append(2)
            else:
                raise Exception('Wrong string input for multiple air-gap type')

        return core_inner_diameter_list, window_h_list, window_w_list, litz_conductor_r, litz_strand_r, litz_strand_n, \
               litz_fill_factor, mu_rel, mult_air_gap_type_list

    def filter_geometry(self, data_matrix):
        """
        Filter out design cases which are not physical possible based on no_of_turns and winding area

        param data_matrix: Matrix containing the design parameters
        :type data_matrix: ndarray
        """

        final_data_matrix1 = np.zeros((len(data_matrix), len(data_matrix[0]) + 6))
        final_data_matrix2 = np.zeros((len(data_matrix), len(data_matrix[0]) + 6))

        for i in range(len(self.litz_conductor_r)):
            temp_var1 = np.full((len(data_matrix), 1), self.litz_conductor_r[i])
            temp_var2 = np.full((len(data_matrix), 1), self.litz_strand_r[i])
            temp_var3 = np.full((len(data_matrix), 1), self.litz_strand_n[i])
            temp_var4 = np.full((len(data_matrix), 1), self.litz_fill_factor[i])
            temp_var5 = np.full((len(data_matrix), 1), np.nan)
            temp_var6 = np.full((len(data_matrix), 1), self.litz_conductor_r[i])

            temp_data_matrix = self.add_column_to_data_matrix(data_matrix, temp_var1, 'litz_conductor_r')
            temp_data_matrix = self.add_column_to_data_matrix(temp_data_matrix, temp_var2, 'litz_strand_r')
            temp_data_matrix = self.add_column_to_data_matrix(temp_data_matrix, temp_var3, 'litz_strand_n')
            temp_data_matrix = self.add_column_to_data_matrix(temp_data_matrix, temp_var4, 'litz_fill_factor')
            temp_data_matrix = self.add_column_to_data_matrix(temp_data_matrix, temp_var5, 'solid_conductor_r')
            temp_data_matrix = self.add_column_to_data_matrix(temp_data_matrix, temp_var6, 'conductor_radius')

            if i == 0:
                final_data_matrix1 = temp_data_matrix
            else:
                final_data_matrix1 = np.concatenate((final_data_matrix1, temp_data_matrix), axis=0)

        for j in range(len(self.solid_conductor_r)):
            temp_var1 = np.full((len(data_matrix), 1), np.nan)
            temp_var2 = np.full((len(data_matrix), 1), np.nan)
            temp_var3 = np.full((len(data_matrix), 1), np.nan)
            temp_var4 = np.full((len(data_matrix), 1), np.nan)
            temp_var5 = np.full((len(data_matrix), 1), self.solid_conductor_r[j])
            temp_var6 = np.full((len(data_matrix), 1), self.solid_conductor_r[j])

            temp_data_matrix = self.add_column_to_data_matrix(data_matrix, temp_var1, 'litz_conductor_r')       #20
            temp_data_matrix = self.add_column_to_data_matrix(temp_data_matrix, temp_var2, 'litz_strand_r')     #21
            temp_data_matrix = self.add_column_to_data_matrix(temp_data_matrix, temp_var3, 'litz_strand_n')     #22
            temp_data_matrix = self.add_column_to_data_matrix(temp_data_matrix, temp_var4, 'litz_fill_factor')  #23
            temp_data_matrix = self.add_column_to_data_matrix(temp_data_matrix, temp_var5, 'solid_conductor_r') #24
            temp_data_matrix = self.add_column_to_data_matrix(temp_data_matrix, temp_var6, 'conductor_radius')  #25

            if j == 0:
                final_data_matrix2 = temp_data_matrix
            else:
                final_data_matrix2 = np.concatenate((final_data_matrix2, temp_data_matrix), axis=0)

        if not len(self.litz_conductor_r) and len(self.solid_conductor_r):
            final_data_matrix = final_data_matrix2
        elif not len(self.solid_conductor_r) and len(self.litz_conductor_r):
            final_data_matrix = final_data_matrix1
        elif len(self.litz_conductor_r) and len(self.solid_conductor_r):
            final_data_matrix = np.concatenate((final_data_matrix1, final_data_matrix2), axis=0)
        else:
            raise Exception("Please input at least one conductor type")

        window_area = final_data_matrix[:, self.param["window_h"]] * final_data_matrix[:, self.param["window_w"]]
        insulation_area = ((self.left_core_insulation + self.right_core_insulation) *
                           final_data_matrix[:, self.param["window_h"]]) + \
                          ((self.top_core_insulation + self.bot_core_insulation) *
                           (final_data_matrix[:, self.param["window_w"]] -
                            (self.left_core_insulation + self.right_core_insulation)))

        data_matrix = final_data_matrix[
            np.where((final_data_matrix[:, self.param["no_of_turns"]] * np.pi * final_data_matrix[:, self.param["conductor_radius"]] ** 2)
                     < (self.winding_factor * (window_area - insulation_area)))]

        return data_matrix

    def filter_goal_inductance(self, data_matrix):
        """
        Filter out design cases which are in between the given goal inductance tolerance band

        param data_matrix: Matrix containing the design parameters
        :type data_matrix: array
        """
        data_matrix = data_matrix[np.where((data_matrix[:, self.param["inductance"]] >
                                            ((100 - self.goal_inductance_percent_tolerance) / 100) *
                                            self.goal_inductance) &
                                           (data_matrix[:, self.param["inductance"]] <
                                            ((100 + self.goal_inductance_percent_tolerance) / 100) *
                                            self.goal_inductance))]
        return data_matrix

    def filter_flux_saturation(self, data_matrix):
        """
        Filter out design cases based on the maximum magnetic flux allowed in the magnetic core.

        param data_matrix: Matrix containing the design parameters
        :type data_matrix: array
        """
        b_sat_dict = {}
        for material_name in self.core_material:
            b_sat_key = material_db.get_material_property(material_name=material_name, property="initial_permeability")
            b_sat_dict[b_sat_key] = material_db.get_material_property(material_name=material_name,
                                                                      property="max_flux_density")
            self.core_material_dict[b_sat_key] = material_name
        print(self.core_material_dict)

        # Creating B_saturated array corresponding to the material type
        b_sat = np.zeros((len(data_matrix), 1))
        for index in range(len(data_matrix)):
            b_sat[index] = b_sat_dict[data_matrix[index, self.param["mu_rel"]]]

        # flux_max = L * i_max / N
        total_flux_max = (data_matrix[:, self.param["inductance"]] * self.peak_current) / data_matrix[:,
                                                                                         self.param["no_of_turns"]]
        b_max_center = total_flux_max / data_matrix[:, self.param["center_leg_area"]]
        b_max_middle = total_flux_max / (
                np.pi * data_matrix[:, self.param["core_inner_diameter"]] * data_matrix[:, self.param["core_h_middle"]])
        b_max_outer = total_flux_max / data_matrix[:, self.param["outer_leg_area"]]

        data_matrix = self.add_column_to_data_matrix(data_matrix, total_flux_max, 'total_flux_max')     # 16
        data_matrix = self.add_column_to_data_matrix(data_matrix, b_max_center, 'b_max_center')         # 17
        data_matrix = self.add_column_to_data_matrix(data_matrix, b_max_middle, 'b_max_middle')         # 18
        data_matrix = self.add_column_to_data_matrix(data_matrix, b_max_outer, 'b_max_outer')           # 19

        data_matrix_temp = np.zeros((0, len(data_matrix[0])))
        for index in range(len(data_matrix)):
            if (data_matrix[index, self.param["b_max_center"]] < (self.percent_of_b_sat / 100) * b_sat[index]) & \
                    (data_matrix[index, self.param["b_max_center"]] < (self.percent_of_b_sat / 100) * b_sat[index]) & \
                    (data_matrix[index, self.param["b_max_outer"]] < (self.percent_of_b_sat / 100) * b_sat[index]):
                data_matrix_temp = np.vstack([data_matrix_temp, data_matrix[index, :]])

        return data_matrix_temp

    def filter_losses(self, data_matrix):
        """
       Filter out design cases based on the calculated hysteresis and DC loss

       param data_matrix: Matrix containing the design parameters
       :type data_matrix: ndarray
        """
        # Filter out data-matrix according to calculated hysteresis loss + DC winding loss
        # Volume chosen as per "Masterthesis_Till_Piepenbrock" pg-45
        mu_imag_dict = {}
        counter = 0
        for material_name in self.core_material:
            mu_imag_key = material_db.get_material_property(material_name=material_name,
                                                            property="initial_permeability")
            mu_imag_dict[mu_imag_key] = counter
            counter = counter + 1

        material_data_list = [material_db.permeability_data_to_pro_file(self.temperature, self.frequency,
                                                                        material_name, "manufacturer_datasheet")
                              for material_name in self.core_material]
        print(material_data_list)

        mu_imag_interpol_func = [interp1d(material_data_list[i][0], material_data_list[i][1], kind="cubic") for i in
                                 range(len(self.core_material))]

        # Creating mu_imag array corresponding to the material type
        mu_imag = np.zeros(len(data_matrix))
        for index in range(len(data_matrix)):
            mu_imag[index] = mu_imag_interpol_func[mu_imag_dict[data_matrix[index, self.param["mu_rel"]]]] \
                (data_matrix[index, self.param["b_max_center"]])

        volume_center = (np.pi * (data_matrix[:, self.param["core_inner_diameter"]] / 2) ** 2) * \
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
                                    (data_matrix[:, self.param["r_inner"]] * 2) / data_matrix[:, self.param["core_inner_diameter"]])
        P_hyst_density_outer = P_hyst_outer * volume_outer
        total_hyst_loss = P_hyst_density_center + (2 * P_hyst_density_middle) + P_hyst_density_outer

        # Winding loss (only DC loss)
        Resistance = (data_matrix[:, self.param["no_of_turns"]] * 2 * np.pi *
                      (data_matrix[:, self.param["core_inner_diameter"]] / 2 + data_matrix[:, self.param["conductor_radius"]])) / \
                     (self.copper_conductivity * (np.pi * (data_matrix[:, self.param["conductor_radius"]] ** 2)))

        DC_loss = ((self.peak_current ** 2) / 2) * Resistance       # Assuming sinusoidal current waveform

        total_loss = DC_loss + total_hyst_loss
        max_total_loss = max(total_loss)
        normalized_total_loss = total_loss / max_total_loss

        total_volume = np.pi * (data_matrix[:, self.param["r_outer"]] ** 2) * data_matrix[:, self.param["core_h"]]
        max_volume = max(total_volume)
        normalized_total_volume = total_volume / max_volume

        data_matrix = self.add_column_to_data_matrix(data_matrix, total_hyst_loss, 'total_hyst_loss')               # 26
        data_matrix = self.add_column_to_data_matrix(data_matrix, DC_loss, 'DC_loss')                               # 27
        data_matrix = self.add_column_to_data_matrix(data_matrix, total_loss, 'total_loss')                         # 28
        data_matrix = self.add_column_to_data_matrix(data_matrix, normalized_total_loss, 'normalized_total_loss')   # 29
        data_matrix = self.add_column_to_data_matrix(data_matrix, total_volume, 'total_volume')                     # 30
        data_matrix = self.add_column_to_data_matrix(data_matrix, normalized_total_volume, 'normalized_total_volume')   # 31

        data_matrix = data_matrix[data_matrix[:, self.param["total_loss"]].argsort()]
        data_matrix = data_matrix[0:int((self.percent_of_total_loss / 100) * len(data_matrix)), :]

        # test_var = data_matrix[np.where(data_matrix[:, self.param["total_volume"]] >= 0.0003574349)]
        # print(test_var)
        # test_var_2 = test_var[np.where(test_var[:, self.param["total_volume"]] <= 0.0003574350)]
        # print(test_var_2)

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
        for i in range(len(self.data_matrix_fem)):

            # MagneticComponent class object
            geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor,
                                        working_directory=working_directory, silent=True)

            core = fmt.Core(core_inner_diameter=self.data_matrix_fem[i, self.param["core_inner_diameter"]],
                            window_w=self.data_matrix_fem[i, self.param["window_w"]],
                            window_h=self.data_matrix_fem[i, self.param["window_h"]],
                            material=self.core_material_dict[self.data_matrix_fem[i, self.param["mu_rel"]]],
                            temperature=self.temperature, frequency=self.frequency, datasource="manufacturer_datasheet")
                            # material="95_100")
                            # mu_rel=self.data_matrix_fem[i, self.param["mu_rel"]], phi_mu_deg=10,
                            # sigma=0.5)
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
            insulation.add_winding_insulations([self.inner_winding_insulation], 0.0001)
            geo.set_insulation(insulation)

            # 5. create winding window and virtual winding windows (vww)
            winding_window = fmt.WindingWindow(core, insulation)
            vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

            # 6. create conductor and set parameters: use solid wires
            winding = fmt.Conductor(0, fmt.Conductivity.Copper)
            if np.isnan(self.data_matrix_fem[i, self.param["solid_conductor_r"]]):

                winding.set_litz_round_conductor(conductor_radius=None if np.isnan(self.data_matrix_fem[i, self.param["litz_conductor_r"]]) else self.data_matrix_fem[i, self.param["litz_conductor_r"]],
                                                 number_strands=None if np.isnan(self.data_matrix_fem[i, self.param["litz_strand_n"]]) else self.data_matrix_fem[i, self.param["litz_strand_n"]],
                                                 strand_radius=None if np.isnan(self.data_matrix_fem[i, self.param["litz_strand_r"]]) else self.data_matrix_fem[i, self.param["litz_strand_r"]],
                                                 fill_factor=None if np.isnan(self.data_matrix_fem[i, self.param["litz_fill_factor"]]) else self.data_matrix_fem[i, self.param["litz_fill_factor"]],
                                                 conductor_arrangement=self.winding_scheme_dict[self.winding_scheme][1])
            else:
                winding.set_solid_round_conductor(conductor_radius=self.data_matrix_fem[i, self.param["solid_conductor_r"]],
                                                  conductor_arrangement=self.winding_scheme_dict[self.winding_scheme][1])

            # 7. add conductor to vww and add winding window to MagneticComponent
            vww.set_winding(winding, int(self.data_matrix_fem[i, self.param["no_of_turns"]]), None)
            geo.set_winding_window(winding_window)

            try:
                # 5. create the model
                geo.create_model(freq=self.frequency, visualize_before=False, save_png=False)

                # 6. start simulation
                geo.single_simulation(freq=self.frequency, current=[self.peak_current], show_results=False)

                shutil.copy2(src_path, data_folder)
                old_filename = os.path.join(data_folder, "log_electro_magnetic.json")
                new_filename = os.path.join(data_folder, f"case{i}.json")
                os.rename(old_filename, new_filename)
                data_files.append(new_filename)
                file_names.append(f"case{i}")
                print(f"Case {i} of {len(self.data_matrix_fem)} completed")

            except (Exception,) as e:
                print("next iteration")
                logging.exception(e)

    def add_column_to_data_matrix(self, data_matrix, column_value, column_name: str):
        """
        Adds column to the given matrix

        param data_matrix: Matrix containing the design parameters
       :type data_matrix: ndarray
       :param column_value: Column to be added
       :type column_value: ndarray
       :param column_name: Identifier of the column
       :type column_name: str
        """
        size = len(data_matrix[0])
        data_matrix = np.hstack((data_matrix, np.reshape(column_value, (len(column_value), 1))))
        self.param[column_name] = size

        return data_matrix

    # def plot_volume_loss(self, x_value: list, y_value: list, x_label: str, y_label: str,
    #         title: str, plot_color: str):
    #     """
    #     Plots estimated normalised volume vs loss graph from reluctance model results
    #
    #     param data_matrix: Matrix containing the design parameters
    #    :type data_matrix: array
    #     """
    #     fig, ax = fmt.plt.subplots()  # Create a figure containing a single axes.
    #     # fmt.plt.title("Normalised volume vs Normalised losses")
    #     # fmt.plt.xlabel("Normalised Volume")
    #     # fmt.plt.ylabel("Normalised Losses")
    #     # ax.plot(data_matrix[:, self.param['normalized_total_volume']],
    #     #         data_matrix[:, self.param['normalized_total_loss']], 'o')
    #     fmt.plt.title(title)
    #     fmt.plt.xlabel(x_label)
    #     fmt.plt.ylabel(y_label)
    #     ax.plot(x_value, y_value, 'o', color='#%02x%02x%02x' % fmt.colors_femmt_default[plot_color])
    #     ax.grid()
    #     fmt.plt.show()


if __name__ == '__main__':
    ad = AutomatedDesign(working_directory='D:/Personal_data/MS_Paderborn/Sem4/Project_2/2022-11-26_fem_simulation_data',
                         magnetic_component='inductor',
                         goal_inductance=120 * 1e-6,
                         frequency=100000,
                         goal_inductance_percent_tolerance=10,
                         winding_scheme='Square',
                         peak_current=8,
                         percent_of_b_sat=70,
                         percent_of_total_loss=30,
                         database_core_names=[],
                         database_litz_names=['1.5x105x0.1', "1.4x200x0.071"],
                         solid_conductor_r=[],  # 0.0013
                         manual_core_inner_diameter=list(np.linspace(0.005, 0.05, 10)),
                         manual_window_h=list(np.linspace(0.01, 0.08, 5)),
                         manual_window_w=list(np.linspace(0.005, 0.04, 10)),
                         no_of_turns=[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
                         n_air_gaps=[1, 2],
                         air_gap_height=list(np.linspace(0.0001, 0.0005, 5)),
                         air_gap_position=list(np.linspace(20, 80, 2)),
                         core_material=['N95'],
                         mult_air_gap_type=['center_distributed'],
                         top_core_insulation=0.001,
                         bot_core_insulation=0.001,
                         left_core_insulation=0.001,
                         right_core_insulation=0.001,
                         inner_winding_insulation=0.0005,
                         temperature=100.0,
                         manual_litz_conductor_r=[],
                         manual_litz_strand_r=[],
                         manual_litz_strand_n=[],
                         manual_litz_fill_factor=[])

    plot_2d(x_value=ad.data_matrix_fem[:, ad.param["total_volume"]],
            y_value=ad.data_matrix_fem[:, ad.param["total_loss"]],
            x_label='Volume / m\u00b3', y_label='Loss / W', title='Volume vs Loss', plot_color='red')

    ad.fem_simulation()

    inductance, total_loss, total_volume, total_cost, annotation_list = load_design \
        (working_directory='D:/Personal_data/MS_Paderborn/Sem4/Project_2/2022-11-26_fem_simulation_data')

    plot_2d(x_value=total_volume, y_value=total_loss, x_label='Volume / m\u00b3', y_label='Loss / W',
            title='Volume vs Loss', annotations=annotation_list, plot_color='red')
    # plot_2d(x_value=total_volume, y_value=total_cost, x_label='Volume / m\u00b3', y_label='Cost / \u20ac',
    #         title='Volume vs Cost', annotations=annotation_list, plot_color='red')
    # plot_3d(x_value=total_volume, y_value=total_loss, z_value=total_cost, x_label='Volume / m\u00b3',
    #         y_label='Loss / W', z_label='Cost / \u20ac', title='Volume vs Loss vs Cost',
    #         annotations=annotation_list, plot_color='red')

    # load_from_single_file(working_directory='D:/Personal_data/MS_Paderborn/Sem4/Project_2/2022-11-20_fem_simulation_data',
    #                file_name='case1922.json')

    # plot_2d(x_value=total_volume, y_value=inductance, x_label='Volume / m\u00b3', y_label='Inductance / H',
    #         title='Volume vs Inductance', annotations=annotation_list, plot_color='red')


