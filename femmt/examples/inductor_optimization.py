"""Optimize an inductor using the reluctance model and FEM simulations."""
# python libraries
import json
import csv
import re
import os
import shutil
from itertools import product
# import logging
import inspect
import time

# 3rd party libraries
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import proj3d
from scipy.interpolate import interp1d
import materialdatabase as mdb

# femmt libraries
import femmt as fmt


material_db = mdb.MaterialDatabase()


def load_from_single_file(working_directory: str, file_name: str):
    """
    Load from a single FEM simulation case for checking the result in detail.

    param working_directory: Working directory where all the simulated cases have been saved
    from the automated design
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
    for _, data in log_parser.data.items():
        frequency = data.sweeps[0].frequency
        current = data.sweeps[0].windings[0].current.real

    geo = fmt.MagneticComponent.decode_settings_from_log(file_path, working_directory)

    geo.create_model(freq=frequency, pre_visualize_geometry=False, save_png=False)

    geo.single_simulation(freq=frequency, current=[current], show_fem_simulation_results=True)


def plot_2d(x_value: list, y_value: list, x_label: str, y_label: str, title: str, plot_color: str, z_value: list = None,
            z_label: str = None, inductance_value: list = None, annotations: list = None):
    """
    Visualize data in 2d plot with popover next to mouse position.

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
    :param inductance_value: Data points for inductance value corresponding to the (x, y, z): (Optional)
    :type inductance_value: list
    :param annotations: Annotations corresponding to the 3D points
    :type annotations: list
    :param plot_color: Color of the plot (the colors are based on 'fmt.colors_femmt_default')
    :type annotations: str
    """
    if annotations is None:
        names = [str(x) for x in list(range(len(x_value)))]
    else:
        temp_var = [int(x) for x in annotations]
        names = [str(x) for x in temp_var]

    if inductance_value is not None:
        l_label = 'L / H'

    if z_value is not None:
        z_value_str = [str(round(z, 3)) for z in z_value]

    if inductance_value is not None:
        l_value_str = [str(round(i_inductance, 6)) for i_inductance in inductance_value]

    x_value_str = [str(round(x, 6)) for x in x_value]
    y_value_str = [str(round(y, 3)) for y in y_value]

    fig, ax = plt.subplots()
    fmt.plt.title(title)
    fmt.plt.xlabel(x_label)
    fmt.plt.ylabel(y_label)

    # c = np.random.randint(1, 5, size=len(y_value))
    # norm = plt.Normalize(1, 4)
    # cmap = plt.cm.RdYlGn

    # sc = plt.scatter(x_value, y_value, c=c, s=50, cmap=cmap, norm=norm)
    if z_value is None:
        sc = plt.scatter(x_value, y_value, c='#%02x%02x%02x' % fmt.colors_femmt_default[plot_color])
    else:
        sc = plt.scatter(x_value, y_value, c=z_value, cmap=plot_color)
        cbar = plt.colorbar(sc)
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel(z_label, rotation=270)

    annot = ax.annotate("", xy=(0, 0), xytext=(20, 20), textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)

    def update_annot(ind):
        """Create popover annotations in 2d plot."""
        pos = sc.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        text = ""
        if z_label is None and inductance_value is None:
            text = "case: {}\n{}: {}\n{}:{}".\
                format(" ".join([names[n] for n in ind["ind"]]),
                       x_label, " ".join([x_value_str[n] for n in ind["ind"]]),
                       y_label, " ".join([y_value_str[n] for n in ind["ind"]]))
        elif z_label is not None and inductance_value is None:
            text = "case: {}\n{}: {}\n{}:{}\n{}:{}". \
                format(" ".join([names[n] for n in ind["ind"]]),
                       x_label, " ".join([x_value_str[n] for n in ind["ind"]]),
                       y_label, " ".join([y_value_str[n] for n in ind["ind"]]),
                       z_label, " ".join([z_value_str[n] for n in ind["ind"]]))
        elif z_label is None and inductance_value is not None:
            text = "case: {}\n{}: {}\n{}:{}\n{}:{}". \
                format(" ".join([names[n] for n in ind["ind"]]),
                       x_label, " ".join([x_value_str[n] for n in ind["ind"]]),
                       y_label, " ".join([y_value_str[n] for n in ind["ind"]]),
                       l_label, " ".join([l_value_str[n] for n in ind["ind"]]))
        else:
            text = "case: {}\n{}: {}\n{}:{}\n{}:{}\n{}:{}".\
                format(" ".join([names[n] for n in ind["ind"]]),
                       x_label, " ".join([x_value_str[n] for n in ind["ind"]]),
                       y_label, " ".join([y_value_str[n] for n in ind["ind"]]),
                       z_label, " ".join([z_value_str[n] for n in ind["ind"]]),
                       l_label, " ".join([l_value_str[n] for n in ind["ind"]]))
        annot.set_text(text)
        # annot.get_bbox_patch().set_facecolor(cmap(norm(c[ind["ind"][0]])))
        annot.get_bbox_patch().set_alpha(0.8)

    def hover(event):
        """Event that is triggered when mouse is hovered. Shows text annotation over data point closest to mouse."""
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
    # print(num)
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
        """Return distance between mouse position and given data point.

        Args:
        ----
            point (np.array): np.array of shape (3,), with x,y,z in data coordinates
            event (MouseEvent): mouse event (which contains mouse position in .x and .xdata)

        Returns:
        -------
            distance (np.float64): distance (in screen coordinates) between mouse pos and data point
        """
        assert point.shape == (3,), "distance: point.shape is wrong: %s, must be (3,)" % point.shape

        # Project 3d data space to 2d data space
        x2, y2, _ = proj3d.proj_transform(point[0], point[1], point[2], plt.gca().get_proj())
        # Convert 2d data space to 2d screen space
        x3, y3 = ax.transData.transform((x2, y2))

        return np.sqrt((x3 - event.x) ** 2 + (y3 - event.y) ** 2)

    def calc_closest_datapoint(array_of_points, mouse_event):
        """Calculate which data point is closest to the mouse position.

        Args:
        ----
        array_of_points (np.array): array of points, of shape (numPoints, 3)
        mouse_event (MouseEvent): mouse event (containing mouse position)

        Returns:
        -------
        smallestIndex (int) - the index (into the array of points X) of the element closest to the mouse position
        """
        X_modified = array_of_points[np.all(~np.isnan(array_of_points), axis=1), 0:3]
        distances = [distance(X_modified[i, 0:3], mouse_event) for i in range(X_modified.shape[0])]
        return np.argmin(distances)

    def annotate_plot(array_of_points, index):
        """Create popover label in 3d chart.

        Args:
        ----
        array_of_points (np.array): array of points, of shape (numPoints, 3)
        index (int): index (into points array X) of item which should be printed
        Returns:
        -------
        None
        """
        # If we have previously displayed another label, remove it first
        if hasattr(annotate_plot, 'label'):
            annotate_plot.label.remove()
        # Get data point from array of points X, at position index
        x2, y2, _ = proj3d.proj_transform(array_of_points[index, 0], array_of_points[index, 1], array_of_points[index, 2], ax.get_proj())
        annotate_plot.label = plt.annotate(
            f'case{index}\nVolume:{x_value_str[index]}\nLoss:{y_value_str[index]}\nCost:{z_value_str[index]}',
            xy=(x2, y2), xytext=(-20, 20), textcoords='offset points', ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
            arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
        fig.canvas.draw()

    def on_mouse_motion(event):
        """Event that is triggered when mouse is moved. Shows text annotation over data point closest to mouse."""
        closestIndex = calc_closest_datapoint(X, event)
        annotate_plot(X, closestIndex)

    fig.canvas.mpl_connect('motion_notify_event', on_mouse_motion)  # on mouse motion
    plt.show()


def filter_after_fem(inductance: list, total_loss: list, total_volume: list, total_cost: list,
                     annotation_list: list, goal_inductance: float, percent_tolerance: int):
    """
    Filter the FEM simulated data based on the inductance tolerance.

    param inductance: Inductance read from FEM simulation cases[in Henry]
    :type inductance: list
    :param total_loss: total_loss (hysteresis + winding) read from FEM simulation cases[in Watt]
    :type total_loss: list
    :param total_volume: total_volume of core read from FEM simulation cases [in cubic meter]
    :type total_volume: list
    :param total_cost: total_cost (core material + winding material) read form FEM simulation cases [in Euro]
    :type total_cost: list
    :param annotation_list: String corresponding to the Case number which displays during plot
    when mouse is hovered over the point
    :type annotation_list: list
    :param goal_inductance: Sets goal inductance for final design [in Henry]
    :type goal_inductance: float
    :param percent_tolerance: percent tolerance with respect to goal inductance [in percent]
    :type percent_tolerance: int
    """
    inductance = np.array(inductance)
    total_loss = np.array(total_loss)
    total_volume = np.array(total_volume)
    total_cost = np.array(total_cost)

    names = np.array(annotation_list)
    num = [re.findall(r'\d+', item) for item in names]
    case_num_list = [int(item[0]) for item in num]

    data_array = np.column_stack((inductance, total_volume, total_loss, total_cost, case_num_list))
    data_array = data_array[np.where((data_array[:, 0] > ((100 - percent_tolerance) / 100) * goal_inductance) & \
                                     (data_array[:, 0] < ((100 + percent_tolerance) / 100) * goal_inductance))]

    return data_array


def load_fem_simulation_results(working_directory: str):
    """
    Load FEM simulation results from given working directory.

    param working_directory: Sets the working directory
    :type fem_simulation_results_directory: str
    """
    working_directories = []
    labels = []
    fem_simulation_results_directory = os.path.join(working_directory, 'fem_simulation_results')
    print("##########################")
    print(f"{fem_simulation_results_directory =}")
    print("##########################")
    file_names = [f for f in os.listdir(fem_simulation_results_directory) if os.path.isfile(os.path.join(fem_simulation_results_directory, f))]

    counter = 0
    for name in file_names:
        temp_var = os.path.join(fem_simulation_results_directory, name)
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
    for _, data in log_parser.data.items():
        inductivities.append(data.sweeps[0].windings[0].flux_over_current)
        total_loss.append(data.total_core_losses + data.total_winding_losses)
        total_volume.append(data.core_2daxi_total_volume)
        total_cost.append(data.total_cost)

    real_inductance = []
    for i in range(len(total_loss)):
        real_inductance.append(inductivities[i].real)

    # read target values
    automated_design_settings_file = os.path.join(working_directory, "automated_design_settings.json")
    with open(automated_design_settings_file, "r") as fd:
        automated_design_settings = json.loads(fd.read())

    return real_inductance, total_loss, total_volume, total_cost, labels, automated_design_settings


class AutomatedDesign:
    """
    AutomatedDesign class implements brute force optimization for magnetic component design.

    It consists of input parameters sweep, filtration, FEM simulation and plotting of the relevant results.
    """

    # copper conductivity (sigma) @ 20 degree celsius (in siemens/meter)
    copper_conductivity = 5.96 * 1e7
    winding_scheme_dict = {'Square': [0.785, fmt.ConductorArrangement.Square],
                           'Hexagonal': [0.907, fmt.ConductorArrangement.Hexagonal]}
    component_type_dict = {'inductor': fmt.ComponentType.Inductor,
                           'integrated_transformer': fmt.ComponentType.IntegratedTransformer}

    def __init__(self, working_directory: str,
                 magnetic_component: str,
                 target_inductance: float,
                 frequency: float,
                 target_inductance_percent_tolerance: int,
                 winding_scheme: str,
                 peak_current: float,
                 percent_of_flux_density_saturation: int,
                 percent_of_total_loss: int,
                 database_core_names: list,
                 database_litz_names: list,
                 solid_conductor_r: list,
                 manual_core_inner_diameter: list,
                 manual_window_h: list,
                 manual_window_w: list,
                 no_of_turns: list,
                 n_air_gaps: list,
                 air_gap_height: list,
                 air_gap_position: list,
                 core_material: list,
                 mult_air_gap_type: list,
                 top_core_insulation: float,
                 bot_core_insulation: float,
                 left_core_insulation: float,
                 right_core_insulation: float,
                 inner_winding_insulation: float,
                 temperature: float,
                 manual_litz_conductor_r: list,
                 manual_litz_strand_r: list,
                 manual_litz_strand_n: list,
                 manual_litz_fill_factor: list):
        """
        Initialize the automated design.

        :param working_directory: Sets the working directory
        :type working_directory: str
        :param magnetic_component: Sets magnetic component: 'inductor', 'integrated transformer'
        :type magnetic_component: str
        :param target_inductance: Sets goal inductance for design [in Henry]
        :type target_inductance: float
        :param frequency: Operating frequency [in Hz]
        :type frequency: float
        :param target_inductance_percent_tolerance: percent tolerance with respect to goal inductance [in percent]
        :type target_inductance_percent_tolerance: float
        :param winding_scheme: Winding scheme: 'Square' or 'Hexagonal'
        :type winding_scheme: str
        :param peak_current: Max current amplitude with assumption of sinusoidal current waveform [in Ampere]
        :type peak_current: float
        :param percent_of_flux_density_saturation: percent of saturation of magnetic flux density [in percent]
        :type percent_of_flux_density_saturation: int
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
        self.set_up_folder_structure(working_directory)

        self.magnetic_component = magnetic_component

        self.goal_inductance = target_inductance
        self.goal_inductance_percent_tolerance = target_inductance_percent_tolerance
        self.peak_current = peak_current
        self.percent_of_b_sat = percent_of_flux_density_saturation
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

        # Set winding and insulation data
        self.winding_scheme = winding_scheme
        self.winding_factor = self.winding_scheme_dict[winding_scheme][0]
        self.top_core_insulation = top_core_insulation
        self.bot_core_insulation = bot_core_insulation
        self.left_core_insulation = left_core_insulation
        self.right_core_insulation = right_core_insulation
        self.inner_winding_insulation = inner_winding_insulation

        # Pre-process function to prepare inputs for Reluctance model
        self.core_inner_diameter_list, self.window_h_list, self.window_w_list, self.litz_conductor_r, \
            self.litz_strand_r, self.litz_strand_n, self.litz_fill_factor, self.mu_r_abs, \
            self.mult_air_gap_type_list = self.input_pre_process()

        # Call to Reluctance model (Class MagneticCircuit)
        mc = fmt.MagneticCircuit(core_inner_diameter=self.core_inner_diameter_list, window_h=self.window_h_list, window_w=self.window_w_list,
                                 no_of_turns=self.no_of_turns, n_air_gaps=self.n_air_gaps,
                                 air_gap_h=self.air_gap_height, air_gap_position=self.air_gap_position,
                                 mu_r_abs=self.mu_r_abs, mult_air_gap_type=self.mult_air_gap_type_list,
                                 air_gap_method='Percent', component_type=self.magnetic_component, sim_type='sweep')
        self.param = mc.get_parameters_position_dict()

        # Filtration of the design cases which are not important
        self.data_matrix_0 = mc.data_matrix
        self.data_matrix_1 = self.filter_reluctance_target_inductance(self.data_matrix_0)
        self.data_matrix_2 = self.filter_reluctance_flux_saturation(self.data_matrix_1)
        self.data_matrix_3 = self.filter_reluctance_winding_window(self.data_matrix_2)
        self.data_matrix_4 = self.filter_reluctance_losses(self.data_matrix_3)
        self.data_matrix_5 = self.filter_reluctance_pareto_front_tolerance(self.data_matrix_4)
        self.data_matrix_fem = self.data_matrix_5

    def set_up_folder_structure(self, working_directory):
        """Set up the folder structure for the inductor optimization."""
        if working_directory is None:
            caller_filename = inspect.stack()[1].filename
            working_directory = os.path.join(os.path.dirname(caller_filename),
                                             "inductor_optimization")

        # generate new and empty working directory
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # set up folders for optimization
        self.optimization_working_directory = working_directory
        self.femmt_working_directory = os.path.join(self.optimization_working_directory, "femmt_simulation")
        self.inductor_reluctance_model_results_directory = os.path.join(
            self.optimization_working_directory, "reluctance_model_results")
        self.inductor_fem_simulations_results_directory = os.path.join(
            self.optimization_working_directory, "fem_simulation_results")
        self.inductor_optimization_input_parameters_file = os.path.join(
            self.optimization_working_directory, "optimization_input_parameters.json")

        fmt.create_folders(self.optimization_working_directory, self.inductor_reluctance_model_results_directory,
                           self.inductor_fem_simulations_results_directory, self.femmt_working_directory)

    def input_pre_process(self):
        """Pre-process the user input to prepare lists for reluctance model."""
        # Creating all possible combinations from the given manual geometry parameters
        all_manual_combinations = list(product(self.manual_core_inner_diameter, self.manual_window_h, self.manual_window_w))
        manual_core_inner_diameter = [item[0] for item in all_manual_combinations]
        manual_window_h = [item[1] for item in all_manual_combinations]
        manual_window_w = [item[2] for item in all_manual_combinations]

        # Segregating core geometry parameters from the core database and saving in lists
        core_db = fmt.core_database()
        db_core_inner_diameter = [core_db[core_name]["core_inner_diameter"] for core_name in self.database_core_names]
        db_window_h = [core_db[core_name]["window_h"] for core_name in self.database_core_names]
        db_window_w = [core_db[core_name]["window_w"] for core_name in self.database_core_names]

        core_inner_diameter_list = db_core_inner_diameter + manual_core_inner_diameter
        window_h_list = db_window_h + manual_window_h
        window_w_list = db_window_w + manual_window_w

        # Segregating core geometry parameters from the litz database and saving in lists
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

        # List of initial permeability extracted from core material database
        mu_r_abs = [material_db.get_material_attribute(material_name=material_name, attribute="initial_permeability")
                    for material_name in self.core_material]

        # Mapping mult_air_gap_type string to float value (To include it in numpy float array easily)
        mult_air_gap_type_list = []
        for item in self.mult_air_gap_type:
            if item == 'edge_distributed':
                mult_air_gap_type_list.append(1)
            elif item == 'center_distributed':
                mult_air_gap_type_list.append(2)
            else:
                raise Exception('Wrong string input for multiple air-gap type')

        return core_inner_diameter_list, window_h_list, window_w_list, litz_conductor_r, litz_strand_r, litz_strand_n, \
            litz_fill_factor, mu_r_abs, mult_air_gap_type_list

    def filter_reluctance_target_inductance(self, data_matrix):
        """
        Filter out design cases which are in between the given goal inductance tolerance band.

        :param data_matrix: Matrix containing the design parameters
        :type data_matrix: ndarray
        """
        data_matrix = data_matrix[np.where((data_matrix[:, self.param["inductance"]] > ((100 - self.goal_inductance_percent_tolerance) / 100) * \
                                            self.goal_inductance) & (data_matrix[:, self.param["inductance"]] < \
                                                                     ((100 + self.goal_inductance_percent_tolerance) / 100) * self.goal_inductance))]
        return data_matrix

    def filter_reluctance_flux_saturation(self, data_matrix):
        """
        Filter out design cases based on the maximum magnetic flux allowed in the magnetic core.

        param data_matrix: Matrix containing the design parameters
        :type data_matrix: ndarray
        """
        # Dictionary to store {initial_permeability: 'Material_name'} to map material_name during FEM iteration
        b_sat_dict = {}
        for material_name in self.core_material:
            b_sat_key = material_db.get_material_attribute(material_name=material_name, attribute="initial_permeability")
            b_sat_dict[b_sat_key] = material_db.get_material_attribute(material_name=material_name, attribute="max_flux_density")
            self.core_material_dict[b_sat_key] = material_name
        # print(self.core_material_dict)

        # Creating B_saturated array corresponding to the material type
        b_sat = np.zeros((len(data_matrix), 1))
        for index in range(len(data_matrix)):
            b_sat[index] = b_sat_dict[data_matrix[index, self.param["mu_r_abs"]]]

        # flux_max = L * i_max / N
        total_flux_max = (data_matrix[:, self.param["inductance"]] * self.peak_current) / data_matrix[:, self.param["no_of_turns"]]
        b_max_center = total_flux_max / data_matrix[:, self.param["center_leg_area"]]
        b_max_middle = total_flux_max / (np.pi * data_matrix[:, self.param["core_inner_diameter"]] * data_matrix[:, self.param["core_h_middle"]])
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

    def filter_reluctance_winding_window(self, data_matrix):
        """
        Filter out design cases which are not physical possible based on no_of_turns and winding area.

        :param data_matrix: Matrix containing the design parameters
        :type data_matrix: ndarray
        """
        final_data_matrix1 = np.zeros((len(data_matrix), len(data_matrix[0]) + 6))
        final_data_matrix2 = np.zeros((len(data_matrix), len(data_matrix[0]) + 6))

        # Adds litz and solid core details to data_matrix with all combinations possible and prepare data matrix for
        # FEM simulation
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

            temp_data_matrix = self.add_column_to_data_matrix(data_matrix, temp_var1, 'litz_conductor_r')        # 20
            temp_data_matrix = self.add_column_to_data_matrix(temp_data_matrix, temp_var2, 'litz_strand_r')      # 21
            temp_data_matrix = self.add_column_to_data_matrix(temp_data_matrix, temp_var3, 'litz_strand_n')      # 22
            temp_data_matrix = self.add_column_to_data_matrix(temp_data_matrix, temp_var4, 'litz_fill_factor')   # 23
            temp_data_matrix = self.add_column_to_data_matrix(temp_data_matrix, temp_var5, 'solid_conductor_r')  # 24
            temp_data_matrix = self.add_column_to_data_matrix(temp_data_matrix, temp_var6, 'conductor_radius')   # 25

            if j == 0:
                final_data_matrix2 = temp_data_matrix
            else:
                final_data_matrix2 = np.concatenate((final_data_matrix2, temp_data_matrix), axis=0)

        # Handles all the different cases of litz or/and solid conductor inputs
        if not len(self.litz_conductor_r) and len(self.solid_conductor_r):
            final_data_matrix = final_data_matrix2
        elif not len(self.solid_conductor_r) and len(self.litz_conductor_r):
            final_data_matrix = final_data_matrix1
        elif len(self.litz_conductor_r) and len(self.solid_conductor_r):
            final_data_matrix = np.concatenate((final_data_matrix1, final_data_matrix2), axis=0)
        else:
            raise Exception("Please input at least one conductor type")

        # Filter based on the geometry
        window_area = final_data_matrix[:, self.param["window_h"]] * final_data_matrix[:, self.param["window_w"]]
        insulation_area = ((self.left_core_insulation + self.right_core_insulation) * final_data_matrix[:, self.param["window_h"]]) + \
                          ((self.top_core_insulation + self.bot_core_insulation) * \
                           (final_data_matrix[:, self.param["window_w"]] - (self.left_core_insulation + self.right_core_insulation)))

        data_matrix = final_data_matrix[
            np.where((final_data_matrix[:, self.param["no_of_turns"]] * np.pi * final_data_matrix[:, self.param["conductor_radius"]] ** 2) < \
                     (self.winding_factor * (window_area - insulation_area)))]

        return data_matrix

    def filter_reluctance_losses(self, data_matrix):
        """
        Filter out design cases based on the calculated hysteresis and DC loss.

        :param data_matrix: Matrix containing the design parameters
        :type data_matrix: ndarray
        """
        # Dictionary to store {initial_permeability: counter}
        mu_r_imag_dict = {}
        counter = 0
        for material_name in self.core_material:
            mu_r_imag_key = material_db.get_material_attribute(material_name=material_name, attribute="initial_permeability")
            mu_r_imag_dict[mu_r_imag_key] = counter
            counter = counter + 1

        material_data_list = [material_db.permeability_data_to_pro_file(material_name=material_name,
                                                                        temperature=self.temperature,
                                                                        frequency=self.frequency,
                                                                        datatype="complex_permeability",
                                                                        datasource=mdb.MaterialDataSource.ManufacturerDatasheet
                                                                        )
                              for material_name in self.core_material]
        # print(material_data_list)

        # Creating interpolation function between mu_imaginary and magnetic flux density
        mu_r_imag_interpol_func = [interp1d(material_data_list[i][0], material_data_list[i][1], kind="cubic") for i in range(len(self.core_material))]

        # Creating mu_imag array corresponding to the material type
        mu_r_imag = np.zeros(len(data_matrix))
        for index in range(len(data_matrix)):
            mu_r_imag[index] = mu_r_imag_interpol_func[mu_r_imag_dict[data_matrix[index, self.param["mu_r_abs"]]]](
                data_matrix[index, self.param["b_max_center"]])

        # Volume chosen as per "Masterthesis_Till_Piepenbrock" pg-45
        volume_center = (np.pi * (data_matrix[:, self.param["core_inner_diameter"]] / 2) ** 2) * \
                        (data_matrix[:, self.param["window_h"]] + data_matrix[:, self.param["core_h_middle"]] - \
                         (data_matrix[:, self.param["n_air_gaps"]] * data_matrix[:, self.param["air_gap_h"]]))
        volume_outer = (np.pi * ((data_matrix[:, self.param["r_outer"]] ** 2) - (data_matrix[:, self.param["r_inner"]] ** 2))) * \
                       (data_matrix[:, self.param["window_h"]] + data_matrix[:, self.param["core_h_middle"]])
        # TODO: Confirm which volume to use

        # volume_center = (np.pi * (data_matrix[:, self.param["core_inner_diameter"]] / 2) ** 2) * \
        #                 (data_matrix[:, self.param["window_h"]])
        # volume_outer = (np.pi * ((data_matrix[:, self.param["r_outer"]] ** 2) -
        #                          (data_matrix[:, self.param["r_inner"]] ** 2))) * \
        #                (data_matrix[:, self.param["window_h"]])

        p_hyst_center = 0.5 * (2 * np.pi * self.frequency) * fmt.mu_0 * mu_r_imag * ((data_matrix[:,
                                                                                      self.param["b_max_center"]] / \
                                                                                      (fmt.mu_0 * data_matrix[:, self.param["mu_r_abs"]])) ** 2)

        p_hyst_outer = 0.5 * (2 * np.pi * self.frequency) * mu_r_imag * fmt.mu_0 * ((data_matrix[:, self.param["b_max_outer"]] / \
                                                                                     (fmt.mu_0 * data_matrix[:, self.param["mu_r_abs"]])) ** 2)

        p_hyst_density_center = p_hyst_center * volume_center
        p_hyst_density_middle = \
            0.5 * (2 * np.pi * self.frequency) * mu_r_imag * fmt.mu_0 * \
            ((data_matrix[:, self.param["total_flux_max"]] / (fmt.mu_0 * data_matrix[:, self.param["mu_r_abs"]])) ** 2) * \
            (1 / (2 * np.pi * data_matrix[:, self.param["core_h_middle"]])) * np.log((data_matrix[:, self.param["r_inner"]] * 2) / \
                                                                                     data_matrix[:, self.param["core_inner_diameter"]])
        p_hyst_density_outer = p_hyst_outer * volume_outer
        total_hyst_loss = p_hyst_density_center + (2 * p_hyst_density_middle) + p_hyst_density_outer

        # Winding loss (only DC loss)
        dc_resistance_wire = (data_matrix[:, self.param["no_of_turns"]] * 2 * np.pi * (data_matrix[:,
                                                                                       self.param["core_inner_diameter"]] / 2 + \
                                                                                       data_matrix[:, self.param["conductor_radius"]])) / \
            (self.copper_conductivity * (np.pi * (data_matrix[:, self.param["conductor_radius"]] ** 2)))

        # I^2 * R loss. Calculation from PEAK current, so division by 2 is needed
        dc_wire_loss = ((self.peak_current ** 2) / 2) * dc_resistance_wire       # Assuming sinusoidal current waveform

        total_hyst_dc_loss = dc_wire_loss + total_hyst_loss
        max_total_loss = max(total_hyst_dc_loss)
        normalized_total_loss = total_hyst_dc_loss / max_total_loss

        core_2daxi_volume = np.pi * (data_matrix[:, self.param["r_outer"]] ** 2) * data_matrix[:, self.param["core_h"]]
        max_core_2daxi_volume = max(core_2daxi_volume)
        normalized_core_2daxi_volume = core_2daxi_volume / max_core_2daxi_volume

        data_matrix = self.add_column_to_data_matrix(data_matrix, total_hyst_loss, 'total_hyst_loss')               # 26
        data_matrix = self.add_column_to_data_matrix(data_matrix, dc_wire_loss, 'dc_wire_loss')                               # 27
        data_matrix = self.add_column_to_data_matrix(data_matrix, total_hyst_dc_loss, 'total_loss')                         # 28
        data_matrix = self.add_column_to_data_matrix(data_matrix, normalized_total_loss, 'normalized_total_loss')   # 29
        data_matrix = self.add_column_to_data_matrix(data_matrix, core_2daxi_volume, 'total_volume')                     # 30
        data_matrix = self.add_column_to_data_matrix(data_matrix, normalized_core_2daxi_volume, 'normalized_total_volume')   # 31

        data_matrix = data_matrix[data_matrix[:, self.param["total_loss"]].argsort()]
        data_matrix = data_matrix[0:int((self.percent_of_total_loss / 100) * len(data_matrix)), :]

        return data_matrix

    def pareto_front_from_data_matrix(self, data_matrix):
        """Get the pareto front from the data matrix."""
        total_loss_vec = data_matrix[:, self.param["total_loss"]]
        core_2daxi_total_volume_vec = data_matrix[:, self.param["total_volume"]]

        tuple_list = np.array(list(zip(total_loss_vec, core_2daxi_total_volume_vec)))

        pareto_tuple_mask_vec = fmt.is_pareto_efficient(tuple_list)

        x_pareto_vec = []
        y_pareto_vec = []

        for count_mask, mask in enumerate(pareto_tuple_mask_vec):
            if mask:
                x_pareto_vec.append(core_2daxi_total_volume_vec[count_mask])
                y_pareto_vec.append(total_loss_vec[count_mask])

        print(f"{len(x_pareto_vec) = }")

        vector_to_sort = np.array([x_pareto_vec, y_pareto_vec])

        # sorting 2d array by 1st row
        # https://stackoverflow.com/questions/49374253/sort-a-numpy-2d-array-by-1st-row-maintaining-columns
        sorted_vector = vector_to_sort[:, vector_to_sort[0].argsort()]
        x_pareto_vec = sorted_vector[0]
        y_pareto_vec = sorted_vector[1]

        return np.array(x_pareto_vec), np.array(y_pareto_vec)

    def filter_reluctance_pareto_front_tolerance(self, data_matrix, factor_min_dc_losses: float = 1):
        """Filter all reluctance calculations to see the pareto front."""
        total_loss_vec = data_matrix[:, self.param["total_loss"]]
        total_volume_vec = data_matrix[:, self.param["total_volume"]]

        min_total_loss = np.min(total_loss_vec)
        print(f"{min_total_loss = }")

        pareto_x_vec, pareto_y_vec = self.pareto_front_from_data_matrix(data_matrix)
        print(f"{pareto_x_vec = }")
        print(f"{pareto_y_vec = }")

        loss_offset = factor_min_dc_losses * min_total_loss

        ref_loss_vec = np.interp(total_volume_vec, pareto_x_vec, pareto_y_vec) + loss_offset

        return_data_matrix = data_matrix[np.where(total_loss_vec < ref_loss_vec)]

        return return_data_matrix

    def fem_simulation(self):
        """Perform FEM simulation of the design cases. Save the result in the given working directory for later analysis."""
        start_time = time.time()

        data_files = []
        file_names = []
        successful_sim_counter = 0
        for count, _ in enumerate(self.data_matrix_fem):

            # MagneticComponent class object
            geo = fmt.MagneticComponent(component_type=self.component_type_dict[self.magnetic_component],
                                        working_directory=self.femmt_working_directory, verbosity=fmt.Verbosity.Silent)

            core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=self.data_matrix_fem[count, self.param["core_inner_diameter"]],
                                                            window_w=self.data_matrix_fem[count, self.param["window_w"]],
                                                            window_h=self.data_matrix_fem[count, self.param["window_h"]],
                                                            core_h=self.data_matrix_fem[count, self.param["core_h"]])

            core = fmt.Core(core_type=fmt.CoreType.Single,
                            core_dimensions=core_dimensions,
                            material=self.core_material_dict[self.data_matrix_fem[count, self.param["mu_r_abs"]]],
                            temperature=self.temperature, frequency=self.frequency,
                            permeability_datasource=fmt.MaterialDataSource.ManufacturerDatasheet,
                            permittivity_datasource=fmt.MaterialDataSource.ManufacturerDatasheet)
            geo.set_core(core)

            # 3. set air gap parameters
            air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
            if int(self.data_matrix_fem[count, self.param["n_air_gaps"]]) == 1:
                air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg,
                                     self.data_matrix_fem[count, self.param["air_gap_h"]],
                                     self.data_matrix_fem[count, self.param["air_gap_position"]])
            else:
                if int(self.data_matrix_fem[count, self.param["mult_air_gap_type"]]) == 1:
                    position_list = list(
                        np.linspace(0, 100, int(self.data_matrix_fem[count, self.param["n_air_gaps"]])))
                    for position in position_list:
                        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg,
                                             self.data_matrix_fem[count, self.param["air_gap_h"]],
                                             position)

                elif int(self.data_matrix_fem[count, self.param["mult_air_gap_type"]]) == 2:
                    position_list = list(
                        np.linspace(0, 100, int(self.data_matrix_fem[count, self.param["n_air_gaps"]]) + 2))
                    position_list.remove(0.0)
                    position_list.remove(100.0)
                    for position in position_list:
                        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg,
                                             self.data_matrix_fem[count, self.param["air_gap_h"]],
                                             position)
            geo.set_air_gaps(air_gaps)

            # 4. set insulations
            insulation = fmt.Insulation()
            insulation.add_core_insulations(self.top_core_insulation, self.bot_core_insulation,
                                            self.left_core_insulation, self.right_core_insulation)
            # insulation.add_winding_insulations([self.inner_winding_insulation], 0.0001)
            insulation.add_winding_insulations([[self.inner_winding_insulation]])
            geo.set_insulation(insulation)

            # 5. create winding window and virtual winding windows (vww)
            winding_window = fmt.WindingWindow(core, insulation)
            vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

            # 6. create conductor and set parameters: use solid wires
            winding = fmt.Conductor(0, fmt.Conductivity.Copper)
            if np.isnan(self.data_matrix_fem[count, self.param["solid_conductor_r"]]):
                winding.set_litz_round_conductor(
                    conductor_radius=None if np.isnan(
                        self.data_matrix_fem[count, self.param["litz_conductor_r"]]) else self.data_matrix_fem[count, self.param["litz_conductor_r"]],
                    number_strands=None if np.isnan(
                        self.data_matrix_fem[count, self.param["litz_strand_n"]]) else self.data_matrix_fem[count, self.param["litz_strand_n"]],
                    strand_radius=None if np.isnan(
                        self.data_matrix_fem[count, self.param["litz_strand_r"]]) else self.data_matrix_fem[count, self.param["litz_strand_r"]],
                    fill_factor=None if np.isnan(
                        self.data_matrix_fem[count, self.param["litz_fill_factor"]]) else self.data_matrix_fem[count, self.param["litz_fill_factor"]],
                    conductor_arrangement=self.winding_scheme_dict[self.winding_scheme][1])
            else:
                winding.set_solid_round_conductor(conductor_radius=self.data_matrix_fem[count, self.param["solid_conductor_r"]],
                                                  conductor_arrangement=self.winding_scheme_dict[self.winding_scheme][1])

            # 7. add conductor to vww and add winding window to MagneticComponent
            vww.set_winding(winding, int(self.data_matrix_fem[count, self.param["no_of_turns"]]), None)
            geo.set_winding_windows([winding_window])

            try:
                # 5. create the model
                geo.create_model(freq=self.frequency, pre_visualize_geometry=False, save_png=False)

                # 6. start simulation
                geo.single_simulation(freq=self.frequency, current=[self.peak_current], show_fem_simulation_results=False)

                source_json_file = os.path.join(self.femmt_working_directory, "results", "log_electro_magnetic.json")
                destination_json_file = os.path.join(self.inductor_fem_simulations_results_directory, f'case_{count}.json')

                shutil.copy(source_json_file, destination_json_file)

                data_files.append(destination_json_file)
                file_names.append(f"case{count}")
                print(f"Case {count} of {len(self.data_matrix_fem)} completed")
                successful_sim_counter = successful_sim_counter + 1

                end_time = time.time()

                time_difference_seconds = end_time - start_time
                time_difference_minutes = time_difference_seconds / 60
                time_difference_hours = time_difference_minutes / 60
                print(f"{time_difference_seconds = }")
                print(f"{time_difference_minutes = }")
                print(f"{time_difference_hours = }")

            except Exception as e:
                print("next iteration")
                # logging.exception(e)
        print(f"Successful FEM simulations: {successful_sim_counter} out of total cases: {len(self.data_matrix_fem)}")

    def add_column_to_data_matrix(self, data_matrix, column_value, column_name: str):
        """
        Add columns to the given matrix.

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

    def write_data_matrix_fem_to_csv(self):
        """Write the data_matrix_fem to csv file for later review of the data."""
        header = list(self.param.keys())
        header.insert(0, 'Case_no.')
        data = self.data_matrix_fem
        a = np.array(range(len(self.data_matrix_fem)))
        data = np.insert(data, 0, a, axis=1)
        file_name = self.working_directory + '/data_matrix_fem.csv'
        with open(file_name, 'w', encoding='UTF8', newline='') as f:
            writer = csv.writer(f)

            # write the header
            writer.writerow(header)

            # write multiple rows
            writer.writerows(data)

    def save_automated_design_settings(self):
        """Create a json file, write all settings that are used to run AutomatedDesign in that particular project."""
        dictionary = {
            "working_directory": self.working_directory,
            "magnetic_component": self.magnetic_component,
            "frequency": self.frequency,
            "temperature": self.temperature,
            "goal_inductance": self.goal_inductance,
            "goal_inductance_percent_tolerance": self.goal_inductance_percent_tolerance,
            "winding_scheme": self.winding_scheme,
            "peak_current": self.peak_current,
            "percent_of_b_sat": self.percent_of_b_sat,
            "percent_of_total_loss": self.percent_of_total_loss,
            "database_core_names": self.database_core_names,
            "database_litz_names": self.database_litz_names,
            "solid_conductor_r": self.solid_conductor_r,
            "manual_core_inner_diameter": self.manual_core_inner_diameter,
            "manual_window_h": self.manual_window_h,
            "manual_window_w": self.manual_window_w,
            "no_of_turns": self.no_of_turns,
            "n_air_gaps": self.n_air_gaps,
            "air_gap_height": self.air_gap_height,
            "air_gap_position": self.air_gap_position,
            "core_material": self.core_material,
            "mult_air_gap_type": self.mult_air_gap_type,
            "top_core_insulation": self.top_core_insulation,
            "bot_core_insulation": self.bot_core_insulation,
            "right_core_insulation": self.right_core_insulation,
            "left_core_insulation": self.left_core_insulation,
            "inner_winding_insulation": self.inner_winding_insulation,
            "manual_litz_conductor_r": self.manual_litz_conductor_r,
            "manual_litz_strand_r": self.manual_litz_strand_r,
            "manual_litz_strand_n": self.manual_litz_strand_n,
            "manual_litz_fill_factor": self.manual_litz_fill_factor
        }

        with open(os.path.join(self.working_directory, "automated_design_settings.json"), "w+", encoding='utf-8') as outfile:
            json.dump(dictionary, outfile, indent=2, ensure_ascii=False)


if __name__ == '__main__':
    task = 'simulate'
    # task = 'load'

    # Input parameters for the Automated Design
    example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
    if not os.path.exists(example_results_folder):
        os.mkdir(example_results_folder)

    # Working directory can be set arbitrarily
    working_directory = os.path.join(example_results_folder, "inductor_optimization")
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    if task == 'simulate':

        core_database = fmt.core_database()
        pq3230 = core_database["PQ 32/30"]
        pq4040 = core_database["PQ 40/40"]
        pq5050 = core_database["PQ 50/50"]
        pq6560 = core_database["PQ 65/60"]

        min_core = pq3230
        max_core = pq6560

        ad = AutomatedDesign(working_directory=working_directory,
                             magnetic_component='inductor',
                             target_inductance=120 * 1e-6,
                             frequency=115000,
                             target_inductance_percent_tolerance=10,
                             winding_scheme='Square',
                             peak_current=8,
                             percent_of_flux_density_saturation=70,
                             percent_of_total_loss=30,
                             database_core_names=[],
                             database_litz_names=['1.5x105x0.1', "1.4x200x0.071", "2.0x405x0.071", "2.0x800x0.05"],
                             solid_conductor_r=[],  # 0.0013
                             manual_core_inner_diameter=list(np.linspace(min_core['core_inner_diameter'], max_core['core_inner_diameter'], 3)),
                             manual_window_h=list(np.linspace(min_core['window_h'], max_core['window_h'], 3)),
                             manual_window_w=list(np.linspace(min_core['window_w'], max_core['window_w'], 3)),
                             no_of_turns=np.arange(1, 80).tolist(),
                             n_air_gaps=[1],
                             air_gap_height=list(np.linspace(0.0001, 0.0005, 20)),
                             air_gap_position=[50],
                             core_material=[fmt.Material.N95],
                             mult_air_gap_type=['center_distributed'],
                             top_core_insulation=0.002,
                             bot_core_insulation=0.002,
                             left_core_insulation=0.0013,
                             right_core_insulation=0.001,
                             inner_winding_insulation=0.0005,
                             temperature=60.0,
                             manual_litz_conductor_r=[],
                             manual_litz_strand_r=[],
                             manual_litz_strand_n=[],
                             manual_litz_fill_factor=[])

        # Create csv file of data_matrix_fem which consist of all the fem simulation cases details
        ad.write_data_matrix_fem_to_csv()

        # Plot of volume vs loss calculated using Reluctance Model
        plot_2d(x_value=ad.data_matrix_fem[:, ad.param["total_volume"]],
                y_value=ad.data_matrix_fem[:, ad.param["total_loss"]],
                x_label='Volume / m\u00b3', y_label='Loss / W', title='Volume vs Loss', plot_color='blue')

        print(f"Total number of cases to be simulated: {len(ad.data_matrix_fem)}")
        print(f"Estimated time of completion of FEM simulations (5 sec/case): {5 * len(ad.data_matrix_fem) / 60 / 60} hours")
        print("##########################################################################################################")
        choice = int(input("Press 1 to run FEM simulations as per the given inputs or any other number to load previous design as per given directory:"))

        if choice == 1:
            # Save simulation settings in json file for later review
            ad.save_automated_design_settings()

            # Run FEM simulation of "self.data_matrix_fem"
            ad.fem_simulation()
    elif task == 'load':

        working_directory = '/home/nikolasf/Dokumente/01_git/30_Python/FEMMT/femmt/examples/example_results/2023-02-28_inductor_optimization_N95_360u_5A'

        # Load design and plot various plots for analysis
        inductance, total_loss, total_volume, total_cost, annotation_list, automated_design_settings = load_fem_simulation_results(
            working_directory=working_directory)

        plot_data = filter_after_fem(inductance=inductance, total_loss=total_loss, total_volume=total_volume, total_cost=total_cost,
                                     annotation_list=annotation_list, goal_inductance=automated_design_settings["goal_inductance"], percent_tolerance=20)

        plot_2d(x_value=plot_data[:, 1], y_value=plot_data[:, 2], z_value=plot_data[:, 3],
                x_label='Volume / m\u00b3', y_label='Loss / W', z_label='Cost / \u20ac', title='Volume vs Loss',
                annotations=plot_data[:, 4], plot_color='RdYlGn_r', inductance_value=plot_data[:, 0])

        plot_2d(x_value=plot_data[:, 1], y_value=plot_data[:, 3], z_value=plot_data[:, 2],
                x_label='Volume / m\u00b3', y_label='Cost / \u20ac', z_label='Loss / W', title='Volume vs Cost',
                annotations=plot_data[:, 4], plot_color='RdYlGn_r', inductance_value=plot_data[:, 0])

    # plot_2d(x_value=data_array[:, 1], y_value=data_array[:, 3], z_value=data_array[:, 2], x_label='Volume / m\u00b3',
    #     y_label='Cost / \u20ac', z_label='Loss / W',
    #         title='Volume vs Cost', plot_color='red')

    # plot_2d(x_value=total_volume, y_value=total_loss, x_label='Volume / m\u00b3', y_label='Loss / W',
    #         title='Volume vs Loss', annotations=annotation_list, plot_color='red')
    # plot_2d(x_value=total_volume, y_value=total_cost, x_label='Volume / m\u00b3', y_label='Cost / \u20ac',
    #         title='Volume vs Cost', annotations=annotation_list, plot_color='red')
    # plot_2d(x_value=total_volume, y_value=inductance, x_label='Volume / m\u00b3', y_label='Inductance / H',
    #         title='Volume vs Inductance', annotations=annotation_list, plot_color='red')
    # plot_3d(x_value=total_volume, y_value=total_loss, z_value=total_cost, x_label='Volume / m\u00b3',
    #         y_label='Loss / W', z_label='Cost / \u20ac', title='Volume vs Loss vs Cost',
    #         annotations=annotation_list, plot_color='red')

    # load_from_single_file(working_directory='D:/Personal_data/MS_Paderborn/Sem4/Project_2/2022-11-27_fem_simulation_data',
    #                file_name='case4.json')
