# python libraries
from typing import List, Dict

# 3rd party libraries
from matplotlib import pyplot as plt
import numpy as np

# femmt libraries
import femmt.functions as ff
from femmt.optimization.ito_dtos import *

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
        l_value_str = [str(round(inductance, 6)) for inductance in inductance_value]

    x_value_str = [str(round(x, 6)) for x in x_value]
    y_value_str = [str(round(y, 3)) for y in y_value]

    fig, ax = plt.subplots()
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    # c = np.random.randint(1, 5, size=len(y_value))
    # norm = plt.Normalize(1, 4)
    # cmap = plt.cm.RdYlGn

    # sc = plt.scatter(x_value, y_value, c=c, s=50, cmap=cmap, norm=norm)
    if z_value is None:
        sc = plt.scatter(x_value, y_value, c='#%02x%02x%02x' % ff.colors_femmt_default[plot_color])
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
        """Create popover annotations in 2d plot"""

        pos = sc.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        text = ""
        if z_label is None and inductance_value is None:
            text = "case: {}\n{}: {}\n{}: {}".\
                format(" ".join([names[n] for n in ind["ind"]]),
                       x_label, " ".join([x_value_str[n] for n in ind["ind"]]),
                       y_label, " ".join([y_value_str[n] for n in ind["ind"]]))
        elif z_label is not None and inductance_value is None:
            text = "case: {}\n{}: {}\n{}: {}\n{}: {}". \
                format(" ".join([names[n] for n in ind["ind"]]),
                       x_label, " ".join([x_value_str[n] for n in ind["ind"]]),
                       y_label, " ".join([y_value_str[n] for n in ind["ind"]]),
                       z_label, " ".join([z_value_str[n] for n in ind["ind"]]))
        elif z_label is None and inductance_value is not None:
            text = "case: {}\n{}: {}\n{}: {}\n{}: {}". \
                format(" ".join([names[n] for n in ind["ind"]]),
                       x_label, " ".join([x_value_str[n] for n in ind["ind"]]),
                       y_label, " ".join([y_value_str[n] for n in ind["ind"]]),
                       l_label, " ".join([l_value_str[n] for n in ind["ind"]]))
        else:
            text = "case: {}\n{}: {}\n{}: {}\n{}: {}\n{}: {}".\
                format(" ".join([names[n] for n in ind["ind"]]),
                       x_label, " ".join([x_value_str[n] for n in ind["ind"]]),
                       y_label, " ".join([y_value_str[n] for n in ind["ind"]]),
                       z_label, " ".join([z_value_str[n] for n in ind["ind"]]),
                       l_label, " ".join([l_value_str[n] for n in ind["ind"]]))
        annot.set_text(text)
        # annot.get_bbox_patch().set_facecolor(cmap(norm(c[ind["ind"][0]])))
        annot.get_bbox_patch().set_alpha(0.8)

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

def pareto_front_from_dtos(dto_list: List[ItoSingleResultFile]) -> tuple:
    """
    Calculates the Pareto front from a list of ItoSingleResultFiles.

    :param dto_list: List of ItoSingleResultFiles
    :type dto_list: List[ItoSingleResultFiles]
    :return: x-Pareto vector, y-Pareto vector
    :rtype: tuple
    """
    x_vec = np.array([])
    y_vec = np.array([])
    tuple_vec = []

    for dto in dto_list:
        x_vec = np.append(x_vec, dto.core_2daxi_total_volume)
        y_vec = np.append(y_vec, dto.total_loss)
        tuple_vec.append((dto.core_2daxi_total_volume, dto.total_loss))

    tuple_vec = np.array(tuple_vec)

    pareto_tuple_mask_vec = is_pareto_efficient(tuple_vec)

    x_pareto_vec = []
    y_pareto_vec = []

    for count_mask, mask in enumerate(pareto_tuple_mask_vec):
        if mask:
            x_pareto_vec.append(x_vec[count_mask])
            y_pareto_vec.append(y_vec[count_mask])

    print(f"{len(x_pareto_vec) = }")

    return np.array(x_pareto_vec), np.array(y_pareto_vec)


def pareto_front_from_result_dicts(result_dict_list: List[Dict]) -> tuple:
    """
    Calculates the Pareto front from a list of result log dictionaries.

    :param result_dict_list: List of result log dictionaries
    :type result_dict_list: List[Dict]
    :return: x-Pareto vector, y-Pareto vector
    :rtype: tuple
    """
    x_vec = np.array([])
    y_vec = np.array([])
    tuple_vec = []

    for result_dict in result_dict_list:
        x_vec = np.append(x_vec, result_dict["misc"]["core_2daxi_total_volume"])
        y_vec = np.append(y_vec, result_dict["total_losses"]["all_windings"] + result_dict["total_losses"]["core"])
        tuple_vec.append((result_dict["misc"]["core_2daxi_total_volume"], result_dict["total_losses"]["all_windings"] + result_dict["total_losses"]["core"]))

    tuple_vec = np.array(tuple_vec)

    pareto_tuple_mask_vec = is_pareto_efficient(tuple_vec)

    x_pareto_vec = []
    y_pareto_vec = []

    for count_mask, mask in enumerate(pareto_tuple_mask_vec):
        if mask:
            x_pareto_vec.append(x_vec[count_mask])
            y_pareto_vec.append(y_vec[count_mask])

    print(f"{len(x_pareto_vec) = }")

    return np.array(x_pareto_vec), np.array(y_pareto_vec)
