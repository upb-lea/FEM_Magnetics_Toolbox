"""Thermal simulation functions."""
# Python standard libraries
import re
import json
import numpy as np
import os
from typing import Dict, List

# Third parry libraries
import gmsh
from onelab import onelab

# Local libraries
import femmt.thermal.thermal_functions as thermal_f
from femmt.thermal.thermal_classes import ConstraintPro, FunctionPro, GroupPro, ParametersPro, PostOperationPro
from femmt.data import FileData

def create_case(boundary_regions, boundary_physical_groups, boundary_temperatures, boundary_flags, k_case,
                function_pro: FunctionPro, parameters_pro: ParametersPro, group_pro: GroupPro,
                constraint_pro: ConstraintPro):
    """
    Set boundary conditions and material parameters for the case around the core.

    TODO Set docstring

    """
    group_pro.add_regions(boundary_regions)
    parameters_pro.add_to_parameters(boundary_temperatures)
    parameters_pro.add_to_parameters(boundary_flags)
    constraint_pro.add_boundary_constraint(
        [x for x in zip(boundary_flags.keys(), boundary_regions.keys(), boundary_temperatures.keys())])

    k = {
        "case_top": k_case["top"],
        "case_top_right": k_case["top_right"],
        "case_right": k_case["right"],
        "case_bot_right": k_case["bot_right"],
        "case_bot": k_case["bot"]
    }

    function_pro.add_dicts(k, None)
    group_pro.add_regions({
        "case_top": boundary_physical_groups["top"],
        "case_top_right": boundary_physical_groups["top_right"],
        "case_right": boundary_physical_groups["right"],
        "case_bot_right": boundary_physical_groups["bot_right"],
        "case_bot": boundary_physical_groups["bot"]
    })


def create_insulation(insulation_tag, k_iso, function_pro: FunctionPro, group_pro: GroupPro):
    """Create insulations for thermal simulation."""
    k_iso = {"insulation": k_iso}

    function_pro.add_dicts(k_iso, None)
    group_pro.add_regions({"insulation": insulation_tag})


def create_background(background_tag, k_air, function_pro: FunctionPro, group_pro: GroupPro):
    """Create simulation background."""
    k_air = {"air": k_air}

    function_pro.add_dicts(k_air, None)
    group_pro.add_regions({"air": background_tag})


def create_core_and_air_gaps(core_tags, k_core, core_area, core_parts_losses, air_gaps_tag, k_air_gaps,
                             function_pro: FunctionPro, group_pro: GroupPro):
    """
    Create core and air gaps configurations with associated thermal properties.

    :param core_tags: Tags identifying different parts of the core.
    :param k_core: Thermal conductivity of the core.
    :param core_area: Areas of different core parts.
    :param core_parts_losses: losses in different core parts.
    :param air_gaps_tag: Tag identifying the air gaps.
    :param k_air_gaps: Thermal conductivity of the air gaps.
    :param function_pro: Object to manage function properties.
    :param group_pro: Object to manage group properties.
    """
    # Initialize dictionaries for volumetric heat flux, thermal conductivity, and region
    q_vol = {}  # Dictionary to hold volumetric heat flux for each core part
    k = {}  # Dictionary to hold thermal conductivity for each core part
    regions = {}  # Dictionary to hold regions for each core part
    core_total_str = "{"

    # Convert core tags to a list of lists for consistency
    core_tags = [[tag] for tag in core_tags]
    for index, tag in enumerate(core_tags):
        # Create a unique name for each core part based on its index
        name = f"core_part_{index + 1}"  # Create a name for the core part based on its position
        core_total_str += f"{name}, "

        # Compute the heat flux for the core part using its loss and area
        q_vol[name] = core_parts_losses[index] / core_area[index]  # TODO Core_area should be calculated for every core_part

        # Store the thermal conductivity of the core part
        k[name] = k_core

        # Associate the core part name with its tag
        regions[name] = tag[0]  # as Tag is list of list

    # If air gaps are defined, assign thermal conductivity and tag
    if air_gaps_tag is not None:
        k["air_gaps"] = k_air_gaps
        regions["air_gaps"] = air_gaps_tag

    regions["core_total"] = core_total_str[:-2] + "}"
    # Update the function_pro and group_pro objects with the newly computed dictionaries
    group_pro.add_regions(regions)
    function_pro.add_dicts(k, q_vol)

def create_windings(winding_tags, k_windings, winding_losses, conductor_radii, wire_distances,
                    function_pro: FunctionPro, group_pro: GroupPro):
    """Create windings for the thermal simulation."""
    q_vol = {}
    k = {}
    regions = {}
    windings_total_str = "{"
    tag_counters = {}

    for winding_index, winding in enumerate(winding_tags):
        if winding is not None and len(winding) > 0:
            if winding is not None and len(winding) > 0:
                original_winding = winding[:]

                if len(original_winding) > len(winding_losses[winding_index]):
                    winding_losses[winding_index] = [winding_losses[winding_index][0] / len(original_winding) for _ in range(len(original_winding))]

            for index, tag in enumerate(winding):
                name = f"winding_{winding_index}_{index}"
                windings_total_str += f"{name}, "
                q_vol[name] = thermal_f.calculate_heat_flux_round_wire(winding_losses[winding_index][index],
                                                                       conductor_radii[winding_index],
                                                                       wire_distances[winding_index][index])

                print(q_vol[name])
                k[name] = k_windings
                if tag not in tag_counters:  # The counter is needed here to create a single region for every turn in case of parallel windings
                    tag_counters[tag] = 0
                else:
                    tag_counters[tag] += 1
                regions[name] = tag + tag_counters[tag]

    # Needs to be added. [:-2] removes the last ', '
    regions["windings_total"] = windings_total_str[:-2] + "}"

    function_pro.add_dicts(k, q_vol)
    group_pro.add_regions(regions)


def create_post_operation(thermal_file_path, thermal_influx_file_path, thermal_material_file_path, sensor_points_file,
                          core_file, insulation_file, winding_file, windings, core_parts, print_sensor_values,
                          post_operation_pro: PostOperationPro, flag_insulation):
    """
    Configure post-operation properties for a thermal simulation.

    :param ...: (Include detailed description for each parameter as needed)
    :param post_operation_pro: Object to manage post-operation properties.
    :param flag_insulation: Flag indicating the presence of insulation.
    """
    # Add pos file generation
    post_operation_pro.add_on_elements_of_statement("T", "Total", thermal_file_path)
    post_operation_pro.add_on_elements_of_statement("influx", "Warm", thermal_influx_file_path)
    post_operation_pro.add_on_elements_of_statement("material", "Total", thermal_material_file_path)

    # Add sensor points file generation
    if print_sensor_values:
        post_operation_pro.add_on_point_statement("T", 0.0084, -0.0114, "GmshParsed", sensor_points_file,
                                                  "first_bottom")
        post_operation_pro.add_on_point_statement("T", 0.0084, 0.0002, "GmshParsed", sensor_points_file, "first_middle",
                                                  True)
        post_operation_pro.add_on_point_statement("T", 0.0084, 0.0072, "GmshParsed", sensor_points_file, "first_top",
                                                  True)
        post_operation_pro.add_on_point_statement("T", 0.011, -0.0114, "GmshParsed", sensor_points_file,
                                                  "second_bottom", True)
        post_operation_pro.add_on_point_statement("T", 0.011, 0.0002, "GmshParsed", sensor_points_file, "second_middle",
                                                  True)
        post_operation_pro.add_on_point_statement("T", 0.011, 0.0072, "GmshParsed", sensor_points_file, "second_top",
                                                  True)
        post_operation_pro.add_on_point_statement("T", 0.0132, -0.0089, "GmshParsed", sensor_points_file,
                                                  "third_bottom", True)
        post_operation_pro.add_on_point_statement("T", 0.0036, 0.0011, "GmshParsed", sensor_points_file,
                                                  "air_gap_upper", True)
        post_operation_pro.add_on_point_statement("T", 0.0036, -0.0011, "GmshParsed", sensor_points_file,
                                                  "air_gap_lower", True)

    # Add regions
    # core file will be GmshParsed file as we have many core parts
    # Adding temperature post-operation statements for each core part if they are defined
    core_parts = [[part] for part in core_parts]
    core_append = False  # Initialize the append flag for core parts
    for core_index, core_part in enumerate(core_parts):
        if core_part is not None and len(core_part) > 0:
            # core_part_value = core_part[0]

            name = f"core_part_{core_index + 1}"
            post_operation_pro.add_on_elements_of_statement("T", name, core_file, "GmshParsed", 0, name, core_append)
            if not core_append:
                core_append = True

    # Adding temperature post-operation statement for insulation if insulation is present
    if flag_insulation:
        post_operation_pro.add_on_elements_of_statement("T", "insulation", insulation_file, "SimpleTable", 0)

    # Adding temperature post-operation statements for each winding, considering the case of parallel windings
    append = False
    for winding_index, winding in enumerate(windings):
        if winding is not None and len(winding) > 0:
            # Remove duplicates from region of parallel windings
            for index, _ in enumerate(winding):
                name = f"winding_{winding_index}_{index}"
                post_operation_pro.add_on_elements_of_statement("T", name, winding_file, "GmshParsed", 0, name, append)
                if not append:
                    append = True


def simulate(onelab_folder_path: str, mesh_file: str, solver_file: str, silent: bool) -> None:
    """
    Start the thermal simulation.

    :param onelab_folder_path: onelab folder path
    :type onelab_folder_path: str
    :param mesh_file: file path to mesh
    :type mesh_file: str
    :param solver_file: file path to solver file
    :type solver_file: str
    :param silent: Set to True to prevent terminal outputs
    :type silent: bool
    """
    onelab_client = onelab.client(__file__)

    if silent:
        verbose = "-verbose 1"
    else:
        verbose = "-verbose 5"

    # Run simulations as sub clients (non-blocking??)
    mygetdp = os.path.join(onelab_folder_path, "getdp")
    onelab_client.runSubClient("myGetDP",
                               mygetdp + " " + solver_file + " -msh " + mesh_file + " -solve analysis -v2 " + verbose)


def parse_simple_table(file_path: str) -> np.array:
    """
    Read a thermal simulation output file (e.g. core.txt or insulation.txt) and translates it to a numpy array.

    The output contains temperatures.

    :param file_path: filepath to core.txt or insulation.txt
    :type file_path: str
    :return: Numpy array including temperatures
    :rtype: np.array
    """
    with open(file_path, "r") as fd:
        lines = fd.readlines()
        # print("last line", lines[-1])
        np_array = np.zeros(len(lines))
        for i, line in enumerate(lines):
            np_array[i] = float(line.split(" ")[5])

        return np_array


def parse_gmsh_parsed(file_path: str):
    """Parse gmsh output to a python dict.

    :param file_path: file path
    :type file_path: str
    """
    regex_view_line = "View \"(?P<key>\w+)\" \{\n"
    regex_SP_line = "SP\(-?\d+\.\d+(e-\d+)?,-?\d+\.\d+(e-\d+)?,0\)\{(?P<value>-?\d+\.\d+)\};\n"

    value_dict = {}

    with open(file_path, "r") as fd:
        lines = fd.readlines()
        current_values = []
        current_key = None
        for line in lines:
            if line.startswith("View"):
                if current_values:
                    if len(current_values) == 1:
                        value_dict[current_key] = current_values[0]
                    else:
                        value_dict[current_key] = np.array(current_values)
                    current_values = []
                current_key = re.search(regex_view_line, line).groupdict()["key"]
            elif line.startswith("SP"):
                if current_key is None:
                    raise Exception("Invalid file format: A 'View'-line must be read before a 'SP'-line")
                current_values.append(float(re.search(regex_SP_line, line).groupdict()["value"]))
            elif line.startswith("}"):
                continue
            else:
                raise Exception(f"Unknown line: {line}")

        # Append the last set of current_values to value_dict after the loop ends
        if current_values:
            if len(current_values) == 1:
                value_dict[current_key] = current_values[0]
            else:
                value_dict[current_key] = np.array(current_values)

    return value_dict

def post_operation(case_volume: float, output_file: str, sensor_points_file: str, core_file: str, insulation_file: str,
                   winding_file: str, flag_insulation: bool = True):
    """Post operations after performing the thermal simulation.

    Calculates minimum, maximum and mean temperatures of core,
    conductors and insulations.

    :param case_volume: case volume in m³
    :type case_volume: float
    :param output_file: filepath path to the thermal output log file
    :type output_file: str
    :param sensor_points_file: filepath with all mesh-cell temperatures of the sensors
    :type sensor_points_file: str
    :param core_file: filepath with all mesh-cell temperatures of the core
    :type core_file: str
    :param insulation_file: filepath with all the mesh-cell temperatures of the insulations
    :type insulation_file: str
    :param winding_file: filepath with all the mesh-cell temperatures of the winding
    :type winding_file: str
    """
    # Extract sensor_points
    sensor_point_values = None
    if sensor_points_file is not None:
        sensor_point_values = parse_gmsh_parsed(sensor_points_file)

    # Extract min/max/averages from core, insulations and windings (and air?)
    # core
    core_values = parse_gmsh_parsed(core_file)
    core_parts = {}

    core_part_min = float('inf')
    core_part_max = -float('inf')
    mean_sum = 0

    for core_items, core_value in core_values.items():
        current_min = core_value.min()
        current_max = core_value.max()
        current_mean = core_value.mean()
        core_parts[core_items] = {
            "min": current_min,
            "max": current_max,
            "mean": current_mean
        }
        if current_min < core_part_min:
            core_part_min = current_min

        if current_max > core_part_max:
            core_part_max = current_max

        mean_sum += current_mean
    core_parts["total"] = {
        "min": core_part_min,
        "max": core_part_max,
        "mean": mean_sum / len(core_values.keys())
    }
    print(len(core_values.keys()))
    # windings
    winding_values = parse_gmsh_parsed(winding_file)
    windings = {}

    winding_min = float('inf')
    winding_max = -float('inf')
    mean_sum = 0

    for winding_name, winding_value in winding_values.items():
        current_min = winding_value.min()
        current_max = winding_value.max()
        current_mean = winding_value.mean()

        windings[winding_name] = {
            "min": current_min,
            "max": current_max,
            "mean": current_mean
        }

        if current_min < winding_min:
            winding_min = current_min

        if current_max > winding_max:
            winding_max = current_max

        mean_sum += current_mean

    windings["total"] = {
        "min": winding_min,
        "max": winding_max,
        "mean": mean_sum / len(winding_values.keys())
    }
    print(len(winding_values.keys()))
    misc = {
        "case_volume": case_volume,
        "case_weight": -1,
    }
    # Fill data for json file
    data = {
        "core_parts": core_parts,
        "windings": windings,
        "misc": misc
    }

    if flag_insulation and insulation_file is not None:  # The flag is to check if the thermal simulation with insulation or without
        # insulations
        insulation_values = parse_simple_table(insulation_file)
        insulation_min = insulation_values.min()
        insulation_max = insulation_values.max()
        insulation_mean = insulation_values.mean()

        data["insulations"] = {
            "min": insulation_min,
            "max": insulation_max,
            "mean": insulation_mean
        }

    if sensor_point_values is not None:
        data["sensor_points"] = sensor_point_values

    with open(output_file, "w") as fd:
        json.dump(data, fd, indent=2)

def run_thermal(file_data: FileData, tags_dict: Dict, thermal_conductivity_dict: Dict, boundary_temperatures: Dict,
                boundary_flags: Dict, boundary_physical_groups: Dict, core_area: List, conductor_radii: float,
                wire_distances: float, case_volume: float,
                show_thermal_fem_results: bool, print_sensor_values: bool, silent: bool, flag_insulation: bool = True):
    """
    Run a thermal simulation.

    :param file_data: Contains paths to every folder and file needed in femmt.
    :param tags_dict: Dictionary containing the tags of the case, air, core and windings
    :param boundary_flags: Dictionary containing the boundary flags (0: Neumann / 1: Dirichlet)
    :param boundary_physical_groups:
    :param boundary_temperatures: dict with boundary temperature
    :param thermal_conductivity_dict: Dictionary containing the thermal conductivity material parameters for case, air, core and windings
    :param core_area: Area of the cross-section of the core
    :param conductor_radii: List of the radius for each winding
    :param wire_distances: List of the outer radius for each winding
    :param show_thermal_fem_results: Boolean - Set true when the results shall be shown in a gmsh window
    :type show_thermal_fem_results: bool
    :param print_sensor_values:
    :param silent: True for silent mode (no terminal outputs)
    :param case_volume: volume of the case in m³

    :return: -
    """
    # Get paths
    onelab_folder_path = file_data.onelab_folder_path
    results_folder_path = file_data.results_folder_path
    model_mesh_file_path = file_data.thermal_mesh_file
    results_log_file_path = file_data.e_m_results_log_path

    # Initial Clearing of gmsh data
    gmsh.clear()

    losses = thermal_f.read_results_log(results_log_file_path)

    # Relative paths
    map_pos_file = os.path.join(results_folder_path, "thermal.pos")
    influx_pos_file = os.path.join(results_folder_path, "thermal_influx.pos")
    material_pos_file = os.path.join(results_folder_path, "thermal_material.pos")
    solver_folder_path = os.path.join(os.path.dirname(__file__), "solver")
    thermal_template_file = os.path.join(solver_folder_path, "Thermal.pro")
    parameters_file = os.path.join(solver_folder_path, "Parameters.pro")
    function_file = os.path.join(solver_folder_path, "Function.pro")
    group_file = os.path.join(solver_folder_path, "Group.pro")
    constraint_file = os.path.join(solver_folder_path, "Constraint.pro")
    post_operation_file = os.path.join(solver_folder_path, "PostOperation.pro")
    sensor_points_file = os.path.join(results_folder_path, "sensor_points.txt") if print_sensor_values else None
    core_file = os.path.join(results_folder_path, "core.txt")
    insulation_file = os.path.join(results_folder_path, "insulation.txt") if flag_insulation else None
    winding_file = os.path.join(results_folder_path, "winding.txt")
    output_file = os.path.join(results_folder_path, "results_thermal.json")

    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.open(model_mesh_file_path)

    # Create file wrappers
    parameters_pro = ParametersPro()
    function_pro = FunctionPro()
    group_pro = GroupPro()
    constraint_pro = ConstraintPro()
    post_operation_pro = PostOperationPro()

    # Extract losses
    winding_losses = []

    # Get number of windings
    winding_keys = [key for key in losses.keys() if key.startswith("winding")]
    num_windings = len(winding_keys)
    # for i in range(1, 4):
    for i in range(1, num_windings + 1):
        key = f"winding{i}"
        inner_winding_list = []
        if key in losses:
            for winding in losses[key]["turns"]:
                inner_winding_list.append(winding)
        winding_losses.append(inner_winding_list)

    # Extract core_parts losses as a list
    core_parts_losses = []

    core_part_keys = [key for key in losses.keys() if key.startswith("total_core_part_")]
    for key in core_part_keys:
        core_parts_losses.append(losses[key])

    # TODO All those pro classes could be used as global variables
    create_case(tags_dict["boundary_regions"], boundary_physical_groups, boundary_temperatures, boundary_flags,
                thermal_conductivity_dict["case"], function_pro, parameters_pro, group_pro, constraint_pro)
    create_background(tags_dict["background_tag"], thermal_conductivity_dict["air"], function_pro, group_pro)

    create_core_and_air_gaps(tags_dict["core_tags"], thermal_conductivity_dict["core"], core_area, core_parts_losses,
                             tags_dict["air_gaps_tag"],
                             thermal_conductivity_dict["air_gaps"], function_pro, group_pro)
    create_windings(tags_dict["winding_tags"], thermal_conductivity_dict["winding"], winding_losses, conductor_radii,
                    wire_distances, function_pro, group_pro)

    create_insulation(tags_dict["insulations_tag"], thermal_conductivity_dict["insulation"], function_pro, group_pro) if flag_insulation else None

    create_post_operation(map_pos_file, influx_pos_file, material_pos_file, sensor_points_file, core_file,
                          insulation_file, winding_file, tags_dict["winding_tags"], tags_dict["core_tags"], print_sensor_values,
                          post_operation_pro, flag_insulation)

    # Create files
    parameters_pro.create_file(parameters_file)
    function_pro.create_file(function_file)

    group_pro.create_file(group_file, tags_dict["air_gaps_tag"] is not None, tags_dict["insulations_tag"] is not None)
    constraint_pro.create_file(constraint_file)
    post_operation_pro.create_file(post_operation_file)

    simulate(onelab_folder_path, model_mesh_file_path, thermal_template_file, silent)

    post_operation(case_volume, output_file, sensor_points_file, core_file, insulation_file, winding_file, flag_insulation)

    if show_thermal_fem_results:
        gmsh.open(map_pos_file)
        gmsh.fltk.run()
