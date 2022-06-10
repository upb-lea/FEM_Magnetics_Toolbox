import gmsh
import re
import json
from os import path
from onelab import onelab
from .thermal_functions import *
from .thermal_classes import ConstraintPro, FunctionPro, GroupPro, ParametersPro, PostOperationPro
from ..femmt_functions import *

def create_case(boundary_regions, boundary_physical_groups, boundary_temperatures, boundary_flags, k_case, function_pro: FunctionPro, parameters_pro: ParametersPro, group_pro: GroupPro, constraint_pro: ConstraintPro):
    """
    Sets boundary conditions and material parameters for the case around the core.

    TODO Set docstring
    
    """
    group_pro.add_regions(boundary_regions)
    parameters_pro.add_to_parameters(boundary_temperatures)
    parameters_pro.add_to_parameters(boundary_flags)
    constraint_pro.add_boundary_constraint([x for x in zip(boundary_flags.keys(), boundary_regions.keys(), boundary_temperatures.keys())])
    
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

def create_isolation(isolation_tag, k_iso, function_pro: FunctionPro, group_pro: GroupPro):
    k_iso = {"isolation": k_iso}

    function_pro.add_dicts(k_iso, None)
    group_pro.add_regions({"isolation": isolation_tag})

def create_background(background_tag, k_air, function_pro: FunctionPro, group_pro: GroupPro):
    k_air = {"air": k_air}

    function_pro.add_dicts(k_air, None)
    group_pro.add_regions({"air": background_tag})

def create_core_and_air_gaps(core_tag, k_core, core_area, core_losses, air_gaps_tag, k_air_gaps, function_pro: FunctionPro, group_pro: GroupPro):
    heat_flux = core_losses/core_area
    print(heat_flux)
    if air_gaps_tag is not None:
        k = {
            "core": k_core,
            "air_gaps": k_air_gaps 
        }
        q_vol = {
            "core": heat_flux
        }
        group_pro.add_regions({
                "core": core_tag,
                "air_gaps": air_gaps_tag
        })
        function_pro.add_dicts(k, q_vol)
    else:
        k = {"core": k_core}
        q_vol = {"core": heat_flux}
        group_pro.add_regions({"core": core_tag})
        function_pro.add_dicts(k, q_vol)

def create_windings(winding_tags, k_windings, winding_losses, conductor_radii, wire_distances, function_pro: FunctionPro, group_pro: GroupPro):
    q_vol = {}
    k = {}
    regions = {}
    windings_total_str = "{"

    for winding_index, winding in enumerate(winding_tags):
        if winding is not None and len(winding) > 0:
            for index, tag in enumerate(winding):
                name = f"winding_{winding_index}_{index}"
                windings_total_str += f"{name}, "
                q_vol[name] = calculate_heat_flux_round_wire(winding_losses[winding_index][index], conductor_radii[winding_index], wire_distances[winding_index][index])
                print(q_vol[name])
                k[name] = k_windings
                regions[name] = tag
                
    # Needs to be added. [:-2] removes the last ', '
    regions["windings_total"] = windings_total_str[:-2] + "}"

    function_pro.add_dicts(k, q_vol)
    group_pro.add_regions(regions)

def create_post_operation(thermal_file_path, thermal_influx_file_path, thermal_material_file_path, sensor_points_file, core_file, isolation_file, winding_file, windings, post_operation_pro: PostOperationPro):

    # Add pos file generation
    post_operation_pro.add_on_elements_of_statement("T", "Total", thermal_file_path)
    post_operation_pro.add_on_elements_of_statement("influx", "Warm", thermal_influx_file_path)
    post_operation_pro.add_on_elements_of_statement("material", "Total", thermal_material_file_path)

    # Add sensor points file generation
    post_operation_pro.add_on_point_statement("T", 0.0084, -0.0114, "GmshParsed", sensor_points_file, "first_bottom")
    post_operation_pro.add_on_point_statement("T", 0.0084, 0.0002, "GmshParsed", sensor_points_file, "first_middle", True)
    post_operation_pro.add_on_point_statement("T", 0.0084, 0.0072, "GmshParsed", sensor_points_file, "first_top", True)
    post_operation_pro.add_on_point_statement("T", 0.011, -0.0114, "GmshParsed", sensor_points_file, "second_bottom", True)
    post_operation_pro.add_on_point_statement("T", 0.011, 0.0002, "GmshParsed", sensor_points_file, "second_middle", True)
    post_operation_pro.add_on_point_statement("T", 0.011, 0.0072, "GmshParsed", sensor_points_file, "second_top", True)
    post_operation_pro.add_on_point_statement("T", 0.0132, -0.0089, "GmshParsed", sensor_points_file, "third_bottom", True)
    post_operation_pro.add_on_point_statement("T", 0.0036, 0.0011, "GmshParsed", sensor_points_file, "air_gap_upper", True)
    post_operation_pro.add_on_point_statement("T", 0.0036, -0.0011, "GmshParsed", sensor_points_file, "air_gap_lower", True)

    # Add regions
    post_operation_pro.add_on_elements_of_statement("T", "core", core_file, "SimpleTable", 0)
    post_operation_pro.add_on_elements_of_statement("T", "isolation", isolation_file, "SimpleTable", 0)

    append = False
    for winding_index, winding in enumerate(windings):
        if winding is not None and len(winding) > 0:
            for index, tag in enumerate(winding):
                name = f"winding_{winding_index}_{index}"
                post_operation_pro.add_on_elements_of_statement("T", name, winding_file, "GmshParsed", 0, name, append)
                if not append:
                    append = True

def simulate(onelab_folder_path, mesh_file, solver_file):
    c = onelab.client(__file__)

    # Run simulations as sub clients (non blocking??)
    mygetdp = path.join(onelab_folder_path, "getdp")
    c.runSubClient("myGetDP", mygetdp + " " + solver_file + " -msh " + mesh_file + " -solve analysis -v2")

def parse_simple_table(file_path):
    with open(file_path, "r") as fd:
        lines = fd.readlines()
        print("lastline", lines[-1])
        np_array = np.zeros(len(lines))
        for i, line in enumerate(lines):
            np_array[i] = float(line.split(" ")[5])

        return np_array

def parse_gmsh_parsed(file_path):
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

    return value_dict

def post_operation(output_file, sensor_points_file, core_file, isolation_file, winding_file):
    # Extract sensor_points
    sensor_point_values = parse_gmsh_parsed(sensor_points_file)

    # Extract min/max/averages from core, isolations and windings (and air?)
    # core
    core_values = parse_simple_table(core_file)
    core_min = core_values.min()
    core_max = core_values.max()
    core_mean = core_values.mean()

    # isolations
    isolation_values = parse_simple_table(isolation_file)
    isolation_min = isolation_values.min()
    isolation_max = isolation_values.max()
    isolation_mean = isolation_values.mean()

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
        "mean": mean_sum/len(winding_values.keys())
    }

    # fill data for json file
    data = {
        "sensor_points": sensor_point_values,
        "core": {
            "min": core_min,
            "max": core_max,
            "mean": core_mean
        },
        "isolations": {
            "min": isolation_min,
            "max": isolation_max,
            "mean": isolation_mean
        },
        "windings": windings,
    }

    with open(output_file, "w") as fd:
        json.dump(data, fd, indent=2)

def run_thermal(onelab_folder_path, results_folder_path, model_mesh_file_path, results_log_file_path, 
    tags_dict, thermal_conductivity_dict, boundary_temperatures, 
    boundary_flags, boundary_physical_groups, core_area, conductor_radii, wire_distances,
    show_results: bool):
    """
    Runs a thermal simulation.
    
    :param onelab_folder_path: Path to the onelab directory
    :param model_mesh_file_path: Path to the .msh file generated by the magnetic simulation
    :oaran results_log_file_path: Path to the results log file generated by the magnetic simulation
    :param tags_dict: Dictionary containing the tags of the case, air, core and windings
    :param thermal_conductivity_dict: Dictionary containing the thermal conductivity material parameters for case, air, core and windings
    :param mesh_size: Settings the mesh size for the case wihich will be constructed around the core
    :param core_area: Area of the cross-section of the core
    :param conductor_radii: List of the radius for each winding 
    :param wire_distances: List of the outer radius for each winding
    :param show_results: Boolean - Set true when the results shall be shown in a gmsh window
    :type show_results: bool

    :param return: -
    """

    # Initial Clearing of gmsh data
    gmsh.clear()
    
    losses = read_results_log(results_log_file_path)

    # Relative paths
    map_pos_file = path.join(results_folder_path, "thermal.pos")
    influx_pos_file = path.join(results_folder_path, "thermal_influx.pos")
    material_pos_file = path.join(results_folder_path, "thermal_material.pos")
    solver_folder_path = path.join(os.path.dirname(__file__), "solver")
    thermal_template_file = path.join(solver_folder_path, "Thermal.pro")
    parameters_file = path.join(solver_folder_path, "Parameters.pro")
    function_file = path.join(solver_folder_path, "Function.pro")
    group_file = path.join(solver_folder_path, "Group.pro")
    constraint_file = path.join(solver_folder_path, "Constraint.pro")
    post_operation_file = path.join(solver_folder_path, "PostOperation.pro")
    sensor_points_file = path.join(results_folder_path, "sensor_points.txt")
    core_file = path.join(results_folder_path, "core.txt")
    isolation_file = path.join(results_folder_path, "isolation.txt")
    winding_file = path.join(results_folder_path, "winding.txt")
    output_file = path.join(results_folder_path, "results_thermal.json")

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
    for i in range(1, 3):
        key = f"winding{i}"
        inner_winding_list = []
        if key in losses:
            for winding in losses[key]["turns"]:
                inner_winding_list.append(winding)
        winding_losses.append(inner_winding_list)

    core_losses = losses["core"]

    # TODO All those pro classes could be used as global variables
    create_case(tags_dict["boundary_regions"], boundary_physical_groups, boundary_temperatures, boundary_flags,
         thermal_conductivity_dict["case"], function_pro, parameters_pro, group_pro, constraint_pro)
    create_background(tags_dict["background_tag"], thermal_conductivity_dict["air"], function_pro, group_pro)
    create_core_and_air_gaps(tags_dict["core_tag"], thermal_conductivity_dict["core"], core_area, core_losses, tags_dict["air_gaps_tag"], 
        thermal_conductivity_dict["air_gaps"], function_pro, group_pro)
    create_windings(tags_dict["winding_tags"], thermal_conductivity_dict["winding"], winding_losses, conductor_radii, wire_distances, function_pro, group_pro)
    create_isolation(tags_dict["isolations_tag"], thermal_conductivity_dict["isolation"], function_pro, group_pro)
    #create_post_operation(map_pos_file.replace("\\", "/"), influx_pos_file.replace("\\", "/"), material_pos_file.replace("\\", "/"), sensor_points_file, core_file,
    #    isolation_file, winding_file, tags_dict["winding_tags"], post_operation_pro)
    create_post_operation(map_pos_file, influx_pos_file, material_pos_file, sensor_points_file, core_file,
        isolation_file, winding_file, tags_dict["winding_tags"], post_operation_pro)

    # Create files
    parameters_pro.create_file(parameters_file)
    function_pro.create_file(function_file)
    group_pro.create_file(group_file, tags_dict["air_gaps_tag"] != None)
    constraint_pro.create_file(constraint_file)
    post_operation_pro.create_file(post_operation_file)

    simulate(onelab_folder_path, model_mesh_file_path, thermal_template_file)

    post_operation(output_file, sensor_points_file, core_file, isolation_file, winding_file)

    if show_results:
        gmsh.open(map_pos_file)
        gmsh.fltk.run()