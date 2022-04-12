import gmsh
from os import path
from onelab import onelab
from .thermal_functions import *
from .thermal_classes import ConstraintPro, FunctionPro, GroupPro, ParametersPro

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
    q_vol = {
        "case_top": 0,
        "case_top_right": 0,
        "case_right": 0,
        "case_bot_right": 0,
        "case_bot": 0 
    }

    function_pro.add_dicts(k, q_vol)
    group_pro.add_regions({
        "case_top": boundary_physical_groups["top"],
        "case_top_right": boundary_physical_groups["top_right"],
        "case_right": boundary_physical_groups["right"],
        "case_bot_right": boundary_physical_groups["bot_right"],
        "case_bot": boundary_physical_groups["bot"] 
    })

    dim_tags = []
    for key in boundary_physical_groups:
        for tag in gmsh.model.getEntitiesForPhysicalGroup(2, boundary_physical_groups[key]):
            dim_tags.append([2, tag])

    return dim_tags

def create_isolation(isolation_tag, k_iso, function_pro: FunctionPro, group_pro: GroupPro):
    k_iso = {"isolation": k_iso}
    q_vol_iso = {"isolation": 0}

    function_pro.add_dicts(k_iso, q_vol_iso)
    group_pro.add_regions({"isolation": isolation_tag})

    return [[2, tag] for tag in gmsh.model.getEntitiesForPhysicalGroup(2, isolation_tag)]

def create_background(background_tag, k_air, function_pro: FunctionPro, group_pro: GroupPro):
    k_air = {"air": k_air}
    q_vol_air = {"air": 0}

    function_pro.add_dicts(k_air, q_vol_air)
    group_pro.add_regions({"air": background_tag})
    
    return [[2, tag] for tag in gmsh.model.getEntitiesForPhysicalGroup(2, background_tag)]

def create_core_and_air_gaps(core_tag, k_core, core_area, core_losses, air_gaps_tag, k_air_gaps, function_pro: FunctionPro, group_pro: GroupPro):
    heat_flux = core_losses/core_area
    
    if air_gaps_tag is not None:
        k = {
            "core": k_core,
            "air_gaps": k_air_gaps 
        }
        q_vol = {
            "core": heat_flux,
            "air_gaps": 0
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
        
    

    return [[2, tag] for tag in gmsh.model.getEntitiesForPhysicalGroup(2, core_tag)], [[2, tag] for tag in gmsh.model.getEntitiesForPhysicalGroup(2, air_gaps_tag)] if air_gaps_tag is not None else None

def create_windings(winding_tags, k_windings, winding_losses, conductor_radii, wire_distances, function_pro: FunctionPro, group_pro: GroupPro):
    q_vol = {}
    k = {}
    regions = {}
    windings_total_str = "{"
    entities = []

    for winding_index, winding in enumerate(winding_tags):
        if winding is not None and len(winding) > 0:
            for index, tag in enumerate(winding):
                name = f"winding_{winding_index}_{index}"
                windings_total_str += f"{name}, "
                q_vol[name] = calculate_heat_flux_round_wire(winding_losses[winding_index][index], conductor_radii[winding_index], wire_distances[winding_index][index])
                k[name] = k_windings
                regions[name] = tag
                for entity in gmsh.model.getEntitiesForPhysicalGroup(2, tag):
                    entities.append(entity)
                
    # Needs to be added. [:-2] removes the last ', '
    regions["windings_total"] = windings_total_str[:-2] + "}"

    function_pro.add_dicts(k, q_vol)
    group_pro.add_regions(regions)

    return [[2, tag] for tag in entities]
    
def simulate(onelab_folder_path, mesh_file, solver_file):
    c = onelab.client(__file__)

    # Run simulations as sub clients (non blocking??)
    mygetdp = path.join(onelab_folder_path, "getdp")
    c.runSubClient("myGetDP", mygetdp + " " + solver_file + " -msh " + mesh_file + " -solve analysis -v2")

def run_thermal(onelab_folder_path, results_folder_path, model_mesh_file_path, results_log_file_path, 
    tags_dict, thermal_conductivity_dict, boundary_temperatures, 
    boundary_flags, boundary_physical_groups, core_area, conductor_radii, wire_distances,
    show_results, pretty_colors = False, show_before_simulation = False):
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
    :param pretty_colors: Boolean - Set true if a specified colorization should be applied
    :param show_before_simulation: -  Set true if the mesh should be shown before running the thermal simulation (e.g. to see the colorization)

    :param return: -
    """

    # Initial Clearing of gmsh data
    gmsh.clear()
    
    losses = read_results_log(results_log_file_path)

    # Relative paths
    map_pos_file = path.join(results_folder_path, "thermal.pos")
    influx_pos_file = path.join(results_folder_path, "thermal_influx.pos")
    solver_folder_path = path.join(os.path.dirname(__file__), "solver")
    thermal_template_file = path.join(solver_folder_path, "Thermal.pro")
    parameters_file = path.join(solver_folder_path, "Parameters.pro")
    function_file = path.join(solver_folder_path, "Function.pro")
    group_file = path.join(solver_folder_path, "Group.pro")
    constraint_file = path.join(solver_folder_path, "Constraint.pro")

    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.open(model_mesh_file_path)
    
    # Create file wrappers
    parameters_pro = ParametersPro()
    function_pro = FunctionPro()
    group_pro = GroupPro()
    constraint_pro = ConstraintPro()

    parameters_pro.add_to_parameters({
        "thermal_file": map_pos_file.replace("\\", "/"),
        "thermal_influx_file": influx_pos_file.replace("\\", "/")
    })

    # Extract losses
    winding_losses = []
    for i in range(1, 3):
        key = f"Winding_{i}"
        inner_winding_list = []
        if key in losses:
            for loss in losses[key]["Turns"]:
                inner_winding_list.append(loss)
        winding_losses.append(inner_winding_list)

    core_losses = losses["Core_Eddy_Current"] + losses["Core_Hysteresis"]

    # TODO All those pro classes could be used as global variables
    case_dim_tags = create_case(tags_dict["boundary_regions"], boundary_physical_groups, boundary_temperatures, boundary_flags, thermal_conductivity_dict["case"], function_pro, parameters_pro, group_pro, constraint_pro)
    background_dim_tags = create_background(tags_dict["background_tag"], thermal_conductivity_dict["air"], function_pro, group_pro)
    core_dim_tags, air_gaps_dim_tags = create_core_and_air_gaps(tags_dict["core_tag"], thermal_conductivity_dict["core"], core_area, core_losses, tags_dict["air_gaps_tag"], thermal_conductivity_dict["air_gaps"], function_pro, group_pro)
    windings_dim_tags = create_windings(tags_dict["winding_tags"], thermal_conductivity_dict["winding"], winding_losses, conductor_radii, wire_distances, function_pro, group_pro)
    isolation_dim_tags = create_isolation(tags_dict["isolations_tag"], thermal_conductivity_dict["isolation"], function_pro, group_pro)

    gmsh.model.geo.synchronize()

    # Colorize
    if pretty_colors:
        gmsh.model.setColor(case_dim_tags, 50, 50, 50)
        gmsh.model.setColor(background_dim_tags, 255, 255, 255)
        gmsh.model.setColor(core_dim_tags, 110, 110, 110)
        gmsh.model.setColor(windings_dim_tags, 192, 28, 40)
        if air_gaps_dim_tags is not None:
            gmsh.model.setColor(air_gaps_dim_tags, 203, 156, 190)
        gmsh.model.setColor(isolation_dim_tags, 59, 59, 59)
        
    gmsh.model.mesh.generate()
    gmsh.write(model_mesh_file_path)

    if show_before_simulation:
        gmsh.fltk.run()

    # Create files
    parameters_pro.create_file(parameters_file)
    function_pro.create_file(function_file)
    group_pro.create_file(group_file, tags_dict["air_gaps_tag"] != None)
    constraint_pro.create_file(constraint_file)

    simulate(onelab_folder_path, model_mesh_file_path, thermal_template_file)

    if show_results:
        gmsh.open(map_pos_file)
        gmsh.fltk.run()