import gmsh
from os import path
from onelab import onelab
from thermal.thermal_functions import *
from thermal.thermal_classes import ConstraintPro, FunctionPro, GroupPro, ParametersPro

def get_entities_from_physical_group_list(dim, ids):
    """
    Merges every entity from every physicalGroup id in ids into one list and returns it.
    """
    entity_tags = []
    
    for current_id in ids:
        entity_tags += list(gmsh.model.getEntitiesForPhysicalGroup(dim, current_id))
        
    return entity_tags

def create_physical_group(dim, entities, name):
    tag = gmsh.model.addPhysicalGroup(dim, entities)
    gmsh.model.setPhysicalName(dim, tag, name)

    return tag

def create_case(core_point_tags, k_case, function_pro: FunctionPro, parameters_pro: ParametersPro, group_pro: GroupPro, constraint_pro: ConstraintPro, mesh_size):
    """
    Creates a case around the core and applied the boundary conditions to it.

    core_point_tags: List of corner points of the core. Order: top left, top right, bottom right, bottom left
    
    """
    
    # Tags for the points for each corner from the core
    tl_point = core_point_tags[0] # Top left - default 5
    tr_point = core_point_tags[1] # Top right - default 4
    br_point = core_point_tags[2] # Bottom right - default 3
    bl_point = core_point_tags[3] # Bottom left - default 2


    top_line = 4
    right_line = 3
    bottom_line = 2

    # Get positions from points
    tl_point_pos  = gmsh.model.getValue(0, tl_point, [])
    tr_point_pos  = gmsh.model.getValue(0, tr_point, [])
    br_point_pos  = gmsh.model.getValue(0, br_point, [])
    bl_point_pos  = gmsh.model.getValue(0, bl_point, [])

    # Case size
    case_gap_top = 0.0015
    case_gap_right = 0.0025
    case_gap_bot = 0.002

    # Create 5 new areas: top, top right, right, bottom right, bottom
    # top
    top_case_left_point = gmsh.model.geo.addPoint(tl_point_pos[0], tl_point_pos[1] + case_gap_top, tl_point_pos[2], mesh_size)
    top_case_right_point = gmsh.model.geo.addPoint(tr_point_pos[0], tr_point_pos[1] + case_gap_top, tr_point_pos[2], mesh_size)
    top_case_left_line = gmsh.model.geo.addLine(tl_point, top_case_left_point)
    top_case_top_line = gmsh.model.geo.addLine(top_case_left_point, top_case_right_point)
    top_case_right_line = gmsh.model.geo.addLine(top_case_right_point, tr_point)
    top_case_curve_loop = gmsh.model.geo.addCurveLoop([top_case_left_line, top_case_top_line, top_case_right_line, top_line])
    top_case_surface = gmsh.model.geo.addPlaneSurface([top_case_curve_loop])

    # top right
    top_right_case_top_right_point = gmsh.model.geo.addPoint(tr_point_pos[0] + case_gap_right, tr_point_pos[1] + case_gap_top, tr_point_pos[2], mesh_size)
    top_right_case_right_point = gmsh.model.geo.addPoint(tr_point_pos[0] + case_gap_right, tr_point_pos[1], tr_point_pos[2], mesh_size)
    top_right_case_bottom_line = gmsh.model.geo.addLine(tr_point, top_right_case_right_point)
    top_right_case_right_line = gmsh.model.geo.addLine(top_right_case_right_point, top_right_case_top_right_point)
    top_right_case_top_line = gmsh.model.geo.addLine(top_right_case_top_right_point, top_case_right_point)
    top_right_case_curve_loop = gmsh.model.geo.addCurveLoop([top_case_right_line, top_right_case_bottom_line, top_right_case_right_line, top_right_case_top_line])
    top_right_case_surface = gmsh.model.geo.addPlaneSurface([top_right_case_curve_loop])

    # right
    right_case_bottom_point = gmsh.model.geo.addPoint(br_point_pos[0] + case_gap_right, br_point_pos[1], br_point_pos[2], mesh_size)
    right_case_right_line = gmsh.model.geo.addLine(top_right_case_right_point, right_case_bottom_point)
    right_case_bottom_line = gmsh.model.geo.addLine(right_case_bottom_point, br_point)
    right_case_curve_loop = gmsh.model.geo.addCurveLoop([top_right_case_bottom_line, right_case_right_line, right_case_bottom_line, right_line])
    right_case_surface = gmsh.model.geo.addPlaneSurface([right_case_curve_loop])

    # bottom right
    bottom_right_case_bottom_right_point = gmsh.model.geo.addPoint(br_point_pos[0] + case_gap_right, br_point_pos[1] - case_gap_bot, br_point_pos[2], mesh_size)
    bottom_right_case_bottom_point = gmsh.model.geo.addPoint(br_point_pos[0], br_point_pos[1] - case_gap_bot, br_point_pos[2], mesh_size)
    bottom_right_case_left_line = gmsh.model.geo.addLine(br_point, bottom_right_case_bottom_point)
    bottom_right_case_bottom_line = gmsh.model.geo.addLine(bottom_right_case_bottom_point, bottom_right_case_bottom_right_point)
    bottom_right_case_right_line = gmsh.model.geo.addLine(bottom_right_case_bottom_right_point, right_case_bottom_point)
    bottom_right_case_curve_loop = gmsh.model.geo.addCurveLoop([right_case_bottom_line, bottom_right_case_left_line, bottom_right_case_bottom_line, bottom_right_case_right_line])
    bottom_right_case_surface = gmsh.model.geo.addPlaneSurface([bottom_right_case_curve_loop])

    # bottom
    bottom_case_bottom_left_point = gmsh.model.geo.addPoint(bl_point_pos[0], bl_point_pos[1] - case_gap_bot, bl_point_pos[2], mesh_size)
    bottom_case_bottom_line = gmsh.model.geo.addLine(bottom_right_case_bottom_point, bottom_case_bottom_left_point)
    bottom_case_left_line = gmsh.model.geo.addLine(bottom_case_bottom_left_point, bl_point)
    bottom_case_curve_loop = gmsh.model.geo.addCurveLoop([bottom_case_bottom_line, bottom_case_left_line, bottom_line, bottom_right_case_left_line])
    bottom_case_surface = gmsh.model.geo.addPlaneSurface([bottom_case_curve_loop])

    gmsh.model.geo.synchronize()

    # Create boundary
    boundary_regions = {
        "region_boundary_top"            : top_case_top_line,
        "region_boundary_top_right"      : top_right_case_top_line,
        "region_boundary_right_top"      : top_right_case_right_line,
        "region_boundary_right"          : right_case_right_line,
        "region_boundary_right_bottom"   : bottom_right_case_right_line,
        "region_boundary_bottom_right"   : bottom_right_case_bottom_line,
        "region_boundary_bottom"         : bottom_case_bottom_line
    }

    boundary_temperatures = {
        "value_boundary_top"             : 273,
        "value_boundary_top_right"       : 273,
        "value_boundary_right_top"       : 273,
        "value_boundary_right"           : 273,
        "value_boundary_right_bottom"    : 273,
        "value_boundary_bottom_right"    : 273,
        "value_boundary_bottom"          : 273
    }

    boundary_flags = {
        "flag_boundary_top"             : 1,
        "flag_boundary_top_right"       : 1,
        "flag_boundary_right_top"       : 1,
        "flag_boundary_right"           : 1,
        "flag_boundary_right_bottom"    : 1,
        "flag_boundary_bottom_right"    : 1,
        "flag_boundary_bottom"          : 1
    }

    for key in boundary_regions:
        boundary_regions[key] = create_physical_group(1, [boundary_regions[key]], key)

    group_pro.add_regions(boundary_regions)
    parameters_pro.add_to_parameters(boundary_temperatures)
    parameters_pro.add_to_parameters(boundary_flags)
    constraint_pro.add_boundary_constraint([x for x in zip(boundary_flags.keys(), boundary_regions.keys(), boundary_temperatures.keys())])
    
    # Add surface physical groups
    # INFO: The physical groups are not created in the createRectWithPhysicalGroup because it causes a bug with the index counter when
    # 1D physical groups (lines) are added after 2D physical groups (surfaces)
    top_surface_physical_group = create_physical_group(2, [top_case_surface], "TopCase")
    top_right_surface_physical_group = create_physical_group(2, [top_right_case_surface], "TopRightCase")
    right_surface_physical_group = create_physical_group(2, [right_case_surface], "RightCase")
    bottom_right_surface_physical_group = create_physical_group(2, [bottom_right_case_surface], "BottomRightCase")
    bottom_surface_physical_group = create_physical_group(2, [bottom_case_surface], "BottomCase")

    k = {"case": k_case}
    q_vol = {"case": 0}

    function_pro.add_dicts(k, q_vol)
    all_surfaces_region = f"{{{top_surface_physical_group}, {top_right_surface_physical_group}, {right_surface_physical_group}, {bottom_right_surface_physical_group}, {bottom_surface_physical_group}}}"
    group_pro.add_regions({"case": all_surfaces_region})

    return [[2, top_case_surface], [2, top_right_case_surface], [2, right_case_surface], [2, bottom_right_case_surface], [2, bottom_case_surface]]

def create_background(background_tag, k_air, function_pro: FunctionPro, group_pro: GroupPro):

    k_air = {"air": k_air}
    q_vol_air = {"air": 0}

    # Airgap with material, currently not used
    # TODO Add airgap material when mesh is split
    q_vol_gap = {"air_gap_1": 0}
    k_gap = {"air_gap_1": 0}

    function_pro.add_dicts(k_air, q_vol_air)
    group_pro.add_regions({"air": background_tag})

    return [[2, tag] for tag in gmsh.model.getEntitiesForPhysicalGroup(2, background_tag)]

def create_core(core_tag, k_core, core_area, core_losses, function_pro: FunctionPro, group_pro: GroupPro):
    heat_flux = core_losses/core_area
    k = {"core": k_core}
    q_vol = {"core": heat_flux}

    function_pro.add_dicts(k, q_vol)
    group_pro.add_regions({"core": core_tag})

    return [[2, tag] for tag in gmsh.model.getEntitiesForPhysicalGroup(2, core_tag)]

def create_windings(winding_tags, k_windings, winding_losses, conductor_radii, function_pro: FunctionPro, group_pro: GroupPro):
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
                q_vol[name] = calculate_heat_flux_round_wire(winding_losses[winding_index][index], conductor_radii[winding_index])
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

def thermal(onelab_folder_path, model_mesh_file_path, results_log_file_path, tags_dict, thermal_conductivity_dict, mesh_size, core_area, conductor_radii, show_results, pretty_colors = False, show_before_simulation = False):
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
    :param show_results: Boolean - Set true when the results shall be shown in a gmsh window
    :param pretty_colors: Boolean - Set true if a specified colorization should be applied
    :param show_before_simulation: -  Set true if the mesh should be shown before running the thermal simulation (e.g. to see the colorization)

    :param return: -
    """
    
    losses = read_results_log(results_log_file_path)

    # Realtive paths
    solver_folder = path.join(path.dirname(__file__), "solver")
    map_pos_file = path.join(solver_folder, "map.pos")
    thermal_template_file = path.join(solver_folder, "Thermal.pro")
    parameters_file = path.join(solver_folder, "Parameters.pro")
    function_file = path.join(solver_folder, "Function.pro")
    group_file = path.join(solver_folder, "Group.pro")
    constraint_file = path.join(solver_folder, "Constraint.pro")
    
    gmsh.initialize() # TODO Is this needed?
    gmsh.open(model_mesh_file_path)
    
    # Create file wrappers
    parameters_pro = ParametersPro()
    function_pro = FunctionPro()
    group_pro = GroupPro()
    constraint_pro = ConstraintPro()

    # Extract losses
    winding_losses = []
    for i in range(1, 3):
        key = f"Winding_{i}"
        inner_winding_list = []
        if key in losses:
            for loss in losses[key]["Turns"]:
                inner_winding_list.append(loss)
        winding_losses.append(inner_winding_list)

    core_losses = losses["Core_Eddy_Current"]

    case_dim_tags = create_case(tags_dict["core_point_tags"], thermal_conductivity_dict["case"], function_pro, parameters_pro, group_pro, constraint_pro, mesh_size)
    background_dim_tags = create_background(tags_dict["background_tag"], thermal_conductivity_dict["air"], function_pro, group_pro)
    core_dim_tags = create_core(tags_dict["core_tag"], thermal_conductivity_dict["core"], core_area, core_losses, function_pro, group_pro)
    windings_dim_tags = create_windings(tags_dict["winding_tags"], thermal_conductivity_dict["winding"], winding_losses, conductor_radii, function_pro, group_pro)

    gmsh.model.geo.synchronize()

    # Colorize
    if pretty_colors:
        gmsh.model.setColor(case_dim_tags, 50, 50, 50)
        gmsh.model.setColor(background_dim_tags, 255, 255, 255)
        gmsh.model.setColor(core_dim_tags, 110, 110, 110)
        gmsh.model.setColor(windings_dim_tags, 192, 28, 40)
        
    gmsh.model.mesh.generate()
    gmsh.write(model_mesh_file_path)

    if show_before_simulation:
        gmsh.fltk.run()

    gmsh.finalize() # TODO Is this needed?

    # Create files
    parameters_pro.create_file(parameters_file)
    function_pro.create_file(function_file)
    group_pro.create_file(group_file)
    constraint_pro.create_file(constraint_file)

    simulate(onelab_folder_path, model_mesh_file_path, thermal_template_file)

    if show_results:
        gmsh.initialize()
        gmsh.open(map_pos_file)
        gmsh.fltk.run()
        gmsh.finalize()

def run_with_config(config):
    simulation = config["simulation"]
    tags = config["tags"]
    thermal_conductivity_dict = config["thermal_conductivities"]

    pretty_colors = simulation["pretty_colors"] == "True"
    show_before_simulation = simulation["show_before_simulation"] == "True"
    onelab_folder_path = simulation["onelab_folder_path"]
    results_log_file_path = simulation["results_log_file_path"]
    model_mesh_file_path = simulation["model_mesh_file_path"]
    mesh_size = simulation["mesh_size"]
    core_area = simulation["core_area"]
    conductor_radii = list(simulation["conductor_radii"])
    show_results = simulation["show_results"] == "True"

    thermal(onelab_folder_path, model_mesh_file_path, results_log_file_path, tags, thermal_conductivity_dict, mesh_size, core_area, conductor_radii, show_results, pretty_colors, show_before_simulation)

if __name__ == "__main__":
    read_config = True

    if read_config:
        config_file_path = "thermal_config.json"
        with open(config_file_path, "r") as fd:
            config = json.loads(fd.read())
            run_with_config(config)
    else:
        # Simulation parameters
        pretty_colors = True
        show_before_simulation = False

        # Currently absolute paths
        onelab_folder_path = r"C:\Uni\Bachelorarbeit\onelab-Windows64"
        results_log_file_path = r"C:\Uni\Bachelorarbeit\github\FEM_Magnetics_Toolbox\femmt\results\result_log.json"
        model_mesh_file_path = r"C:\Uni\Bachelorarbeit\github\FEM_Magnetics_Toolbox\femmt\thermal\thermal_mesh.msh"
        winding_tags = [list(range(4000, 4036)), list(range(7000, 7011)), None]
        tags = {
            "core_tag": 2000,
            "background_tag": 1000,
            "winding_tags": winding_tags,
            "core_point_tags": [5, 4, 3, 2] # Order: top left, top right, bottom right, bottom left
        }
        thermal_conductivity_dict = {
            "air": 0.0263,
            "case": 0.3,
            "core": 5,
            "winding": 400 
        }
        mesh_size = 0.001
        core_area = 0.00077
        conductor_radii = [0.0011, 0.0011]
        show_results = True

        thermal(onelab_folder_path, model_mesh_file_path, results_log_file_path, tags, thermal_conductivity_dict, mesh_size, core_area, conductor_radii, show_results, pretty_colors, show_before_simulation)