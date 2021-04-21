from functions import inner_points, min_max_inner_points, call_for_path
from onelab import onelab
import gmsh, pathlib, sys, os
import FEMMT as geo
import numpy as np


# Initialization
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.model.add("geometry")


# ------------------------------------------ Geometry -------------------------------------------
# Core generation
if geo.core_type == "EI":
    # --------------------------------------- Points --------------------------------------------
    if geo.y_symmetric == 1:
        if geo.axi_symmetric == 1:

            # Find points of air gaps (used later)
            if geo.n_air_gaps > 0:
                center_right = min_max_inner_points(geo.p_window[4], geo.p_window[6], geo.p_air_gaps)
                island_right = inner_points(geo.p_window[4], geo.p_window[6], geo.p_air_gaps)

            # Pre-Definitions
            # Points
            p_core = []
            p_island = []
            p_cond = []
            # Curves
            l_bound_core = []
            l_bound_air = []
            l_core_air = []
            l_cond = []
            curve_loop_cond = []
            # Curve Loops
            curve_loop_island = []
            curve_loop_air = []
            curve_loop_bound = []
            # Plane Surfaces
            plane_surface_core = []
            plane_surface_cond = []
            plane_surface_air = []

            # =====================
            # Main Core
            # Points of Main Core (index refers to sketch)
            if geo.n_air_gaps > 0:
                p_core.append(gmsh.model.geo.addPoint(0, center_right[0][1], center_right[0][2], center_right[0][3]))
            if geo.n_air_gaps == 0:
                p_core.append(None)  # dummy filled for no air gap special case
            p_core.append(gmsh.model.geo.addPoint(0, geo.p_outer[1][1], geo.p_outer[1][2], geo.p_outer[1][3]))
            p_core.append(gmsh.model.geo.addPoint(geo.p_outer[1][0], geo.p_outer[1][1], geo.p_outer[1][2], geo.p_outer[1][3]))
            p_core.append(gmsh.model.geo.addPoint(geo.p_outer[3][0], geo.p_outer[3][1], geo.p_outer[3][2], geo.p_outer[3][3]))
            p_core.append(gmsh.model.geo.addPoint(0, geo.p_outer[3][1], geo.p_outer[3][2], geo.p_outer[3][3]))
            if geo.n_air_gaps > 0:
                p_core.append(gmsh.model.geo.addPoint(0, center_right[1][1], center_right[1][2], center_right[1][3]))
                p_core.append(gmsh.model.geo.addPoint(center_right[1][0], center_right[1][1], center_right[1][2], center_right[1][3]))
            if geo.n_air_gaps == 0:
                p_core.append(None)  # dummy filled for no air gap special case
                p_core.append(None)  # dummy filled for no air gap special case
            p_core.append(gmsh.model.geo.addPoint(geo.p_window[6][0], geo.p_window[6][1], geo.p_window[6][2], geo.p_window[6][3]))
            p_core.append(gmsh.model.geo.addPoint(geo.p_window[7][0], geo.p_window[7][1], geo.p_window[7][2], geo.p_window[7][3]))
            p_core.append(gmsh.model.geo.addPoint(geo.p_window[5][0], geo.p_window[5][1], geo.p_window[5][2], geo.p_window[5][3]))
            p_core.append(gmsh.model.geo.addPoint(geo.p_window[4][0], geo.p_window[4][1], geo.p_window[4][2], geo.p_window[4][3]))
            if geo.n_air_gaps > 0:
                p_core.append(gmsh.model.geo.addPoint(center_right[0][0], center_right[0][1], center_right[0][2], center_right[0][3]))
            if geo.n_air_gaps == 0:
                p_core.append(None)  # dummy filled for no air gap special case
            # Curves of Main Core (index refers to sketch)
# To be added: Case Air Gaps directly on outer leg
            # Curves: Boundary - Core
            if geo.n_air_gaps > 0:
                l_bound_core.append(gmsh.model.geo.addLine(p_core[0], p_core[1]))
            l_bound_core.append(gmsh.model.geo.addLine(p_core[1], p_core[2]))
            l_bound_core.append(gmsh.model.geo.addLine(p_core[2], p_core[3]))
            l_bound_core.append(gmsh.model.geo.addLine(p_core[3], p_core[4]))
            if geo.n_air_gaps > 0:
                l_bound_core.append(gmsh.model.geo.addLine(p_core[4], p_core[5]))
            if geo.n_air_gaps == 0:
                l_bound_core.append(gmsh.model.geo.addLine(p_core[4], p_core[1]))
            # Curves: Core - Air
            if geo.n_air_gaps > 0:
                l_core_air.append(gmsh.model.geo.addLine(p_core[5], p_core[6]))
                l_core_air.append(gmsh.model.geo.addLine(p_core[6], p_core[7]))
            l_core_air.append(gmsh.model.geo.addLine(p_core[7], p_core[8]))
            l_core_air.append(gmsh.model.geo.addLine(p_core[8], p_core[9]))
            l_core_air.append(gmsh.model.geo.addLine(p_core[9], p_core[10]))
            if geo.n_air_gaps > 0:
                l_core_air.append(gmsh.model.geo.addLine(p_core[10], p_core[11]))
                l_core_air.append(gmsh.model.geo.addLine(p_core[11], p_core[0]))
            if geo.n_air_gaps == 0:
                l_core_air.append(gmsh.model.geo.addLine(p_core[10], p_core[7]))
            # Plane: Main Core --> plane_surface_core[0]
            if geo.n_air_gaps > 0:
                curve_loop_core = gmsh.model.geo.addCurveLoop(l_bound_core + l_core_air)
                plane_surface_core.append(gmsh.model.geo.addPlaneSurface([curve_loop_core]))
            if geo.n_air_gaps == 0:
                curve_loop_bound_core = gmsh.model.geo.addCurveLoop(l_bound_core)
                curve_loop_core_air = gmsh.model.geo.addCurveLoop(l_core_air)
                plane_surface_core.append(gmsh.model.geo.addPlaneSurface([curve_loop_bound_core, curve_loop_core_air]))
            # =====================
            # Core Islands (Core parts between Air Gaps)
            # Points of Core Islands (index refers to sketch)
            if geo.n_air_gaps != 0:
                while island_right.shape[0] > 0:
                    min_island_right = np.argmin(island_right[:, 1])
                    p_island.append(gmsh.model.geo.addPoint(0, island_right[min_island_right, 1], island_right[min_island_right, 2], island_right[min_island_right, 3]))
                    p_island.append(gmsh.model.geo.addPoint(island_right[min_island_right, 0], island_right[min_island_right, 1], island_right[min_island_right, 2], island_right[min_island_right, 3]))
                    island_right = np.delete(island_right, min_island_right, 0)
                # Curves of Core Islands (index refers to sketch)
                for i in range(0, int(len(p_island) / 4)):
                    l_core_air.append(gmsh.model.geo.addLine(p_island[4 * i + 0], p_island[4 * i + 1]))
                    l_core_air.append(gmsh.model.geo.addLine(p_island[4 * i + 1], p_island[4 * i + 3]))
                    l_core_air.append(gmsh.model.geo.addLine(p_island[4 * i + 3], p_island[4 * i + 2]))
                    l_bound_core.append(gmsh.model.geo.addLine(p_island[4 * i + 2], p_island[4 * i + 0]))
                    # Iterative plane creation
                    curve_loop_island.append(gmsh.model.geo.addCurveLoop([l_core_air[-3], l_core_air[-2], l_core_air[-1], l_bound_core[-1]]))
                    plane_surface_core.append(gmsh.model.geo.addPlaneSurface([curve_loop_island[-1]]))

            # Curves: Boundary - Air
            if geo.n_air_gaps == 1:
                l_bound_air.append(gmsh.model.geo.addLine(p_core[0], p_core[5]))
            else:
                for i in range(0, int(len(p_island) / 4)):
                    if i == 0:  # First Line
                        l_bound_air.append(gmsh.model.geo.addLine(p_core[0], p_island[0]))
                    else:  # Middle Lines
                        l_bound_air.append(gmsh.model.geo.addLine(p_island[4*(i-1)+2], p_island[4*i+0]))
                    if i == int(len(p_island) / 4) - 1:  # Last Line
                        l_bound_air.append(gmsh.model.geo.addLine(p_island[-2], p_core[5]))
            # =====================
            # Conductors
            # Points of Conductors
            for i in range(0, geo.p_conductor.shape[0]):
                p_cond.append(gmsh.model.geo.addPoint(geo.p_conductor[i][0], geo.p_conductor[i][1], 0, geo.c_conductor))
            # Curves of Conductors
            for i in range(0, int(len(p_cond) / 4)):
                l_cond.append(gmsh.model.geo.addLine(p_cond[4 * i + 0], p_cond[4 * i + 2]))
                l_cond.append(gmsh.model.geo.addLine(p_cond[4 * i + 2], p_cond[4 * i + 3]))
                l_cond.append(gmsh.model.geo.addLine(p_cond[4 * i + 3], p_cond[4 * i + 1]))
                l_cond.append(gmsh.model.geo.addLine(p_cond[4 * i + 1], p_cond[4 * i + 0]))
                # Iterative plane creation
                curve_loop_cond.append(gmsh.model.geo.addCurveLoop(
                    [l_cond[i * 4 + 0], l_cond[i * 4 + 1], l_cond[i * 4 + 2], l_cond[i * 4 + 3]]))
                plane_surface_cond.append(gmsh.model.geo.addPlaneSurface([curve_loop_cond[i]]))
            # =====================
            # Air (Points are partwise double designated)
            l_air_tmp = l_core_air[:7]
            for i in range(0, len(l_bound_air)):
                l_air_tmp.append(l_bound_air[i])
                if i < len(l_bound_air)-1:
                    l_air_tmp.append(l_core_air[7+3*i])
                    l_air_tmp.append(l_core_air[7+3*i+1])
                    l_air_tmp.append(l_core_air[7+3*i+2])
            curve_loop_air.append(gmsh.model.geo.addCurveLoop(l_air_tmp))
            plane_surface_air.append(gmsh.model.geo.addPlaneSurface(curve_loop_air + curve_loop_cond))

            # =====================
            # Bound
            l_bound_tmp = l_bound_core[:5]
            for i in range(0, len(l_bound_air)):
                l_bound_tmp.append(l_bound_air[-i-1])
                if i != len(l_bound_air)-1:  # last run
                    l_bound_tmp.append(l_bound_core[-i-1])
            # curve_loop_bound.append(gmsh.model.geo.addCurveLoop(l_bound_tmp, reorient=True))


# Define physical Surfaces and Curves
# Core
ps_core = gmsh.model.geo.addPhysicalGroup(2, plane_surface_core, tag=2000)
# Conductors
ps_cond = []
if geo.n_conductors == 2 and geo.conductor_type == "stacked":
    ps_cond[0] = gmsh.model.geo.addPhysicalGroup(2, [plane_surface_cond[0]], tag=4000)
    ps_cond[1] = gmsh.model.geo.addPhysicalGroup(2, [plane_surface_cond[1]], tag=4001)
if geo.conductor_type == "foil":
        for i in range(0, geo.n_conductors):
            ps_cond.append(gmsh.model.geo.addPhysicalGroup(2, [plane_surface_cond[i]], tag=4000+i))
# Air
ps_air = gmsh.model.geo.addPhysicalGroup(2, plane_surface_air, tag=1000)
# Boundary
pc_bound = gmsh.model.geo.addPhysicalGroup(1, l_bound_tmp, tag=1111)

gmsh.model.setPhysicalName(2, ps_core, "CORE")
#gmsh.model.setPhysicalName(2, ps_cond, "COND")
gmsh.model.setPhysicalName(2, ps_air, "AIR")
gmsh.model.setPhysicalName(1, pc_bound, "BOUND")


# Synchronize again
gmsh.model.geo.synchronize()

# Output .msh file
gmsh.option.setNumber("Mesh.SaveAll", 0)
gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
gmsh.option.setNumber("Mesh.SurfaceFaces", 1)


# Colors
for i in range(0, len(plane_surface_core)):
    gmsh.model.setColor([(2, plane_surface_core[i])], 50, 50, 50)
gmsh.model.setColor([(2, plane_surface_air[0])], 0, 0, 230)
for i in range(0, len(plane_surface_cond)):
    gmsh.model.setColor([(2, plane_surface_cond[i])], 150, 150, 0)


# Mesh generation
gmsh.model.mesh.generate(2)

# Parent folder path
path = str(pathlib.Path(__file__).parent.absolute())

# Check operating system
if sys.platform == "linux" or sys.platform == "linux2":
    gmsh.write(path + "/geometry.msh")
elif sys.platform == "darwin":
    # OS X
    gmsh.write(path + "/geometry.msh")
elif sys.platform == "win32":
    gmsh.write(path + "/geometry.msh")  # Win10 can handle slash

# Open gmsh GUI for visualization
#gmsh.fltk.run()

# Terminate gmsh
gmsh.finalize()


# ------------------------------------- File Communication ----------------------------------
# All shared control variables and parameters are passed to a Prolog file
text_file = open("Parameter.pro", "w")


# -- Control Flags --
if geo.flag_excitation_type == 'current':
    text_file.write(f"Flag_ImposedVoltage = 0;\n")
if geo.flag_excitation_type == 'voltage':
    text_file.write(f"Flag_ImposedVoltage = 1;\n")


# -- Geometry --
# Number of conductors
text_file.write("NbrCond = %s;\n" % geo.n_conductors)
# Coordinates of the rectangular winding window
if geo.axi_symmetric == 1:
    text_file.write("Xw1 = %s;\n" % geo.p_window[4, 0])
    text_file.write("Xw2 = %s;\n" % geo.p_window[5, 0])
else:
    raise NotImplementedError("Only axi symmetric case implemented :(")


# -- Materials --

# Nature Constants
text_file.write(f"mu0 = 4.e-7 * Pi;\n")
text_file.write(f"nu0 = 1 / mu0;\n")

# Material Properties
# Conductor Material
text_file.write(f"SigmaCu = 6e7;\n")

# Core Material
if geo.frequency == 0:
    if geo.flag_non_linear_core == 1:
        text_file.write(f"Flag_NL = 1;\n")
        text_file.write(f"Core_Material = {geo.core_material};\n")
    else:
        text_file.write(f"Flag_NL = 0;\n")
        text_file.write(f"mur = {geo.mu_rel};\n")
if geo.frequency != 0:
    text_file.write(f"Flag_NL = 0;\n")
    text_file.write(f"mur = {geo.mu_rel};\n")


# -- Excitation --
# Imposed current, current density or voltage
if geo.flag_excitation_type == 'current':
    text_file.write(f"Val_EE = {geo.current};\n")
if geo.flag_excitation_type == 'current_density':
    text_file.write(f"Val_EE = {geo.current_density};\n")
if geo.flag_excitation_type == 'voltage':
    text_file.write(f"Val_EE = {geo.voltage};\n")

# Frequency and reduced Frequency
text_file.write("Rc = Sqrt[1/Pi]*1e-3;\n")
text_file.write("Flag_imposedRr = %s;\n" % geo.flag_imposed_reduced_frequency)
if geo.flag_imposed_reduced_frequency == 1:
    text_file.write("Rr = %s;\n" % geo.red_freq)
    text_file.write("delta = Rc/Rr;\n")
    text_file.write("Freq  = 1/(delta*delta*mu0*SigmaCu*Pi);\n")
else:
    text_file.write("Freq = %s;\n" % geo.frequency)
    text_file.write("delta = 1/Sqrt[mu0*SigmaCu*Freq*Pi];\n")
    text_file.write("Rr = Rc/delta;\n")


text_file.close()


# ---------------------------------------- Simulation ---------------------------------------
# create a new onelab client
c = onelab.client(__file__)


if os.path.isfile(path + "/config.py"):
    import config
    with open('config.py') as f:
        if 'mygetdp' in f.read():
            if os.path.isfile(config.mygetdp + ".exe") or os.path.isfile(config.mygetdp + ".sh"):
                mygetdp = config.mygetdp
        else:
            mygetdp = call_for_path("mygetdp")
else:
    mygetdp = call_for_path("mygetdp")


# create a onelab variable for the model name
inductor = c.defineString('Inductor model', value='inductor')


# get model file names with correct path
msh_file = c.getPath('geometry.msh')
solver = c.getPath('ind_axi_python_controlled' + '.pro')

# Run simulations as sub clients (non blocking??)
c.runSubClient('myGetDP', mygetdp + ' ' + solver + ' -msh ' + msh_file + ' -solve Analysis -v2')


# ---------------------------------------- Visualization in gmsh ---------------------------------------
gmsh.initialize()
epsilon = 1e-9
# Mesh
gmsh.option.setNumber("Mesh.SurfaceEdges", 0)

# Ohmic losses (weightend effective value of current density)
gmsh.open("res/j2F.pos")
gmsh.option.setNumber("View[0].ScaleType", 2)
gmsh.option.setNumber("View[0].RangeType", 2)
gmsh.option.setNumber("View[0].SaturateValues", 1)
gmsh.option.setNumber("View[0].CustomMin", gmsh.option.getNumber("View[0].Min") + epsilon)
gmsh.option.setNumber("View[0].CustomMax", gmsh.option.getNumber("View[0].Max"))
gmsh.option.setNumber("View[0].ColormapNumber", 1)
gmsh.option.setNumber("View[0].IntervalsType", 2)
gmsh.option.setNumber("View[0].NbIso", 40)

# Magnetic flux density
gmsh.open("res/Magb.pos")
gmsh.option.setNumber("View[1].ScaleType", 1)
gmsh.option.setNumber("View[1].RangeType", 1)
gmsh.option.setNumber("View[1].CustomMin", gmsh.option.getNumber("View[1].Min") + epsilon)
gmsh.option.setNumber("View[1].CustomMax", gmsh.option.getNumber("View[1].Max"))
gmsh.option.setNumber("View[1].ColormapNumber", 1)
gmsh.option.setNumber("View[1].IntervalsType", 2)
gmsh.option.setNumber("View[1].NbIso", 40)

print(gmsh.option.getNumber("View[0].Max"))
gmsh.fltk.run()
gmsh.finalize()