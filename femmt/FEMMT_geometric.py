from onelab import onelab
import gmsh, pathlib, sys
import FEMMT as geo
import numpy as np

# Initialization
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.model.add("geometry")


#---------------------------------- Functions ----------------------------------------
def inner_points(a, b, input_points):
    # This function returns the input points that have a common coordinate as the two
    # interval borders a and b
    [min, max] = [None, None]
    output = input_points
    dim = None
    # Find equal dimension
    for i in range(0, 3):
        if a[i] == b[i] and a[i] != 0:
            dim = i
    if dim == None:
        raise Exception("Given points do not have a common dimension")
    n = 0
    while n < output.shape[0]:
        if a[dim] != output[n, dim]:
            output = np.delete(output, n, 0)
        else:
            n += 1
    if output.shape[0] == 0:
        raise Exception("Not implemented Error: No air gaps between interval borders")
    if output.shape[0] % 2 == 1:
        raise Exception("Odd number of input points")
    if dim == 2:
        raise Exception("Not implemented Error: Only 2D is implemeted")
    dim2 = (dim+1) % 2
    if output.shape[0] >= 2:
        argmax = np.argmax(output[:, dim2])
        output = np.delete(output, argmax, 0)
        argmin = np.argmin(output[:, dim2])
        output = np.delete(output, argmin, 0)
        if output.shape[0] == 0:
            print("Only one air gap in this leg. No island needed.")
    return output

def min_max_inner_points(a, b, input_points):
    # This function returns the input points that have a common coordinate and 
    # the minimum distance from the interval borders.
    [min, max] = [None, None]
    buffer = input_points
    dim = None
    # Find equal dimension
    for i in range(0, 3):
        if a[i] == b[i] and a[i] != 0:
            dim = i
    if dim == None:
        raise Exception("Given points do not have a common dimension")
    n = 0
    while n < buffer.shape[0]:
        if a[dim] != buffer[n, dim]:
            buffer = np.delete(buffer, n, 0)
        else:
            n += 1
    if buffer.shape[0] == 0:
        print("No air gaps between interval borders")
    if buffer.shape[0] % 2 == 1:
        raise Exception("Odd number of input points")
    if dim == 2:
        raise Exception("Not implemented Error: Only 2D is implemeted")
    dim2 = (dim+1) % 2
    if buffer.shape[0] >= 2:
        argmax = np.argmax(buffer[:, 1])
        max = buffer[argmax]
        argmin = np.argmin(buffer[:, 1])
        min = buffer[argmin]
    return [min, max]

#------------------------------------------ Geometry -------------------------------------------
# Core generation
if geo.core_type == "EI":
    #--------------------------------------- Points --------------------------------------------
    if geo.y_symmetric == 1:
        if geo.axi_symmetric == 1:

            # Find points of air gaps (used later)
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
            curve_loop_core = []
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
            p_core.append(gmsh.model.geo.addPoint(0, center_right[0][1], center_right[0][2], center_right[0][3]))
            p_core.append(gmsh.model.geo.addPoint(0, geo.p_outer[1][1], geo.p_outer[1][2], geo.p_outer[1][3]))
            p_core.append(gmsh.model.geo.addPoint(geo.p_outer[1][0], geo.p_outer[1][1], geo.p_outer[1][2], geo.p_outer[1][3]))
            p_core.append(gmsh.model.geo.addPoint(geo.p_outer[3][0], geo.p_outer[3][1], geo.p_outer[3][2], geo.p_outer[3][3]))
            p_core.append(gmsh.model.geo.addPoint(0, geo.p_outer[3][1], geo.p_outer[3][2], geo.p_outer[3][3]))
            p_core.append(gmsh.model.geo.addPoint(0, center_right[1][1], center_right[1][2], center_right[1][3]))
            p_core.append(gmsh.model.geo.addPoint(center_right[1][0], center_right[1][1], center_right[1][2], center_right[1][3]))
            p_core.append(gmsh.model.geo.addPoint(geo.p_window[6][0], geo.p_window[6][1], geo.p_window[6][2], geo.p_window[6][3]))
            p_core.append(gmsh.model.geo.addPoint(geo.p_window[7][0], geo.p_window[7][1], geo.p_window[7][2], geo.p_window[7][3]))
            p_core.append(gmsh.model.geo.addPoint(geo.p_window[5][0], geo.p_window[5][1], geo.p_window[5][2], geo.p_window[5][3]))
            p_core.append(gmsh.model.geo.addPoint(geo.p_window[4][0], geo.p_window[4][1], geo.p_window[4][2], geo.p_window[4][3]))
            p_core.append(gmsh.model.geo.addPoint(center_right[0][0], center_right[0][1], center_right[0][2], center_right[0][3]))
            # Curves of Main Core (index refers to sketch)
# To be added: Case Air Gaps directly on outer leg
            # Curves: Boundary - Core
            l_bound_core.append(gmsh.model.geo.addLine(p_core[0], p_core[1]))
            l_bound_core.append(gmsh.model.geo.addLine(p_core[1], p_core[2]))
            l_bound_core.append(gmsh.model.geo.addLine(p_core[2], p_core[3]))
            l_bound_core.append(gmsh.model.geo.addLine(p_core[3], p_core[4]))
            l_bound_core.append(gmsh.model.geo.addLine(p_core[4], p_core[5]))
            # Curves: Core - Air
            l_core_air.append(gmsh.model.geo.addLine(p_core[5], p_core[6]))
            l_core_air.append(gmsh.model.geo.addLine(p_core[6], p_core[7]))
            l_core_air.append(gmsh.model.geo.addLine(p_core[7], p_core[8]))
            l_core_air.append(gmsh.model.geo.addLine(p_core[8], p_core[9]))
            l_core_air.append(gmsh.model.geo.addLine(p_core[9], p_core[10]))
            l_core_air.append(gmsh.model.geo.addLine(p_core[10], p_core[11]))
            l_core_air.append(gmsh.model.geo.addLine(p_core[11], p_core[0]))
            # Plane: Main Core --> plane_surface_core[0]
            curve_loop_core = gmsh.model.geo.addCurveLoop(l_bound_core + l_core_air)
            plane_surface_core.append(gmsh.model.geo.addPlaneSurface([curve_loop_core]))

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
            #curve_loop_bound.append(gmsh.model.geo.addCurveLoop(l_bound_tmp, reorient=True))



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

# Colors
#for i in range(0, len(plane_surface_core)):
#    gmsh.model.setColor([(2, plane_surface_core[i])], 50, 50, 50)
#gmsh.model.setColor([(2, plane_surface_air)], 0, 0, 150)
#gmsh.model.setColor([(1, l_bound[0]), (1, l_bound[1]), (1, l_bound[2]), (1, l_bound[3])], 200, 0, 0)


gmsh.model.mesh.generate(2)


path = str(pathlib.Path(__file__).parent.absolute())
gmsh.write(path + "\\geometry.msh")
gmsh.fltk.run()
gmsh.finalize()



#----------------------------------------Simulation---------------------------------------
# create a new onelab client
c = onelab.client(__file__)

# create a onelab variable for the model name
simple_inductor = c.defineString('Simple Inductor Model', value='simple_inductor')

# get model file names with correct path
mygmsh = 'C:/onelab-Windows64/gmsh'
mygetdp = 'C:/onelab-Windows64/getdp'

for s in range(9):
   n = c.getString('Solver.Name' + str(s))
   if(n == 'GetDP'):
      mygetdp = c.getString('Solver.Executable' + str(s))
      break
if(not len(mygetdp)):
   c.sendError('This appears to be the first time you are trying to run GetDP')
   c.sendError('Please run a GetDP model interactively once with Gmsh to ' +
               'initialize the solver location')
   exit(0)

c.sendInfo('Will use gmsh={0} and getdp={1}'.format(mygmsh, mygetdp))


# create a onelab variable for the model name
inductor = c.defineString('Inductor model', value='inductor')

# we're done if we don't do the actual calculation
if c.action == 'check' :
   exit(0)

# get model file names with correct path
#inductor_geo = c.getPath(inductor + '.geo') #not needed if mesh is directly created with python
msh_file = c.getPath('geometry.msh')
solver = c.getPath('ind_axi_python_controlled' + '.pro')


# run gmsh as a subclient
# not needed if mesh is directly created with python
#c.runSubClient('myGmsh', mygmsh + ' ' + inductor_geo + ' -2 -v 2')

# run getdp as a subclient
c.runSubClient('myGetDP',  mygetdp + ' ' + solver + ' -msh ' + msh_file + ' -solve Analysis -v2')

