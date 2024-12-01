import gmsh
import matplotlib.pyplot as plt
import numpy as np
import re
from matplotlib import pyplot
import time

PLOT = True

gmsh.initialize()
start_tot = time.time()
start = time.time()
# path of gmsh-file needed
# gmsh.open("C:/Users/schacht/PycharmProjects/FEM_Magnetics_Toolbox/femmt/examples/example_results/basic_inductor/mesh/electro_magnetic.msh")
gmsh.open("C:/Users/sebas/Desktop/Python/FEM_Magnetics_Toolbox/femmt/examples/example_results/basic_inductor/mesh/electro_magnetic.msh")

# path of Magb.pos needed (file containing the data of the magnetic flux density for the mesh)
# path = 'C:/Users/schacht/PycharmProjects/FEM_Magnetics_Toolbox/femmt/examples/example_results/basic_inductor/results/fields/Magb.pos'
path = 'C:/Users/sebas/Desktop/Python/FEM_Magnetics_Toolbox/femmt/examples/example_results/basic_inductor/results/fields/Magb.pos'

# read in Magb.pos <- data of magnetic flux density
with open(path, 'r') as file:
    content = file.read()  # type of content: string

print("1", time.time()-start)
# start = time.time()

# find both occurrence of "Magnitude B-Field / T" in Magb.pos
indexes_B_Field = [m.start() for m in re.finditer('"Magnitude B-Field / T"', content)]
indexes_Elements = [m.start() for m in re.finditer('Elements', content)]
# filter magnetic flux density out of Magb.pos(data is between the both occurrences of "Magnitude B-Field / T") and put data into list
data_B_Field = [float(x) for x in content[indexes_B_Field[0]:indexes_B_Field[1]].split()[11:-3]]
data_Elements = [float(x) for x in content[indexes_Elements[0]:indexes_Elements[1]].split()[2:-1]]
# divide long list into lists for every triangle/mesh-cell
chunks_B_Field = [data_B_Field[x:x+5] for x in range(0, len(data_B_Field), 5)]
chunks_Elements = [data_Elements[x:x+8] for x in range(0, len(data_Elements), 8)]

# print("2", time.time()-start)
# start = time.time()

# get all triangles with there tags and node and put them in right format -> [TriangleTag, NodeTag1, NodeTag2, NodeTag3]
TriangleTags, TriangleNodeTags = gmsh.model.mesh.getElementsByType(elementType=2)
TriangleNodeTags = TriangleNodeTags.reshape((-1, 3))
TriangleNodeTags = [[TriangleTags[index], value[0], value[1], value[2]] for index, value in enumerate(TriangleNodeTags)]

# append the average of the magnetic flux density on every triangle
Triangle_Core = []
for index, TriangleNode in enumerate(TriangleNodeTags):
    TriangleNodeTags[index].append(np.mean(chunks_B_Field[index][2:]))
    TriangleNodeTags[index].append(chunks_Elements[index][3])

Triangle_Core = []
for TriangleNode in TriangleNodeTags:
    if TriangleNode[-1] // 10000 == 12:
        Triangle_Core.append(TriangleNode.copy())
# get the coordinates for all nodes
NodeTags, coords, parametricCoord = gmsh.model.mesh.getNodes()
coords = coords.reshape((-1, 3))
NodeCoords = [[NodeTags[index], value[0], value[1], value[2]] for index, value in enumerate(coords)]
xyz = coords.reshape(-1, 3)

surface_area = 0
for index, TriangleNode in enumerate(Triangle_Core):

    xi = xyz[int(TriangleNode[1:4][0]-1), :]
    xj = xyz[int(TriangleNode[1:4][1]-1), :]
    xk = xyz[int(TriangleNode[1:4][2]-1), :]

    a = np.sqrt(np.dot(xi - xj, xi - xj))
    b = np.sqrt(np.dot(xi - xk, xi - xk))
    c = np.sqrt(np.dot(xj - xk, xj - xk))

    s = (a+b+c)/2
    dA = np.sqrt(s*(s-a)*(s-b)*(s-c))

    Triangle_Core[index].append(dA)

    surface_area += dA
# print("area of core:", surface_area)
# print("3", time.time()-start)
# print(Triangle_Core)

b_field_avg = np.sum([x[-3]*x[-1] for x in Triangle_Core])/surface_area

print("tot", time.time()-start_tot)
print(b_field_avg)

if PLOT:
    Triangle_Core_All_Nodes = []
    for TriangleNode in TriangleNodeTags:
        if TriangleNode[-1] // 10000 == 12:
            Triangle_Core_All_Nodes.append(TriangleNode[-5:-2])
    # ALL NODES FOR TRIANGLES IN CORE
    Triangle_Core_All_Nodes = set(x for xs in Triangle_Core_All_Nodes for x in xs)

    NodeTags, coords, parametricCoord = gmsh.model.mesh.getNodes()
    coords = coords.reshape((-1, 3))
    NodeCoords = {}
    for index, value in enumerate(coords):
        NodeCoords[NodeTags[index]] = [value[0], value[1], value[2]]
    NodeMagneticFlux = {}
    TriangleMagneticFluxDensity = [[x[0], x[2], x[3], x[4]] for x in chunks_B_Field]
    for index, node in enumerate(TriangleNodeTags):
        NodeMagneticFlux[node[1]] = TriangleMagneticFluxDensity[index][1]
        NodeMagneticFlux[node[2]] = TriangleMagneticFluxDensity[index][2]
        NodeMagneticFlux[node[3]] = TriangleMagneticFluxDensity[index][3]
    Triangle_Core_B_And_Coords = []
    for TriangleNode in Triangle_Core_All_Nodes:
        a = [TriangleNode, NodeCoords[TriangleNode], NodeMagneticFlux[TriangleNode]]
        Triangle_Core_B_And_Coords.append(a)
    x_coords_core = [x[1][0] for x in Triangle_Core_B_And_Coords]
    y_coords_core = [x[1][1] for x in Triangle_Core_B_And_Coords]
    flux_density = [x[-1] for x in Triangle_Core_B_And_Coords]

    fig, ax = plt.subplots(1, 1)
    # # ax.plot(x_coords_core, y_coords_core, marker="o", markersize=2, linestyle="None")
    plt.scatter(x_coords_core, y_coords_core, c=flux_density, cmap='viridis', linewidths=0.1)
    plt.colorbar()
    plt.grid()
    plt.show(block=False)
    # # plt.draw()

gmsh.fltk.run()
gmsh.finalize()

