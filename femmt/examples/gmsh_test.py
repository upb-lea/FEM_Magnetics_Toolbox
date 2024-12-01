import gmsh
import matplotlib.pyplot as plt
import numpy as np
import re
from matplotlib import pyplot
import time

# nubmer of air gasp
no_of_air_gaps = 4
# number of points that are necessary to describe the outer line of the core (the points are the corners of the core)
points_of_core = 10 + 4 * no_of_air_gaps

gmsh.initialize()
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
start = time.time()
start_tot = time.time()

# find both occurrence of "Magnitude B-Field / T" in Magb.pos
indexes = [m.start() for m in re.finditer('"Magnitude B-Field / T"', content)]
indexes_Elements = [m.start() for m in re.finditer('Elements', content)]
# filter magnetic flux density out of Magb.pos(data is between the both occurrences of "Magnitude B-Field / T") and put data into list
data = [float(x) for x in content[indexes[0]:indexes[1]].split()[11:-3]]
data_Elements = [float(x) for x in content[indexes_Elements[0]:indexes_Elements[1]].split()[2:-1]]
# divide long list into lists for every triangle/mesh-cell
chunks = [data[x:x+5] for x in range(0, len(data), 5)]
chunks_Elements = [data_Elements[x:x+8] for x in range(0, len(data_Elements), 8)]
print(len(chunks))
print(len(chunks_Elements))
# put lists in right format -> [TriangleTag, MagFluxDenOfNode1, MagFluxDenOfNode2, MagFluxDenOfNode3]
# before [TriangleTag, 3, MagFluxDenOfNode1, MagFluxDenOfNode2, MagFluxDenOfNode3]
TriangleMagneticFluxDensity = [[x[0], x[2], x[3], x[4]] for x in chunks]
# print(TriangleMagneticFluxDensity)
# print(len(TriangleMagneticFluxDensity))

# get all triangles with there tags and node and put them in right format -> [TriangleTag, NodeTag1, NodeTag2, NodeTag3]
TriangleTags, TriangleNodeTags = gmsh.model.mesh.getElementsByType(elementType=2)
TriangleNodeTags = TriangleNodeTags.reshape((-1, 3))
TriangleNodeTags = [[TriangleTags[index], value[0], value[1], value[2]] for index, value in enumerate(TriangleNodeTags)]

# append the average of the magnetic flux density on every triangle
for index, TriangleNode in enumerate(TriangleNodeTags):
    TriangleNodeTags[index].append(np.mean(chunks[index][2:]))
    TriangleNodeTags[index].append(chunks_Elements[index][3])
# print(TriangleNodeTags[:100])
# print(len(TriangleNodeTags))

Triangle_Core = []
Triangle_Core_All_Nodes = []
for TriangleNode in TriangleNodeTags:
    if TriangleNode[-1] // 10000 == 12:
        Triangle_Core.append(TriangleNode)
        Triangle_Core_All_Nodes.append(TriangleNode[-5:-2])

Triangle_Core_All_Nodes = set(x for xs in Triangle_Core_All_Nodes for x in xs)
print(Triangle_Core_All_Nodes)
print(len(Triangle_Core_All_Nodes))


# print(Triangle_Core[0:100])
# print(Triangle_Core[-100:])
print(len(Triangle_Core))
print("2", time.time()-start)
start = time.time()

# create a dict with the NodeTag as key and the corresponding value of the magnetic flux density of the node as value e.g. {755: 0.00547812 ...}
NodeMagneticFlux = {}
for index, node in enumerate(TriangleNodeTags):
    NodeMagneticFlux[node[1]] = TriangleMagneticFluxDensity[index][1]
    NodeMagneticFlux[node[2]] = TriangleMagneticFluxDensity[index][2]
    NodeMagneticFlux[node[3]] = TriangleMagneticFluxDensity[index][3]
# print(NodeMagneticFlux)
# print(len(NodeMagneticFlux))

# get the coordinates for all nodes
NodeTags, coords, parametricCoord = gmsh.model.mesh.getNodes()
coords = coords.reshape((-1, 3))
NodeCoords = [[NodeTags[index], value[0], value[1], value[2]] for index, value in enumerate(coords)]

# get the coordinates of the corners of the core
xyz = coords.reshape(-1, 3)
coords_of_core = xyz[0:points_of_core]
x_coords = sorted(set([a[0] for a in coords_of_core]))
y_coords = sorted(set([a[1] for a in coords_of_core]))
# print(x_coords)
# print(y_coords)

print("3", time.time()-start)
start = time.time()
core_nodes = NodeCoords.copy()
# get all nodes that are not positioned in the core or on the boundaries of the core
# nodes_outside_of_core = []
for node in NodeCoords:
    for i in range(no_of_air_gaps+1):
        if i == 0:
            if (x_coords[1] < node[1] < x_coords[2]) and (y_coords[1] < node[2] < y_coords[-2]):  # x1 < x < x2 AND -y1 < y < y1
                core_nodes.remove(node)
                # nodes_outside_of_core.append(node)
        else:
            if not (i*2+3 == no_of_air_gaps):
                if (0 <= node[1] <= x_coords[1]) and (y_coords[i*2] < node[2] < y_coords[i*2+1]):  # 0 <= x <= x1 AND yx < y < yx+1
                    core_nodes.remove(node)
                    # nodes_outside_of_core.append(node)
# print(nodes_outside_of_core)

print("4", time.time()-start)
start = time.time()

# get all nodes that are positioned in the core
# core_nodes = NodeCoords.copy()
# for node in nodes_outside_of_core:
#     core_nodes.remove(node)

print("5", time.time()-start)
start = time.time()

# append the magnetic flux density on the core_nodes [NodeTag, x-coordinate, y-coordinate, z-coordinate, magnetic flux density]
for index, node in enumerate(core_nodes):
    core_nodes[index].append(NodeMagneticFlux[node[0]])
# print(len(nodes_outside_of_core))
# print(len(core_nodes))
# print(len(core_nodes) + len(nodes_outside_of_core))
# print(len(NodeCoords))

print("6", time.time()-start)
start = time.time()

# get all triangle that are inside the core
core_node_Tags = [x[0] for x in core_nodes]
core_triangles = []
for TriangleNode in TriangleNodeTags:
    if TriangleNode[1] in core_node_Tags:
        if TriangleNode[2] in core_node_Tags:
            if TriangleNode[3] in core_node_Tags:
                core_triangles.append(TriangleNode)
# print(core_triangles[0:100])
# print(core_triangles[-100:])
# print([core[0] for core in core_triangles])
# print(len(core_triangles))

print("7", time.time()-start)
print("total", time.time()-start_tot)

surface_area = 0
for index, TriangleNode in enumerate(core_triangles):

    xi = xyz[int(TriangleNode[1:4][0]-1), :]
    xj = xyz[int(TriangleNode[1:4][1]-1), :]
    xk = xyz[int(TriangleNode[1:4][2]-1), :]

    a = np.sqrt(np.dot(xi - xj, xi - xj))
    b = np.sqrt(np.dot(xi - xk, xi - xk))
    c = np.sqrt(np.dot(xj - xk, xj - xk))

    s = (a+b+c)/2
    dA = np.sqrt(s*(s-a)*(s-b)*(s-c))

    core_triangles[index].append(dA)

    surface_area += dA
print("area of core:", surface_area)
# print(core_triangles)
print(len(core_triangles))

x_coords_core = [x[1] for x in core_nodes]
y_coords_core = [x[2] for x in core_nodes]
flux_density = [x[4] for x in core_nodes]


# for tag, xyz_e in zip(nodeTags, xyz):
#     print(f"Node # {tag} is at {xyz_e}")
# print("Physical gmsh: ", gmsh.model.getPhysicalGroups())

fig, ax = plt.subplots(1, 1)
# # ax.plot(x_coords_core, y_coords_core, marker="o", markersize=2, linestyle="None")
plt.scatter(x_coords_core, y_coords_core, c=flux_density, cmap='viridis', linewidths=0.1)
plt.colorbar()
plt.grid()
plt.show(block=False)
# # plt.draw()


surface_area = 0.0
TriangleTags, TriangleNodeTags = gmsh.model.mesh.getElementsByType(elementType=2)
TriangleNodeTags = TriangleNodeTags.reshape((-1, 3))
for tag, nodes in zip(TriangleTags, TriangleNodeTags):
    # print(f"Element # {tag} has nodes {nodes}")

    xi = xyz[int(nodes[0]-1), :]
    xj = xyz[int(nodes[1]-1), :]
    xk = xyz[int(nodes[2]-1), :]

    a = np.sqrt(np.dot(xi - xj, xi - xj))
    b = np.sqrt(np.dot(xi - xk, xi - xk))
    c = np.sqrt(np.dot(xj - xk, xj - xk))

    s = (a+b+c)/2
    dA = np.sqrt(s*(s-a)*(s-b)*(s-c))

    surface_area += dA
print(surface_area)

# defined_element_type = gmsh.model.mesh.getElementTypes()
# print(f"{defined_element_type=}")

# Points
# eleTags, nodeTags = gmsh.model.mesh.getElementsByType(elementType=1)
# print(eleTags)
# print(nodeTags)

# Lines
# eleTags, nodeTags = gmsh.model.mesh.getElementsByType(elementType=15)
# print(eleTags)
# print(nodeTags)

gmsh.fltk.run()
gmsh.finalize()
