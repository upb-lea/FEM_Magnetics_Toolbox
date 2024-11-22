import gmsh
import numpy as np
import re

no_of_air_gaps = 1
points_of_core = 10 + 4 * no_of_air_gaps

gmsh.initialize()
gmsh.open("C:/Users/schacht/PycharmProjects/FEM_Magnetics_Toolbox/femmt/examples/example_results/basic_inductor/mesh/electro_magnetic.msh")

dim = -1
tag = -1

nodeTags, coords, parametricCoord = gmsh.model.mesh.getNodes(dim, tag)

coords = coords.reshape((-1, 3))

# print(nodeTags)
# print(coords)
xyz = coords.reshape(-1, 3)
coords_of_core = xyz[0:points_of_core]

x_coords = [a[0] for a in coords_of_core]
y_coords = [a[1] for a in coords_of_core]
print(sorted(set(x_coords)))
print(sorted(set(y_coords)))
print(sorted(set(abs(np.array(y_coords)))))
print(xyz[0:points_of_core+1])
print(xyz[755])
print(xyz[832])
print(xyz[931])
print(xyz[5598])
print(xyz[5630])
print("test")
print(xyz[5630])
print(list(gmsh.model.mesh.getNode(5631)))
# for tag, xyz_e in zip(nodeTags, xyz):
#     print(f"Node # {tag} is at {xyz_e}")

defined_element_type = gmsh.model.mesh.getElementTypes()
print(f"{defined_element_type=}")

eletype = 2
tag = -1
eleTags, nodeTags = gmsh.model.mesh.getElementsByType(eletype, tag)

nodeTags = nodeTags.reshape((-1, 3))

print(eleTags)
print(nodeTags)

surface_area = 0.0

for tag, nodes in zip(eleTags, nodeTags):
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

eletype = 15
tag = -1
eleTags, nodeTags = gmsh.model.mesh.getElementsByType(eletype, tag)

# print(eleTags)

eletype = 1
tag = -1
eleTags, nodeTags = gmsh.model.mesh.getElementsByType(eletype, tag)

# print(eleTags)
# print(nodeTags)

gmsh.fltk.run()

gmsh.finalize()


path = 'C:/Users/schacht/PycharmProjects/FEM_Magnetics_Toolbox/femmt/examples/example_results/basic_inductor/results/fields/Magb.pos'

with open(path, 'r') as file:
    content = file.read()

# print(content)
# print(type(content))
# print(len(content))
# print([m.start() for m in re.finditer('"Magnitude B-Field / T"', content)])
indexes = [m.start() for m in re.finditer('"Magnitude B-Field / T"', content)]
# print(content[indexes[0]:indexes[1]].split()[11:-3])
print(len(content[indexes[0]:indexes[1]].split()[11:-3]))
data = [float(x) for x in content[indexes[0]:indexes[1]].split()[11:-3]]
chunks = [data[x:x+5] for x in range(0, len(data), 5)]
print(chunks)
print(len(chunks))

# print([i for i, x in enumerate(content[2793144:6028848].split()[11:-3]) if x == "3"])
# print(content.find('"Magnitude B-Field / T"'))
