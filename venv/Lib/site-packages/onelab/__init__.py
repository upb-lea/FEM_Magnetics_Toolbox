name = "onelab"

import logging
import os
import numpy as np
import pandas as pd


class Mesh():
    """a mesh container"""

    def __init__(self):
        """create a Mesh instance with nodes and elements"""

        self.elements = pd.DataFrame(columns=['pid', 'elid', 'type', 'n_nodes', 'nodes', 'nidxs'])
        self.nodes = pd.DataFrame(columns=['nid', 'x', 'y', 'z'])

    @classmethod
    def from_gmsh(cls, gmsh):
        """create mesh from gmsh"""

        elements = pd.DataFrame(columns=['pid', 'elid', 'type', 'n_nodes', 'nodes'])
        etypes, elids, enids = gmsh.model.mesh.getElements()
        typedict = dict(zip(etypes, range(len(etypes))))
        logging.debug(typedict)
        # quads
        idx = typedict.get(3)
        if idx is not None:
            quad_nids = elids[idx]
            quad_nodes = np.array(enids[idx])
            quad_nodes = quad_nodes.reshape(len(quad_nodes) // 4, 4).tolist()
            quads = pd.DataFrame()
            quads["elid"] = quad_nids
            quads["n_nodes"] = 4
            quads["nodes"] = quad_nodes
            quads["nidxs"] = quad_nodes
            quads["type"] = "shell4"
            quads["pid"] = 1
            elements = elements.append(quads, sort=False)
        # tria
        idx = typedict.get(2)
        if idx is not None:
            tria_nids = elids[idx]
            tria_nodes = np.array(enids[idx])
            tria_nodes = tria_nodes.reshape(len(tria_nodes) // 3, 3).tolist()
            tria_nodes
            trias = pd.DataFrame()
            trias["elid"] = tria_nids
            trias["n_nodes"] = 3
            trias["nodes"] = tria_nodes
            trias["nidxs"] = tria_nodes
            trias["type"] = "shell3"
            trias["pid"] = 1
            elements = elements.append(trias, sort=False)

        nids, coord, parametric_coord = gmsh.model.mesh.getNodes()
        coord = np.array(coord)
        coord = coord.reshape(len(coord) // 3, 3)
        nodes = pd.DataFrame({'nid': nids, 'x': coord[:, 0], 'y': coord[:, 1], 'z': coord[:, 2]})
        nodes.index = nodes.nid.values.copy()

        elements.reset_index(inplace=True, drop=True)

        self = cls()
        self.elements = elements
        self.nodes = nodes
        return self

    def __str__(self):
        return "(Mesh nodes:{} elements:{})".format(len(self.elements), len(self.nodes))

    __repr__ = __str__
