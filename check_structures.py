
import qml
from chemhelp import cheminfo
import numpy as np
from rdkit import Chem

import sys

args = sys.argv[1:]

filename = args[0]

molobjs = cheminfo.read_sdffile(filename)

for i, molobj in enumerate(molobjs):

    molobj = next(molobjs)


    # stat = cheminfo.molobj_optimize(molobj)
    # print(stat)

    dist = Chem.rdmolops.Get3DDistanceMatrix(molobj)
    np.fill_diagonal(dist, 10.0)
    min_dist = np.min(dist)

    if min_dist < 0.01:
        print(i, min_dist)
        smi = cheminfo.molobj_to_smiles(molobj)
        molobj = cheminfo.conformationalsearch(smi)

        dist = Chem.rdmolops.Get3DDistanceMatrix(molobj)
        np.fill_diagonal(dist, 10.0)
        min_dist = np.min(dist)

        print(smi)
        print(min_dist)

    # atoms, coord = cheminfo.molobj_to_xyz(molobj)

    # atoms = list(atoms)
    # many_atoms = [atoms]
    # mbtypes = qml.representations.get_slatm_mbtypes(many_atoms)

    # rep = qml.representations.generate_slatm(coord, atoms, mbtypes)
    # print(cheminfo.molobj_to_smiles(molobj))
    # print(rep.mean())



