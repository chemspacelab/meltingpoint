
import gzip

import multiprocessing as mp
import matplotlib.pyplot as plt
import numpy as np
import rdkit.Chem as Chem

import misc
from chemhelp import cheminfo

ALLOWED_ATOMS = [
    1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 27, 35, 53
]


def is_mol_allowed(atoms, allowed_atoms=ALLOWED_ATOMS):

    atoms = np.unique(atoms)

    for atom in atoms:
        if atom not in ALLOWED_ATOMS:
            return False

    return True


def main():

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--scratch', action='store', help='', metavar="DIR", default="_tmp_")
    parser.add_argument('--json', action='store', help='', metavar="FILE")
    parser.add_argument('-j', '--procs', action='store', help='pararallize', metavar="int", default=0, type=int)

    args = parser.parse_args()

    if args.scratch[-1] != "/":
        args.scratch += "/"

    data = misc.load_json(args.json)

    keys = data.keys()
    keys = list(keys)

    canonical_data = {}

    for key in keys:

        molobj, status = cheminfo.smiles_to_molobj(key)

        if molobj is None:
            print("error none mol:", key)
            continue

        smiles = cheminfo.molobj_to_smiles(molobj, remove_hs=True)

        if "." in smiles:
            print("error multi mol:", smiles)
            continue

        atoms = cheminfo.molobj_to_atoms(molobj)

        if not is_mol_allowed(atoms):
            print("error heavy mol:", smiles)
            continue

        canonical_data[smiles] = data[key]


    misc.save_json(args.scratch + "molecule_data", canonical_data)
    misc.save_obj(args.scratch + "molecule_data", canonical_data)

    return


if __name__ == '__main__':
    main()
