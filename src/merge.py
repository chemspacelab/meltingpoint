
try:
    from matplotlib_venn import venn2, venn3
except:
    venn2 = None
    venn3 = None

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

def is_allowed_atoms(atoms, allowed_atoms=ALLOWED_ATOMS):

    atoms = np.unique(atoms)

    for atom in atoms:
        if atom not in allowed_atoms:
            return False

    return True


def filter_value(values):

    kelvin = 273.15

    line = np.array(values)
    mean = np.mean(line)
    dev = np.abs(line - mean)
    dev_max = np.max(dev)

    if dev_max > 10:
        return False

    if mean > 250+kelvin:
        return False

    if mean < 50+kelvin:
        return False

    return True


def filter_molobj(molobj):

    # GetRingInfo
    info = molobj.GetRingInfo()
    n_rings = info.NumRings()

    # if n_rings == 0:
    #     return False

    # if n_rings > 2:
    #     return False

    atoms = cheminfo.molobj_to_atoms(molobj)
    if not is_allowed_atoms(atoms):
        return False

    n_atoms = len(atoms)
    n_heavy_atoms, = np.where(atoms > 1)
    n_heavy_atoms = len(n_heavy_atoms)

    # # no long chains
    # aromatic_atoms = molobj.GetAromaticAtoms()
    # aromatic_atoms = [atom for atom in aromatic_atoms]
    # aromatic_atoms = [atom.GetAtomicNum() for atom in aromatic_atoms]
    # n_atomatic_atoms = len(aromatic_atoms)
    #
    # n_non_aromatic_atoms = n_heavy_atoms - n_atomatic_atoms 
    #
    # if n_non_aromatic_atoms > 7:
    #     return False

    if n_heavy_atoms < 10:
        return False

    if n_heavy_atoms > 20:
        return False

    if n_atoms > 40:
        return False

    return True


def filter_dict(molecules):

    keys = molecules.keys()
    keys = list(keys)

    max_atoms = 0

    for key in keys:

        molobj, status = cheminfo.smiles_to_molobj(key)

        if molobj is None:
            continue

        status = filter_molobj(molobj)


        if not status:
            del molecules[key]
            print(key, status)
            continue

        status = filter_value(molecules[key])

        if not status:
            print(status, key, molecules[key])
            del molecules[key]
            continue

        # Report
        atoms = cheminfo.molobj_to_atoms(molobj)
        n_atoms = len(atoms)

        if n_atoms > max_atoms:
            max_atoms = n_atoms

        continue


    print("max atoms: ", max_atoms)


    return molecules



def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--scratch', action='store', help='', metavar="DIR", default="_tmp_")
    parser.add_argument('--sdf', action='store', help='', metavar="FILE", nargs="+", default=[])
    parser.add_argument('--dict', action='store', help='', metavar="FILE", nargs="+", default=[])
    parser.add_argument('--name', action='store', help='', metavar="STR", nargs="+")
    parser.add_argument('--filename', action='store', help='', metavar="STR")
    parser.add_argument('--filter', action='store_true', help='')
    parser.add_argument('-j', '--procs', action='store', help='pararallize', metavar="int", default=0, type=int)

    args = parser.parse_args()

    if args.scratch[-1] != "/":
        args.scratch += "/"

    print()
    databases_set = []
    databases_dict = []

    for sdf in args.sdf:
        molobjs = cheminfo.read_sdffile(sdf)
        molobjs = list(molobjs)
        smiles = [cheminfo.molobj_to_smiles(molobj, remove_hs=True) for molobj in molobjs]
        smiles = set(smiles)
        databases_set.append(smiles)
        print(sdf, len(smiles))

    for filename in args.dict:
        data = misc.load_obj(filename)
        smiles = data.keys()
        smiles = set(smiles)
        databases_set.append(smiles)
        databases_dict.append(data)
        print(filename, len(smiles))


    if args.scratch is not None:

        # Merge databases
        everything = {}

        for data in databases_dict:

            keys = data.keys()

            for key in keys:

                if key not in everything:
                    everything[key] = []

                everything[key] += data[key]


        if args.filter:
            everything = filter_dict(everything)

        keys = everything.keys()
        print("n items", len(keys))

        # Save
        misc.save_json(args.scratch + "molecule_data", everything)
        misc.save_obj(args.scratch + "molecule_data", everything)


    if args.name is not None:

        n_db = len(databases_set)

        if n_db == 2:
            venn2(databases_set, set_labels=args.name)
        elif n_db == 3:
            venn3(databases_set, set_labels=args.name)

        plt.savefig(args.filename)


    return


if __name__ == '__main__':
    main()

