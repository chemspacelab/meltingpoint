
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



def molobjfilter(molobj):

    # GetRingInfo
    info = molobj.GetRingInfo()
    n_rings = info.NumRings()

    # if n_rings == 0:
    #     return False

    # if n_rings > 2:
    #     return False

    # really_small_space = [1, 5, 6, 7, 8]
    really_small_space = [1, 6]
    atoms = cheminfo.molobj_to_atoms(molobj)
    if not is_allowed_atoms(atoms,  allowed_atoms=really_small_space):
        return False

    # n_atoms = len(atoms)
    # n_heavy_atoms, = np.where(atoms > 1)
    # n_heavy_atoms = len(n_heavy_atoms)
    #
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

    return True


def valuefilter(line):

    line = line.split()
    line = line[1:]

    if len(line) == 1:
        return True

    line = [float(x) for x in line]
    line = np.array(line)
    mean = np.mean(line)
    dev = np.abs(line - mean)
    dev_max = np.max(dev)

    if dev_max > 20:
        return False

    median = np.median(line)

    if median > 700:
        return False

    if median < 100:
        return False

    return True


def parse():


    return


def main():

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--scratch', action='store', help='', metavar="DIR", default="_tmp_")
    parser.add_argument('--sdf', action='store', help='', metavar="FILE") #, nargs="+", default=[])
    parser.add_argument('--properties', action='store', help='', metavar="FILE") #, nargs="+", default=[])
    parser.add_argument('-j', '--procs', action='store', help='pararallize', metavar="int", default=0, type=int)

    args = parser.parse_args()

    if args.scratch[-1] != "/":
        args.scratch += "/"

    fsdf = gzip.open(args.scratch + "structures.sdf.gz", 'w')
    fprop = open(args.scratch + "properties.csv", 'w')

    molecules = cheminfo.read_sdffile(args.sdf)
    properties = open(args.properties, 'r')

    for molobj, line in zip(molecules, properties):

        status = molobjfilter(molobj)

        if not status:
            continue

        status = valuefilter(line)

        if not status:
            continue

        smiles = cheminfo.molobj_to_smiles(molobj, remove_hs=True)

        print(smiles)

        sdfstr = cheminfo.molobj_to_sdfstr(molobj)
        sdfstr += "$$$$\n"
        fsdf.write(sdfstr.encode())

        fprop.write(line)

    fsdf.close()
    fprop.close()

    properties.close()

    return


if __name__ == '__main__':
    main()

