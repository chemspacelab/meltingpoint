
from matplotlib_venn import venn2, venn3

import gzip

import multiprocessing as mp
import matplotlib.pyplot as plt
import numpy as np
import rdkit.Chem as Chem

import misc
from chemhelp import cheminfo


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--scratch', action='store', help='', metavar="DIR", default="_tmp_")
    parser.add_argument('--sdf', action='store', help='', metavar="FILE", nargs="+")
    parser.add_argument('--name', action='store', help='', metavar="STR", nargs="+")
    parser.add_argument('--filename', action='store', help='', metavar="STR")
    parser.add_argument('--dict', action='store', help='', metavar="FILE", nargs="+")
    parser.add_argument('-j', '--procs', action='store', help='pararallize', metavar="int", default=0, type=int)

    args = parser.parse_args()

    databases = []

    for sdf in args.sdf:
        molobjs = cheminfo.read_sdffile(sdf)
        molobjs = list(molobjs)
        smiles = [cheminfo.molobj_to_smiles(molobj, remove_hs=True) for molobj in molobjs]
        smiles = set(smiles)
        databases.append(smiles)

        print(sdf, len(smiles))

    for filename in args.dict:
        data = misc.load_obj(filename)
        smiles = data.keys()
        smiles = set(smiles)
        databases.append(smiles)

        print(filename, len(smiles))


    n_db = len(databases)

    if n_db == 2:
        venn2(databases, set_labels=args.name)
    elif n_db == 3:
        venn3(databases, set_labels=args.name)

    plt.savefig(args.filename)


    return


if __name__ == '__main__':
    main()

