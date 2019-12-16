
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


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--scratch', action='store', help='', metavar="DIR", default="_tmp_")
    parser.add_argument('--sdf', action='store', help='', metavar="FILE", nargs="+", default=[])
    parser.add_argument('--dict', action='store', help='', metavar="FILE", nargs="+", default=[])
    parser.add_argument('--name', action='store', help='', metavar="STR", nargs="+")
    parser.add_argument('--filename', action='store', help='', metavar="STR")
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

        keys = everything.keys()

        print()
        print(args.scratch, len(keys))

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

