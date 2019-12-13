
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

    for key in keys:

        molobj, status = cheminfo.smiles_to_molobj(key)

        if molobj is None:
            print("delete", key)
            del data[key]


    misc.save_json(args.scratch + "molecule_data", data)
    misc.save_obj(args.scratch + "molecule_data", data)

    return


if __name__ == '__main__':
    main()
