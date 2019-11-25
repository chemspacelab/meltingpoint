
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle

import misc
from chemhelp import cheminfo


def read_csv(filename, read_header=True, sep=None):

    with open(filename, 'r') as f:

        if read_header:
            header = next(f)

        lines = []

        for line in f:

            line = line.replace("*", "")
            line = line.strip()
            line = line.split(sep)
            lines.append(line)

    if read_header:
        return header, lines

    return lines


# @misc.memory.cache
def main():

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', action='store', help='', metavar='FILE')
    parser.add_argument('-j', '--procs', action='store', help='', type=int, metavar='int', default=0)
    args = parser.parse_args()

    # data = pd.read_csv(args.data, sep=" ")

    header, data = read_csv(args.data, read_header=True)

    data = clean_data(data)

    misc.save_obj(args.data.replace(".txt", ""), data)

    return


def clean_data(listdata):

    data = {}

    for row in listdata:

        idx = row[0]
        smi = row[1]
        value = row[3]
        value = float(value)

        molobj, status = cheminfo.smiles_to_molobj(smi)

        if molobj is None:
            print("error:", smi)
            continue

        smi = cheminfo.molobj_to_smiles(molobj, remove_hs=True)

        if smi not in data:
            data[smi] = []

        data[smi].append(value)


    return data



if __name__ == "__main__":
    main()

