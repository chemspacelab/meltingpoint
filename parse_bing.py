
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle

import misc
from chemhelp import cheminfo


def read_csv(filename, read_header=True, sep=" "):


    with open(filename, 'r') as f:

        if read_header:
            header = next(f)

        lines = []

        for line in f:

            line = line.strip()
            line = line.split(sep)
            lines.append(line)

    if read_header:
        return header, lines

    return lines


@misc.memory.cache
def main():

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', action='store', help='', metavar='FILE')
    parser.add_argument('-j', '--procs', action='store', help='', type=int, metavar='int', default=0)
    args = parser.parse_args()

    # data = pd.read_csv(args.data, sep=" ")

    data = read_csv(args.data)

    print(data)

    return data


def clean_data(df):

    smiles = df.iloc[1]

    data = {}

    for index, row in df.iterrows():

        print(row)
        row = list(row)
        print(row)
        quit()

        smi = row.iloc[1]
        # value = row.mpC + 273.15

        molobj, status = cheminfo.smiles_to_molobj(smi)

        if molobj is None:
            print("error:", smi)
            continue

        smi = cheminfo.molobj_to_smiles(molobj, remove_hs=True)

        if smi not in data:
            data[smi] = []

        data[smi].append(value)


    return



if __name__ == "__main__":
    pd = main()
    data = clean_data(pd)

