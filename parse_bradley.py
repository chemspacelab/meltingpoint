
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle

import misc
from chemhelp import cheminfo

ALLOWED_ATOMS = [
    1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 27, 35, 53
]

@misc.memory.cache
def main():


    xsl = pd.read_excel("data/BradleyMeltingPointDataset.xlsx")
    col_names = [col for col in xsl]

    data = xsl.loc[xsl["donotuse"] != "x"]

    data = data.drop(columns="donotuse")
    data = data.drop(columns="donotusebecause")
    data = data.drop(columns="name")
    data = data.drop(columns="csid")
    data = data.drop(columns="link")
    data = data.drop(columns="source")

    data.to_csv("data/melting_bradley.csv", sep=",", index=False)

    return data


def clean_data(df, scratch):

    smiles = df.iloc[1]

    data = {}

    atom_types = []

    for index, row in df.iterrows():

        smi = row.smiles
        value = row.mpC + 273.15

        molobj, status = cheminfo.smiles_to_molobj(smi)

        if molobj is None:
            print("error:", smi)
            continue

        smi = cheminfo.molobj_to_smiles(molobj, remove_hs=True)

        # Atoms
        atoms = cheminfo.molobj_to_atoms(molobj)
        atom_types += list(atoms)

        if smi not in data:
            data[smi] = []

        data[smi].append(value)

    atom_types, counts = np.unique(atom_types, return_counts=True)

    for atom, count in zip(atom_types, counts):
        print(atom, count)

    misc.save_obj(scratch + "molecule_data", data)
    misc.save_json(scratch + "molecule_data", data)

    return

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--scratch', action='store', help='', metavar="DIR", default="_tmp_")
    parser.add_argument('-j', '--procs', action='store', help='pararallize', metavar="int", default=0, type=int)

    args = parser.parse_args()

    if args.scratch[-1] != "/":
        args.scratch += "/"

    pd = main()
    data = clean_data(pd, args.scratch)

