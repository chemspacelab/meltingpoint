

import rdkit
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import rdkit.Chem.ChemicalForceFields as ChemicalForceFields
import rdkit.Chem.rdMolDescriptors as rdMolDescriptors

from chemhelp import cheminfo
import glob

import matplotlib.pyplot as plt
import numpy as np

import misc
import views


def canonical(smiles):
    """
    Translate smiles into a canonical form
    """

    molobj = Chem.MolFromSmiles(smiles)
    smiles = Chem.MolToSmiles(molobj, canonical=True)

    return smiles


def stoichiometry(smiles, typ="str", include_hydrogen=True):
    """


    """

    molobj = Chem.MolFromSmiles(smiles)

    if include_hydrogen:
        molobj = Chem.AddHs(molobj)

    atoms = [atom for atom in molobj.GetAtoms()]

    if typ == "str":
        atoms = [atom.GetSymbol() for atom in atoms]
        atoms = np.array(atoms)
        atomkey, atomcount = np.unique(atoms, return_counts=True)
        idxs = np.argsort(atomkey)
        rtnstoi = ""
        for idx in idxs:
            rtnstoi += atomkey[idx]
            rtnstoi += str(atomcount[idx])

    else:
        atoms = [atom.GetAtomicNum() for atom in atoms]
        # not needed
        rtnstoi = atoms

    return rtnstoi


@misc.memory.cache
def split_dict(filename, load_func):

    data = load_func(filename)

    xvalues = []
    yvalues = []

    for key in data.keys():

        value = data[key]

        if isinstance(value, list):
            value = np.mean(value)

        if isinstance(value, dict):
            value = value["K"]
            value = np.mean(value)

        smiles = key

        stoi = stoichiometry(smiles, include_hydrogen=False, typ="int")
        N = len(stoi)

        xvalues.append(N)
        yvalues.append(value)

        del stoi

    xvalues = np.array(xvalues, dtype=int)
    yvalues = np.array(yvalues)

    return xvalues, yvalues


def view_values_molecules(xvalues, yvalues, filename):
    """
    xvalues - no of atoms
    yvalues - phase transistion
    """

    print("Max atoms:", xvalues.max())
    print("Min atoms:", xvalues.min())
    print("Mean and std atoms:", xvalues.mean(), xvalues.std())

    y_mean = yvalues.mean()
    y_std = yvalues.std()

    print("Max value", yvalues.max())
    print("Min value", yvalues.min())
    print("Mean value", yvalues.mean(), yvalues.std())

    # Filter
    max_atoms = 90
    idxs, = np.where(xvalues < max_atoms)
    xvalues = xvalues[idxs]
    yvalues = yvalues[idxs]

    # Fillter outliers
    idxs, = np.where(yvalues < y_mean+y_std*4)
    xvalues = xvalues[idxs]
    yvalues = yvalues[idxs]

    n_items = xvalues.shape[0]

    print(f"Total {n_items} items for {filename}")

    views.histogram_2d_with_kde(xvalues, yvalues, filename=filename+"_overview")

    return


def main():

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--dict', action='store', help='', metavar="FILE")
    parser.add_argument('--json', action='store', help='', metavar="FILE")

    args = parser.parse_args()

    if args.dict:
        xvalues, yvalues = split_dict(args.dict, misc.load_obj)
        filename = args.dict

    if args.json:
        xvalues, yvalues = split_dict(args.json, misc.load_obj)
        filename = args.json

    view_values_molecules(xvalues, yvalues, filename)

    return

if __name__ == '__main__':
    main()
