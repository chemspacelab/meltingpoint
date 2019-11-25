

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




def view_values_molecules(filename):

    data = misc.load_obj(filename)


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

    # Filter
    max_atoms = 90
    idxs, = np.where(xvalues < max_atoms)

    xvalues = xvalues[idxs]
    yvalues = yvalues[idxs]


    n_items = xvalues.shape[0]

    print(f"Total {n_items} items for {filename}")

    views.histogram_2d_with_kde(xvalues, yvalues, filename=filename+"_overview")

    return



def main():

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--dict', action='store', help='', metavar="dir", default="_tmp_")

    args = parser.parse_args()

    view_values_molecules(args.dict)

    return

if __name__ == '__main__':
    main()
