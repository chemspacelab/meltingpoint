
import numpy as np

from chemhelp import cheminfo

import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem

import qml
from qml.kernels import get_atomic_local_kernel, get_local_kernel
from qml.kernels import gaussian_kernel, laplacian_kernel
from qml.kernels import kpca
from qml.math import svd_solve
from qml.representations import generate_fchl_acsf

import matplotlib.pyplot as plt
import misc


def search_similar_molecules(smiles_list, refsmi=None):
    """

    https://arxiv.org/abs/1903.11372
    https://jcheminf.biomedcentral.com/articles/10.1186/s13321-015-0069-3

    """

    fp1 = Chem.RDKFingerprint(mol1)
    fp2 = Chem.RDKFingerprint(mol2)
    # print "RDK fingerprint: ",DataStructs.TanimotoSimilarity(fp1,fp2)


    return


def overview_properties_pca():

    elements = []

    with open('data/sdf/subset_properties.csv', 'r') as f:
        properties = f.readlines()
        properties = [float(x) for x in properties]
        properties = np.array(properties)

    representations = []
    molobjs = cheminfo.read_sdffile("data/sdf/subset_structures.sdf")

    mols_atoms = []
    mols_coord = []

    n_atoms = 0
    n_items = 500

    for i, molobj in enumerate(molobjs):

        atoms, coord = cheminfo.molobj_to_xyz(molobj)

        mols_atoms.append(atoms)
        mols_coord.append(coord)

        elements += list(np.unique(atoms))
        elements = list(np.unique(elements))

        if len(atoms) > n_atoms:
            n_atoms = len(atoms)

        i += 1
        if i == n_items:
            break

    properties = properties[:n_items]

    print(elements)
    print(n_atoms)
    print(len(mols_atoms))

    distance_cut = 20.0
    parameters = {
        "pad": n_atoms,
	'nRs2': 22,
	'nRs3': 17,
	'eta2': 0.41,
	'eta3': 0.97,
	'three_body_weight': 45.83,
	'three_body_decay': 2.39,
	'two_body_decay': 2.39,
        "rcut": distance_cut,
        "acut": distance_cut,
        "elements": elements
    }

    for atoms, coord in zip(mols_atoms, mols_coord):
        representation = generate_fchl_acsf(atoms, coord, **parameters)
        representations.append(representation)

    representations = np.array(representations)

    sigma = 10.

    kernel = qml.kernels.get_local_kernel(representations, representations, mols_atoms, mols_atoms, sigma)

    print(kernel.shape)

    pca = kpca(kernel, n=2)

    fig, axs = plt.subplots(2, 1, figsize=(5,10))
    sc = axs[0].scatter(*pca, c=properties)
    fig.colorbar(sc, ax=axs[0])
    im = axs[1].imshow(kernel)
    fig.colorbar(im, ax=axs[1])
    fig.savefig("_tmp_pca_prop.png")

    return



def search_molcules(mollist, proplist, conf_scr="_tmp_ensemble_"):

    sublist_mol = []
    sublist_prop = []
    sublist_idxs = []

    for idx, (molobj, prop) in enumerate(zip(mollist, proplist)):

        atoms = molobj.GetAtoms()
        atoms = [atom.GetSymbol() for atom in atoms]

        atoms = np.array(atoms)
        uatm, counts = np.unique(atoms, return_counts=True)

        if 'C' not in uatm:
            continue

        c_idx, = np.where(uatm == 'C')
        c_idx = c_idx[0]

        if counts[c_idx] > 8:
            continue

        if counts[c_idx] < 3:
            continue

        h_idx, = np.where(uatm == 'C')
        h_idx = h_idx[0]

        counts[h_idx] = 0

        N = sum(counts)

        if N > 10:
            continue

        smi = cheminfo.molobj_to_smiles(molobj, remove_hs=True)

        try:

            value, stddev = prop.strip().split()

            value = float(value)
            stddev = float(stddev)

        except:

            value = prop.strip()
            value = float(value)

        sublist_mol.append(molobj)
        sublist_prop.append(value)
        sublist_idxs.append(idx)

    return sublist_mol, sublist_prop, sublist_idxs


def main():

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--scratch', action='store', help='', metavar="dir", default="tmp2/")
    parser.add_argument('--randomseed', action='store', help='random seed', metavar="int", default=666)
    parser.add_argument('--sdf', action='store', help='', metavar="file")
    parser.add_argument('-j', '--cpu', action='store', help='pararallize', metavar="int", default=0)

    args = parser.parse_args()


    molecules = cheminfo.read_sdffile('data/sdf/structures.sdf.gz')
    properties = open('data/sdf/properties.csv', 'r')

    sub_mol, sub_prop, idxs = search_molcules(molecules, properties)

    properties.close()


    fm = open('data/sdf/subset_structures.sdf', 'w')
    fp = open('data/sdf/subset_properties.csv', 'w')

    for mol, prop in zip(sub_mol, sub_prop):

        sdf = cheminfo.molobj_to_sdfstr(mol)
        fm.write(sdf)
        fm.write("$$$$\n")
        fp.write(str(prop) + "\n")

    fm.close()
    fp.close()

    for i, idx in enumerate(idxs):
        from_dir = "_tmp_ensemble_/"
        to_dir = "_tmp_subset_/conformers/"
        cmd = "cp {:}{:}.sdf {:}{:}.sdf".format(from_dir, str(idx), to_dir, str(i))
        list(misc.shell(cmd))
        cmd = "cp {:}{:}.energies.npy {:}{:}.energies.npy".format(from_dir, str(idx), to_dir, str(i))
        list(misc.shell(cmd))

        print(cmd)

    return

if __name__ == "__main__":
    # overview_properties_pca()
    main()

