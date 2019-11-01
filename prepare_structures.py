
import gzip

import numpy as np
import matplotlib as pyplot

import rdkit
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import rdkit.Chem.ChemicalForceFields as ChemicalForceFields
import rdkit.Chem.rdMolDescriptors as rdMolDescriptors

from chemhelp import cheminfo
import sys

import pickle

import misc


# TODO Do parallel
# TODO Save all conformers. save energy as property
# TODO systematic?
# TODO optimize
# TODO one file per molecule
#       sdf/structures.11965.sdf.gz



def generate_conformers(molobj, max_conf=100, min_conf=10):

    status = AllChem.EmbedMolecule(molobj)
    status = AllChem.UFFOptimizeMolecule(molobj)

    rot_bond = rdMolDescriptors.CalcNumRotatableBonds(molobj)

    confs = min(1 + 3*rot_bond, max_conf)
    confs = max(confs, min_conf)

    AllChem.EmbedMultipleConfs(molobj, numConfs=confs,
                useExpTorsionAnglePrefs=True,
                useBasicKnowledge=True)

    res = AllChem.MMFFOptimizeMoleculeConfs(molobj)
    res = np.array(res)

    status = res[:,0]
    energies = res[:,1]

    return energies


def get_conformations(idx, **kwargs):

    smi = Chem.MolToSmiles(molecule)
    energies = generate_conformers(molecule)

    misc.save_npy(scr + str(im) + ".energies", energies)

    txtsdf = cheminfo.molobj_to_sdfstr(molecule)

    fsdf = open(scr + str(im) + ".sdf", 'w')
    fsdf.write(txtsdf)
    fsdf.close()

    print(im, "{:} {:5.2f} {:5.2f}".format(smi, energies.mean(), energies.std()))

    return


def conformation(filename):

    scr = "_tmp_ensemble_/"

    molecules = cheminfo.read_sdffile(filename)

    for im, molecule in enumerate(molecules):

        smi = Chem.MolToSmiles(molecule)
        energies = generate_conformers(molecule)

        misc.save_npy(scr + str(im) + ".energies", energies)

        txtsdf = cheminfo.molobj_to_sdfstr(molecule)

        fsdf = open(scr + str(im) + ".sdf", 'w')
        fsdf.write(txtsdf)
        fsdf.close()

        print(im, "{:} {:5.2f} {:5.2f}".format(smi, energies.mean(), energies.std()))


    return


def main():

    db = misc.load_obj("data/melting_bradley")

    keys = db.keys()

    xaxis = []
    yaxis = []

    fullsdf = ""
    fsdf = gzip.open("data/sdf/structures.sdf.gz", 'w')
    fprop = open("data/sdf/properties.csv", 'w')

    for i, key in enumerate(keys):

        smi = key
        idx = db[key]['idx'][0] # First for
        kelvin = db[key]['K']
        kelvin = np.array(kelvin)

        # TODO histogram of sigma
        standard_deviation = np.std(kelvin)
        mean = np.mean(kelvin)

        # print(idx, mean, standard_deviation)

        # TODO Do conformation search
        # TODO Ignore very small and very big molecules

        molobj = Chem.MolFromSmiles(smi)

        n_atoms = len(molobj.GetAtoms())

        # TODO This is a choice
        if n_atoms > 50: continue

        molobj = Chem.AddHs(molobj)

        if molobj is None: continue

        # conformers = cheminfo.molobj_conformers(molobj, 10)

        molobj = cheminfo.conformationalsearch(smi)

        if molobj is None: continue

        # if len(conformers) == 0: continue

        # Find lowest conformer
        # energies = cheminfo.get_conformer_energies(molobj)
        # idx = np.argmin(energies)
        # conformer = conformers[idx]

        # sdfstr = cheminfo.molobj_to_sdfstr(molobj, conf_idx=idx)
        sdfstr = cheminfo.molobj_to_sdfstr(molobj)

        print(i, mean, standard_deviation)

        fsdf.write(sdfstr.encode())

        propstr = "{:} {:} {:}\n".format(i, mean, standard_deviation)
        fprop.write(propstr)

    fsdf.close()
    fprop.close()

    return


if __name__ == "__main__":

    conformation("data/sdf/subset_structures.sdf")


