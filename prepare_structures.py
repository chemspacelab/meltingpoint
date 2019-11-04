
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
import time

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


def get_conformations(line, scr="_tmp_ensemble_/", **kwargs):

    im, molecule = line

    # smi = Chem.MolToSmiles(molecule)
    energies = generate_conformers(molecule)

    misc.save_npy(scr + str(im) + ".energies", energies)

    txtsdf = cheminfo.molobj_to_sdfstr(molecule)

    fsdf = open(scr + str(im) + ".sdf", 'w')
    fsdf.write(txtsdf)
    fsdf.close()

    print(im, "{:} {:5.2f} {:5.2f}".format("smi", energies.mean(), energies.std()))

    return


def conformation(filename, procs=0):

    scr = "_tmp_ensemble_/"

    molecules = cheminfo.read_sdffile(filename)

    if procs == 0:
        for im, molecule in enumerate(molecules):
            get_conformations((im,molecule), scr=scr)

    else:
        def workpackages():
            for im, molecule in enumerate(molecules):
                yield im, molecule

        lines = workpackages()

        results = misc.parallel(lines, get_conformations, [], {"scr":scr}, procs=procs)
        for result in results:
            pass

        # misc.parallel(lines)#, get_conformations, [], {"scr":scr}, procs=procs)

    return


def prepare_sdf_and_csv_procs(line, **kwargs):

    result = prepare_sdf_and_csv(*line)

    return result

def prepare_sdf_and_csv(smi, values, **kwargs):

    kelvin = np.array(values)

    # TODO histogram of sigma
    standard_deviation = np.std(kelvin)
    mean = np.mean(kelvin)

    molobj = Chem.MolFromSmiles(smi)

    n_atoms = len(molobj.GetAtoms())

    # NOTE This is a choice
    if n_atoms > 50: return None
    if n_atoms < 4: return None

    molobj = Chem.AddHs(molobj)

    if molobj is None: return None

    molobj = cheminfo.conformationalsearch(smi)

    if molobj is None: return None

    sdfstr = cheminfo.molobj_to_sdfstr(molobj)

    # fsdf.write(sdfstr.encode())
    #
    # propstr = "{:} {:} {:}\n".format(i, mean, standard_deviation)
    # fprop.write(propstr)
    #
    # other_end = time.time()
    print("{:4.1f}".format(mean), "{:1.2f}".format(standard_deviation))

    return molobj, mean, standard_deviation


def main(datafile, procs=0):

    db = misc.load_obj(datafile)

    keys = db.keys()

    print("total keys:", len(keys))

    xaxis = []
    yaxis = []

    fullsdf = ""
    fsdf = gzip.open("data/sdf/structures.sdf.gz", 'w')
    fprop = open("data/sdf/properties.csv", 'w')


    # TODO Do parallel

    if procs == 0:
        for i, key in enumerate(keys):

            smi = key
            kelvin = db[key]

            results = prepare_sdf_and_csv(smi, kelvin)

            if results is None: continue

            molobj, mean, standard_deviation = results

            sdfstr = cheminfo.molobj_to_sdfstr(molobj)
            fsdf.write(sdfstr.encode())

            propstr = "{:} {:} {:}\n".format(i, mean, standard_deviation)
            fprop.write(propstr)


    else:

        def workpackages():
            for i, key in enumerate(keys):

                # if i > 5000: break

                smi = key
                kelvin = db[key]
                yield smi, kelvin

        lines = workpackages()

        results = misc.parallel(lines, prepare_sdf_and_csv_procs, [], {}, procs=procs)

        print("streaming results")

        for i, result in enumerate(results):

            if result is None: continue

            molobj, mean, standard_deviation = result

            sdfstr = cheminfo.molobj_to_sdfstr(molobj)
            fsdf.write(sdfstr.encode())

            propstr = "{:} {:}\n".format(mean, standard_deviation)
            fprop.write(propstr)


    fsdf.close()
    fprop.close()

    return



if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', action='store', help='', metavar='FILE')
    parser.add_argument('--sdf', action='store', help='', metavar='FILE')
    parser.add_argument('-j', '--procs', action='store', help='', type=int, metavar='int', default=0)

    args = parser.parse_args()

    if args.data:
        main(args.data, procs=args.procs)

    if args.sdf:
        conformation("data/sdf/structures.sdf.gz", procs=args.procs)


