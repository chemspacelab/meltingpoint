
import gzip

import numpy as np
import matplotlib as pyplot

import rdkit
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import rdkit.Chem.ChemicalForceFields as ChemicalForceFields
import rdkit.Chem.rdMolDescriptors as rdMolDescriptors

from multiprocessing import Pool

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


ALLOWED_ATOMS = [
    1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 27, 35, 53
]

def is_allowed_atoms(atoms, allowed_atoms=ALLOWED_ATOMS):

    atoms = np.unique(atoms)

    for atom in atoms:
        if atom not in ALLOWED_ATOMS:
            return False

    return True


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

    try:
        result = prepare_sdf_and_csv(*line, **kwargs)
        return result
    except:
        print("procs exception made")
        return None


def prepare_sdf_and_csv(smi, values, debug=True, **kwargs):

    kelvin = np.array(values)

    # 
    standard_deviation = np.std(kelvin)
    mean = np.mean(kelvin)

    # Load molecule information
    molobj = Chem.MolFromSmiles(smi)
    atoms = cheminfo.molobj_to_atoms(molobj)
    n_atoms = len(atoms)

    # NOTE This is a choice
    # NOTE Filter organic chemistry
    if n_atoms > 50: return None
    if n_atoms < 4: return None
    if not is_allowed_atoms(atoms): return None

    molobj = Chem.AddHs(molobj)

    if molobj is None: return None

    molobj = cheminfo.conformationalsearch(smi)

    if molobj is None: return None

    # sdfstr = cheminfo.molobj_to_sdfstr(molobj)

    if debug:
        print("{:4.1f}".format(mean), "{:1.2f}".format(standard_deviation))

    # return molobj, mean, standard_deviation, values
    return molobj, values


def main(datafile, procs=0, scr="_tmp_"):

    db = misc.load_obj(datafile)

    keys = db.keys()

    print("total keys:", len(keys))

    xaxis = []
    yaxis = []

    if procs == 0:
        def get_results():

            for i, key in enumerate(keys):

                smi = key
                kelvin = db[key]
                result = prepare_sdf_and_csv(smi, kelvin)
                if result is None: continue

                yield result

        results = get_results()

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

    # Write results

    fullsdf = ""
    fsdf = gzip.open("data/sdf/structures.sdf.gz", 'w')
    fprop = open("data/sdf/properties.csv", 'w')

    for i, result in enumerate(results):

        if result is None: continue

        molobj, values = result

        sdfstr = cheminfo.molobj_to_sdfstr(molobj)
        fsdf.write(sdfstr.encode())

        valuesstr = " ".join(values)
        # propstr = "{:} {:}\n".format(mean, standard_deviation)
        propstr = f"{i} " + valuestr
        fprop.write(propstr)

    fsdf.close()
    fprop.close()

    return


def set_structures(datadict, scratch, procs=0):
    """
    take dict of smiles->value and generate sdf from smiles.
    Put in scratch/structures.sdf.gz
    Put values in scratch/properties.{txt,npy}

    """


    keys = datadict.keys()
    results = []

    # no mp
    if procs == 0:

        def get_results():
            values = []
            for key in keys:
                values.append(datadict[key])

            for smi, value in zip(keys, values):
                result = prepare_sdf_and_csv(smi, value)
                yield result

        results = get_results()

    # scale it out
    elif procs > 0:

        def workpackages():
            for i, key in enumerate(keys):

                smi = key
                kelvin = datadict[key]
                yield smi, kelvin

        lines = workpackages()

        import multiprocessing.util as util
        util.log_to_stderr(util.SUBDEBUG)

        p = Pool(procs)
        results = p.map(prepare_sdf_and_csv_procs, lines)


    print("wating for results")
    fsdf = gzip.open(scratch + "structures.sdf.gz", 'w')
    fprop = open(scratch + "properties.csv", 'w')

    for i, result in enumerate(results):

        if result is None: continue

        molobj, values = result

        mean = np.mean(values)

        prtstr = np.round(values, decimals=1)

        print("save {:4.2f}".format(mean), "-", prtstr)

        sdfstr = cheminfo.molobj_to_sdfstr(molobj)
        sdfstr += "$$$$\n"
        fsdf.write(sdfstr.encode())

        valuesstr = " ".join([str(x) for x in values])
        # propstr = "{:} {:}\n".format(mean, standard_deviation)
        propstr = f"{i} " + valuesstr + "\n"
        fprop.write(propstr)

    fsdf.close()
    fprop.close()

    return


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--datadict', action='store', help='', metavar='FILE')
    parser.add_argument('--data', action='store', help='', metavar='FILE')
    parser.add_argument('--sdf', action='store', help='', metavar='FILE')
    parser.add_argument('--scratch', action='store', help='', metavar='DIR')
    parser.add_argument('-j', '--procs', action='store', help='', type=int, metavar='int', default=0)

    args = parser.parse_args()

    if args.scratch[-1] != "/":
        args.scratch += "/"

    if args.datadict:
        data = misc.load_obj(args.datadict)
        set_structures(data, args.scratch, procs=args.procs)

    if args.data:
        main(args.data, procs=args.procs)

    if args.sdf:
        conformation("data/sdf/structures.sdf.gz", procs=args.procs)


