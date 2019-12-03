
import matplotlib.pyplot as plt

import gzip
import misc
import numpy as np
from chemhelp import cheminfo

from scipy.spatial import ConvexHull, distance

import rdkit.Chem as Chem


def parse_molobj(molobj):

    # TODO Prepare for pool


    return molobj, value


def parse_ochem(filename, procs=0, debug=False):

    # molobjs = cheminfo.read_sdffile("_tmp_ochem/meltingpoints_weight100_200.sdf.gz")
    molobjs = cheminfo.read_sdffile(filename)

    properties = []

    for molobj in molobjs:

        if molobj is None: continue

        mol_smi = cheminfo.molobj_to_smiles(molobj)
        props = molobj.GetPropsAsDict()
        keys = props.keys()

        if "SMILES" not in keys:
            continue

        prop_smiles = props["SMILES"]

        # Ignore multi molecules
        if "." in prop_smiles:
            if debug:
                print(f"ignore: {prop_smiles}")
            continue

        # Count
        atoms = cheminfo.molobj_to_atoms(molobj)

        if len(atoms) < 3:
            if debug:
                print("ignore small", props)
            continue

        if len(atoms) > 40:
            if debug:
                print("ignore large", props)
            continue

        atoms_carbons, = np.where(atoms == 6)
        if len(atoms_carbons) < 1:
            if debug:
                print("ignore non-org", props)
            continue

        # Add hydrogens and optimize structure
        molobj = cheminfo.molobj_add_hydrogens(molobj)
        status = cheminfo.molobj_optimize(molobj)

        # if unconverged
        if status != 0:

            # try the smiles
            molobj, status = cheminfo.smiles_to_molobj(prop_smiles)
            if molobj is None:
                print("error", props)
                continue

            molobj = cheminfo.molobj_add_hydrogens(molobj)
            status = cheminfo.molobj_optimize(molobj)

            if status != 0:
                print("error", props)
                continue


        idx_value = [key for key in keys if "measured, converted" in key]
        idx_value = idx_value[0]

        idx_unit = [key for key in keys if "UNIT" in key]
        idx_unit = [key for key in idx_unit if "Point" in key]
        idx_unit = idx_unit[0]

        prop_unit = props[idx_unit]
        prop_value = props[idx_value]

        if prop_unit == "Celsius":
            prop_value += 273.15
        elif prop_unit == "K":
            pass
        else:
            print("error unknown unit", prop_unit, props)
            continue

        properties.append(prop_value)

    return molobjs, properties


def main():

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--scratch', action='store', help='', metavar="DIR", default="_tmp_")
    parser.add_argument('--sdf', action='store', help='', metavar="FILE")
    parser.add_argument('-j', '--procs', action='store', help='pararallize', metavar="int", default=0)

    args = parser.parse_args()

    if args.scratch[-1] != "/":
        args.scratch += "/"

    molobjs, values = parse_ochem(args.sdf, procs=0, debug=True)

    print("wating for results")
    fsdf = gzip.open(args.scratch + "structures.sdf.gz", 'w')
    fprop = open(args.scratch + "properties.csv", 'w')

    for molobj, value in zip(molobjs, values):

        sdfstr = cheminfo.molobj_to_sdfstr(molobj)
        sdfstr += "$$$$\n"

        propstr = "{:} {:}\n".format(value, 0.0)
        fprop.write(propstr)

    fsdf.close()
    fprop.close()

    return


if __name__ == '__main__':
    main()
