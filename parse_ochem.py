
import gzip

import multiprocessing as mp
import matplotlib.pyplot as plt
import numpy as np
import rdkit.Chem as Chem
from scipy.spatial import ConvexHull, distance

import misc
from chemhelp import cheminfo


def parse_molobj(molobj, debug=False, **kwargs):

    if molobj is None:
        return None, None

    mol_smi = cheminfo.molobj_to_smiles(molobj)
    props = molobj.GetPropsAsDict()
    keys = props.keys()

    result = parse_molandprop(molobj, props)

    return result


def parse_molandprop(*args, debug=False, **kwargs):

    if len(args) > 1:
        molobj = args[0]
        props = args[1]
    else:
        molobj, props = args[0]

    if molobj is None:
        return None, None

    keys = props.keys()

    if "SMILES" not in keys:
        return None, None

    prop_smiles = props["SMILES"]

    # Ignore multi molecules
    if "." in prop_smiles:
        if debug:
            print(f"ignore: {prop_smiles}")
        return None, None

    # Count
    atoms = cheminfo.molobj_to_atoms(molobj)

    if len(atoms) < 3:
        if debug:
            print("ignore small", props)
        return None, None

    # if len(atoms) > 40:
    #     if debug:
    #         print("ignore large", props)
    #     return None, None

    atoms_carbons, = np.where(atoms == 6)
    if len(atoms_carbons) < 1:
        if debug:
            print("ignore non-org", props)
        return None, None

    # Add hydrogens and optimize structure
    molobj = cheminfo.molobj_add_hydrogens(molobj)
    status = cheminfo.molobj_optimize(molobj)

    # if unconverged
    if status != 0:

        # try the smiles
        molobj, status = cheminfo.smiles_to_molobj(prop_smiles)
        if molobj is None:
            print("error", props)
            return None, None

        molobj = cheminfo.molobj_add_hydrogens(molobj)
        status = cheminfo.molobj_optimize(molobj)

        if status != 0:
            print("error", props)
            return None, None


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
        return None, None

    return molobj, prop_value


def parse_ochem(filename, procs=0, debug=False):

    molobjs = cheminfo.read_sdffile(filename)
    # molobjs = list(molobjs)

    success_properties = []
    success_molobjs = []

    if procs > 0:

        def generate_input(molobjs):
            for molobj in molobjs:
                if molobj is None: continue
                props = molobj.GetPropsAsDict()
                yield molobj, props

        pool = mp.Pool(processes=procs)
        results = pool.map(parse_molandprop, generate_input(molobjs))

        for result in results:
            molobj, value = result
            if molobj is None:
                continue

            success_properties.append(value)
            success_molobjs.append(molobj)

    else:

        for molobj in molobjs:
            molobj, value = parse_molobj(molobj, debug=debug)

            if molobj is None:
                continue

            success_properties.append(value)
            success_molobjs.append(molobj)

    return success_molobjs, success_properties


def main():

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--scratch', action='store', help='', metavar="DIR", default="_tmp_")
    parser.add_argument('--sdf', action='store', help='', metavar="FILE", nargs="+")
    parser.add_argument('-j', '--procs', action='store', help='pararallize', metavar="int", default=0, type=int)

    args = parser.parse_args()

    if args.scratch[-1] != "/":
        args.scratch += "/"

    # fsdf = gzip.open(args.scratch + "structures.sdf.gz", 'w')
    # fprop = open(args.scratch + "properties.csv", 'w')
    mol_val_dict = {}

    for sdf in args.sdf:

        print("reading", sdf)

        molobjs, values = parse_ochem(sdf, debug=True, procs=args.procs)

        for molobj, value in zip(molobjs, values):

            smiles = cheminfo.molobj_to_smiles(molobj, remove_hs=True)

            if "smiles" not in mol_val_dict:
                mol_val_dict[smiles] = []
            else:
                print("duplicate", smiles)

            mol_val_dict[smiles].append(value)

            # sdfstr = cheminfo.molobj_to_sdfstr(molobj)
            # sdfstr += "$$$$\n"
            #
            # propstr = "{:} {:}\n".format(value, 0.0)
            # fprop.write(propstr)

    # fsdf.close()
    # fprop.close()

    misc.save_json(args.scratch + "molecule_data", mol_val_dict)
    misc.save_obj(args.scratch + "molecule_data", mol_val_dict)

    return


if __name__ == '__main__':
    main()
