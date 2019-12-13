import os
import misc
import numpy as np
import qml
from chemhelp import cheminfo
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem


from functools import partial
from multiprocessing import Process, Pool
import multiprocessing.managers


class MyManager(multiprocessing.managers.BaseManager):
    pass
MyManager.register('np_zeros', np.zeros, multiprocessing.managers.ArrayProxy)



def procs_representation_slatm(args, **kwargs):

    coord = args[0]
    atoms = args[1]
    mbtypes = kwargs['mbtypes']

    rep = qml.representations.generate_slatm(coord, atoms, mbtypes)

    return rep


def get_representations_slatm(atoms,
    structures, scr="_tmp_/",
    mbtypes=None,
    debug=True,
    procs=0,
    **kwargs):
    """
    atoms -- list of molecule atoms

    """

    if mbtypes is None:

        filename_mbtypes = scr + "slatm.mbtypes"

        try: mbtypes = misc.load_obj(filename_mbtypes)
        except FileNotFoundError:

            print("Generate slatm mbtypes")
            mbtypes = qml.representations.get_slatm_mbtypes(atoms)
            misc.save_obj(filename_mbtypes, mbtypes)


    if debug:
        print("Generate slatm representations")

    replist = []

    # Set OMP
    if procs > 1:
        os.environ["OMP_NUM_THREADS"] = "1"

        workargs = zip(structures, atoms)
        workargs = list(workargs)

        pool = Pool(processes=procs)
        funcname = partial(procs_representation_slatm, mbtypes=mbtypes)
        replist = pool.map(funcname, workargs)

    else:
        for i, (coord, atom) in enumerate(zip(structures, atoms)):
            rep = qml.representations.generate_slatm(coord, atom, mbtypes)
            replist.append(rep)

    # replist = [qml.representations.generate_slatm(coordinate, atom, mbtypes) for coordinate, atom in zip(structures, atoms)]
    replist = np.array(replist)

    # for i, rep in enumerate(replist):
    #     m = rep.mean()
    #     if np.isnan(m):
    #         print(i, rep.mean())
    # print(replist.mean())

    return replist


def get_representations_fchl(atoms, structures, max_atoms=35, cut_distance=10**6, **kwargs):

    print("Generate fchl18 representations")
    replist = [qml.fchl.generate_representation(coordinate, atom, max_size=max_atoms, cut_distance=cut_distance) for coordinate, atom in zip(structures, atoms)]
    replist = np.array(replist)

    return replist


def get_representations_fchl19(atoms, structures, max_atoms=35, cut_distances=8.0, **kwargs):

    print("Generate fchl19 representations")

    elements = []
    for atom in atoms:
        elements += list(np.unique(atom))
        elements = list(np.unique(elements))

    parameters = {
        "pad": max_atoms,
	'nRs2': 22,
	'nRs3': 17,
	'eta2': 0.41,
	'eta3': 0.97,
	'three_body_weight': 45.83,
	'three_body_decay': 2.39,
	'two_body_decay': 2.39,
        "rcut": cut_distances,
        "acut": cut_distances,
        "elements": elements
    }

    replist = [qml.representations.generate_fchl_acsf(atom, coord, **parameters) for atom,coord in zip(atoms, structures)]
    replist = np.array(replist)

    return replist


def get_representations_cm(atoms, structures, max_atoms=23, **kwargs):

    print("Generate coulomb matrix representations")
    replist = [qml.representations.generate_coulomb_matrix(nuclear_charges, coordinates, size=max_atoms) for coordinates, nuclear_charges in zip(structures, atoms)]
    replist = np.array(replist)

    return replist


def get_representations_bob(atoms, structures, max_atoms=23, asize=None, **kwargs):

    if asize is None:
        print("Generate atypes")
        asize = {}

        for charges in atoms:
            uc, count = np.unique(charges, return_counts=True)

            for atom, N, in zip(uc, count):

                atom = cheminfo.atom_str(atom)

                if atom not in asize:
                    asize[atom] = N
                else:
                    if asize[atom] < N:
                        asize[atom] = N

    print("Generate bag-of-bond representations")

    replist = [qml.representations.generate_bob(charges, coordinates, [], size=max_atoms, asize=asize) for charges, coordinates in zip(atoms, structures)]
    replist = np.asarray(replist)

    return replist


def molobjs_to_fingerprints(molobjs, procs=0):

    # TODO do this in parallel

    fps = []

    for molobj in molobjs:
        fp1 = Chem.RDKFingerprint(molobj)
        fps.append(fp1)

    return fps


def molobjs_to_xyzs(molobjs):

    mol_atoms = []
    mol_coord = []

    for molobj in molobjs:
        atoms, coord = cheminfo.molobj_to_xyz(molobj)
        mol_atoms.append(atoms)
        mol_coord.append(coord)

    return mol_atoms, mol_coord


def xyzs_to_representations(mol_atoms, mol_coord, name="cm", **kwargs):

    if name == "slatm":

        reprs = get_representations_slatm(mol_atoms, mol_coord, **kwargs)

    elif name == "fchl" or name == "fchl18":

        reprs = get_representations_fchl(mol_atoms, mol_coord, **kwargs)

    elif name == "fchl19":

        reprs = get_representations_fchl19(mol_atoms, mol_coord, **kwargs)

    elif name == "cm":

        reprs = get_representations_cm(mol_atoms, mol_coord, **kwargs)

    elif name == "bob":

        reprs = get_representations_bob(mol_atoms, mol_coord, **kwargs)

    else:
        quit("error representation unknown:", name)

    return reprs


def print_test(idx, **kwargs):

    print(idx, "wat")

    return 0, np.array([1.0])


def get_avg_repr(idx, scr="_tmp_ensemble_/", **kwargs):

    name = "slatm"

    energies = misc.load_npy(scr + str(idx) + ".energies")
    molobjs = cheminfo.read_sdffile(scr + str(idx) + ".sdf")
    molobjs = [mol for mol in molobjs]

    xyzs = molobjs_to_xyzs(molobjs)
    reprs = xyzs_to_representations(*xyzs, **kwargs)

    # Boltzmann factors
    factors = np.exp(-energies)
    factors /= np.sum(factors)

    length = reprs.shape[1]
    avgrep = np.zeros(length)

    for rep, factor in zip(reprs, factors):
        avgrep += factor*rep

    print(idx, avgrep.shape)

    if "array" in kwargs:
        results = kwargs["array"]
        results[idx,:] = avgrep

    else:
        return idx, avgrep


def generate_conformer_representation(scr="_tmp_ensemble_/", procs=0):

    names = ["cm", "slatm", "bob"]
    name = "slatm"

    mbtypes = misc.load_npy(scr +"slatm.mbtypes")

    # TODO Calculate max_size
    mol_atoms = misc.load_obj(scr + "atoms")
    max_atoms = [len(atoms) for atoms in mol_atoms]
    max_atoms = max(max_atoms)

    kwargs = {
        "name": name,
        "mbtypes": mbtypes,
        "debug": False,
        "max_atoms": max_atoms,
    }

    # n_total = 1285
    n_total = 3456
    idxs = range(n_total)

    avgreps = [0]*n_total

    if procs == 0:

        for idx in idxs:

            idx, avgrep = get_avg_repr(idx, **kwargs)
            avgreps[idx] = avgrep

    else:

        idx, rep = get_avg_repr(0, **kwargs)
        rep_size = rep.shape[0]
        print("rep size", rep_size)

        m = MyManager()
        m.start()

        results = m.np_zeros((n_total, rep_size))

        pool = Pool(32)

        kwargs["array"] = results
        func = partial(get_avg_repr, **kwargs)
        pool.map(func, idxs)
        avgreps = results

        # results = misc.parallel(idxs, get_avg_repr, [], kwargs, procs=nprocs)
        #
        # for result in results:
        #     idx, avgrep = result
        #     avgreps[idx] = avgrep
        #     print(idx, avgrep.mean())

    avgreps = np.array(avgreps)
    misc.save_npy(scr + "repr.avgslatm", avgreps)

    return


def main():

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--scratch', action='store', help='', metavar="dir", default="_tmp_")
    parser.add_argument('--conformers', action='store_true', help='')
    parser.add_argument('--sdf', action='store', help='', metavar="file")
    parser.add_argument('-j', '--procs', action='store', help='pararallize', metavar="int", default=0, type=int)

    args = parser.parse_args()

    if args.scratch[-1] != "/":
        args.scratch += "/"

    if not args.conformers:
        molobjs = cheminfo.read_sdffile(args.sdf)
        molobjs = [mol for mol in molobjs]
        # molobjs = molobjs[:10]

        # representation_names = ["cm", "fchl18", "fchl19", "slatm", "bob"]
        # representation_names = ["fchl18"]
        # representation_names = ["bob"]
        representation_names = ["slatm", "bob", "cm"]

        xyzs = molobjs_to_xyzs(molobjs)

        mol_atoms, mol_coords = xyzs
        misc.save_obj(args.scratch + "atoms", mol_atoms)

        # Print unique atoms
        unique_atoms = []
        for atoms in mol_atoms:
            unique_atoms += list(np.unique(atoms))

        unique_atoms = np.array(unique_atoms)
        unique_atoms = unique_atoms.flatten()
        unique_atoms = np.unique(unique_atoms)

        # Calculate max_size
        max_atoms = [len(atoms) for atoms in mol_atoms]
        max_atoms = max(max_atoms)

        print("total mols:", len(mol_coords))
        print("atom types:", unique_atoms)
        print("max atoms: ", max_atoms)

        quit()

        # Gas phase
        for name in representation_names:
            representations = xyzs_to_representations(*xyzs, name=name, scr=args.scratch, max_atoms=max_atoms, procs=args.procs)
            print(representations.shape)
            misc.save_npy(args.scratch + "repr." + name, representations)
            del representations

        # # fingerprints
        print("Generate fingerprints")
        fingerprints = molobjs_to_fingerprints(molobjs)
        misc.save_obj(args.scratch + "repr." + "fp", fingerprints)


    # Ensemble
    if args.conformers:
        generate_conformer_representation(scr=args.scratch, procs=args.procs)

    return


if __name__ == '__main__':
    main()
