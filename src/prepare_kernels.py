
import os
from qml.kernels.distance import l2_distance
import qml

import time
import rdkit
import sklearn
import copy
import numpy as np
import misc
import views

import fingerprints

def generate_l2_distances(representations):
    l2dist = l2_distance(representations, representations)
    return l2dist


def get_fchl18_kernels(reps, sigmas=None, return_sigmas=False):

    if sigmas is None:
        sigmas = [0.1*2**(i) for i in range(10)]

    parameters = {
        "sigma": sigmas,
        "lambda": [0.0]
    }

    do_alchemy = "off"

    kernel_args = {
        "sigma": parameters["sigma"],
        "alchemy": do_alchemy,
    }

    kernels = qml.fchl.get_local_symmetric_kernels(reps, kernel_args=kernel_args)

    if return_sigmas:
        return sigmas, kernels

    return kernels


def get_fchl19_kernels(reps, atoms, sigmas=None, return_sigmas=False):

    if sigmas is None:
        sigmas = [2**(i) for i in range(1,10)]


    kernels = qml.kernels.get_local_symmetric_kernels(reps, atoms, sigmas)

    # kernel = qml.kernels.get_local_kernel(reps, reps, atoms, atoms, sigma)
    # kernels = qml.kernels.get_local_kernels(reps, reps, atoms, atoms, sigmas)

    if return_sigmas:
        return sigmas, kernels

    return kernels


def dump_distances_and_kernels(scr, name, procs=0):

    # TODO Properties should be read by scr!!
    # properties
    # print("Saving properties")
    # with open(scr + 'properties.csv', 'r') as f:
    #     properties = f.readlines()
    #     properties = [x.split()[0] for x in properties]
    #     properties = [float(x) for x in properties]
    #     properties = np.array(properties)

    # print(properties.shape)
    # misc.save_npy(scr + "properties", properties)

    representation_names_coordbased = ["cm", "slatm", "bob"]
    representation_names_molbased = ["morgan", "rdkitfp"]

    if procs != 0:
        os.environ["OMP_NUM_THREADS"] = str(procs)

    # Prepare fchl kernels
    if name == "fclh18":
        print("Generating fchl18 kernel")
        start = time.time()
        reps = misc.load_npy(scr + "repr." + "fchl18")
        print("shape:", reps.shape)
        sigmas, kernels = get_fchl18_kernels(reps, return_sigmas=True)
        end = time.time()
        print("time:", end-start)
        misc.save_npy(scr + "fchl18." + "sigmas", sigmas)
        misc.save_npy(scr + "kernels." + "fchl18", kernels)

        reps = None
        del reps
        kernels = None
        del kernels

    elif name == "fchl19":
        print("Generating fchl19 kernel")
        reps = misc.load_npy(scr + "repr." + "fchl19")
        print("shape:", reps.shape)
        atoms = misc.load_obj(scr + "atoms")
        start = time.time()
        sigmas, kernels = get_fchl19_kernels(reps, atoms, return_sigmas=True)
        end = time.time()
        print("time:", end-start)
        misc.save_npy(scr + "fchl19." + "sigmas", sigmas)
        misc.save_npy(scr + "kernels." + "fchl19", kernels)

    elif name in representation_names_coordbased:
        print("Distance", name)
        representations = misc.load_npy(scr + "repr." + name)
        print(representations.shape)
        dist = generate_l2_distances(representations)
        misc.save_npy(scr + "dist." + name, dist)

        dist = None
        del dist

    elif name == "rdkitfp" or name == "morgan":

        print("Generating fingerprint kernel", name)
        representations_fp = misc.load_npy(scr + "repr." + name)
        representations_fp = np.asarray(representations_fp, dtype=np.float)
        kernel = fingerprints.bitmap_jaccard_kernel(representations_fp)
        misc.save_npy(scr + "kernel." + name, kernel)

        # kernel = fingerprints.fingerprints_to_kernel(representations_fp, representations_fp, procs=procs)

    else:
        print("error: unknown representation", name)
        quit()

    return


def main():

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--scratch', action='store', help='', metavar="dir", default="_tmp_")
    parser.add_argument('-j', '--procs', action='store', help='pararallize', metavar="int", default=0, type=int)
    parser.add_argument('--get-kernels', action='store_true', help='')

    parser.add_argument('-r', '--representations', action='store', help='', metavar="STR", nargs="+")

    args = parser.parse_args()

    if args.scratch[-1] != "/":
        args.scratch += "/"

    if args.procs == -1:
        args.procs = int(os.cpu_count())
        print("set procs", args.procs)

    representation_names_coordbased = ["cm", "fchl18", "fchl19", "slatm", "bob"]
    representation_names_molbased = ["morgan", "rdkitfp"]

    if args.representations is None:
        representation_names = ["slatm", "bob", "cm", "rdkitfp", "morgan"]
    else:
        representation_names = args.representations

    for name in representation_names:
        print("goto", name)
        dump_distances_and_kernels(args.scratch, name, procs=args.procs)

    return

if __name__ == '__main__':
    main()
