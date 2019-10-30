
from qml.kernels.distance import l2_distance
import qml

import rdkit
import sklearn
import copy
import numpy as np
import misc
import views


def score_kernel(kernel, properties, idxs_train, idxs_test):

    kernel_train = copy.deepcopy(kernel[np.ix_(idxs_train, idxs_train)])
    kernel_test  = copy.deepcopy(kernel[np.ix_(idxs_test, idxs_train)])

    properties_train = properties[idxs_train]
    properties_test = properties[idxs_train]

    alpha = math.cho_solve(kernel_train, properties_train)

    predictions = np.dot(kernel_test, alpha)

    diff = predictions - properties_test
    diff = diff**2
    rmse = np.sqrt(diff.mean())

    return rmse


def cross_validation_score(kernel, properties, score_func=score_kernel):

    fold_five = sklearn.model_selection.KFold(n_splits=5, random_state=42, shuffle=True)
    n_items = kernel.shape[0]
    X = list(range(n_items))

    scores = []

    for idxs_train, idxs_test in fold_five.split(X):

        score = score_func(kernel, properties, idxs_train, idxs_test)
        scores.append(score)

        quit()

    return scores


def generate_l2_distances(representations):
    l2dist = l2_distance(representations, representations)
    return l2dist


def get_kernels_l2distance(l2distance, parameters):

    kernel_copy = copy.deepcopy(l2distance)
    kernel_copy = np.square(kernel_copy)
    kernel_copy *= -1

    for k_sigma in parameters["sigma"]:

        kernel = copy.deepcopy(kernel_copy)
        kernel /= (2*k_sigma**2)
        kernel = np.exp(kernel)

        diag_kernel = kernel[np.diag_indices_from(kernel)]

        for k_lambda in parameters["lambda"]:
            kernel[np.diag_indices_from(kernel)] = diag_kernel + k_lambda

        # TODO
        yield kernel


def get_fchl18_kernels(reps):

    parameters = {
        "sigma": [1.0, 2.5, 5, 10],
        "lambda": [float(10)**-8]
    }

    kernel_args = {
        "sigma": parameters["sigma"]
    }

    kernels = qml.fchl.get_local_symmetric_kernels(reps, kernel_args=kernel_args)

    return kernels


def get_fchl19_kernels(reps, atoms):

    sigma = 2.0

    kernels = qml.kernels.get_local_kernel(reps, reps, atoms, atoms, sigma)

    return kernels


def get_fp_kernel(reps):

    n_items = len(reps)

    kernel = np.zeros((n_items, n_items))

    for i, repi in enumerate(reps):
        for j, repj in enumerate(reps):
            if i > j: continue
            kernel[i,j] = rdkit.DataStructs.FingerprintSimilarity(repi, repj)
            kernel[j,i] = kernel[i,j]

    return kernel


def dump_distances_and_kernels(scr):

    # properties
    # print("Saving properties")
    # misc.save_npy(scr + "properties", properties)
    # with open('data/sdf/subset_properties.csv', 'r') as f:
    #     properties = f.readlines()
    #     properties = [float(x) for x in properties]
    #     properties = np.array(properties)

    # Prepare distances
    # representation_names = ["cm", "bob", "slatm"]
    # for name in representation_names:
    #     print("Distance", name)
    #     representations = misc.load_npy(args.scratch + "repr." + name)
    #     dist = generate_l2_distances(representations)
    #     misc.save_npy(scr + "dist." + name, dist)

    # Prepare fchl kernels
    # print("Generating fchl18 kernel")
    # reps = misc.load_npy(scr + "repr." + "fchl18")
    # kernels = get_fchl18_kernels(reps)
    # misc.save_npy(scr + "kernel." + "fchl18", kernels)

    print("Generating fchl19 kernel")
    reps = misc.load_npy(scr + "repr." + "fchl19")
    atoms = misc.load_obj(scr + "atoms")
    kernels = get_fchl19_kernels(reps, atoms)
    misc.save_npy(scr + "kernel." + "fchl19", kernels)

    # print("Generating fingerprint kernel")
    # representations_fp = misc.load_obj(scr + "repr.fp")
    # kernel = get_fp_kernel(representations_fp)
    # misc.save_npy(scr + "kernel.fp", kernel)


    return


def main():

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--scratch', action='store', help='', metavar="dir", default="_tmp_")
    parser.add_argument('--randomseed', action='store', help='random seed', metavar="int", default=1)
    parser.add_argument('-j', '--cpu', action='store', help='pararallize', metavar="int", default=0)

    args = parser.parse_args()

    if args.scratch[-1] != "/":
        args.scratch += "/"


    dump_distances_and_kernels(args.scratch)


    return


if __name__ == '__main__':
    main()


