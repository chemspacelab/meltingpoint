
from qml.kernels.distance import l2_distance
import qml

import rdkit
import sklearn
import copy
import numpy as np
import misc
import views

import matplotlib.pyplot as plt

def score_kernel(kernel, properties, idxs_train, idxs_test):

    kernel_train = copy.deepcopy(kernel[np.ix_(idxs_train, idxs_train)])
    kernel_test  = copy.deepcopy(kernel[np.ix_(idxs_test, idxs_train)])

    properties_train = properties[idxs_train]
    properties_test = properties[idxs_test]

    alpha = qml.math.cho_solve(kernel_train, properties_train, l2reg=1e-6)

    predictions = np.dot(kernel_test, alpha)

    diff = predictions - properties_test
    diff = diff**2
    rmse = np.sqrt(diff.mean())

    return rmse


def cross_validation_score(kernel, properties, score_func=score_kernel, **kwargs):

    fold_five = sklearn.model_selection.KFold(n_splits=5, random_state=42, shuffle=True)
    n_items = kernel.shape[0]
    X = list(range(n_items))

    scores = []

    for idxs_train, idxs_test in fold_five.split(X):

        score = cross_validated_learning_curve(kernel, properties, idxs_train, idxs_test, **kwargs)
        # score = score_func(kernel, properties, idxs_train, idxs_test)
        scores.append(score)

    scores = np.array(scores)
    scores = scores.T

    return scores


def cross_validated_learning_curve(kernel, properties, idxs_train, idxs_test,
    n_trains=[2**x for x in range(4, 15)],
    check_len=True,
    **kwargs):

    n_items = len(idxs_train)

    scores = []

    for n in n_trains:

        if n > n_items and check_len: break

        idxs = idxs_train[:n]

        score = score_kernel(kernel, properties, idxs, idxs_test, **kwargs)
        scores.append(score)

    return scores


def learning_curves(scr):

    # Define n_training
    n_trains=[2**x for x in range(4, 12)]
    misc.save_npy(scr + "n_train", n_trains)


    # Load properties
    properties = misc.load_npy(scr + "properties")

    # TODO Load done kernels:
    names = ["fchl19", "fp"]
    for name in names:
        kernel = misc.load_npy(scr + "kernel." + name)
        scores = cross_validation_score(kernel, properties, n_trains=n_trains)
        misc.save_npy(scr + "score."+name, scores)
        scores = np.around(np.mean(scores, axis=1), decimals=2)
        print(name, scores)


    # TODO Load multi kernels
    names = ["fchl18"]
    for name in names:
        kernel = misc.load_npy(scr + "kernel." + name)[0]
        scores = cross_validation_score(kernel, properties, n_trains=n_trains)
        misc.save_npy(scr + "score."+name, scores)
        scores = np.around(np.mean(scores, axis=1), decimals=2)
        print(name, scores)


    # Load distance kernels
    models = []
    parameters = {
        "name": "slatm",
        "sigma": [5000.0],
        "lambda": [0.0]
    }
    models.append(parameters)
    parameters = {
        "name": "cm",
        "sigma": [200.0],
        "lambda": [0.0]
    }
    models.append(parameters)
    parameters = {
        "name": "bob",
        "sigma": [200.0],
        "lambda": [0.0]
    }
    models.append(parameters)
    parameters = {
        "name": "avgslatm",
        "sigma": [5000.0],
        "lambda": [0.0]
    }
    models.append(parameters)

    for model in models:
        name = model["name"]
        parameters = model
        dist = misc.load_npy(scr + "dist." + name)
        kernels = get_kernels_l2distance(dist, parameters)
        kernel = next(kernels)
        scores = cross_validation_score(kernel, properties, n_trains=n_trains)
        misc.save_npy(scr + "score."+name, scores)
        scores = np.around(np.mean(scores, axis=1), decimals=2)
        print(name, scores)


    return



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


def get_fchl18_kernels(reps, sigmas=None, return_sigmas=False):

    if sigmas is None:
        sigmas = [0.1*2**(i) for i in range(10)]

    parameters = {
        "sigma": sigmas,
        "lambda": [0.0]
    }

    kernel_args = {
        "sigma": parameters["sigma"]
    }

    kernels = qml.fchl.get_local_symmetric_kernels(reps, kernel_args=kernel_args)

    if return_sigmas:
        sigmas, kernels

    return kernels


def get_fchl19_kernels(reps, atoms, sigmas=None, return_sigmas=False):

    if sigmas is None:
        sigmas = [2**(i) for i in range(1,10)]

    # kernel = qml.kernels.get_local_kernel(reps, reps, atoms, atoms, sigma)
    kernels = qml.kernels.get_local_kernel_gaussian(reps, reps, atoms, atoms, sigmas)

    if return_sigmas:
        return sigmas, kernels

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


def training_all():

    # properties
    properties = misc.load_npy(args.scratch + "properties")

    # fchls
    kernel = misc.load_npy(args.scratch + "kernel." + "fchl18")


    return


def dump_distances_and_kernels(scr):

    # properties
    print("Saving properties")
    with open('data/sdf/subset_properties.csv', 'r') as f:
        properties = f.readlines()
        properties = [float(x) for x in properties]
        properties = np.array(properties)

    misc.save_npy(scr + "properties", properties)

    # Prepare distances
    representation_names = ["cm", "bob", "slatm", "avgslatm"]
    for name in representation_names:
        break
        print("Distance", name)
        representations = misc.load_npy(scr + "repr." + name)
        print(representations.shape)
        dist = generate_l2_distances(representations)
        misc.save_npy(scr + "dist." + name, dist)

    # Prepare fchl kernels
    print("Generating fchl18 kernel")
    reps = misc.load_npy(scr + "repr." + "fchl18")
    sigmas, kernels = get_fchl18_kernels(reps, return_sigmas=True)
    misc.save_npy(scr + "fchl18." + "sigmas", sigmas)
    misc.save_npy(scr + "kernels." + "fchl18", kernels)

    print("Generating fchl19 kernel")
    reps = misc.load_npy(scr + "repr." + "fchl19")
    atoms = misc.load_obj(scr + "atoms")
    sigmas, kernels = get_fchl19_kernels(reps, atoms, return_sigmas=True)
    misc.save_npy(scr + "fchl19." + "sigmas", sigmas)
    misc.save_npy(scr + "kernels." + "fchl19", kernels)

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
    parser.add_argument('--plot', action='store_true')
    parser.add_argument('-j', '--cpu', action='store', help='pararallize', metavar="int", default=0)
    parser.add_argument('--get-kernels', action='store_true', help='')
    parser.add_argument('--get-learning-curves', action='store_true', help='')

    args = parser.parse_args()

    if args.scratch[-1] != "/":
        args.scratch += "/"

    np.random.seed(args.randomseed)

    plot_errors(args.scratch)

    if args.get_kernels:
        dump_distances_and_kernels(args.scratch)

    if args.get_learning_curves:
        learning_curves(args.scratch)

    return


if __name__ == '__main__':
    main()


