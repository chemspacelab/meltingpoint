
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

    # TODO Define n_training

    n_trains=[2**x for x in range(4, 15)]

    misc.save_npy(scr + "n_train", n_trains)


    properties = misc.load_npy(scr + "properties")

    # TODO Load done kernels:

    # names = ["fchl19", "fp"]
    #
    #
    # for name in names:
    #     break
    #     kernel = misc.load_npy(scr + "kernel." + name)
    #     scores = cross_validation_score(kernel, properties, n_trains=n_trains)
    #
    #     print(name)
    #     misc.save_npy(scr + "score."+name, scores)
    #
    #     scores = np.around(np.mean(scores, axis=1), decimals=2)
    #     print(scores)


    # TODO Load multi kernels
    # names = ["fchl18"]
    # for name in names:
    #     break
    #     kernel = misc.load_npy(scr + "kernel." + name)[0]
    #     scores = cross_validation_score(kernel, properties, n_trains=n_trains)
    #
    #     print(name)
    #     misc.save_npy(scr + "score."+name, scores)
    #
    #     scores = np.around(np.mean(scores, axis=1), decimals=2)
    #     print(scores)


    # Load distance kernels
    # names = ["slatm", "cm", "bob"]
    # names = ["cm"]
    names = ["slatm"]
    parameters = {
        "sigma": [5000.0],
        # "sigma": [200.0],
        "lambda": [0.0]
    }

    for name in names:
        dist = misc.load_npy(scr + "dist." + name)
        kernels = get_kernels_l2distance(dist, parameters)
        kernel = next(kernels)

        scores = cross_validation_score(kernel, properties, n_trains=n_trains)

        print(name)
        misc.save_npy(scr + "score."+name, scores)

        scores = np.around(np.mean(scores, axis=1), decimals=2)
        print(scores)


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





def training_all():

    # properties
    properties = misc.load_npy(args.scratch + "properties")

    # fchls
    kernel = misc.load_npy(args.scratch + "kernel." + "fchl18")


    return


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

    # print("Generating fchl19 kernel")
    # reps = misc.load_npy(scr + "repr." + "fchl19")
    # atoms = misc.load_obj(scr + "atoms")
    # kernels = get_fchl19_kernels(reps, atoms)
    # misc.save_npy(scr + "kernel." + "fchl19", kernels)

    # print("Generating fingerprint kernel")
    # representations_fp = misc.load_obj(scr + "repr.fp")
    # kernel = get_fp_kernel(representations_fp)
    # misc.save_npy(scr + "kernel.fp", kernel)


    # Ensemble
    name = "slatm"
    representations = misc.load_npy("_tmp_ensemble_/" + "repr." + name)
    dist = generate_l2_distances(representations)
    misc.save_npy(scr + "dist." + name, dist)


    return


def plot_errors(scr):

    fig, axes = plt.subplots(1, 1, figsize=(4,4))
    ax = axes

    n_trains=[2**x for x in range(4, 4+7)]
    names = ["cm", "bob", "fchl18", "fchl19", "fp", "slatm"]

    for name in names:

        scores = misc.load_npy(scr + "score."+name)
        mean = scores.mean(axis=1)
        std = scores.std(axis=1)

        ax.errorbar(n_trains, mean, std,
            fmt='-o',
            # color="k",
            capsize=3,
            markersize=4, label=name.upper())

    ykeys = [300, 150, 75, 40]
    xkeys = n_trains

    views.learning_curve_error(ax, xkeys, ykeys,
        x_range=(10, 1100),
        y_range=(35, 350))

    leg = ax.legend(ncol=2, frameon=False)
    # leg.get_frame().set_edgecolor('b') # change color
    # leg.get_frame().set_linewidth(0.0) # remove box

    ax.set_xlabel('Training set size', fontweight='medium', fontsize=13)
    ax.set_ylabel('RMSE [K]', fontweight='medium', fontsize=13)

    plt.savefig("learning_melt.png", bbox_inches="tight")
    plt.savefig("learning_melt.pdf", bbox_inches="tight")

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

    np.random.seed(args.randomseed)

    # plot_errors(args.scratch)

    # dump_distances_and_kernels(args.scratch)

    learning_curves(args.scratch)

    return


if __name__ == '__main__':
    main()


