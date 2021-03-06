
from qml.kernels.distance import l2_distance
import qml

import gc
import time
import rdkit
import sklearn
import copy
import numpy as np
import misc
import views
from tqdm import tqdm

import matplotlib.pyplot as plt
import functools

def score_kernels(kernels, properties, idxs_train, idxs_test, l2reg=1e-6):


    return


def score_kernel(kernel, properties, idxs_train, idxs_test, l2reg=1e-6):

    kernel_train = copy.deepcopy(kernel[np.ix_(idxs_train, idxs_train)])
    kernel_test  = copy.deepcopy(kernel[np.ix_(idxs_test, idxs_train)])

    properties_train = properties[idxs_train]
    properties_test = properties[idxs_test]

    alpha = qml.math.cho_solve(kernel_train, properties_train, l2reg=l2reg)

    predictions = np.dot(kernel_test, alpha)

    diff = predictions - properties_test
    diff = diff**2
    rmse = np.sqrt(diff.mean())

    return rmse


def get_predictions(kernel, properties, idxs_train, idxs_test, l2reg=0.0,
    rtn_train=False):

    kernel_train = copy.deepcopy(kernel[np.ix_(idxs_train, idxs_train)])
    kernel_test  = copy.deepcopy(kernel[np.ix_(idxs_test, idxs_train)])

    properties_train = properties[idxs_train]
    alpha = qml.math.cho_solve(kernel_train, properties_train, l2reg=l2reg)

    test_predictions = np.dot(kernel_test, alpha)

    if not rtn_train:
        del kernel_train
        del kernel_test
        gc.collect()
        return test_predictions

    train_predictions = np.dot(kernel_train, alpha)

    del kernel_train
    del kernel_test
    gc.collect()
    return test_predictions, train_predictions


def score_rmse(kernel, properties, idxs_train, idxs_test, l2reg=1e-6):

    predictions = get_predictions(kernel, properties, idxs_train, idxs_test, l2reg=l2reg)

    # properties_train = properties[idxs_train]
    properties_test = properties[idxs_test]

    diff = predictions - properties_test
    diff = diff**2
    rmse = np.sqrt(diff.mean())

    if np.isnan(rmse):
        rmse = np.inf

    return rmse


def score_rmse_testtrain(kernel, properties, idxs_train, idxs_test, l2reg=0.0):

    predictions_test, predictions_train = get_predictions(kernel, properties, idxs_train, idxs_test, l2reg=l2reg)

    properties_train = properties[idxs_train]
    properties_test = properties[idxs_test]

    diff = predictions_test - properties_test
    diff = diff**2
    rmse_test = np.sqrt(diff.mean())

    diff = predictions_train - properties_train
    diff = diff**2
    rmse_train = np.sqrt(diff.mean())

    return rmse_test, rmse_train


def learning_curves(kernel, properties, idxs_train, idxs_test,
    score_func=score_rmse,
    training_points=None,
    check_len=True,
    **kwargs):
    """
    """

    n_items = len(idxs_train)

    if training_points is None:
        training_points = [2**x for x in range(4, 15)]

    scores = []
    for n in training_points:

        if n > n_items and check_len: break

        idxs = idxs_train[:n]

        score = score_func(kernel, properties, idxs, idxs_test, **kwargs)
        scores.append(score)

    return scores


def cross_validation(kernels, properties,
    score_func=score_rmse,
    training_points=None,
    n_splits=5):
    """

    iterate over kernels and select best index for each n_learn

    """

    fold_five = sklearn.model_selection.KFold(n_splits=n_splits, random_state=45, shuffle=True)

    n_points = len(training_points)

    kernel_scores = []
    mean_scores = []
    winner_idxs = np.zeros(n_points, dtype=int)
    winner_idxs -= 1
    winner_mean = np.zeros(n_points)
    winner_mean += np.inf

    for idxk, kernel in enumerate(kernels):

        n_items = kernel.shape[0]
        X = list(range(n_items))

        kernel_score = []

        for i, (idxs_train, idxs_test) in enumerate(fold_five.split(X)):

            # un-ordered idxs_train
            np.random.seed(45+i)
            np.random.shuffle(idxs_train)

            training_scores = learning_curves(kernel, properties, idxs_train, idxs_test,
                score_func=score_rmse,
                training_points=training_points)

            kernel_score.append(training_scores)

        kernel_score = np.array(kernel_score)
        kernel_score = kernel_score.T
        kernel_scores.append(kernel_score)

        mean_score = np.mean(kernel_score, axis=1)

        for i in range(n_points):
            this_mean = mean_score[i]
            that_mean = winner_mean[i]

            if this_mean < that_mean:
                winner_mean[i] = this_mean
                winner_idxs[i] = idxk

    scores = []
    for i, idx in enumerate(winner_idxs):
        score = kernel_scores[idx][i,:]
        scores.append(kernel_scores[idx][i])

    return winner_idxs, scores


# def cross_validation_score(kernel, properties, score_func=score_kernel,
#     n_trains=[2**x for x in range(4, 15)],
#     parameters=None,
#     **kwargs):
#     """
#
#     """
#
#     # fold-it
#     fold_five = sklearn.model_selection.KFold(n_splits=5, random_state=45, shuffle=True)
#     n_items = kernel.shape[0]
#     X = list(range(n_items))
#
#     # for hp opt
#     l2regs = [10**-x for x in range(1,10,2)] + [0.0]
#
#     reg_winners = []
#     reg_scores = []
#
#     for l2reg in l2regs:
#
#         scores = []
#
#         for idxs_train, idxs_test in fold_five.split(X):
#
#             learn_score = cross_validated_learning_curve(kernel, properties, idxs_train, idxs_test,
#                 l2reg=l2reg,
#                 n_trains=n_trains,
#                 **kwargs)
#             scores.append(learn_score)
#
#         scores = np.array(scores)
#         scores = scores.T
#
#         reg_scores.append(scores)
#
#         print(scores)
#         quit()
#
#     # for i in len():
#
#
#     quit()
#
#     return


# def cross_validated_learning_curve(kernel, properties, idxs_train, idxs_test,
#     n_trains=[2**x for x in range(4, 15)],
#     check_len=True,
#     l2reg=1e-6,
#     **kwargs):
#
#     n_items = len(idxs_train)
#
#     score_parameters = []
#
#     scores = []
#     for n in n_trains:
#
#         if n > n_items and check_len: break
#
#         idxs = idxs_train[:n]
#
#         score = score_kernel(kernel, properties, idxs, idxs_test, l2reg=l2reg, **kwargs)
#         scores.append(score)
#
#     return scores


def dump_kernel_scores(scr, names=[]):

    # Predefined reg
    l2regs = [10**-x for x in range(1, 6, 2)] + [0.0]
    n_l2regs = len(l2regs)

    # Define n_training
    # n_trains=[2**x for x in range(4, 12)]
    n_trains=[2**x for x in range(4, 17)]
    n_trains = np.array(n_trains, dtype=int)
    n_items = misc.load_txt(scr + "n_items")

    n_train_idx, = np.where(n_trains < n_items*4.0/5.0)
    n_trains = n_trains[n_train_idx]
    n_trains = list(n_trains) # + [-1]

    print("Assume total items", n_items,
            "N train", "{:5.1f}".format(np.floor(n_items*4/5)),
            "N test", "{:5.1f}".format(np.ceil(n_items*1/5)))
    print("Training:", list(n_trains))
    misc.save_npy(scr + "n_train", n_trains)

    # Load properties
    try:
        properties = misc.load_npy(scr + "properties")
    except:
        with open(scr + "properties.csv", 'r') as f:
            lines = f.readlines()
            properties = []
            for line in lines:

                values = [float(x) for x in line.split()]
                values = values[1:]
                value = np.median(values)
                properties.append(value)

            properties = np.array(properties)
            misc.save_npy(scr + "properties", properties)


    print(n_items, "==", len(properties))
    assert n_items == len(properties)

    # Load done kernel
    this_names = ["rdkitfp", "morgan"]
    for name in names:

        break

        if name not in this_names:
            continue

        print("scoring", name)

        now = time.time()

        print("load kernel", name)
        kernel = misc.load_npy(scr + "kernel." + name)

        n_len = kernel.shape[0]
        diaidx = np.diag_indices(n_len)

        def scan_kernels(debug=True):
            kernel[diaidx] += l2regs[0]
            yield kernel
            # for i in tqdm.tqdm(range(1, n_l2regs), ncols=47, ascii=True, desc=name):
            for i in range(1, n_l2regs):
                kernel[diaidx] += -l2regs[i-1] +l2regs[i]
                yield kernel

        generator = functools.partial(tqdm, scan_kernels(), ncols=75, ascii=True, desc=name+ " kernels", total=n_l2regs)

        print("scan kernels", name)
        idx_winners, scores = cross_validation(generator(), properties, training_points=n_trains)
        misc.save_npy(scr + "score."+name, scores)
        scores = np.around(np.mean(scores, axis=1), decimals=2)

        # Save parameters
        winner_parameters = {}
        for ni, index in enumerate(idx_winners):

            n = n_trains[ni]
            l2reg = l2regs[index]

            parameters = {
                "reg": l2reg,
            }

            winner_parameters[str(n)] = parameters

        nower = time.time()

        print("time: {:10.1f}s".format(nower-now))
        print(name, list(scores))

        misc.save_json(scr + "parameters."+name, winner_parameters)

        print("saved")

        kernel = None
        del kernel

    # Load multi kernels (reg search)
    this_names = ["fchl19", "fchl18"]
    for name in names:
        break
        kernels = misc.load_npy(scr + "kernels." + name)

        n_l2regs = len(l2regs)
        n_kernels = kernels.shape[0]
        n_len = kernels[0].shape[0]

        diaidx = np.diag_indices(n_len)

        def scan_kernels():
            for kernel in kernels:
                kernel[diaidx] += l2regs[0]
                yield kernel
                for i in range(1, n_l2regs):
                    kernel[diaidx] += -l2regs[i-1] +l2regs[i]
                    yield kernel

        idx_winners, scores = cross_validation(scan_kernels(), properties, training_points=n_trains)
        misc.save_npy(scr + "score."+name, scores)
        scores = np.around(np.mean(scores, axis=1), decimals=2)

        # Clean
        kernels = None
        del kernels

        # Save parameters
        winner_parameters = {}
        for ni, index in enumerate(idx_winners):

            # convert linear index to multi-dimensions
            idx_parameters = np.unravel_index([index], (n_kernels, n_l2regs))
            i, j = idx_parameters
            i = int(i[0])
            j = int(j[0])

            n = n_trains[ni]
            sigma = i
            l2reg = l2regs[j]

            parameters = {
                "sigma": sigma,
                "reg": l2reg,
            }

            winner_parameters[str(n)] = parameters

        misc.save_json(scr + "parameters."+name, winner_parameters)

        print(name, scores)


    # Load distance kernels
    models = []
    parameters = {
        "name": "rdkitfp",
        "sigma": [2**x for x in range(1, 12, 2)],
        # "sigma": [2**x for x in np.arange(20, 40, 0.5)],
        # "lambda": l2regs,
        # "lambda":  [10.0**-x for x in np.arange(1, 10, 1)]
        "lambda":  [10.0**-6],
    }
    models.append(parameters)
    parameters = {
        "name": "slatm",
        "sigma": [2**x for x in range(1, 12, 2)],
        # "sigma": [2**x for x in np.arange(20, 40, 0.5)],
        # "lambda": l2regs,
        # "lambda":  [10.0**-x for x in np.arange(1, 10, 1)]
        "lambda":  [10.0**-6],
    }
    models.append(parameters)
    parameters = {
        "name": "cm",
        "sigma": [2**x for x in range(1, 12, 2)],
        "lambda": l2regs,
    }
    models.append(parameters)
    parameters = {
        "name": "bob",
        "sigma": [2**x for x in range(1, 12, 2)],
        "lambda": l2regs,
    }
    models.append(parameters)
    parameters = {
        "name": "avgslatm",
        "sigma": [2**x for x in range(1, 20, 2)],
        "lambda": l2regs,
    }
    # models.append(parameters)

    for model in models:
        name = model["name"]

        if name not in names:
            continue

        print("scoring", name)

        parameters = model

        n_sigma = len(parameters["sigma"])
        n_lambda = len(parameters["lambda"])

        print("parameter range")
        print("sigma", min(parameters["sigma"]), max(parameters["sigma"]))

        dist = misc.load_npy(scr + "dist." + name)
        kernels = get_kernels_l2distance(dist, parameters)

        # Cross validate
        idx_winners, scores = cross_validation(kernels, properties, training_points=n_trains)

        # Save scores
        misc.save_npy(scr + "score."+name, scores)
        scores = np.around(np.mean(scores, axis=1), decimals=2)

        # Save parameters
        winner_parameters = {}
        for ni, index in enumerate(idx_winners):

            # convert linear index to multi-dimensions
            idx_parameters = np.unravel_index([index], (n_sigma, n_lambda))
            i, j = idx_parameters
            i = int(i[0])
            j = int(j[0])

            n = n_trains[ni]
            sigma = parameters["sigma"][i]
            l2reg = parameters["lambda"][j]

            this_parameters = {
                "sigma": str(sigma),
                "reg": str(l2reg),
            }

            winner_parameters[str(n)] = this_parameters


        print(name, scores)
        misc.save_json(scr + "parameters."+name, winner_parameters)



    quit()

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

            yield kernel


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


def get_fp_kernel(reps):

    # TODO Do this in parallel

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

    # TODO Properties should be read by scr!!

    # properties
    print("Saving properties")
    with open(scr + 'properties.csv', 'r') as f:
        properties = f.readlines()
        properties = [x.split()[0] for x in properties]
        properties = [float(x) for x in properties]
        properties = np.array(properties)

    print("properties", properties.shape)

    misc.save_npy(scr + "properties", properties)

    # Prepare distances
    representation_names = ["cm", "bob", "slatm"] # + ["avgslatm"]
    for name in representation_names:
        print("Distance", name)
        representations = misc.load_npy(scr + "repr." + name)
        print(representations.shape)
        dist = generate_l2_distances(representations)
        misc.save_npy(scr + "dist." + name, dist)

        dist = None
        del dist

    # Prepare fchl kernels
    if False:
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

    if False:
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

    if True:
        print("Generating fingerprint kernel")
        representations_fp = misc.load_obj(scr + "repr.fp")
        kernel = get_fp_kernel(representations_fp)
        misc.save_npy(scr + "kernel.fp", kernel)

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
    parser.add_argument('--names', nargs="+", help='', default=[])

    args = parser.parse_args()

    if args.scratch[-1] != "/":
        args.scratch += "/"

    np.random.seed(args.randomseed)

    if args.get_kernels:
        dump_distances_and_kernels(args.scratch)

    if args.get_learning_curves:
        dump_kernel_scores(args.scratch, names=args.names)

    return


if __name__ == '__main__':
    main()


