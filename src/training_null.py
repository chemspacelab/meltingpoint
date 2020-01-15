#!/usr/bin/env python

import sklearn.model_selection
import gzip
import os
import numpy as np

from chemhelp import cheminfo
import misc


def main():

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--scratch', action='store', help='', metavar="dir", default="_tmp_")
    parser.add_argument('--randomseed', action='store', help='random seed', metavar="int", default=1)
    parser.add_argument('-j', '--procs', action='store', help='pararallize', type=int, metavar="int", default=0)
    args = parser.parse_args()

    if args.scratch[-1] != "/":
        args.scratch += "/"

    # Not that random
    # np.random.seed(args.randomseed)

    # Get properties
    properties = misc.load_npy(args.scratch + "properties")
    # molobjs = cheminfo.read_sdffile(args.scratch + "structures.sdf.gz")

    n_items = len(properties)
    X = np.arange(n_items)

    # Train
    n_splits = 5
    n_train = misc.load_npy(args.scratch + "n_train")

    fold_five = sklearn.model_selection.KFold(
        n_splits=n_splits,
        random_state=45,
        shuffle=True)

    scores = []

    for i, (idxs_train, idxs_test) in enumerate(fold_five.split(X)):

        # un-ordered idxs_train
        np.random.seed(45+i)
        np.random.shuffle(idxs_train)

        learning_curve = []

        for n in n_train:
            idxs = idxs_train[:n]

            train = properties[idxs]
            model = train.mean()

            test = properties[idxs_test]

            # predict
            sign_diff = model - test

            # rmse
            diff = sign_diff**2
            rmse_test = np.sqrt(diff.mean())

            # save
            learning_curve.append(rmse_test)

        scores.append(learning_curve)

    scores = np.array(scores)
    scores = scores.T

    mean_score = np.mean(scores, axis=1)
    print(mean_score)
    misc.save_npy(args.scratch + "score.null", scores)

    return


if __name__ == '__main__':
    main()
