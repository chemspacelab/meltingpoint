
import cv
import fingerprints
import itertools
from chemhelp import cheminfo
import misc
import numpy as np
import scipy

from sklearn.model_selection import cross_val_score, GridSearchCV
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import MinMaxScaler



def get_best_rfr(X, y):

    # TODO
    param_grid={
        'max_depth': range(3,7),
        'n_estimators': (10, 50, 100, 1000),
    }
    gsc = GridSearchCV(
        estimator=RandomForestRegressor(),
        param_grid=param_grid,
        cv=5,
        scoring='neg_mean_squared_error', verbose=10,
        n_jobs=2)
    grid_result = gsc.fit(X, y)
    best_params = grid_result.best_params_
    print(best_params)

    quit()

    rfrargs = {
        "random_state": 0,
        "n_estimators": 1000,
        "max_depth": 5,
        "criterion": "mse",
    }

    clf = RandomForestRegressor(**rfrargs)
    clf.fit(X, y)

    return clf


def main():

    # L. Breiman, "Random Forests", Machine Learning, 45(1), 5-32, 2001.

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--scratch', action='store', help='', metavar="dir", default="_tmp_")
    parser.add_argument('--randomseed', action='store', help='random seed', metavar="int", default=1)
    parser.add_argument('-j', '--procs', action='store', help='pararallize', type=int, metavar="int", default=0)
    args = parser.parse_args()

    if args.scratch[-1] != "/":
        args.scratch += "/"

    # Not that random
    np.random.seed(args.randomseed)

    # Get properties
    properties = misc.load_npy(args.scratch + "properties")
    molobjs = cheminfo.read_sdffile(args.scratch + "structures.sdf.gz")

    X = []

    try:
        X = misc.load_npy(args.scratch + "repr.rdkitfp")
        print("loaded")
    except:
        for molobj in molobjs:
            bitmap = fingerprints.get_rdkitfp(molobj)
            X.append(bitmap)

    X = np.asarray(X)
    y = properties

    # load predefined training points
    n_train = misc.load_npy(args.scratch + "n_train")


    # CV
    idxs = np.array(list(range(len(properties))), dtype=int)
    scores = []

    for idxs_train, idxs_test in cv.cross_view(idxs):

        learning_curve = []

        for n in n_train:
            idxs = idxs_train[:n]

            clf = get_best_rfr(X[idxs], y[idxs])

            # training error
            # predictions = clf.predict(X)

            # predictions
            predictions = clf.predict(X[idxs_test])
            diff = predictions-y[idxs_test]
            diff = diff**2
            rmse_test = np.sqrt(diff.mean())
            learning_curve.append(rmse_test)
            print(n, rmse_test)

        scores.append(learning_curve)

    scores = np.array(scores)
    scores = scores.T

    mean_score = np.mean(scores, axis=1)
    print(mean_score)
    misc.save_npy(args.scratch + "score.rfr", scores)


if __name__ == '__main__':
    main()
