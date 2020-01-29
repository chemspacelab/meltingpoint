
import numpy as np
import sklearn.model_selection

def cross_view(items, n_splits=5, randomseed=45):

    fold_obj = sklearn.model_selection.KFold(
        n_splits=n_splits,
        random_state=randomseed,
        shuffle=True)

    for i, (idxs_train, idxs_test) in enumerate(fold_obj.split(items)):

        # un-ordered idxs_train
        np.random.seed(randomseed+i)
        np.random.shuffle(idxs_train)

        yield idxs_train, idxs_test


