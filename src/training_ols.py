#!/usr/bin/env python

import sklearn.model_selection
from functools import partial
from multiprocessing import Process, Pool
import gzip
import os
import numpy as np
import pandas as pd
import patsy
import rdkit
import rdkit.Chem.Lipinski
import rdkit.Chem.rdPartialCharges as rcr
import scipy.spatial
import statsmodels.api as sm
import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Descriptors3D#, PandasTools

from chemhelp import cheminfo
import misc

def read_dataset(basepath):
    properties = pd.read_csv(basepath + "properties.csv", sep=" ", names="prop stddev".split())
    fh = gzip.open(basepath + 'structures.sdf.gz')
    sdf = Chem.ForwardSDMolSupplier(fh, removeHs=False, sanitize=True)
    mols = [_ for _ in sdf]
    return properties, mols

def get_dipole(mol):
    def nuclei_nuclei(coordinates, charges):
        shift = coordinates - coordinates.mean(axis=0)
        return np.sum(shift.T * charges, axis=1)
    coords = mol.GetConformer(0).GetPositions()
    charges = np.zeros(mol.GetNumAtoms())
    charges = [float(_.GetAtomicNum()) for _ in mol.GetAtoms()]
    return np.linalg.norm(nuclei_nuclei(coords, charges))


def extract_feature(*args):

    if isinstance(args[0], tuple):
        row = args[0][0]
        mol = args[0][1]
    else:
        row = args[0]
        mol = args[1]

    fs = dict()
    fs['prop'] = row
    # fs['prop'] = row[1].prop
    # fs['stddev'] = row[1].stddev

    # counts
    fs['natoms'] = mol.GetNumAtoms()
    fs['naromatic'] = len(mol.GetAromaticAtoms())
    fs['nheavy'] = mol.GetNumHeavyAtoms()
    fs['nbonds'] = mol.GetNumBonds()
    fs['molwt'] = Descriptors.MolWt(mol)

    # 3d properties
    coords = mol.GetConformer(0).GetPositions()
    hull = scipy.spatial.ConvexHull(coords, qhull_options='QJ')
    fs['convexarea'] = hull.area
    fs['convexvolume'] = hull.volume

    try:
        fs['dipole'] = get_dipole(mol)
    except:
        return False

    for kind in ['RadiusOfGyration',]:
        fs[kind.lower()] = getattr(Descriptors3D, kind)(mol)
    fs['volume'] = AllChem.ComputeMolVolume(mol)    

    # environments
    for atom in mol.GetAtoms():
        this = atom.GetSymbol()
        neighbors = sorted([_.GetSymbol() for _ in atom.GetNeighbors()])
        countkey(fs, 'atom_%s' % this)
        countkey(fs, 'atomenv_%s%s' % (this, ''.join(neighbors)))

    for bond in mol.GetBonds():
        symbols = bond.GetBeginAtom().GetSymbol(), bond.GetEndAtom().GetSymbol()
        border = int(bond.GetBondType())
        bondkey = 'bond_%s%d' % (''.join(sorted(symbols)), border)
        countkey(fs, bondkey)

    # Lipinski
    fs['HBaccept'] = rdkit.Chem.Lipinski.NumHAcceptors(mol)
    fs['HBdonate'] = rdkit.Chem.Lipinski.NumHDonors(mol)
    fs['HBnhoh'] = rdkit.Chem.Lipinski.NHOHCount(mol)
    fs['HBno'] = rdkit.Chem.Lipinski.NOCount(mol)
    fs['rotatablebonds'] = rdkit.Chem.Lipinski.NumRotatableBonds(mol) 

    return fs

def extract_features(p, sdf, procs=0):
    res = []

    if procs > 0:

        def generate_input():
            for row, mol in zip(p, sdf):
                yield (row, mol)

        workargs = generate_input()
        funcname = extract_feature
        pool = Pool(processes=procs)
        # results = list(tqdm.tqdm(pool.imap(funcname, workargs), total=len(p)))
        res = pool.map(funcname, workargs)

    else:
        for row, mol in tqdm.tqdm(zip(p, sdf), total=len(p)):

            fs = extract_feature(row, mol)
            print(mol)

            # Finalise
            res.append(fs)


    return res

def countkey(d, key):
    if key not in d:
        d[key] = 0
    d[key] += 1
    return d

def regression(df, columns, powers, train_idx, test_idx):
    df = df.copy()
    for column in columns:
        try:
            df[column] = df[column] ** powers[column.split('_')[0]]
        except:
            continue
    
    y, X = patsy.dmatrices('prop ~ %s' % (' + '.join(columns)), data=df, return_type="matrix")
    mod = sm.OLS(y[train_idx], X[train_idx])
    result = mod.fit()
    
    residuals = result.predict(X[test_idx]) - y[test_idx, 0]
    return residuals
    
def fit_model(table, train_idx, test_idx):
    # select columns to include
    columns = []
    for colname in table.columns:
        if colname == 'prop':
            continue
        nonzero = np.sum(table[colname].values != 0)
        if nonzero < len(table)/100:
            continue
        columns.append(colname)

    # powers (small, but noticeable improvement)
    powers = {'volume': 0.25, 'radiusofgyration': -0.5, 'nheavy': 0.8, 'nbonds': 1., 'natoms': 1.25, 'naromatic': 1., 'molwt': -0.5, 'dipole': 0.2, 'atom': 0.85, 'atomenv': 0.85, 'bond': 1., 'convexarea': 0.85, 'convexvolume': 0.85,}
    
    # regression
    residuals = regression(table, columns, powers, train_idx, test_idx)
    return residuals


def example():
    basepath = "/mnt/c/Users/guido/workcopies/avl-notebooks/guido/melting/data/phase_transistions/bing_mp/"

    table = pd.DataFrame(extract_features(*read_dataset(basepath)))

    # important: remove NaN values
    table = table.fillna(0)

    # optional: selection 
    table = table.query("prop < 700 & natoms <40 & stddev < 1")

    # poor man's cross validation
    for i in range(5):
        table = table.sample(frac=1).reset_index(drop=True)
        print (np.abs(fit_model(table, list(range(1000)), list(range(1000, len(table))))).mean())


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
    np.random.seed(args.randomseed)

    # Get properties
    properties = misc.load_npy(args.scratch + "properties")
    molobjs = cheminfo.read_sdffile(args.scratch + "structures.sdf.gz")

    # Get features
    filename = "repr.ols"
    if os.path.exists(args.scratch + filename + ".pkl"):
        features = misc.load_obj(args.scratch + filename)

    else:
        features = extract_features(properties, molobjs, procs=args.procs)
        features = pd.DataFrame(features)
        features = features.fillna(0)
        misc.save_obj(args.scratch + filename, features)


    n_items = len(features)
    X = np.arange(n_items)

    # Train
    n_splits = 5
    n_train = misc.load_npy(args.scratch + "n_train")

    fold_five = sklearn.model_selection.KFold(
        n_splits=n_splits,
        random_state=45,
        shuffle=True)

    scores = []

    for idxs_train, idxs_test in fold_five.split(X):

        learning_curve = []

        for n in n_train:
            idxs = idxs_train[:n]

            # signed difference
            sign_diff = fit_model(features, idxs, idxs_test)

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
    misc.save_npy(args.scratch + "score.ols", scores)

    return


if __name__ == '__main__':
    main()
