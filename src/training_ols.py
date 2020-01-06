#!/usr/bin/env python
import statsmodels.api as sm
import numpy as np
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, Descriptors3D
from rdkit.Chem import AllChem
import scipy.spatial
import gzip
import patsy
import tqdm
import rdkit.Chem.rdPartialCharges as rcr
from rdkit.Chem import PandasTools
import rdkit.Chem.Lipinski

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

def extract_features(p, sdf):
    res = []
    for row, mol in tqdm.tqdm(zip(p.iterrows(), sdf), total=len(p)):
        fs = dict()
        fs['prop'] = row[1].prop
        fs['stddev'] = row[1].stddev
        
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
            continue
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

if __name__ == '__main__':
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