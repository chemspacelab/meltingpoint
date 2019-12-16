"""

List of Available Fingerprints in rdkit:
Fingerprint Type    Notes   Language
RDKit   a Daylight-like fingerprint based on hashing molecular subgraphs    C++
Atom Pairs  JCICS 25:64-73 (1985)   C++
Topological Torsions    JCICS 27:82-5 (1987)    C++
MACCS keys  Using the 166 public keys implemented as SMARTS C++
Morgan/Circular Fingerprints based on the Morgan algorithm, similar to the ECFP/FCFP fingerprints JCIM 50:742-54 (2010).    C++
2D Pharmacophore    Uses topological distances between pharmacophoric points.   C++
Pattern a topological fingerprint optimized for substructure screening  C++
Extended Reduced Graphs Derived from the ErG fingerprint published by Stiefl et al. in JCIM 46:208–20 (2006).
NOTE: these functions return an array of floats, not the usual fingerprint types


The Tanimoto coefficent is determined by looking at the number of chemical
features that are common to both molecules (the intersection of the data
strings) compared to the number of chemical features that are in either (the
union of the data strings). The Dice coefficient also compares these values but
using a slightly different weighting. 

The Tanimoto coefficient is the ratio of the number of features common to both
molecules to the total number of features, i.e.

( A intersect B ) / ( A + B - ( A intersect B ) ) 

The range is 0 to 1 inclusive. 

The Dice coefficient is the number of features in common to both molecules
relative to the average size of the total number of features present, i.e.

( A intersect B ) / 0.5 ( A + B ) 

The weighting factor comes from the 0.5 in the denominator. The range is 0 to 1.


"""

import multiprocessing.managers
import os
from functools import partial
from multiprocessing import Pool, Process

import numpy as np
import rdkit
import rdkit.Chem as Chem

from chemhelp import cheminfo

class MyManager(multiprocessing.managers.BaseManager):
    pass
MyManager.register('np_zeros', np.zeros, multiprocessing.managers.ArrayProxy)

def tanimoto_similarity(a, b):
    return rdkit.DataStructs.FingerprintSimilarity(a,b)

def dice_similarity(a, b):
    return rdkit.DataStructs.DiceSimilarity(a,b)

def procs_similarity(args, similarity=rdkit.DataStructs.FingerprintSimilarity, kernel=None, symmetric=True):

    idx = list(args[0])
    i, j = idx
    reps = args[1]

    value = similarity(*reps)

    if kernel is not None:

        kernel[i,j] = value

        if i != j and symmetric:
            kernel[j,i] = value

        return

    else:

        return idx, value


def fingerprints_to_kernel(fps1, fps2, procs=0, similarity=rdkit.DataStructs.FingerprintSimilarity):

    if id(fps1) == id(fps2):
        symmetric = True
    else:
        symmetric = False

    n1_items = len(fps1)
    n2_items = len(fps2)

    def server(a, b):
        for i, aval in enumerate(a):
            for j, bval in enumerate(b):
                if i > j and symmetric: continue
                workpackage = (i,j), (aval,bval)
                yield workpackage

    if procs == 0:

        kernel = np.zeros((n1_items, n2_items))
        for idx, reps in server(fps1, fps2):
            value = similarity(*reps)
            kernel[idx] = value
            if symmetric:
                ridx = list(reversed(idx))
                kernel[ridx] = value

    else:

        m = MyManager()
        m.start()
        # NOTE call this with [i,j] and NOT [i][j]
        results = m.np_zeros((n1_items, n2_items))

        kwargs = {}
        kwargs["kernel"] = results
        kwargs["similarity"] = similarity
        kwargs["symmetric"] = symmetric
        func = partial(procs_similarity, **kwargs)
        pool = Pool(procs)
        pool.map(func, server(fps1, fps2))

        kernel = np.array(results)

    return kernel


def get_morgan(molobj, radius=2, **kwargs):
    """

    This family of fingerprints, better known as circular fingerprints [5], is
    built by applying the Morgan algorithm to a set of user-supplied atom
    invariants. When generating Morgan fingerprints, the radius of the
    fingerprint must also be provided :

    The default atom invariants use connectivity information similar to those
    used for the well known ECFP family of fingerprints. Feature-based
    invariants, similar to those used for the FCFP fingerprints, can also be
    used. The feature definitions used are defined in the section Feature
    Definitions Used in the Morgan Fingerprints. At times this can lead to
    quite different similarity scores:



    Feature Definitions Used in the Morgan Fingerprints
    These are adapted from the definitions in Gobbi, A. & Poppinger, D.
    “Genetic optimization of combinatorial libraries.” Biotechnology and
    Bioengineering 61, 47-54 (1998).

    Feature SMARTS
    Donor   [$([N;!H0;v3,v4&+1]),$([O,S;H1;+0]),n&H1&+0]
    Acceptor    [$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=[O,N,P,S])]),n&H0&+0,$([o,s;+0;!$([o,s]:n);!$([o,s]:c:n)])]
    Aromatic    [a]
    Halogen [F,Cl,Br,I]
    Basic   [#7;+,$([N;H2&+0][$([C,a]);!$([C,a](=O))]),$([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))]),$([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))])]
    Acidic  [$([C,S](=[O,S,P])-[O;H1,-1])]


    >>> m1 = Chem.MolFromSmiles('c1ccccn1')
    >>> m2 = Chem.MolFromSmiles('c1ccco1')
    >>> fp1 = AllChem.GetMorganFingerprint(m1, 2)
    >>> fp2 = AllChem.GetMorganFingerprint(m2, 2)
    >>> ffp1 = AllChem.GetMorganFingerprint(m1, 2, useFeatures=True)
    >>> ffp2 = AllChem.GetMorganFingerprint(m2, 2, useFeatures=True)
    >>> DataStructs.DiceSimilarity(fp1,fp2)
    0.36...
    >>> DataStructs.DiceSimilarity(ffp1,ffp2)
    0.90...
    """

    fp = Chem.rdMolDescriptors.GetMorganFingerprint(molobj, radius)

    return fp


def get_rdkit(molobj, **kwargs):

    fp1 = Chem.RDKFingerprint(molobj)

    return fp1


def molobjs_to_fps(molobjs, procs=0, fingerfunc=get_rdkit, **kwargs):

    results = []

    if procs == 0:
        for molobj in molobjs:
            results.append(fingerfunc(molobj, **kwargs))

    else:

        os.environ["OMP_NUM_THREADS"] = "1"
        workargs = list(molobjs)
        pool = Pool(processes=procs)
        funcname = partial(fingerfunc, **kwargs)
        results = pool.map(funcname, workargs)

    return results


def main():

    # smiles_list = ['c1ccccn1', 'c1ccco1']*10
    # molobjs = [cheminfo.smiles_to_molobj(smiles)[0] for smiles in smiles_list]

    molobjs = cheminfo.read_sdffile("_tmp_bing_bp_/structures.sdf.gz")
    molobjs = list(molobjs)
    molobjs = molobjs[:500]

    fingerprints = molobjs_to_fps(molobjs, procs=2)
    kernel = fingerprints_to_kernel(fingerprints, fingerprints, procs=2, similarity=dice_similarity)

    print(kernel)

    return

if __name__ == '__main__':
    main()
