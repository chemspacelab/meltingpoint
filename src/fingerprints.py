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
import time
import itertools
import multiprocessing.managers
import os
from functools import partial
from multiprocessing import Pool, Process

import numpy as np
import rdkit
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem

import misc
from chemhelp import cheminfo

import numba
import bitmap_kernels

class MyManager(multiprocessing.managers.BaseManager):
    pass
MyManager.register('np_zeros', np.zeros, multiprocessing.managers.ArrayProxy)


def tanimoto_similarity(a, b):
    return rdkit.DataStructs.FingerprintSimilarity(a,b)


def dice_similarity(a, b):
    return rdkit.DataStructs.DiceSimilarity(a,b)


def procs_similarity_index(args,
    arraya=None,
    arrayb=None,
    similarity=rdkit.DataStructs.FingerprintSimilarity,
    kernel=None,
    dim=None,
    symmetric=True,
    **kwargs):

    i = args
    # i = args[0]
    # j = args[1]

    values = np.zeros(dim)

    for j in range(dim):
        value = similarity(arraya[i], arrayb[j])
        values[j] = value
        # value = 5
        # kernel[i,j] = value

    # if i != j and symmetric:
    #     kernel[j,i] = value

    return i, values


def procs_similarity(args, similarity=rdkit.DataStructs.FingerprintSimilarity, kernel=None, symmetric=True, **kwargs):

    idx = list(args[0])
    i, j = idx
    reps = args[1]

    value = similarity(*reps)

    # print(i,j,value)

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

    #
    print("allocate kernel")
    kernel = np.zeros((n1_items, n2_items))

    if procs == 0:

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

        # print("allocating kernel", n1_items, n2_items)
        # results = m.np_zeros((n1_items, n2_items))
        # print("done allocating")
        #
        # results[0,0] = 600

        kwargs = {}
        # kwargs["kernel"] = results
        kwargs["similarity"] = similarity
        kwargs["symmetric"] = symmetric
        kwargs["arraya"] = fps1
        kwargs["arrayb"] = fps2
        kwargs["dim"] = n2_items
        func = partial(procs_similarity_index, **kwargs)

        # lines = server(fps1, fps2)
        # lines = list(lines)

        idxa = np.arange(n1_items ,dtype=int)
        idxb = np.arange(n2_items ,dtype=int)

        # lines = itertools.product(idxa, idxb)
        # lines = list(lines)

        print("starting pool", procs)
        # results = misc.parallel(lines, func, [], {}, procs=procs, maxsize=2*procs)
        # results = list(results)
        # print("done pool")

        pool = Pool(procs)
        print("pooling kernel")
        results = pool.map(func, idxa)
        pool.close()
        pool.join()

        for i, values in results:
            kernel[i,:] = values[:]

            # TODO if symmetric

        # kernel = np.array(results)


    return kernel


def get_morgan(molobj, radius=2, bitsize=3*1024, bits=True, **kwargs):
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

    # TODO useBondTypes=False
    # TODO useFeatures=True

    if bits:
        fp = AllChem.GetMorganFingerprintAsBitVect(molobj,2,
            nBits=1024*5,
            useFeatures=True)
        fp = fp_to_bitmap(fp)

    else:
        fp = Chem.rdMolDescriptors.GetMorganFingerprint(molobj, radius)

    return fp


def get_rdkitfp(molobj, bits=True, **kwargs):

    fp1 = Chem.RDKFingerprint(molobj)

    if bits:
        fp1 = fp_to_bitmap(fp1)

    return fp1


def molobjs_to_fps(molobjs, procs=0, fingerfunc=get_rdkitfp, bits=True, **kwargs):

    kwargs["bits"] = True

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
        pool.close()
        pool.join()

    if bits:
        results = np.array(results, dtype=np.float)

    return results

def fp_to_bitmap(fp, dtype=np.float):

    bit_size = fp.GetNumBits()
    bit_str = fp.ToBitString()
    # bit_size = fp.GetLength()
    # bit_bin = fp.ToBinary()

    bitmap = np.zeros(bit_size, dtype=dtype)

    for i, char in enumerate(bit_str):
        bitmap[i] = int(char)

    return bitmap


def jaccard_index(bm1, bm2):

    dot = np.dot(bm1,bm2)
    length1 = np.sum(bm1)
    length2 = np.sum(bm2)
    length3 = np.dot(bm1, bm1)

    if dot == 0:
        return 0.0

    s = float(dot) / (length1 + length2 - dot)

    return s


def dice_coefficient(vec1, vec2):

    dot = np.dot(vec1,vec2)

    if dot == 0:
        return 0.0

    length1 = np.sum(vec1)
    length2 = np.sum(vec2)
    s = 2.0*dot / (length1+length2)

    return s


@numba.njit(parallel=True)
def bitmap_jaccard_kernel(vectors):

    lengths = np.sum(vectors, axis=1)

    shape = lengths.shape[0]
    lengths = lengths.reshape((shape, 1))

    kernel = np.dot(vectors, vectors.T)
    kernel /= (lengths + lengths.T - kernel)

    return kernel


def test_kernel():

    smiles = ['c1ccccn1']
    smiles += ['c1ccco1']
    smiles += ['Oc1ccccc1']
    smiles += ['Nc1ccccc1']
    smiles += ['CCO']
    smiles += ['CCN']
    molobjs = [cheminfo.smiles_to_molobj(x)[0] for x in smiles]

    molobjs = cheminfo.read_sdffile("_tmp_bing_bp_/structures.sdf.gz")
    molobjs = [next(molobjs) for _ in range(5000)]

    init = time.time()
    vectors = molobjs_to_fps(molobjs)

    print("init", time.time()-init)

    time_pykernel = time.time()
    kernel = bitmap_jaccard_kernel(vectors)
    print("pykernel", time.time()- time_pykernel)
    print(kernel)

    del kernel

    n_items = vectors.shape[0]
    # kernel = np.zeros((n_items, n_items))

    vectors = vectors.T
    vectors = np.array(vectors, dtype=int)

    # help(bitmap_kernels)
    time_fkernel = time.time()
    kernel = bitmap_kernels.symmetric_jaccard_kernel(n_items, vectors)
    print("fokernel", time.time()-time_fkernel)
    print(kernel)

    return


def test_dot():

    smiles = "Oc1ccccc1"
    molobj, status = cheminfo.smiles_to_molobj(smiles)
    fp1 = get_rdkitfp(molobj)
    bm = fp_to_bitmap(fp1)

    print(list(bm))

    # hello = np.array([0, 1, 0,0,0,0,0,0,0,0,0,1])
    # res = np.dot(hello, hello)

    bm = np.array(bm, dtype=int)

    s = np.sum(bm)
    other = np.dot(bm,bm)

    print(s, other)

    return

def main():

    # smiles_list = ['c1ccccn1', 'c1ccco1']*10
    # molobjs = [cheminfo.smiles_to_molobj(smiles)[0] for smiles in smiles_list]

    smiles1 = 'c1ccccn1'
    smiles2 = 'c1ccco1'
    smiles1 = 'Oc1ccccc1'
    smiles2 = 'Nc1ccccc1'
    # smiles1 = 'CCO'
    # smiles2 = 'CCN'
    molobj1, status = cheminfo.smiles_to_molobj(smiles1)
    molobj2, status = cheminfo.smiles_to_molobj(smiles2)

    fp1 = get_rdkitfp(molobj1)
    fp2 = get_rdkitfp(molobj2)
    bm1 = fp_to_bitmap(fp1)
    bm2 = fp_to_bitmap(fp2)

    print(bm1)

    print()

    sim = rdkit.DataStructs.FingerprintSimilarity(fp1,fp2)
    print(sim)

    sim = jaccard_index(bm1, bm2)
    print(sim)

    sim = dice_coefficient(bm1, bm2)
    print(sim)

    print()


    fp1 = AllChem.GetMorganFingerprintAsBitVect(molobj1,2,nBits=1024*5,useFeatures=True)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(molobj2,2,nBits=1024*5,useFeatures=True)
    bm1 = fp_to_bitmap(fp1)
    bm2 = fp_to_bitmap(fp2)

    sim = jaccard_index(bm1, bm2)
    print(sim)
    sim = rdkit.DataStructs.FingerprintSimilarity(fp1,fp2)
    print(sim)

    fp1 = get_morgan(molobj1)
    fp2 = get_morgan(molobj2)
    sim = AllChem.DataStructs.DiceSimilarity(fp1,fp2)
    print(sim)


    # molobjs = cheminfo.read_sdffile("_tmp_bing_bp_/structures.sdf.gz")
    # molobjs = [next(molobjs) for _ in range(20)]
    #
    # fingerprints = molobjs_to_fps(molobjs, procs=2)
    # kernel = fingerprints_to_kernel(fingerprints, fingerprints, procs=2, similarity=dice_similarity)
    #
    # print(kernel)

    return

if __name__ == '__main__':
    test_kernel()
