
import matplotlib.pyplot as plt

import misc
import numpy as np
from chemhelp import cheminfo

from scipy.spatial import ConvexHull, distance

def rmse(X, Y):
    """
    Root-Mean-Square Error

    Lower Error = RMSE \left( 1- \sqrt{ 1- \frac{1.96\sqrt{2}}{\sqrt{N-1}} }  \right )
    Upper Error = RMSE \left(    \sqrt{ 1+ \frac{1.96\sqrt{2}}{\sqrt{N-1}} } - 1 \right )

    This only works for N >= 8.6832, otherwise the lower error will be
    imaginary.

    Parameters:
    X -- One dimensional Numpy array of floats
    Y -- One dimensional Numpy array of floats

    Returns:
    rmse -- Root-mean-square error between X and Y
    le -- Lower error on the RMSE value
    ue -- Upper error on the RMSE value
    """

    N, = X.shape

    if N < 9:
        print("Not enough points. {} datapoints given. At least 9 is required".format(N))
        return

    diff = X - Y
    diff = diff**2
    rmse = np.sqrt(diff.mean())

    le = rmse * (1.0 - np.sqrt(1-1.96*np.sqrt(2.0)/np.sqrt(N-1)))
    ue = rmse * (np.sqrt(1 + 1.96*np.sqrt(2.0)/np.sqrt(N-1))-1)

    return rmse, le, ue


def mae(X, Y):
    """
    Mean Absolute Error (MAE)

    Lower Error =  MAE_X \left( 1- \sqrt{ 1- \frac{1.96\sqrt{2}}{\sqrt{N-1}} }  \right )
    Upper Error =  MAE_X \left(  \sqrt{ 1+ \frac{1.96\sqrt{2}}{\sqrt{N-1}} }-1  \right )

    Parameters:
    X -- One dimensional Numpy array of floats
    Y -- One dimensional Numpy array of floats

    Returns:
    mae -- Mean-absolute error between X and Y
    le -- Lower error on the MAE value
    ue -- Upper error on the MAE value
    """

    N, = X.shape

    mae = np.abs(X - Y)
    mae = mae.mean()

    le =  mae * (1 - np.sqrt(1 - 1.96*np.sqrt(2)/np.sqrt(N-1) ) )
    ue =  mae * (    np.sqrt(1 + 1.96*np.sqrt(2)/np.sqrt(N-1) ) -1 )

    return mae, le, ue



def main():

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--scratch', action='store', help='', metavar="dir", default="_tmp_")
    parser.add_argument('-j', '--procs', action='store', help='pararallize', metavar="int", default=0, type=int)

    args = parser.parse_args()

    if args.scratch[-1] != "/":
        args.scratch += "/"


    properties = misc.load_npy(args.scratch + "properties")
    molecules = cheminfo.read_sdffile(args.scratch + "structures.sdf.gz")


    heavy_atoms = []
    distances = []
    volumes = []

    for mol in molecules:

        # atoms = cheminfo.molobj_to_atoms(mol)
        atoms, coord = cheminfo.molobj_to_xyz(mol)

        idx = np.where(atoms != 1)
        atoms = atoms[idx]
        N = len(atoms)
        heavy_atoms.append(N)

        hull = ConvexHull(coord, qhull_options="QJ")

        vol = hull.volume
        volumes.append(vol)

        avgdist = distance.pdist(coord)
        avgdist = np.mean(avgdist)

        distances.append(avgdist)


    heavy_atoms = np.array(heavy_atoms)
    volumes = np.array(volumes)
    distances = np.array(distances)

    #
    #
    #

    representation = distances

    # linear fit
    p = np.polyfit(representation, properties, 3)
    p = np.poly1d(p)

    results = p(representation)
    rmse_error = rmse(results, properties)

    print(rmse_error)

    plt.scatter(representation, properties, c=heavy_atoms, s=0.8)
    x_prop = np.linspace(min(representation), max(representation), 80)
    plt.plot(x_prop, p(x_prop), "k-")

    plt.savefig("i_can_member_it")
    plt.clf()

    return



if __name__ == '__main__':
    main()
