"""
Caleb Bell (2016). thermo: Chemical properties component of Chemical Engineering Design Library (ChEDL)
https://github.com/CalebBell/thermo.

pip install thermo

"""


import matplotlib.pyplot as plt

import misc
import numpy as np
from chemhelp import cheminfo


import thermo

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


    # Read properties
    properties = misc.load_npy(args.scratch + "properties")
    molecules = cheminfo.read_sdffile(args.scratch + "structures.sdf.gz")
    molecules = list(molecules)


    heavy_atoms = []
    predictions = []
    errors = []

    for mol, prop in zip(molecules, properties):

        smi = cheminfo.molobj_to_smiles(mol, remove_hs=True)
        J = thermo.joback.Joback(smi)
        # J = thermo.joback.Joback('CC(=O)C')
        # J = thermo.joback.Joback('CCC(=O)OC(=O)CC')

        status = J.status

        atoms, coord = cheminfo.molobj_to_xyz(mol)
        idx = np.where(atoms != 1)
        atoms = atoms[idx]
        N = len(atoms)
        heavy_atoms.append(N)

        if "Did not match all atoms present" in status:
            errors.append(1)
            predictions.append(float("nan"))
            continue

        try:
            estimate = J.estimate()
        except TypeError:
            errors.append(1)
            predictions.append(float("nan"))
            continue


        errors.append(0)

        T_b = estimate["Tb"]
        T_m = estimate["Tm"]

        predictions.append(T_m)



    errors = np.array(errors, dtype=int)

    idx_success, = np.where(errors == 0)

    heavy_atoms = np.array(heavy_atoms)
    predictions = np.array(predictions)
    properties = np.array(properties)

    predictions = predictions[idx_success]
    properties = properties[idx_success]
    heavy_atoms = heavy_atoms[idx_success]

    print("total", errors.shape[0], "filter", idx_success.shape[0])
    print()
    print(rmse(properties, predictions))

    plt.plot(properties, properties, "-k")
    plt.scatter(properties, predictions, s=0.95, alpha=0.8, c=heavy_atoms)


    plt.xlabel("True")
    plt.ylabel("Predicted")

    plt.savefig("_fig_joback")
    plt.clf()


    return



if __name__ == '__main__':
    main()
