
import statsmodels.api as sm

from scipy import stats
from scipy.stats import halfnorm
from scipy.stats import gaussian_kde
from scipy.stats import norm

import rdkit
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import rdkit.Chem.ChemicalForceFields as ChemicalForceFields
import rdkit.Chem.rdMolDescriptors as rdMolDescriptors

from chemhelp import cheminfo
import glob

import matplotlib.pyplot as plt
import numpy as np

import misc
import views

import matplotlib as mpl

plt.rc('font', size=14)

fontname = "Fira Sans"
fontweight = "bold"
# color_std = "#d81b6a"
# color_hl = "#800031"

plt.rc('legend', fontsize=15)
mpl.rcParams['font.sans-serif'] = fontname
mpl.rcParams['font.family'] = "sans-serif"
mpl.rcParams['font.weight'] = fontweight

fontargs = {
    "fontweight": fontweight,
    "fontname": fontname
}


def canonical(smiles):
    """
    Translate smiles into a canonical form
    """

    molobj = Chem.MolFromSmiles(smiles)
    smiles = Chem.MolToSmiles(molobj, canonical=True)

    return smiles


def stoichiometry(smiles, typ="str", include_hydrogen=True):
    """


    """

    molobj = Chem.MolFromSmiles(smiles)

    if include_hydrogen:
        molobj = Chem.AddHs(molobj)

    atoms = [atom for atom in molobj.GetAtoms()]

    if typ == "str":
        atoms = [atom.GetSymbol() for atom in atoms]
        atoms = np.array(atoms)
        atomkey, atomcount = np.unique(atoms, return_counts=True)
        idxs = np.argsort(atomkey)
        rtnstoi = ""
        for idx in idxs:
            rtnstoi += atomkey[idx]
            rtnstoi += str(atomcount[idx])

    else:
        atoms = [atom.GetAtomicNum() for atom in atoms]
        # not needed
        rtnstoi = atoms

    return rtnstoi


@misc.memory.cache
def split_dict(filename, load_func):

    data = load_func(filename)

    xvalues = []
    yvalues = []
    yvalues_n = []
    yvalues_std = []
    yvalues_diff = []

    yvalues_distances = []

    for key in data.keys():

        value = data[key]

        if isinstance(value, list):

            uvalues = np.round(value, decimals=0)
            uvalues = np.unique(uvalues)

            n_values = len(uvalues)
            if n_values < 2:
                std_values = 400
            else:
                # std_values = np.std(uvalues, ddof=1.5) / 0.8
                std_values = np.std(uvalues)

            for i in range(n_values):
                for j in range(i+1, n_values):
                    diff = abs(uvalues[i]-uvalues[j])
                    yvalues_distances.append(diff)

            yvalues_n.append(n_values)
            yvalues_std.append(std_values)

            value = np.mean(value)

        if isinstance(value, dict):
            value = value["K"]
            value = np.mean(value)

        yvalues.append(value)

        smiles = key
        stoi = stoichiometry(smiles, include_hydrogen=False, typ="int")
        N = len(stoi)
        xvalues.append(N)

    print(xvalues)

    xvalues = np.array(xvalues, dtype=int)
    yvalues = np.array(yvalues)
    yvalues_std = np.array(yvalues_std)
    yvalues_n = np.array(yvalues_n)
    yvalues_distances = np.array(yvalues_distances)

    return xvalues, yvalues, yvalues_n, yvalues_std, yvalues_distances


def view_std_values(yvalstd, filename):

    # std of zero is not interesting
    # idx, = np.where(yvalstd > 1.0)
    # yvalstd = yvalstd[idx]

    # HACKED
    n_points = yvalstd.shape[0]
    yvalstd_copy = np.zeros(2*n_points)
    yvalstd_copy[:n_points] = yvalstd
    yvalstd_copy[n_points:] = -yvalstd
    # mu, std = norm.fit(yvalstd_copy)


    # Half-norm fit
    # mu, std = halfnorm.fit(yvalstd, floc=0, scale=1.0)
    # print("gaussian on data", what)

    # Trunc-norm fit
    mu = 0.0
    std = stats.truncnorm.fit(yvalstd, 1.0, 100.0, floc=0, scale=1.0)
    std = std[-1]
    std = np.sqrt(std)/np.sqrt(2)

    print("std:", std)

    fig, axs = plt.subplots(1, 1, figsize=(6,6))

    ax = axs

    min_val = np.min(yvalstd)
    max_val = np.max(yvalstd)

    if False:
        bins = np.linspace(min_val, max_val, 300)
        gaussian_kernel = gaussian_kde(yvalstd)
        values = gaussian_kernel(bins)
        ax.plot(bins, values, "k", linewidth=1.0)

    else:

        n, bins, patches = ax.hist(yvalstd, bins=30, histtype='stepfilled', color="k", density=False)

        hist, bins = np.histogram(yvalstd, density=True, bins=30)
        # figqq = sm.qqplot(yvalstd, norm, fit=True, line='45')
        # figqq.savefig(filename + "_qq")
        # bincentres = [(bins[i]+bins[i+1])/2. for i in range(len(bins)-1)]
        # barplot = ax.plot(bincentres, hist, 'b.')
        # barplt = ax.bar(bin_edges, hist, align="edge")


    # x = np.linspace(min_val, max_val, 1000)
    # p = norm.pdf(x, mu, std)
    # # p *= 4000
    # ax.plot(x, p, 'r-')

    xticks = ax.get_xticks()[1:]
    yticks = ax.get_yticks()

    tick_width = xticks[1] - xticks[0]
    idx_max, = np.where(xticks < max_val+tick_width)
    xticks = xticks[idx_max]

    views.set_border(ax, xticks, yticks)

    # inset ax
    right_inset_ax = fig.add_axes([.4, .4, .4, .4])

    max_val = bins[0] + bins[1]-bins[0]
    idx, = np.where(yvalstd > max_val)
    zoom_yval = yvalstd[idx]

    idx, = np.where(zoom_yval < 250)
    zoom_yval = zoom_yval[idx]

    n, bins, patches = right_inset_ax.hist(zoom_yval, bins=30, histtype='stepfilled', color="k")
    xticks = right_inset_ax.get_xticks()
    idx, = np.where(xticks > 0)
    xticks = [0] + list(xticks[idx])
    yticks = right_inset_ax.get_yticks()
    views.set_border(right_inset_ax, xticks, yticks)
    right_inset_ax.set_xlim( xticks[0]-0.5*(xticks[1]-xticks[0]), xticks[-1] )
    right_inset_ax.set_xlabel(f"First bin removed")

    # Save
    views.save(filename + "_std", fig=fig)

    return


def view_values_molecules(xvalues, yvalues, filename):
    """
    xvalues - no of atoms
    yvalues - phase transistion
    """

    print("Max atoms:", xvalues.max())
    print("Min atoms:", xvalues.min())
    print("Mean and std atoms:", xvalues.mean(), xvalues.std())

    y_mean = yvalues.mean()
    y_std = yvalues.std()

    print("Max value", yvalues.max())
    print("Min value", yvalues.min())
    print("Mean value", yvalues.mean(), yvalues.std())

    # Filter
    max_atoms = 90
    idxs, = np.where(xvalues < max_atoms)
    xvalues = xvalues[idxs]
    yvalues = yvalues[idxs]

    # Fillter outliers
    idxs, = np.where(yvalues < y_mean+y_std*4)
    xvalues = xvalues[idxs]
    yvalues = yvalues[idxs]

    n_items = xvalues.shape[0]

    print(f"Total {n_items} items for {filename}")

    views.histogram_2d_with_kde(xvalues, yvalues, filename=filename+"_overview")

    return


def main():

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--dict', action='store', help='', metavar="FILE")
    parser.add_argument('--json', action='store', help='', metavar="FILE")

    args = parser.parse_args()

    if args.dict:
        xvalues, yvalues, yvaln, yvalstd, yvaldist = split_dict(args.dict, misc.load_obj)
        filename = args.dict

    if args.json:
        xvalues, yvalues = split_dict(args.json, misc.load_obj)
        filename = args.json

    # 2d histogram
    view_values_molecules(xvalues, yvalues, filename)

    if yvaln.shape[0] == 0:
        return

    print("no. unique distances", len(yvaldist))

    idx, = np.where(yvaln >= 2)
    print("no. of values>2", len(idx))
    # yvalstd = yvalstd[idx]

    view_std_values(yvaldist, filename)

    return

if __name__ == '__main__':
    main()
