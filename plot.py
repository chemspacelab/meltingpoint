import glob

import matplotlib.pyplot as plt
import numpy as np

import misc
import views


def plot_errors(scr):

    fig, axes = plt.subplots(1, 1, figsize=(8,4))
    ax = axes

    # n_trains=[2**x for x in range(4, 4+7)]
    n_trains = misc.load_npy(scr + "n_train")

    names = ["cm", "bob", "fchl18", "fchl19", "fp", "slatm"]
    names = glob.glob(scr + "score.*")

    fix_name = lambda x: x.replace(scr, "").replace(".npy", "").replace("score.", "")

    names = [fix_name(x) for x in names]

    lines = []
    last_points = []

    for name in names:

        scores = misc.load_npy(scr + "score."+name)
        mean = scores.mean(axis=1)
        std = scores.std(axis=1)

        line = ax.errorbar(n_trains, mean, std,
            fmt='-o',
            # color="k",
            capsize=3,
            lw=1,
            markersize=4, label=name.upper())

        lines.append(line)
        last_points.append(mean[-1])

    ykeys = [300, 150, 75, 40]
    xkeys = n_trains

    views.learning_curve_error(ax, xkeys, ykeys,
        x_range=(10, max(n_trains)*5),
        y_range=(35, 350))


    # views.legend_colorcoded(ax, lines, names)

    # learning legends

    idxs = np.argsort(last_points)
    idxs = np.flip(idxs, axis=0)
    offset = 0.06

    for n, idx in enumerate(idxs):

        name = names[idx]
        point = last_points[idx]
        color = plt.getp(lines[idx][0], 'color')

        ax.text(0.8, 0.46-offset*n, name.upper(),
            fontweight='bold',
            color=color,
            transform=ax.transAxes)
    #


    ax.set_xlabel('Training set size', fontweight='medium', fontsize=11)
    ax.set_ylabel('RMSE [K]', fontweight='medium', fontsize=11)

    plt.savefig(scr + "learning_curves.png", bbox_inches="tight")
    plt.savefig(scr + "learning_curves.pdf", bbox_inches="tight")



    return


def main():

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--scratch', action='store', help='', metavar="dir", default="_tmp_")

    args = parser.parse_args()

    if args.scratch[-1] != "/":
        args.scratch += "/"

    plot_errors(args.scratch)

    return


if __name__ == "__main__":
    main()
