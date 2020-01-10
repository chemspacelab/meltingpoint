import glob

import matplotlib.pyplot as plt
import numpy as np

import misc
import views


def plot_errors(scr):

    fig, axes = plt.subplots(1, 1, figsize=(8,4))
    ax = axes

    # n_trains=[2**x for x in range(4, 4+7)]
    try:
        n_trains = misc.load_npy(scr + "n_train")
    except FileNotFoundError:
        n_trains = misc.load_txt(scr + "n_train")


    names = ["cm", "bob", "fchl18", "fchl19", "fp", "slatm"]
    names = glob.glob(scr + "score.*")

    fix_name = lambda x: x.replace(scr, "").replace(".npy", "").replace("score.", "")

    names = [fix_name(x) for x in names]

    lines = []
    last_points = []

    y_min = np.inf
    y_max = -np.inf

    for name in names:

        scores = misc.load_npy(scr + "score."+name)
        mean = scores.mean(axis=1)
        std = scores.std(axis=1)

        valid_scores, = np.where(mean < 500)
        x_mean = n_trains[valid_scores]
        mean = mean[valid_scores]
        std = std[valid_scores]

        line = ax.errorbar(x_mean, mean, std,
            fmt='-o',
            # color="k",
            capsize=3,
            lw=1,
            markersize=4, label=name.upper())

        lines.append(line)
        last_points.append(mean[-1])

        max_mean = max(mean) + max(std)
        if max_mean > y_max:
            y_max = max_mean

        min_mean = min(mean) - max(std)
        if min_mean < y_min:
            y_min = min_mean

    y_min = np.floor(y_min)
    y_min = int(np.floor(y_min / 10.0)) * 10
    y_max = int(np.ceil(y_max) / 10.0) * 10

    ykeys = []

    print("y", y_min, y_max)

    diff = y_max- y_min
    if diff < 50:
        y_min -= 40

    if y_min < 0.0:
        y_min = 50

    if y_max > 120:
        y_max = 120

    # ykeys = np.arange(y_min, y_max, 30)
    ykeys = np.geomspace(y_min, y_max, num=5)

    ykeys = [int(np.ceil(y) / 5.0) * 5 for y in ykeys]

    # ykeys = [40 +10*x for x in range(0, 12, 2)]
    xkeys = n_trains

    print("x", n_trains)

    views.learning_curve_error(ax, xkeys, ykeys,
        x_range=(10, max(n_trains)*1.3),
        y_range=(y_min*0.95, y_max*1.12))


    views.legend_colorcoded(ax, lines, names)

    # learning legends

    # idxs = np.argsort(last_points)
    # idxs = np.flip(idxs, axis=0)
    # offset = 0.06
    #
    # for n, idx in enumerate(idxs):
    #
    #     name = names[idx]
    #     point = last_points[idx]
    #     color = plt.getp(lines[idx][0], 'color')
    #
    #     ax.text(0.8, 0.46-offset*n, name.upper(),
    #         fontweight='bold',
    #         color=color,
    #         transform=ax.transAxes)
    #

    # help(ax.grid)
    # ax.grid( linestyle='-', linewidth=.5, axis="x")

    # ax.grid(True)

    ax.set_xlabel('Training set size', fontweight='medium', fontsize=11)
    ax.set_ylabel('RMSE [Kelvin]', fontweight='medium', fontsize=11)

    plt.savefig(scr + "learning_curves.png", bbox_inches="tight")
    plt.savefig(scr + "learning_curves.pdf", bbox_inches="tight")

    print(scr + "learning_curves.png")

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
