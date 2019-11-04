
def plot_errors(scr):

    fig, axes = plt.subplots(1, 1, figsize=(4,4))
    ax = axes

    n_trains=[2**x for x in range(4, 4+7)]
    names = ["cm", "bob", "fchl18", "fchl19", "fp", "slatm"]

    for name in names:

        scores = misc.load_npy(scr + "score."+name)
        mean = scores.mean(axis=1)
        std = scores.std(axis=1)

        ax.errorbar(n_trains, mean, std,
            fmt='-o',
            capsize=3,
            markersize=4, label=name.upper())

    ykeys = [300, 150, 75, 40]
    xkeys = n_trains

    views.learning_curve_error(ax, xkeys, ykeys,
        x_range=(10, 1100),
        y_range=(35, 350))

    leg = ax.legend(ncol=2, frameon=False)

    ax.set_xlabel('Training set size', fontweight='medium', fontsize=13)
    ax.set_ylabel('RMSE [K]', fontweight='medium', fontsize=13)

    plt.savefig("learning_melt.png", bbox_inches="tight")
    plt.savefig("learning_melt.pdf", bbox_inches="tight")

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

