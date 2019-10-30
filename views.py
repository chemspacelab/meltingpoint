
from qml.kernels import kpca
import numpy as np
import matplotlib.pyplot as plt


def kde_2d():

    np.random.seed(20)

    data = np.random.multivariate_normal((0, 0), [[0.8, 0.05], [0.05, 0.7]], 100)
    x = data[:, 0]
    y = data[:, 1]
    xmin, xmax = -3, 3
    ymin, ymax = -3, 3

    # Peform the kernel density estimate
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([x, y])
    kernel = st.gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)

    fig = plt.figure()
    ax = fig.gca()
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    # Contourf plot
    cfset = ax.contourf(xx, yy, f, cmap='Blues')
    ## Or kernel density estimate plot instead of the contourf plot
    #ax.imshow(np.rot90(f), cmap='Blues', extent=[xmin, xmax, ymin, ymax])
    # Contour plot
    cset = ax.contour(xx, yy, f, colors='k', linewidths=0.7)
    # Label plot
    ax.clabel(cset, inline=1, fontsize=8)
    # ax.set_xlabel('Y1')
    # ax.set_ylabel('Y0')

    plt.savefig("fig_kde")


    return


def pca_with_properties(kernel, properties, filename):

    pca = kpca(kernel, n=2)

    fig, axs = plt.subplots(2, 1, figsize=(5,10))
    sc = axs[0].scatter(*pca, c=properties)
    fig.colorbar(sc, ax=axs[0])
    im = axs[1].imshow(kernel)
    fig.colorbar(im, ax=axs[1])
    fig.savefig(filename)

    return


def histogram_2d_with_kde(xvalues, yvalues,
        xlabel="# Heavy atoms",
        ylabel="Melting point [Kelvin]",
        filename="overview_scathis",
        debug=False):


    from matplotlib import rc

    # plt.rc('text', usetex=True)
    # plt.rc('font', family='serif')
    # plt.rc('font', size=18)
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


    nullfmt = NullFormatter()

    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.1]
    rect_histy = [left_h, bottom, 0.1, height]

    plt.figure(1, figsize=(8, 8))

    ax_scatter = plt.axes(rect_scatter)
    ax_histx = plt.axes(rect_histx)
    ax_histy = plt.axes(rect_histy)


    # scatter plot
    # ax_scatter.scatter(xvalues, yvalues, color="k", alpha=0.4)

    lineswidth=0.0 # white lines
    lineswidth=0.2 # perfect fit
    colormap = 'Greys'
    colormap = 'PuRd'

    hb = ax_scatter.hexbin(xvalues, yvalues, gridsize=60, bins='log', cmap='PuRd', linewidths=0.3, mincnt=1)

    # define binwidth
    x_max = np.max(xvalues)
    x_min = np.min(xvalues)
    x_binwidth = (abs(x_min) + x_max) / 30.0
    x_binwidth = int(x_binwidth)
    x_binwidth = 1
    x_bins = np.arange(x_min, x_max+x_binwidth, x_binwidth)

    y_max = np.max(yvalues)
    y_min = np.min(yvalues)
    y_binwidth = (abs(y_min) + y_max) / 50.0
    y_binwidth = int(y_binwidth)
    y_bins = np.arange(y_min, y_max+y_binwidth, y_binwidth)

    # xlim = (x_min-x_binwidth, x_max+x_binwidth)
    # ylim = (y_min-y_binwidth, y_max+y_binwidth)

    # Set limits and ticks of scatter
    xlim = (0, 70)
    ax_scatter.set_xlim(xlim)
    # ax_scatter.set_ylim(ylim)

    xkeys = np.arange(10, x_max+x_binwidth*2, 10)
    xkeys = [1] + list(xkeys)
    ykeys = np.arange(0, y_max+y_binwidth, 100)

    # Histogram

    bins = np.linspace(x_min, x_max, 200)
    gaussian_kernel = gaussian_kde(xvalues)
    values = gaussian_kernel(bins)
    ax_histx.plot(bins, values, "k", linewidth=1.0)

    bins = np.linspace(y_min, y_max, 200)
    gaussian_kernel = gaussian_kde(yvalues)
    values = gaussian_kernel(bins)
    ax_histy.plot(values, bins, "k", linewidth=1.0)

    # ax_histx.hist(xvalues, bins=x_bins, histtype='step', color="k")
    # ax_histy.hist(yvalues, bins=y_bins, orientation='horizontal', histtype='step', color="k")

    ax_histx.set_xlim(ax_scatter.get_xlim())
    ax_histy.set_ylim(ax_scatter.get_ylim())

    # pretty
    if not debug:
        ax_histx.xaxis.set_major_formatter(nullfmt)
        ax_histy.yaxis.set_major_formatter(nullfmt)

        set_border(ax_scatter, xkeys, ykeys)
        set_border(ax_histx, [], [], border=[False, False, False, False])
        set_border(ax_histy, [], [], border=[False, False, False, False])

        ax_histx.set_xticks([], [])
        ax_histy.set_yticks([], [])


    ax_scatter.set_xlabel(xlabel, **fontargs)
    ax_scatter.set_ylabel(ylabel, **fontargs)

    plt.savefig(filename, bbox_inches="tight")
    plt.savefig(filename + ".pdf", bbox_inches="tight")

    plt.clf()

    return


