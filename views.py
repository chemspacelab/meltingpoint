
from qml.kernels import kpca
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import NullFormatter

from matplotlib import ticker
from scipy.stats import gaussian_kde

def set_border(ax, xkeys, ykeys,
    border=[False, False, True, True]):
    """
    Set border CSS style
    """

    ax.set_xticks(xkeys)
    ax.set_yticks(ykeys)

    spines = ax.spines.items()

    for direction, spine in spines:

        # spine.set_linewidth(1.2)

        if direction == "top":
            spine.set_visible(border[0])

        if direction == "right":
            spine.set_visible(border[1])

        if direction == "bottom":
            spine.set_visible(border[2])

            if border[2]:
                spine.set_bounds(min(xkeys), max(xkeys))

        if direction == "left":
            spine.set_visible(border[3])

            if border[3]:
                spine.set_bounds(min(ykeys), max(ykeys))

    return



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
        xlabel="Heavy atoms",
        ylabel="Kelvin",
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

    # Hack to make MPL hide the overlap of hexacons
    lineswidth=0.0 # white lines
    lineswidth=0.2 # perfect fit
    lineswidth=0.3 # fit for pngs
    lineswidth=0.4 # fit for pngs

    colormap = 'Greys'
    colormap = 'PuRd'

    hex_density = 50

    hexbinpar = {
        'gridsize': hex_density,
        'cmap': colormap,
        'linewidths': lineswidth,
        'mincnt': 1,
        'bins': 'log',
    }

    hb = ax_scatter.hexbin(xvalues, yvalues, **hexbinpar)

    # define binwidth
    x_max = np.max(xvalues)
    x_min = np.min(xvalues)
    x_binwidth = (abs(x_min) + x_max) / 30.0
    x_binwidth = int(x_binwidth)
    x_binwidth = 1.0
    x_bins = np.arange(x_min, x_max+x_binwidth, x_binwidth)

    y_max = np.max(yvalues)
    y_min = np.min(yvalues)
    y_binwidth = (abs(y_min) + y_max) / 50.0
    y_binwidth = int(y_binwidth)
    y_bins = np.arange(y_min, y_max+y_binwidth, y_binwidth)


    # Set limits and ticks of scatter
    xlim = (x_min-x_binwidth*2, x_max+x_binwidth*2)
    ylim = (0-y_binwidth*2, y_max+y_binwidth*2)
    ax_scatter.set_xlim(xlim)
    ax_scatter.set_ylim(ylim)

    xkeys = np.arange(10, x_max+x_binwidth*2, 10)
    xkeys = [1] + list(xkeys)
    ykeys = np.arange(0, y_max+y_binwidth, 100)

    # Histogram

    if True:
        bins = np.linspace(min(xkeys), max(xkeys), 300)
        gaussian_kernel = gaussian_kde(xvalues)
        values = gaussian_kernel(bins)
        ax_histx.plot(bins, values, "k", linewidth=1.0)

        bins = np.linspace(min(ykeys), max(ykeys), 300)
        gaussian_kernel = gaussian_kde(yvalues)
        values = gaussian_kernel(bins)
        ax_histy.plot(values, bins, "k", linewidth=1.0)

    else:
        ax_histx.hist(xvalues, bins=x_bins, histtype='step', color="k")
        ax_histy.hist(yvalues, bins=y_bins, orientation='horizontal', histtype='step', color="k")

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



def formatter_int(x, pos):
    return "%i" % x

def formatter_float(x, pos):
    return "%4.2f" % x

def formatter_off(x, pos):
    return ""

def formatter_notrail(x, pos):
    """
    remove trailing zeros
    """

    if x.is_integer():
        return formatter_int(x, pos)
    else:
        return '{0:g}'.format(x)

    return


def hex2color(s):
    """
    Function from MPL lib.

    Take a hex string *s* and return the corresponding rgb 3-tuple
    Example: #efefef -> (0.93725, 0.93725, 0.93725)
    """
    hexColorPattern = re.compile("\A#[a-fA-F0-9]{6}\Z")
    # if not isinstance(s, basestring):
    #     raise TypeError('hex2color requires a string argument')
    if hexColorPattern.match(s) is None:
        raise ValueError('invalid hex color string "%s"' % s)
    return tuple([int(n, 16)/255.0 for n in (s[1:3], s[3:5], s[5:7])])


def set_global_font(font="Fira Sans", fontsize=14):
    """
    """

    plt.rc('legend', fontsize=fontsize)
    mpl.rcParams['font.sans-serif'] = font
    mpl.rcParams['font.family'] = "sans-serif"
    mpl.rcParams['font.weight'] = "medium" # font.weight         : medium


    mpl.rcParams.update({'figure.autolayout': True})

    return


def set_custom_color():

    from cycler import cycler
    mpl.colors.ColorConverter.colors['r'] = hex2color('#e41a1c')
    mpl.colors.ColorConverter.colors['b'] = hex2color('#377eb8')
    mpl.colors.ColorConverter.colors['g'] = hex2color('#4daf4a')
    mpl.colors.ColorConverter.colors['p'] = hex2color('#984ea3')
    mpl.colors.ColorConverter.colors['y'] = hex2color('#ff7f00')

    return


def set_custom_lines():

    plt.rc('lines', antialiased=True) # render lines in antialised (no jaggies)
    plt.rc('xtick.major', size=6)      # major tick size in points
    plt.rc('xtick.minor', size=6)      # minor tick size in points
    plt.rc('xtick.major', pad=6)       # distance to major tick label in points
    plt.rc('xtick.minor', pad=6)       # distance to the minor tick label in points
    # plt.rc('xtick', color='111111')    # color of the tick labels
    # plt.rc('xtick', direction='out')    # direction: in or out
    #
    plt.rc('ytick.major', size=6)      # major tick size in points
    plt.rc('ytick.minor', size=6)      # minor tick size in points
    plt.rc('ytick.major', pad=6)       # distance to major tick label in points
    plt.rc('ytick.minor', pad=6)       # distance to the minor tick label in points
    # plt.rc('ytick', color='111111')    # color of the tick labels
    # plt.rc('ytick', direction='in')    # direction: in or out

    return


def plot_learning_curve(ax, xkeys, ykeys,
                        border=[True, False, True, False],
                        loglog=True,
                        show_legend=True):
    """

    """

    if loglog:
        ax.set_xscale('log')
        ax.set_yscale('log')

    if show_legend:
        leg = ax.legend(loc="best", borderaxespad=0., framealpha=1.0, fancybox=False, borderpad=1)
        leg.get_frame().set_linewidth(0.0)
        leg.get_frame().set_facecolor('#ffffff')

    ax.yaxis.grid(True, zorder=0)

    if border:
        # I like the css standard
        spines = ax.spines.items()
        for direction, spine in spines:
            if direction == "top": spine.set_visible(border[0])
            if direction == "right": spine.set_visible(border[1])
            if direction == "bottom": spine.set_visible(border[2])
            if direction == "left": spine.set_visible(border[3])


        # spines[0].set_visible(False) # left
        # spines[0].set_visible(border[3]) # left
        # spines[1].set_visible(border[1]) # right
        # spines[2].set_visible(border[2]) # bottom
        # spines[3].set_visible(border[0]) # top

    ax.set_xticks(xkeys)
    ax.set_xlim((min(xkeys)*(1-0.1), max(xkeys)*(1+0.1)))

    ax.set_yticks(ykeys)
    ax.set_ylim((min(ykeys), max(ykeys)))

    ax.yaxis.set_major_formatter(ticker.FuncFormatter(yformatter))
    ax.yaxis.set_minor_formatter(ticker.FuncFormatter(off_formatter))

    ax.xaxis.set_major_formatter(ticker.FuncFormatter(xformatter))
    ax.xaxis.set_minor_formatter(ticker.FuncFormatter(off_formatter))

    return None


def learning_curve_error(ax, xkeys, ykeys,
                         x_range=None,
                         y_range=None,
                         border=[False, False, True, True],
                         loglog=True,
                         show_legend=True):
    """

    """

    if loglog:
        ax.set_xscale('log')
        ax.set_yscale('log')


    ax.set_xticks(xkeys)
    ax.set_yticks(ykeys)

    if x_range is None:
        ax.set_xlim((min(xkeys)*(1-0.1), max(xkeys)*(1+0.1)))
    else:
        ax.set_xlim(tuple(x_range))

    if y_range is None:
        ax.set_ylim((min(ykeys), max(ykeys)))
    else:
        ax.set_ylim(tuple(y_range))


    ax.yaxis.set_major_formatter(ticker.FuncFormatter(formatter_int))
    ax.yaxis.set_minor_formatter(ticker.FuncFormatter(formatter_off))

    ax.xaxis.set_major_formatter(ticker.FuncFormatter(formatter_int))
    ax.xaxis.set_minor_formatter(ticker.FuncFormatter(formatter_off))


    # ax.xaxis.set_tick_params(width=1.2)
    # ax.yaxis.set_tick_params(width=1.2)

    if border:
        # I like the css standard
        spines = ax.spines.items()
        for direction, spine in spines:

            # spine.set_linewidth(1.2)

            if direction == "top":
                spine.set_visible(border[0])

            if direction == "right":
                spine.set_visible(border[1])

            if direction == "bottom":
                spine.set_visible(border[2])
                spine.set_bounds(min(xkeys), max(xkeys))

            if direction == "left":
                spine.set_visible(border[3])
                spine.set_bounds(min(ykeys), max(ykeys))


    # Remove the small minor ticks in log-log plot
    ax.minorticks_off()

    return



def legend_colorcoded(ax, lines, names):

    from matplotlib.legend_handler import HandlerLine2D
    import matplotlib.lines as mlines

    legkwargs = {
        "ncol": 6,
        "frameon": False,
        "columnspacing": 0.5,
        "handletextpad": -0.5,
        "handlelength": 1,
        # "markerfirst":False,
        # "markerscale": 0,
    }

    handles = []
    for line, name in zip(lines, names):
        color = plt.getp(line[0], 'color')
        handle = mlines.Line2D([], [], color=color, marker='None', linestyle='None',
            markersize=4, label=name.upper())
        handles.append(handle)

    # leg = ax.legend(**legkwargs)
    leg = ax.legend(handles=handles, **legkwargs)

    for text, line in zip(leg.get_texts(), lines):
        color = plt.getp(line[0], 'color')
        plt.setp(text, color=color)

    # leg.get_frame().set_edgecolor('b') # change color
    # leg.get_frame().set_linewidth(0.0) # remove box


    return


