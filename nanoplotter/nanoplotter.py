# wdecoster
'''
This module provides functions for plotting data extracted from Oxford Nanopore sequencing
reads and alignments, but some of it's functions can also be used for other applications.


FUNCTIONS
* Check if a specified color is a valid matplotlib color
checkvalidColor(color)
* Check if a specified output format is valid
checkvalidFormat(format)
* Create a bivariate plot with dots, hexbins and/or kernel density estimates.
Also arguments for specifying axis names, color and xlim/ylim
scatter(x, y, names, path, color, format, plots, stat=None, log=False, minvalx=0, minvaly=0)
* Create cumulative yield plot and evaluate read length and quality over time
timePlots(df, path, color, format)
* Create length distribution histogram and density curve
lengthPlots(array, name, path, n50, color, format, log=False)
* Create flowcell physical layout in numpy array
makeLayout()
* Present the activity (number of reads) per channel on the flowcell as a heatmap
spatialHeatmap(array, title, path, color, format)

'''


from __future__ import division
import logging
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import seaborn as sns


def check_valid_color(color):
    '''
    Check if the color provided by the user is valid
    If color is invalid the default is returned.
    '''
    if color in list(mcolors.CSS4_COLORS.keys()) + ["#4CB391"]:
        logging.info("Nanoplotter: Valid color {}.".format(color))
        return color
    else:
        logging.info("Nanoplotter: Invalid color {}, using default.".format(color))
        sys.stderr.write("Invalid color {}, using default.\n".format(color))
        return "#4CB391"


def check_valid_format(figformat):
    '''
    Check if the specified format is valid.
    If format is invalid the default is returned.
    Probably installation-dependent
    '''
    fig = plt.figure()
    if figformat in list(fig.canvas.get_supported_filetypes().keys()):
        logging.info("Nanoplotter: valid output format {}".format(figformat))
        return figformat
    else:
        logging.info("Nanoplotter: invalid output format {}".format(figformat))
        sys.stderr.write("Invalid format {}, using default.\n".format(figformat))
        return "png"


def scatter(x, y, names, path, color, figformat, plots, stat=None, log=False, minvalx=0, minvaly=0):
    '''
    Plotting functionq
    Create three types of joint plots of length vs quality, containing marginal summaries
    -A scatter plot with histograms on axes
    -A hexagonal binned plot with histograms on axes
    -A kernel density plot with density curves on axes, subsampled to 10000 reads if required
    '''
    logging.info("Nanoplotter: Creating {} vs {} plots using statistics from {} reads.".format(
        names[0], names[1], x.size))
    sns.set(style="ticks")
    maxvalx = np.amax(x)
    maxvaly = np.amax(y)
    if plots["hex"]:
        plot = sns.jointplot(
            x=x,
            y=y,
            kind="hex",
            color=color,
            stat_func=stat,
            space=0,
            xlim=(minvalx, maxvalx),
            ylim=(minvaly, maxvaly),
            size=10)
        plot.set_axis_labels(names[0], names[1])
        if log:
            ticks = [10**i for i in range(10) if not 10**i > 10 * (10**maxvalx)]
            plot.ax_joint.set_xticks(np.log10(ticks))
            plot.ax_joint.set_xticklabels(ticks)
        plot.savefig(path + "_hex." + figformat, format=figformat, dpi=100)
    sns.set(style="darkgrid")
    if plots["dot"]:
        plot = sns.jointplot(
            x=x,
            y=y,
            kind="scatter",
            color=color,
            stat_func=stat,
            xlim=(minvalx, maxvalx),
            ylim=(minvaly, maxvaly),
            space=0,
            size=10,
            joint_kws={"s": 1})
        plot.set_axis_labels(names[0], names[1])
        if log:
            ticks = [10**i for i in range(10) if not 10**i > 10 * (10**maxvalx)]
            plot.ax_joint.set_xticks(np.log10(ticks))
            plot.ax_joint.set_xticklabels(ticks)
        plot.savefig(path + "_dot." + figformat, format=figformat, dpi=100)
    if plots["kde"]:
        plot = sns.jointplot(
            x=x,
            y=y,
            kind="kde",
            clip=((0, np.Inf), (0, np.Inf)),
            xlim=(minvalx, maxvalx),
            ylim=(minvaly, maxvaly),
            space=0,
            color=color,
            stat_func=stat,
            shade_lowest=False,
            size=10)
        plot.set_axis_labels(names[0], names[1])
        if log:
            ticks = [10**i for i in range(10) if not 10**i > 10 * (10**maxvalx)]
            plot.ax_joint.set_xticks(np.log10(ticks))
            plot.ax_joint.set_xticklabels(ticks)
        plot.savefig(path + "_kde." + figformat, format=figformat, dpi=1000)
    plt.close("all")


def check_valid_time(times):
    '''
    Check if the data contains reads created within the same 96-hours timeframe
    if not, return false and warn the user that time plots are invalid and not created
    '''
    timediff = (times.max() - times.min()).days
    if timediff < 7:
        return True
    else:
        sys.stderr.write("\nWarning: data generated is from more than 7 days.\n")
        sys.stderr.write("Likely this indicates you are combining multiple runs.\n")
        sys.stderr.write("As such plots based on time are invalid and therefore skipped.\n\n")
        logging.warning("Time plots not created: invalid timespan: {} days".format(str(timediff)))
        return False


def time_plots(df, path, color, figformat):
    '''
    Plotting function
    Making plots of time vs read length, time vs quality and cumulative yield
    '''
    if check_valid_time(df["start_time"]):
        logging.info("Nanoplotter: Creating timeplots.")
        dfs = df.sort_values("start_time")
        dfs["cumyield_gb"] = dfs["lengths"].cumsum() / 10**9
        dfs_sparse = dfs.sample(min(2000, len(df.index)))
        dfs_sparse["time"] = dfs_sparse["start_time"].astype('timedelta64[s]')
        maxtime = dfs_sparse.time.max()
        ticks = [int(i) for i in range(0, 168, 4) if not i > (maxtime / 3600)]

        if "quals" in df:
            g = sns.JointGrid(
                x='time',
                y='quals',
                data=dfs_sparse,
                space=0,
                size=10,
                xlim=(0, maxtime))
            g.plot_joint(plt.scatter, color=color)
            g.ax_joint.set_xticks([i * 3600 for i in ticks])
            g.ax_joint.set_xticklabels(ticks)
            g.ax_marg_y.hist(dfs_sparse['quals'].dropna(), orientation="horizontal", color=color)
            g.set_axis_labels('Run time (hours)', 'Median average basecall quality')
            g.savefig(path + "TimeQualityScatterPlot." + figformat, format=figformat, dpi=100)

        g = sns.JointGrid(
            x='time',
            y="lengths",
            data=dfs_sparse,
            space=0,
            size=10,
            xlim=(0, maxtime))
        g.plot_joint(plt.scatter, color=color)
        g.ax_joint.set_xticks([i * 3600 for i in ticks])
        g.ax_joint.set_xticklabels(ticks)
        g.ax_marg_y.hist(dfs_sparse["lengths"].dropna(), orientation="horizontal", color=color)
        g.set_axis_labels('Run time (hours)', 'Median read length')
        g.savefig(path + "TimeLengthScatterPlot." + figformat, format=figformat, dpi=100)
        plt.close("all")

        ax = sns.regplot(
            x='time',
            y="cumyield_gb",
            data=dfs_sparse,
            x_ci=None,
            fit_reg=False,
            color=color,
            scatter_kws={"s": 5})
        ax.set(
            xticks=[i * 3600 for i in ticks],
            xticklabels=ticks,
            xlabel='Run time (hours)',
            ylabel='Cumulative yield in gigabase')
        fig = ax.get_figure()
        fig.savefig(path + "CumulativeYieldPlot." + figformat, format=figformat, dpi=100)
        plt.close("all")


def length_plots(array, name, path, n50, color, figformat, log=False):
    '''
    Plotting function
    Create density plot and histogram based on a numpy array
    containing read lengths or transformed read lengths
    '''
    logging.info("Nanoplotter: Creating length plots for {}.".format(name))
    logging.info("Nanoplotter: Using {} reads with read length N50 of {}.".format(array.size, n50))
    maxvalx = np.amax(array)

    ax = sns.distplot(
        a=array,
        kde=True,
        hist=False,
        bins=round(int(maxvalx) / 500),
        color=color,
        kde_kws={"label": name, "clip": (0, maxvalx)})
    if log:
        ticks = [10**i for i in range(10) if not 10**i > 10 * (10**maxvalx)]
        ax.set(xticks=np.log10(ticks), xticklabels=ticks)
    fig = ax.get_figure()
    fig.savefig(path + "DensityCurve" + name.replace(' ', '') +
                "." + figformat, format=figformat, dpi=100)
    plt.close("all")

    ax = sns.distplot(
        a=array,
        kde=False,
        hist=True,
        bins=round(int(maxvalx) / 500),
        color=color)
    if log:
        ax.set(
            xticks=np.log10(ticks),
            xticklabels=ticks)
        plt.axvline(np.log10(n50))
        plt.annotate('N50', xy=(np.log10(n50), np.amax(
            [h.get_height() for h in ax.patches])), size=8)
    else:
        plt.axvline(n50)
        plt.annotate('N50', xy=(n50, np.amax([h.get_height() for h in ax.patches])), size=8)
    ax.set(xlabel='Read length', ylabel='Number of reads')
    fig = ax.get_figure()
    fig.savefig(path + "Histogram" + name.replace(' ', '') +
                "." + figformat, format=figformat, dpi=100)
    plt.close("all")


def make_layout():
    '''
    Make the physical layout of the MinION flowcell
    based on https://bioinformatics.stackexchange.com/a/749/681
    returned as a numpy array
    '''
    layoutlist = []
    for i, j in zip([33, 481, 417, 353, 289, 225, 161, 97], [8, 456, 392, 328, 264, 200, 136, 72]):
        for n in range(4):
            layoutlist.append(list(range(i + n * 8, (i + n * 8) + 8, 1)) +
                              list(range(j + n * 8, (j + n * 8) - 8, -1)))
    return np.array(layoutlist).transpose()


def spatial_heatmap(array, title, path, color, figformat):
    '''
    Plotting function
    Taking channel information and creating post run channel activity plots
    '''
    logging.info("Nanoplotter: Creating activity map for {} using statistics from {} reads.".format(
        title.lower(), array.size))
    layout = make_layout()
    activityData = np.zeros((16, 32))
    valueCounts = pd.value_counts(pd.Series(array))
    for entry in valueCounts.keys():
        activityData[np.where(layout == entry)] = valueCounts[entry]

    plt.figure()
    ax = sns.heatmap(
        data=activityData,
        xticklabels=range(1, 33),
        yticklabels=range(1, 17),
        square=True,
        cbar_kws={"orientation": "horizontal"},
        cmap=color,
        linewidths=0.20)
    ax.set_title(title)
    fig = ax.get_figure()
    fig.savefig(path + "." + figformat, format=figformat, dpi=100)
    plt.close("all")


def violin_or_box_plot(df, y, figformat, path, violin=True, log=False):
    '''
    Plotting function
    Create a violinplot from the received DataFrame
    The x-axis should be divided based on the 'dataset' column,
    the y-axis is specified in the arguments
    '''
    if violin:
        logging.info("Nanoplotter: Creating violin plot for {}.".format(y))
        ax = sns.violinplot(x="dataset", y=y, data=df, inner=None, cut=0)
    else:
        logging.info("Nanoplotter: Creating box plot for {}.".format(y))
        ax = sns.boxplot(x="dataset", y=y, data=df)
    if log:
        ticks = [10**i for i in range(10) if not 10**i > 10 * (10**np.amax(df[y]))]
        ax.set(yticks=np.log10(ticks), yticklabels=ticks)
    fig = ax.get_figure()
    fig.savefig(
        path + "NanoComp_" + y.replace(' ', '_') + '.' + figformat,
        format=figformat,
        dpi=100)
    plt.close("all")


checkvalidColor = check_valid_color
checkvalidFormat = check_valid_format
spatialHeatmap = spatial_heatmap
lengthPlots = length_plots
timePlots = time_plots
