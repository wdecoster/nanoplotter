# wdecoster
"""
This module provides functions for plotting data extracted from Oxford Nanopore sequencing
reads and alignments, but some of it's functions can also be used for other applications.


FUNCTIONS
* Check if a specified color is a valid matplotlib color
check_valid_color(color)
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

"""


import logging
import sys
from datetime import timedelta
import pandas as pd
import numpy as np
import base64
from math import ceil
import io
import urllib
from collections import namedtuple
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import seaborn as sns
from pauvre.marginplot import margin_plot
import plotly
import plotly.graph_objs as go


class Plot(object):
    """A Plot object is defined by a path to the output file and the title of the plot."""

    def __init__(self, path, title):
        self.path = path
        self.title = title
        self.fig = None
        self.html = None

    def encode(self):
        if self.html:
            return self.html
        elif self.fig:
            return self.encode2()
        else:
            return self.encode1()

    def encode1(self):
        """Return the base64 encoding of the figure file and insert in html image tag."""
        data_uri = base64.b64encode(open(self.path, 'rb').read()).decode('utf-8').replace('\n', '')
        return '<img src="data:image/png;base64,{0}">'.format(data_uri)

    def encode2(self):
        """Return the base64 encoding of the fig attribute and insert in html image tag."""
        buf = io.BytesIO()
        self.fig.savefig(buf, format='png', bbox_inches='tight', dpi=100)
        buf.seek(0)
        string = base64.b64encode(buf.read())
        return '<img src="data:image/png;base64,{0}">'.format(urllib.parse.quote(string))


class Layout(object):
    def __init__(self, structure, template, xticks, yticks):
        self.structure = structure
        self.template = template
        self.xticks = xticks
        self.yticks = yticks


def check_valid_color(color):
    """Check if the color provided by the user is valid.

    If color is invalid the default is returned.
    """
    if color in list(mcolors.CSS4_COLORS.keys()) + ["#4CB391"]:
        logging.info("Nanoplotter: Valid color {}.".format(color))
        return color
    else:
        logging.info("Nanoplotter: Invalid color {}, using default.".format(color))
        sys.stderr.write("Invalid color {}, using default.\n".format(color))
        return "#4CB391"


def check_valid_format(figformat):
    """Check if the specified figure format is valid.

    If format is invalid the default is returned.
    Probably installation-dependent
    """
    fig = plt.figure()
    if figformat in list(fig.canvas.get_supported_filetypes().keys()):
        logging.info("Nanoplotter: valid output format {}".format(figformat))
        return figformat
    else:
        logging.info("Nanoplotter: invalid output format {}".format(figformat))
        sys.stderr.write("Invalid format {}, using default.\n".format(figformat))
        return "png"


def plot_settings(plot_settings, dpi):
    sns.set(**plot_settings)
    mpl.rcParams['savefig.dpi'] = dpi


def scatter(x, y, names, path, plots, color="#4CB391", figformat="png",
            stat=None, log=False, minvalx=0, minvaly=0, title=None, plot_settings=None):
    """Create bivariate plots.

    Create four types of bivariate plots of x vs y, containing marginal summaries
    -A scatter plot with histograms on axes
    -A hexagonal binned plot with histograms on axes
    -A kernel density plot with density curves on axes
    -A pauvre-style plot using code from https://github.com/conchoecia/pauvre
    """
    logging.info("Nanoplotter: Creating {} vs {} plots using statistics from {} reads.".format(
        names[0], names[1], x.size))
    sns.set(style="ticks", **plot_settings)
    maxvalx = np.amax(x)
    maxvaly = np.amax(y)

    plots_made = []

    if plots["hex"]:
        hex_plot = Plot(
            path=path + "_hex." + figformat,
            title="{} vs {} plot using hexagonal bins".format(names[0], names[1]))
        plot = sns.jointplot(
            x=x,
            y=y,
            kind="hex",
            color=color,
            stat_func=stat,
            space=0,
            xlim=(minvalx, maxvalx),
            ylim=(minvaly, maxvaly),
            height=10)
        plot.set_axis_labels(names[0], names[1])
        if log:
            hex_plot.title = hex_plot.title + " after log transformation of read lengths"
            ticks = [10**i for i in range(10) if not 10**i > 10 * (10**maxvalx)]
            plot.ax_joint.set_xticks(np.log10(ticks))
            plot.ax_marg_x.set_xticks(np.log10(ticks))
            plot.ax_joint.set_xticklabels(ticks)
        plt.subplots_adjust(top=0.90)
        plot.fig.suptitle(title or "{} vs {} plot".format(names[0], names[1]), fontsize=25)
        hex_plot.fig = plot
        plot.savefig(hex_plot.path, format=figformat, bbox_inches="tight")
        plots_made.append(hex_plot)

    sns.set(style="darkgrid", **plot_settings)
    if plots["dot"]:
        dot_plot = Plot(
            path=path + "_dot." + figformat,
            title="{} vs {} plot using dots".format(names[0], names[1]))
        plot = sns.jointplot(
            x=x,
            y=y,
            kind="scatter",
            color=color,
            stat_func=stat,
            xlim=(minvalx, maxvalx),
            ylim=(minvaly, maxvaly),
            space=0,
            height=10,
            joint_kws={"s": 1})
        plot.set_axis_labels(names[0], names[1])
        if log:
            dot_plot.title = dot_plot.title + " after log transformation of read lengths"
            ticks = [10**i for i in range(10) if not 10**i > 10 * (10**maxvalx)]
            plot.ax_joint.set_xticks(np.log10(ticks))
            plot.ax_marg_x.set_xticks(np.log10(ticks))
            plot.ax_joint.set_xticklabels(ticks)
        plt.subplots_adjust(top=0.90)
        plot.fig.suptitle(title or "{} vs {} plot".format(names[0], names[1]), fontsize=25)
        dot_plot.fig = plot
        plot.savefig(dot_plot.path, format=figformat, bbox_inches="tight")
        plots_made.append(dot_plot)

    if plots["kde"]:
        idx = np.random.choice(x.index, min(2000, len(x)), replace=False)
        kde_plot = Plot(
            path=path + "_kde." + figformat,
            title="{} vs {} plot using a kernel density estimation".format(names[0], names[1]))
        plot = sns.jointplot(
            x=x[idx],
            y=y[idx],
            kind="kde",
            clip=((0, np.Inf), (0, np.Inf)),
            xlim=(minvalx, maxvalx),
            ylim=(minvaly, maxvaly),
            space=0,
            color=color,
            stat_func=stat,
            shade_lowest=False,
            height=10)
        plot.set_axis_labels(names[0], names[1])
        if log:
            kde_plot.title = kde_plot.title + " after log transformation of read lengths"
            ticks = [10**i for i in range(10) if not 10**i > 10 * (10**maxvalx)]
            plot.ax_joint.set_xticks(np.log10(ticks))
            plot.ax_marg_x.set_xticks(np.log10(ticks))
            plot.ax_joint.set_xticklabels(ticks)
        plt.subplots_adjust(top=0.90)
        plot.fig.suptitle(title or "{} vs {} plot".format(names[0], names[1]), fontsize=25)
        kde_plot.fig = plot
        plot.savefig(kde_plot.path, format=figformat, bbox_inches="tight")
        plots_made.append(kde_plot)

    if plots["pauvre"] and names == ['Read lengths', 'Average read quality'] and log is False:
        pauvre_plot = Plot(
            path=path + "_pauvre." + figformat,
            title="{} vs {} plot using pauvre-style @conchoecia".format(names[0], names[1]))
        sns.set(style="white", **plot_settings)
        margin_plot(df=pd.DataFrame({"length": x, "meanQual": y}),
                    Y_AXES=False,
                    title=title or "Length vs Quality in Pauvre-style",
                    plot_maxlen=None,
                    plot_minlen=0,
                    plot_maxqual=None,
                    plot_minqual=0,
                    lengthbin=None,
                    qualbin=None,
                    BASENAME="whatever",
                    path=pauvre_plot.path,
                    fileform=[figformat],
                    dpi=600,
                    TRANSPARENT=True,
                    QUIET=True)
        plots_made.append(pauvre_plot)
    plt.close("all")
    return plots_made


def check_valid_time_and_sort(df, timescol, days=5, warning=True):
    """Check if the data contains reads created within the same `days` timeframe.

    if not, print warning and only return part of the data which is within `days` days
    Resetting the index twice to get also an "index" column for plotting the cum_yield_reads plot
    """
    timediff = (df[timescol].max() - df[timescol].min()).days
    if timediff < days:
        return df.sort_values(timescol).reset_index(drop=True).reset_index()
    else:
        if warning:
            sys.stderr.write(
                "\nWarning: data generated is from more than {} days.\n".format(str(days)))
            sys.stderr.write("Likely this indicates you are combining multiple runs.\n")
            sys.stderr.write(
                "Plots based on time are invalid and therefore truncated to first {} days.\n\n"
                .format(str(days)))
            logging.warning("Time plots truncated to first {} days: invalid timespan: {} days"
                            .format(str(days), str(timediff)))
        return df[df[timescol] < timedelta(days=days)] \
            .sort_values(timescol) \
            .reset_index(drop=True) \
            .reset_index()


def time_plots(df, path, title=None, color="#4CB391", figformat="png",
               log_length=False, plot_settings=None):
    """Making plots of time vs read length, time vs quality and cumulative yield."""
    dfs = check_valid_time_and_sort(df, "start_time")
    logging.info("Nanoplotter: Creating timeplots using {} reads.".format(len(dfs)))
    cumyields = cumulative_yield(dfs=dfs.set_index("start_time"),
                                 path=path,
                                 figformat=figformat,
                                 title=title,
                                 color=color)
    reads_pores_over_time = plot_over_time(dfs=dfs.set_index("start_time"),
                                           path=path,
                                           figformat=figformat,
                                           title=title,
                                           color=color)
    violins = violin_plots_over_time(dfs=dfs,
                                     path=path,
                                     figformat=figformat,
                                     title=title,
                                     log_length=log_length,
                                     plot_settings=plot_settings)
    return cumyields + reads_pores_over_time + violins


def violin_plots_over_time(dfs, path, figformat, title,
                           log_length=False, plot_settings=None):
    maxtime = dfs["start_time"].max().total_seconds()
    time_length = Plot(
        path=path + "TimeLengthViolinPlot." + figformat,
        title="Violin plot of read lengths over time")
    sns.set(style="white", **plot_settings)
    labels = [str(i) + "-" + str(i + 3) for i in range(0, 168, 3) if not i > (maxtime / 3600)]
    if log_length:
        length_column = "log_lengths"
    else:
        length_column = "lengths"

    dfs['timebin'] = pd.cut(
        x=dfs["start_time"],
        bins=ceil((maxtime / 3600) / 3),
        labels=labels)
    ax = sns.violinplot(
        x="timebin",
        y=length_column,
        data=dfs,
        inner=None,
        cut=0,
        linewidth=0)
    ax.set(
        xlabel='Interval (hours)',
        ylabel="Read length",
        title=title or time_length.title)
    if log_length:
        ticks = [10**i for i in range(10) if not 10**i > 10 * np.amax(dfs["lengths"])]
        ax.set(
            yticks=np.log10(ticks),
            yticklabels=ticks)
    plt.xticks(rotation=45, ha='center', fontsize=8)
    fig = ax.get_figure()
    time_length.fig = fig
    fig.savefig(
        fname=time_length.path,
        format=figformat,

        bbox_inches='tight')
    plt.close("all")

    plots = [time_length]

    if "quals" in dfs:
        time_qual = Plot(
            path=path + "TimeQualityViolinPlot." + figformat,
            title="Violin plot of quality over time")
        sns.set(style="white", **plot_settings)
        ax = sns.violinplot(
            x="timebin",
            y="quals",
            data=dfs,
            inner=None,
            cut=0,
            linewidth=0)
        ax.set(
            xlabel='Interval (hours)',
            ylabel="Basecall quality",
            title=title or time_qual.title)
        plt.xticks(rotation=45, ha='center', fontsize=8)
        fig = ax.get_figure()
        time_qual.fig = fig
        fig.savefig(
            fname=time_qual.path,
            format=figformat,
            bbox_inches='tight')
        plots.append(time_qual)
        plt.close("all")

    if "duration" in dfs:
        time_duration = Plot(
            path=path + "TimeSequencingSpeed_ViolinPlot." + figformat,
            title="Violin plot of sequencing speed over time")
        labels = [str(i) + "-" + str(i + 3) for i in range(0, 168, 3) if not i > (maxtime / 3600)]
        dfs['timebin'] = pd.cut(
            x=dfs["start_time"],
            bins=ceil((maxtime / 3600) / 3),
            labels=labels)
        ax = sns.violinplot(
            x=dfs["timebin"],
            y=dfs["lengths"] / dfs["duration"],
            inner=None,
            cut=0,
            linewidth=0)
        ax.set(
            xlabel='Interval (hours)',
            ylabel="Sequencing speed (nucleotides/second)",
            title=title or time_duration.title)
        plt.xticks(rotation=45, ha='center', fontsize=8)
        fig = ax.get_figure()
        time_duration.fig = fig
        fig.savefig(
            fname=time_duration.path,
            format=figformat,
            bbox_inches='tight')
        plots.append(time_duration)
        plt.close("all")
    return plots


def plot_over_time(dfs, path, figformat, title, color):
    num_reads = Plot(
        path=path + "NumberOfReads_Over_Time." + figformat,
        title="Number of reads over time")
    s = dfs.loc[:, "lengths"].resample('10T').count()
    ax = sns.regplot(
        x=s.index.total_seconds() / 3600,
        y=s,
        x_ci=None,
        fit_reg=False,
        color=color,
        scatter_kws={"s": 3})
    ax.set(
        xlabel='Run time (hours)',
        ylabel='Number of reads per 10 minutes',
        title=title or num_reads.title)
    fig = ax.get_figure()
    num_reads.fig = fig
    fig.savefig(num_reads.path, format=figformat, bbox_inches="tight")
    plt.close("all")
    plots = [num_reads]

    if "channelIDs" in dfs:
        pores_over_time = Plot(
            path=path + "ActivePores_Over_Time." + figformat,
            title="Number of active pores over time")
        s = dfs.loc[:, "channelIDs"].resample('10T').nunique()
        ax = sns.regplot(
            x=s.index.total_seconds() / 3600,
            y=s,
            x_ci=None,
            fit_reg=False,
            color=color,
            scatter_kws={"s": 3})
        ax.set(
            xlabel='Run time (hours)',
            ylabel='Active pores per 10 minutes',
            title=title or pores_over_time.title)
        fig = ax.get_figure()
        pores_over_time.fig = fig
        fig.savefig(pores_over_time.path, format=figformat, bbox_inches="tight")
        plt.close("all")
        plots.append(pores_over_time)
    return plots


def cumulative_yield(dfs, path, figformat, title, color):
    cum_yield_gb = Plot(
        path=path + "CumulativeYieldPlot_Gigabases." + figformat,
        title="Cumulative yield")
    s = dfs.loc[:, "lengths"].cumsum().resample('10T').max() / 1e9
    ax = sns.regplot(
        x=s.index.total_seconds() / 3600,
        y=s,
        x_ci=None,
        fit_reg=False,
        color=color,
        scatter_kws={"s": 3})
    ax.set(
        xlabel='Run time (hours)',
        ylabel='Cumulative yield in gigabase',
        title=title or cum_yield_gb.title)
    fig = ax.get_figure()
    cum_yield_gb.fig = fig
    fig.savefig(cum_yield_gb.path, format=figformat, bbox_inches="tight")
    plt.close("all")

    cum_yield_reads = Plot(
        path=path + "CumulativeYieldPlot_NumberOfReads." + figformat,
        title="Cumulative yield")
    s = dfs.loc[:, "lengths"].resample('10T').count().cumsum()
    ax = sns.regplot(
        x=s.index.total_seconds() / 3600,
        y=s,
        x_ci=None,
        fit_reg=False,
        color=color,
        scatter_kws={"s": 3})
    ax.set(
        xlabel='Run time (hours)',
        ylabel='Cumulative yield in number of reads',
        title=title or cum_yield_reads.title)
    fig = ax.get_figure()
    cum_yield_reads.fig = fig
    fig.savefig(cum_yield_reads.path, format=figformat, bbox_inches="tight")
    plt.close("all")
    return [cum_yield_gb, cum_yield_reads]


def length_plots(array, name, path, title=None, n50=None, color="#4CB391", figformat="png"):
    """Create histogram of normal and log transformed read lengths."""
    logging.info("Nanoplotter: Creating length plots for {}.".format(name))
    maxvalx = np.amax(array)
    if n50:
        logging.info("Nanoplotter: Using {} reads with read length N50 of {}bp and maximum of {}bp."
                     .format(array.size, n50, maxvalx))
    else:
        logging.info("Nanoplotter: Using {} reads maximum of {}bp.".format(array.size, maxvalx))

    HistType = namedtuple('HistType', 'weight name ylabel')
    plots = []
    for h_type in [HistType(None, "", "Number of reads"),
                   HistType(array, "Weighted ", "Number of bases")]:
        histogram = Plot(
            path=path + h_type.name.replace(" ", "_") + "Histogram" +
            name.replace(' ', '') + "." + figformat,
            title=h_type.name + "Histogram of read lengths")
        ax = sns.distplot(
            a=array,
            kde=False,
            hist=True,
            bins=max(round(int(maxvalx) / 500), 10),
            color=color,
            hist_kws=dict(weights=h_type.weight,
                          edgecolor=color,
                          linewidth=0.2,
                          alpha=0.8))
        if n50:
            plt.axvline(n50)
            plt.annotate('N50', xy=(n50, np.amax([h.get_height() for h in ax.patches])), size=8)
        ax.set(
            xlabel='Read length',
            ylabel=h_type.ylabel,
            title=title or histogram.title)
        plt.ticklabel_format(style='plain', axis='y')
        fig = ax.get_figure()
        histogram.fig = fig
        fig.savefig(histogram.path, format=figformat, bbox_inches="tight")
        plt.close("all")

        log_histogram = Plot(
            path=path + h_type.name.replace(" ", "_") + "LogTransformed_Histogram" +
            name.replace(' ', '') + "." + figformat,
            title=h_type.name + "Histogram of read lengths after log transformation")
        ax = sns.distplot(
            a=np.log10(array),
            kde=False,
            hist=True,
            color=color,
            hist_kws=dict(weights=h_type.weight,
                          edgecolor=color,
                          linewidth=0.2,
                          alpha=0.8))
        ticks = [10**i for i in range(10) if not 10**i > 10 * maxvalx]
        ax.set(
            xticks=np.log10(ticks),
            xticklabels=ticks,
            xlabel='Read length',
            ylabel=h_type.ylabel,
            title=title or log_histogram.title)
        if n50:
            plt.axvline(np.log10(n50))
            plt.annotate('N50', xy=(np.log10(n50), np.amax(
                [h.get_height() for h in ax.patches])), size=8)
        plt.ticklabel_format(style='plain', axis='y')
        fig = ax.get_figure()
        log_histogram.fig = fig
        fig.savefig(log_histogram.path, format=figformat, bbox_inches="tight")
        plt.close("all")
        plots.extend([histogram, log_histogram])
    plots.append(yield_by_minimal_length_plot(array, name, path, title=None,
                                              n50=None, color="#4CB391", figformat=figformat))
    return plots


def yield_by_minimal_length_plot(array, name, path,
                                 title=None, n50=None, color="#4CB391", figformat="png"):
    df = pd.DataFrame(data={"lengths": np.sort(array)[::-1]})
    df["cumyield_gb"] = df["lengths"].cumsum() / 10**9
    yield_by_length = Plot(
        path=path + "Yield_By_Length." + figformat,
        title="Yield by length")
    ax = sns.regplot(
        x='lengths',
        y="cumyield_gb",
        data=df,
        x_ci=None,
        fit_reg=False,
        color=color,
        scatter_kws={"s": 3})
    ax.set(
        xlabel='Read length',
        ylabel='Cumulative yield for minimal length',
        title=title or yield_by_length.title)
    fig = ax.get_figure()
    yield_by_length.fig = fig
    fig.savefig(yield_by_length.path, format=figformat, bbox_inches="tight")
    plt.close("all")
    return yield_by_length


def make_layout(maxval):
    """Make the physical layout of the MinION flowcell.
    based on https://bioinformatics.stackexchange.com/a/749/681
    returned as a numpy array
    """
    if maxval > 512:
        return Layout(
            structure=np.concatenate([np.array([list(range(10 * i + 1, i * 10 + 11))
                                                for i in range(25)]) + j
                                      for j in range(0, 3000, 250)],
                                     axis=1),
            template=np.zeros((25, 120)),
            xticks=range(1, 121),
            yticks=range(1, 26))
    else:
        layoutlist = []
        for i, j in zip(
                [33, 481, 417, 353, 289, 225, 161, 97],
                [8, 456, 392, 328, 264, 200, 136, 72]):
            for n in range(4):
                layoutlist.append(list(range(i + n * 8, (i + n * 8) + 8, 1)) +
                                  list(range(j + n * 8, (j + n * 8) - 8, -1)))
        return Layout(
            structure=np.array(layoutlist).transpose(),
            template=np.zeros((16, 32)),
            xticks=range(1, 33),
            yticks=range(1, 17))


def spatial_heatmap(array, path, title=None, color="Greens", figformat="png"):
    """Taking channel information and creating post run channel activity plots."""
    logging.info("Nanoplotter: Creating heatmap of reads per channel using {} reads."
                 .format(array.size))
    activity_map = Plot(
        path=path + "." + figformat,
        title="Number of reads generated per channel")
    layout = make_layout(maxval=np.amax(array))
    valueCounts = pd.value_counts(pd.Series(array))
    for entry in valueCounts.keys():
        layout.template[np.where(layout.structure == entry)] = valueCounts[entry]
    plt.figure()
    ax = sns.heatmap(
        data=pd.DataFrame(layout.template, index=layout.yticks, columns=layout.xticks),
        xticklabels="auto",
        yticklabels="auto",
        square=True,
        cbar_kws={"orientation": "horizontal"},
        cmap=color,
        linewidths=0.20)
    ax.set_title(title or activity_map.title)
    fig = ax.get_figure()
    activity_map.fig = fig
    fig.savefig(activity_map.path, format=figformat)
    plt.close("all")
    return [activity_map]


def violin_or_box_plot(df, y, figformat, path, y_name,
                       title=None, violin=True, log=False, palette=None):
    """Create a violin or boxplot from the received DataFrame.

    The x-axis should be divided based on the 'dataset' column,
    the y-axis is specified in the arguments
    """
    violin_comp = Plot(
        path=path + "NanoComp_" + y.replace(' ', '_') + '.' + figformat,
        title="Comparing {}".format(y))
    if y == "quals":
        violin_comp.title = "Comparing base call quality scores"
    if violin:
        logging.info("Nanoplotter: Creating violin plot for {}.".format(y))
        ax = sns.violinplot(
            x="dataset",
            y=y,
            data=df,
            inner=None,
            cut=0,
            palette=palette,
            linewidth=0)
    else:
        logging.info("Nanoplotter: Creating box plot for {}.".format(y))
        ax = sns.boxplot(
            x="dataset",
            y=y,
            data=df,
            palette=palette)
    if log:
        ticks = [10**i for i in range(10) if not 10**i > 10 * (10**np.amax(df[y]))]
        ax.set(
            yticks=np.log10(ticks),
            yticklabels=ticks)
    ax.set(title=title or violin_comp.title,
           ylabel=y_name)
    plt.xticks(rotation=30, ha='center')
    fig = ax.get_figure()
    violin_comp.fig = fig
    fig.savefig(
        fname=violin_comp.path,
        format=figformat,
        bbox_inches='tight')
    plt.close("all")
    return [violin_comp]


def output_barplot(df, figformat, path, title=None, palette=None):
    """Create barplots based on number of reads and total sum of nucleotides sequenced."""
    logging.info("Nanoplotter: Creating barplots for number of reads and total throughput.")
    read_count = Plot(
        path=path + "NanoComp_number_of_reads." + figformat,
        title="Comparing number of reads")
    ax = sns.countplot(
        x="dataset",
        data=df,
        palette=palette)
    ax.set(
        ylabel='Number of reads',
        title=title or read_count.title)
    plt.xticks(rotation=30, ha='center')
    fig = ax.get_figure()
    read_count.fig = fig
    fig.savefig(
        fname=read_count.path,
        format=figformat,
        bbox_inches='tight')
    plt.close("all")

    throughput_bases = Plot(
        path=path + "NanoComp_total_throughput." + figformat,
        title="Comparing throughput in gigabases sequenced")
    throughput = df.groupby('dataset')['lengths'].sum()
    ax = sns.barplot(
        x=list(throughput.index),
        y=throughput / 1e9,
        palette=palette,
        order=df["dataset"].unique())
    ax.set(
        ylabel='Total gigabase sequenced',
        title=title or throughput_bases.title)
    plt.xticks(rotation=30, ha='center')
    fig = ax.get_figure()
    throughput_bases.fig = fig
    fig.savefig(
        fname=throughput_bases.path,
        format=figformat,
        bbox_inches='tight')
    plt.close("all")
    return read_count, throughput_bases


def compare_cumulative_yields(df, path, palette=None, title=None):
    if palette is None:
        palette = plotly.colors.DEFAULT_PLOTLY_COLORS * 5
    dfs = check_valid_time_and_sort(df, "start_time").set_index("start_time")

    logging.info("Nanoplotter: Creating cumulative yield plots using {} reads.".format(len(dfs)))
    cum_yield_gb = Plot(
        path=path + "NanoComp_CumulativeYieldPlot_Gigabases.html",
        title="Cumulative yield")
    data = []
    for d, c in zip(df.dataset.unique(), palette):
        s = dfs.loc[dfs.dataset == d, "lengths"].cumsum().resample('10T').max() / 1e9
        data.append(go.Scatter(x=s.index.total_seconds() / 3600,
                               y=s,
                               opacity=0.75,
                               name=d,
                               marker=dict(color=c))
                    )
    cum_yield_gb.html = plotly.offline.plot({
        "data": data,
        "layout": go.Layout(barmode='overlay',
                            title=title or cum_yield_gb.title,
                            xaxis=dict(title="Time (hours)"),
                            yaxis=dict(title="Yield (gigabase)"),
                            )},
        output_type="div",
        show_link=False)
    with open(cum_yield_gb.path, 'w') as html_out:
        html_out.write(cum_yield_gb.html)
    return [cum_yield_gb]


def overlay_histogram(df, path, palette=None):
    """
    Use plotly to create an overlay of length histograms
    Return html code

    Only has 10 colors, which get recycled up to 5 times.
    """
    if palette is None:
        palette = plotly.colors.DEFAULT_PLOTLY_COLORS * 5

    overlay_hist = Plot(
        path=path + "NanoComp_OverlayHistogram.html",
        title="Histogram of read lengths")
    overlay_hist.html = plot_overlay_histogram(
        df, palette, title=overlay_hist.title, histnorm="")
    with open(overlay_hist.path, 'w') as html_out:
        html_out.write(overlay_hist.html)

    overlay_hist_normalized = Plot(
        path=path + "NanoComp_OverlayHistogram_Normalized.html",
        title="Normalized histogram of read lengths")
    overlay_hist_normalized.html = plot_overlay_histogram(
        df, palette, title=overlay_hist_normalized.title, histnorm="probability")
    with open(overlay_hist_normalized.path, 'w') as html_out:
        html_out.write(overlay_hist_normalized.html)

    overlay_log_hist = Plot(
        path=path + "NanoComp_OverlayLogHistogram.html",
        title="Histogram of log transformed read lengths")
    overlay_log_hist.html = plot_overlay_log_histogram(
        df, palette, title=overlay_log_hist.title, histnorm="")
    with open(overlay_log_hist.path, 'w') as html_out:
        html_out.write(overlay_log_hist.html)

    overlay_log_hist_normalized = Plot(
        path=path + "NanoComp_OverlayLogHistogram_Normalized.html",
        title="Normalized histogram of log transformed read lengths")
    overlay_log_hist_normalized.html = plot_overlay_log_histogram(
        df, palette, title=overlay_log_hist_normalized.title, histnorm="probability")
    with open(overlay_log_hist_normalized.path, 'w') as html_out:
        html_out.write(overlay_log_hist_normalized.html)

    return [overlay_hist, overlay_hist_normalized, overlay_log_hist, overlay_log_hist_normalized]


def plot_overlay_histogram(df, palette, title, histnorm):
    data = [go.Histogram(x=df.loc[df.dataset == d, "lengths"],
                         opacity=0.4,
                         name=d,
                         histnorm=histnorm,
                         marker=dict(color=c))
            for d, c in zip(df.dataset.unique(), palette)]

    return plotly.offline.plot(
        {"data": data,
         "layout": go.Layout(barmode='overlay',
                             title=title)},
        output_type="div",
        show_link=False)


def plot_overlay_log_histogram(df, palette, title, histnorm):
    data = [go.Histogram(x=np.log10(df.loc[df.dataset == d, "lengths"]),
                         opacity=0.4,
                         name=d,
                         histnorm=histnorm,
                         marker=dict(color=c))
            for d, c in zip(df.dataset.unique(), palette)]
    xtickvals = [10**i for i in range(10) if not 10**i > 10 * np.amax(df["lengths"])]
    return plotly.offline.plot(
        {"data": data,
         "layout": go.Layout(barmode='overlay',
                             title=title,
                             xaxis=dict(tickvals=np.log10(xtickvals),
                                        ticktext=xtickvals)
                             )
         },
        output_type="div",
        show_link=False)


def run_tests():
    import pickle
    df = pickle.load(open("nanotest/sequencing_summary.pickle", "rb"))
    scatter(
        x=df["lengths"],
        y=df["quals"],
        names=['Read lengths', 'Average read quality'],
        path="LengthvsQualityScatterPlot",
        plots={'dot': 1, 'kde': 1, 'hex': 1, 'pauvre': 1},
        plot_settings=dict(font_scale=1))
    time_plots(
        df=df,
        path=".",
        color="#4CB391",
        plot_settings=dict(font_scale=1))
    length_plots(
        array=df["lengths"],
        name="lengths",
        path=".")
    spatial_heatmap(
        array=df["channelIDs"],
        title="Number of reads generated per channel",
        path="ActivityMap_ReadsPerChannel")


checkvalidColor = check_valid_color
checkvalidFormat = check_valid_format
spatialHeatmap = spatial_heatmap
lengthPlots = length_plots
timePlots = time_plots


if __name__ == "__main__":
    run_tests()
