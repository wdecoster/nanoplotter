import sys
import logging
from nanoplotter.plot import Plot
from datetime import timedelta
import seaborn as sns
import matplotlib.pyplot as plt
from math import ceil
import pandas as pd
import numpy as np


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
