from nanoplotter.plot import Plot
from nanoplotter.timeplots import check_valid_time_and_sort, add_time_bins
import logging
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import plotly
import plotly.graph_objs as go


def violin_or_box_plot(df, y, figformat, path, y_name,
                       title=None, violin=True, log=False, palette=None):
    """Create a violin or boxplot from the received DataFrame.

    The x-axis should be divided based on the 'dataset' column,
    the y-axis is specified in the arguments
    """
    violin_comp = Plot(path=path + "NanoComp_" + y.replace(' ', '_') + '.' + figformat,
                       title="Comparing {}".format(y))
    if y == "quals":
        violin_comp.title = "Comparing base call quality scores"
    if violin:
        logging.info("Nanoplotter: Creating violin plot for {}.".format(y))
        ax = sns.violinplot(x="dataset",
                            y=y,
                            data=df,
                            inner=None,
                            cut=0,
                            palette=palette,
                            linewidth=0)
    else:
        logging.info("Nanoplotter: Creating box plot for {}.".format(y))
        ax = sns.boxplot(x="dataset",
                         y=y,
                         data=df,
                         palette=palette)
    if log:
        ticks = [10**i for i in range(10) if not 10**i > 10 * (10**np.amax(df[y]))]
        ax.set(yticks=np.log10(ticks),
               yticklabels=ticks)
    ax.set(title=title or violin_comp.title,
           ylabel=y_name)
    plt.xticks(rotation=30, ha='center')
    violin_comp.fig = ax.get_figure()
    violin_comp.save(format=figformat)
    plt.close("all")
    return [violin_comp]


def output_barplot(df, figformat, path, title=None, palette=None):
    """Create barplots based on number of reads and total sum of nucleotides sequenced."""
    logging.info("Nanoplotter: Creating barplots for number of reads and total throughput.")
    read_count = Plot(path=path + "NanoComp_number_of_reads." + figformat,
                      title="Comparing number of reads")
    ax = sns.countplot(x="dataset",
                       data=df,
                       palette=palette)
    ax.set(ylabel='Number of reads',
           title=title or read_count.title)
    plt.xticks(rotation=30, ha='center')
    read_count.fig = ax.get_figure()
    read_count.save(format=figformat)
    plt.close("all")

    throughput_bases = Plot(path=path + "NanoComp_total_throughput." + figformat,
                            title="Comparing throughput in gigabases sequenced")
    throughput = df.groupby('dataset')['lengths'].sum()
    ax = sns.barplot(x=list(throughput.index),
                     y=throughput / 1e9,
                     palette=palette,
                     order=df["dataset"].unique())
    ax.set(ylabel='Total gigabase sequenced',
           title=title or throughput_bases.title)
    plt.xticks(rotation=30, ha='center')
    throughput_bases.fig = ax.get_figure()
    throughput_bases.save(format=figformat)
    plt.close("all")
    return read_count, throughput_bases


def compare_sequencing_speed(df, figformat, path, title=None, palette=None):
    logging.info("Nanoplotter: creating comparison of sequencing speed over time.")
    seq_speed = Plot(path=path + "NanoComp_sequencing_speed_over_time." + figformat,
                     title="Sequencing speed over time")
    dfs = check_valid_time_and_sort(df, "start_time")
    dfs['timebin'] = add_time_bins(dfs)
    ax = sns.violinplot(x=dfs["timebin"],
                        y=dfs["lengths"] / dfs["duration"],
                        hue=dfs["dataset"],
                        inner=None,
                        cut=0,
                        linewidth=0)
    ax.set(xlabel='Interval (hours)',
           ylabel="Sequencing speed (nucleotides/second)")
    plt.xticks(rotation=45, ha='center', fontsize=8)
    seq_speed.fig = ax.get_figure()
    seq_speed.save(format=figformat)
    plt.close("all")
    return [seq_speed]


def compare_cumulative_yields(df, path, palette=None, title=None):
    if palette is None:
        palette = plotly.colors.DEFAULT_PLOTLY_COLORS * 5
    dfs = check_valid_time_and_sort(df, "start_time").set_index("start_time")

    logging.info("Nanoplotter: Creating cumulative yield plots using {} reads.".format(len(dfs)))
    cum_yield_gb = Plot(path=path + "NanoComp_CumulativeYieldPlot_Gigabases.html",
                        title="Cumulative yield")
    data = []
    for d, c in zip(df["dataset"].unique(), palette):
        s = dfs.loc[dfs["dataset"] == d, "lengths"].cumsum().resample('10T').max() / 1e9
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

    cum_yield_gb.fig = go.Figure({
        "data": data,
        "layout": go.Layout(barmode='overlay',
                            title=title or cum_yield_gb.title,
                            xaxis=dict(title="Time (hours)"),
                            yaxis=dict(title="Yield (gigabase)"),
                            )})
    cum_yield_gb.save()
    return [cum_yield_gb]


def overlay_histogram(df, path, palette=None):
    """
    Use plotly to create an overlay of length histograms
    Return html code, but also save as png

    Only has 10 colors, which get recycled up to 5 times.
    """
    if palette is None:
        palette = plotly.colors.DEFAULT_PLOTLY_COLORS * 5

    hist = Plot(path=path + "NanoComp_OverlayHistogram.html",
                title="Histogram of read lengths")
    hist.html, hist.fig = plot_overlay_histogram(df, palette, title=hist.title)
    hist.save()

    hist_norm = Plot(path=path + "NanoComp_OverlayHistogram_Normalized.html",
                     title="Normalized histogram of read lengths")
    hist_norm.html, hist_norm.fig = plot_overlay_histogram(
        df, palette, title=hist_norm.title, histnorm="probability")
    hist_norm.save()

    log_hist = Plot(path=path + "NanoComp_OverlayLogHistogram.html",
                    title="Histogram of log transformed read lengths")
    log_hist.html, log_hist.fig = plot_log_histogram(df, palette, title=log_hist.title)
    log_hist.save()

    log_hist_norm = Plot(path=path + "NanoComp_OverlayLogHistogram_Normalized.html",
                         title="Normalized histogram of log transformed read lengths")
    log_hist_norm.html, log_hist_norm.fig = plot_log_histogram(
        df, palette, title=log_hist_norm.title, histnorm="probability")
    log_hist_norm.save()

    return [hist, hist_norm, log_hist, log_hist_norm]


def plot_overlay_histogram(df, palette, title, histnorm=""):
    data = [go.Histogram(x=df.loc[df["dataset"] == d, "lengths"],
                         opacity=0.4,
                         name=d,
                         histnorm=histnorm,
                         marker=dict(color=c))
            for d, c in zip(df["dataset"].unique(), palette)]
    html = plotly.offline.plot(
        {"data": data,
         "layout": go.Layout(barmode='overlay',
                             title=title)},
        output_type="div",
        show_link=False)
    fig = go.Figure(
        {"data": data,
         "layout": go.Layout(barmode='overlay',
                             title=title)})
    return html, fig


def plot_log_histogram(df, palette, title, histnorm=""):
    """
    Plot overlaying histograms with log transformation of length
    Return both html and fig for png
    """
    data = [go.Histogram(x=np.log10(df.loc[df["dataset"] == d, "lengths"]),
                         opacity=0.4,
                         name=d,
                         histnorm=histnorm,
                         marker=dict(color=c))
            for d, c in zip(df["dataset"].unique(), palette)]
    xtickvals = [10**i for i in range(10) if not 10**i > 10 * np.amax(df["lengths"])]
    html = plotly.offline.plot(
        {"data": data,
         "layout": go.Layout(barmode='overlay',
                             title=title,
                             xaxis=dict(tickvals=np.log10(xtickvals),
                                        ticktext=xtickvals))},
        output_type="div",
        show_link=False)
    fig = go.Figure(
        {"data": data,
         "layout": go.Layout(barmode='overlay',
                             title=title,
                             xaxis=dict(tickvals=np.log10(xtickvals),
                                        ticktext=xtickvals))})
    return html, fig
