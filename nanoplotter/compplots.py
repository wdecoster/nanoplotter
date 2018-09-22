from nanoplotter.plot import Plot
from nanoplotter.timeplots import check_valid_time_and_sort
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
