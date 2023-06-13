import numpy as np
import seaborn as sb
from matplotlib import pyplot as plt


def plot_metrics(df, col_name='method', batch_metrics=None, bio_metrics=None, palette=None):

    sb.set_context("paper")

    if palette is None:
        palette = "viridis_r"
        # sb.color_palette("ch:s=.25,rot=-.25", as_cmap=True)

    if batch_metrics is None:
        batch_metrics = ["ASW_label/batch", "PCR_batch", "graph_conn", "kBET", "iLISI"]

    if bio_metrics is None:
        bio_metrics = [
            "NMI_cluster/label",
            "ARI_cluster/label",
            "ASW_label",
            "cell_cycle_conservation",
            "isolated_label_F1",
            "isolated_label_silhouette",
            "cLISI",
            "hvg_overlap",
            "trajectory",
        ]

    df = pd.melt(df, id_vars=col_name, var_name="metric", value_name="value")
    df['metric'] = df[col_name].copy()

    conditions = [(df['metric'].isin(batch_metrics)), (df['metric'].isin(bio_metrics))]
    metric_type = ["Batch Correction", "Biological Conservation"]

    df["metric_type"] = np.select(conditions, metric_type)
    df["metric"] = df["metric"].str.replace("_", " ")
    df["rank"] = df.groupby("metric")["value"].rank(ascending=False)

    dims = df[["metric_type", "metric"]].drop_duplicates()["metric_type"].value_counts()
    n_metrics = dims.sum()
    n_methods = df["method"].nunique()
    dim_x = (n_metrics + dims.shape[0]) * 0.48
    dim_y = np.max([2.5, n_methods])

    # Build plot
    fig, axs = plt.subplots(
        nrows=1,
        ncols=dims.shape[0],
        figsize=(dim_x, dim_y),
        sharey=True,
        gridspec_kw=dict(width_ratios=list(dims)),
    )

    for i, metric_type in enumerate(dims.index):
        sb.despine(bottom=True, left=True)
        sb.scatterplot(
            data=df.query(f'metric_type == "{metric_type}"'),
            x="metric",
            y="method",
            hue="rank", 
            palette=palette,
            size="value",
            sizes=(0, 100),
            ax=axs[i],
        )
        axs[i].set(title=metric_type, xlabel=None, ylabel=None)
        axs[i].tick_params(axis="x", rotation=90)
        axs[i].legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)

        for t in axs[i].legend_.texts:
            t.set_text(t.get_text()[:5])

    fig.tight_layout()