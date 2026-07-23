from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.ticker import PercentFormatter


def plot_jp_heatmap(
    similarity_csv: Path,
    output_path: Path,
    style=None,
) -> None:
    dpi = style.dpi if style else 200
    label_fs = style.axis_label_fontsize if style else 12
    label_fw = style.axis_label_fontweight if style else "bold"
    tick_fs = style.tick_fontsize if style else 8
    tick_fw = style.tick_fontweight if style else "bold"

    df = pd.read_csv(similarity_csv)
    df = df.dropna(subset=["path_jaccard_core", "path_jaccard_branch"])
    n_total = len(df)
    if n_total == 0:
        raise ValueError("No data after dropping NaN path_jaccard_core / path_jaccard_branch.")

    fig, ax = plt.subplots(figsize=(5, 4))

    hb = ax.hexbin(
        df["path_jaccard_core"].values,
        df["path_jaccard_branch"].values,
        gridsize=15,
        cmap="Reds",
        mincnt=1,
    )

    counts = hb.get_array()
    fractions = counts / float(n_total)
    hb.set_array(fractions)
    hb.set_norm(Normalize(vmin=0.0, vmax=fractions.max() if fractions.size else 1.0))

    cbar = fig.colorbar(hb, ax=ax, format=PercentFormatter(xmax=1))
    cbar.set_label("\nfraction of all comparisons", fontsize=label_fs, fontweight=label_fw)

    ax.set_xlabel("\npath Jaccard (core)", fontsize=label_fs, fontweight=label_fw)
    ax.set_ylabel("\npath Jaccard (branch)", fontsize=label_fs, fontweight=label_fw)

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(tick_fs)
        label.set_fontweight(tick_fw)

    fig.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
