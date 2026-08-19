from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.ticker import PercentFormatter


def plot_branch_core_heatmap(similarity_csv: Path, output_path: Path, style=None) -> None:
    dpi = style.dpi if style else 200
    label_fs = style.axis_label_fontsize if style else 12
    label_fw = style.axis_label_fontweight if style else "bold"
    tick_fs = style.tick_fontsize if style else 8
    tick_fw = style.tick_fontweight if style else "bold"

    sim_long_df = pd.read_csv(similarity_csv)

    sim_long_df = sim_long_df.dropna(subset=["path_jaccard_core", "path_jaccard_total"])
    n_total = len(sim_long_df)
    if n_total == 0:
        raise ValueError("No data points after dropping NaNs in 'path_jaccard_core' and 'path_jaccard_total'.")

    plt.figure(figsize=(5, 4))

    hb = plt.hexbin(
        sim_long_df["path_jaccard_core"].values,
        sim_long_df["path_jaccard_total"].values,
        gridsize=15,
        cmap="Reds",
        mincnt=1,
    )

    counts = hb.get_array()
    fractions = counts / float(n_total)

    hb.set_array(fractions)
    hb.set_norm(Normalize(vmin=0.0, vmax=fractions.max() if fractions.size else 1.0))

    cbar = plt.colorbar(hb, format=PercentFormatter(xmax=1))
    cbar.set_label("\nFraction of all comparisons", fontsize=label_fs, fontweight=label_fw)

    plt.xlabel("\npath Jaccard (core)", fontsize=label_fs, fontweight=label_fw)
    plt.ylabel("\npath Jaccard (core + branch)", fontsize=label_fs, fontweight=label_fw)
    plt.title("", fontsize=label_fs, fontweight=label_fw)

    plt.axvline(x=0.9, color="gray", linestyle="--", linewidth=1.5)

    ax = plt.gca()
    for label in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
        label.set_fontweight(tick_fw)
        label.set_fontsize(tick_fs)

    plt.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=dpi)
    plt.close()
