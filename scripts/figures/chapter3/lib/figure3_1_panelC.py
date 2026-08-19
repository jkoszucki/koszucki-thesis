from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap


def _freq_colormap():
    """light gray (low %) → medium gray (moderate) → red (high %)."""
    return LinearSegmentedColormap.from_list(
        "freq_map",
        [(0.0, "#dddddd"), (0.5, "#888888"), (1.0, "#cc0000")],
    )


def plot_figure3_1_panelC(
    modifications_csv: Path,
    output_path: Path,
    style=None,
) -> None:
    dpi = style.dpi if style else 200
    label_fs = style.axis_label_fontsize - 3 if style else 12
    label_fw = style.axis_label_fontweight if style else "bold"
    tick_fs = style.tick_fontsize if style else 8
    tick_fw = style.tick_fontweight if style else "bold"

    pyr_color = style.pyruvylation_color if style else "#4f81bd"
    ace_color = style.acetylation_color if style else "#f2c94c"
    gray_color = style.gray_color if style else "#bfbfbf"

    df = pd.read_csv(modifications_csv)

    counts = df.groupby("modification").size().reset_index()
    counts.columns = ["modification", "n_total"]
    counts = counts.sort_values("n_total", ascending=False)
    mod_order = counts["modification"].tolist()

    def get_bar_color(mod):
        if mod == "pyruvylation":
            return pyr_color
        if mod == "acetylation":
            return ace_color
        return gray_color

    cmap = _freq_colormap()
    n_bins = 10
    bin_edges = np.linspace(0, 100, n_bins + 1)       # [0, 10, 20, ..., 100]
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2  # [5, 15, ..., 95]

    fig, (ax_bar, ax_na, ax_grid) = plt.subplots(
        1, 3, figsize=(8, 3), sharey=True,
        gridspec_kw={"wspace": 0.08, "width_ratios": [4, 0.5, 4]},
    )

    # ---- Left panel: horizontal bar chart ----
    y_pos = list(range(len(mod_order)))
    ax_bar.barh(
        y_pos,
        counts["n_total"].tolist(),
        color=[get_bar_color(m) for m in mod_order],
        edgecolor="black",
        linewidth=0.5,
    )
    ax_bar.set_yticks(y_pos)
    ax_bar.set_yticklabels(mod_order, fontsize=tick_fs, fontweight=tick_fw)
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel("\ncount", fontsize=label_fs, fontweight=label_fw)
    for label in ax_bar.get_xticklabels():
        label.set_fontsize(tick_fs)
        label.set_fontweight(tick_fw)
    for i, (_, row) in enumerate(counts.iterrows()):
        ax_bar.text(
            row["n_total"] + 0.3, i, f"n={int(row['n_total'])}",
            va="center", fontsize=tick_fs,
        )
    ax_bar.set_xlim(0, counts["n_total"].max() + 5)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)

    # ---- Middle and right panels: compute counts (shared global max) ----
    cell_h = 0.65
    cell_w = 1.0

    bin_counts = {}
    unknown_counts = {}
    for mod in mod_order:
        mod_df = df[df["modification"] == mod]["mean_frequency_percent"]
        freq_vals = mod_df.dropna().values
        unknown_counts[mod] = int(mod_df.isna().sum())
        for j in range(n_bins):
            lo, hi = bin_edges[j], bin_edges[j + 1]
            mask = (freq_vals >= lo) & (freq_vals <= hi if j == n_bins - 1 else freq_vals < hi)
            bin_counts[(mod, j)] = int(mask.sum())

    all_counts = list(bin_counts.values()) + list(unknown_counts.values())
    global_max = max(all_counts) or 1

    def draw_cell(ax, x, y, count):
        color = cmap(count / global_max) if count > 0 else "white"
        ax.add_patch(mpatches.Rectangle(
            (x, y - cell_h / 2), cell_w, cell_h,
            facecolor=color, edgecolor="#aaaaaa", linewidth=0.5,
        ))
        if count > 0:
            text_color = "white" if count / global_max > 0.4 else "#333333"
            ax.text(
                x + cell_w / 2, y, str(count),
                ha="center", va="center", fontsize=tick_fs, color=text_color,
            )

    # ---- Middle panel: unknown frequency column ----
    for i, mod in enumerate(mod_order):
        draw_cell(ax_na, 0, i, unknown_counts[mod])

    ax_na.set_xlim(0, cell_w)
    ax_na.set_xticks([cell_w / 2])
    ax_na.set_xticklabels(["unknown\nfrequency"], fontsize=tick_fs)
    ax_na.set_xlabel("", fontsize=label_fs, fontweight=label_fw)
    ax_na.spines["top"].set_visible(False)
    ax_na.spines["right"].set_visible(False)
    ax_na.spines["left"].set_visible(False)
    ax_na.tick_params(left=False, labelrotation=90)

    # ---- Right panel: 10-bin frequency grid ----
    for i, mod in enumerate(mod_order):
        for j in range(n_bins):
            draw_cell(ax_grid, j * cell_w, i, bin_counts[(mod, j)])

    ax_grid.set_xlim(0, n_bins * cell_w)
    ax_grid.set_xticks([0, 5, 10])
    ax_grid.set_xticklabels(["0 %", "50 %", "100 %"], fontsize=tick_fs)
    ax_grid.set_xlabel("\nmean frequency per\nCPS repeating unit (%)", fontsize=label_fs, fontweight=label_fw)
    ax_grid.spines["top"].set_visible(False)
    ax_grid.spines["right"].set_visible(False)
    ax_grid.spines["left"].set_visible(False)
    ax_grid.tick_params(left=False)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
