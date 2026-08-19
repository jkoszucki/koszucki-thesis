from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text


def plot_figure3_2_panelA(
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

    x = df["path_jaccard_core"].values
    y = df["path_jaccard_branch"].values
    labels = (df["structure_id_1"] + "/" + df["structure_id_2"]).values

    fig, ax = plt.subplots(figsize=(5, 4))

    ax.scatter(x, y, s=18, alpha=0.55, color="#888888", linewidths=0)

    # Threshold lines: sim_core (x=0.6), sim_branch (y=0.6)
    ax.axvline(0.6, linestyle="--", linewidth=0.8, color="#bbbbbb")
    ax.axhline(0.6, linestyle="--", linewidth=0.8, color="#bbbbbb")

    # Pearson r
    if len(x) >= 2 and np.std(x) > 0 and np.std(y) > 0:
        r_val = np.corrcoef(x, y)[0, 1]
        r_label = f"r = {r_val:.2f}"
    else:
        r_label = "r = NA"

    ax.set_title(r_label)

    # Label points where pj_branch >= 0.6 or pj_core >= 0.6
    mask = (y >= 0.6) | (x >= 0.6)

    ax.set_xlim(-0.02, 1.05)
    ax.set_ylim(-0.02, 1.05)
    ax.set_xlabel("\npath Jaccard (core)", fontsize=label_fs, fontweight=label_fw)
    ax.set_ylabel("\npath Jaccard (branch)", fontsize=label_fs, fontweight=label_fw)

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(tick_fs)
        label.set_fontweight(tick_fw)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()

    # texts = []
    # for xi, yi, lab, is_out in zip(x, y, labels, mask):
    #     if is_out:
    #         texts.append(ax.text(xi, yi, str(lab), fontsize=3.5, fontweight="bold"))
    # if texts:
    #     adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="#aaaaaa", lw=0.6))

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
