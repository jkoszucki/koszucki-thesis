"""
Figure 4.2, Panel C — TM-score matrix of the structural representatives.

An alternative to the representative network (`figure4_2_panelC_representatives.py`) that
shows the same comparison without a threshold. The network binarises a continuum: it draws
an edge above a cut and nothing below, and the cut has to be defended. The matrix shows
every pairwise TM-score on a continuous scale, so the block structure — two folds with
almost no similarity between them — is read directly off the picture rather than asserted
by an edge rule.

Rows and columns are the representatives of the TM >= 0.75 clusters (19 at present), ordered by the
dendrogram of those representatives so structurally similar ones sit adjacent and the
blocks fall on the diagonal. The dendrogram is drawn alongside; its 0.5 cut is the same
same-fold threshold that defines the two clusters ordering figure 4.2B's x-axis, so the two
panels can be read against each other.

Colour: a single hue, light to dark, over the full TM range in the data (0.2-1.0). One hue
is what a magnitude scale takes — a multi-hue or rainbow ramp would imply categories that
do not exist here, and a single hue is inherently safe for colour-vision deficiency. Cells
are not individually labelled; the colourbar carries the magnitude and 324 numbers would
bury the block structure that is the point of the figure.

Reference proteins keep the panel B label colours — red = active (experimental evidence),
black = predicted — so the same protein reads the same way across panels. Cluster size is
appended to each label, since one row can stand for up to 13 proteins.

Outputs:
    plots_dir/figure4_2-panelC.png / .pdf
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

from figure4_2_panelC_tmalign import _node_type

mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

_CMAP        = "Blues"     # single hue, light -> dark; magnitude, not category
_VMIN, _VMAX = 0.20, 1.00
_SAME_FOLD   = 0.50        # dendrogram reference line; the panel 3.2B cluster cut

_REFERENCE_LABEL_COLOR  = "#b2182b"   # characterised enzyme  — matches figure4_2_panelB
_DEFAULT_LABEL_COLOR    = "#000000"   # everything predicted: GWAS + K-locus candidates


def _label_color(protein_id: str) -> str:
    """Red marks experimental evidence; black marks a prediction.

    The GWAS predictor used to take a pale red of its own, which read as a third
    category when it is simply another prediction (user, 2026-07-28).
    """
    if _node_type(protein_id) == "EXPERIMENTAL":
        return _REFERENCE_LABEL_COLOR
    return _DEFAULT_LABEL_COLOR


def _matrix(tmscore_tsv: Path, ids: list[str]) -> np.ndarray:
    idx = {k: i for i, k in enumerate(ids)}
    sim = np.ones((len(ids), len(ids)))
    for row in pd.read_csv(tmscore_tsv, sep="\t").itertuples(index=False):
        if row.id1 in idx and row.id2 in idx:
            sim[idx[row.id1], idx[row.id2]] = sim[idx[row.id2], idx[row.id1]] = row.tm_max
    return sim


def plot_figure4_2_panelC_heatmap(
    tmscore_tsv: Path,
    representatives_node_tsv: Path,
    plots_dir: Path,
    style=None,
    show_dendrogram: bool = True,
) -> None:
    """
    Draw the representative-vs-representative TM-score matrix.

    Args:
        tmscore_tsv:              cps_acetylases/tmalign/acetylases_kloci_tmscore.tsv
        representatives_node_tsv: plots/figure4_2-panelC/representatives/node.tsv —
                                  defines which proteins are representatives and how many
                                  proteins each one stands for
        plots_dir:                scripts/figures/chapter4/plots
        style:                    cfg.style (dpi, fonts)
        show_dendrogram:          draw the linkage tree beside the matrix
    """
    tick_fs = getattr(style, "tick_fontsize", 8)
    label_fs = getattr(style, "axis_label_fontsize", 12)
    dpi = getattr(style, "dpi", 300)

    nodes = pd.read_csv(representatives_node_tsv, sep="\t")
    ids = nodes["name"].tolist()
    sizes = dict(zip(nodes["name"], nodes["cluster_size"]))

    sim = _matrix(tmscore_tsv, ids)
    dist = 1.0 - sim
    np.fill_diagonal(dist, 0.0)
    z = linkage(squareform(dist, checks=False), method="average")
    order = dendrogram(z, no_plot=True)["leaves"]

    ordered_ids = [ids[i] for i in order]
    m = sim[np.ix_(order, order)]
    labels = [f"{_short(p)}  (n={sizes[p]})" for p in ordered_ids]

    # Dendrogram gets a narrow column; the matrix stays square.
    if show_dendrogram:
        fig, (ax_dend, ax) = plt.subplots(
            1, 2, figsize=(9.4, 7.0), gridspec_kw={"width_ratios": [1, 4], "wspace": 0.02},
        )
        # One neutral colour for the tree: the matrix already shows the block structure,
        # and scipy's default per-cluster colouring would add a second, competing
        # categorical scheme to a figure whose only colour encoding is magnitude.
        dendrogram(z, ax=ax_dend, orientation="left", no_labels=True,
                   link_color_func=lambda _: "#4d4d4d")
        ax_dend.axvline(1.0 - _SAME_FOLD, color="black", linestyle="--", linewidth=0.8)
        # scipy places leaf i at y = 10i + 5, so the tree spans [0, 10n]. Pinning the axis
        # to exactly that range — and letting the matrix fill its own axis (aspect="auto")
        # — makes leaf i sit on matrix row i. Without this the tree is drawn over the full
        # figure height while the square matrix is not, and the two silently disagree.
        ax_dend.set_ylim(10 * len(ordered_ids), 0)
        ax_dend.set_xlabel("1 − TM-score", fontsize=tick_fs)
        ax_dend.tick_params(axis="x", labelsize=tick_fs - 1)
        ax_dend.set_yticks([])
        for spine in ("top", "right", "left"):
            ax_dend.spines[spine].set_visible(False)
    else:
        fig, ax = plt.subplots(figsize=(7.6, 7.0))

    im = ax.imshow(m, cmap=_CMAP, vmin=_VMIN, vmax=_VMAX, aspect="auto")

    ax.set_xticks(range(len(ordered_ids)))
    ax.set_yticks(range(len(ordered_ids)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=tick_fs, fontweight="bold")
    ax.set_yticklabels(labels, fontsize=tick_fs, fontweight="bold")
    # Row labels go on the right: the dendrogram occupies the left margin, and labels
    # there would be drawn over the tree.
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    for tick, pid in zip(ax.get_xticklabels(), ordered_ids):
        tick.set_color(_label_color(pid))
    for tick, pid in zip(ax.get_yticklabels(), ordered_ids):
        tick.set_color(_label_color(pid))

    # Recessive cell separators — a 2px surface gap so adjacent cells stay distinct.
    ax.set_xticks(np.arange(-0.5, len(ordered_ids), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(ordered_ids), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.2)
    ax.tick_params(which="minor", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Steal space from both axes, not just the matrix: taking it from `ax` alone shortens
    # the matrix while the dendrogram keeps full height, and the leaves stop lining up
    # with the rows.
    cbar_axes = [ax_dend, ax] if show_dendrogram else [ax]
    cbar = fig.colorbar(im, ax=cbar_axes, orientation="horizontal",
                        fraction=0.035, pad=0.22, shrink=0.45, anchor=(1.0, 1.0))
    cbar.set_label("TM-score", fontsize=label_fs - 2, fontweight="bold")
    cbar.ax.tick_params(labelsize=tick_fs)
    cbar.outline.set_visible(False)

    plots_dir = Path(plots_dir)
    plots_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        out = plots_dir / f"figure4_2-panelC.{ext}"
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        print(f"  [figure4_2-panelC] → {out.name}")
    plt.close(fig)

    off = m[~np.eye(len(m), dtype=bool)]
    print(f"  {len(ordered_ids)} representatives covering {int(sum(sizes.values()))} proteins")
    print(f"  off-diagonal TM: min {off.min():.2f}  median {np.median(off):.2f}  max {off.max():.2f}")


def _short(protein_id: str) -> str:
    """Row label: K-locus for candidates, active_AT_/gwas_AT_ for reference proteins."""
    from figure4_2_panelC_tmalign import _label
    return _label(protein_id, _node_type(protein_id))
