"""
Figure 4.2, Panel A — best SGNH hydrolase predictors per K-locus (precision ≥ 0.50).

Horizontal dot plot:
  - Y-axis : F1-score
  - X-axis : K-locus (all 35 shown; those without a predictor are greyed out)
  - Marker shape   : clustering level (PCI80C80 = *, PCI50C50 = D, …)
  - Marker size    : PC_abundance (3 size classes: 1–10, 11–25, >25)
  - Edge linewidth : precision ≥ 0.8 → thick (good); precision < 0.8 → thin (likely)
  - Colour         : sgnh_domain_color

Reads:
    sgnh_dir/gwas_sgnh_best.tsv  (12 loci: KL2, KL6, KL8, KL11, KL16, KL24, KL30, KL35, KL55, KL64, KL111, KL125)

Outputs:
    plots_dir/figure4_2-panelA.png / .pdf
    plots_dir/legends/figure4_2-panelA-legend.png / .pdf
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pandas as pd

_MARKERS = {
    "PCI80C80": "*",
    "PCI50C50": "D",
    "PCI80C50": "^",
    "PCI50C80": "s",
    "PCI00C80": "P",
    "PCI00C50": "o",
}

_SIZE_BINS   = [(1, 10), (11, 25), (26, float("inf"))]
_SIZE_VALUES = [60, 180, 450]
_SIZE_LABELS = ["1–10", "11–25", ">25"]

_BEST_PERFORMING_COLOR = "#d62728"


def _kl_sort_key(locus: str) -> int:
    try:
        return int(locus.replace("KL", ""))
    except ValueError:
        return 9999


def _abundance_to_size(ab: float) -> float:
    for (lo, hi), s in zip(_SIZE_BINS, _SIZE_VALUES):
        if lo <= ab <= hi:
            return s
    return _SIZE_VALUES[-1]


def plot_figure4_2_panelA(
    gwas_sgnh_best_tsv: Path,
    best_predictors_tsv: Path,
    pyseer_hits_tsv: Path,
    plots_dir: Path,
    style=None,
) -> None:
    """
    Produce Figure 4.2 Panel A: horizontal dot plot of best SGNH predictor per K-locus.

    Args:
        gwas_sgnh_best_tsv:  sgnh-hydrolases/gwas_sgnh_best.tsv
        best_predictors_tsv: rbp_best_predictors/best_predictors_gwas.tsv
        pyseer_hits_tsv:     input_dir/gwas/3_GWAS/3_PROCESSING/pyseer_hits_filtered.tsv; source of the 35-locus x-axis
        plots_dir:           output directory for plots
        style:               cfg.style (optional)
    """
    axis_label_fs = getattr(style, "axis_label_fontsize",   12)
    axis_label_fw = getattr(style, "axis_label_fontweight", "bold")
    tick_fs       = getattr(style, "tick_fontsize",          8)
    tick_fw       = getattr(style, "tick_fontweight",        "bold")
    dpi           = getattr(style, "dpi",                    300)
    fill_color    = getattr(style, "sgnh_domain_color",      "#c9a227")

    df = pd.read_csv(gwas_sgnh_best_tsv, sep="\t")

    bp = pd.read_csv(best_predictors_tsv, sep="\t")
    best_sgnh_loci = set(bp.loc[bp["ecod_folder"] == "sgnh-ecod-reported-topology", "locus"])

    pyseer = pd.read_csv(pyseer_hits_tsv, sep="\t")
    pyseer = pyseer[pyseer["mode"] == "lasso"]

    loci_present = set(df["locus"].unique())
    all_loci     = sorted(pyseer["locus"].unique(), key=_kl_sort_key)
    x_pos        = {locus: i for i, locus in enumerate(all_loci)}

    fig, ax = plt.subplots(figsize=(10, 4))

    for _, row in df.iterrows():
        marker    = _MARKERS.get(row["clustering_level"], "o")
        size      = _abundance_to_size(row["PC_abundance"]) * (2.5 if row["locus"] == "KL16" else 1.0)
        linewidth = 2.0 if row["precision"] >= 0.8 else 0.5
        color     = _BEST_PERFORMING_COLOR if row["locus"] in best_sgnh_loci else fill_color

        ax.scatter(
            x_pos[row["locus"]], row["F1_score"],
            marker=marker,
            s=size,
            color=color,
            edgecolors="black",
            linewidths=linewidth,
            zorder=3,
            alpha=0.9,
        )

    ax.set_xlim(-0.5, len(all_loci) - 0.5)
    ax.set_xticks(range(len(all_loci)))
    ax.set_xticklabels(all_loci, rotation=45, ha="right")
    ax.set_ylim(0, 1)
    ax.set_ylabel("F1-score", fontsize=axis_label_fs, fontweight=axis_label_fw)
    ax.set_xlabel("K-locus",  fontsize=axis_label_fs, fontweight=axis_label_fw)

    ax.tick_params(axis="x", labelsize=tick_fs)
    ax.tick_params(axis="y", labelsize=tick_fs)
    for lbl in ax.get_xticklabels():
        if lbl.get_text() in loci_present:
            lbl.set_fontweight(tick_fw)
        else:
            lbl.set_fontweight("light")
            lbl.set_color("gray")
    for lbl in ax.get_yticklabels():
        lbl.set_fontweight(tick_fw)

    ax.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.4, color="gray")
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plots_dir.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    for ext in ("png", "pdf"):
        out = plots_dir / f"figure4_2-panelA.{ext}"
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        print(f"  [figure4_2-panelA] → {out.name}")
    plt.close(fig)

    # -----------------------------------------------------------------------
    # Legend — three groups in one column, saved to legends/
    # -----------------------------------------------------------------------
    versions_present = sorted(df["clustering_level"].unique())

    shape_handles = [
        plt.scatter([], [], marker=_MARKERS.get(v, "o"), s=80,
                    color=fill_color, edgecolors="black", linewidths=0.8,
                    label=v)
        for v in versions_present
    ]
    size_handles = [
        plt.scatter([], [], marker="o", s=s, color=fill_color,
                    edgecolors="black", linewidths=0.8, label=label)
        for s, label in zip(_SIZE_VALUES, _SIZE_LABELS)
    ]
    edge_handles = [
        mlines.Line2D([], [], marker="o", color="w", markerfacecolor=fill_color,
                      markeredgecolor="black", markeredgewidth=2.0,
                      markersize=8, label="precision ≥ 0.8"),
        mlines.Line2D([], [], marker="o", color="w", markerfacecolor=fill_color,
                      markeredgecolor="black", markeredgewidth=0.5,
                      markersize=8, label="precision < 0.8"),
    ]
    color_handles = [
        plt.scatter([], [], marker="o", s=80, color=_BEST_PERFORMING_COLOR,
                    edgecolors="black", linewidths=0.8,
                    label="best-performing predictor"),
        plt.scatter([], [], marker="o", s=80, color=fill_color,
                    edgecolors="black", linewidths=0.8,
                    label="additional predictor"),
    ]

    def _header(text):
        latex = text.replace(" ", "~")
        return mlines.Line2D([], [], color="none", label=f"$\\bf{{{latex}}}$")

    def _spacer():
        return mlines.Line2D([], [], color="none", label=" ")

    all_handles = (
        [_header("Clustering level")]               + shape_handles +
        [_spacer(), _header("Number of sequences")] + size_handles  +
        [_spacer(), _header("Precision")]            + edge_handles  +
        [_spacer(), _header("Predictor class")]      + color_handles
    )
    fig_leg, ax_leg = plt.subplots(figsize=(2.8, 5.5))
    ax_leg.axis("off")
    ax_leg.legend(
        handles=all_handles, fontsize=tick_fs,
        loc="center", ncol=1,
        frameon=False, handletextpad=0.5, labelspacing=0.4,
    )

    legends_dir = plots_dir / "legends"
    legends_dir.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    for ext in ("png", "pdf"):
        out = legends_dir / f"figure4_2-panelA-legend.{ext}"
        fig_leg.savefig(out, dpi=dpi, bbox_inches="tight")
        print(f"  [figure4_2-panelA-legend] → legends/{out.name}")
    plt.close(fig_leg)
