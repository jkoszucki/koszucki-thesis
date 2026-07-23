"""
Chapter 2 — SSRBH depolymerase predictors per K-locus (Figure 2.3, Panel A).

Horizontal dot plot:
  - X-axis : K-locus labels (all 35), sorted numerically
  - Y-axis : F1 score
  - Marker shape   : prediction_strength → strong = square (s), likely = triangle (^)
  - Marker size    : PC_abundance (3 size classes: 1–10, 11–25, >25)
  - Edge linewidth : prediction_strength → strong=1.5, likely=0.5
  - Colour         : #0E470E (dark green)

Multiple points per K-locus are possible (>1 depolymerase per capsule type).

Reads:
    depolymerases_gwas_tsv — output_dir/rbp_depolymerases/depolymerases_gwas.tsv
    pyseer_hits_tsv        — input_dir/gwas/3_GWAS/3_PROCESSING/pyseer_hits_filtered.tsv
    best_predictors_tsv    — output_dir/rbp_best_predictors/best_predictors_gwas.tsv

Outputs:
    plots_dir/figure2_3-panelA.png / .pdf
    plots_dir/legends/figure2_3-panelA-legend.png / .pdf
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pandas as pd

_STRENGTH_MARKERS = {
    "strong": "s",
    "likely": "^",
}

_SIZE_BINS   = [(1, 10), (11, 25), (26, float("inf"))]
_SIZE_VALUES = [60, 180, 450]
_SIZE_LABELS = ["1–10", "11–25", ">25"]

_STRENGTH_COLORS = {
    "strong": "#0E470E",
    "likely": "#d7ead9",
}
_BEST_PERFORMING_COLOR  = "#d62728"


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


def plot_figure2_3_panelA(
    depolymerases_gwas_tsv: Path,
    pyseer_hits_tsv: Path,
    best_predictors_tsv: Path,
    plots_dir: Path,
    style=None,
) -> None:
    """
    Produce depolymerase predictor dot plot: F1 score per K-locus.

    Args:
        depolymerases_gwas_tsv: output_dir/rbp_depolymerases/depolymerases_gwas.tsv; columns locus, PC, prediction_strength
        pyseer_hits_tsv:        input_dir/gwas/3_GWAS/3_PROCESSING/pyseer_hits_filtered.tsv; all lasso GWAS hits with F1, recall, PC_abundance
        best_predictors_tsv:    output_dir/rbp_best_predictors/best_predictors_gwas.tsv; columns locus, ecod_type, ...
        plots_dir:              output directory for plots
        style:                  cfg.style (optional)
    """
    axis_label_fs = getattr(style, "axis_label_fontsize",   12)
    axis_label_fw = getattr(style, "axis_label_fontweight", "bold")
    tick_fs       = getattr(style, "tick_fontsize",          8)
    tick_fw       = getattr(style, "tick_fontweight",        "bold")
    dpi           = getattr(style, "dpi",                    300)

    # ---- Load & join ----
    s3 = pd.read_csv(depolymerases_gwas_tsv, sep="\t")[["locus", "PC", "prediction_strength"]]
    pyseer = pd.read_csv(pyseer_hits_tsv, sep="\t")
    pyseer = pyseer[pyseer["mode"] == "lasso"]
    best_predictors = pd.read_csv(best_predictors_tsv, sep="\t")
    best_ssrbh_loci = set(best_predictors.loc[best_predictors["ecod_type"] == "ssrbh-ecod", "locus"])

    df = s3.merge(
        pyseer[["PC", "locus", "F1_score", "recall", "PC_abundance"]],
        on=["PC", "locus"],
        how="inner",
    )
    df = df.sort_values("F1_score", ascending=False).drop_duplicates(subset=["locus", "PC"])

    # Flag the single best-performing PC per best-performing locus (highest F1)
    best_pcs = set()
    for locus in best_ssrbh_loci:
        sub = df[df["locus"] == locus]
        if not sub.empty:
            best_pcs.add((locus, sub.loc[sub["F1_score"].idxmax(), "PC"]))
    df["best_performing"] = df.apply(lambda r: (r["locus"], r["PC"]) in best_pcs, axis=1)

    if df.empty:
        print("  [figure2_3-panelA] WARNING: no rows after join — check PC/locus keys")
        return

    print(f"  [figure2_3-panelA] {len(df)} points from {df['locus'].nunique()} K-loci")

    loci_present = set(df["locus"].unique())
    all_loci     = sorted(pyseer["locus"].unique(), key=_kl_sort_key)
    x_pos        = {locus: i for i, locus in enumerate(all_loci)}

    fig, ax = plt.subplots(figsize=(10, 4))

    for _, row in df.iterrows():
        if row["locus"] not in x_pos:
            continue
        marker    = _STRENGTH_MARKERS.get(row["prediction_strength"], "o")
        size      = _abundance_to_size(row["PC_abundance"])
        linewidth = 1.5 if row["prediction_strength"] == "strong" else 0.5
        color     = _BEST_PERFORMING_COLOR if row["best_performing"] else _STRENGTH_COLORS.get(row["prediction_strength"], "#0E470E")

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
    ax.set_xticklabels(all_loci, rotation=45, ha="center")
    ax.set_ylim(0, 1)
    ax.set_ylabel("F1 score", fontsize=axis_label_fs, fontweight=axis_label_fw)
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
        out = plots_dir / f"figure2_3-panelA.{ext}"
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        print(f"  [figure2_3-panelA] → {out.name}")
    plt.close(fig)

    # -----------------------------------------------------------------------
    # Legend
    # -----------------------------------------------------------------------
    shape_handles = [
        plt.scatter([], [], marker=_STRENGTH_MARKERS["strong"], s=80,
                    color=_BEST_PERFORMING_COLOR, edgecolors="black", linewidths=1.5,
                    label="strong (best-performing predictor)"),
        plt.scatter([], [], marker=_STRENGTH_MARKERS["strong"], s=80,
                    color=_STRENGTH_COLORS["strong"], edgecolors="black", linewidths=1.5,
                    label="strong"),
        plt.scatter([], [], marker=_STRENGTH_MARKERS["likely"], s=80,
                    color=_STRENGTH_COLORS["likely"], edgecolors="black", linewidths=0.5,
                    label="likely"),
    ]
    size_handles = [
        plt.scatter([], [], marker="o", s=s, color=_STRENGTH_COLORS["strong"],
                    edgecolors="black", linewidths=0.8, label=label)
        for s, label in zip(_SIZE_VALUES, _SIZE_LABELS)
    ]

    def _header(text):
        latex = text.replace(" ", "~")
        return mlines.Line2D([], [], color="none", label=f"$\\bf{{{latex}}}$")

    def _spacer():
        return mlines.Line2D([], [], color="none", label=" ")

    all_handles = (
        [_header("Prediction strength")]            + shape_handles +
        [_spacer(), _header("Number of sequences")] + size_handles
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
        out = legends_dir / f"figure2_3-panelA-legend.{ext}"
        fig_leg.savefig(out, dpi=dpi, bbox_inches="tight")
        print(f"  [figure2_3-panelA-legend] → legends/{out.name}")
    plt.close(fig_leg)
