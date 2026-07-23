"""
Chapter 2 — SSRBH depolymerase predictors: precision vs recall (Figure 2.3, Panel B).

Scatter plot:
  - X-axis : precision
  - Y-axis : recall
  - Marker shape   : prediction_strength → strong = square (s), likely = triangle (^)
  - Marker size    : large, fixed
  - Alpha          : 0.65 (semi-transparent to show overlap)
  - Colour         : #0E470E (dark green)

Uses the same data join as figure2_3_panelA:
    depolymerases_gwas.tsv × pyseer_hits_filtered.tsv on (PC, locus).

Reads:
    depolymerases_gwas_tsv — output_dir/rbp_depolymerases/depolymerases_gwas.tsv
    pyseer_hits_tsv        — input_dir/gwas/3_GWAS/3_PROCESSING/pyseer_hits_filtered.tsv
    best_predictors_tsv    — output_dir/rbp_best_predictors/best_predictors_gwas.tsv

Outputs:
    plots_dir/figure2_3-panelB.png / .pdf
    plots_dir/legends/figure2_3-panelB-legend.png / .pdf
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

_STRENGTH_COLORS = {
    "strong": "#0E470E",
    "likely": "#d7ead9",
}
_BEST_PERFORMING_COLOR  = "#d62728"
_SIZE_BINS   = [(1, 10), (11, 25), (26, float("inf"))]
_SIZE_VALUES = [60, 180, 450]
_SIZE_LABELS = ["1–10", "11–25", ">25"]
_POINT_ALPHA = 0.65


def _abundance_to_size(ab: float) -> float:
    for (lo, hi), s in zip(_SIZE_BINS, _SIZE_VALUES):
        if lo <= ab <= hi:
            return s
    return _SIZE_VALUES[-1]


def _load_data(depolymerases_gwas_tsv: Path, pyseer_hits_tsv: Path, best_predictors_tsv: Path) -> pd.DataFrame:
    s3 = pd.read_csv(depolymerases_gwas_tsv, sep="\t")[["locus", "PC", "prediction_strength"]]
    pyseer = pd.read_csv(pyseer_hits_tsv, sep="\t")
    pyseer = pyseer[pyseer["mode"] == "lasso"]
    best_predictors = pd.read_csv(best_predictors_tsv, sep="\t")
    best_ssrbh_loci = set(best_predictors.loc[best_predictors["ecod_type"] == "ssrbh-ecod", "locus"])
    df = s3.merge(
        pyseer[["PC", "locus", "precision", "recall", "F1_score", "PC_abundance"]],
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
    return df


def plot_figure2_3_panelB(
    depolymerases_gwas_tsv: Path,
    pyseer_hits_tsv: Path,
    best_predictors_tsv: Path,
    plots_dir: Path,
    style=None,
) -> None:
    """
    Produce Figure 2.3B: precision vs recall scatter for SSRBH depolymerase predictors.

    Args:
        depolymerases_gwas_tsv: output_dir/rbp_depolymerases/depolymerases_gwas.tsv
        pyseer_hits_tsv:        input_dir/gwas/3_GWAS/3_PROCESSING/pyseer_hits_filtered.tsv
        best_predictors_tsv:    output_dir/rbp_best_predictors/best_predictors_gwas.tsv
        plots_dir:              output directory for plots
        style:                  cfg.style (optional)
    """
    axis_label_fs = getattr(style, "axis_label_fontsize",   12)
    axis_label_fw = getattr(style, "axis_label_fontweight", "bold")
    tick_fs       = getattr(style, "tick_fontsize",          8)
    tick_fw       = getattr(style, "tick_fontweight",        "bold")
    dpi           = getattr(style, "dpi",                    300)

    df = _load_data(depolymerases_gwas_tsv, pyseer_hits_tsv, best_predictors_tsv)

    if df.empty:
        print("  [figure2_3-panelB] WARNING: no rows after join — check PC/locus keys")
        return

    print(f"  [figure2_3-panelB] {len(df)} points from {df['locus'].nunique()} K-loci")

    fig, ax = plt.subplots(figsize=(5, 5))

    for _, row in df.iterrows():
        marker = _STRENGTH_MARKERS.get(row["prediction_strength"], "o")
        color  = _BEST_PERFORMING_COLOR if row["best_performing"] else _STRENGTH_COLORS.get(row["prediction_strength"], "#0E470E")
        ax.scatter(
            row["precision"], row["recall"],
            marker=marker,
            s=_abundance_to_size(row["PC_abundance"]),
            color=color,
            edgecolors="black",
            linewidths=0.6,
            alpha=_POINT_ALPHA,
            zorder=3,
        )

    ax.set_xlim(0, 1.1)
    ax.set_ylim(0, 1.1)
    ax.set_xlabel("Precision", fontsize=axis_label_fs, fontweight=axis_label_fw)
    ax.set_ylabel("Recall",    fontsize=axis_label_fs, fontweight=axis_label_fw)

    ax.tick_params(axis="both", labelsize=tick_fs)
    for lbl in ax.get_xticklabels() + ax.get_yticklabels():
        lbl.set_fontweight(tick_fw)

    ax.grid(linestyle="--", linewidth=0.5, alpha=0.4, color="gray")
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plots_dir.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    for ext in ("png", "pdf"):
        out = plots_dir / f"figure2_3-panelB.{ext}"
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        print(f"  [figure2_3-panelB] → {out.name}")
    plt.close(fig)

    # -----------------------------------------------------------------------
    # Legend
    # -----------------------------------------------------------------------
    shape_handles = [
        plt.scatter([], [], marker=_STRENGTH_MARKERS["strong"], s=120,
                    color=_BEST_PERFORMING_COLOR, edgecolors="black", linewidths=0.6,
                    alpha=_POINT_ALPHA, label="strong (best-performing predictor)"),
        plt.scatter([], [], marker=_STRENGTH_MARKERS["strong"], s=120,
                    color=_STRENGTH_COLORS["strong"], edgecolors="black", linewidths=0.6,
                    alpha=_POINT_ALPHA, label="strong"),
        plt.scatter([], [], marker=_STRENGTH_MARKERS["likely"], s=120,
                    color=_STRENGTH_COLORS["likely"], edgecolors="black", linewidths=0.6,
                    alpha=_POINT_ALPHA, label="likely"),
    ]

    def _header(text):
        latex = text.replace(" ", "~")
        return mlines.Line2D([], [], color="none", label=f"$\\bf{{{latex}}}$")

    size_handles = [
        plt.scatter([], [], marker="o", s=s, color=_STRENGTH_COLORS["strong"],
                    edgecolors="black", linewidths=0.6, alpha=_POINT_ALPHA, label=label)
        for s, label in zip(_SIZE_VALUES, _SIZE_LABELS)
    ]

    def _spacer():
        return mlines.Line2D([], [], color="none", label=" ")

    all_handles = (
        [_header("Prediction strength")]            + shape_handles +
        [_spacer(), _header("Number of sequences")] + size_handles
    )

    fig_leg, ax_leg = plt.subplots(figsize=(2.8, 3.8))
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
        out = legends_dir / f"figure2_3-panelB-legend.{ext}"
        fig_leg.savefig(out, dpi=dpi, bbox_inches="tight")
        print(f"  [figure2_3-panelB-legend] → legends/{out.name}")
    plt.close(fig_leg)
