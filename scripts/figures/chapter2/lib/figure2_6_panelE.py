"""
Chapter 2 — Figure 2.6 Panel E: Alignment length vs protein length (scatter).

For each active depolymerase pair with bidirectional coverage >= 50 %, plots
alignment length (X) against protein length (Y).  Each alignment contributes
two points — one for the shorter protein (square) and one for the longer
protein (circle).  Marker fill indicates whether the alignment spans the
C-terminus: filled = spans the C-terminal; hollow = does not span the C-terminal.

Reads:
    blastp_all_hits_tsv  — analysis_dir/blastp_all_hits.tsv

Writes:
    plots_dir/figure2_6-panelE.png / .pdf
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import pandas as pd

_ACTIVE_SOURCES = {"virulent", "prophage"}

_COLOR_SHORTER = "#1f77b4"
_COLOR_LONGER  = "#d62728"


def plot_figure2_6_panelE(
    blastp_all_hits_tsv: Path,
    plots_dir: Path,
    style=None,
) -> None:
    axis_label_fs = getattr(style, "axis_label_fontsize",   12) * 1.2
    axis_label_fw = getattr(style, "axis_label_fontweight", "bold")
    tick_fs       = getattr(style, "tick_fontsize",           8) * 1.2
    tick_fw       = getattr(style, "tick_fontweight",        "bold")
    dpi           = getattr(style, "dpi",                    300)

    hits = pd.read_csv(blastp_all_hits_tsv, sep="\t")

    active = hits[
        hits["q_source"].isin(_ACTIVE_SOURCES) &
        hits["s_source"].isin(_ACTIVE_SOURCES) &
        (hits["bidir_cov_min"] >= 0.50)
    ].copy()

    # Deduplicate symmetric pairs — keep canonical direction (query < target)
    active["_pair_key"] = active.apply(
        lambda r: "_||_".join(sorted([str(r["query"]), str(r["target"])])), axis=1
    )
    active = active.drop_duplicates(subset="_pair_key").copy()

    active["aln_len"]    = active["qend"] - active["qstart"] + 1
    active["shorter_len"] = active[["qlen", "slen"]].min(axis=1)
    active["longer_len"]  = active[["qlen", "slen"]].max(axis=1)

    spans   = active[active["reaches_end"] == True]
    nospans = active[active["reaches_end"] == False]

    print(f"  [figure2_6-panelE] Unique pairs: {len(active)}")
    print(f"  [figure2_6-panelE] Spans C-terminal:        n={len(spans)}")
    print(f"  [figure2_6-panelE] Does not span C-terminal: n={len(nospans)}")

    ms  = 35   # marker size
    lw  = 1.2  # edge linewidth for hollow markers
    tfs = tick_fs * 0.8 * 1.5

    ax_min = 400
    upper  = max(active["longer_len"].max(), active["aln_len"].max()) * 1.02

    rows = [
        ("shorter protein", "s", _COLOR_SHORTER, "shorter_len"),
        ("longer protein",  "o", _COLOR_LONGER,  "longer_len"),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(7, 3.5),
                             sharex=True, sharey=True)

    legend_handles = [
        mlines.Line2D([], [], marker="o", color="w",
                      markerfacecolor="gray", markeredgecolor="gray",
                      markersize=6, label="spans the C-terminal"),
        mlines.Line2D([], [], marker="o", color="w",
                      markerfacecolor="none", markeredgecolor="gray",
                      markersize=6, label="does not span the C-terminal"),
    ]

    for ax, (row_label, marker, color, len_col) in zip(axes, rows):
        ax.scatter(spans["aln_len"], spans[len_col],
                   marker=marker, color=color, s=ms, alpha=0.75,
                   edgecolors="none", zorder=3)
        ax.scatter(nospans["aln_len"], nospans[len_col],
                   marker=marker, facecolors="none", edgecolors=color,
                   linewidths=lw, s=ms, alpha=0.85, zorder=3)

        ax.plot([ax_min, upper], [ax_min, upper], color="gray", linewidth=0.8,
                linestyle="--", alpha=0.5, zorder=2)

        ax.set_xlim(ax_min, upper)
        ax.set_ylim(ax_min, upper)

        ax.set_title(row_label, fontsize=tfs, fontweight=axis_label_fw,
                     loc="left", pad=3, color=color)
        ax.set_xlabel("Alignment length (residues)", fontsize=tfs,
                      fontweight=axis_label_fw)
        ax.tick_params(labelsize=tfs * 0.7)
        for lbl in ax.get_xticklabels() + ax.get_yticklabels():
            lbl.set_fontweight(tick_fw)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(linestyle="--", linewidth=0.5, alpha=0.4, color="gray")
        ax.set_axisbelow(True)

    axes[0].set_ylabel("Protein length\n(residues)", fontsize=tfs,
                       fontweight=axis_label_fw)
    axes[0].legend(handles=legend_handles, fontsize=(tfs - 1) * 0.7,
                   frameon=False, loc="upper left")

    plots_dir = Path(plots_dir)
    plots_dir.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    for ext in ("png", "pdf"):
        out = plots_dir / f"figure2_6-panelE.{ext}"
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        print(f"  [figure2_6-panelE] → {out.name}")
    plt.close(fig)
