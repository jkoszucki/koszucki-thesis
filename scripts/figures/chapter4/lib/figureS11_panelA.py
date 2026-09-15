"""
S11 Figure, Panel A — Functional annotation of putative acetyltransferases.

For each candidate acetyltransferase, reads its FoldSeek-vs-AFDB50
(AlphaFold3/UniProt v6) result table, keeps hits with min. probability >= 0.7, and
assigns each remaining hit to exactly one functional category by matching its TARGET
annotation (first matching category wins, in priority order). Two subplots apply two
different Tmscore cutoffs to the same hits.

Categories (keyword matched against TARGET, case-insensitive; priority order = list order):
    Acetyltransferase, Acyltransferase, OpgC protein, Hexapeptide / Trimeric LpxA, Other

"Hexapeptide" and "Trimeric LpxA" are one category: the hexapeptide repeat *is* the
left-handed beta-helix, and "Trimeric LpxA-like enzymes" is simply the ECOD topology name
for proteins built from it. Splitting them separated two names for the same structural
feature, so they are counted together under one colour.

Every filtered hit is counted: anything matching none of the four named patterns falls into
"Other" (grey) rather than being discarded. "Other" is ~7-10% of hits and is dominated by
named O-acetylases whose annotation carries no 'acetyltransferase' string -- chiefly
Peptidoglycan/LPS O-acetylase OafA/YrhL, Fucose 4-O-acetylase, NodL and ExoZ -- alongside
genuinely unannotated entries ("Transferase", "Uncharacterized protein"). Folding those
named O-acetylases into an acetyltransferase category would be a functional claim this
figure does not make, so they stay in "Other".

Note on OpgC: the pattern matches the literal string "OpgC" only. OpgC's UniProt name is
"Glucans biosynthesis protein C", so hits annotated that way land in "Other" and the OpgC
bar is correspondingly small -- a deliberate choice, not an oversight. Because that leaves
it at 0.0%/0.1% of each cluster's hits, too small to see, it is now listed in
`_HIDDEN_CATEGORIES`: still matched and printed, but folded into "Other" for the plot so it
takes no legend entry. The chapter text's OpgC prevalence (14 of 37 long-cluster proteins)
counts the UniProt name too and is computed outside this figure.

Each protein is one stacked bar: segment height = raw count of that protein's hits in each
category (not normalised), so taller bars mean more categorised hits. Both subplots share
the same y-axis, so the TM-score >= 0.50 subplot's systematically taller bars visually show
that the lower threshold recovers far more hits than TM-score >= 0.75. Categories keep the
same colour across both proteins and subplots.

Protein set (73 plotted): the 69 K-locus candidates plus four reference acetyltransferases
from S1_Table.xlsx -- PROTEIN01-03 (characterised: NeuO from E. coli, and the two
Klebsiella enzymes WcsU/K2 and orf13-A1142/K57) and PROTEIN05 (best-performing GWAS
predictor for KL111). PROTEIN04 (KL30) is present in S3_Data but excluded here; see
`_EXCLUDED_PROTEINS`. Everything else `foldseek_dir` contains is plotted, so the set grows
simply by filing another result table there.

Reference proteins are marked by x-tick label colour rather than a separate track:
red (`_REFERENCE_LABEL_COLOR`) for the enzymes with experimental evidence of activity,
black for everything predicted — the GWAS predictor and the K-locus candidates alike. The
point of putting them on the same axis is to let the reader see whether a characterised
enzyme's annotation profile resembles the candidates that cluster with it structurally.

Note on PROTEIN01 (NeuO): until 2026-07-27 both S1_Table and S4_Data held a 170-residue
sequence/model for it, which was in fact KL1_11_wcsT — so this panel plotted KL1_11_wcsT's
FoldSeek table twice and the TM-align ordering collapsed the two as near-identical
(TM 0.950). The real 307-residue NeuO is now in place everywhere; it stays in structural
cluster 1 but no longer pairs with KL1_11_wcsT above 0.73.

X-axis order: by default the proteins are sorted by structural similarity, using the
TM-align dendrogram leaf order written by processing/acetyl-proc (see
`structure_order.py`). Ordering by structure -- an axis independent of the FoldSeek
annotations being plotted -- lets the reader judge whether annotation composition tracks
fold. A dashed black divider separates the structural clusters. Proteins with no AF3
model cannot be placed by fold; they are appended at the right behind a dotted grey
divider, which is deliberately a different line style because they are unplaced, not a
third cluster. Without the order file the proteins fall back to name order.

Reads:
    foldseek_dir/*.xlsx  (one per candidate protein; column TARGET holds the
                           functional annotation of each hit, PROB. and Tmscore the filters)
    structure_order_tsv  (optional) cps_acetylases/acetylases_kloci_structure_order.tsv

Outputs:
    plots_dir/supplementary/figureS11-panelA.png / .pdf
    plots_dir/supplementary/legends/figureS11-panelA-legend.png / .pdf
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

PROB_MIN = 0.7
TM_SCORE_HIGH = 0.75
TM_SCORE_LOW  = 0.5

# Horizontal slack per bar beyond the bare 1.4x-font-size needed by 45-degree labels.
# Raise it if labels ever collide; lower it to tighten the panel.
_LABEL_PACKING = 1.15

# Combined height of the two stacked subplots, in inches. Width is derived from the label
# size and bar count (see below), so this is the only free dimension.
_PANEL_HEIGHT = 4.94

_OTHER_LABEL = "Other"

# Proteins whose FoldSeek result file exists in S3_Data but is excluded from this panel.
# S3_Data is deposited supplementary data and keeps the file; the exclusion is a figure-level
# decision, recorded here rather than by deleting evidence.
_EXCLUDED_PROTEINS = {
    "PROTEIN04_GWAS_AC_K30":
        "false-positive acetyltransferase call: FoldSeek returns 234 hits of which exactly "
        "one passes PROB. >= 0.7, and that hit is the protein itself (AF-A0A483L4A7, 100% "
        "identity); the remainder are <5% identity coiled-coil matches to chemotaxis / "
        "methyl-accepting transducer proteins. All-vs-all TM-align agrees: its best score "
        "against the other 73 proteins is 0.347, so it shares a fold with none of them. "
        "Re-running the search as a monomer reproduced the trimer result exactly, ruling "
        "out a search artefact.",
}

_LBH_LABEL = "Hexapeptide / Trimeric LpxA"

# Categories still matched and reported on stdout, but folded into "Other" for the plot so
# they take no legend entry (user, 2026-07-28). OpgC is here because the pattern matches the
# literal string only (see the note above): it catches 5 of the 14 long-cluster proteins
# that have an OpgC hit under any name, and 0.0%/0.1% of each cluster's hits, so the legend
# entry promised a category the reader could not see. The count is still printed, and the
# figure loses nothing but an invisible sliver.
_HIDDEN_CATEGORIES = {"OpgC protein"}

_CATEGORY_PATTERNS = {
    "Acetyltransferase": re.compile(r"acetyl\s*transferase", re.IGNORECASE),
    "Acyltransferase":   re.compile(r"acyl\s*transferase", re.IGNORECASE),
    "OpgC protein":      re.compile(r"opgc", re.IGNORECASE),
    # One category: the hexapeptide repeat is the left-handed beta-helix, and "Trimeric
    # LpxA-like" is the ECOD topology name for that same fold.
    _LBH_LABEL:          re.compile(r"hexapeptide|trimeric\s*lpxa", re.IGNORECASE),
}

_CATEGORY_COLORS = {
    # Purple is the one hue figure 4.2 does not otherwise use: panel A spends green, gold,
    # blue, red and salmon; panel B's other categories are orange, aqua and pink; panel C
    # is a blue sequential ramp; panel D is red/blue. The category was previously
    # #5090dd, which collided with both panel A's sslbh blue and panel C's ramp.
    "Acetyltransferase": "#8c6bb1",  # purple
    "Acyltransferase":   "#ee8358",  # pastel orange
    "OpgC protein":      "#44bd91",  # pastel aqua
    # Deeper pink than the other pastels: this category is only ~2% of hits, so a pale
    # fill was hard to pick out in a thin band. Kept magenta-leaning to stay clearly
    # distinct from the red used on the x-tick labels of the characterised enzymes.
    _LBH_LABEL:          "#d9569b",  # pink — the merged left-handed beta-helix category
    _OTHER_LABEL:        "#bfbfbf",  # grey — matched no named pattern
}

# Reference proteins from S1_Table sit on the same axis as the 69 K-locus candidates; their
# x-tick labels are coloured red so the reader can pick them out without a separate track.
# The axis carries one distinction only: red = active (experimental evidence of activity),
# black = predicted. The GWAS predictor previously took a pale red of its own, which read as
# a third category when it is a prediction exactly like the 69 K-locus candidates it sits
# among (user, 2026-07-28); it is now black like them.
# Red carries no other meaning in these figures (green = SSRBH, gold = SGNH, blue = SSLBH)
# and collides with none of the bar categories above.
_REFERENCE_LABEL_COLOR = "#b2182b"   # active — experimentally characterised
_DEFAULT_LABEL_COLOR   = "#000000"   # predicted — GWAS predictor + K-locus candidates

_EXPERIMENTAL_PATTERN = re.compile(r"_MOD_AC_", re.IGNORECASE)


def _label_color(protein: str) -> str:
    if _EXPERIMENTAL_PATTERN.search(protein):
        return _REFERENCE_LABEL_COLOR
    return _DEFAULT_LABEL_COLOR


# The S1_Table reference proteins are plotted next to K-locus candidates named by locus and
# gene, so `PROTEIN02_MOD_AC_K2` says nothing to a reader of the figure. On the x-axis they
# are relabelled by what they are and which serotype they act on; the internal protein ids
# are untouched everywhere else, so this is display only.
_DISPLAY_LABEL_PATTERNS = (
    (re.compile(r"^PROTEIN\d+_MOD_AC_(?P<serotype>.+)$",  re.IGNORECASE), "active_AT_{serotype}"),
    (re.compile(r"^PROTEIN\d+_GWAS_AC_(?P<serotype>.+)$", re.IGNORECASE), "gwas_AT_{serotype}"),
)


def _display_label(protein: str) -> str:
    """X-axis label for a protein; K-locus candidates keep their id unchanged."""
    for pattern, template in _DISPLAY_LABEL_PATTERNS:
        m = pattern.match(protein)
        if m:
            return template.format(serotype=m.group("serotype"))
    return protein


def _all_categories() -> list[str]:
    """Categories drawn and shown in the legend, in stacking order."""
    return [c for c in _CATEGORY_PATTERNS if c not in _HIDDEN_CATEGORIES] + [_OTHER_LABEL]


def _load_hits(foldseek_dir: Path) -> dict[str, pd.DataFrame]:
    hits = {}
    for path in sorted(foldseek_dir.glob("*.xlsx")):
        if path.name.startswith("~$"):
            continue
        if path.stem in _EXCLUDED_PROTEINS:
            print(f"  [excluded] {path.stem}: {_EXCLUDED_PROTEINS[path.stem]}")
            continue
        hits[path.stem] = pd.read_excel(path)
    return hits


def _assign_category(target: str) -> str:
    for cat, pattern in _CATEGORY_PATTERNS.items():
        if pattern.search(target):
            return cat
    return _OTHER_LABEL


def _protein_category_counts(hits: dict[str, pd.DataFrame], tm_score_min: float) -> pd.DataFrame:
    """Per-protein hit counts per plotted category; hidden categories fold into "Other"."""
    all_categories = list(_CATEGORY_PATTERNS) + [_OTHER_LABEL]
    rows = {}
    for protein, df in hits.items():
        filtered = df[(df["PROB."] >= PROB_MIN) & (df["Tmscore"] >= tm_score_min)]
        value_counts = filtered["TARGET"].apply(_assign_category).value_counts()
        rows[protein] = pd.Series({cat: value_counts.get(cat, 0) for cat in all_categories})
    counts = pd.DataFrame(rows).T

    for hidden in _HIDDEN_CATEGORIES:
        n = counts[hidden].sum()
        if n:
            print(f"  [{hidden}] {n} hits ({(counts[hidden] > 0).sum()} proteins) "
                  f"folded into '{_OTHER_LABEL}' — no legend entry")
        counts[_OTHER_LABEL] += counts[hidden]
    return counts.drop(columns=list(_HIDDEN_CATEGORIES))


def _resolve_order(hits: dict[str, pd.DataFrame], structure_order_tsv: Path | None):
    """Return (ordered protein list, cluster boundaries, start of the unplaced block).

    Proteins with no AF3 model are absent from the structural order and cannot be placed
    by fold, so they are appended after every ordered protein. They are separated by a
    distinct divider rather than a cluster boundary — they are not a third cluster, they
    are simply unplaced, and the difference matters for reading the panel.
    """
    if structure_order_tsv is None or not Path(structure_order_tsv).is_file():
        if structure_order_tsv is not None:
            print(f"  [warn] {Path(structure_order_tsv).name} not found — falling back to name order")
        return sorted(hits), [], None

    order_df = pd.read_csv(structure_order_tsv, sep="\t").sort_values("order_index")
    ordered = [p for p in order_df["protein_id"] if p in hits]

    clusters = order_df.set_index("protein_id")["structure_cluster"]
    boundaries = [
        i - 0.5
        for i in range(1, len(ordered))
        if clusters.get(ordered[i]) != clusters.get(ordered[i - 1])
    ]

    unplaced_start = None
    dropped = sorted(set(hits) - set(ordered))
    if dropped:
        print(f"  [warn] no structural placement (no AF3 model), appended at the right: "
              f"{', '.join(dropped)}")
        unplaced_start = len(ordered) - 0.5
        ordered += dropped

    print(f"  ordered by structural similarity ({clusters.nunique()} clusters)")
    return ordered, boundaries, unplaced_start


def plot_figureS11_panelA(
    foldseek_dir: Path,
    plots_dir: Path,
    structure_order_tsv: Path | None = None,
    style=None,
) -> None:
    """
    Produce S11 Figure Panel A: functional annotation of putative acetyltransferases.

    Args:
        foldseek_dir:        S3_Data (69 FoldSeek-vs-AFDB50 .xlsx files, one per candidate protein)
        plots_dir:           output directory for plots
        structure_order_tsv: cps_acetylases/acetylases_kloci_structure_order.tsv; orders the
                             x-axis by TM-align structural similarity. None → name order.
        style:               cfg.style (optional)
    """
    axis_label_fs = getattr(style, "axis_label_fontsize",   12)
    axis_label_fw = getattr(style, "axis_label_fontweight", "bold")
    tick_fs       = getattr(style, "tick_fontsize",          8)
    tick_fw       = getattr(style, "tick_fontweight",        "bold")
    dpi           = getattr(style, "dpi",                    300)

    hits = _load_hits(foldseek_dir)
    print(f"  {len(hits)} candidate acetyltransferases loaded from {foldseek_dir.name}")

    proteins, cluster_boundaries, unplaced_start = _resolve_order(hits, structure_order_tsv)
    n_proteins = len(proteins)

    counts_high = _protein_category_counts(hits, TM_SCORE_HIGH).reindex(proteins)
    counts_low  = _protein_category_counts(hits, TM_SCORE_LOW).reindex(proteins)

    categories = _all_categories()
    x_pos = range(n_proteins)
    y_max = max(counts_high.sum(axis=1).max(), counts_low.sum(axis=1).max())

    for label, counts in (("TM≥0.75", counts_high), ("TM≥0.50", counts_low)):
        totals = counts.sum()
        n_all = totals.sum()
        other_pct = totals[_OTHER_LABEL] / n_all
        print(f"  {label}: {n_all} hits, {totals[_OTHER_LABEL]} in '{_OTHER_LABEL}' ({other_pct:.1%})")

    # 45-degree tick labels need ~1.4x their font size of horizontal room per bar, so at this
    # many bars the labels set the figure width, not the plot content. Width is therefore
    # derived from the label size rather than fixed: raising the font without widening the
    # figure would just collide the labels.
    tick_label_fs = axis_label_fs - 2.5        # x-tick labels (protein names)
    fig_width = n_proteins * tick_label_fs * 1.4 / 72 * _LABEL_PACKING
    fig, (ax_top, ax_bottom) = plt.subplots(2, 1, figsize=(fig_width, _PANEL_HEIGHT),
                                            sharex=True, sharey=True)

    for ax, counts, tm_score in (
        (ax_top,    counts_high, TM_SCORE_HIGH),
        (ax_bottom, counts_low,  TM_SCORE_LOW),
    ):
        bottom = pd.Series(0.0, index=proteins)
        for cat in categories:
            ax.bar(
                x_pos, counts[cat].values, bottom=bottom.values,
                color=_CATEGORY_COLORS[cat], edgecolor="black", linewidth=0.3,
                width=0.8,
            )
            bottom += counts[cat]

        ax.set_ylim(0, y_max * 1.05)
        ax.set_yticks([0, 250, 500, 750, 1000])
        ax.text(
            1.01, 0.5, f"≥ {tm_score:.2f}",
            transform=ax.transAxes, va="center", ha="left",
            fontsize=axis_label_fs, fontweight=tick_fw,
        )
        ax.tick_params(axis="y", labelsize=tick_fs)
        for lbl in ax.get_yticklabels():
            lbl.set_fontweight(tick_fw)
        for boundary in cluster_boundaries:
            ax.axvline(boundary, color="black", linestyle="--", linewidth=1.0, zorder=4)
        if unplaced_start is not None:
            ax.axvline(unplaced_start, color="#808080", linestyle=":", linewidth=1.0, zorder=4)

        ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.4, color="gray")
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # "TM-score" heads the column of per-subplot threshold values on the right
    ax_top.text(
        1.01, 1.04, "TM-score",
        transform=ax_top.transAxes, va="bottom", ha="left",
        fontsize=axis_label_fs, fontweight=axis_label_fw,
    )

    ax_bottom.set_xlim(-0.5, n_proteins - 0.5)
    ax_bottom.set_xticks(list(x_pos))
    ax_bottom.set_xticklabels([_display_label(p) for p in proteins], rotation=45, ha="right",
                              fontsize=tick_label_fs, fontweight="bold")
    for tick_label, protein in zip(ax_bottom.get_xticklabels(), proteins):
        tick_label.set_color(_label_color(protein))

    supylabel = fig.supylabel("Number of hits", fontsize=axis_label_fs + 3,
                              fontweight=axis_label_fw)

    out_dir = plots_dir / "supplementary"
    out_dir.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()

    # tight_layout centres supylabel on the figure, which the tall x-tick margin drags
    # below the axes; re-centre it on the two plots instead.
    supylabel.set_y((ax_top.get_position().y1 + ax_bottom.get_position().y0) / 2)
    for ext in ("png", "pdf"):
        out = out_dir / f"figureS11-panelA.{ext}"
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        print(f"  [figureS11-panelA] → {out.name}")
    plt.close(fig)

    # -----------------------------------------------------------------------
    # Legend — category colours, shared across both subplots
    # -----------------------------------------------------------------------
    handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=_CATEGORY_COLORS[cat], edgecolor="black", linewidth=0.3, label=cat)
        for cat in categories
    ]
    # Label-colour key, each entry shown only when such a protein is on the axis. The label
    # text is the distinction itself — the id prefixes it used to name in brackets are
    # visible on the axis already.
    label_key = [
        (_REFERENCE_LABEL_COLOR, "active"),
        (_DEFAULT_LABEL_COLOR,   "predicted"),
    ]
    shown_key = [
        (color, text) for color, text in label_key
        if any(_label_color(p) == color for p in proteins)
    ]
    handles += [
        plt.Line2D([], [], linestyle="none", marker="$\\mathbf{Aa}$", markersize=9,
                   color=color, label=text)
        for color, text in shown_key
    ]

    fig_leg, ax_leg = plt.subplots(figsize=(3.2, 1.8 + 0.25 * len(shown_key)))
    ax_leg.axis("off")
    ax_leg.legend(
        handles=handles, fontsize=tick_fs,
        loc="center", ncol=1,
        frameon=False, handletextpad=0.5, labelspacing=0.4,
    )

    legends_dir = plots_dir / "supplementary" / "legends"
    legends_dir.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    for ext in ("png", "pdf"):
        out = legends_dir / f"figureS11-panelA-legend.{ext}"
        fig_leg.savefig(out, dpi=dpi, bbox_inches="tight")
        print(f"  [figureS11-panelA-legend] → supplementary/legends/{out.name}")
    plt.close(fig_leg)
