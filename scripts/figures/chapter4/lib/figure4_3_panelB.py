"""
Figure 4.3, Panel B — acetyltransferase detection tier vs SGNH deacetylase predictor.

A three-row table asking whether the co-occurrence reported in Figure 4.3A survives when
the acetyltransferase call is made on progressively weaker evidence. Each row is one
detection tier, tested against the same SGNH deacetylase predictor split:

    BLASTp           source_blastp   — homology to a characterised acetyltransferase
    keyword search   source_keyword  — Kaptive product annotation
    all candidates   union           — any method, i.e. the Figure 4.3A row 2 call

Universe: the 35 GWAS K-loci of Figure 4.3A (`pyseer_hits_filtered.tsv`, `mode == lasso`),
split 12 with an SGNH deacetylase predictor vs 23 without. This is a gene-presence test on
both axes, so all 35 loci are used — the 29-locus solved-structure restriction applies only
to tests involving acetylation itself. See `docs/chapters/CHAPTER4.md`, "Association tests".

Statistics: two-sided Fisher's exact test; odds ratios are the **conditional-MLE** OR, the
convention the chapter's association-test table uses. The bottom row reproduces test 2 of
that table (9/12 vs 9/23, OR 4.46, P = 0.075).

Nothing here is hard-coded — every count is computed from the pipeline TSVs and checked
against `_EXPECTED` before drawing, so a change upstream fails loudly instead of silently
redrawing a table the text no longer matches.

Outputs:
    plots_dir/figure4_3-panelB.png / .pdf
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio

mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

# Detection tiers, top to bottom: weakest evidence first, the union last.
_TIERS = [
    ("BLASTp",          "source_blastp"),
    ("keyword search",  "source_keyword"),
    ("all candidates",  None),                 # None = union of every method
]

# Author's table, 2026-09-04. Guards the four counts per row that the Results text
# quotes; OR and P are recomputed and compared to 2 d.p. / 3 d.p. respectively.
_EXPECTED = {
    "BLASTp":         dict(proteins=28, loci=9,  with_=4, without=5, or_=1.77, p=0.685),
    "keyword search": dict(proteins=58, loci=15, with_=8, without=7, or_=4.36, p=0.071),
    "all candidates": dict(proteins=69, loci=18, with_=9, without=9, or_=4.46, p=0.075),
}

# Header wrapping is manual: these are long noun phrases and matplotlib does no wrapping
# of its own, so the line breaks are chosen to keep each column narrow enough to fit.
_HEADERS = [
    "Acetyltransferase\ndetection method",
    "# candidate\nacetyltransferases",
    "# GWAS K-loci with\nacetyltransferase",
    "# GWAS K-loci with\nacetyltransferase and\ndeacetylase predictor",
    "# GWAS K-loci with\nacetyltransferase and\nwithout deacetylase\npredictor",
    "OR",
    "$\\it{p}$",
]
_COL_W = [0.200, 0.130, 0.150, 0.175, 0.185, 0.080, 0.080]

_C_HEADER_BG = "#3a3a3a"
_C_ROW_BG    = "#ffffff"
_C_HIGHLIGHT = "#dbeaf4"   # tint of sslbh_color — marks the union row
_C_BORDER    = "#999999"


def _locus_of(protein_id: str) -> str | None:
    """'KL61_12_wcqZ' -> 'KL61'. Reference proteins (PROTEIN0x) have no locus."""
    m = re.match(r"(KL\d+)", str(protein_id))
    return m.group(1) if m else None


def _truthy(series: pd.Series) -> pd.Series:
    """The source_* columns are written as real booleans, but survive a round-trip
    through TSV as the strings 'True'/'False'. Accept both."""
    if series.dtype == bool:
        return series
    return series.astype(str).str.strip().str.lower().isin({"true", "yes", "1", "1.0", "y"})


def _locus_num(locus: str) -> int:
    return int(locus[2:])


def _build_rows(acetylases_kloci_tsv: Path,
                deacetylases_gwas_tsv: Path,
                pyseer_hits_tsv: Path) -> list[dict]:
    pyseer = pd.read_csv(pyseer_hits_tsv, sep="\t")
    pyseer = pyseer[pyseer["mode"] == "lasso"]
    universe = set(pyseer["locus"].unique())

    kl = pd.read_csv(acetylases_kloci_tsv, sep="\t")
    kl["locus"] = kl["protein_id"].map(_locus_of)

    deac = pd.read_csv(deacetylases_gwas_tsv, sep="\t")
    predictor_loci = set(deac["locus"]) & universe
    n_pred = len(predictor_loci)
    n_nopred = len(universe) - n_pred

    rows = []
    for label, column in _TIERS:
        mask = pd.Series(True, index=kl.index) if column is None else _truthy(kl[column])
        sub = kl[mask]
        called = set(sub["locus"].dropna()) & universe

        a = len(called & predictor_loci)          # predictor +, acetyltransferase +
        b = n_pred - a                            # predictor +, acetyltransferase -
        c = len(called - predictor_loci)          # predictor -, acetyltransferase +
        d = n_nopred - c                          # predictor -, acetyltransferase -

        _, p = fisher_exact([[a, b], [c, d]])
        res = odds_ratio([[a, b], [c, d]])
        ci = res.confidence_interval(0.95)

        rows.append(dict(
            label=label, proteins=len(sub), loci=len(called),
            a=a, n_pred=n_pred, c=c, n_nopred=n_nopred,
            or_=res.statistic, ci_low=ci.low, ci_high=ci.high, p=p,
            highlight=(column is None),
        ))
    return rows


def _verify(rows: list[dict]) -> None:
    """Fail loudly if the pipeline no longer reproduces the author's table.

    Figure 3.2C has the same guard for the same reason: a table panel whose numbers
    drift away from the Results text is worse than no panel, because nothing about the
    rendered image says it has gone stale.
    """
    problems = []
    for row in rows:
        exp = _EXPECTED.get(row["label"])
        if exp is None:
            problems.append(f"{row['label']}: no expected values recorded")
            continue
        for key, got, want, tol in [
            ("proteins", row["proteins"], exp["proteins"], 0),
            ("K-loci",   row["loci"],     exp["loci"],     0),
            ("with",     row["a"],        exp["with_"],    0),
            ("without",  row["c"],        exp["without"],  0),
            ("OR",       row["or_"],      exp["or_"],      0.005),
            ("P",        row["p"],        exp["p"],        0.0005),
        ]:
            if abs(got - want) > tol:
                problems.append(f"{row['label']} {key}: expected {want}, computed {got}")
    if problems:
        raise ValueError(
            "Figure 4.3B no longer reproduces the recorded detection-tier table "
            f"({'; '.join(problems)}). Reconcile the Results text and _EXPECTED with the "
            "pipeline before regenerating the figure."
        )


def plot_figure4_3_panelB(
    acetylases_kloci_tsv: Path,
    deacetylases_gwas_tsv: Path,
    pyseer_hits_tsv: Path,
    plots_dir: Path,
    style=None,
) -> None:
    """
    Produce Figure 4.3 Panel B: acetyltransferase detection tier vs SGNH predictor.

    Args:
        acetylases_kloci_tsv:  cps_acetylases/acetylases_kloci.tsv (69 candidates, source_* flags)
        deacetylases_gwas_tsv: rbp_deacetylases/deacetylases_gwas.tsv (12 SGNH predictor loci)
        pyseer_hits_tsv:       input_dir/data-gwas/3_GWAS/3_PROCESSING/pyseer_hits_filtered.tsv;
                               source of the 35-locus universe, as in panel A
        plots_dir:             output directory for plots
        style:                 cfg.style (optional)
    """
    dpi = getattr(style, "dpi", 300)

    rows = _build_rows(acetylases_kloci_tsv, deacetylases_gwas_tsv, pyseer_hits_tsv)
    _verify(rows)

    header_fs, cell_fs = 7.0, 8.0
    fig_w = 7.5
    row_h, header_h = 0.30, 0.68
    fig_h = header_h + len(rows) * row_h + 0.10

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    total_h = header_h + len(rows) * row_h
    x_starts, x = [], 0.0
    for w in _COL_W:
        x_starts.append(x)
        x += w

    def _cell(x0, y0, w, h, facecolor):
        ax.add_patch(mpatches.Rectangle(
            (x0, y0), w, h, facecolor=facecolor,
            edgecolor=_C_BORDER, linewidth=0.5, clip_on=False))

    # --- header band ---
    y_head = 1.0 - header_h / total_h
    for label, x0, w in zip(_HEADERS, x_starts, _COL_W):
        _cell(x0, y_head, w, header_h / total_h, _C_HEADER_BG)
        ax.text(x0 + w / 2, y_head + header_h / total_h / 2, label,
                ha="center", va="center", fontsize=header_fs,
                fontweight="bold", color="white", linespacing=1.35)

    # --- data rows ---
    for i, row in enumerate(rows):
        y0 = 1.0 - (header_h + (i + 1) * row_h) / total_h
        h = row_h / total_h
        bg = _C_HIGHLIGHT if row["highlight"] else _C_ROW_BG
        pct_with = 100 * row["a"] / row["n_pred"]
        pct_without = 100 * row["c"] / row["n_nopred"]
        values = [
            row["label"],
            f"{row['proteins']}",
            f"{row['loci']}",
            f"{row['a']}/{row['n_pred']} ({pct_with:.0f}%)",
            f"{row['c']}/{row['n_nopred']} ({pct_without:.0f}%)",
            f"{row['or_']:.1f}",
            f"{row['p']:.2f}" if row["p"] >= 0.1 else f"{row['p']:.3f}",
        ]
        for j, (val, x0, w) in enumerate(zip(values, x_starts, _COL_W)):
            _cell(x0, y0, w, h, bg)
            ax.text(x0 + (0.012 if j == 0 else w / 2), y0 + h / 2, val,
                    ha="left" if j == 0 else "center", va="center",
                    fontsize=cell_fs, color="#222222",
                    fontweight="bold" if row["highlight"] and j == 0 else "normal")

    plots_dir.mkdir(parents=True, exist_ok=True)
    plt.tight_layout(pad=0.1)
    for ext in ("png", "pdf"):
        out = plots_dir / f"figure4_3-panelB.{ext}"
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        print(f"  [figure4_3-panelB] → {out.name}")
    plt.close(fig)

    for row in rows:
        print(f"    {row['label']:19s} n={row['proteins']:2d} loci={row['loci']:2d} "
              f"{row['a']}/{row['n_pred']} vs {row['c']}/{row['n_nopred']} "
              f"OR={row['or_']:.2f} (95% CI {row['ci_low']:.2f}–{row['ci_high']:.2f}) "
              f"P={row['p']:.3f}")
