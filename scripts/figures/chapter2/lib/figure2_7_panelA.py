"""
Chapter 2 — Figure 2.7 Panel A: Active depolymerase pairwise BLASTP comparison.

Compares pairwise amino-acid sequence identity between active depolymerases
from virulent phages (virulent_active TSV) and temperate phages (manualsearch +
gwas_active TSVs). Points are classified into six categories by specificity match
and source pairing.

Two subplots (rows):
  - Row 1: bidirectional coverage >= 80%
  - Row 2: bidirectional coverage >= 50%
Both rows use evalue <= 1e-3, pident >= 20%.

Reads:
    virulent_active_tsv         — output_dir/rbp_depolymerases/depolymerases_virulent_active.tsv
    gwas_tsv                    — output_dir/rbp_depolymerases/depolymerases_gwas.tsv
    manualsearch_active_tsv     — output_dir/rbp_depolymerases/depolymerases_manualsearch_active.tsv
    manualsearch_inactive_tsv   — output_dir/rbp_depolymerases/depolymerases_manualsearch_inactive.tsv
    manualsearch_notproduced_tsv — output_dir/rbp_depolymerases/depolymerases_manualsearch_notproduced.tsv
    gwas_active_tsv             — output_dir/rbp_depolymerases/depolymerases_gwas_active.tsv

Writes:
    analysis_dir/active_enzymes.tsv      — cleaned combined table
    analysis_dir/blastp_all_hits.tsv     — all annotated BLASTP hits
    plots_dir/figure2_7-panelA.png / .pdf
    plots_dir/legends/figure2_7-panelA-legend.png / .pdf
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

EVALUE_CUTOFF           = 1e-3
PIDENT_CUTOFF           = 20.0
REACHES_END_TOLERANCE   = 50

_ACTIVE_SOURCES = {"virulent", "prophage"}

# (label, specificity_type, source_type)
_CATEGORIES = [
    ("same-specificity\n(all)",                      "same",     "all"),
    ("same-specificity\n(prophage)",                 "same",     "prophage-prophage"),
    ("same-specificity\n(virulent)",                 "same",     "virulent-virulent"),
    ("same-specificity\n(prophage-virulent)",        "same",     "prophage-virulent"),
    ("same-specificity\n(prediction)",               "same",     "prediction"),
    ("same-specificity\n(prophage\nnot produced)",   "same",     "prophage_notproduced"),
    ("distinct-specificity\n(all)",                      "distinct", "all"),
    ("distinct-specificity\n(prophage)",                 "distinct", "prophage-prophage"),
    ("distinct-specificity\n(virulent)",                 "distinct", "virulent-virulent"),
    ("distinct-specificity\n(prophage-virulent)",        "distinct", "prophage-virulent"),
    ("distinct-specificity\n(prediction)",               "distinct", "prediction"),
    ("distinct-specificity\n(prophage\nnot produced)",   "distinct", "prophage_notproduced"),
]

_SAME_LABELS     = [c[0] for c in _CATEGORIES if c[1] == "same"]
_DISTINCT_LABELS = [c[0] for c in _CATEGORIES if c[1] == "distinct"]
_CAT_LABELS = [c[0] for c in _CATEGORIES]

# (cov_thresh, id_min, id_max, label)
_COV_ROWS = [
    (0.50, 0.20, 1.00, "bidirectional coverage ≥50%, sequence identity ≥20%"),
]

_DEFAULT_PROPHAGE_COLOR     = "#9900aa"   # dark magenta
_DEFAULT_VIRULENT_COLOR     = "#45087a"   # dark purple
_DEFAULT_GRAY_COLOR         = "#991b1b"   # dark red (cross-source / both)
_DEFAULT_INACTIVE_COLOR     = "#5c0000"   # dark deep red (prophage inactive)
_DEFAULT_NOTPRODUCED_COLOR  = "#4d4d4d"   # dark gray (prophage not produced)
_DEFAULT_PREDICTION_COLOR   = "#1a7a1a"   # dark green (GWAS predictions)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _normalize_klocus(val) -> str | None:
    if pd.isna(val) or str(val).strip() in ("", "—", "-"):
        return None
    val = str(val).strip()
    if val.upper().startswith("KL") or val.upper().startswith("KN"):
        return val
    if val.upper().startswith("K"):
        return "KL" + val[1:]
    return val


def _load_virulent(path):
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={"protein_seq": "seq"})
    df["source"] = "virulent"
    df["K_locus_specificity"] = df["K_locus_specificity"].map(_normalize_klocus)
    df = df.dropna(subset=["K_locus_specificity"])
    return df[["proteinID", "source", "K_locus_specificity", "seq"]].drop_duplicates(subset=["proteinID"]).reset_index(drop=True)


def _load_gwas_predictions(path):
    df = pd.read_csv(path, sep="\t")
    df["source"] = "prediction"
    df["K_locus_specificity"] = df["locus"].map(_normalize_klocus)
    df = df.dropna(subset=["K_locus_specificity"])
    return df[["proteinID", "source", "K_locus_specificity", "seq"]].drop_duplicates(subset=["proteinID"]).reset_index(drop=True)


def _load_prophage_from_tsv(path, source_label, use_host_as_specificity=False):
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={"protein seq": "seq"})
    df["source"] = source_label
    spec_col = "K_locus_host" if use_host_as_specificity else "K_locus_specificity"
    df["K_locus_specificity"] = df[spec_col].map(_normalize_klocus)
    df = df.dropna(subset=["K_locus_specificity"])
    return df[["proteinID", "source", "K_locus_specificity", "seq"]].drop_duplicates(subset=["proteinID"]).reset_index(drop=True)


# ---------------------------------------------------------------------------
# BLASTP
# ---------------------------------------------------------------------------

def _run_blastp(combined: pd.DataFrame, tmp_dir: Path) -> pd.DataFrame:
    tmp_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = tmp_dir / "active_depo.fasta"
    blast_out  = tmp_dir / "active_depo_blastp.tsv"

    with open(fasta_path, "w") as fh:
        for _, row in combined.iterrows():
            fh.write(f">{row['proteinID']}\n{row['seq']}\n")

    subprocess.run(
        ["makeblastdb", "-in", str(fasta_path), "-dbtype", "prot"],
        check=True, capture_output=True,
    )
    subprocess.run(
        [
            "blastp",
            "-query", str(fasta_path),
            "-db",    str(fasta_path),
            "-out",   str(blast_out),
            "-outfmt", "6 qseqid sseqid evalue pident bitscore qlen qstart qend slen sstart send",
            "-evalue", "10",
        ],
        check=True, capture_output=True,
    )

    cols = ["query", "target", "evalue", "pident", "bitscore",
            "qlen", "qstart", "qend", "slen", "sstart", "send"]
    return pd.read_csv(blast_out, sep="\t", names=cols, header=None)


# ---------------------------------------------------------------------------
# Annotation
# ---------------------------------------------------------------------------

def _annotate_hits(raw: pd.DataFrame, combined: pd.DataFrame) -> pd.DataFrame:
    spec_map   = dict(zip(combined["proteinID"], combined["K_locus_specificity"]))
    source_map = dict(zip(combined["proteinID"], combined["source"]))

    hits = raw[
        (raw["evalue"] <= EVALUE_CUTOFF) &
        (raw["pident"] >= PIDENT_CUTOFF) &
        (raw["query"] != raw["target"])
    ].copy()

    # Merge adjacent HSPs for the same (query, target) where gap < 100 aa in both sequences
    merged_rows = []
    for (q, t), grp in hits.groupby(["query", "target"]):
        grp = grp.sort_values("qstart").reset_index(drop=True)
        current = grp.loc[0].to_dict()
        for i in range(1, len(grp)):
            nxt = grp.loc[i].to_dict()
            q_gap = nxt["qstart"] - current["qend"] - 1
            s_gap = nxt["sstart"] - current["send"] - 1
            if q_gap < 100 and s_gap < 100:
                len1 = current["qend"] - current["qstart"] + 1
                len2 = nxt["qend"]    - nxt["qstart"]    + 1
                current["pident"]   = (current["pident"] * len1 + nxt["pident"] * len2) / (len1 + len2)
                current["qend"]     = nxt["qend"]
                current["send"]     = nxt["send"]
                current["bitscore"] = current["bitscore"] + nxt["bitscore"]
                current["evalue"]   = min(current["evalue"], nxt["evalue"])
            else:
                merged_rows.append(current)
                current = nxt
        merged_rows.append(current)
    hits = pd.DataFrame(merged_rows)

    # Keep best-bitscore hit per unordered pair
    hits["pair_key"] = [
        frozenset([q, t]) for q, t in zip(hits["query"], hits["target"])
    ]
    hits = hits.sort_values("bitscore", ascending=False)
    hits = hits.drop_duplicates(subset=["pair_key"])

    hits["qcov"]         = (hits["qend"] - hits["qstart"] + 1) / hits["qlen"]
    hits["scov"]         = (hits["send"] - hits["sstart"] + 1) / hits["slen"]
    hits["bidir_cov_min"] = hits[["qcov", "scov"]].min(axis=1)
    hits["reaches_end"]  = (
        (hits["qend"] >= hits["qlen"] - REACHES_END_TOLERANCE) &
        (hits["send"] >= hits["slen"] - REACHES_END_TOLERANCE)
    )

    hits["q_specificity"]    = hits["query"].map(spec_map)
    hits["s_specificity"]    = hits["target"].map(spec_map)
    hits["q_source"]         = hits["query"].map(source_map)
    hits["s_source"]         = hits["target"].map(source_map)
    def _spec_overlap(a, b) -> bool:
        a_parts = {_normalize_klocus(p.strip()) for p in str(a).split("/")}
        b_parts = {_normalize_klocus(p.strip()) for p in str(b).split("/")}
        return bool(a_parts & b_parts)

    hits["specificity_match"] = hits.apply(
        lambda r: _spec_overlap(r["q_specificity"], r["s_specificity"]), axis=1
    )
    def _source_pair(q, s):
        if q == "prophage_inactive" and s in _ACTIVE_SOURCES:
            return "prophage_inactive"
        if s == "prophage_inactive" and q in _ACTIVE_SOURCES:
            return "prophage_inactive"
        if q == "prophage_notproduced" and s in _ACTIVE_SOURCES:
            return "prophage_notproduced"
        if s == "prophage_notproduced" and q in _ACTIVE_SOURCES:
            return "prophage_notproduced"
        if q == "prediction" and s in _ACTIVE_SOURCES:
            return "prediction"
        if s == "prediction" and q in _ACTIVE_SOURCES:
            return "prediction"
        if q == "prophage" and s == "prophage":
            return "prophage-prophage"
        if q == "virulent" and s == "virulent":
            return "virulent-virulent"
        if q in _ACTIVE_SOURCES and s in _ACTIVE_SOURCES:
            return "prophage-virulent"
        return "other"

    hits["source_pair"] = hits.apply(lambda r: _source_pair(r["q_source"], r["s_source"]), axis=1)
    hits["pident_frac"] = hits["pident"] / 100.0
    return hits.drop(columns=["pair_key"]).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Category assignment — one row per (hit × category) it belongs to
# ---------------------------------------------------------------------------

def _assign_categories(hits: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, hit in hits.iterrows():
        for label, spec_type, src_type in _CATEGORIES:
            if spec_type == "same"     and not hit["specificity_match"]:
                continue
            if spec_type == "distinct" and hit["specificity_match"]:
                continue
            if src_type == "all":
                # "all" covers only active-active pairs (unchanged from original)
                if hit["q_source"] not in _ACTIVE_SOURCES or hit["s_source"] not in _ACTIVE_SOURCES:
                    continue
            elif hit["source_pair"] != src_type:
                continue
            d = hit.to_dict()
            d["category"] = label
            rows.append(d)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def _plot(cat_hits: pd.DataFrame, plots_dir: Path, style=None) -> None:
    axis_label_fs  = getattr(style, "axis_label_fontsize",   12) * 1.2
    axis_label_fw  = getattr(style, "axis_label_fontweight", "bold")
    tick_fs        = getattr(style, "tick_fontsize",           8) * 1.2
    tick_fw        = getattr(style, "tick_fontweight",        "bold")
    dpi            = getattr(style, "dpi",                    300)
    prophage_color    = _DEFAULT_PROPHAGE_COLOR
    virulent_color    = _DEFAULT_VIRULENT_COLOR
    gray_color        = _DEFAULT_GRAY_COLOR
    notproduced_color = _DEFAULT_NOTPRODUCED_COLOR
    prediction_color  = _DEFAULT_PREDICTION_COLOR

    color_map = {
        "prophage-prophage":    prophage_color,
        "virulent-virulent":    virulent_color,
        "prophage-virulent":    gray_color,
        "prediction":           prediction_color,
        "prophage_notproduced": notproduced_color,
    }

    cov_thresh, id_min, id_max, cov_label = _COV_ROWS[0]
    subset = cat_hits[
        (cat_hits["bidir_cov_min"] >= cov_thresh) &
        (cat_hits["pident_frac"]   >= id_min) &
        (cat_hits["pident_frac"]   <= id_max)
    ]

    np.random.seed(42)

    # Left column: active pairs (all / prophage / virulent / prophage-virulent)
    # Right column: uncharacterised proteins (prediction / prophage not produced)
    left_tick_labels  = ["all", "prophage", "virulent", "prophage-\nvirulent"]
    right_tick_labels = ["prediction", "prophage\nnot produced"]

    spec_rows = [
        ("same capsule type\nspecificity",       _SAME_LABELS[:4],     _SAME_LABELS[4:]),
        ("distinct capsule type\nspecificities", _DISTINCT_LABELS[:4], _DISTINCT_LABELS[4:]),
    ]

    fig, axes = plt.subplots(
        2, 2, figsize=(11, 4.5),
        gridspec_kw={"width_ratios": [4, 2]},
        sharex="col", sharey="row",
    )

    def _draw_panel(ax, cat_labels):
        x_positions = {lbl: i for i, lbl in enumerate(cat_labels)}
        box_data, box_positions = [], []
        n_per_cat = {}
        for cat_label in cat_labels:
            cat_data = subset[subset["category"] == cat_label]
            n_per_cat[cat_label] = len(cat_data)
            if cat_data.empty:
                continue
            x_base  = x_positions[cat_label]
            x_vals  = x_base + np.random.uniform(-0.25, 0.25, size=len(cat_data))
            y_vals  = cat_data["pident_frac"].values
            colors  = [color_map[sp] for sp in cat_data["source_pair"]]
            markers = ["o" if re else "^" for re in cat_data["reaches_end"]]
            for x, y, c, m in zip(x_vals, y_vals, colors, markers):
                ax.scatter(x, y, color=c, marker=m, s=60, alpha=0.7,
                           edgecolors="none", zorder=3)
            box_data.append(y_vals)
            box_positions.append(x_base)
        if box_data:
            ax.boxplot(
                box_data, positions=box_positions, widths=0.45,
                patch_artist=True, showfliers=False, zorder=4,
                boxprops=dict(facecolor="none", edgecolor="gray", linewidth=1.2),
                whiskerprops=dict(color="gray", linewidth=1.0),
                capprops=dict(color="gray", linewidth=1.2),
                medianprops=dict(color="gray", linewidth=2.0),
            )
        ax.set_ylim(0, 1.05)
        ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax.tick_params(axis="y", labelsize=tick_fs)
        for lbl in ax.get_yticklabels():
            lbl.set_fontweight(tick_fw)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.4, color="gray")
        ax.set_axisbelow(True)
        return n_per_cat

    # n_counts[row_idx][col_idx] = {cat_label: n}
    n_counts = {}
    for row_idx, (row_title, left_labels, right_labels) in enumerate(spec_rows):
        n_left  = _draw_panel(axes[row_idx, 0], left_labels)
        n_right = _draw_panel(axes[row_idx, 1], right_labels)
        n_counts[row_idx] = {0: n_left, 1: n_right}
        axes[row_idx, 1].yaxis.set_label_position("right")
        axes[row_idx, 1].set_ylabel(row_title, fontsize=tick_fs + 1, fontweight="bold",
                                    rotation=0, labelpad=60, va="center")

    # cat_labels[row_idx][col_idx] → list of category label strings in x order
    cat_labels_by_row_col = {
        0: {0: _SAME_LABELS[:4],      1: _SAME_LABELS[4:]},
        1: {0: _DISTINCT_LABELS[:4],  1: _DISTINCT_LABELS[4:]},
    }

    for col_idx, base_tick_labels in enumerate([left_tick_labels, right_tick_labels]):
        # Bottom row: original bold category labels (unchanged)
        axes[1, col_idx].set_xticks(range(len(base_tick_labels)))
        axes[1, col_idx].set_xticklabels(base_tick_labels, fontsize=tick_fs, rotation=0, ha="center")
        for lbl in axes[1, col_idx].get_xticklabels():
            lbl.set_fontweight(tick_fw)

        # All rows: gray n= annotations just above the top of each panel
        for row_idx in range(2):
            cats = cat_labels_by_row_col[row_idx][col_idx]
            for i, cat in enumerate(cats):
                n = n_counts[row_idx][col_idx].get(cat, 0)
                axes[row_idx, col_idx].text(
                    i, 1.07, f"n={n}", ha="center", va="bottom",
                    fontsize=tick_fs - 2, color="#888888", clip_on=False,
                )

    axes[1, 0].set_xlabel("active depolymerases", fontsize=tick_fs + 1, fontweight="bold", labelpad=20)
    axes[1, 1].set_xlabel("putative depolymerases", fontsize=tick_fs + 1, fontweight="bold", labelpad=20)

    plots_dir.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    fig.text(0.04, 0.5, "receptor-binding domain\nsequence identity",
             va="center", ha="center", rotation=90,
             fontsize=axis_label_fs, fontweight=axis_label_fw)
    plt.subplots_adjust(left=0.13)
    for ext in ("png", "pdf"):
        out = plots_dir / f"figure2_7-panelA.{ext}"
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        print(f"  [figure2_7-panelA] → {out.name}")
    plt.close(fig)

    # Legend
    legend_handles = [
        plt.scatter([], [], marker="o", s=60, color=prophage_color,
                    alpha=0.8, edgecolors="none", label="Prophage–prophage pair"),
        plt.scatter([], [], marker="o", s=60, color=virulent_color,
                    alpha=0.8, edgecolors="none", label="Virulent–virulent pair"),
        plt.scatter([], [], marker="o", s=60, color=gray_color,
                    alpha=0.8, edgecolors="none", label="Prophage–virulent pair"),
        plt.scatter([], [], marker="o", s=60, color=prediction_color,
                    alpha=0.8, edgecolors="none", label="GWAS prediction vs active pair"),
        plt.scatter([], [], marker="o", s=60, color=notproduced_color,
                    alpha=0.8, edgecolors="none", label="Prophage not produced vs active pair"),
        plt.scatter([], [], marker="o", s=60, color="gray",
                    alpha=0.8, edgecolors="none", label="Alignment reaches C-terminal"),
        plt.scatter([], [], marker="^", s=60, color="gray",
                    alpha=0.8, edgecolors="none", label="Alignment does not reach C-terminal"),
    ]

    def _header(text):
        latex = text.replace(" ", "~")
        return mlines.Line2D([], [], color="none", label=f"$\\bf{{{latex}}}$")

    all_handles = (
        [_header("Source")] + legend_handles[:5] +
        [_header("C-terminal")] + legend_handles[5:]
    )
    fig_leg, ax_leg = plt.subplots(figsize=(4, 4.5))
    ax_leg.axis("off")
    ax_leg.legend(handles=all_handles, fontsize=tick_fs, loc="center",
                  ncol=1, frameon=False, handletextpad=0.5, labelspacing=0.4)

    legends_dir = plots_dir / "legends"
    legends_dir.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    for ext in ("png", "pdf"):
        out = legends_dir / f"figure2_7-panelA-legend.{ext}"
        fig_leg.savefig(out, dpi=dpi, bbox_inches="tight")
        print(f"  [figure2_7-panelA-legend] → legends/{out.name}")
    plt.close(fig_leg)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def plot_figure2_7_panelA(
    virulent_active_tsv: Path,
    gwas_tsv: Path,
    manualsearch_active_tsv: Path,
    manualsearch_inactive_tsv: Path,
    manualsearch_notproduced_tsv: Path,
    gwas_active_tsv: Path,
    plots_dir: Path,
    analysis_dir: Path,
    tmp_dir: Path,
    style=None,
) -> None:
    """
    Compare active depolymerases from virulent phages and temperate phages,
    plus GWAS predictions and inactive / not-produced prophage proteins.

    Args:
        virulent_active_tsv:          output_dir/rbp_depolymerases/depolymerases_virulent_active.tsv
        gwas_tsv:                     output_dir/rbp_depolymerases/depolymerases_gwas.tsv
        manualsearch_active_tsv:      output_dir/rbp_depolymerases/depolymerases_manualsearch_active.tsv
        manualsearch_inactive_tsv:    output_dir/rbp_depolymerases/depolymerases_manualsearch_inactive.tsv
        manualsearch_notproduced_tsv: output_dir/rbp_depolymerases/depolymerases_manualsearch_notproduced.tsv
        gwas_active_tsv:              output_dir/rbp_depolymerases/depolymerases_gwas_active.tsv
        plots_dir:                    output directory for plots
        analysis_dir:                 output directory for analysis TSVs
        tmp_dir:                      scratch directory for BLAST db and raw output
        style:                        cfg.style (optional)
    """
    analysis_dir = Path(analysis_dir)
    analysis_dir.mkdir(parents=True, exist_ok=True)

    print("  [figure2_7-panelA] Loading virulent active ...")
    s1 = _load_virulent(Path(virulent_active_tsv))
    print(f"    {len(s1)} proteins")

    print("  [figure2_7-panelA] Loading GWAS predictions ...")
    s3 = _load_gwas_predictions(Path(gwas_tsv))
    print(f"    {len(s3)} proteins")

    print("  [figure2_7-panelA] Loading prophage active ...")
    s4_manual = _load_prophage_from_tsv(Path(manualsearch_active_tsv), "prophage")
    s4_gwas   = _load_prophage_from_tsv(Path(gwas_active_tsv), "prophage")
    s4 = pd.concat([s4_manual, s4_gwas], ignore_index=True).drop_duplicates(subset=["proteinID"]).reset_index(drop=True)
    print(f"    {len(s4)} proteins")

    print("  [figure2_7-panelA] Loading prophage inactive + not produced ...")
    s4_inactive    = _load_prophage_from_tsv(Path(manualsearch_inactive_tsv), "prophage_inactive", use_host_as_specificity=True)
    s4_notproduced = _load_prophage_from_tsv(Path(manualsearch_notproduced_tsv), "prophage_notproduced", use_host_as_specificity=True)
    print(f"    inactive: {len(s4_inactive)}, not produced: {len(s4_notproduced)}")

    combined = pd.concat([s1, s3, s4, s4_inactive, s4_notproduced], ignore_index=True)
    combined_path = analysis_dir / "active_enzymes.tsv"
    combined.to_csv(combined_path, sep="\t", index=False)
    print(f"  [figure2_7-panelA] Combined: {len(combined)} proteins → {combined_path.name}")

    print("  [figure2_7-panelA] Running BLASTP all-vs-all ...")
    raw = _run_blastp(combined, Path(tmp_dir))
    print(f"    Raw BLAST rows: {len(raw)}")

    hits = _annotate_hits(raw, combined)
    hits_path = analysis_dir / "blastp_all_hits.tsv"
    hits.to_csv(hits_path, sep="\t", index=False)
    print(f"  [figure2_7-panelA] Annotated pairs: {len(hits)} → {hits_path.name}")

    cat_hits = _assign_categories(hits)

    # Clean table: one row per (pair × category), coverage_category = highest qualifying threshold
    clean = cat_hits[cat_hits["bidir_cov_min"] >= 0.50].copy()
    clean_out = pd.DataFrame({
        "protein1":             clean["query"].values,
        "protein2":             clean["target"].values,
        "protein1_specificity": clean["q_specificity"].values,
        "protein2_specificity": clean["s_specificity"].values,
        "seq_similarity":       clean["pident_frac"].values,
        "coverage_category":    ">=50%",
        "C-terminus":           ["included" if r else "not_included" for r in clean["reaches_end"]],
        "category":             clean["category"].str.replace("\n", " ", regex=False).values,
    })
    clean_path = analysis_dir / "figure2_7_panelA_table.tsv"
    clean_out.to_csv(clean_path, sep="\t", index=False)
    print(f"  [figure2_7-panelA] Clean table: {len(clean_out)} rows → {clean_path.name}")

    _plot(cat_hits, Path(plots_dir), style=style)
