"""
Render phandango-style cladogram + heatmap PNGs for all GWAS PCs.

Reads phandango files written by gwas-proc Step 2h:
    gwas-data/{ecod_folder}/{locus}/{clustering_level}/{PC}/phandango/
        subtree.nwk
        variants.csv    (columns: genomeID, {locus}, {locus}:colour, {PC}, {PC}:colour)

Renders:
    gwas-data/{ecod_folder}/{locus}/{clustering_level}/{PC}/phandango/{locus}.png

Checkpoint: skips any PC whose PNG already exists.

K-locus confidence colours (from bacteria_metadata.tsv):
    Perfect / Very high  →  kpam_color  (#1f77b4)
    High                 →  #5a9fd4
    Good                 →  #8ec5e8
    Low                  →  #e8a0a0
    Absent / unknown     →  #FDFEFE
"""

from __future__ import annotations

from io import StringIO
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import pandas as pd
from Bio import Phylo


_ABSENT_COLOR = "#FDFEFE"

_CONF_COLOR: dict[str, str] = {
    "Perfect":   "#1f77b4",
    "Very high": "#1f77b4",
    "High":      "#5a9fd4",
    "Good":      "#8ec5e8",
    "Low":       "#e8a0a0",
}

# Only genomes with these confidence levels are retained in the plot
_ABOVE_GOOD = {"Perfect", "Very high", "High"}

_TP_COLOR = "#2ca02c"   # True positive  — locus present, PC detected
_FN_COLOR = "#d62728"   # False negative — locus present, PC absent
# False positive uses ecod-type color passed at render time

_ECOD_PC_COLOR_ATTRS = {
    "sgnh-ecod":  ("sgnh_domain_color", "#c9a227"),
    "ssrbh-ecod": ("ssrbh_color",       "#2ca02c"),
    "other-ecod": ("gray_color",        "#bfbfbf"),
    "no-ecod":    ("gray_color",        "#bfbfbf"),
}


# ---------------------------------------------------------------------------
# Tree layout
# ---------------------------------------------------------------------------

def _get_leaf_order(root_clade) -> list[str]:
    leaves: list[str] = []
    def _dfs(clade):
        if clade.is_terminal():
            leaves.append(clade.name)
        else:
            for child in clade.clades:
                _dfs(child)
    _dfs(root_clade)
    return leaves


def _compute_positions(root_clade) -> dict[Any, tuple[float, float]]:
    def _max_depth(clade, d):
        if clade.is_terminal():
            return d
        return max(_max_depth(c, d + 1) for c in clade.clades)
    max_d = _max_depth(root_clade, 0)

    pos: dict[Any, tuple[float, float]] = {}
    counter = [0]

    def _assign(clade, depth: int):
        if clade.is_terminal():
            pos[id(clade)] = (max_d, counter[0])
            counter[0] += 1
        else:
            for child in clade.clades:
                _assign(child, depth + 1)
            ys = [pos[id(c)][1] for c in clade.clades]
            pos[id(clade)] = (depth, (min(ys) + max(ys)) / 2.0)

    _assign(root_clade, 0)
    return pos


def _draw_tree(ax: plt.Axes, root_clade, positions: dict, n_leaves: int) -> None:
    def _draw(clade):
        px, py = positions[id(clade)]
        for child in clade.clades:
            cx, cy = positions[id(child)]
            ax.plot([px, cx], [cy, cy], color="black", lw=0.5)
            _draw(child)
        if not clade.is_terminal():
            child_ys = [positions[id(c)][1] for c in clade.clades]
            ax.plot([px, px], [min(child_ys), max(child_ys)], color="black", lw=0.5)
        else:
            ax.text(px + 0.05, py, clade.name, va="center", ha="left",
                    fontsize=2.3, color="black")
    _draw(root_clade)
    ax.axis("off")
    ax.set_ylim(n_leaves - 0.5, -0.5)
    max_d = max(p[0] for p in positions.values())
    ax.set_xlim(0, max_d * 1.4)


def _draw_heatmap(
    ax: plt.Axes,
    leaf_order: list[str],
    values: dict[str, float],
    color_present: str,
    conf_colors: dict[str, str] | None = None,
) -> None:
    n = len(leaf_order)
    for i, gid in enumerate(leaf_order):
        val = values.get(gid, 0)
        if conf_colors is not None:
            color = conf_colors.get(gid, _ABSENT_COLOR)
        elif val:
            color = color_present
        else:
            color = _ABSENT_COLOR
        ax.add_patch(plt.Rectangle((0, n - 1 - i), 1, 1, color=color, linewidth=0))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, n)
    ax.axis("off")


def _draw_tp_fn_fp_heatmap(
    ax: plt.Axes,
    leaf_order: list[str],
    keep_genomes: set[str],
    pc_vals: dict[str, float],
    fp_color: str,
) -> None:
    """
    One filled rectangle per genome, coloured by TP/FN/FP classification.
    "Locus present" = genome in keep_genomes (>Good confidence only).
        locus=1, pc=1  → TP  (green)
        locus=1, pc=0  → FN  (red)
        locus=0, pc=1  → FP  (teal)
        locus=0, pc=0  → absent (white)
    """
    n = len(leaf_order)
    for i, gid in enumerate(leaf_order):
        locus = gid in keep_genomes
        pc    = bool(pc_vals.get(gid, 0))
        if locus and pc:
            color = _TP_COLOR
        elif locus and not pc:
            color = _FN_COLOR
        elif not locus and pc:
            color = fp_color
        else:
            color = _ABSENT_COLOR
        ax.add_patch(plt.Rectangle((0, n - 1 - i), 1, 1, color=color, linewidth=0))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, n)
    ax.axis("off")


# ---------------------------------------------------------------------------
# Single-PC render
# ---------------------------------------------------------------------------

def _ecod_pc_color(ecod_type: str, style) -> str:
    attr, default = _ECOD_PC_COLOR_ATTRS.get(ecod_type, ("gray_color", "#bfbfbf"))
    return getattr(style, attr, default) if style is not None else default


def _render_one(
    locus: str,
    pc: str,
    ecod_type: str,
    phandango_dir: Path,
    conf_colors_locus: dict[str, str],
    dpi: int,
    style=None,
) -> None:
    nwk_path = phandango_dir / "subtree.nwk"
    csv_path = phandango_dir / "variants.csv"
    out_png  = phandango_dir / f"{locus}.png"

    tree      = Phylo.read(StringIO(nwk_path.read_text()), "newick")
    root      = tree.root
    leaves    = _get_leaf_order(root)
    positions = _compute_positions(root)
    n         = len(leaves)

    # keep_genomes: >Good confidence only (used for K-locus column colouring and TP/FN/FP)
    keep_genomes = set(conf_colors_locus.keys())

    variants = pd.read_csv(csv_path)
    variants["genomeID"] = variants["genomeID"].astype(str)
    locus_vals = dict(zip(variants["genomeID"],
                          variants.get(locus, pd.Series(dtype=float)).fillna(0)))
    pc_vals    = dict(zip(variants["genomeID"],
                          variants.get(pc,    pd.Series(dtype=float)).fillna(0)))

    locus_conf = {gid: conf_colors_locus.get(gid, _ABSENT_COLOR) for gid in leaves}
    pc_color   = _ecod_pc_color(ecod_type, style)
    fp_color   = "#17becf"

    fig = plt.figure(figsize=(2.8, max(4, n * 0.046)))
    gs  = gridspec.GridSpec(
        2, 4,
        width_ratios=[7, 1.5, 1.5, 1.5],
        height_ratios=[0.15, 1],
        hspace=0.01, wspace=0.12,
        left=0.01, right=0.99, top=0.97, bottom=0.01,
    )
    ax_h_locus = fig.add_subplot(gs[0, 1])
    ax_h_pc    = fig.add_subplot(gs[0, 2])
    ax_h_tp    = fig.add_subplot(gs[0, 3])
    ax_tree    = fig.add_subplot(gs[1, 0])
    ax_locus   = fig.add_subplot(gs[1, 1])
    ax_pc      = fig.add_subplot(gs[1, 2])
    ax_tp      = fig.add_subplot(gs[1, 3])

    _draw_tree(ax_tree, root, positions, n)
    _draw_heatmap(ax_locus, leaves, locus_vals, "#1f77b4", conf_colors=locus_conf)
    _draw_heatmap(ax_pc,    leaves, pc_vals,    pc_color)
    _draw_tp_fn_fp_heatmap(ax_tp, leaves, keep_genomes, pc_vals, fp_color=fp_color)

    for ax_h, label, color in [
        (ax_h_locus, locus,      "#1f77b4"),
        (ax_h_pc,    pc,         pc_color),
        (ax_h_tp,    "TP/FN/FP", "#555555"),
    ]:
        ax_h.set_xlim(0, 1); ax_h.set_ylim(0, 1); ax_h.axis("off")
        ax_h.text(0.5, 0.5, label, ha="center", va="center",
                  fontsize=12, fontweight="bold", color="black",
                  rotation=90, transform=ax_h.transAxes)

    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def render_all_phandango(
    gwas_hits_tsv: Path,
    gwas_data_dir: Path,
    meta_tsv: Path,
    dpi: int = 150,
    style=None,
) -> None:
    """
    Render phandango PNG for every PC in gwas_hits.tsv.

    Skips any PC whose PNG already exists (checkpoint).

    Args:
        gwas_hits_tsv: processing/gwas-data/gwas_hits.tsv
        gwas_data_dir: output_dir/processing/gwas-data
        meta_tsv:      bacteria_metadata.tsv (for Kaptive confidence colours)
        dpi:           output resolution (default 150)
    """
    hits = pd.read_csv(gwas_hits_tsv, sep="\t")
    rows = hits.drop_duplicates(["locus", "clustering_level", "PC", "ecod_folder"])

    # Build confidence colour map once — only >Good confidence genomes retained
    # conf_map: {locus: {genomeID: hex}}  (only Perfect/Very high/High)
    conf_map: dict[str, dict[str, str]] = {}
    if meta_tsv.exists():
        meta = pd.read_csv(meta_tsv, sep="\t",
                           usecols=["genomeID", "K_locus", "K_locus_confidence"])
        meta["genomeID"] = meta["genomeID"].astype(str)
        for locus, grp in meta.groupby("K_locus"):
            conf_map[locus] = {
                row["genomeID"]: _CONF_COLOR[row["K_locus_confidence"]]
                for _, row in grp.iterrows()
                if row["K_locus_confidence"] in _ABOVE_GOOD
            }

    total = len(rows)
    ok = skipped = errors = 0

    for _, row in rows.iterrows():
        locus, cl, pc = row["locus"], row["clustering_level"], row["PC"]
        ecod_folder   = row["ecod_folder"]
        ecod_type     = row["ecod_type"]

        phandango_dir = gwas_data_dir / ecod_folder / locus / cl / pc / "phandango"
        out_png       = phandango_dir / f"{locus}.png"

        if out_png.exists():
            skipped += 1
            continue

        if not (phandango_dir / "subtree.nwk").exists():
            errors += 1
            print(f"  [missing] {locus}/{cl}/{pc}: phandango files not exported yet")
            continue

        try:
            _render_one(locus, pc, ecod_type, phandango_dir,
                        conf_map.get(locus, {}), dpi, style)
            ok += 1
        except Exception as exc:
            errors += 1
            print(f"  [error] {locus}/{cl}/{pc}: {exc}")

    print(f"  Phandango PNGs: {ok} rendered, {skipped} skipped (checkpoint), "
          f"{errors} errors — {total} total PCs")
