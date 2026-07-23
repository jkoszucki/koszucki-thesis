"""
Export phandango input files (subtree + variants matrix) for each best SGNH predictor.

Source:
    gwas_path/3_GWAS/4_ANALYZE/1_PER_LOCUS/{locus}/lasso/phandango/{clustering_level}/
        _subtree.nwk     — phylogenetic subtree for this K-locus
        _variants.csv    — presence/absence for all significant PCs in this locus

Destination (one directory per best predictor):
    gwas-data/sgnh-ecod-reported-topology/{locus}/{clustering_level}/{PC}/phandango/
        subtree.nwk      — copied as-is
        variants.csv     — filtered: genomeID, {locus}, {locus}:colour, {PC}, {PC}:colour
                           leaves aligned to tree; absent cells → #FDFEFE

If the best-predictor PC is absent from _variants.csv, the binary matrix at
{clustering_level}/3_binary_matrix.tsv is used as a fallback.

Colour scheme:
    K-locus present  →  style.kpam_color        (default #1f77b4)
    PC present       →  style.sgnh_domain_color  (default #c9a227)
    absent           →  #FDFEFE
"""

from __future__ import annotations

import shutil
from io import StringIO
from pathlib import Path

import pandas as pd
from Bio import Phylo


_ABSENT_COLOR = "#FDFEFE"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _tree_leaves(nwk_path: Path) -> set[str]:
    tree = Phylo.read(StringIO(nwk_path.read_text()), "newick")
    return {c.name for c in tree.get_terminals()}


def _pc_from_binary_matrix(mmseqs_dir: Path, clustering_level: str, pc: str) -> pd.Series:
    """Return Series (index=genomeID) with 0/1 presence from binary matrix."""
    matrix_path = mmseqs_dir / clustering_level / "3_binary_matrix.tsv"
    row = pd.read_csv(matrix_path, sep="\t", index_col=0).loc[pc]
    row.name = pc
    return row


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def export_phandango(
    gwas_sgnh_best_tsv: Path,
    phandango_src_root: Path,
    mmseqs_dir: Path,
    gwas_data_dir: Path,
    style=None,
) -> None:
    """
    Export phandango subtree + variants for each best SGNH predictor.

    Args:
        gwas_sgnh_best_tsv:  sgnh-hydrolases/gwas_sgnh_best.tsv
        phandango_src_root:  gwas_path/3_GWAS/4_ANALYZE/1_PER_LOCUS
        mmseqs_dir:          gwas_path/3_GWAS/1_INTERMEDIATE/2_MMSEQS
        gwas_data_dir:       output_dir/processing/gwas-data
        style:               cfg.style (optional)
    """
    locus_color = getattr(style, "kpam_color",        "#1f77b4")
    pc_color    = getattr(style, "sgnh_domain_color", "#c9a227")

    best = pd.read_csv(gwas_sgnh_best_tsv, sep="\t")
    ok, skipped = 0, []

    for _, row in best.iterrows():
        locus, cl, pc = row["locus"], row["clustering_level"], row["PC"]

        src_dir  = phandango_src_root / locus / "lasso" / "phandango" / cl
        nwk_src  = src_dir / "_subtree.nwk"
        var_src  = src_dir / "_variants.csv"

        if not nwk_src.exists() or not var_src.exists():
            skipped.append(f"{locus}/{cl}/{pc}: source files missing in {src_dir}")
            continue

        # destination
        dest_dir = (
            gwas_data_dir
            / "sgnh-ecod-reported-topology"
            / locus / cl / pc / "phandango"
        )
        dest_dir.mkdir(parents=True, exist_ok=True)

        # --- subtree: copy as-is ---
        shutil.copy(nwk_src, dest_dir / "subtree.nwk")

        # --- variants: filter to genomeID + locus + PC columns ---
        leaves   = _tree_leaves(nwk_src)
        variants = pd.read_csv(var_src)
        variants["genomeID"] = variants["genomeID"].astype(str)
        variants = variants[variants["genomeID"].isin(leaves)]

        if locus not in variants.columns:
            skipped.append(f"{locus}/{cl}/{pc}: locus column '{locus}' missing in variants")
            continue

        out = variants[["genomeID", locus]].copy()
        out[f"{locus}:colour"] = out[locus].map({1: locus_color, 0: _ABSENT_COLOR})

        if pc in variants.columns:
            out[pc] = variants[pc]
        else:
            try:
                pc_series = _pc_from_binary_matrix(mmseqs_dir, cl, pc)
            except (KeyError, FileNotFoundError) as exc:
                skipped.append(f"{locus}/{cl}/{pc}: binary matrix fallback failed: {exc}")
                continue
            out[pc] = (
                pc_series.reindex(out["genomeID"].values).fillna(0).astype(int).values
            )

        out[f"{pc}:colour"] = out[pc].map({1: pc_color, 0: _ABSENT_COLOR})
        out.to_csv(dest_dir / "variants.csv", index=False)
        ok += 1
        print(f"    {locus}/{cl}/{pc} → phandango/ ({len(out)} leaves)")

    print(f"  Phandango export: {ok}/{len(best)} predictors written")
    if skipped:
        for s in skipped:
            print(f"  [skip] {s}")
