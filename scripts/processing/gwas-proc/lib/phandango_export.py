"""
Export phandango input files (subtree + variants matrix) for every PC in gwas_hits.tsv.

Source:
    gwas_path/3_GWAS/4_ANALYZE/1_PER_LOCUS/{locus}/lasso/phandango/{clustering_level}/
        _subtree.nwk     — phylogenetic subtree for this K-locus
        _variants.csv    — presence/absence for all significant PCs in this locus

Destination (one directory per unique locus × clustering_level × PC):
    gwas-data/{ecod_folder}/{locus}/{clustering_level}/{PC}/phandango/
        subtree.nwk      — copied as-is
        variants.csv     — filtered: genomeID, {locus}, {locus}:colour, {PC}, {PC}:colour
                           leaves aligned to tree; absent cells → #FDFEFE

Checkpoint: skips any destination already containing variants.csv.

If the PC is absent from _variants.csv, falls back to binary matrix at
{clustering_level}/3_binary_matrix.tsv.

Colour scheme:
    K-locus present  →  style.kpam_color        (default #1f77b4)
    PC present       →  style.sgnh_domain_color  (default #c9a227)  [SGNH]
                     →  style.ssrbh_color         (default #1f77b4)  [SSRBH]
                     →  style.sslbh_color         (default #9467bd)  [other-ecod / no-ecod SSLBH]
                     →  style.gray_color          (default #bfbfbf)  [other-ecod]
                     →  style.no_ecod_color       (default #ffffff)  [no-ecod]
    absent           →  #FDFEFE
"""

from __future__ import annotations

import shutil
from io import StringIO
from pathlib import Path

import pandas as pd
from Bio import Phylo


_ABSENT_COLOR = "#FDFEFE"

_ECOD_PC_COLORS = {
    "sgnh-ecod":  "sgnh_domain_color",
    "ssrbh-ecod": "ssrbh_color",
    "other-ecod": "gray_color",
    "no-ecod":    "no_ecod_color",
}
_ECOD_PC_DEFAULTS = {
    "sgnh-ecod":  "#c9a227",
    "ssrbh-ecod": "#1f77b4",
    "other-ecod": "#bfbfbf",
    "no-ecod":    "#ffffff",
}


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


def _pc_color(ecod_type: str, style) -> str:
    attr    = _ECOD_PC_COLORS.get(ecod_type, "gray_color")
    default = _ECOD_PC_DEFAULTS.get(ecod_type, "#bfbfbf")
    return getattr(style, attr, default)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def export_phandango(
    gwas_hits_tsv: Path,
    phandango_src_root: Path,
    mmseqs_dir: Path,
    gwas_data_dir: Path,
    style=None,
) -> None:
    """
    Export phandango subtree + variants for every PC in gwas_hits.tsv.

    Args:
        gwas_hits_tsv:      gwas-data/gwas_hits.tsv
        phandango_src_root: gwas_path/3_GWAS/4_ANALYZE/1_PER_LOCUS
        mmseqs_dir:         gwas_path/3_GWAS/1_INTERMEDIATE/2_MMSEQS
        gwas_data_dir:      output_dir/processing/gwas-data
        style:              cfg.style (optional)
    """
    locus_color = getattr(style, "kpam_color", "#1f77b4")

    hits = pd.read_csv(gwas_hits_tsv, sep="\t")
    rows = hits.drop_duplicates(["locus", "clustering_level", "PC", "ecod_folder"])

    total  = len(rows)
    ok     = 0
    skipped_missing = 0
    skipped_checkpoint = 0

    for _, row in rows.iterrows():
        locus, cl, pc = row["locus"], row["clustering_level"], row["PC"]
        ecod_folder   = row["ecod_folder"]
        ecod_type     = row["ecod_type"]

        dest_dir = gwas_data_dir / ecod_folder / locus / cl / pc / "phandango"

        # checkpoint
        if (dest_dir / "variants.csv").exists():
            skipped_checkpoint += 1
            continue

        src_dir = phandango_src_root / locus / "lasso" / "phandango" / cl
        nwk_src = src_dir / "_subtree.nwk"
        var_src = src_dir / "_variants.csv"

        if not nwk_src.exists() or not var_src.exists():
            skipped_missing += 1
            print(f"  [missing] {locus}/{cl}/{pc}: source not found in {src_dir}")
            continue

        dest_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy(nwk_src, dest_dir / "subtree.nwk")

        leaves   = _tree_leaves(nwk_src)
        variants = pd.read_csv(var_src)
        variants["genomeID"] = variants["genomeID"].astype(str)
        variants = variants[variants["genomeID"].isin(leaves)]

        if locus not in variants.columns:
            print(f"  [skip] {locus}/{cl}/{pc}: locus column missing in variants")
            skipped_missing += 1
            continue

        pc_color = _pc_color(ecod_type, style)

        out = variants[["genomeID", locus]].copy()
        out[f"{locus}:colour"] = out[locus].map({1: locus_color, 0: _ABSENT_COLOR})

        if pc in variants.columns:
            out[pc] = variants[pc]
        else:
            try:
                pc_series = _pc_from_binary_matrix(mmseqs_dir, cl, pc)
            except (KeyError, FileNotFoundError) as exc:
                print(f"  [skip] {locus}/{cl}/{pc}: binary matrix fallback failed: {exc}")
                skipped_missing += 1
                continue
            out[pc] = (
                pc_series.reindex(out["genomeID"].values).fillna(0).astype(int).values
            )

        out[f"{pc}:colour"] = out[pc].map({1: pc_color, 0: _ABSENT_COLOR})
        out.to_csv(dest_dir / "variants.csv", index=False)
        ok += 1

    print(f"  Phandango export: {ok} written, "
          f"{skipped_checkpoint} skipped (checkpoint), "
          f"{skipped_missing} skipped (source missing) — {total} total PCs")
