"""
TM-align all-vs-all structural comparison for CIF files.

Requires the `tmtools` conda environment:
    conda run -n tmtools python script.py

Usage
-----
    from tmalign import tm_align_all_vs_all

    paths = list(Path("structures/").glob("*.cif"))
    df = tm_align_all_vs_all(paths)
    df.to_csv("tmscore.tsv", sep="\\t", index=False)

Output columns
--------------
    id1, id2            — filename stems
    tm_norm_1           — TM-score normalised by length of structure 1
    tm_norm_2           — TM-score normalised by length of structure 2
    tm_max              — max(tm_norm_1, tm_norm_2)
    rmsd                — RMSD of aligned residues
"""

from __future__ import annotations

from itertools import combinations
from pathlib import Path

import pandas as pd
from Bio.PDB import MMCIFParser
from tmtools import tm_align
from tmtools.io import get_residue_data


def _load_chain(cif_path: Path) -> tuple[str, object]:
    """Return (stem_id, first_chain) for a CIF file."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(cif_path.stem, str(cif_path))
    chain = next(structure.get_chains())
    return cif_path.stem, chain


def tm_align_all_vs_all(cif_paths: list[Path]) -> pd.DataFrame:
    """
    Run TM-align for every pair of CIF structures and return results as a DataFrame.

    Each unordered pair appears once. TM-scores for both normalisation directions
    are reported.

    Args:
        cif_paths: list of paths to CIF structure files.

    Returns:
        DataFrame with columns: id1, id2, tm_norm_1, tm_norm_2, rmsd.
    """
    # Load all structures once
    print(f"  Loading {len(cif_paths)} structures …")
    structures: list[tuple[str, object]] = []
    for path in cif_paths:
        stem, chain = _load_chain(path)
        structures.append((stem, chain))
        print(f"    {stem}")

    rows = []
    pairs = list(combinations(range(len(structures)), 2))
    print(f"  Running TM-align on {len(pairs)} pairs …")

    for i, j in pairs:
        id1, chain1 = structures[i]
        id2, chain2 = structures[j]

        coords1, seq1 = get_residue_data(chain1)
        coords2, seq2 = get_residue_data(chain2)

        result = tm_align(coords1, coords2, seq1, seq2)
        tm1 = round(result.tm_norm_chain1, 4)
        tm2 = round(result.tm_norm_chain2, 4)
        rows.append({
            "id1":       id1,
            "id2":       id2,
            "tm_norm_1": tm1,
            "tm_norm_2": tm2,
            "tm_max":    max(tm1, tm2),
            "rmsd":      round(result.rmsd, 4),
        })
        print(f"    {id1} vs {id2}: TM1={result.tm_norm_chain1:.4f}  TM2={result.tm_norm_chain2:.4f}  RMSD={result.rmsd:.4f}")

    return pd.DataFrame(rows, columns=["id1", "id2", "tm_norm_1", "tm_norm_2", "tm_max", "rmsd"])
