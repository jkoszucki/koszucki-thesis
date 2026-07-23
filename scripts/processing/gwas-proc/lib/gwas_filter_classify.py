"""
Step 1: Filter pyseer GWAS hits and classify by ECOD topology.

Quality filters applied to pyseer_hits_filtered.tsv:
    pvalue_corr ≤ 0.05   (skipped with warning if column absent)
    precision   ≥ 0.50
    F1_score    ≥ 0.50
    MCC         ≥ 0.50

ECOD classification uses reported_topology_PC from clusters_functions_best_all.tsv,
joined by (clustering_level, PC). Multiple rows per (PC, version) always carry the
same reported_topology_PC, so a single representative row is taken per pair.
    Topology → ecod_type:
        "sgnh-ecod"   — SGNH hydrolase
        "ssrbh-ecod"  — Pectin-lyase like (right-handed beta-helix; SSRBH fold)
        "other-ecod"  — any other annotated topology
        "no-ecod"     — no hit (reported_topology_PC absent, empty, or "no hit")

Output:
    gwas_root/gwas_hits.tsv — one row per (locus, ecod_type, clustering_level, PC)
    Includes ecod_folder column: e.g. "sgnh-ecod-reported-topology"
"""

from __future__ import annotations

import pandas as pd
from pathlib import Path

PVALUE_CORR_MAX    = 0.05
GWAS_PRECISION_MIN = 0.50
GWAS_F1_MIN        = 0.50
GWAS_MCC_MIN       = 0.50

_ECOD_FOLDER = {
    "sgnh-ecod":  "sgnh-ecod-reported-topology",
    "ssrbh-ecod": "ssrbh-ecod-reported-topology",
    "other-ecod": "other-ecod-reported-topology",
    "no-ecod":    "no-ecod-reported-topology",
}


def _classify_ecod_type(topology: str | None) -> str:
    """Map reported_topology_PC to one of four ecod_type labels."""
    if not topology or pd.isna(topology) or str(topology).lower() == "no hit":
        return "no-ecod"
    t = topology.lower()
    if "sgnh" in t:
        return "sgnh-ecod"
    if "pectin" in t:
        return "ssrbh-ecod"
    return "other-ecod"


def _load_topology_map(functions_path: Path) -> pd.DataFrame:
    """
    Load clusters_functions_best_all.tsv and return one row per (clustering_level, PC)
    with column reported_topology_PC.

    Multiple rows per (PC, version) always carry the same reported_topology_PC,
    so duplicates are dropped after selecting the key columns.
    """
    df = pd.read_csv(functions_path, sep="\t", usecols=["PC", "version", "reported_topology_PC"])
    df = df.rename(columns={"version": "clustering_level"})
    df = df.drop_duplicates(subset=["clustering_level", "PC"])
    print(f"  {len(df):,} (clustering_level, PC) topology entries loaded.")
    return df


def filter_and_classify(
    table_path: Path,
    functions_path: Path,
    gwas_root: Path,
) -> pd.DataFrame:
    """
    Filter pyseer hits by quality thresholds and classify by ECOD topology.

    reported_topology_PC is joined from clusters_functions_best_all.tsv
    by (clustering_level, PC). Hits with no matching entry → no-ecod.

    Writes gwas_root/gwas_hits.tsv with all original pyseer columns plus:
        clustering_level, reported_topology_PC, ecod_type, ecod_folder

    Args:
        table_path:     pyseer_hits_filtered.tsv
        functions_path: clusters_functions_best_all.tsv
        gwas_root:      output root; gwas_hits.tsv written here

    Returns:
        gwas_hits DataFrame
    """
    gwas_root.mkdir(parents=True, exist_ok=True)

    print(f"  Loading {table_path.name} …")
    df = pd.read_csv(table_path, sep="\t")
    print(f"  {len(df):,} rows loaded.")
    df = df.rename(columns={"version": "clustering_level"})

    flt = (
        (df["mode"]      == "lasso") &
        (df["precision"] >= GWAS_PRECISION_MIN) &
        (df["F1_score"]  >= GWAS_F1_MIN) &
        (df["MCC"]       >= GWAS_MCC_MIN)
    )
    if "pvalue_corr" in df.columns:
        flt &= (df["pvalue_corr"] <= PVALUE_CORR_MAX)
    else:
        print("  [warn] pvalue_corr column absent — skipping pvalue filter.")
    df = df[flt].copy()
    print(f"  {len(df):,} rows after quality filters.")

    topology_map = _load_topology_map(functions_path)
    df = df.merge(topology_map, on=["clustering_level", "PC"], how="left")
    df["reported_topology_PC"] = df["reported_topology_PC"].fillna("")

    df["ecod_type"]   = df["reported_topology_PC"].apply(_classify_ecod_type)
    df["ecod_folder"] = df["ecod_type"].map(_ECOD_FOLDER)

    print("  ECOD-type breakdown:")
    for ecod_type, grp in df.groupby("ecod_type"):
        print(f"    {ecod_type}: {len(grp):,} rows, {grp['locus'].nunique()} unique loci")

    out = gwas_root / "gwas_hits.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"  Saved → {out}  ({len(df):,} rows)")
    return df
