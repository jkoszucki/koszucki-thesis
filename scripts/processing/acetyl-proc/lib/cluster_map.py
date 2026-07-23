"""
Flag pyseer GWAS hits with SSLBH / acetyltransferase evidence.

Approach:
  1. Load pyseer_hits_filtered.tsv; rename version→clustering_level; filter pvalue_corr ≤ 0.05
  2. For each clustering level: scan alignment FASTAs to find which GWAS PC contains
     each SSLBH protein (identified at PC80 level)
  3. Add three columns to the pyseer table:
       acetyl_hhsearch  — True / False
       acetyl_hhsearch_pc80   — comma-separated PC80 IDs (empty if not detected)
       acetyl_hhsearch_db     — comma-separated evidence sources (ECOD, PHROGs, PFAM)
  4. Write acetyltransferase/acetyl-gwas/pyseer_hits_sslbh_pvalcor005.tsv
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

PVALUE_CORR_MAX = 0.05


def flag_pyseer_with_sslbh(
    pyseer_tsv: Path,
    protein_to_pc80: dict[str, dict],
    mmseqs_dir: Path,
    acetyl_dir: Path,
) -> pd.DataFrame:
    """
    Load pyseer_hits_filtered.tsv, filter pvalue_corr ≤ 0.05, flag SSLBH-positive rows.

    Args:
        pyseer_tsv:      input_dir/.../pyseer_hits_filtered.tsv
        protein_to_pc80: {protein_id: {"pc80": str, "detected_by": str}}
                         flat mapping built from detected PC80s
        mmseqs_dir:      input_dir/.../2_MMSEQS; alignment FASTAs per clustering level
        acetyl_dir:      output_dir/acetyltransferase/

    Returns:
        Flagged DataFrame (all rows, pvalue_corr ≤ 0.05 filter applied).
    """
    # ---- Load & filter ----
    df = pd.read_csv(pyseer_tsv, sep="\t")
    df = df.rename(columns={"version": "clustering_level"})
    n_raw = len(df)

    if "mode" in df.columns:
        df = df[df["mode"] == "lasso"].copy()
    else:
        print("  [warn] mode column absent — no mode filter applied")
    if "pvalue_corr" in df.columns:
        df = df[df["pvalue_corr"] <= PVALUE_CORR_MAX].copy()
    else:
        print("  [warn] pvalue_corr column absent — no p-value filter applied")
    print(f"  pyseer: {n_raw:,} rows → {len(df):,} after mode==lasso, pvalue_corr ≤ {PVALUE_CORR_MAX}")

    # ---- Build GWAS PC → SSLBH info mapping ----
    target_proteins = set(protein_to_pc80)
    cl_pc_to_hits: dict[tuple[str, str], list[dict]] = {}

    for cl_name in df["clustering_level"].unique():
        protein_to_pc = _scan_fastas_for_proteins(
            mmseqs_dir / cl_name, target_proteins
        )
        print(f"  {cl_name}: {len(protein_to_pc)} / {len(target_proteins)} SSLBH proteins found in FASTAs")

        for pid, pc_name in protein_to_pc.items():
            key = (cl_name, pc_name)
            info = protein_to_pc80[pid]
            cl_pc_to_hits.setdefault(key, []).append(info)

    n_flagged = sum(
        1 for _, row in df.iterrows()
        if (row["clustering_level"], row["PC"]) in cl_pc_to_hits
    )
    print(f"  SSLBH-positive GWAS rows: {n_flagged} / {len(df)}")

    # ---- Add flag columns ----
    acetyl_hhsearch   = []
    acetyl_hhsearch_pc80    = []
    acetyl_hhsearch_db      = []

    for _, row in df.iterrows():
        hits = cl_pc_to_hits.get((row["clustering_level"], row["PC"]), [])
        if hits:
            acetyl_hhsearch.append(True)
            acetyl_hhsearch_pc80.append(",".join(sorted({h["pc80"] for h in hits})))
            acetyl_hhsearch_db.append(",".join(sorted({src for h in hits for src in h["detected_by"].split(",")})))
        else:
            acetyl_hhsearch.append(False)
            acetyl_hhsearch_pc80.append("")
            acetyl_hhsearch_db.append("")

    df["acetyl_hhsearch"] = acetyl_hhsearch
    df["acetyl_hhsearch_pc80"]  = acetyl_hhsearch_pc80
    df["acetyl_hhsearch_db"]    = acetyl_hhsearch_db

    # ---- Write ----
    out_dir = acetyl_dir / "acetyl-gwas"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "pyseer_hits_sslbh_pvalcor005.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"  → {out}  ({len(df):,} rows, {n_flagged} acetyl_hhsearch=True)")
    return df


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _scan_fastas_for_proteins(
    cl_dir: Path,
    target_proteins: set[str],
) -> dict[str, str]:
    """
    Scan alignment FASTAs in cl_dir/alignments/ to build {protein_id: pc_name}
    for any protein in target_proteins. Reads only FASTA header lines.
    Stops early once all targets are found.
    """
    alignments_dir = cl_dir / "alignments"
    if not alignments_dir.exists():
        print(f"  [warn] alignments dir not found: {alignments_dir}")
        return {}

    protein_to_pc: dict[str, str] = {}
    remaining = set(target_proteins)

    for fasta_path in sorted(alignments_dir.glob("*.fasta")):
        if not remaining:
            break
        pc_name = fasta_path.stem
        with open(fasta_path) as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                pid = line[1:].strip()
                if pid in remaining:
                    protein_to_pc[pid] = pc_name
                    remaining.discard(pid)

    return protein_to_pc
