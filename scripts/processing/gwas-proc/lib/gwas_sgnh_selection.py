"""
Select best GWAS predictors per K-locus.

Best SGNH predictor per K-locus:
    Reads gwas_hits.tsv; filters ecod_type == "sgnh-ecod".
    Selects highest F1_score per locus across all clustering levels.
    Maps representative_protein_id and representative_sequence
    (first record from MMseqs2 alignment FASTA) for each selected predictor.
    → gwas_root/gwas_sgnh_hits.tsv
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def _read_first_fasta(fasta_path: Path) -> tuple[str, str]:
    """Return (protein_id, sequence) for the first record in a FASTA file."""
    protein_id, seq_lines, found = "", [], False
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if found:
                    break
                protein_id = line[1:]
                found = True
            elif found:
                seq_lines.append(line)
    return protein_id, "".join(seq_lines)


def _map_representative_sequences(best: pd.DataFrame, mmseqs_dir: Path) -> pd.DataFrame:
    """Add representative_protein_id and representative_sequence columns."""
    pids, seqs = [], []
    missing = 0
    for _, row in best.iterrows():
        fasta = mmseqs_dir / row["clustering_level"] / "alignments" / f"{row['PC']}.fasta"
        if not fasta.exists():
            print(f"  [warn] alignment FASTA not found: {fasta}")
            pids.append("")
            seqs.append("")
            missing += 1
        else:
            pid, seq = _read_first_fasta(fasta)
            pids.append(pid)
            seqs.append(seq)
    if missing:
        print(f"  {missing} alignment FASTAs missing.")
    best = best.copy()
    best["representative_protein_id"] = pids
    best["representative_sequence"]   = seqs
    return best


def select_sgnh_predictors(
    gwas_hits_tsv: Path,
    mmseqs_dir: Path,
    gwas_root: Path,
) -> pd.DataFrame:
    """
    Best SGNH predictor per K-locus: highest F1_score from sgnh-ecod hits.

    Args:
        gwas_hits_tsv: gwas_root / "gwas_hits.tsv"
        mmseqs_dir:    input_dir/.../2_MMSEQS (for representative sequences)
        gwas_root:     output root; gwas_sgnh_hits.tsv written here

    Returns:
        DataFrame with one row per K-locus (best SGNH predictor).
    """
    df = pd.read_csv(gwas_hits_tsv, sep="\t")
    sgnh = df[df["ecod_type"] == "sgnh-ecod"].copy()
    print(f"  {len(sgnh):,} sgnh-ecod hits across {sgnh['locus'].nunique()} loci")

    out = gwas_root / "gwas_sgnh_hits.tsv"
    if sgnh.empty:
        print("  No SGNH hits found.")
        pd.DataFrame().to_csv(out, sep="\t", index=False)
        return pd.DataFrame()

    best = (sgnh.loc[sgnh.groupby("locus")["F1_score"].idxmax()]
               .sort_values("locus")
               .reset_index(drop=True))
    print(f"  {len(best)} best SGNH predictors selected.")

    best = _map_representative_sequences(best, mmseqs_dir)
    best.to_csv(out, sep="\t", index=False)
    print(f"  Saved → {out}")
    return best


