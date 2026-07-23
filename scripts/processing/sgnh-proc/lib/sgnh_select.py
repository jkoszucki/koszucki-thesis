"""
Select SGNH hydrolase hits from gwas_hits.tsv and organise outputs.

Steps:
  1. Filter gwas_hits.tsv for ecod_type == "sgnh-ecod"
     → sgnh-hydrolases/gwas_sgnh_hits.tsv   (all SGNH hits)
  2. Select best predictor per K-locus (highest F1_score); map representative
     protein ID and sequence from MMseqs2 alignment FASTA
     → sgnh-hydrolases/gwas_sgnh_best.tsv
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# Representative sequence helpers (ported from gwas_sgnh_selection.py)
# ---------------------------------------------------------------------------

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
        print(f"  {missing} alignment FASTA(s) missing — representative_* columns will be empty for those rows.")
    best = best.copy()
    best["representative_protein_id"] = pids
    best["representative_sequence"]   = seqs
    return best


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def select_sgnh(
    gwas_hits_tsv: Path,
    mmseqs_dir: Path,
    sgnh_dir: Path,
    precision_threshold: float = 0.80,
) -> None:
    """
    Filter SGNH hits, write all-hits and best-per-locus TSVs.

    Args:
        gwas_hits_tsv:       path to gwas_hits.tsv (output of gwas-proc)
        mmseqs_dir:          input MMseqs2 directory (for representative sequences)
        sgnh_dir:            output_dir/sgnh-hydrolases (all outputs written here)
        precision_threshold: minimum precision for best-predictor selection (default 0.80)
    """
    sgnh_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Step 1 — all SGNH hits
    # ------------------------------------------------------------------
    df   = pd.read_csv(gwas_hits_tsv, sep="\t")
    sgnh = df[df["ecod_type"] == "sgnh-ecod"].copy()
    print(f"  {len(sgnh):,} sgnh-ecod rows across {sgnh['locus'].nunique()} loci")

    all_out = sgnh_dir / "gwas_sgnh_hits.tsv"
    sgnh.to_csv(all_out, sep="\t", index=False)
    print(f"  Saved → {all_out.name}")

    if sgnh.empty:
        print("  No SGNH hits — skipping best-predictor selection.")
        return

    # ------------------------------------------------------------------
    # Step 2 — best predictor per K-locus (precision >= threshold)
    # ------------------------------------------------------------------
    eligible = sgnh[sgnh["precision"] >= precision_threshold]
    print(f"  {len(eligible):,} rows after precision ≥ {precision_threshold} filter "
          f"({eligible['locus'].nunique()} loci)")
    best = (
        eligible.loc[eligible.groupby("locus")["F1_score"].idxmax()]
                .sort_values("locus")
                .reset_index(drop=True)
    )
    print(f"  {len(best)} best SGNH predictors selected (one per K-locus)")
    best = _map_representative_sequences(best, mmseqs_dir)

    best_out = sgnh_dir / "gwas_sgnh_best.tsv"
    best.to_csv(best_out, sep="\t", index=False)
    print(f"  Saved → {best_out.name}")
