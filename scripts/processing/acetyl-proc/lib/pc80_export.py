"""
Export PC80 acetyltransferase annotation.

Writes to acetyltransferase/acetyl-annot-pc80/:

    acetyl-hhsearch-pc80.tsv — raw_hhsuite.tsv rows for detected PC80s (all triggering hits)
    acetyl-pc80.tsv       — one row per protein: pc80, detected_by,
                            proteinID, prophageID, genomeID, sequence

    {pc80}/
        sequence.fasta         — representative protein (first in pc2proteins order)
        against-prophages/
            raw_blast.tsv      — BLASTP vs prophage DB (checkpoint: skip if exists)
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pandas as pd

_PROPHAGE_ID_PATTERN = re.compile(r"_PROTEIN_\d+$", re.IGNORECASE)

BLAST_COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "qcovhsp", "qlen", "slen",
]


def export_pc80_annotation(
    detected: pd.DataFrame,
    raw_hhsuite_tsv: Path,
    pc2proteins_tsv: Path,
    prophages_metadata_tsv: Path,
    seq_index: dict[str, str],
    db_path: Path,
    acetyl_dir: Path,
    run_blast: bool = True,
) -> None:
    """
    Write acetyl-annot-pc80/ contents for all detected PC80s.

    Args:
        detected:              output of detect_sslbh_pc80() — pc80, detected_by, protein_ids
        raw_hhsuite_tsv:       raw HHsearch results (for acetyl-hhsearch-pc80.tsv)
        pc2proteins_tsv:       PC → proteinID mapping
        prophages_metadata_tsv: prophageID → genomeID
        seq_index:             protein_id → sequence (from prepare_prophage_db)
        db_path:               prophage BLAST DB path
        acetyl_dir:            output_dir/acetyltransferase/
        run_blast:             if False, skip BLASTP
    """
    out_root = acetyl_dir / "acetyl-annot-pc80"
    out_root.mkdir(parents=True, exist_ok=True)

    pc80_set = set(detected["pc80"])

    # ---- acetyl-hhsearch-pc80.tsv ----
    _write_hhsearch_table(raw_hhsuite_tsv, pc80_set, out_root)

    # ---- Build pc80 annotation table ----
    prophage_to_genome = _load_prophage_to_genome(prophages_metadata_tsv)
    pc2prot_df = _load_pc2proteins_full(pc2proteins_tsv)

    annotation_rows = []
    for _, det_row in detected.iterrows():
        pc80       = det_row["pc80"]
        detected_by = det_row["detected_by"]
        protein_ids = det_row["protein_ids"]

        for pid in protein_ids:
            phage_id  = _PROPHAGE_ID_PATTERN.sub("", pid)
            genome_id = prophage_to_genome.get(phage_id, "")
            seq       = seq_index.get(pid, "")
            annotation_rows.append({
                "pc80":        pc80,
                "detected_by": detected_by,
                "proteinID":   pid,
                "prophageID":  phage_id,
                "genomeID":    genome_id,
                "sequence":    seq,
            })

    annot_df = pd.DataFrame(annotation_rows)
    annot_tsv = out_root / "acetyl-pc80.tsv"
    annot_df.to_csv(annot_tsv, sep="\t", index=False)
    print(f"  → {annot_tsv}  ({len(annot_df)} protein rows across {len(pc80_set)} PC80s)")

    # ---- Per-PC80 files ----
    for _, det_row in detected.iterrows():
        pc80       = det_row["pc80"]
        protein_ids = det_row["protein_ids"]

        # Representative = first protein with a resolvable sequence
        rep_id, rep_seq = _pick_representative(protein_ids, seq_index)
        if not rep_id:
            print(f"  [warn] {pc80}: no sequence found for any protein — skipping per-PC export")
            continue

        pc_dir = out_root / pc80
        pc_dir.mkdir(exist_ok=True)

        # sequence.fasta
        with open(pc_dir / "sequence.fasta", "w") as fh:
            fh.write(f">{rep_id}\n{rep_seq}\n")

        # BLASTP
        _blast_one(
            query_id=rep_id,
            seq=rep_seq,
            db_path=db_path,
            seq_index=seq_index,
            tmp_dir=pc_dir,
            protein_dir=pc_dir,
            run_blast=run_blast,
        )

    print(f"  Per-PC export complete → {out_root}")


# ---------------------------------------------------------------------------
# Helpers — data loading
# ---------------------------------------------------------------------------

def _write_hhsearch_table(
    raw_hhsuite_tsv: Path,
    pc80_set: set[str],
    out_root: Path,
) -> None:
    df = pd.read_csv(raw_hhsuite_tsv, sep="\t")
    df.columns = [c.lower() for c in df.columns]
    query_col = "query" if "query" in df.columns else df.columns[0]
    filtered = df[df[query_col].isin(pc80_set)].copy()
    out = out_root / "acetyl-hhsearch-pc80.tsv"
    filtered.to_csv(out, sep="\t", index=False)
    print(f"  → {out}  ({len(filtered)} HHsearch rows for {len(pc80_set)} PC80s)")


def _load_prophage_to_genome(prophages_metadata_tsv: Path) -> dict[str, str]:
    df = pd.read_csv(prophages_metadata_tsv, sep="\t")
    df.columns = [c.lower() for c in df.columns]
    return dict(zip(df["prophageid"], df["genomeid"]))


def _load_pc2proteins_full(pc2proteins_tsv: Path) -> pd.DataFrame:
    df = pd.read_csv(pc2proteins_tsv, sep="\t", dtype=str)
    df.columns = [c.lower() for c in df.columns]
    return df


def _pick_representative(protein_ids: list[str], seq_index: dict[str, str]) -> tuple[str, str]:
    """Return (protein_id, sequence) for the first protein with a known sequence."""
    for pid in protein_ids:
        seq = seq_index.get(pid, "")
        if seq:
            return pid, seq
    return "", ""


# ---------------------------------------------------------------------------
# Helpers — BLASTP
# ---------------------------------------------------------------------------

def _blast_one(
    query_id: str,
    seq: str,
    db_path: Path,
    seq_index: dict[str, str],
    tmp_dir: Path,
    protein_dir: Path,
    run_blast: bool = True,
    evalue: float = 1e-5,
) -> None:
    against_dir = protein_dir / "against-prophages"
    against_dir.mkdir(parents=True, exist_ok=True)
    raw_tsv = against_dir / "raw_blast.tsv"

    if raw_tsv.exists():
        print(f"  [{query_id}] raw_blast.tsv exists, skipping", flush=True)
        return
    if not run_blast:
        return

    print(f"  BLASTing {query_id} …", flush=True)
    q_fasta = tmp_dir / f"_query_{query_id}.faa"
    q_fasta.write_text(f">{query_id}\n{seq}\n")

    result = subprocess.run(
        [
            "blastp",
            "-query", str(q_fasta),
            "-db", str(db_path),
            "-outfmt",
            "6 qseqid sseqid pident length mismatch gapopen "
            "qstart qend sstart send evalue bitscore qcovhsp qlen slen",
            "-evalue", str(evalue),
            "-max_hsps", "1",
        ],
        capture_output=True, text=True, check=True,
    )
    q_fasta.unlink(missing_ok=True)

    rows = [line.strip().split("\t") for line in result.stdout.splitlines() if line.strip()]
    blast_df = pd.DataFrame(rows, columns=BLAST_COLS)
    blast_df["sseq"] = blast_df["sseqid"].map(seq_index)
    blast_df.to_csv(raw_tsv, sep="\t", index=False)
    print(f"    {len(blast_df)} hits → {raw_tsv.name}", flush=True)
