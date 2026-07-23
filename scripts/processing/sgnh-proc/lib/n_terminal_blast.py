"""
BLAST the shared N-terminal fragment against KPSC prophage proteins.

The query is extracted automatically from the PC0262 alignment FASTA
(first sequence, first 120 amino acids — the conserved N-terminal anchor).

Output: n_terminal_dir/raw_blast.tsv
        n_terminal_dir/query.fasta
Format: qseqid, sseqid, pident, length, mismatch, gapopen,
        qstart, qend, sstart, send, evalue, bitscore, qcovhsp, qlen, slen, sseq
Checkpoint: skips if raw_blast.tsv already exists.
"""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import pandas as pd

N_TERMINAL_LENGTH = 120

BLAST_COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp", "qlen", "slen",
]


def _extract_query(pc_fasta: Path, length: int = N_TERMINAL_LENGTH) -> tuple[str, str]:
    """
    Return (query_id, n_terminal_seq) from the first record of pc_fasta.

    query_id is built from the protein ID and pc_fasta path components:
      {protein_id}_{locus}_{pc}_{clustering_level}_1-{length}
    e.g. KPN_B10_PHAGE037_M_PROTEIN_49_KL24_PC0262_PCI80C50_1-120
    """
    from Bio import SeqIO

    record = next(SeqIO.parse(pc_fasta, "fasta"))
    seq = str(record.seq).replace("-", "")[:length]

    # Derive locus / clustering_level / pc from path:
    # …/KL24/lasso/alignments/PCI80C50/PC0262.fasta
    parts  = pc_fasta.parts
    pc     = pc_fasta.stem                         # PC0262
    cl     = pc_fasta.parent.name                  # PCI80C50
    locus  = parts[parts.index("1_PER_LOCUS") + 1] # KL24

    query_id = f"{record.id}_{locus}_{pc}_{cl}_1-{length}"
    return query_id, seq


def blast_n_terminal(
    prophage_faa_pattern: str,
    blast_db_dir: Path,
    n_terminal_dir: Path,
    pc_fasta: Path,
) -> None:
    """
    BLAST the N-terminal anchor against the KPSC prophage DB.

    The query sequence is extracted from the first record of pc_fasta
    (residues 1–120 after stripping alignment gaps).

    Args:
        prophage_faa_pattern: glob for prophage FASTAs (passed to prepare_prophage_db)
        blast_db_dir:         directory where the prophage BLAST DB lives
        n_terminal_dir:       output directory; receives query.fasta + raw_blast.tsv
        pc_fasta:             PC alignment FASTA to extract query from
    """
    from prophage_db import prepare_prophage_db

    n_terminal_dir.mkdir(parents=True, exist_ok=True)

    query_id, query_seq = _extract_query(pc_fasta)
    query_fasta = n_terminal_dir / "query.fasta"
    query_fasta.write_text(f">{query_id}\n{query_seq}\n")
    print(f"  query: {query_id} ({len(query_seq)} aa)")

    raw_tsv = n_terminal_dir / "raw_blast.tsv"
    if raw_tsv.exists():
        print(f"  raw_blast.tsv exists — skipping")
        return

    db_path, seq_index = prepare_prophage_db(prophage_faa_pattern, blast_db_dir)
    print(f"  {len(seq_index)} prophage proteins in DB")

    with tempfile.TemporaryDirectory() as tmp_str:
        q_fasta = Path(tmp_str) / "query.faa"
        q_fasta.write_text(f">{query_id}\n{query_seq}\n")

        result = subprocess.run(
            [
                "blastp",
                "-query",  str(q_fasta),
                "-db",     str(db_path),
                "-outfmt", "6 qseqid sseqid pident length mismatch gapopen "
                           "qstart qend sstart send evalue bitscore qcovhsp qlen slen",
                "-evalue", "1e-5",
                "-max_hsps", "1",
            ],
            capture_output=True, text=True, check=True,
        )

    rows = [line.split("\t") for line in result.stdout.strip().splitlines() if line]
    blast_df = pd.DataFrame(rows, columns=BLAST_COLS)

    for col in ("pident", "evalue", "bitscore", "qcovhsp"):
        blast_df[col] = pd.to_numeric(blast_df[col], errors="coerce")
    for col in ("length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "qlen", "slen"):
        blast_df[col] = pd.to_numeric(blast_df[col], errors="coerce").astype("Int64")

    blast_df["sseq"] = blast_df["sseqid"].map(seq_index)
    blast_df.to_csv(raw_tsv, sep="\t", index=False)
    print(f"  {len(blast_df)} hits → {raw_tsv}")
