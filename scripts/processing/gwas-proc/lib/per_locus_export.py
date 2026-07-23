"""
Step 2: Export per-PC data into the 4D output hierarchy.

For each hit in gwas_hits.tsv:
    gwas_root/{ecod_folder}/{locus}/{clustering_level}/{PC}/
        protein/
            pc.fasta               — MMseqs2 alignment FASTA (full cluster)
            sequence.fasta         — representative sequence (first from alignment)
            against-prophages/
                raw_blast.tsv      — BLASTP hits vs all prophage proteins; columns:
                                     qseqid, sseqid, pident, length, mismatch, gapopen,
                                     qstart, qend, sstart, send, evalue, bitscore,
                                     qcovhsp, qlen, slen, sseq


"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from prophage_db import prepare_prophage_db

import pandas as pd


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

def _pc_dir(gwas_root: Path, row) -> Path:
    """Return the 4D base directory for a gwas_hits row."""
    return gwas_root / row["ecod_folder"] / row["locus"] / row["clustering_level"] / row["PC"]


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



# ---------------------------------------------------------------------------
# Public: alignment FASTAs
# ---------------------------------------------------------------------------

def export_per_locus(
    gwas_hits_tsv: Path,
    mmseqs_dir: Path,
    gwas_root: Path,
) -> None:
    """
    Copy alignment FASTAs (→ pc.fasta) into the 4D output hierarchy.

    Args:
        gwas_hits_tsv: gwas_root / "gwas_hits.tsv"
        mmseqs_dir:    input_dir/.../2_MMSEQS
        gwas_root:     output root
    """
    hits = pd.read_csv(gwas_hits_tsv, sep="\t")

    # --- Alignments → pc.fasta ---
    print("  Alignments (pc.fasta) …")
    align_ok, align_missing = 0, []
    for _, row in hits.iterrows():
        src = mmseqs_dir / row["clustering_level"] / "alignments" / f"{row['PC']}.fasta"
        if not src.exists():
            align_missing.append(f"{row['locus']}/{row['clustering_level']}/{row['PC']}")
            continue
        dest_dir = _pc_dir(gwas_root, row) / "protein"
        dest_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest_dir / "pc.fasta")
        align_ok += 1
    print(f"    Copied {align_ok}/{len(hits)}")
    for m in align_missing[:10]:
        print(f"    missing: {m}")
    if len(align_missing) > 10:
        print(f"    … and {len(align_missing) - 10} more")


# ---------------------------------------------------------------------------
# Public: sequence FASTAs
# ---------------------------------------------------------------------------

def export_per_protein(
    gwas_hits_tsv: Path,
    mmseqs_dir: Path,
    gwas_root: Path,
) -> None:
    """
    Write sequence.fasta for each predictor.
    Representative sequence = first record from the MMseqs2 alignment FASTA.
    → gwas_root/{ecod_folder}/{locus}/{cl}/{PC}/protein/sequence.fasta
    """
    hits = pd.read_csv(gwas_hits_tsv, sep="\t")

    for _, row in hits.iterrows():
        locus, pc, cl = row["locus"], row["PC"], row["clustering_level"]
        fasta_path = mmseqs_dir / cl / "alignments" / f"{pc}.fasta"
        if not fasta_path.exists():
            print(f"  [warn] {locus}/{cl}/{pc}: alignment FASTA not found")
            continue
        protein_id, seq = _read_first_fasta(fasta_path)
        if not seq:
            print(f"  [warn] {locus}/{cl}/{pc}: empty sequence")
            continue

        protein_dir = _pc_dir(gwas_root, row) / "protein"
        protein_dir.mkdir(parents=True, exist_ok=True)
        header = protein_id if protein_id else f"{locus}_{pc}"
        with open(protein_dir / "sequence.fasta", "w") as fh:
            fh.write(f">{header}\n{seq}\n")


# ---------------------------------------------------------------------------
# Public: BLAST vs prophage proteins
# ---------------------------------------------------------------------------


def _run_blastp(
    query_id: str, query_seq: str, db_path: Path, tmp_dir: Path, evalue: float = 1e-5
) -> list[tuple[str, ...]]:
    """Run blastp; return list of raw hit rows (outfmt 6 fields)."""
    q_fasta = tmp_dir / f"{query_id}_q.faa"
    q_fasta.write_text(f">{query_id}\n{query_seq}\n")
    result = subprocess.run(
        ["blastp", "-query", str(q_fasta), "-db", str(db_path),
         "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp qlen slen",
         "-evalue", str(evalue), "-max_hsps", "1"],
        capture_output=True, text=True, check=True,
    )
    return [
        tuple(line.strip().split("\t"))
        for line in result.stdout.splitlines()
        if line.strip()
    ]


def _blast_one_protein(
    query_id: str,
    seq: str,
    db_path: Path,
    seq_index: dict[str, str],
    db_dir: Path,
    protein_dir: Path,
    run_blast: bool = True,
) -> None:
    """
    BLAST one protein against the prophage DB.
    Writes protein_dir/against-prophages/raw_blast.tsv with subject sequences appended.
    Skips if raw_blast.tsv already exists.
    """
    blast_cols = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp", "qlen", "slen",
    ]
    against_dir = protein_dir / "against-prophages"
    against_dir.mkdir(parents=True, exist_ok=True)

    raw_tsv = against_dir / "raw_blast.tsv"
    if raw_tsv.exists():
        print(f"  [{query_id}] raw_blast.tsv exists, skipping", flush=True)
        return
    if not run_blast:
        return

    print(f"  BLASTing {query_id} …", flush=True)
    blast_rows = _run_blastp(query_id, seq, db_path, db_dir)
    blast_df = pd.DataFrame(blast_rows, columns=blast_cols)
    blast_df["sseq"] = blast_df["sseqid"].map(seq_index)
    blast_df.to_csv(raw_tsv, sep="\t", index=False)
    print(f"    {len(blast_df)} hits → {raw_tsv.name}", flush=True)


def blast_repr_vs_prophage(
    gwas_hits_tsv: Path,
    mmseqs_dir: Path,
    prophage_faa_pattern: str,
    gwas_root: Path,
    blast_db_dir: Path,
    run_blast: bool = True,
    max_proteins: int | None = None,
) -> None:
    """
    BLAST all predictor proteins against prophage proteins.
    Representative sequence per predictor = first record from alignment FASTA.

    Args:
        max_proteins: if set, stop after this many proteins are BLASTed (for testing).

    Output: gwas_root/{ecod_folder}/{locus}/{cl}/{PC}/protein/against-prophages/raw_blast.tsv
    """
    print("  Prophage BLAST DB …")
    db_path, seq_index = prepare_prophage_db(prophage_faa_pattern, blast_db_dir)
    print(f"  {len(seq_index)} prophage proteins indexed")

    n_blasted = 0
    hits = pd.read_csv(gwas_hits_tsv, sep="\t")
    for _, row in hits.iterrows():
        if max_proteins is not None and n_blasted >= max_proteins:
            break
        locus, pc, cl = row["locus"], row["PC"], row["clustering_level"]
        fasta_path = mmseqs_dir / cl / "alignments" / f"{pc}.fasta"
        if not fasta_path.exists():
            print(f"  [{locus}/{cl}/{pc}] alignment FASTA not found — skipping BLAST")
            continue
        protein_id, rep_seq = _read_first_fasta(fasta_path)
        if not rep_seq:
            print(f"  [{locus}/{cl}/{pc}] empty sequence — skipping")
            continue
        _blast_one_protein(
            query_id=protein_id or f"{locus}_{pc}", seq=rep_seq,
            db_path=db_path, seq_index=seq_index, db_dir=blast_db_dir,
            protein_dir=_pc_dir(gwas_root, row) / "protein",
            run_blast=run_blast,
        )
        n_blasted += 1
