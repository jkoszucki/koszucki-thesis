"""
scripts/helpers/prophage_db.py

Prepare the central prophage protein FASTA and BLASTP database.

Outputs (output_dir/other/prophage_blastp_db/):
    prophage_proteins.faa   — all prophage proteins with curated headers
    prophage_db.*           — blastp DB files (makeblastdb)

Checkpoint: if prophage_db.* files already exist the script exits immediately
without overwriting anything. Delete the DB files to force a rebuild.

Usage:
    conda run -n jkoszucki python prophage_db.py
"""

from __future__ import annotations

import glob as _glob
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import Config

from Bio import SeqIO

cfg = Config()

PROPHAGE_FAA_PATTERN = str(cfg.gwas_path / "2_PROPHAGES/4_FASTA_CDS_AA/*.faa")
DB_DIR               = cfg.output_dir / "other" / "prophage_blastp_db"


# ---------------------------------------------------------------------------
# Prophage DB preparation
# ---------------------------------------------------------------------------

def prepare_prophage_db(
    prophage_faa_pattern: str,
    db_dir: Path,
) -> tuple[Path, dict[str, str]]:
    """
    1. Load all prophage FASTAs; extract protein name from ' ||| '-delimited header.
       Format: >CDS_1 ||| KPN_B1_PHAGE002_M_PROTEIN_1 ||| ...
    2. Write prophage_proteins.faa with clean >{protein_name} headers.
    3. Build blastp DB from prophage_proteins.faa.

    Checkpoint: if prophage_db.* files already exist, skip steps 2–3.
    Returns (db_path, seq_index {protein_name: sequence}).
    """
    db_dir.mkdir(parents=True, exist_ok=True)
    db_path  = db_dir / "prophage_db"
    faa_out  = db_dir / "prophage_proteins.faa"

    faa_files = sorted(_glob.glob(prophage_faa_pattern))
    if not faa_files:
        raise FileNotFoundError(f"No prophage FASTAs matched: {prophage_faa_pattern}")

    # Load all records with curated headers
    records: list[tuple[str, str]] = []
    for faa in faa_files:
        for rec in SeqIO.parse(faa, "fasta"):
            parts = rec.description.split(" ||| ")
            label = parts[1] if len(parts) >= 2 else rec.id
            records.append((label, str(rec.seq)))

    seq_index = dict(records)

    # Checkpoint — skip if DB already exists
    if list(db_dir.glob("prophage_db.*")):
        print(f"  Prophage DB already exists ({len(seq_index)} proteins) — skipping")
        return db_path, seq_index

    # Step 1: write prophage_proteins.faa
    print(f"  Writing {faa_out.name} ({len(records)} sequences) …")
    with open(faa_out, "w") as out:
        for label, seq in records:
            out.write(f">{label}\n{seq}\n")

    # Step 2: build BLAST DB
    print("  Building BLAST DB …")
    subprocess.run(
        ["makeblastdb", "-in", str(faa_out), "-dbtype", "prot", "-out", str(db_path)],
        check=True,
    )
    print(f"  Done → {db_dir}")
    return db_path, seq_index


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    db_path, seq_index = prepare_prophage_db(PROPHAGE_FAA_PATTERN, DB_DIR)
    print(f"\n  DB: {db_path}")
    print(f"  Proteins: {len(seq_index)}")
