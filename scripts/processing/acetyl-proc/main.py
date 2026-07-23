"""
processing/acetyl-proc — acetyltransferase detection and annotation.

Pipeline:
  1. detect_sslbh_pc80()         — detect acetyltransferase PC80s from raw_hhsuite
                                    (ECOD SSLBH / PHROGs keywords / PFAM PF01757)
  2. flag_pyseer_with_sslbh()    — flag GWAS PCs in pyseer_hits_filtered (pvalue_corr ≤ 0.05)
                                    by tracing PC80 proteins → GWAS clustering levels
                                    → acetyltransferase/acetyl-gwas/pyseer_hits_sslbh_pvalcor005.tsv
  3. export_pc80_annotation()    — write PC80-level annotation tables and per-PC files
                                    → acetyltransferase/acetyl-annot-pc80/

Output root: cfg.output_dir / "processing" / "acetyltransferase"
Feeds figures/chapter3/.

TM-align screening against no-ecod GWAS predictors is retired — it depended
on the per-PC/protein-fold AF3 structure attachment in gwas-proc, which is
no longer populated (see gwas-proc/main.py). See other/unused/tmalign_screen.py
and other/unused/merge_tmalign.py.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))

from config import Config
from sslbh_detect import detect_sslbh_pc80
from cluster_map import flag_pyseer_with_sslbh
from pc80_export import export_pc80_annotation
from prophage_db import prepare_prophage_db
from build_kloci import build_acetylases_kloci
import pandas as pd

cfg = Config()

# ---------------------------------------------------------------------------
# Run flags
# ---------------------------------------------------------------------------
BUILD_ACETYLASES_KLOCI = True   # True → copy acetylases_kloci.tsv from input/cps/
RUN_PROPHAGE_BLAST = True   # True → BLASTP per PC80; skips if raw_blast.tsv exists

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
raw_hhsuite_tsv        = cfg.gwas_path / "2_PROPHAGES/raw_hhsuite.tsv"
pc2proteins_tsv        = cfg.gwas_path / "2_PROPHAGES/pcs2proteins.tsv"
prophages_metadata_tsv = cfg.gwas_path / "2_PROPHAGES/prophages_metadata.tsv"
pyseer_tsv             = cfg.gwas_path / "3_GWAS/3_PROCESSING/pyseer_hits_filtered.tsv"
mmseqs_dir             = cfg.gwas_path / "3_GWAS/1_INTERMEDIATE/2_MMSEQS"
prophage_faa_pattern   = str(cfg.gwas_path / "2_PROPHAGES/4_FASTA_CDS_AA/*.faa")

blast_db_dir     = cfg.output_dir / "other" / "prophage_blastp_db"

acetyl_dir         = cfg.output_dir / "processing" / "acetyltransferase"
cps_acetylases_dir = cfg.output_dir / "cps_acetylases"
s3_table_xlsx                = cfg.input_dir / "supplementary-thesis/supplementary-tables/S3_Table.xlsx"
enzymes_xlsx                = cfg.input_dir / "supplementary-thesis/supplementary-tables/S1_Table.xlsx"

# ---------------------------------------------------------------------------
# Step 0 — Copy acetylases_kloci.tsv from input
# ---------------------------------------------------------------------------
if BUILD_ACETYLASES_KLOCI:
    print("Step 0: Building acetylases_kloci.tsv from S3_Table.xlsx …")
    cps_acetylases_dir.mkdir(parents=True, exist_ok=True)
    kloci_df = build_acetylases_kloci(
        s3_table_xlsx = s3_table_xlsx,
        enzymes_xlsx  = enzymes_xlsx,
    )
    dest = cps_acetylases_dir / "acetylases_kloci.tsv"
    kloci_df.to_csv(dest, sep="\t", index=False)
    print(f"  acetylases_kloci.tsv — {len(kloci_df)} rows → {dest}")
else:
    print("[skip] Step 0: acetylases_kloci.tsv (BUILD_ACETYLASES_KLOCI = False)")

# ---------------------------------------------------------------------------
# Step 1 — Detect SSLBH / acetyltransferase PC80s
# ---------------------------------------------------------------------------
print("Step 1: Detecting SSLBH/acetyltransferase PC80s from raw_hhsuite …")
detected = detect_sslbh_pc80(raw_hhsuite_tsv, pc2proteins_tsv)
print(f"  {len(detected)} PC80s detected")

# Build flat protein → PC80 info map for downstream steps
protein_to_pc80: dict[str, dict] = {}
for _, row in detected.iterrows():
    info = {"pc80": row["pc80"], "detected_by": row["detected_by"]}
    for pid in row["protein_ids"]:
        protein_to_pc80[pid] = info
print(f"  {len(protein_to_pc80)} SSLBH protein IDs indexed")

# ---------------------------------------------------------------------------
# Step 2 — Flag pyseer GWAS hits with SSLBH evidence
# ---------------------------------------------------------------------------
print("\nStep 2: Flagging pyseer GWAS hits with SSLBH evidence …")
flag_pyseer_with_sslbh(
    pyseer_tsv      = pyseer_tsv,
    protein_to_pc80 = protein_to_pc80,
    mmseqs_dir      = mmseqs_dir,
    acetyl_dir      = acetyl_dir,
)

# Write acetylases_gwas.tsv — SSLBH-positive subset → cps_acetylases/
_flagged_tsv = acetyl_dir / "acetyl-gwas" / "pyseer_hits_sslbh_pvalcor005.tsv"
if _flagged_tsv.exists():
    _gwas_df = pd.read_csv(_flagged_tsv, sep="\t")
    _acetyl_hits = _gwas_df[_gwas_df["acetyl_hhsearch"] == True].copy()
    cps_acetylases_dir.mkdir(parents=True, exist_ok=True)
    _dst = cps_acetylases_dir / "acetylases_gwas.tsv"
    _acetyl_hits.to_csv(_dst, sep="\t", index=False)
    print(f"  acetylases_gwas.tsv — {len(_acetyl_hits)} rows ({_acetyl_hits['locus'].nunique()} loci) → {_dst}")

# ---------------------------------------------------------------------------
# Step 3 — Export PC80 annotation
# ---------------------------------------------------------------------------
print("\nStep 3: Exporting PC80 annotation …")
blast_db_path, seq_index = prepare_prophage_db(prophage_faa_pattern, blast_db_dir)
print(f"  {len(seq_index)} prophage proteins in index")

export_pc80_annotation(
    detected               = detected,
    raw_hhsuite_tsv        = raw_hhsuite_tsv,
    pc2proteins_tsv        = pc2proteins_tsv,
    prophages_metadata_tsv = prophages_metadata_tsv,
    seq_index              = seq_index,
    db_path                = blast_db_path,
    acetyl_dir             = acetyl_dir,
    run_blast              = RUN_PROPHAGE_BLAST,
)

print("\nDone.")
