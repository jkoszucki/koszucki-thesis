"""
processing/gwas-proc — shared data preparation for chapters 2 and 3.

Pipeline:
  1. filter_and_classify()           → gwas_hits.tsv
  2. per-PC file export
       2a. export_per_locus()            → {ecod_folder}/{locus}/{cl}/{PC}/protein/pc.fasta
       2b. export_per_protein()          → protein/sequence.fasta
       2d. prepare_prophage_db()         → other/prophage_blastp_db/prophage_proteins.faa, prophage_db.*
       2e. blast_repr_vs_prophage()      → protein/against-prophages/raw_blast.tsv

SGNH predictor selection moved to processing/sgnh-proc/.

AF3 structures are no longer attached per-PC here — that per-fold-group AF3
attachment (via manual-outputs/manual-upload) is retired; the live AF3
pipeline is other/alphafold3/ (see scripts/helpers/build_af3_index.py).

Output root: cfg.output_dir / "processing" / "gwas-data"
Feeds figures/chapter2/ (sgnh-ecod entries, via sgnh-proc symlinks)
and  figures/chapter3/ (no-ecod, ssrbh-ecod entries).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))

from config import Config
from gwas_filter_classify import filter_and_classify
from per_locus_export import (
    export_per_locus,
    export_per_protein,
    blast_repr_vs_prophage,
)
from phandango_export import export_phandango
from phandango_render import render_all_phandango
from prophage_db import prepare_prophage_db
from depolymerases_tables import build_depolymerases_tables
from best_predictors import build_best_predictors

cfg = Config()

# ---------------------------------------------------------------------------
# Run flags — toggle steps without modifying config.yml
# When ON: step runs with per-PC checkpointing — skips any PC whose output
# file already exists (raw_blast.tsv / grr.tsv). Delete files to recompute.
# When OFF: step is skipped entirely.
# ---------------------------------------------------------------------------
BUILD_DEPOLYMERASES_TABLES = True   # True → write rbp_depolymerases/ TSVs from PLOS S1/S3/S4 tables
BUILD_BEST_PREDICTORS      = True   # True → write rbp_best_predictors/best_predictors_gwas.tsv (precision≥0.8)
RUN_PROPHAGE_BLAST       = False    # True → blastp per PC; skips if raw_blast.tsv exists
RUN_PHANDANGO_EXPORT     = False    # True → export subtree + variants per PC; skips if variants.csv exists
RUN_PHANDANGO_RENDER     = False    # True → render PNG per PC; skips if PNG exists

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
table_path     = cfg.gwas_path / "3_GWAS/3_PROCESSING/pyseer_hits_filtered.tsv"
functions_path = cfg.gwas_path / "3_GWAS/1_INTERMEDIATE/3_FUNCTIONS/clusters_functions_best_all.tsv"
mmseqs_dir     = cfg.gwas_path / "3_GWAS/1_INTERMEDIATE/2_MMSEQS"
gwas_root      = cfg.output_dir / "processing" / "gwas-data"

blast_db_dir        = cfg.output_dir / "other" / "prophage_blastp_db"

plos_tables_dir       = cfg.input_dir / "supplementary-otwinowska" / "supplementary-tables"
rbp_depolymerases_dir = cfg.output_dir / "rbp_depolymerases"
rbp_best_dir          = cfg.output_dir / "rbp_best_predictors"

# ---------------------------------------------------------------------------
# Step 0 — build rbp_depolymerases/ from PLOS supplement tables
# ---------------------------------------------------------------------------
if BUILD_DEPOLYMERASES_TABLES:
    print("Step 0: Building rbp_depolymerases/ tables from PLOS supplement tables …")
    build_depolymerases_tables(plos_tables_dir, rbp_depolymerases_dir)
else:
    print("\n[skip] Step 0: depolymerases tables (BUILD_DEPOLYMERASES_TABLES = False)")

# ---------------------------------------------------------------------------
# Step 1 — filter, classify, write gwas_hits.tsv
# ---------------------------------------------------------------------------
print("Step 1: Filtering and classifying GWAS hits …")
filter_and_classify(table_path, functions_path, gwas_root)

gwas_hits_tsv = gwas_root / "gwas_hits.tsv"

if BUILD_BEST_PREDICTORS:
    print("\nStep 1b: Building best_predictors_gwas.tsv (precision≥0.8, F1≥0.5, MCC≥0.5) …")
    build_best_predictors(gwas_hits_tsv, rbp_best_dir)
else:
    print("\n[skip] Step 1b: best predictors (BUILD_BEST_PREDICTORS = False)")

# ---------------------------------------------------------------------------
# Step 2 — per-PC file export
# ---------------------------------------------------------------------------
print("\nStep 2a: Exporting per-PC alignment FASTAs …")
export_per_locus(
    gwas_hits_tsv = gwas_hits_tsv,
    mmseqs_dir    = mmseqs_dir,
    gwas_root     = gwas_root,
)

print("\nStep 2b: Writing representative sequence FASTAs …")
export_per_protein(
    gwas_hits_tsv = gwas_hits_tsv,
    mmseqs_dir    = mmseqs_dir,
    gwas_root     = gwas_root,
)

print("\nStep 2d: Building prophage BLAST DB …")
prepare_prophage_db(
    prophage_faa_pattern = str(cfg.gwas_path / "2_PROPHAGES/4_FASTA_CDS_AA/*.faa"),
    db_dir               = blast_db_dir,
)

print("\nStep 2e: BLASTP vs prophage proteins …")
blast_repr_vs_prophage(
    gwas_hits_tsv        = gwas_hits_tsv,
    mmseqs_dir           = mmseqs_dir,
    prophage_faa_pattern = str(cfg.gwas_path / "2_PROPHAGES/4_FASTA_CDS_AA/*.faa"),
    gwas_root            = gwas_root,
    blast_db_dir         = blast_db_dir,
    run_blast            = RUN_PROPHAGE_BLAST,
)

if RUN_PHANDANGO_EXPORT:
    print("\nStep 2h: Exporting phandango files …")
    export_phandango(
        gwas_hits_tsv      = gwas_hits_tsv,
        phandango_src_root = cfg.gwas_path / "3_GWAS/4_ANALYZE/1_PER_LOCUS",
        mmseqs_dir         = mmseqs_dir,
        gwas_data_dir      = gwas_root,
        style              = cfg.style,
    )
else:
    print("\n[skip] Step 2h: phandango export (RUN_PHANDANGO_EXPORT = False)")

if RUN_PHANDANGO_RENDER:
    print("\nStep 2i: Rendering phandango PNGs …")
    render_all_phandango(
        gwas_hits_tsv = gwas_hits_tsv,
        gwas_data_dir = gwas_root,
        meta_tsv      = cfg.gwas_path / "1_BACTERIA/bacteria_metadata.tsv",
        dpi           = getattr(cfg.style, "dpi", 150),
        style         = cfg.style,
    )
else:
    print("\n[skip] Step 2i: phandango render (RUN_PHANDANGO_RENDER = False)")

print("\nDone.")
# SGNH predictor selection → run processing/sgnh-proc/ next.
