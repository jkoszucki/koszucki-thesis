"""
processing/sgnh-proc — SGNH hydrolase predictor selection and organisation.

Reads gwas_hits.tsv produced by processing/gwas-proc/.
Filters for ecod_type == "sgnh-ecod" and writes:
  1. gwas_sgnh_hits.tsv  — all SGNH hits
  2. gwas_sgnh_best.tsv  — best predictor per K-locus (highest F1, with representative sequence)

Output root: cfg.output_dir / "processing" / "sgnh-hydrolases"
Feeds figures/chapter2/.

per-pc-hits/ and per-pc-best/ symlink trees (and the MAFFT alignment step that
depended on them) are retired — unused in the manuscript.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))

from config import Config
from sgnh_select import select_sgnh
from n_terminal_blast import blast_n_terminal
from n_terminal_analysis import analyse_n_terminal
from rbps import build_rbp_folder

cfg = Config()

gwas_hits_tsv    = cfg.output_dir / "processing" / "gwas-data" / "gwas_hits.tsv"
mmseqs_dir       = cfg.gwas_path / "3_GWAS/1_INTERMEDIATE/2_MMSEQS"
sgnh_dir         = cfg.output_dir / "processing" / "sgnh-hydrolases"
rbp_deacetylases_dir = cfg.output_dir / "rbp_deacetylases"

# N-terminal anchor analysis — standalone output under other/n-terminal/
n_terminal_dir   = cfg.output_dir / "other" / "n-terminal"
n_terminal_pc_fasta = (
    cfg.gwas_path / "3_GWAS/4_ANALYZE/1_PER_LOCUS/KL24/lasso/alignments/PCI80C50/PC0262.fasta"
)
n_terminal_genbank_dir = cfg.gwas_path / "2_PROPHAGES/1_GENBANK_GENOMES"

# ---------------------------------------------------------------------------
# Run flags
# ---------------------------------------------------------------------------
RUN_N_TERMINAL_BLAST    = True    # True → BLASTP N-terminal fragment vs KPSC prophages; skips if raw_blast.tsv exists
RUN_N_TERMINAL_ANALYSIS = True    # True → filter hits, map clustering, select representatives; skips if 3_one_per_cluster.tsv exists
RUN_RBP_FOLDER          = True    # True → build n-terminal/rbps/ (fasta, json, cif symlinks, png renders)
RENDER_RBP_STRUCTURES   = True    # True → render highlighted PNGs (pymol env); False → skip renders only

# ---------------------------------------------------------------------------
# Step 1 — select SGNH hits, write TSVs, recreate symlinks
# ---------------------------------------------------------------------------
print("Step 1: Selecting SGNH hydrolase predictors …")
select_sgnh(
    gwas_hits_tsv       = gwas_hits_tsv,
    mmseqs_dir          = mmseqs_dir,
    sgnh_dir            = sgnh_dir,
    precision_threshold = 0.50,
)

# Write deacetylases_gwas.tsv to rbp_deacetylases/
import shutil, pandas as pd
rbp_deacetylases_dir.mkdir(parents=True, exist_ok=True)
_src = sgnh_dir / "gwas_sgnh_best.tsv"
_dst = rbp_deacetylases_dir / "deacetylases_gwas.tsv"
shutil.copy2(_src, _dst)
print(f"  deacetylases_gwas.tsv — {len(pd.read_csv(_dst, sep=chr(9)))} rows → {_dst}")

# Phandango files are exported by processing/gwas-proc/ (Step 2h)
# and are accessible under processing/gwas-data/ directly.

# ---------------------------------------------------------------------------
# Step 2 — BLAST N-terminal fragment vs KPSC prophage proteins
# ---------------------------------------------------------------------------
if RUN_N_TERMINAL_BLAST:
    print("\nStep 2: BLASTing N-terminal fragment vs KPSC prophage proteins …")
    blast_n_terminal(
        prophage_faa_pattern = str(cfg.gwas_path / "2_PROPHAGES/4_FASTA_CDS_AA/*.faa"),
        blast_db_dir         = cfg.output_dir / "other" / "prophage_blastp_db",
        n_terminal_dir       = n_terminal_dir,
        pc_fasta             = n_terminal_pc_fasta,
    )
else:
    print("\n[skip] Step 2: N-terminal BLAST (RUN_N_TERMINAL_BLAST = False)")

# ---------------------------------------------------------------------------
# Step 3 — N-terminal BLAST analysis: filter, cluster mapping, representative selection
# ---------------------------------------------------------------------------
if RUN_N_TERMINAL_ANALYSIS:
    print("\nStep 3: Analysing N-terminal BLAST hits …")
    analyse_n_terminal(
        prophage_metadata_tsv = cfg.gwas_path / "2_PROPHAGES/prophages_metadata.tsv",
        bacteria_metadata_tsv = cfg.gwas_path / "1_BACTERIA/bacteria_metadata.tsv",
        sgnh_dir              = n_terminal_dir.parent,   # other/ → n-terminal/ appended inside
        genbank_dir           = n_terminal_genbank_dir,
    )
else:
    print("\n[skip] Step 3: N-terminal analysis (RUN_N_TERMINAL_ANALYSIS = False)")

# ---------------------------------------------------------------------------
# Step 4 — build n-terminal/prophages_rbps/ folder
# ---------------------------------------------------------------------------
if RUN_RBP_FOLDER:
    print("\nStep 4: Building n-terminal/prophages_rbps/ folder …")
    build_rbp_folder(
        prophage_faa_pattern = str(cfg.gwas_path / "2_PROPHAGES/4_FASTA_CDS_AA/*.faa"),
        blast_db_dir         = cfg.output_dir / "other" / "prophage_blastp_db",
        af3_raw_dir          = cfg.input_dir / "supplementary-thesis" / "supplementary-data" / "S4_Data" / "1_INPUT" / "1_RAW_ALPHAFOLD3",
        n_terminal_dir       = n_terminal_dir,
        render_structures    = RENDER_RBP_STRUCTURES,
    )
else:
    print("\n[skip] Step 4: RBP folder (RUN_RBP_FOLDER = False)")

print("\nDone.")
