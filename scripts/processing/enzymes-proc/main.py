"""
processing/enzymes-proc — experimental enzyme data preparation.

Pipeline:
  1. build_literature_tables() → rbp_deacetylases/deacetylases_literature.tsv,
                                  cps_acetylases/acetylases_literature_active.tsv

Both tables already carry a `sequence` column, and the corresponding AF3
structures are indexed under output_dir/other/alphafold3/ (built + maintained by
scripts/helpers/build_af3_index.py) — so there is no separate per-protein
output tree here.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))

from config import Config
from literature_tables import build_literature_tables

cfg = Config()

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
enzymes_xlsx = cfg.input_dir / "supplementary-thesis" / "supplementary-tables" / "S1_Table.xlsx"

# ---------------------------------------------------------------------------
# Step 1 — literature tables from S1_Table.xlsx
# ---------------------------------------------------------------------------
print("Step 1: Building literature tables from S1_Table.xlsx …")
build_literature_tables(
    enzymes_xlsx         = enzymes_xlsx,
    rbp_deacetylases_dir = cfg.output_dir / "rbp_deacetylases",
    cps_acetylases_dir   = cfg.output_dir / "cps_acetylases",
)

print("\nDone.")
