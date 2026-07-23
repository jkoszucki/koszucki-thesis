"""
Render AF3 best models from 2_BEST_MODELS → 1_DRAWN.

Reads:  input_dir/supplementary-thesis/supplementary-data/S4_Data/2_BEST_MODELS/*.cif  (symlinks; read-only)
        output_dir/other/alphafold3/summary.tsv                                        (assigned_set + ecod_type for color)
Writes: output_dir/other/alphafold3/1_DRAWN/{stem}.png

Color priority per canonical_name (from status_per_protein.tsv):
  1. Experimentally active set (ACTIVE_SETS) → dark color from af3_set_colors (−30% lightness)
  2. best_predictors_gwas → light color by ecod_type (sgnh→gold, ssrbh→green)
  3. Any other set → light color from af3_set_colors
  If a protein appears in both an active and a non-active set, the active color wins.

Run flags
---------
RENDER_MODE    — "fast" for quick checks, "accurate" for thesis
ORIENT_NTERM   — rotate each structure so chain A N-terminal CA is at top centre
N_STRUCTURES   — limit to first N structures (None = all); useful for spot-checks
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import Config
from structure_renderer import StructureRenderer

cfg = Config()

PATH2    = cfg.input_dir  / "supplementary-thesis" / "supplementary-data" / "S4_Data" / "2_BEST_MODELS"
PATH3    = cfg.output_dir / "other" / "alphafold3" / "1_DRAWN"
SUMMARY  = cfg.output_dir / "other" / "alphafold3" / "status_per_protein.tsv"

# Sets that contain experimentally validated active proteins — rendered darker
ACTIVE_SETS = frozenset({
    "depolymerases_virulent_active",
    "depolymerases_manualsearch_active",
    "depolymerases_gwas_active",
    "deacetylases_literature",
    "deacetylases_gwas",
    "acetylases_literature_active",
})

# ---------------------------------------------------------------------------
# Run flags
# ---------------------------------------------------------------------------
RENDER_MODE   = "accurate"
ORIENT_NTERM  = True
N_STRUCTURES  = None  # None = all; set to e.g. 3 for a spot-check

# ---------------------------------------------------------------------------
# Build canonical_name → color from summary.tsv
# ---------------------------------------------------------------------------
summary = pd.read_csv(SUMMARY, sep="\t")

name_to_color: dict[str, str] = {}
for canonical, group in summary.groupby("canonical_name"):
    active = group[group["status"] != "EXCLUDED"]
    if active.empty:
        name_to_color[canonical] = cfg.style.gray_color
        continue

    # Priority 1: experimentally active sets → dark color from af3_set_colors
    active_exp = active[active["assigned_set"].isin(ACTIVE_SETS)]
    if not active_exp.empty:
        chosen = active_exp.iloc[0]["assigned_set"]
        name_to_color[canonical] = cfg.style.af3_set_colors.get(chosen, cfg.style.gray_color)
        continue

    # Priority 2: best_predictors_gwas → light color resolved by ecod_type
    bp = active[active["assigned_set"] == "best_predictors_gwas"]
    if not bp.empty:
        ecod = str(bp.iloc[0].get("ecod_type", "")).strip()
        if ecod == "sgnh-ecod":
            name_to_color[canonical] = cfg.style.sgnh_domain_color
        else:
            name_to_color[canonical] = cfg.style.ssrbh_color
        continue

    # Priority 3: any other set → light color from af3_set_colors
    chosen = active.iloc[0]["assigned_set"]
    name_to_color[canonical] = cfg.style.af3_set_colors.get(chosen, cfg.style.gray_color)

# ---------------------------------------------------------------------------
# Render
# ---------------------------------------------------------------------------
cifs = sorted(PATH2.glob("*.cif"))
if N_STRUCTURES is not None:
    cifs = cifs[:N_STRUCTURES]

if not cifs:
    print(f"No CIF symlinks found in {PATH2}")
    sys.exit(0)

# Group by colour so each batch shares the same renderer settings
groups: dict[str, list[Path]] = defaultdict(list)
for cif in cifs:
    canonical = cif.stem.removesuffix("_model_0")
    color = name_to_color.get(canonical, cfg.style.gray_color)
    groups[color].append(cif)

for color, group_cifs in groups.items():
    renderer = StructureRenderer(
        render_mode  = RENDER_MODE,
        color_mode   = "gray",
        gray_color   = color,
        orient_nterm = ORIENT_NTERM,
    )
    renderer.render_to_dir(group_cifs, out_dir=PATH3)
