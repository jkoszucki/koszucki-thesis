"""
figures/chapter4 — Chapter 4 plots: SGNH hydrolase deacetylases + acetyltransferases.

Reads from:
  - output_dir/processing/sgnh-hydrolases/, output_dir/processing/gwas-data/
  - output_dir/rbp_best_predictors/, rbp_depolymerases/, rbp_deacetylases/, cps_acetylases/
  - input_dir/supplementary-thesis/supplementary-tables/S1_Table.xlsx + S2_Table.xlsx
  - input_dir/supplementary-thesis/supplementary-data/S3_Data/
Writes plots to scripts/figures/chapter4/plots/.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))

from config import Config
from figure4_1_panelA import plot_figure4_1_panelA
from figure4_1_panelB import plot_figure4_1_panelB
from figure4_2_panelA import plot_figure4_2_panelA
from figure4_2_panelB import plot_figure4_2_panelB
from figure4_2_panelC_heatmap import plot_figure4_2_panelC_heatmap
from figure4_2_panelC_representatives import plot_figure4_2_panelC_representatives
from figure4_2_panelC_tmalign import plot_figure4_2_panelC_tmalign

cfg = Config()

plots_dir       = Path(__file__).resolve().parent / "plots"
sgnh_dir        = cfg.output_dir / "processing" / "sgnh-hydrolases"
enzymes_xlsx    = cfg.input_dir / "supplementary-thesis" / "supplementary-tables" / "S1_Table.xlsx"
pyseer_hits_tsv = cfg.gwas_path / "3_GWAS" / "3_PROCESSING" / "pyseer_hits_filtered.tsv"

# ---------------------------------------------------------------------------
# Run flags
# ---------------------------------------------------------------------------
PLOT_FIGURE4_1_PANELA = True   # SGNH predictor recall dot plot (12 K-loci, precision >= 0.50)
PLOT_FIGURE4_1_PANELB = True   # SGNH similarity network: 12 predictions + 4 experimental (Cytoscape)
PLOT_FIGURE4_2_PANELA = True   # GWAS prediction + O-acetylation overview table (35 K-loci)
PLOT_FIGURE4_2_PANELB = True   # Functional annotation of acetyltransferases (Foldseek hits; 69 K-locus + 4 reference)
PLOT_FIGURE4_2_PANELC = True   # Acetyltransferase TM-score matrix: 73 proteins collapsed to 19 representatives
PLOT_FIGURE4_2_PANELC_TMALIGN = True   # Supporting all-vs-all structural similarity network (Cytoscape, 73 proteins)

# ---------------------------------------------------------------------------
# Figure 4.1A — SGNH predictor recall dot plot
# ---------------------------------------------------------------------------
if PLOT_FIGURE4_1_PANELA:
    print("\nFigure 4.1A: SGNH predictor recall dot plot …")
    plot_figure4_1_panelA(
        gwas_sgnh_best_tsv  = sgnh_dir / "gwas_sgnh_best.tsv",
        best_predictors_tsv = cfg.output_dir / "rbp_best_predictors" / "best_predictors_gwas.tsv",
        pyseer_hits_tsv     = pyseer_hits_tsv,
        plots_dir           = plots_dir,
        style               = cfg.style,
    )
else:
    print("[skip] Figure 4.1A")

# ---------------------------------------------------------------------------
# Figure 4.1B — SGNH similarity network (Cytoscape)
# ---------------------------------------------------------------------------
if PLOT_FIGURE4_1_PANELB:
    print("\nFigure 4.1B: SGNH similarity network (Cytoscape) …")
    plot_figure4_1_panelB(
        gwas_sgnh_best_tsv = sgnh_dir / "gwas_sgnh_best.tsv",
        enzymes_xlsx       = enzymes_xlsx,
        plots_dir          = plots_dir,
        style              = cfg.style,
    )
else:
    print("[skip] Figure 4.1B")

# ---------------------------------------------------------------------------
# Figure 4.2X — GWAS prediction + O-acetylation overview table
# ---------------------------------------------------------------------------
if PLOT_FIGURE4_2_PANELA:
    print("\nFigure 4.2A: GWAS prediction + O-acetylation overview table …")
    plot_figure4_2_panelA(
        best_predictors_tsv    = cfg.output_dir / "rbp_best_predictors" / "best_predictors_gwas.tsv",
        depolymerases_gwas_tsv = cfg.output_dir / "rbp_depolymerases"  / "depolymerases_gwas.tsv",
        deacetylases_gwas_tsv  = cfg.output_dir / "rbp_deacetylases"   / "deacetylases_gwas.tsv",
        acetylases_gwas_tsv    = cfg.output_dir / "cps_acetylases"     / "acetylases_gwas.tsv",
        acetylases_kloci_tsv   = cfg.output_dir / "cps_acetylases"     / "acetylases_kloci.tsv",
        cps_xlsx               = cfg.input_dir  / "supplementary-thesis" / "supplementary-tables" / "S2_Table.xlsx",
        pyseer_hits_tsv        = pyseer_hits_tsv,
        plots_dir              = plots_dir,
        style                  = cfg.style,
    )
else:
    print("[skip] Figure 4.2A")

# ---------------------------------------------------------------------------
# Figure 4.2B — Functional annotation of putative acetyltransferases
# ---------------------------------------------------------------------------
if PLOT_FIGURE4_2_PANELB:
    print("\nFigure 4.2B: Functional annotation of putative acetyltransferases …")
    plot_figure4_2_panelB(
        foldseek_dir        = cfg.input_dir / "supplementary-thesis" / "supplementary-data" / "S3_Data",
        structure_order_tsv = cfg.output_dir / "cps_acetylases" / "tmalign" / "acetylases_kloci_structure_order.tsv",
        plots_dir           = plots_dir,
        style               = cfg.style,
    )
else:
    print("[skip] Figure 4.2B")

# ---------------------------------------------------------------------------
# Figure 4.2C — TM-score matrix of the 18 structural representatives
#
# Two steps, in order: the representative selection writes the node table that names the
# 18 rows/columns, and the heatmap reads it back. Running the heatmap alone against a
# stale node table would plot the wrong protein set, so both sit behind one flag.
# ---------------------------------------------------------------------------
if PLOT_FIGURE4_2_PANELC:
    _tmalign_dir = cfg.output_dir / "cps_acetylases" / "tmalign"
    print("\nFigure 4.2C: structural representatives (TM >= 0.75 clusters) …")
    plot_figure4_2_panelC_representatives(
        tmscore_tsv         = _tmalign_dir / "acetylases_kloci_tmscore.tsv",
        structure_order_tsv = _tmalign_dir / "acetylases_kloci_structure_order.tsv",
        plots_dir           = plots_dir,
    )
    print("\nFigure 4.2C: representative TM-score matrix …")
    plot_figure4_2_panelC_heatmap(
        tmscore_tsv              = _tmalign_dir / "acetylases_kloci_tmscore.tsv",
        representatives_node_tsv = plots_dir / "figure4_2-panelC" / "representatives" / "node.tsv",
        plots_dir                = plots_dir,
        style                    = cfg.style,
    )
else:
    print("[skip] Figure 4.2C")

# ---------------------------------------------------------------------------
# Figure 4.2C (structure half) — TM-align structural similarity network
# ---------------------------------------------------------------------------
if PLOT_FIGURE4_2_PANELC_TMALIGN:
    print("\nFigure 4.2C: Acetyltransferase structural similarity network (TM-align) …")
    _tmalign_dir = cfg.output_dir / "cps_acetylases" / "tmalign"
    plot_figure4_2_panelC_tmalign(
        tmscore_tsv          = _tmalign_dir / "acetylases_kloci_tmscore.tsv",
        structure_order_tsv  = _tmalign_dir / "acetylases_kloci_structure_order.tsv",
        acetylases_kloci_tsv = cfg.output_dir / "cps_acetylases" / "acetylases_kloci.tsv",
        enzymes_xlsx         = enzymes_xlsx,
        plots_dir            = plots_dir,
    )
else:
    print("[skip] Figure 4.2C (TM-align)")

print("\nDone.")
