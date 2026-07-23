"""
figures/chapter3 — Chapter 3 plots: SGNH hydrolase deacetylases + acetyltransferases.

Reads from:
  - output_dir/processing/sgnh-hydrolases/
  - output_dir/enzymes/
Writes plots to scripts/figures/chapter3/plots/.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))

from config import Config
from figure3_1_panelA import plot_figure3_1_panelA
from figure3_1_panelB import plot_figure3_1_panelB
from figure3_2_panelA import plot_figure3_2_panelA
from figure3_2_panelB import plot_figure3_2_panelB

cfg = Config()

plots_dir       = Path(__file__).resolve().parent / "plots"
sgnh_dir        = cfg.output_dir / "processing" / "sgnh-hydrolases"
enzymes_xlsx    = cfg.input_dir / "supplementary-thesis" / "supplementary-tables" / "S1_Table.xlsx"
pyseer_hits_tsv = cfg.gwas_path / "3_GWAS" / "3_PROCESSING" / "pyseer_hits_filtered.tsv"

# ---------------------------------------------------------------------------
# Run flags
# ---------------------------------------------------------------------------
PLOT_FIGURE3_1_PANELA = True   # SGNH predictor recall dot plot (12 K-loci, precision >= 0.50)
PLOT_FIGURE3_1_PANELB = True   # SGNH similarity network: 12 predictions + 4 experimental (Cytoscape)
PLOT_FIGURE3_2_PANELA = True   # GWAS prediction + O-acetylation overview table (35 K-loci)
PLOT_FIGURE3_2_PANELB = True   # Acetyltransferase SSN (69 K-locus + 3 experimental + 2 GWAS best)

# ---------------------------------------------------------------------------
# Figure 3.1A — SGNH predictor recall dot plot
# ---------------------------------------------------------------------------
if PLOT_FIGURE3_1_PANELA:
    print("\nFigure 3.1A: SGNH predictor recall dot plot …")
    plot_figure3_1_panelA(
        gwas_sgnh_best_tsv  = sgnh_dir / "gwas_sgnh_best.tsv",
        best_predictors_tsv = cfg.output_dir / "rbp_best_predictors" / "best_predictors_gwas.tsv",
        pyseer_hits_tsv     = pyseer_hits_tsv,
        plots_dir           = plots_dir,
        style               = cfg.style,
    )
else:
    print("[skip] Figure 3.1A")

# ---------------------------------------------------------------------------
# Figure 3.1B — SGNH similarity network (Cytoscape)
# ---------------------------------------------------------------------------
if PLOT_FIGURE3_1_PANELB:
    print("\nFigure 3.1B: SGNH similarity network (Cytoscape) …")
    plot_figure3_1_panelB(
        gwas_sgnh_best_tsv = sgnh_dir / "gwas_sgnh_best.tsv",
        enzymes_xlsx       = enzymes_xlsx,
        plots_dir          = plots_dir,
        style              = cfg.style,
    )
else:
    print("[skip] Figure 3.1B")

# ---------------------------------------------------------------------------
# Figure 3.2X — GWAS prediction + O-acetylation overview table
# ---------------------------------------------------------------------------
if PLOT_FIGURE3_2_PANELA:
    print("\nFigure 3.2A: GWAS prediction + O-acetylation overview table …")
    plot_figure3_2_panelA(
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
    print("[skip] Figure 3.2A")

# ---------------------------------------------------------------------------
# Figure 3.2B — Acetyltransferase similarity network (Cytoscape)
# ---------------------------------------------------------------------------
if PLOT_FIGURE3_2_PANELB:
    print("\nFigure 3.2B: Acetyltransferase similarity network (Cytoscape) …")
    _no_ecod_dir = cfg.output_dir / "processing" / "gwas-data" / "no-ecod-reported-topology"
    plot_figure3_2_panelB(
        acetylases_kloci_tsv      = cfg.output_dir / "cps_acetylases" / "acetylases_kloci.tsv",
        acetylases_literature_tsv = cfg.output_dir / "cps_acetylases" / "acetylases_literature_active.tsv",
        gwas_best_at_fastas       = {
            "KL30":  _no_ecod_dir / "KL30"  / "PCI80C80" / "PC0675" / "protein" / "sequence.fasta",
            "KL111": _no_ecod_dir / "KL111" / "PCI80C80" / "PC0915" / "protein" / "sequence.fasta",
        },
        plots_dir = plots_dir,
        style     = cfg.style,
    )
else:
    print("[skip] Figure 3.2B")

print("\nDone.")
