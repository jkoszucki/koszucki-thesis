"""
figures/chapter2 — Chapter 2 plots: SSRBH depolymerase GWAS predictors.

Reads from:
  - output_dir/processing/sgnh-hydrolases/
  - output_dir/rbp_depolymerases/
  - output_dir/rbp_best_predictors/
  - output_dir/processing/gwas-data/
Writes plots to scripts/figures/chapter2/plots/.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))

from config import Config
from figure2_3_panelA import plot_figure2_3_panelA
from figure2_3_panelB import plot_figure2_3_panelB

cfg = Config()

plots_dir        = Path(__file__).resolve().parent / "plots"
sgnh_dir         = cfg.output_dir / "processing" / "sgnh-hydrolases"
pyseer_hits_tsv  = cfg.gwas_path / "3_GWAS" / "3_PROCESSING" / "pyseer_hits_filtered.tsv"

# ---------------------------------------------------------------------------
# Run flags
# ---------------------------------------------------------------------------
# Figure 2.7B (N-terminal anchor sharing) is produced by processing/sgnh-proc
# (n_terminal_blast.py, n_terminal_analysis.py, rbps.py) + manual AI assembly —
# see docs/AGENT_CONTEXT.md. The old figure2_7_panelB.py here was superseded
# and referenced removed data (output_dir/enzymes, sgnh-hydrolases/per-pc-best);
# moved to other/unused/.
PLOT_FIGURE2_3_PANELA        = True    # SSRBH depolymerase predictors per K-locus (F1 dot plot)
PLOT_FIGURE2_3_PANELB        = True    # SSRBH depolymerase precision vs recall scatter
PLOT_FIGURE2_6_PANELA        = True    # active depolymerase pairwise BLASTP comparison (virulent vs prophage)
PLOT_FIGURE2_6_PANELB        = True    # Cytoscape network of distinct-specificity depolymerase pairs
PLOT_FIGURE2_6_PANELE        = True    # scatter — alignment length (X) vs protein length (Y), shorter/longer + spans fill; matches manuscript caption

# ---------------------------------------------------------------------------
# Figure 2.3A — SSRBH depolymerase predictor dot plot
# ---------------------------------------------------------------------------
if PLOT_FIGURE2_3_PANELA:
    print("\nFigure 2.3A: SSRBH depolymerase predictor dot plot …")
    plot_figure2_3_panelA(
        depolymerases_gwas_tsv = cfg.output_dir / "rbp_depolymerases" / "depolymerases_gwas.tsv",
        pyseer_hits_tsv        = pyseer_hits_tsv,
        best_predictors_tsv    = cfg.output_dir / "rbp_best_predictors" / "best_predictors_gwas.tsv",
        plots_dir              = plots_dir,
        style                  = cfg.style,
    )
else:
    print("[skip] Figure 2.3A")

# ---------------------------------------------------------------------------
# Figure 2.3B — SSRBH depolymerase precision vs recall scatter
# ---------------------------------------------------------------------------
if PLOT_FIGURE2_3_PANELB:
    print("\nFigure 2.3B: SSRBH depolymerase precision vs recall …")
    plot_figure2_3_panelB(
        depolymerases_gwas_tsv = cfg.output_dir / "rbp_depolymerases" / "depolymerases_gwas.tsv",
        pyseer_hits_tsv        = pyseer_hits_tsv,
        best_predictors_tsv    = cfg.output_dir / "rbp_best_predictors" / "best_predictors_gwas.tsv",
        plots_dir              = plots_dir,
        style                  = cfg.style,
    )
else:
    print("[skip] Figure 2.3B")

# ---------------------------------------------------------------------------
# Figure 2.6A — active depolymerase pairwise BLASTP comparison
# ---------------------------------------------------------------------------
if PLOT_FIGURE2_6_PANELA:
    print("\nFigure 2.6A: active depolymerase pairwise BLASTP …")
    from figure2_6_panelA import plot_figure2_6_panelA
    plot_figure2_6_panelA(
        virulent_active_tsv          = cfg.output_dir / "rbp_depolymerases" / "depolymerases_virulent_active.tsv",
        gwas_tsv                     = cfg.output_dir / "rbp_depolymerases" / "depolymerases_gwas.tsv",
        manualsearch_active_tsv      = cfg.output_dir / "rbp_depolymerases" / "depolymerases_manualsearch_active.tsv",
        manualsearch_inactive_tsv    = cfg.output_dir / "rbp_depolymerases" / "depolymerases_manualsearch_inactive.tsv",
        manualsearch_notproduced_tsv = cfg.output_dir / "rbp_depolymerases" / "depolymerases_manualsearch_notproduced.tsv",
        gwas_active_tsv              = cfg.output_dir / "rbp_depolymerases" / "depolymerases_gwas_active.tsv",
        plots_dir    = plots_dir,
        analysis_dir = plots_dir / "figure2_6_panelA",
        tmp_dir      = plots_dir / "figure2_6_panelA" / "blastp_tmp",
        style        = cfg.style,
    )
else:
    print("[skip] Figure 2.6A")

# ---------------------------------------------------------------------------
# Figure 2.6E — scatter: alignment length (X) vs protein length (Y)
# (matches manuscript caption; the earlier strip-plot design and an older
# unversioned draft are archived in other/unused/ as figure2_6_panelE_v1.py
# and figure2_6_panelE.py)
# ---------------------------------------------------------------------------
if PLOT_FIGURE2_6_PANELE:
    print("\nFigure 2.6E: alignment length vs protein length (scatter) …")
    from figure2_6_panelE import plot_figure2_6_panelE
    plot_figure2_6_panelE(
        blastp_all_hits_tsv = plots_dir / "figure2_6_panelA" / "blastp_all_hits.tsv",
        plots_dir           = plots_dir,
        style               = cfg.style,
    )
else:
    print("[skip] Figure 2.6E")

# ---------------------------------------------------------------------------
# Figure 2.6B — distinct-specificity depolymerase network (Cytoscape)
# ---------------------------------------------------------------------------
if PLOT_FIGURE2_6_PANELB:
    print("\nFigure 2.6B: distinct-specificity depolymerase network …")
    from figure2_6_panelB import plot_figure2_6_panelB
    plot_figure2_6_panelB(
        active_enzymes_tsv = plots_dir / "figure2_6_panelA" / "active_enzymes.tsv",
        hits_tsv           = plots_dir / "figure2_6_panelA" / "figure2_6_panelA_table.tsv",
        plots_dir          = plots_dir,
        style              = cfg.style,
    )
else:
    print("[skip] Figure 2.6B")

print("\nDone.")
