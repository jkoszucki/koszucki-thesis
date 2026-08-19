"""
figures/chapter3 — Chapter 3 plots: CPS K-type diversity and O-acetylation.

Reads from:
  - output_dir/cps_structures/  (built by processing/cps-proc/main.py — run that first)
Writes plots to scripts/figures/chapter3/plots/.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))

from config import Config
from figure3_1 import KTypePlotAPI
from figure3_1_panelC import plot_figure3_1_panelC
from figure3_2_panelA import plot_figure3_2_panelA
from figure3_2_panelB import plot_figure3_2_panelB
from figure3_2_panelC import plot_figure3_2_panelC

cfg = Config()

analysis_output_dir = cfg.output_dir / "cps_structures"
plots_dir            = Path(__file__).resolve().parent / "plots"

# ---------------------------------------------------------------------------
# Run flags — toggle steps without modifying config.yml
# ---------------------------------------------------------------------------
PLOT_FIGURE3_1_PANELA = True   # cumulative NMR structures over time
PLOT_FIGURE3_1_PANELC = True   # modification frequency grid
PLOT_FIGURE3_1_PANELD = True   # OAc/OPy distribution per monosaccharide
PLOT_FIGURE3_2_PANELA = True   # J_core vs J_branch scatter
PLOT_FIGURE3_2_PANELB = True   # CPS structural similarity network (Cytoscape)
PLOT_FIGURE3_2_PANELC = True   # summary table of structurally related CPS pairs

# ---------------------------------------------------------------------------
# Read prepared tables — built by processing/cps-proc/main.py, not here
# ---------------------------------------------------------------------------
plots_dir.mkdir(parents=True, exist_ok=True)

processed_csv     = analysis_output_dir / "ktypes.csv"
similarity_csv    = analysis_output_dir / "ktypes_sim.csv"
modifications_csv = analysis_output_dir / "ktypes_modifications.csv"

for csv_path in (processed_csv, similarity_csv, modifications_csv):
    if not csv_path.exists():
        raise FileNotFoundError(
            f"{csv_path} not found — run `conda run -n jkoszucki python "
            "scripts/processing/cps-proc/main.py` first to build the CPS tables."
        )

# ---------------------------------------------------------------------------
# Figure 3.1A — cumulative NMR structures over time
# ---------------------------------------------------------------------------
if PLOT_FIGURE3_1_PANELA:
    print("\nFigure 3.1A: cumulative CPS structures over time …")
    plot_api = KTypePlotAPI(
        ktypes_csv=processed_csv,
        modifications_csv=modifications_csv,
        style=cfg.style,
    )
    for ext in ("png", "pdf"):
        plot_api.plot_cumulative_structures(
            output_path=plots_dir / f"figure3_1-panelA.{ext}",
            show=False,
        )
else:
    print("[skip] Figure 3.1A")

# ---------------------------------------------------------------------------
# Figure 3.1C — count + mean frequency grid for four modification types
# ---------------------------------------------------------------------------
if PLOT_FIGURE3_1_PANELC:
    print("\nFigure 3.1C: modification frequency grid …")
    for ext in ("png", "pdf"):
        plot_figure3_1_panelC(
            modifications_csv=modifications_csv,
            output_path=plots_dir / f"figure3_1-panelC.{ext}",
            style=cfg.style,
        )
else:
    print("[skip] Figure 3.1C")

# ---------------------------------------------------------------------------
# Figure 3.1D — OAc/OPy distribution per monosaccharide
# ---------------------------------------------------------------------------
if PLOT_FIGURE3_1_PANELD:
    print("\nFigure 3.1D: OAc/OPy distribution per monosaccharide …")
    plot_api = KTypePlotAPI(
        ktypes_csv=processed_csv,
        modifications_csv=modifications_csv,
        style=cfg.style,
    )
    for ext in ("png", "pdf"):
        plot_api.plot_modification_and_monosaccharide_distribution(
            output_path=plots_dir / f"figure3_1-panelD.{ext}",
            show=False,
        )
else:
    print("[skip] Figure 3.1D")

# ---------------------------------------------------------------------------
# Figure 3.2A — J_core vs J_branch scatter (all K-type pairs)
# ---------------------------------------------------------------------------
if PLOT_FIGURE3_2_PANELA:
    print("\nFigure 3.2A: J_core vs J_branch scatter …")
    for ext in ("png", "pdf"):
        plot_figure3_2_panelA(
            similarity_csv=similarity_csv,
            output_path=plots_dir / f"figure3_2-panelA.{ext}",
            style=cfg.style,
        )
else:
    print("[skip] Figure 3.2A")

# ---------------------------------------------------------------------------
# Figure 3.2B — CPS structural similarity network (Cytoscape)
# ---------------------------------------------------------------------------
if PLOT_FIGURE3_2_PANELB:
    print("\nFigure 3.2B: CPS similarity network (Cytoscape) …")
    plot_figure3_2_panelB(
        similarity_csv=similarity_csv,
        ktypes_csv=processed_csv,
        modifications_csv=modifications_csv,
        output_dir=plots_dir / "figure3_2-panelB",
        style=cfg.style,
    )
else:
    print("[skip] Figure 3.2B")

# ---------------------------------------------------------------------------
# Figure 3.2C — summary table of structurally related CPS pairs
# ---------------------------------------------------------------------------
if PLOT_FIGURE3_2_PANELC:
    print("\nFigure 3.2C: CPS pairs summary table …")
    for ext in ("png", "pdf"):
        plot_figure3_2_panelC(
            output_path=plots_dir / f"figure3_2-panelC.{ext}",
            style=cfg.style,
        )
else:
    print("[skip] Figure 3.2C")

print("\nDone.")
