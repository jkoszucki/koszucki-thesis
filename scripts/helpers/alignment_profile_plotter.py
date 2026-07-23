"""
Shared alignment gap-frequency profile plotter — importable by any chapter script.

Computes per-column gap frequencies from one or more protein alignments,
bins them into equal-width position bins, and plots individual traces
(pale gray) + median + IQR ribbon per panel.

Each panel can receive unaligned sequences (MAFFT is called automatically)
or pre-aligned sequences (MAFFT skipped).

Usage (in any figures/chapterX/main.py or lib/)
------------------------------------------------
    sys.path.insert(0, str(Path(__file__).resolve().parents[N] / "helpers"))
    from alignment_profile_plotter import AlignmentProfilePlotter

    plotter = AlignmentProfilePlotter(
        n_bins           = 30,
        gap_chars        = "-.",
        # Individual trace appearance
        trace_color      = "#cccccc",
        trace_alpha      = 0.30,
        trace_linewidth  = 0.8,
        # Median + IQR ribbon appearance
        ribbon_alpha     = 0.30,
        median_linewidth = 2.5,
        # Layout
        figsize          = None,     # None → auto (10 wide, 5 per panel)
        panel_hspace     = 0.35,
        style            = cfg.style,
    )

    plotter.plot(
        panels = [
            {
                "sequences":   fasta_path_or_list_of_seqs,
                "title":       "MMseqs2 alignments",
                "color":       "#1f77b4",     # optional; cycles through defaults if omitted
                "pre_aligned": False,         # optional; True → skip MAFFT
            },
            {
                "sequences":   another_fasta_path,
                "title":       "Prophage BLAST alignments",
                "color":       "#d62728",
                "pre_aligned": True,
            },
        ],
        out_path = plots_dir / "gap_frequency_profile.png",
    )

Panel specification
-------------------
Each panel dict accepts:
    sequences   : list[str] | Path
                  Unaligned sequences (list of amino-acid strings or FASTA path) if
                  pre_aligned=False; already-aligned sequences if pre_aligned=True.
    title       : str        Panel title (appended with n= and bin count).
    color       : str        Hex color for median line and IQR ribbon.
                             Optional; defaults cycle through _DEFAULT_COLORS.
    pre_aligned : bool       If True, skip MAFFT. Default: False.
"""

from __future__ import annotations

import io
import subprocess
import tempfile
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from Bio import SeqIO


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

_DEFAULT_COLORS = ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd"]
_GAP_CHARS_DEFAULT = set("-.")


# ---------------------------------------------------------------------------
# AlignmentProfilePlotter
# ---------------------------------------------------------------------------

class AlignmentProfilePlotter:
    """
    Gap-frequency profile plot from one or more protein alignments.

    Parameters
    ----------
    n_bins : int
        Number of equal-width bins over normalised position 0–100%.
        Default: 30.
    gap_chars : str
        Characters counted as gaps when computing per-column gap frequency.
        Default: "-.".
    trace_color : str
        Hex color of individual per-predictor profile lines. Default: "#cccccc".
    trace_alpha : float
        Opacity of individual traces. Default: 0.30.
    trace_linewidth : float
        Line width of individual traces. Default: 0.8.
    ribbon_alpha : float
        Opacity of the IQR ribbon. Default: 0.30.
    median_linewidth : float
        Line width of the median line. Default: 2.5.
    figsize : tuple[float, float] | None
        Figure dimensions (width, height) in inches. If None, auto-sized to
        (10, 5 × n_panels). Default: None.
    panel_hspace : float
        Vertical spacing between panels (matplotlib hspace). Default: 0.35.
    style : cfg.style
        Provides axis_label_fontsize, axis_label_fontweight, tick_fontsize,
        tick_fontweight, dpi. Falls back to hardcoded defaults if None.
    """

    def __init__(
        self,
        n_bins:           int                      = 30,
        gap_chars:        str                      = "-.",
        trace_color:      str                      = "#cccccc",
        trace_alpha:      float                    = 0.30,
        trace_linewidth:  float                    = 0.8,
        ribbon_alpha:     float                    = 0.30,
        median_linewidth: float                    = 2.5,
        figsize:          tuple[float, float] | None = None,
        panel_hspace:     float                    = 0.35,
        style:            Any                      = None,
    ) -> None:
        self.n_bins           = n_bins
        self.gap_chars        = set(gap_chars)
        self.trace_color      = trace_color
        self.trace_alpha      = trace_alpha
        self.trace_linewidth  = trace_linewidth
        self.ribbon_alpha     = ribbon_alpha
        self.median_linewidth = median_linewidth
        self.figsize          = figsize
        self.panel_hspace     = panel_hspace

        self.label_fs = getattr(style, "axis_label_fontsize",   12)
        self.label_fw = getattr(style, "axis_label_fontweight", "bold")
        self.tick_fs  = getattr(style, "tick_fontsize",          8)
        self.tick_fw  = getattr(style, "tick_fontweight",       "bold")
        self.dpi      = getattr(style, "dpi",                   300)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def plot(self, panels: list[dict], out_path: Path) -> None:
        """
        Compute profiles for each panel and save the figure.

        Args:
            panels:    List of panel specification dicts. Each dict:
                         sequences   — list[str] | Path (FASTA)
                         title       — str
                         color       — str hex (optional)
                         pre_aligned — bool (optional, default False)
            out_path:  Output PNG path.
        """
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)

        computed = []
        for i, panel in enumerate(panels):
            seqs = self._resolve_sequences(panel)
            profiles = self._compute_profiles(seqs, panel.get("pre_aligned", False))
            if not profiles:
                print(f"  [alignment] panel '{panel.get('title', i)}': "
                      "no valid profiles — skipping.")
                continue
            computed.append({
                "profiles": profiles,
                "title":    panel.get("title", f"Panel {i + 1}"),
                "color":    panel.get("color", _DEFAULT_COLORS[i % len(_DEFAULT_COLORS)]),
            })

        if not computed:
            print("  [alignment] no profiles computed — figure not saved.")
            return

        self._draw(computed, out_path)
        print(f"  Saved → {out_path}")

    # ------------------------------------------------------------------
    # Sequence loading
    # ------------------------------------------------------------------

    @staticmethod
    def _resolve_sequences(panel: dict) -> list[str]:
        """
        Return a flat list of amino-acid strings from a panel dict.
        Accepts list[str] or a FASTA Path.
        """
        seqs = panel["sequences"]
        if isinstance(seqs, (str, Path)):
            return [str(r.seq) for r in SeqIO.parse(seqs, "fasta")]
        return list(seqs)

    # ------------------------------------------------------------------
    # Profile computation
    # ------------------------------------------------------------------

    def _compute_profiles(
        self, sequences: list[str], pre_aligned: bool
    ) -> list[np.ndarray]:
        """
        Align (unless pre_aligned) and compute a gap-frequency profile
        for each individual sequence in the batch.

        When sequences share a common alignment (e.g. all from one FASTA),
        treat the whole batch as one profile. When sequences come from multiple
        separate FASTAs, call this once per FASTA and collect the result.

        Here we treat the entire list as a single alignment → one profile.
        To get per-FASTA profiles, call _compute_profiles once per FASTA.
        """
        if len(sequences) < 2:
            print(f"  [alignment] fewer than 2 sequences — skipping.")
            return []

        if pre_aligned:
            aligned = sequences
        else:
            aligned = self._align(sequences)
            if aligned is None:
                return []

        lengths = {len(s) for s in aligned}
        if len(lengths) > 1:
            print("  [alignment] unequal lengths after alignment — skipping.")
            return []

        return [self._gap_profile(aligned)]

    def _align(self, sequences: list[str]) -> list[str] | None:
        """Write sequences to a temp FASTA and run MAFFT --auto."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            for i, seq in enumerate(sequences):
                f.write(f">seq{i}\n{seq}\n")
            tmp_path = f.name
        try:
            result = subprocess.run(
                ["mafft", "--auto", "--quiet", tmp_path],
                capture_output=True, text=True, check=True,
            )
            return [str(r.seq) for r in SeqIO.parse(io.StringIO(result.stdout), "fasta")]
        except subprocess.CalledProcessError as e:
            print(f"  [alignment] MAFFT failed: {e}")
            return None
        finally:
            Path(tmp_path).unlink(missing_ok=True)

    def _gap_profile(self, seqs: list[str]) -> np.ndarray:
        """Bin per-column gap frequencies into n_bins equal-width position bins."""
        aln_len  = len(seqs[0])
        n_seqs   = len(seqs)
        gap_freq = np.array([
            sum(1 for s in seqs if s[i] in self.gap_chars) / n_seqs
            for i in range(aln_len)
        ])
        positions  = np.linspace(0, 100, aln_len)
        bin_edges  = np.linspace(0, 100, self.n_bins + 1)
        bin_means  = []
        for b in range(self.n_bins):
            lo, hi   = bin_edges[b], bin_edges[b + 1]
            is_last  = b == self.n_bins - 1
            mask     = (positions >= lo) & (positions <= hi if is_last else positions < hi)
            bin_means.append(gap_freq[mask].mean() if mask.sum() > 0 else np.nan)
        return np.array(bin_means)

    # ------------------------------------------------------------------
    # Drawing
    # ------------------------------------------------------------------

    def _draw(self, computed: list[dict], out_path: Path) -> None:
        n = len(computed)
        figsize = self.figsize if self.figsize else (10, 5 * n)

        sns.set_style("whitegrid")

        if n == 1:
            fig, axes = plt.subplots(1, 1, figsize=figsize)
            axes = [axes]
        else:
            fig, axes = plt.subplots(
                n, 1, figsize=figsize, sharex=True,
                gridspec_kw={"hspace": self.panel_hspace},
            )

        for ax, panel_data in zip(axes, computed):
            self._draw_panel(ax, **panel_data)

        # X-label only on bottom axis
        axes[-1].set_xlabel(
            "Normalised alignment position (%)",
            fontsize=self.label_fs, fontweight=self.label_fw,
        )

        plt.tight_layout()
        fig.savefig(out_path, dpi=self.dpi, bbox_inches="tight")
        plt.close(fig)

    def _draw_panel(
        self,
        ax,
        profiles: list[np.ndarray],
        title: str,
        color: str,
    ) -> None:
        bin_edges   = np.linspace(0, 100, self.n_bins + 1)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        mat         = np.array(profiles)
        median      = np.nanmedian(mat, axis=0)
        q1          = np.nanpercentile(mat, 25, axis=0)
        q3          = np.nanpercentile(mat, 75, axis=0)

        for profile in profiles:
            ax.plot(bin_centers, profile,
                    color=self.trace_color,
                    linewidth=self.trace_linewidth,
                    alpha=self.trace_alpha)

        ax.fill_between(bin_centers, q1, q3,
                        color=color, alpha=self.ribbon_alpha,
                        label="IQR (25th–75th percentile)")
        ax.plot(bin_centers, median,
                color=color, linewidth=self.median_linewidth,
                linestyle="-", label="Median")

        ax.set_xlim(0, 100)
        ax.set_ylim(0, 1)
        ax.set_xticks([0, 25, 50, 75, 100])
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_ylabel("Gap frequency",
                      fontsize=self.label_fs, fontweight=self.label_fw)
        ax.set_title(f"{title} (n={len(profiles)}, bins={self.n_bins})",
                     fontsize=self.label_fs + 1, fontweight=self.label_fw)
        ax.tick_params(labelsize=self.tick_fs)
        for lbl in ax.get_xticklabels() + ax.get_yticklabels():
            lbl.set_fontweight(self.tick_fw)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(fontsize=self.tick_fs, frameon=True, loc="upper right")
