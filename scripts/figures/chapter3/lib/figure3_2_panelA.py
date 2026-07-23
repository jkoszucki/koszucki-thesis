"""
Figure 3.2, Panel X — GWAS prediction and O-acetylation overview table.

Table with 35 most-abundant K-loci as columns and five evidence rows:

  Row 0: Best GWAS predictor per locus (precision >= 0.80, best F1xMCC)
          coloured by ECOD class:
            sgnh-ecod  → deacetylase gold (strong shade)
            ssrbh-ecod → depolymerase green (strong shade)
            no-ecod    → gray
  Row 1: Predicted prophage depolymerases (PLOS Biology S3, 26 proteins)
          strong shade → ssrbh_color; likely shade → light green
  Row 2: Predicted prophage deacetylases (SGNH best per locus, 12 K-loci)
          strong (precision >= 0.80) → sgnh_domain_color; likely → light gold
  Row 3: K-locus acetyltransferase candidates (acetylases_kloci.tsv, 59 K-loci)
          present → sslbh_color (blue)
  Row 4: O-acetylation reported in CPS repeating unit (modifications sheet)
          kpam_reported YES → dark red; NO → light salmon

All filled cells show a count. Cells without data are white.

Outputs:
    plots_dir/figure3_2-panelX.png / .pdf
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import pandas as pd

mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

_ROW_LABELS = [
    "Best GWAS predictor",
    "Predicted depolymerase",
    "Predicted deacetylase",
    "K-locus acetyltransferase",
    "CPS O-acetylation",
]

# KL111 and KL30 best predictors are acetyltransferases identified via FoldSeek;
# not captured in acetylases_gwas.tsv (HHsearch-based) — manually assigned here.
_FOLDSEEK_ACETYLTRANSFERASES = {"KL111", "KL30"}

# ---------------------------------------------------------------------------
# Colors
# ---------------------------------------------------------------------------
_C_DEPOLY_STRONG  = "#a6d1a6"   # ssrbh_color
_C_DEPOLY_LIKELY  = "#d4ead4"
_C_DEAC_STRONG    = "#c9a227"   # sgnh_domain_color
_C_DEAC_LIKELY    = "#e8d88c"
_C_ACETYL         = "#4393c3"   # sslbh_color
_C_OTHER          = "#bfbfbf"   # gray_color
_C_OAC_KPAM      = "#d62728"   # kpam-reported O-acetylation
_C_OAC_NKPAM     = "#f4a582"   # non-kpam O-acetylation
_C_NA            = "#ffffff"   # n/a / 0 — white with label


def _locus_num(locus: str) -> int:
    """Return the integer ID from a locus string, e.g. 'KL103' -> 103."""
    return int(locus[2:])


def _ecod_to_color(ecod_folder: str) -> str:
    if ecod_folder == "sgnh-ecod-reported-topology":
        return _C_DEAC_STRONG
    if ecod_folder == "ssrbh-ecod-reported-topology":
        return _C_DEPOLY_STRONG
    return _C_OTHER


def _ktype_to_locus(ktype: str) -> str | None:
    """Convert 'K2' -> 'KL2'; returns None for ambiguous types like 'K21b', 'UNK'."""
    s = ktype.strip()
    if not s.startswith("K"):
        return None
    rest = s[1:]
    if not rest.isdigit():
        return None
    return f"KL{rest}"


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _load_best_predictor(
    best_predictors_tsv: Path,
    acetylases_gwas_tsv: Path,
) -> dict[str, tuple[str, int]]:
    """Returns {locus: (color, count=1)}.

    ECOD annotation is primary; if the best predictor PC also has an HHsearch
    acetyltransferase hit in acetylases_gwas (acetyl_hhsearch=True), override
    to acetyltransferase color regardless of ecod_folder.
    """
    df = pd.read_csv(best_predictors_tsv, sep="\t")
    ac = pd.read_csv(acetylases_gwas_tsv, sep="\t")
    ac_positives = set(zip(ac[ac["acetyl_hhsearch"] == True]["locus"],
                           ac[ac["acetyl_hhsearch"] == True]["PC"]))
    result = {}
    for _, row in df.iterrows():
        if row["locus"] in _FOLDSEEK_ACETYLTRANSFERASES:
            color = _C_ACETYL
        elif (row["locus"], row["PC"]) in ac_positives:
            color = _C_ACETYL
        else:
            color = _ecod_to_color(row["ecod_folder"])
        result[row["locus"]] = (color, 1)
    return result


def _load_depolymerases(depoly_tsv: Path) -> dict[str, list[tuple[str, int]]]:
    """Returns {locus: [(color, count), ...]}; split into strong/likely subcells when both present."""
    df = pd.read_csv(depoly_tsv, sep="\t")
    counts: dict[str, dict[str, int]] = {}
    for _, row in df.iterrows():
        locus    = row["locus"]
        strength = row["prediction_strength"]
        if locus not in counts:
            counts[locus] = {"strong": 0, "likely": 0}
        key = "strong" if strength == "strong" else "likely"
        counts[locus][key] += 1
    result = {}
    for locus, c in counts.items():
        subcells = []
        if c["strong"]:
            subcells.append((_C_DEPOLY_STRONG, c["strong"]))
        if c["likely"]:
            subcells.append((_C_DEPOLY_LIKELY, c["likely"]))
        result[locus] = subcells   # len==1: single cell; len==2: split top/bottom
    return result


def _load_deacetylases(deac_tsv: Path) -> dict[str, tuple[str, int]]:
    """Returns {locus: (color, count=1)}; strong/likely from precision >= 0.80."""
    df = pd.read_csv(deac_tsv, sep="\t")
    result = {}
    for _, row in df.iterrows():
        color = _C_DEAC_STRONG if row["precision"] >= 0.80 else _C_DEAC_LIKELY
        result[row["locus"]] = (color, 1)
    return result


def _load_acetylases_kloci(acetylases_kloci_tsv: Path) -> dict[str, tuple[str, int]]:
    """Returns {locus: (color, count)}."""
    df = pd.read_csv(acetylases_kloci_tsv, sep="\t")
    df["locus"] = df["protein_id"].str.extract(r"^(KL\d+)")
    counts = df.groupby("locus").size()
    return {locus: (_C_ACETYL, int(n)) for locus, n in counts.items()}


def _load_o_acetylation(cps_xlsx: Path) -> dict[str, tuple[str, int]]:
    """Returns {locus: (color, count)}; kpam hue.

    Only canonical KPAM-sourced structures are considered. LIT-only structures
    (e.g. CPS78_K2) share a serotype with a canonical K-locus but represent a
    distinct capsule variant and must not be attributed to that K-locus.
    """
    ktypes = pd.read_excel(cps_xlsx, sheet_name="capsule_structures")
    canonical_ids = set(ktypes.loc[ktypes["Source"] == "KPAM", "structure_id"])
    mods = pd.read_excel(cps_xlsx, sheet_name="capsule_modifications")
    oac = mods[
        (mods["modification"] == "acetylation") &
        (mods["structure_id"].isin(canonical_ids))
    ].copy()
    oac["locus"] = oac["K-type"].apply(_ktype_to_locus)
    oac = oac.dropna(subset=["locus"])
    result = {}
    for locus, grp in oac.groupby("locus"):
        has_kpam = (grp["kpam_reported"] == "YES").any()
        color = _C_OAC_KPAM if has_kpam else _C_OAC_NKPAM
        result[locus] = (color, len(grp))
    return result


def _load_acetyl_merged(
    acetylases_kloci_tsv: Path,
    cps_xlsx: Path,
    all_k_loci: list[str],
) -> dict[str, list[tuple[str, object]]]:
    """Returns {locus: [at_subcell, oac_subcell]} — always exactly two subcells.

    Left  — K-locus acetyltransferase: blue with count if present; gray, no label if absent.
    Right — CPS O-acetylation: red with count if documented; gray with 'n/a' for loci >= 99;
            gray, no label if absent (loci < 99).

    A locus is included when AT is present or (OAc is present and locus < 99).
    """
    at_data  = _load_acetylases_kloci(acetylases_kloci_tsv)
    oac_data = _load_o_acetylation(cps_xlsx)

    loci = (
        set(at_data)
        | {l for l in oac_data if _locus_num(l) < 99}
        | {l for l in all_k_loci if _locus_num(l) >= 99}
    )
    result = {}
    for locus in loci:
        if locus in at_data:
            at_subcell = at_data[locus]
        elif _locus_num(locus) >= 99:
            at_subcell = ("#ffffff", None)   # empty for loci >= 99 without AT
        else:
            at_subcell = (_C_NA, 0)          # gray with 0 for loci < 99 without AT
        if _locus_num(locus) >= 99:
            oac_subcell = ("#ffffff", "n/a")
        elif locus in oac_data:
            oac_subcell = oac_data[locus]
        else:
            oac_subcell = (_C_NA, 0)
        result[locus] = [at_subcell, oac_subcell]
    return result


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------

def _draw_cell(ax, x, y, w, h, color, count, cells_fontsize, number_color="#000000"):
    rect = patches.Rectangle((x, y), w, h, facecolor=color, edgecolor="black", linewidth=0.5)
    ax.add_patch(rect)
    if count is not None:
        ax.text(
            x + w / 2, y + h / 2, str(count),
            ha="center", va="center",
            color=number_color, fontsize=cells_fontsize, fontweight="bold",
        )


def _draw_split_cell(ax, x, y, w, h, subcells, cells_fontsize):
    """Draw a cell split into len(subcells) vertical bands; first entry on left."""
    n     = len(subcells)
    sub_w = w / n
    for i, (color, count) in enumerate(subcells):
        _draw_cell(ax, x + i * sub_w, y, sub_w, h, color, count,
                   cells_fontsize - 0.5)


def _draw_empty_cell(ax, x, y, w, h):
    rect = patches.Rectangle((x, y), w, h, facecolor="white", edgecolor="black", linewidth=0.5)
    ax.add_patch(rect)


def _draw_na_cell(ax, x, y, w, h, cells_fontsize):
    rect = patches.Rectangle((x, y), w, h, facecolor=_C_NA, edgecolor="black", linewidth=0.5)
    ax.add_patch(rect)
    ax.text(
        x + w / 2, y + h / 2, "n/a",
        ha="center", va="center",
        color="#000000", fontsize=cells_fontsize, fontweight="bold",
    )


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def plot_figure3_2_panelA(
    best_predictors_tsv: Path,
    depolymerases_gwas_tsv: Path,
    deacetylases_gwas_tsv: Path,
    acetylases_gwas_tsv: Path,
    acetylases_kloci_tsv: Path,
    cps_xlsx: Path,
    pyseer_hits_tsv: Path,
    plots_dir: Path,
    style=None,
    show_best_predictor: bool = False,
) -> None:
    """
    Produce Figure 3.2 Panel X: GWAS prediction overview table.

    Args:
        best_predictors_tsv:    rbp_best_predictors/best_predictors_gwas.tsv
        depolymerases_gwas_tsv: rbp_depolymerases/depolymerases_gwas.tsv
        deacetylases_gwas_tsv:  rbp_deacetylases/deacetylases_gwas.tsv
        acetylases_gwas_tsv:    cps_acetylases/acetylases_gwas.tsv
        acetylases_kloci_tsv:   cps_acetylases/acetylases_kloci.tsv
        cps_xlsx:               input/cps/cps.xlsx
        pyseer_hits_tsv:        input_dir/gwas/3_GWAS/3_PROCESSING/pyseer_hits_filtered.tsv; source of the 35-locus column universe
        plots_dir:              output directory for plots
        style:                  cfg.style (optional)
        show_best_predictor:    include the "Best GWAS predictor" row (default False)
    """
    dpi           = getattr(style, "dpi", 300)
    labels_fs     = 8
    cells_fs      = 5
    fig_w         = 7.5
    margin        = 0.3
    header_space  = 0.1

    pyseer = pd.read_csv(pyseer_hits_tsv, sep="\t")
    pyseer = pyseer[pyseer["mode"] == "lasso"]
    all_k_loci = sorted(pyseer["locus"].unique(), key=_locus_num)

    # --- build row list ---
    all_labels = [
        "Best GWAS predictor",
        "Predicted depolymerase",
        "Predicted deacetylase",
        "Acetyltransferase/acetylation",
    ]
    all_data = [
        _load_best_predictor(best_predictors_tsv, acetylases_gwas_tsv),
        _load_depolymerases(depolymerases_gwas_tsv),
        _load_deacetylases(deacetylases_gwas_tsv),
        _load_acetyl_merged(acetylases_kloci_tsv, cps_xlsx, all_k_loci),
    ]

    if not show_best_predictor:
        all_labels = all_labels[1:]
        all_data   = all_data[1:]

    row_labels = all_labels
    row_data   = all_data

    num_rows = len(row_labels)
    num_cols = len(all_k_loci)
    fig_h    = 0.3 + num_rows * 0.225   # scale height with row count
    cell_w   = (fig_w - 2 * margin) / num_cols
    cell_h   = (fig_h - 2 * margin) / num_rows

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(0, fig_w)
    ax.set_ylim(0, fig_h)
    ax.axis("off")

    # --- column headers ---
    header_y = fig_h - margin + header_space
    for i, locus in enumerate(all_k_loci):
        x = margin + i * cell_w + cell_w / 2
        ax.text(x, header_y, locus, ha="center", va="bottom",
                fontsize=labels_fs, fontweight="bold", rotation=60)

    # --- cells ---
    for row_idx, data in enumerate(row_data):
        for col_idx, locus in enumerate(all_k_loci):
            x = margin + col_idx * cell_w
            y = margin + (num_rows - 1 - row_idx) * cell_h
            if locus in data:
                val = data[locus]
                if isinstance(val, list):
                    if len(val) == 1:
                        _draw_cell(ax, x, y, cell_w, cell_h, val[0][0], val[0][1], cells_fs)
                    else:
                        _draw_split_cell(ax, x, y, cell_w, cell_h, val, cells_fs)
                else:
                    color, count = val
                    _draw_cell(ax, x, y, cell_w, cell_h, color, count, cells_fs)
            else:
                _draw_empty_cell(ax, x, y, cell_w, cell_h)

    # --- row labels (left side) ---
    for row_idx, label in enumerate(row_labels):
        y = margin + (num_rows - 1 - row_idx) * cell_h + cell_h / 2
        ax.text(margin - 0.05, y, label, ha="right", va="center",
                fontsize=labels_fs, fontweight="bold")

    plots_dir.mkdir(parents=True, exist_ok=True)
    plt.tight_layout(pad=0)
    for ext in ("png", "pdf"):
        out = plots_dir / f"figure3_2-panelA.{ext}"
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        print(f"  [figure3_2-panelA] → {out.name}")
    plt.close(fig)

    # -----------------------------------------------------------------------
    # Legend — color swatches, saved to legends/
    # -----------------------------------------------------------------------
    import matplotlib.patches as mpatches

    def _swatch(color, label, edgecolor="black"):
        return mpatches.Patch(facecolor=color, edgecolor=edgecolor, linewidth=0.6, label=label)

    def _header(text):
        import matplotlib.lines as mlines
        latex = text.replace(" ", "~")
        return mlines.Line2D([], [], color="none", label=f"$\\bf{{{latex}}}$")

    def _spacer():
        import matplotlib.lines as mlines
        return mlines.Line2D([], [], color="none", label=" ")

    handles = [
        _header("Predicted depolymerase"),
        _swatch(_C_DEPOLY_STRONG, "strong (precision ≥ 0.8)"),
        _swatch(_C_DEPOLY_LIKELY, "likely (precision < 0.8)"),
        _spacer(),
        _header("Predicted deacetylase"),
        _swatch(_C_DEAC_STRONG, "strong (precision ≥ 0.8)"),
        _swatch(_C_DEAC_LIKELY, "likely (precision < 0.8)"),
        _spacer(),
        _header("Acetyltransferase / acetylation"),
        _swatch(_C_ACETYL,   "K-locus acetyltransferase present"),
        _swatch(_C_OAC_KPAM, "O-acetylation (KPAM-reported)"),
        _swatch(_C_OAC_NKPAM, "O-acetylation (literature only)"),
        _swatch("#ffffff", "absent / n/a", edgecolor="#999999"),
    ]

    fig_leg, ax_leg = plt.subplots(figsize=(2.8, 4.0))
    ax_leg.axis("off")
    ax_leg.legend(
        handles=handles, fontsize=labels_fs,
        loc="center", ncol=1,
        frameon=False, handletextpad=0.5, labelspacing=0.4,
    )

    legends_dir = plots_dir / "legends"
    legends_dir.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    for ext in ("png", "pdf"):
        out = legends_dir / f"figure3_2-panelA-legend.{ext}"
        fig_leg.savefig(out, dpi=dpi, bbox_inches="tight")
        print(f"  [figure3_2-panelA-legend] → legends/{out.name}")
    plt.close(fig_leg)
