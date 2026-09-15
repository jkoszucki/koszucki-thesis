"""
Render the AlphaFold3 models of the figure 3.2B reference acetyltransferases.

The S1_Table reference proteins are called out on the panel's x-axis by label colour; their
structures are placed alongside the panel during manual assembly, so the renders use the
same two colours as the labels rather than the usual pLDDT scale:

    active_AT_*   experimentally characterised enzymes (K1, K2, K57)
                                                        -> _ACTIVE_COLOR    (label red)
    gwas_AT_*     putative GWAS-predicted enzymes (KL111)
                                                        -> _KLOCUS_AT_COLOR (blue)

Everything is rendered from the **monomer** AF3 model, never the homotrimer that KL30 and
KL111 also have under their best-predictor names.

Red marks an enzyme with experimental evidence. KL111 is rendered in the **blue** used for
K-locus acetyltransferase candidates throughout the chapter (`sslbh_color` / the
`acetylases_kloci` AF3 set colour), because it is a prediction like they are, not a
characterised enzyme (user, 2026-07-28). It previously used a pale red, "same hue as the
characterised enzymes but lighter".

The pale red `_PREDICTION_LABEL_COLOR` is gone from both panels as of 2026-07-28: the
x-tick labels in `figures/chapter4/lib/figureS11_panelA.py` and
`figureS11_panelB.py` now carry one distinction only, red = active / black =
predicted, so KL111's label is black there while its structure is rendered blue here. Both
say "prediction", in the vocabulary each figure has available — a label colour cannot use
the K-locus-candidate blue without clashing with the panel B bar categories.

Last validated against the manuscript caption on 2026-07-28 (user-confirmed): panel B reads
"three experimentally characterised enzymes (red) and the GWAS-predicted enzyme of KL111
(black)", panel D "red for the three experimentally characterised enzymes (K1, K2, K57),
blue for the predictions". The earlier "pale red" wording is superseded.

Rendering runs through `helpers/structure_renderer.StructureRenderer` in `color_mode="gray"`,
which applies one uniform colour to the whole model; PyMOL is invoked internally as a
subprocess, so this module is called from `jkoszucki` like any other. Always `accurate`
mode — these are thesis figure outputs.

Checkpointed: a protein whose PNG already exists is skipped, so a routine acetyl-proc rerun
does not re-render (accurate mode costs ~2 min per structure).

Output: other/alphafold3/1_DRAWN/{protein_id_lowercase}_model_0.png — the AF3 render
tree, so the figure colours are what the index and its organised view show.
"""

from __future__ import annotations

from pathlib import Path

from structure_renderer import StructureRenderer

_ACTIVE_COLOR    = "#b2182b"   # characterised enzymes — matches the panel's label red
_KLOCUS_AT_COLOR = "#4393c3"   # sslbh_color — the blue used for K-locus AT candidates

_MODEL_SUFFIX = "_model_0.cif"

# protein_id -> colour. Grouped by what the colour means, not by protein number.
#
# All four are rendered from the **monomer** models in S4_Data (verified single-chain:
# 170/221/326/360 residues). For KL111 that matters — it also has a homotrimer model filed
# under its best_predictors_gwas name (`PCI80C80_PC0915_KL111_model_0.cif`), which is not
# what this render uses.
#
# PROTEIN04_GWAS_AC_K30 is deliberately absent. KL30 is out of figure 3.2 entirely — it is
# a false-positive acetyltransferase call (no FoldSeek homologue beyond itself, no TM-score
# above 0.35 against any of the other 73 proteins) and appears only in the methods text.
# Its AF3 model and S3_Data FoldSeek table remain as deposited data, and it stays indexed
# under `acetylases_gwas_putative`; it simply gets no figure-coloured render.
REFERENCE_STRUCTURE_COLORS = {
    "PROTEIN01_MOD_AC_K1":    _ACTIVE_COLOR,
    "PROTEIN02_MOD_AC_K2":    _ACTIVE_COLOR,
    "PROTEIN03_MOD_AC_K57":   _ACTIVE_COLOR,
    "PROTEIN05_GWAS_AC_K111": _KLOCUS_AT_COLOR,
}


def render_reference_structures(
    s4_data_dir: Path,
    out_dir: Path,
    style=None,
    render_mode: str = "accurate",
) -> None:
    """
    Render each reference acetyltransferase in its assigned colour.

    Args:
        s4_data_dir: input_dir/.../S4_Data (read-only AF3 models)
        out_dir:     other/alphafold3/1_DRAWN — the AF3 render tree. Writing straight
                     there (rather than to a private folder and copying) keeps one copy
                     of each render: `build_af3_index.py` Step F builds
                     2_DRAWN_AND_ORGANISED from it, and Step G skips any protein that
                     already has a PNG, so these colours survive a reindex instead of
                     being replaced by the per-set colour `render_af3.py` would assign.
        style:       cfg.style (supplies dpi)
        render_mode: 'accurate' for thesis output; 'fast' only for quick checks
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = _MODEL_SUFFIX.replace(".cif", "")

    # One renderer per colour — colour is set on the renderer, not per structure.
    by_color: dict[str, list[Path]] = {}
    missing, skipped = [], []

    for protein_id, color in REFERENCE_STRUCTURE_COLORS.items():
        if (out_dir / f"{protein_id.lower()}{stem}.png").exists():
            skipped.append(protein_id)
            continue
        cif = next(
            (p for p in (s4_data_dir / f"{c}{_MODEL_SUFFIX}" for c in (protein_id, protein_id.lower()))
             if p.is_file()),
            None,
        )
        if cif is None:
            missing.append(protein_id)
            continue
        by_color.setdefault(color, []).append(cif)

    if skipped:
        print(f"  {len(skipped)} already rendered, skipping: {', '.join(skipped)}")
    if missing:
        print(f"  [warn] no AF3 model for: {', '.join(missing)}")
    if not by_color:
        print("  Nothing to render.")
        return

    for color, cifs in by_color.items():
        print(f"  rendering {len(cifs)} structure(s) in {color} …")
        renderer = StructureRenderer(
            render_mode = render_mode,
            color_mode  = "gray",      # uniform fill, colour supplied below
            gray_color  = color,
            style       = style,
        )
        renderer.render_to_dir(cifs, out_dir)

    # render_to_dir names PNGs after the CIF stem. S4_Data files are lowercase, which is
    # already the AF3 canonical name, but normalise anything filed under the upper-case
    # protein_id so the tree is consistent.
    for protein_id in REFERENCE_STRUCTURE_COLORS:
        src = out_dir / f"{protein_id}{stem}.png"
        if src.exists():
            src.rename(out_dir / f"{protein_id.lower()}{stem}.png")
    print(f"  → {out_dir}")
