# koszucki-thesis

Complete workflow for PhD thesis analysis and figure reproducibility — phage receptor-binding
proteins and capsular polysaccharide (CPS) diversity in *Klebsiella pneumoniae*.

Covers three analysis threads:

- **GWAS-predicted phage receptor-binding proteins (RBPs)** — depolymerases (SSRBH fold) and
  putative deacetylases (SGNH fold), with experimental validation.
- **SGNH hydrolase diversity and acetyltransferase detection** — sequence/structure comparison
  (BLASTp, HHsearch, FoldSeek) of K-locus- and prophage-encoded acetyltransferases (SSLBH fold).
- **CPS structural diversity and O-acetylation** — pairwise structural similarity across 81
  resolved *Klebsiella* CPS repeating-unit structures, and their chemical modification patterns.

## Repository structure

```
.
├── config/
│   └── config.yml          ← paths, plot style, colour palette
└── scripts/
    ├── processing/         ← reads raw input, prepares data, writes to an external output directory
    │   ├── gwas-proc/          GWAS filtering, classification, per-PC export
    │   ├── sgnh-proc/          SGNH hydrolase predictor selection and organisation
    │   ├── cps-proc/           CPS K-type diversity analysis and structure diagrams
    │   ├── enzymes-proc/       experimental enzyme literature tables
    │   └── acetyl-proc/        acetyltransferase detection and evidence tables
    ├── figures/            ← reads prepared data, writes plots (this is the only output tracked in-repo)
    │   └── chapterX/
    │       ├── main.py         entry point with a run-flags block per panel
    │       ├── lib/             one module per figure panel
    │       ├── plots/           generated PNG/PDF + Cytoscape node/edge tables
    │       └── figures/         manually assembled panels (Adobe Illustrator sources)
    └── helpers/            ← shared Config/Style classes, structure rendering, AF3 indexing
```

Raw input data and full analysis output (large intermediate tables, AlphaFold3 structures, etc.)
live outside this repository and are not tracked here — only the code and the final figure
plots/panels are.

## Pipeline

```
gwas-proc → sgnh-proc   → sgnh-hydrolases/            → figures/chapter2 + figures/chapter3
gwas-proc               → gwas-data/                   → figures/chapter2 + figures/chapter3
           acetyl-proc  → acetyltransferase/            → figures/chapter3
cps-proc                → cps_structures/ + cps_drawn/  → figures/chapter4
enzymes-proc             → rbp_deacetylases/ + cps_acetylases/  (literature tables, shared)

  ↓ (all of the above)
helpers/build_af3_index.py → indexes every canonical protein set against its AlphaFold3 structure
```

Each processing module is checkpointed — safe to rerun; only new or missing outputs are
recomputed.

| Chapter | Topic | Processing modules |
|---------|-------|--------------------|
| 2 | GWAS depolymerases (SSRBH) and putative deacetylases (SGNH); experimental validation | `gwas-proc` → `sgnh-proc` |
| 3 | SGNH hydrolase diversity + acetyltransferase detection (SSLBH, HHsearch/FoldSeek) | `gwas-proc` + `sgnh-proc` + `acetyl-proc` |
| 4 | CPS K-type diversity and O-acetylation | `cps-proc` |

## Environment setup

The primary environment (`jkoszucki`) is exported to `env.yaml`:

```bash
conda env create -f env.yaml
conda activate jkoszucki
```

Structure rendering (PyMOL), K-locus extraction (Kaptive), and structural alignment (TM-align)
each require their own separate environment — `pymol`, `kaptive`, and `tmtools` respectively —
not included in `env.yaml`.

## Running

All scripts run under the `jkoszucki` conda environment unless noted otherwise (PyMOL structure
rendering requires the `pymol` environment; Kaptive requires `kaptive`; TM-align requires
`tmtools`).

```bash
# Full processing pipeline
conda run -n jkoszucki python scripts/run_pipeline.py

# A single processing module
conda run -n jkoszucki python scripts/processing/gwas-proc/main.py

# Figure generation, per chapter
conda run -n jkoszucki python scripts/figures/chapter2/main.py
conda run -n jkoszucki python scripts/figures/chapter3/main.py
conda run -n jkoszucki python scripts/figures/chapter4/main.py
```

Each `main.py` has a run-flags block at the top (plain booleans) to toggle individual panels
without touching `config.yml`.

## Configuration

`config/config.yml` holds all paths (input/output directories, GWAS data root) and shared plot
style (fonts, DPI, ECOD topology colours). `Config` and `Style` (in `scripts/helpers/config.py`)
are loaded once per script and passed down to every `lib/` module — nothing under `lib/` reads
the config file directly.

## Naming convention

Figure artefacts are named from the manuscript figure number, prefixed with the chapter number
(e.g. manuscript Figure 2.3 panel A → `figure2_3_panelA`):

| Artefact | Pattern | Example |
|----------|---------|---------|
| lib script | `lib/figure{C}_{F}_panel{P}.py` | `lib/figure2_3_panelA.py` |
| entry-point function | `plot_figure{C}_{F}_panel{P}` | `plot_figure2_3_panelA` |
| plot output | `plots/figure{C}_{F}-panel{P}.png` / `.pdf` | `plots/figure2_3-panelA.png` |

## Citation

This repository accompanies a PhD thesis building on Otwinowska, Koszucki et al. 2026, *PLOS
Biology* ("Capsular specificity in temperate phages of *Klebsiella pneumoniae* is driven by
diverse receptor-binding enzymes").
