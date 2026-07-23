# koszucki-thesis

Complete workflow for PhD thesis analysis and figure reproducibility — phage receptor-binding
proteins and capsular polysaccharide (CPS) diversity in *Klebsiella pneumoniae*.

Covers three analysis threads:

- **GWAS-predicted phage receptor-binding proteins (RBPs)** — depolymerases and
  putative deacetylases
- **putative prophage deacetylases and K-locus acetyltrasferase diversity** — sequence and structure comparison
- **capsule structural diversity and O-acetylation** — pairwise structural similarity across 81
  NMR-resolved *Klebsiella* capsule repeating-unit structures along with their chemical modification.

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


| Chapter | Topic | Processing modules |
|---------|-------|--------------------|
| 2 | GWAS depolymerases and putative deacetylases; experimental proteins | `gwas-proc` → `sgnh-proc` |
| 3 | SGNH hydrolase diversity + acetyltransferase detection | `gwas-proc` + `sgnh-proc` + `acetyl-proc` |
| 4 | Capsule structures diversity and acetylation | `cps-proc` |

## Environment setup

The primary environment (`jkoszucki`) is exported to `env.yaml`:

```bash
conda env create -f env.yaml
conda activate jkoszucki
```

Structure rendering (PyMOL), K-locus extraction (Kaptive), and structural alignment (TM-align)
each require their own separate environment — `pymol`, `kaptive`, and `tmtools` respectively —
not included in `env.yaml`.

### Data Availability

1. Download the raw input data and supplementary material from Figshare: `https://figshare.com/s/abd6ae6a4af427f66a70`.
2. Unpack it locally, then point `config/config.yml` at it:

```yaml
paths:
  input_dir:  /path/to/downloaded/input
  output_dir: /path/to/analysis/output
  gwas_path:  /path/to/downloaded/input/data-gwas
```

`input_dir` is read-only throughout the pipeline; `output_dir` is written to by every processing
module and may be empty on first run.

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


## AI Statement

This codebase was written with the help of Claude Code (Anthropic).

## Citation

This repository accompanies a PhD thesis building on Otwinowska, Koszucki et al. 2026, *PLOS
Biology*.
