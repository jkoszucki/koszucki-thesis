"""
Build sgnh-hydrolases/n-terminal/rbps/ folder.

Contents:
  rbps.fasta         — all RBP sequences (from prophage protein DB)
  af3_trimers.json   — homotrimer AF3 upload jobs for all proteins
  structures_cif/    — symlinks: {PROTEIN_ID}.cif → model_0 CIF from AF3 raw output
  structures_png/    — rendered PNGs (res 1–130 orange) for N-terminal-sharing proteins only

Protein lists are hardcoded below. Checkpoints per output file/symlink.
"""

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

# ---------------------------------------------------------------------------
# Protein lists — hardcoded
# ---------------------------------------------------------------------------

ALL_PROTEINS = [
    "KPN_B2_PHAGE002_M_PROTEIN_55",
    "KPN_B2_PHAGE002_M_PROTEIN_58",
    "KVV_B3_PHAGE169_L_PROTEIN_71",
    "KVV_B3_PHAGE169_L_PROTEIN_72",
    "KVV_B1_PHAGE242_M_PROTEIN_14",
    "KVV_B1_PHAGE242_M_PROTEIN_13",
    "KVV_B1_PHAGE242_M_PROTEIN_12",
    "KVV_B1_PHAGE242_M_PROTEIN_11",
    "KPN_B18_PHAGE103_M_PROTEIN_13",
    "KPN_B18_PHAGE103_M_PROTEIN_12",
    "KPN_B18_PHAGE103_M_PROTEIN_11",
    "KPN_B18_PHAGE103_M_PROTEIN_10",
]

N_TERMINAL_PROTEINS = [
    "KPN_B2_PHAGE002_M_PROTEIN_55",
    "KPN_B2_PHAGE002_M_PROTEIN_58",
    "KVV_B3_PHAGE169_L_PROTEIN_71",
    "KVV_B3_PHAGE169_L_PROTEIN_72",
    "KVV_B1_PHAGE242_M_PROTEIN_14",
    "KVV_B1_PHAGE242_M_PROTEIN_11",
    "KPN_B18_PHAGE103_M_PROTEIN_13",
]

# Highlighted region
NTERM_RES_START = 1
NTERM_RES_END   = 130
NTERM_COLOR     = "orange"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _find_cif(protein: str, af3_raw_dir: Path) -> Path | None:
    """Walk af3_raw_dir recursively for fold_{protein_lower}_model_0.cif."""
    target = f"fold_{protein.lower()}_model_0.cif"
    for cif in af3_raw_dir.rglob(target):
        return cif
    return None


def _render_gray(cif: Path, out_png: Path) -> bool:
    with tempfile.NamedTemporaryFile(suffix=".pml", mode="w", delete=False) as f:
        f.write(f"""\
reinitialize
load {cif}, protein
hide everything
show cartoon, protein
color gray70, protein
set ray_opaque_background, 1
bg_color white
set ray_trace_mode, 1
reset
orient
zoom all, 40
ray 4500, 3000
png {out_png}, dpi=300
quit
""")
        pml = Path(f.name)
    r = subprocess.run(
        ["conda", "run", "-n", "pymol", "pymol", "-c", str(pml)],
        capture_output=True, text=True,
    )
    pml.unlink(missing_ok=True)
    return r.returncode == 0


def _render_highlighted(cif: Path, out_png: Path,
                        res_start: int, res_end: int, color: str) -> bool:
    with tempfile.NamedTemporaryFile(suffix=".pml", mode="w", delete=False) as f:
        f.write(f"""\
reinitialize
load {cif}, protein
hide everything
show cartoon, protein
color gray70, protein
color {color}, protein and resi {res_start}-{res_end}
set ray_opaque_background, 1
bg_color white
set ray_trace_mode, 1
reset
orient
zoom all, 40
ray 4500, 3000
png {out_png}, dpi=300
quit
""")
        pml = Path(f.name)
    r = subprocess.run(
        ["conda", "run", "-n", "pymol", "pymol", "-c", str(pml)],
        capture_output=True, text=True,
    )
    pml.unlink(missing_ok=True)
    return r.returncode == 0


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def build_rbp_folder(
    prophage_faa_pattern: str,
    blast_db_dir: Path,
    af3_raw_dir: Path,
    n_terminal_dir: Path,
    render_structures: bool = True,
) -> None:
    """
    Build other/n-terminal/prophages_rbps/.

    Args:
        prophage_faa_pattern: glob pattern for prophage FASTAs (seq source)
        blast_db_dir:         prophage BLAST DB dir (used by prepare_prophage_db)
        af3_raw_dir:          input_dir/supplementary-thesis/supplementary-data/S4_Data/1_INPUT/1_RAW_ALPHAFOLD3/ (read-only)
        n_terminal_dir:       output_dir/other/n-terminal/
        render_structures:    if False, skip PNG rendering
    """
    from prophage_db import prepare_prophage_db

    out_dir = n_terminal_dir / "prophages_rbps"
    out_dir.mkdir(parents=True, exist_ok=True)

    # -- sequences --
    _, raw_seq_index = prepare_prophage_db(prophage_faa_pattern, blast_db_dir)
    # Strip X (unknown aa) — present in some prophage DB entries
    seq_index = {k: v.replace("X", "") for k, v in raw_seq_index.items()}

    # 1. rbps.fasta
    fasta_path = out_dir / "rbps.fasta"
    if fasta_path.exists():
        print(f"  rbps.fasta exists — skipping (delete to recompute)")
    else:
        missing = []
        with open(fasta_path, "w") as fh:
            for protein in ALL_PROTEINS:
                seq = seq_index.get(protein)
                if not seq:
                    missing.append(protein)
                    continue
                fh.write(f">{protein}\n{seq}\n")
        if missing:
            print(f"  [warn] sequences not found: {missing}")
        written = len(ALL_PROTEINS) - len(missing)
        print(f"  rbps.fasta — {written}/{len(ALL_PROTEINS)} sequences written")

    # 2. af3_trimers.json — always regenerated; excludes proteins with CIF already downloaded
    af3_json = out_dir / "af3_trimers.json"
    cif_dir  = out_dir / "structures_cif"
    jobs = []
    done = []
    for protein in ALL_PROTEINS:
        if (cif_dir / f"{protein}.cif").exists():
            done.append(protein)
            continue
        seq = seq_index.get(protein)
        if not seq:
            continue
        jobs.append({
            "name": protein.lower(),
            "modelSeeds": [],
            "sequences": [
                {"proteinChain": {"sequence": seq, "count": 3}}
            ],
            "dialect": "alphafoldserver",
            "version": 1,
        })
    af3_json.write_text(json.dumps(jobs, indent=2))
    print(f"  af3_trimers.json — {len(jobs)} pending, {len(done)} already downloaded")

    # 3. structures_cif/ symlinks
    cif_dir = out_dir / "structures_cif"
    cif_dir.mkdir(exist_ok=True)
    ok = skipped = missing_cif = 0
    for protein in ALL_PROTEINS:
        dst = cif_dir / f"{protein}.cif"
        if dst.exists() or dst.is_symlink():
            skipped += 1
            continue
        src = _find_cif(protein, af3_raw_dir)
        if src is None:
            missing_cif += 1
            continue
        dst.symlink_to(src)
        ok += 1
    print(f"  structures_cif/ — {ok} symlinked, {skipped} skipped, {missing_cif} missing")

    # 4. structures_png/ — N-terminal-sharing proteins only
    if not render_structures:
        print(f"  structures_png/ — skipped (render_structures=False)")
        return

    png_dir = out_dir / "structures_png"
    png_dir.mkdir(exist_ok=True)
    ok = skipped = errors = missing_cif = 0

    # N-terminal-sharing: orange res 1–130, rest gray
    for protein in N_TERMINAL_PROTEINS:
        out_png = png_dir / f"{protein}.png"
        if out_png.exists():
            skipped += 1
            continue
        cif = cif_dir / f"{protein}.cif"
        if not cif.exists():
            missing_cif += 1
            continue
        success = _render_highlighted(cif, out_png, NTERM_RES_START, NTERM_RES_END, NTERM_COLOR)
        if success:
            print(f"    → {out_png.name}")
            ok += 1
        else:
            print(f"    [render error] {protein}")
            errors += 1

    # Remaining proteins: plain gray
    remaining = [p for p in ALL_PROTEINS if p not in N_TERMINAL_PROTEINS]
    for protein in remaining:
        out_png = png_dir / f"{protein}.png"
        if out_png.exists():
            skipped += 1
            continue
        cif = cif_dir / f"{protein}.cif"
        if not cif.exists():
            missing_cif += 1
            continue
        success = _render_gray(cif, out_png)
        if success:
            print(f"    → {out_png.name}")
            ok += 1
        else:
            print(f"    [render error] {protein}")
            errors += 1

    print(f"  structures_png/ — {ok} rendered, {skipped} skipped, "
          f"{missing_cif} missing cif, {errors} errors")
