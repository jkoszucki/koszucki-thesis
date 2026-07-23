"""
Render a folder of CIF structures coloured by protein source (virulent / prophage).

Usage:
    conda run -n pymol python render_cif_by_source.py \
        --cif_dir <path> \
        --source_tsv <active_enzymes.tsv> \
        [--virulent_color #b794c5] \
        [--prophage_color #e968a3] \
        [--out_dir <path>]

Rendering parameters match draw_structure.py (supplementary_text_S3):
    ray_trace_mode=1, ray_trace_gain=0.1, ambient=0.6, specular=0.2,
    shininess=10, reflect=0.0, antialias=2, viewport 1024×768, ray 3000×2400, dpi=300.

Background is transparent (ray_opaque_background=0).
Output PNGs are written next to each CIF, or to --out_dir if specified.
"""

import argparse
from pathlib import Path

from pymol import cmd
import pandas as pd


DEFAULT_VIRULENT = "#b794c5"
DEFAULT_PROPHAGE = "#e968a3"


def hex_to_rgb(h: str) -> list[float]:
    h = h.lstrip("#")
    return [int(h[i:i+2], 16) / 255.0 for i in (0, 2, 4)]


def render_cif_by_source(
    cif_dir: Path,
    source_tsv: Path | None = None,
    virulent_color: str = DEFAULT_VIRULENT,
    prophage_color: str = DEFAULT_PROPHAGE,
    out_dir: Path | None = None,
) -> None:
    cif_dir = Path(cif_dir)
    out_dir = Path(out_dir) if out_dir else cif_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # Build protein → source map from TSV if provided
    source_map: dict[str, str] = {}
    if source_tsv:
        df = pd.read_csv(source_tsv, sep="\t")
        source_map = dict(zip(df["proteinID"], df["source"]))

    cmd.set_color("virulent_col", hex_to_rgb(virulent_color))
    cmd.set_color("prophage_col", hex_to_rgb(prophage_color))

    for cif in sorted(cif_dir.glob("*.cif")):
        # Match protein ID: try full stem, then stem without KL suffix
        stem = cif.stem
        source = source_map.get(stem)
        if source is None:
            # try stripping trailing _KLxx / _KNxx
            base = "_".join(stem.split("_")[:-1]) if "_" in stem else stem
            source = source_map.get(base, "virulent")

        color = "virulent_col" if source == "virulent" else "prophage_col"
        out_png = out_dir / f"{stem}.png"

        print(f"Rendering {stem} ({source}) ...")

        cmd.load(str(cif), stem)
        cmd.bg_color("white")
        cmd.hide("everything", stem)
        cmd.show("cartoon", stem)
        cmd.color(color, stem)

        cmd.set("ray_trace_mode", 1)
        cmd.set("ray_trace_gain", 0.1)
        cmd.set("ambient",        0.6)
        cmd.set("specular",       0.2)
        cmd.set("shininess",      10)
        cmd.set("reflect",        0.0)
        cmd.set("antialias",      2)
        cmd.set("ray_opaque_background", 0)

        cmd.orient(stem)
        cmd.center(stem)
        cmd.zoom(stem, buffer=50)

        cmd.viewport(1024, 768)
        cmd.ray(3000, 2400)
        cmd.png(str(out_png), dpi=300)
        cmd.delete(stem)

        print(f"  → {out_png.name}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cif_dir",        required=True)
    parser.add_argument("--source_tsv",     default=None)
    parser.add_argument("--virulent_color", default=DEFAULT_VIRULENT)
    parser.add_argument("--prophage_color", default=DEFAULT_PROPHAGE)
    parser.add_argument("--out_dir",        default=None)
    args = parser.parse_args()

    render_cif_by_source(
        cif_dir        = args.cif_dir,
        source_tsv     = args.source_tsv,
        virulent_color = args.virulent_color,
        prophage_color = args.prophage_color,
        out_dir        = args.out_dir,
    )
