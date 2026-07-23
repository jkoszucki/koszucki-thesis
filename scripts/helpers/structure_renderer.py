"""
Shared AF3 structure renderer — importable by any chapter script.

This file serves two roles:

  1. Imported as a library (jkoszucki env)
     Provides StructureRenderer, which builds job lists, checks for existing
     PNGs, and calls this same file as a subprocess in the pymol env.

  2. Run directly as a PyMOL worker (pymol env)
     Receives rendering parameters as JSON CLI args, renders each CIF.
     PyMOL is imported only inside the worker functions — never at module
     level — so importing this file in the jkoszucki env works without PyMOL.

Usage (in any figures/chapterX/main.py)
----------------------------------------
    sys.path.insert(0, str(Path(__file__).resolve().parents[N] / "helpers"))
    from structure_renderer import StructureRenderer

    renderer = StructureRenderer(
        render_mode   = RENDER_MODE,    # "fast" | "accurate"
        color_mode    = "plddt",        # "plddt" | "domain" | "chain" | "gray"
        domain_colors = [               # only used with color_mode="domain"
            ("1-350",   "#c9a227"),     # (resi_selection, hex_color)
            ("400-500", "#1f77b4"),
        ],
        gray_color    = "#cfcfce",      # base color for unlabelled residues
        background    = "white",
        width         = 3000,
        height        = 2400,
        style         = cfg.style,      # provides dpi
    )

    # Option A — write structure.png next to each CIF (standard for data hierarchy)
    renderer.render_next_to(cif_paths)

    # Option B — write {stem}.png into a named output directory (for plots/)
    renderer.render_to_dir(cif_paths, out_dir=plots_dir / "structures")

Color modes
-----------
  "plddt"   pLDDT (b-factor), standard AlphaFold scheme (default):
                > 90  dark blue   — very high confidence
                > 70  light blue  — confident
                > 50  yellow      — low confidence
                ≤ 50  orange      — very low confidence
  "domain"  Uniform gray_color base; named residue ranges coloured by domain_colors.
  "chain"   Each chain coloured distinctly (PyMOL util.cbc palette).
  "gray"    Uniform gray_color throughout.

Render modes
------------
  "fast"      No ray tracing — fast preview (default).
  "accurate"  Ray tracing enabled — publication quality (slow).
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_DEFAULT_GRAY    = "#cfcfce"
_DEFAULT_BG      = "white"
_DEFAULT_WIDTH   = 3000
_DEFAULT_HEIGHT  = 2400

_PLDDT_COLORS: list[tuple[int, list[float]]] = [
    (90, [  0/255,  83/255, 214/255]),   # very high  > 90  dark blue
    (70, [101/255, 203/255, 243/255]),   # confident  > 70  light blue
    (50, [255/255, 219/255,  19/255]),   # low        > 50  yellow
    ( 0, [255/255, 125/255,  69/255]),   # very low  <= 50  orange
]

_COLOR_MODES   = ("plddt", "domain", "chain", "gray")
_RENDER_MODES  = ("fast", "accurate")


# ---------------------------------------------------------------------------
# PyMOL worker internals — only executed when run as __main__
# ---------------------------------------------------------------------------

def _hex_to_rgb(hex_color: str) -> list[float]:
    h = hex_color.lstrip("#")
    return [int(h[i:i+2], 16) / 255.0 for i in (0, 2, 4)]


def _apply_plddt(obj: str) -> None:
    from pymol import cmd  # noqa: PLC0415 — pymol env only
    for i, (threshold, rgb) in enumerate(reversed(_PLDDT_COLORS)):
        name = f"_plddt_{i}"
        cmd.set_color(name, rgb)
        if threshold == 0:
            cmd.color(name, obj)
        else:
            cmd.color(name, f"({obj}) and b > {threshold}")


def _apply_domain(obj: str, gray_color: str, domain_colors: list[list]) -> None:
    from pymol import cmd  # noqa: PLC0415 — pymol env only
    cmd.set_color("_gray_base", _hex_to_rgb(gray_color))
    cmd.color("_gray_base", obj)
    for idx, (resi_sel, hex_color) in enumerate(domain_colors):
        color_name = f"_domain_{idx}"
        sel_name   = f"_sel_{idx}"
        cmd.set_color(color_name, _hex_to_rgb(hex_color))
        cmd.select(sel_name, f"({obj}) and resi {resi_sel}")
        cmd.color(color_name, sel_name)
        cmd.delete(sel_name)


def _apply_chain(obj: str) -> None:
    from pymol import cmd, util  # noqa: PLC0415 — pymol env only
    util.cbc(obj)


def _apply_gray(obj: str, gray_color: str) -> None:
    from pymol import cmd  # noqa: PLC0415 — pymol env only
    cmd.set_color("_gray_uniform", _hex_to_rgb(gray_color))
    cmd.color("_gray_uniform", obj)


def _mat3_mul(A: list, B: list) -> list:
    """Multiply two 3x3 matrices given as flat 9-element row-major lists."""
    C = [0.0] * 9
    for i in range(3):
        for k in range(3):
            for j in range(3):
                C[i*3+j] += A[i*3+k] * B[k*3+j]
    return C


def _centroid(pts: list) -> tuple:
    n = len(pts)
    return (sum(p[0] for p in pts) / n,
            sum(p[1] for p in pts) / n,
            sum(p[2] for p in pts) / n)


def _orient_nterm_to_top(obj: str, window: int = 200) -> None:
    """Rotate the view in-plane so N-terminal is at top, C-terminal at bottom.

    Uses the centroids of the first and last `window` CA atoms as stable
    reference points for the N→C axis, avoiding noise from disordered termini.
    Uses cmd.turn("z", ...) which rotates around the camera/screen Z axis —
    no depth-axis tilt. A self-correction step flips 180° if the convention
    sign is opposite.
    """
    from pymol import cmd  # noqa: PLC0415

    # Collect all chain A CA atoms, sorted by residue sequence number
    space: dict = {"atoms": []}
    cmd.iterate_state(1, f"({obj}) and chain A and name CA",
                      "atoms.append((int(resi), x, y, z))", space=space)
    if not space["atoms"]:
        cmd.iterate_state(1, f"({obj}) and name CA",
                          "atoms.append((int(resi), x, y, z))", space=space)
    if len(space["atoms"]) < 2:
        return

    atoms = sorted(space["atoms"], key=lambda a: a[0])
    coords = [(a[1], a[2], a[3]) for a in atoms]

    # Window size: at most `window`, but never more than half the chain
    k = min(window, len(coords) // 2)
    if k < 1:
        k = 1

    nterm_ref = _centroid(coords[:k])   # centroid of first 200 CAs
    cterm_ref = _centroid(coords[-k:])  # centroid of last  200 CAs

    # N→C direction vector in world space
    dx = cterm_ref[0] - nterm_ref[0]
    dy = cterm_ref[1] - nterm_ref[1]
    dz = cterm_ref[2] - nterm_ref[2]

    # Project N→C onto screen axes using the current view rotation
    r = cmd.get_view()[:9]
    sx = r[0]*dx + r[1]*dy + r[2]*dz   # screen X of N→C
    sy = r[3]*dx + r[4]*dy + r[5]*dz   # screen Y of N→C

    # Rotate so C→N = (-sx, -sy) points to +Y (top)
    angle = math.degrees(math.atan2(-sx, -sy))

    # cmd.turn("z", ...) rotates around the camera/screen Z — pure in-plane
    cmd.turn("z", angle)

    # Self-correction: verify N-terminal centroid ended up above the COM
    com = cmd.centerofmass(obj)
    r2 = cmd.get_view()[:9]
    vx = nterm_ref[0] - com[0]
    vy = nterm_ref[1] - com[1]
    vz = nterm_ref[2] - com[2]
    nterm_screen_y = r2[3]*vx + r2[4]*vy + r2[5]*vz
    if nterm_screen_y < 0:
        cmd.turn("z", 180)


def _render_one(
    cif_path: str,
    png_path: str,
    render_mode: str,
    color_mode: str,
    domain_colors: list[list],
    gray_color: str,
    background: str,
    width: int,
    height: int,
    dpi: int,
    orient_nterm: bool = False,
) -> None:
    from pymol import cmd  # noqa: PLC0415 — pymol env only

    os.makedirs(os.path.dirname(png_path), exist_ok=True)

    obj = os.path.splitext(os.path.basename(cif_path))[0]
    cmd.reinitialize()
    cmd.load(cif_path, obj)
    cmd.hide("everything", obj)
    cmd.show("cartoon", obj)

    if color_mode == "plddt":
        _apply_plddt(obj)
    elif color_mode == "domain":
        _apply_domain(obj, gray_color, domain_colors)
    elif color_mode == "chain":
        _apply_chain(obj)
    elif color_mode == "gray":
        _apply_gray(obj, gray_color)

    cmd.bg_color(background)
    cmd.set("antialias", 2)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_opaque_background", 0)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_smooth_loops", 1)
    cmd.viewport(width, height)
    cmd.orient(obj)
    cmd.center(obj)

    if orient_nterm:
        _orient_nterm_to_top(obj)

    cmd.zoom(obj, buffer=50)

    ray = 1 if render_mode == "accurate" else 0
    cmd.png(png_path, width=width, height=height, dpi=dpi, ray=ray)
    cmd.delete(obj)


def _worker_main() -> None:
    """Entry point when run as a subprocess in the pymol conda env."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--jobs",          required=True,
                        help='JSON array of {"cif":"...", "out":"..."} objects')
    parser.add_argument("--render-mode",   default="fast",   choices=list(_RENDER_MODES))
    parser.add_argument("--color-mode",    default="plddt",  choices=list(_COLOR_MODES))
    parser.add_argument("--domain-colors", default="[]",
                        help='JSON array of [resi_selection, hex_color] pairs')
    parser.add_argument("--gray-color",    default=_DEFAULT_GRAY)
    parser.add_argument("--background",    default=_DEFAULT_BG)
    parser.add_argument("--width",         type=int, default=_DEFAULT_WIDTH)
    parser.add_argument("--height",        type=int, default=_DEFAULT_HEIGHT)
    parser.add_argument("--dpi",           type=int, default=300)
    parser.add_argument("--orient-nterm",  action="store_true", default=False)
    args = parser.parse_args()

    jobs          = json.loads(args.jobs)
    domain_colors = json.loads(args.domain_colors)
    ok, failed = 0, []
    total      = len(jobs)

    for i, job in enumerate(jobs, 1):
        cif_name = os.path.basename(job["cif"])
        try:
            _render_one(
                cif_path      = job["cif"],
                png_path      = job["out"],
                render_mode   = args.render_mode,
                color_mode    = args.color_mode,
                domain_colors = domain_colors,
                gray_color    = args.gray_color,
                background    = args.background,
                width         = args.width,
                height        = args.height,
                dpi           = args.dpi,
                orient_nterm  = args.orient_nterm,
            )
            print(f"  [{i}/{total}] rendered  {cif_name}  →  {os.path.basename(job['out'])}", flush=True)
            ok += 1
        except Exception as exc:
            print(f"  [{i}/{total}] FAILED    {cif_name}: {exc}", flush=True)
            failed.append(job["cif"])

    print(f"\nDone: {ok}/{total} rendered.")
    if failed:
        sys.exit(1)


# ---------------------------------------------------------------------------
# StructureRenderer — importable class (jkoszucki env)
# ---------------------------------------------------------------------------

class StructureRenderer:
    """
    Render AF3 CIF structures to PNG via PyMOL.

    All rendering parameters are set at construction and applied to every
    render call. Already-rendered PNGs are skipped (checkpoint).

    Parameters
    ----------
    render_mode : "fast" | "accurate"
        "fast"     → no ray tracing; quick preview.
        "accurate" → ray tracing; publication quality (slow).
    color_mode : "plddt" | "domain" | "chain" | "gray"
        "plddt"   → colour by b-factor (pLDDT); standard AlphaFold scheme.
        "domain"  → gray base; named residue ranges coloured by domain_colors.
        "chain"   → each chain coloured distinctly.
        "gray"    → uniform gray throughout.
    domain_colors : list of (resi_selection, hex_color)
        Applied when color_mode="domain". resi_selection is PyMOL resi syntax,
        e.g. "1-350" or "1-50+100-150". Residues not listed keep gray_color.
    gray_color : str
        Hex color used as the base / unlabelled-residue color.
        Default: "#cfcfce".
    background : str
        PyMOL background color name or hex. Default: "white".
    width, height : int
        PNG dimensions in pixels. Default: 3000 × 2400.
    style : cfg.style
        Provides dpi. If None, defaults to 300.
    """

    _SELF = Path(__file__).resolve()

    def __init__(
        self,
        render_mode:   str                           = "fast",
        color_mode:    str                           = "plddt",
        domain_colors: list[tuple[str, str]] | None = None,
        gray_color:    str                           = _DEFAULT_GRAY,
        background:    str                           = _DEFAULT_BG,
        width:         int                           = _DEFAULT_WIDTH,
        height:        int                           = _DEFAULT_HEIGHT,
        orient_nterm:  bool                          = False,
        style                                        = None,
    ) -> None:
        if render_mode not in _RENDER_MODES:
            raise ValueError(f"render_mode must be one of {_RENDER_MODES}, got {render_mode!r}")
        if color_mode not in _COLOR_MODES:
            raise ValueError(f"color_mode must be one of {_COLOR_MODES}, got {color_mode!r}")

        self.render_mode   = render_mode
        self.color_mode    = color_mode
        self.domain_colors = list(domain_colors) if domain_colors else []
        self.gray_color    = gray_color
        self.background    = background
        self.width         = width
        self.height        = height
        self.orient_nterm  = orient_nterm
        self.dpi           = getattr(style, "dpi", 300)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    _SECS_PER_STRUCTURE = {"fast": 5, "accurate": 120}

    def render_all(self, output_dir: Path) -> None:
        """
        Discover every *.cif under output_dir and render it in-place.

        Output filename: {cif_stem}_{render_mode}.png (next to the CIF).
        Already-rendered PNGs are skipped (checkpoint).
        Prints a time estimate before starting and per-structure progress.
        """
        output_dir = Path(output_dir)
        _exclude = {"prev"}
        cif_paths  = [
            p for p in sorted(output_dir.rglob("*.cif"))
            if not any(part in _exclude for part in p.parts)
        ]
        print(f"  Found {len(cif_paths)} CIF file(s) under {output_dir}")

        jobs, skipped = [], []
        for cif in cif_paths:
            png = cif.parent / f"{cif.stem}_{self.render_mode}.png"
            if png.exists():
                skipped.append(cif)
                print(f"  [skip] {cif.name}  →  {png.name} already exists")
            else:
                jobs.append({"cif": str(cif), "out": str(png)})

        print(f"\n  Checkpoint: {len(skipped)} already rendered, {len(jobs)} queued.")

        if not jobs:
            print("  Nothing to render.")
            return

        secs = len(jobs) * self._SECS_PER_STRUCTURE.get(self.render_mode, 60)
        est  = f"~{secs // 60}min {secs % 60}s" if secs >= 60 else f"~{secs}s"
        print(f"  Estimated time ({self.render_mode} mode, ~2 min/structure): {est}")
        print(f"  Rendering {len(jobs)} structure(s) [{self.render_mode}, {self.color_mode}] …\n")

        self._dispatch(jobs)

    def render_next_to(self, cif_paths: list[Path]) -> None:
        """
        Render each CIF; write structure.png in the same directory as the CIF.

        Standard for data-hierarchy PNGs (next to structure.cif):
            gwas-data/{locus}/{ecod_type}/{cl}/{PC}/protein/structure.png
            acetyltransferase/{source}/{id}/structure.png
        """
        jobs, skipped, missing = self._build_jobs_next_to(cif_paths)
        self._report(jobs, skipped, missing)
        self._dispatch(jobs)

    def render_to_dir(self, cif_paths: list[Path], out_dir: Path) -> None:
        """
        Render each CIF; write {stem}.png into out_dir.

        For figure panel PNGs that land in plots/ rather than next to the CIF.
        """
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        jobs, skipped, missing = self._build_jobs_to_dir(cif_paths, out_dir)
        self._report(jobs, skipped, missing)
        self._dispatch(jobs)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _build_jobs_next_to(self, cif_paths):
        jobs, skipped, missing = [], [], []
        for cif in cif_paths:
            cif = Path(cif)
            if not cif.exists():
                missing.append(str(cif))
                continue
            png = cif.parent / "structure.png"
            if png.exists():
                skipped.append(cif.parent.name)
            else:
                jobs.append({"cif": str(cif), "out": str(png)})
        return jobs, skipped, missing

    def _build_jobs_to_dir(self, cif_paths, out_dir):
        jobs, skipped, missing = [], [], []
        for cif in cif_paths:
            cif = Path(cif)
            if not cif.exists():
                missing.append(str(cif))
                continue
            png = out_dir / f"{cif.stem}.png"
            if png.exists():
                skipped.append(cif.stem)
            else:
                jobs.append({"cif": str(cif), "out": str(png)})
        return jobs, skipped, missing

    def _dispatch(self, jobs: list[dict]) -> None:
        if not jobs:
            print("  Nothing to render.")
            return
        print(f"  Rendering {len(jobs)} structure(s) "
              f"[{self.render_mode}, {self.color_mode}] …")
        cmd = [
            "conda", "run", "--no-capture-output", "-n", "pymol", "python",
            str(self._SELF),
            "--jobs",          json.dumps(jobs),
            "--render-mode",   self.render_mode,
            "--color-mode",    self.color_mode,
            "--domain-colors", json.dumps(self.domain_colors),
            "--gray-color",    self.gray_color,
            "--background",    self.background,
            "--width",         str(self.width),
            "--height",        str(self.height),
            "--dpi",           str(self.dpi),
        ]
        if self.orient_nterm:
            cmd.append("--orient-nterm")
        result = subprocess.run(cmd, text=True, capture_output=False)
        if result.returncode != 0:
            print(f"  [warn] renderer exited with code {result.returncode}")

    @staticmethod
    def _report(jobs: list, skipped: list, missing: list) -> None:
        print(f"  Checkpoint: {len(skipped)} already rendered, {len(jobs)} queued.")
        for m in missing:
            print(f"  [warn] CIF not found: {m}")


def _scan_main() -> None:
    """
    Entry point for --scan mode (jkoszucki env).

    Discovers all *.cif files under cfg.output_dir and renders each one
    using 'accurate' mode, writing {cif_stem}_accurate.png next to the CIF.

    Usage:
        conda run -n jkoszucki python structure_renderer.py --scan
    """
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from config import Config  # noqa: PLC0415

    cfg = Config()
    renderer = StructureRenderer(render_mode="accurate", color_mode="plddt")
    renderer.render_all(cfg.output_dir)


if __name__ == "__main__":
    if "--scan" in sys.argv:
        _scan_main()
    else:
        _worker_main()
