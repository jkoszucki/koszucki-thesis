"""
Run the full processing workflow, in dependency order, ending with the AF3 index rebuild.

  1. processing/gwas-proc/main.py
  2. processing/sgnh-proc/main.py
  3. processing/enzymes-proc/main.py
  4. processing/acetyl-proc/main.py
  5. processing/cps-proc/main.py
  6. helpers/build_af3_index.py   ← integral last step: rebuilds summary.tsv, 2_BEST_MODELS
                                     symlinks, missing/ batches, 4_ORGANISED symlinks, and
                                     renders any not-yet-rendered nterminal structures (pymol env)

Each step is run as a subprocess in the current interpreter (run this script itself under
the jkoszucki conda env). A step's failure stops the run — later steps assume earlier
output exists.

Usage:
    conda run -n jkoszucki python scripts/run_pipeline.py
    conda run -n jkoszucki python scripts/run_pipeline.py --skip-af3   # all processing, no AF3 index step
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent

STEPS = [
    ("gwas-proc",    ROOT / "processing" / "gwas-proc"   / "main.py"),
    ("sgnh-proc",    ROOT / "processing" / "sgnh-proc"    / "main.py"),
    ("enzymes-proc", ROOT / "processing" / "enzymes-proc" / "main.py"),
    ("acetyl-proc",  ROOT / "processing" / "acetyl-proc"  / "main.py"),
    ("cps-proc",     ROOT / "processing" / "cps-proc"     / "main.py"),
]

AF3_INDEX_STEP = ("af3-index", ROOT / "helpers" / "build_af3_index.py")


def _run(name: str, script: Path) -> None:
    print(f"\n{'=' * 70}\n{name}\n{'=' * 70}")
    result = subprocess.run([sys.executable, str(script)])
    if result.returncode != 0:
        print(f"\n[FAILED] {name} exited with code {result.returncode} — stopping.")
        sys.exit(result.returncode)


def main() -> None:
    skip_af3 = "--skip-af3" in sys.argv

    for name, script in STEPS:
        _run(name, script)

    if skip_af3:
        print("\n[skip] AF3 index rebuild (--skip-af3)")
    else:
        _run(*AF3_INDEX_STEP)

    print("\nPipeline complete.")


if __name__ == "__main__":
    main()
