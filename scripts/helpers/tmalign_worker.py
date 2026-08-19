"""
TM-align all-vs-all worker — runs in the `tmtools` conda environment.

That environment carries only tmtools, numpy and biopython, so this script sticks to
the standard library for I/O and does no analysis: it writes a long-format TSV of
pairwise TM-scores and nothing else. Callers in the `jkoszucki` environment invoke it
via subprocess and do their own clustering.

    conda run -n tmtools python scripts/helpers/tmalign_worker.py \
        --cif-list paths.txt --out tmscore.tsv

`--cif-list` is a text file with one CIF path per line; blank lines are ignored. Every
unordered pair is compared once.

Output columns
--------------
    id1, id2      — filename stems, with any trailing "_model_0" removed
    tm_norm_1     — TM-score normalised by the length of structure 1
    tm_norm_2     — TM-score normalised by the length of structure 2
    tm_max        — max(tm_norm_1, tm_norm_2)
    rmsd          — RMSD over aligned residues
"""

from __future__ import annotations

import argparse
import csv
import sys
from itertools import combinations
from pathlib import Path

from Bio.PDB import MMCIFParser
from tmtools import tm_align
from tmtools.io import get_residue_data

_MODEL_SUFFIX = "_model_0"


def _structure_id(cif_path: Path) -> str:
    stem = cif_path.stem
    return stem[: -len(_MODEL_SUFFIX)] if stem.endswith(_MODEL_SUFFIX) else stem


def _load(cif_path: Path):
    """Return (id, coords, seq) for the first chain of a CIF file."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(cif_path.stem, str(cif_path))
    chain = next(structure.get_chains())
    coords, seq = get_residue_data(chain)
    return _structure_id(cif_path), coords, seq


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cif-list", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    paths = [
        Path(line.strip())
        for line in args.cif_list.read_text().splitlines()
        if line.strip()
    ]
    missing = [p for p in paths if not p.is_file()]
    if missing:
        print(f"ERROR: {len(missing)} CIF path(s) do not exist, e.g. {missing[0]}", file=sys.stderr)
        return 1

    print(f"  loading {len(paths)} structures …", flush=True)
    loaded = [_load(p) for p in paths]

    pairs = list(combinations(range(len(loaded)), 2))
    print(f"  running TM-align on {len(pairs)} pairs …", flush=True)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["id1", "id2", "tm_norm_1", "tm_norm_2", "tm_max", "rmsd"])

        for n, (i, j) in enumerate(pairs, start=1):
            id1, coords1, seq1 = loaded[i]
            id2, coords2, seq2 = loaded[j]
            res = tm_align(coords1, coords2, seq1, seq2)
            tm1 = round(res.tm_norm_chain1, 4)
            tm2 = round(res.tm_norm_chain2, 4)
            writer.writerow([id1, id2, tm1, tm2, max(tm1, tm2), round(res.rmsd, 4)])

            if n % 250 == 0 or n == len(pairs):
                print(f"    {n}/{len(pairs)} pairs", flush=True)

    print(f"  → {args.out}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
